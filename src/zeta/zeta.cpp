#include <silex/zeta.hpp>

#include "zeta_internal.hpp"

#include <flint/flint.h>

#include <silex/diagnostics.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/acb.hpp>
#include <silex/flint/dirichlet.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/flint/nmod_poly_factor.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/signature.hpp>
#include <silex/unit.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <utility>
#include <vector>

namespace silex {
namespace {

constexpr ulong kBfMinCutoff = UWORD(70);
constexpr ulong kBfDefaultMaxCutoff = UWORD(20000);
constexpr slong kBfMaxTargetBits = WORD(20);
constexpr slong kBfMaxWorkPrecision = WORD(2048);
constexpr slong kBfArrayLinearFactorMaxDegree = WORD(6);
constexpr slong kValidationInitialWorkPrecision = WORD(64);
constexpr slong kValidationReserveBits = WORD(20);

struct BfResidueDegreeContext {
    const fmpq_poly_struct* defining_polynomial = nullptr;
    std::vector<flint::Fmpz> defining_coefficients;
    slong defining_degree = -1;
    flint::Fmpz order_discriminant;
    flint::Fmpz equation_index;
    const FactorBase* residue_degree_base = nullptr;
    bool has_defining_polynomial_route = false;
    bool has_order_discriminant = false;
    bool check_equation_index = false;
    bool has_complete_residue_degree_base = false;
};

struct BfSummandScratch {
    flint::Arb log_norm;
    flint::Arb sqrt_norm_power;
    flint::Arb term1;
    flint::Arb term2;
};

struct BfPrimeTermScratch {
    flint::Arb term;
    flint::Arb weighted_term;
    BfSummandScratch summand;
};

bool bf_residue_degree_cache_lookup(
        std::vector<slong>& out,
        detail::ZetaBfResidueDegreeCache* cache,
        ulong p,
        const DiagnosticsContext* diagnostics) noexcept {
    if (cache == nullptr || cache->entries.empty()) {
        return false;
    }

    std::size_t& hint = cache->lookup_hint;
    if (hint < cache->entries.size() && cache->entries[hint].p == p) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.zeta_bf.term.prime_decomposition_type.cache_hit");
        const detail::ZetaBfResidueDegreeCacheEntry& entry =
                cache->entries[hint];
        out.assign(cache->residue_degrees.begin() +
                           static_cast<std::ptrdiff_t>(entry.offset),
                   cache->residue_degrees.begin() +
                           static_cast<std::ptrdiff_t>(
                                   entry.offset + entry.length));
        ++hint;
        return true;
    }
    if (p > cache->entries.back().p) {
        return false;
    }

    auto it = std::lower_bound(
            cache->entries.begin(), cache->entries.end(), p,
            [](const detail::ZetaBfResidueDegreeCacheEntry& entry,
               ulong prime) {
                return entry.p < prime;
            });
    if (it == cache->entries.end() || it->p != p) {
        return false;
    }

    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.cache_hit");
    out.assign(cache->residue_degrees.begin() +
                       static_cast<std::ptrdiff_t>(it->offset),
               cache->residue_degrees.begin() +
                       static_cast<std::ptrdiff_t>(it->offset + it->length));
    hint = static_cast<std::size_t>(it - cache->entries.begin()) + 1;
    return true;
}

detail::ZetaBfResidueDegreeCacheEntry
bf_residue_degree_cache_entry(detail::ZetaBfResidueDegreeCache* cache,
                              ulong p,
                              const std::vector<slong>& residue_degrees)
        noexcept {
    const std::size_t offset = cache->residue_degrees.size();
    cache->residue_degrees.insert(cache->residue_degrees.end(),
                                  residue_degrees.begin(),
                                  residue_degrees.end());
    return detail::ZetaBfResidueDegreeCacheEntry{
            p, offset, residue_degrees.size()};
}

void bf_residue_degree_cache_store(
        detail::ZetaBfResidueDegreeCache* cache,
        ulong p,
        const std::vector<slong>& residue_degrees,
        const DiagnosticsContext* diagnostics) noexcept {
    if (cache == nullptr || residue_degrees.empty()) {
        return;
    }
    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.cache_store");
    if (cache->entries.empty() || cache->entries.back().p < p) {
        cache->entries.push_back(
                bf_residue_degree_cache_entry(cache, p, residue_degrees));
        return;
    }

    auto it = std::lower_bound(
            cache->entries.begin(), cache->entries.end(), p,
            [](const detail::ZetaBfResidueDegreeCacheEntry& entry,
               ulong prime) {
                return entry.p < prime;
            });
    if (it != cache->entries.end() && it->p == p) {
        *it = bf_residue_degree_cache_entry(cache, p, residue_degrees);
        return;
    }
    cache->entries.insert(
            it, bf_residue_degree_cache_entry(cache, p, residue_degrees));
}

bool supported_order(const Order& order, slong precision) noexcept {
    return precision > 0 && order.has_basis() && order.is_maximal();
}

bool fmpq_poly_is_monic_integral(const fmpq_poly_t polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree < 1) {
        return false;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (fmpz_is_one(fmpq_denref(coeff.raw())) == 0) {
            return false;
        }
    }

    fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, degree);
    return fmpz_is_one(fmpq_numref(coeff.raw())) != 0;
}

bool dirichlet_value_matches_kronecker(const flint::DirichletGroup& group,
                                       const flint::DirichletChar& character,
                                       flint::FmpzConstRef discriminant,
                                       ulong n) noexcept {
    flint::Fmpz zn;
    flint::fmpz_set_ui(flint::FmpzRef(zn), n);
    const int expected = flint::fmpz_kronecker(discriminant,
                                              flint::FmpzConstRef(zn));
    if (expected == 0) {
        return true;
    }

    const ulong got = flint::dirichlet_chi(group, character, n);
    if (expected == 1) {
        return got == 0;
    }
    return expected == -1 && group.exponent() % 2 == 0 &&
           got == group.exponent() / 2;
}

bool find_quadratic_character(flint::DirichletChar& out,
                              const flint::DirichletGroup& group,
                              flint::FmpzConstRef discriminant) noexcept {
    flint::DirichletChar candidate(group);
    for (ulong j = 0; j < group.character_count(); ++j) {
        flint::dirichlet_char_index(candidate, group, j);
        if (!flint::dirichlet_char_is_real(group, candidate) ||
            !flint::dirichlet_char_is_primitive(group, candidate)) {
            continue;
        }

        bool ok = true;
        for (ulong n = 1; n <= group.modulus(); ++n) {
            if (!dirichlet_value_matches_kronecker(
                        group, candidate, discriminant, n)) {
                ok = false;
                break;
            }
        }
        if (ok) {
            flint::dirichlet_char_set(out, group, candidate);
            return true;
        }
    }
    return false;
}

bool quadratic_residue(flint::Arb& out,
                       const Order& order,
                       slong precision) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr ||
        field->backend_kind() != NumberFieldBackendKind::quadratic ||
        precision <= 0) {
        return false;
    }

    flint::Fmpz conductor;
    if (!order.quadratic_conductor(flint::FmpzRef(conductor)) ||
        !flint::fmpz_is_one(conductor)) {
        return false;
    }

    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    if (!flint::fmpz_abs_fits_ui(flint::FmpzConstRef(abs_discriminant))) {
        return false;
    }

    const ulong modulus =
            flint::fmpz_get_ui(flint::FmpzConstRef(abs_discriminant));
    if (modulus == 0) {
        return false;
    }

    flint::DirichletGroup group(modulus);
    if (!group.is_initialized()) {
        return false;
    }
    flint::DirichletChar character(group);
    if (!find_quadratic_character(character, group,
                                  flint::FmpzConstRef(discriminant))) {
        return false;
    }

    flint::Fmpq one;
    flint::Acb value;
    flint::fmpq_one(one);
    flint::acb_dirichlet_l_fmpq(
            value, flint::FmpqConstRef(one), group.raw(), character.raw(),
            FLINT_MAX(precision + 64, 128));
    if (!flint::acb_is_finite(value) ||
        !flint::arb_contains_zero(flint::acb_imag_part(value)) ||
        !flint::arb_is_positive(flint::acb_real_part(value))) {
        return false;
    }

    flint::arb_set(flint::ArbRef(out), flint::acb_real_part(value));
    return true;
}

bool class_regulator_product_from_residue(flint::ArbRef out,
                                          const Order& order,
                                          flint::ArbConstRef residue,
                                          slong precision,
                                          const DiagnosticsContext* diagnostics)
        noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.class_regulator_product_from_residue");
    if (!supported_order(order, precision) ||
        !flint::arb_is_finite(residue) ||
        !flint::arb_is_positive(residue)) {
        return false;
    }

    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    Signature sig;
    flint::Fmpz discriminant;
    flint::Fmpz roots;
    if (!sig.compute(*field) ||
        !root_of_unity_order(flint::FmpzRef(roots), *field) ||
        !order.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }

    flint::Fmpz abs_discriminant;
    flint::Arb product;
    flint::Arb denominator;
    flint::Arb pi;

    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    flint::arb_set_fmpz(product, flint::FmpzConstRef(abs_discriminant));
    flint::arb_sqrt(product, product, precision);
    flint::arb_mul(product, product, residue, precision);
    flint::arb_mul_fmpz(product, product, flint::FmpzConstRef(roots),
                        precision);

    flint::arb_one(denominator);
    flint::arb_mul_2exp_si(denominator, denominator, sig.r1());
    flint::arb_const_pi(pi, precision);
    flint::arb_mul_2exp_si(pi, pi, 1);
    for (slong i = 0; i < sig.r2(); ++i) {
        flint::arb_mul(denominator, denominator, pi, precision);
    }

    flint::arb_div(product, product, denominator, precision);
    if (!flint::arb_is_finite(product) ||
        !flint::arb_is_positive(product)) {
        return false;
    }

    flint::arb_set(out, flint::ArbConstRef(product));
    return true;
}

slong bf_target_bits(slong precision) noexcept {
    slong bits = precision / 32;
    if (bits < 2) {
        bits = 2;
    }
    if (bits > kBfMaxTargetBits) {
        bits = kBfMaxTargetBits;
    }
    return bits;
}

void bf_target(flint::Arb& out, slong precision) noexcept {
    flint::arb_one(out);
    flint::arb_mul_2exp_si(out, out, -bf_target_bits(precision));
}

bool bf_radius_lt(const flint::Arb& value, const flint::Arb& target) noexcept {
    flint::Arb radius;
    flint::arb_get_rad_arb(radius, value);
    return flint::arb_lt(radius, target);
}

ulong bf_adjust_cutoff(ulong cutoff) noexcept {
    while (cutoff % 9 != 0) {
        ++cutoff;
    }
    return cutoff;
}

bool bf_next_cutoff(ulong& next, ulong cutoff, ulong max_cutoff) noexcept {
    if (cutoff >= max_cutoff) {
        next = bf_adjust_cutoff(max_cutoff);
        return false;
    }

    if (cutoff > UWORD_MAX / 2) {
        return false;
    }

    ulong candidate = 2 * cutoff;
    if (candidate > max_cutoff) {
        candidate = max_cutoff;
    }

    next = bf_adjust_cutoff(candidate);
    return next > cutoff;
}

bool bf_power_ui_less_than(ulong& out,
                           ulong base,
                           slong exponent,
                           ulong cutoff) noexcept {
    if (base < 2 || exponent <= 0 || cutoff <= 2) {
        return false;
    }

    const ulong limit = cutoff - 1;
    ulong value = 1;
    for (slong i = 0; i < exponent; ++i) {
        if (value > limit / base) {
            return false;
        }
        value *= base;
    }

    out = value;
    return true;
}

void bf_summand_ui(flint::Arb& out,
                   BfSummandScratch& scratch,
                   ulong norm_power,
                   slong exponent,
                   const flint::Arb& sqrt_cutoff_log_cutoff,
                   slong precision) noexcept {
    flint::arb_sqrt_ui(scratch.sqrt_norm_power, norm_power, precision);

    flint::arb_div_ui(scratch.term1, sqrt_cutoff_log_cutoff, norm_power,
                      precision);
    flint::arb_div_ui(scratch.term1, scratch.term1,
                      static_cast<ulong>(exponent), precision);

    flint::arb_div(scratch.term2, scratch.log_norm,
                   flint::ArbConstRef(scratch.sqrt_norm_power), precision);
    flint::arb_sub(out, scratch.term1, scratch.term2, precision);
}

bool bf_linear_factor_count_from_squarefree_nmod(
        slong& out,
        const flint::NmodPoly& reduced,
        ulong modulus) noexcept {
    const slong degree = nmod_poly_degree(reduced.raw());
    if (degree < 1) {
        return false;
    }
    if (degree == 1) {
        out = 1;
        return true;
    }

    // FLINT 3.5 nmod_poly/powmod_binexp.c and nmod_poly/gcd.c expose the
    // low-level kernels used by their owning-polynomial wrappers.  For the
    // profiled low degrees, fixed scratch storage avoids wrapper allocation
    // while FLINT still computes x^p mod f and gcd(f, x^p - x).
    if (degree <= kBfArrayLinearFactorMaxDegree) {
        constexpr std::size_t kScratchLength =
                static_cast<std::size_t>(
                        kBfArrayLinearFactorMaxDegree + 1);
        std::array<ulong, kScratchLength> x{};
        std::array<ulong, kScratchLength> frobenius{};
        std::array<ulong, kScratchLength> roots{};
        x[1] = 1;
        _nmod_poly_powmod_ui_binexp(
                frobenius.data(), x.data(), modulus,
                reduced.raw()->coeffs, reduced.raw()->length,
                reduced.raw()->mod);
        frobenius[1] =
                nmod_sub(frobenius[1], UWORD(1), reduced.raw()->mod);

        slong frobenius_length = degree;
        while (frobenius_length > 0 &&
               frobenius[static_cast<std::size_t>(
                           frobenius_length - 1)] == 0) {
            --frobenius_length;
        }
        if (frobenius_length == 0) {
            out = degree;
            return true;
        }

        const slong gcd_length = _nmod_poly_gcd(
                roots.data(), reduced.raw()->coeffs, reduced.raw()->length,
                frobenius.data(), frobenius_length, reduced.raw()->mod);
        out = gcd_length - 1;
        return out >= 0 && out <= degree;
    }

    flint::NmodPoly reverse(modulus);
    flint::NmodPoly inverse(modulus);
    flint::NmodPoly frobenius(modulus);
    flint::NmodPoly x(modulus);
    flint::NmodPoly roots(modulus);

    nmod_poly_reverse(reverse.raw(), reduced.raw(), reduced.raw()->length);
    nmod_poly_inv_series(inverse.raw(), reverse.raw(), reduced.raw()->length);
    nmod_poly_powmod_x_ui_preinv(frobenius.raw(), modulus, reduced.raw(),
                                 inverse.raw());

    nmod_poly_set_coeff_ui(x.raw(), 1, 1);
    nmod_poly_sub(frobenius.raw(), frobenius.raw(), x.raw());
    nmod_poly_gcd(roots.raw(), frobenius.raw(), reduced.raw());

    out = nmod_poly_degree(roots.raw());
    return out >= 0 && out <= degree;
}

bool bf_residue_degrees_from_small_squarefree_nmod(
        std::vector<slong>& out,
        const flint::NmodPoly& reduced,
        ulong modulus,
        const DiagnosticsContext* diagnostics) noexcept {
    const slong degree = nmod_poly_degree(reduced.raw());
    if (degree < 1 || degree > 3) {
        return false;
    }

    out.clear();
    out.reserve(static_cast<std::size_t>(degree));
    if (degree == 1) {
        out.push_back(1);
        return true;
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.nmod."
            "small_degree_roots");
    slong root_count = -1;
    if (!bf_linear_factor_count_from_squarefree_nmod(
                root_count, reduced, modulus)) {
        return false;
    }
    if (degree == 2) {
        if (root_count == 0) {
            out.push_back(2);
            return true;
        }
        if (root_count == 2) {
            out.push_back(1);
            out.push_back(1);
            return true;
        }
        out.clear();
        return false;
    }

    if (root_count == 0) {
        out.push_back(3);
        return true;
    }
    if (root_count == 1) {
        out.push_back(1);
        out.push_back(2);
        return true;
    }
    if (root_count == 3) {
        out.push_back(1);
        out.push_back(1);
        out.push_back(1);
        return true;
    }

    out.clear();
    return false;
}

bool bf_reduce_defining_polynomial_nmod(
        flint::NmodPoly& out,
        const BfResidueDegreeContext& context,
        ulong modulus,
        const DiagnosticsContext* diagnostics) noexcept {
    if (context.defining_degree <= 0 ||
        context.defining_coefficients.size() !=
                static_cast<std::size_t>(context.defining_degree + 1) ||
        modulus < 2) {
        return false;
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.nmod."
            "reduce_coefficients");
    for (slong i = 0; i <= context.defining_degree; ++i) {
        nmod_poly_set_coeff_ui(
                out.raw(), i,
                fmpz_fdiv_ui(context.defining_coefficients[
                                     static_cast<std::size_t>(i)]
                                     .raw(),
                             modulus));
    }
    return nmod_poly_degree(out.raw()) > 0;
}

bool bf_residue_degrees_from_equation_order_nmod(
        std::vector<slong>& out,
        const Order& order,
        const BfResidueDegreeContext& context,
        std::vector<slong>& degree_scratch,
        ulong modulus,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.nmod");
    if (degree_scratch.size() <
                static_cast<std::size_t>(context.defining_degree) ||
        modulus < 2) {
        return false;
    }

    flint::NmodPoly reduced(modulus);
    if (!bf_reduce_defining_polynomial_nmod(reduced, context, modulus,
                                            diagnostics)) {
        return false;
    }

    if (bf_residue_degrees_from_small_squarefree_nmod(out, reduced, modulus,
                                                      diagnostics)) {
        return true;
    }

    flint::NmodPolyFactor factorization;
    slong* degree_data = degree_scratch.data();
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.zeta_bf.term.prime_decomposition_type.nmod."
                "factor_distinct_deg");
        nmod_poly_factor_distinct_deg(factorization.raw(), reduced.raw(),
                                      &degree_data);
    }
    if (factorization.raw()->num <= 0) {
        return false;
    }

    out.clear();
    out.reserve(static_cast<std::size_t>(order.degree()));
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        const slong degree = degree_scratch[static_cast<std::size_t>(i)];
        const slong product_degree =
                nmod_poly_degree(factorization.raw()->p + i);
        if (degree <= 0 || product_degree <= 0 ||
            product_degree % degree != 0) {
            out.clear();
            return false;
        }

        const slong count = product_degree / degree;
        for (slong j = 0; j < count; ++j) {
            out.push_back(degree);
        }
    }

    return !out.empty();
}

bool bf_init_residue_degree_context(
        BfResidueDegreeContext& out,
        const Order& order,
        const FactorBase* residue_degree_base = nullptr) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr || !order.is_maximal() ||
        field->raw_flint_field() == nullptr ||
        !fmpq_poly_is_monic_integral(field->raw_flint_field()->pol)) {
        return false;
    }

    BfResidueDegreeContext next;
    next.defining_polynomial = field->raw_flint_field()->pol;
    next.defining_degree = fmpq_poly_degree(next.defining_polynomial);
    next.defining_coefficients.reserve(
            static_cast<std::size_t>(next.defining_degree + 1));
    for (slong i = 0; i <= next.defining_degree; ++i) {
        flint::Fmpz coeff;
        fmpq_poly_get_coeff_fmpz(coeff.raw(), next.defining_polynomial, i);
        next.defining_coefficients.push_back(std::move(coeff));
    }
    if (!order.discriminant(flint::FmpzRef(next.order_discriminant))) {
        return false;
    }
    next.has_order_discriminant = true;
    next.has_defining_polynomial_route = true;

    if (!order.is_equation_order()) {
        Order equation = Order::equation_order(*field);
        if (!equation.is_defined() ||
            !order_index(flint::FmpzRef(next.equation_index), equation,
                         order)) {
            return false;
        }
        next.check_equation_index = true;
    }

    if (residue_degree_base != nullptr &&
        residue_degree_base->is_defined() &&
        residue_degree_base->rational_prime_blocks_are_complete() &&
        same_order_parent(residue_degree_base->parent(), &order)) {
        next.residue_degree_base = residue_degree_base;
        next.has_complete_residue_degree_base = true;
    }

    out = std::move(next);
    return true;
}

bool bf_residue_degree_context_uses_defining_polynomial(
        const BfResidueDegreeContext& context,
        ulong p) noexcept {
    if (!context.has_defining_polynomial_route ||
        context.defining_polynomial == nullptr) {
        return false;
    }
    if (context.has_order_discriminant &&
        fmpz_fdiv_ui(context.order_discriminant.raw(), p) == 0) {
        return false;
    }
    return !context.check_equation_index ||
           fmpz_fdiv_ui(context.equation_index.raw(), p) != 0;
}

bool bf_residue_degrees_from_defining_polynomial_ui(
        std::vector<slong>& out,
        const BfResidueDegreeContext& context,
        const Order& order,
        std::vector<slong>& degree_scratch,
        ulong p,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!bf_residue_degree_context_uses_defining_polynomial(context, p)) {
        return false;
    }

    return bf_residue_degrees_from_equation_order_nmod(
            out, order, context, degree_scratch, p, diagnostics);
}

bool bf_linear_residue_degree_count_from_defining_polynomial_ui(
        slong& out,
        const BfResidueDegreeContext& context,
        ulong p,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!bf_residue_degree_context_uses_defining_polynomial(context, p)) {
        return false;
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.nmod");
    flint::NmodPoly reduced(p);
    if (!bf_reduce_defining_polynomial_nmod(reduced, context, p,
                                            diagnostics)) {
        return false;
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.nmod."
            "cutoff_linear_factors");
    return bf_linear_factor_count_from_squarefree_nmod(out, reduced, p);
}

bool bf_residue_degrees_from_complete_factor_base(
        std::vector<slong>& out,
        const BfResidueDegreeContext& context,
        flint::FmpzConstRef p,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!context.has_complete_residue_degree_base ||
        context.residue_degree_base == nullptr) {
        return false;
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.term.prime_decomposition_type.factor_base");
    const FactorBase& base = *context.residue_degree_base;
    out.clear();
    for (slong block = 0; block < base.rational_prime_block_count();
         ++block) {
        flint::Fmpz rational_prime;
        slong length = 0;
        if (!base.rational_prime_block_data(
                    flint::FmpzRef(rational_prime), length, block)) {
            return false;
        }
        const int cmp = flint::fmpz_cmp(flint::FmpzConstRef(rational_prime),
                                        p);
        if (cmp > 0) {
            return false;
        }
        if (cmp != 0) {
            continue;
        }

        out.reserve(static_cast<std::size_t>(length));
        for (slong offset = 0; offset < length; ++offset) {
            slong index = -1;
            if (!base.rational_prime_block_index(index, block, offset)) {
                out.clear();
                return false;
            }
            const PrimeIdeal* prime = base.prime_at(index);
            if (prime == nullptr || prime->residue_degree() <= 0) {
                out.clear();
                return false;
            }
            out.push_back(prime->residue_degree());
        }
        return !out.empty();
    }
    return false;
}

bool bf_prime_residue_degrees_ui(std::vector<slong>& out,
                                 const BfResidueDegreeContext& context,
                                 const Order& order,
                                 std::vector<slong>& degree_scratch,
                                 ulong p,
                                 const DiagnosticsContext* diagnostics,
                                 detail::ZetaBfResidueDegreeCache*
                                         residue_degree_cache = nullptr)
        noexcept {
    if (p < 2) {
        return false;
    }
    if (bf_residue_degree_cache_lookup(out, residue_degree_cache, p,
                                       diagnostics)) {
        return true;
    }

    bool found = bf_residue_degrees_from_defining_polynomial_ui(
            out, context, order, degree_scratch, p, diagnostics);
    if (!found) {
        flint::Fmpz p_value;
        flint::fmpz_set_ui(flint::FmpzRef(p_value), p);
        found = bf_residue_degrees_from_complete_factor_base(
                out, context, flint::FmpzConstRef(p_value), diagnostics);
        if (!found) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.zeta_bf.term.prime_decomposition_type."
                    "prime_ideal");
            PrimeIdealList primes;
            if (!decompose_prime(primes, order,
                                 flint::FmpzConstRef(p_value))) {
                return false;
            }

            out.clear();
            out.reserve(static_cast<std::size_t>(primes.size()));
            for (slong i = 0; i < primes.size(); ++i) {
                const PrimeIdeal* prime = primes.at(i);
                if (prime == nullptr || prime->residue_degree() <= 0) {
                    out.clear();
                    return false;
                }
                out.push_back(prime->residue_degree());
            }
            found = true;
        }
    }

    if (!found || out.empty()) {
        return false;
    }
    bf_residue_degree_cache_store(residue_degree_cache, p, out, diagnostics);
    return true;
}

constexpr slong kFactorBaseBoundPrecision = WORD(64);

bool factor_base_bound_bdf_right_side(
        flint::Arb& out,
        const Order& order,
        const BfResidueDegreeContext& residue_context,
        std::vector<slong>& residue_degrees,
        std::vector<slong>& degree_scratch,
        detail::ZetaBfResidueDegreeCache& residue_degree_cache,
        slong degree,
        slong real_places,
        double x0,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!std::isfinite(x0) || x0 < 2.0 ||
        x0 >= static_cast<double>(WORD_MAX)) {
        return false;
    }
    const double cutoff_double = std::ceil(x0);
    if (!std::isfinite(cutoff_double) || cutoff_double < 2.0 ||
        cutoff_double >= static_cast<double>(WORD_MAX)) {
        return false;
    }
    const ulong cutoff = static_cast<ulong>(cutoff_double);

    flint::Arb log_cutoff;
    flint::Arb sum;
    flint::Arb log_norm;
    flint::Arb sqrt_norm_power;
    flint::Arb scaled_log;
    flint::Arb ratio;
    flint::Arb weight;
    flint::Arb term;
    flint::arb_log_ui(log_cutoff, cutoff, kFactorBaseBoundPrecision);
    flint::arb_zero(flint::ArbRef(sum));

    for (ulong p = 2; p < cutoff; p = n_nextprime(p, 1)) {
        if (!bf_prime_residue_degrees_ui(
                    residue_degrees, residue_context, order, degree_scratch,
                    p, diagnostics, &residue_degree_cache)) {
            return false;
        }
        for (slong residue_degree : residue_degrees) {
            ulong norm = 0;
            if (!bf_power_ui_less_than(norm, p, residue_degree, cutoff)) {
                continue;
            }

            flint::arb_log_ui(log_norm, norm,
                              kFactorBaseBoundPrecision);
            ulong norm_power = 1;
            for (slong exponent = 1;
                 norm_power <= (cutoff - 1) / norm;
                 ++exponent) {
                norm_power *= norm;
                flint::arb_sqrt_ui(sqrt_norm_power, norm_power,
                                   kFactorBaseBoundPrecision);
                flint::arb_mul_ui(scaled_log, log_norm,
                                  static_cast<ulong>(exponent),
                                  kFactorBaseBoundPrecision);
                flint::arb_div(ratio, scaled_log, log_cutoff,
                               kFactorBaseBoundPrecision);
                flint::arb_one(weight);
                flint::arb_sub(weight, weight, ratio,
                               kFactorBaseBoundPrecision);
                flint::arb_div(term, log_norm, sqrt_norm_power,
                               kFactorBaseBoundPrecision);
                flint::arb_mul(term, term, weight,
                               kFactorBaseBoundPrecision);
                flint::arb_add(sum, sum, term,
                               kFactorBaseBoundPrecision);
            }
        }
    }

    flint::Arb pi;
    flint::Arb pi_squared;
    flint::Arb catalan;
    flint::Arb archimedean;
    flint::Arb doubled_sum;
    flint::arb_const_pi(pi, kFactorBaseBoundPrecision);
    flint::arb_sqr(pi_squared, pi, kFactorBaseBoundPrecision);
    flint::arb_mul_ui(pi_squared, pi_squared,
                      static_cast<ulong>(degree),
                      kFactorBaseBoundPrecision);
    flint::arb_div_ui(pi_squared, pi_squared, UWORD(2),
                      kFactorBaseBoundPrecision);
    ::arb_const_catalan(catalan.raw(), kFactorBaseBoundPrecision);
    flint::arb_mul_ui(catalan, catalan,
                      static_cast<ulong>(4 * real_places),
                      kFactorBaseBoundPrecision);
    flint::arb_add(archimedean, pi_squared, catalan,
                   kFactorBaseBoundPrecision);
    flint::arb_div(archimedean, archimedean, log_cutoff,
                   kFactorBaseBoundPrecision);
    flint::arb_mul_2exp_si(doubled_sum, sum, 1);
    flint::arb_sub(out, doubled_sum, archimedean,
                   kFactorBaseBoundPrecision);
    return flint::arb_is_finite(out);
}

bool factor_base_bound_bdf(flint::FmpzRef out,
                                 const Order& order,
                                 const DiagnosticsContext* diagnostics)
        noexcept {
    if (out.raw() == nullptr || !order.is_maximal() || order.degree() <= 0) {
        return false;
    }
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    Signature signature;
    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    BfResidueDegreeContext residue_context;
    if (!signature.compute(*field) ||
        !order.discriminant(flint::FmpzRef(discriminant)) ||
        !bf_init_residue_degree_context(residue_context, order)) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    if (flint::fmpz_sgn(flint::FmpzConstRef(abs_discriminant)) <= 0) {
        return false;
    }

    flint::Arb target;
    flint::Arb log_discriminant;
    flint::Arb euler;
    flint::Arb pi;
    flint::Arb log_eight_pi;
    flint::Arb degree_term;
    flint::Arb real_term;
    flint::arb_log_fmpz(log_discriminant,
                        flint::FmpzConstRef(abs_discriminant),
                        kFactorBaseBoundPrecision);
    ::arb_const_euler(euler.raw(), kFactorBaseBoundPrecision);
    flint::arb_const_pi(pi, kFactorBaseBoundPrecision);
    flint::arb_mul_ui(log_eight_pi, pi, UWORD(8),
                      kFactorBaseBoundPrecision);
    flint::arb_log(log_eight_pi, log_eight_pi,
                   kFactorBaseBoundPrecision);
    flint::arb_add(degree_term, euler, log_eight_pi,
                   kFactorBaseBoundPrecision);
    flint::arb_mul_ui(degree_term, degree_term,
                      static_cast<ulong>(order.degree()),
                      kFactorBaseBoundPrecision);
    flint::arb_mul_ui(real_term, pi,
                      static_cast<ulong>(signature.r1()),
                      kFactorBaseBoundPrecision);
    flint::arb_div_ui(real_term, real_term, UWORD(2),
                      kFactorBaseBoundPrecision);
    flint::arb_sub(target, log_discriminant, degree_term,
                   kFactorBaseBoundPrecision);
    flint::arb_sub(target, target, real_term,
                   kFactorBaseBoundPrecision);
    if (!flint::arb_is_finite(target)) {
        return false;
    }

    std::vector<slong> residue_degrees;
    residue_degrees.reserve(static_cast<std::size_t>(order.degree()));
    std::vector<slong> degree_scratch(static_cast<std::size_t>(
            FLINT_MAX(order.degree(), WORD(1))));
    detail::ZetaBfResidueDegreeCache residue_degree_cache;
    flint::Arb right_side;
    double x0 = 100.0;
    double x1 = 2.0 * x0;
    if (!factor_base_bound_bdf_right_side(
                right_side, order, residue_context, residue_degrees,
                degree_scratch, residue_degree_cache, order.degree(),
                signature.r1(), x0, diagnostics)) {
        return false;
    }

    while (flint::arb_lt(right_side, target)) {
        x0 = x1;
        x1 = 2.0 * x0;
        if (!std::isfinite(x1) ||
            !factor_base_bound_bdf_right_side(
                    right_side, order, residue_context, residue_degrees,
                    degree_scratch, residue_degree_cache, order.degree(),
                    signature.r1(), x1, diagnostics)) {
            return false;
        }
    }
    if (!flint::arb_gt(right_side, target)) {
        return false;
    }

    double distance = std::fabs(x0 - x1);
    while (!(flint::arb_gt(right_side, target) && distance < 100.0)) {
        if (flint::arb_lt(right_side, target)) {
            x1 = x0 + 3.0 * distance / 2.0;
        } else if (flint::arb_gt(right_side, target)) {
            x1 = x0 - distance / 2.0;
        } else {
            return false;
        }
        distance = std::fabs(x1 - x0);
        x0 = x1;
        if (!factor_base_bound_bdf_right_side(
                    right_side, order, residue_context, residue_degrees,
                    degree_scratch, residue_degree_cache, order.degree(),
                    signature.r1(), x0, diagnostics)) {
            return false;
        }
    }

    const double result = std::ceil(x0);
    if (!std::isfinite(result) || result < 0.0 ||
        result >= static_cast<double>(WORD_MAX)) {
        return false;
    }
    flint::fmpz_set_ui(out, static_cast<ulong>(result));
    return true;
}

bool factor_base_bound_bach(flint::FmpzRef out,
                                  const Order& order) noexcept {
    if (out.raw() == nullptr || !order.is_maximal() || order.degree() <= 0) {
        return false;
    }
    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    if (flint::fmpz_sgn(flint::FmpzConstRef(abs_discriminant)) <= 0) {
        return false;
    }

    flint::Arb bound;
    flint::Arb rounded;
    flint::arb_log_fmpz(bound, flint::FmpzConstRef(abs_discriminant),
                        kFactorBaseBoundPrecision);
    flint::arb_sqr(bound, bound, kFactorBaseBoundPrecision);
    flint::arb_mul_ui(bound, bound,
                      order.degree() == 2 ? UWORD(6) : UWORD(12),
                      kFactorBaseBoundPrecision);
    ::arb_ceil(rounded.raw(), bound.raw(), kFactorBaseBoundPrecision);
    return ::arb_get_unique_fmpz(out.raw(), rounded.raw()) != 0;
}

bool grh_factor_base_bound_impl(
        flint::FmpzRef out,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    flint::Fmpz bdf;
    flint::Fmpz bach;
    if (!factor_base_bound_bdf(flint::FmpzRef(bdf), order,
                                     diagnostics) ||
        !factor_base_bound_bach(flint::FmpzRef(bach), order)) {
        return false;
    }
    if (flint::fmpz_cmp(flint::FmpzConstRef(bdf),
                        flint::FmpzConstRef(bach)) <= 0) {
        flint::fmpz_set(out, flint::FmpzConstRef(bdf));
    } else {
        flint::fmpz_set(out, flint::FmpzConstRef(bach));
    }
    if (flint::fmpz_is_zero(out)) {
        flint::fmpz_one(out);
    }
    return flint::fmpz_sgn(out) > 0;
}

double tail_residue_backtransform(long R1,
                        long R2,
                        double rK,
                        long C,
                        double C2,
                        double C3,
                        double r1K,
                        double r2K,
                        double logC,
                        double logC2,
                        double logC3) noexcept {
    const double rQ = 1.83787706641;
    const double r1Q = 1.98505372441;
    const double r2Q = 1.07991541347;
    return std::fabs(
            (R1 + R2 - 1) *
                    (12 * logC3 + 4 * logC2 - 9 * logC - 6) /
                    (2 * C * logC3) +
            (rK - rQ) * (6 * logC2 + 5 * logC + 2) / (C * logC3) -
            R2 * (6 * logC2 + 11 * logC + 6) / (C2 * logC2) -
            2 * (r1K - r1Q) * (3 * logC2 + 4 * logC + 2) /
                    (C2 * logC3) +
            (R1 + R2 - 1) *
                    (12 * logC3 + 40 * logC2 + 45 * logC + 18) /
                    (6 * C3 * logC3) +
            (r2K - r2Q) * (2 * logC2 + 3 * logC + 2) /
                    (C3 * logC3));
}

double tail_residue(long R1,
                    long R2,
                    double al2K,
                    double rKm,
                    double rKM,
                    double r1Km,
                    double r1KM,
                    double r2Km,
                    double r2KM,
                    double C,
                    long i) noexcept {
    static constexpr double tab[] = {
            0.50409264803,    0.26205336997,    0.14815491171,
            0.08770540561,    0.05347651832,    0.03328934284,
            0.02104510690,    0.01346475900,    0.00869778586,
            0.00566279855,    0.00371111950,    0.00244567837,
            0.00161948049,    0.00107686891,    0.00071868750,
            0.00048119961,    0.00032312188,    0.00021753772,
            0.00014679818,    9.9272855581E-5,  6.7263969995E-5,
            4.5656812967E-5,  3.1041124593E-5,  2.1136011590E-5,
            1.4411645381E-5,  9.8393304088E-6,  6.7257395409E-6,
            4.6025878272E-6,  3.1529719271E-6,  2.1620490021E-6,
            1.4839266071E-6};
    const double logC = std::log(C);
    const double logC2 = logC * logC;
    const double logC3 = logC * logC2;
    const double C2 = C * C;
    const double C3 = C * C2;
    const double E1 = i > 30 ? 0 : tab[i];
    // Preserve the source analytic-tail call order and C conversions.
    const double back1 = tail_residue_backtransform(
            static_cast<long>(rKm), static_cast<long>(r1KM), r2Km,
            static_cast<long>(C), C2, C3, static_cast<double>(R1),
            static_cast<double>(R2), logC, logC2, logC3);
    const double back2 = tail_residue_backtransform(
            static_cast<long>(rKM), static_cast<long>(r1Km), r2KM,
            static_cast<long>(C), C2, C3, static_cast<double>(R1),
            static_cast<double>(R2), logC, logC2, logC3);
    return al2K * ((33 * logC2 + 22 * logC + 8) /
                           (8 * logC3 * std::sqrt(C)) +
                   15 * E1 / 16) +
           std::max(back1, back2) / 2 +
           ((R1 + R2 - 1) * 4 * C + R2) * (C2 + 6 * logC) /
                   (4 * C2 * C2 * logC2);
}

long prime_needed(long N, long R1, long R2, double logD) noexcept {
    const double lim = 0.25;
    const double al2K = 0.3526 * logD - 0.8212 * N + 4.5007;
    const double rKm = -1.0155 * logD + 2.1042 * N - 8.3419;
    const double rKM = -0.5 * logD + 1.2076 * N + 1;
    const double r1Km = -logD + 1.4150 * N;
    const double r1KM = -logD + 1.9851 * N;
    const double r2Km = -logD + 0.9151 * N;
    const double r2KM = -logD + 1.0800 * N;
    long Cmin = 3;
    long Cmax = 3;
    long i = 0;
    while (tail_residue(R1, R2, al2K, rKm, rKM, r1Km, r1KM, r2Km,
                        r2KM, Cmax, i) > lim) {
        Cmin = Cmax;
        if (Cmax > (LONG_MAX / 2)) {
            return 0;
        }
        Cmax *= 2;
        ++i;
    }
    --i;
    while (Cmax - Cmin > 1) {
        const long t = (Cmin + Cmax) / 2;
        if (tail_residue(R1, R2, al2K, rKm, rKM, r1Km, r1KM, r2Km,
                         r2KM, t, i) > lim) {
            Cmin = t;
        } else {
            Cmax = t;
        }
    }
    return Cmax;
}

bool log_inverse_residue(double& out,
                              const Order& order,
                              ulong limc,
                              const DiagnosticsContext* diagnostics,
                              detail::ZetaBfResidueDegreeCache*
                                      residue_degree_cache = nullptr)
        noexcept {
    if (limc < 3) {
        return false;
    }

    const double log_limc = std::log(static_cast<double>(limc));
    const double log_limc2 = log_limc * log_limc;
    if (!std::isfinite(log_limc) || log_limc <= 0 || log_limc2 <= 0) {
        return false;
    }

    double denc = 1 / (std::pow(static_cast<double>(limc), 3.0) *
                       log_limc * log_limc2);
    const double c2 = (log_limc2 + 3 * log_limc / 2 + 1) * denc;
    denc *= limc;
    const double c1 = (3 * log_limc2 + 4 * log_limc + 2) * denc;
    denc *= limc;
    const double c0 = (3 * log_limc2 + 5 * log_limc / 2 + 1) * denc;
    if (!std::isfinite(c0) || !std::isfinite(c1) || !std::isfinite(c2)) {
        return false;
    }

    BfResidueDegreeContext residue_context;
    if (!bf_init_residue_degree_context(residue_context, order)) {
        return false;
    }

    std::vector<slong> residue_degrees;
    residue_degrees.reserve(static_cast<std::size_t>(order.degree()));
    std::vector<slong> degree_scratch(static_cast<std::size_t>(
            FLINT_MAX(order.degree(), WORD(1))));
    std::vector<slong> degree_counts(
            static_cast<std::size_t>(order.degree() + 1));

    double log_inverse_residue_sum = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.inverse_residue.prime_loop");
        for (ulong p = 2;; p = n_nextprime(p, 1)) {
            const double p_double = static_cast<double>(p);
            const double logp = std::log(p_double);
            const long limp = static_cast<long>(log_limc / logp);
            if (limp < 1) {
                break;
            }

            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.inverse_residue.prime_decomposition_type");
                if (!bf_prime_residue_degrees_ui(
                            residue_degrees, residue_context, order,
                            degree_scratch, p, diagnostics,
                            residue_degree_cache)) {
                    return false;
                }
            }

            std::fill(degree_counts.begin(), degree_counts.end(), WORD(0));
            for (slong degree : residue_degrees) {
                if (degree <= 0 ||
                    degree >= static_cast<slong>(degree_counts.size())) {
                    return false;
                }
                ++degree_counts[static_cast<std::size_t>(degree)];
            }

            log_inverse_residue_sum += 1 / p_double;
            double NPk = p_double;
            for (long k = 2; k <= limp; ++k) {
                NPk *= p_double;
                log_inverse_residue_sum += 1 / (k * NPk);
            }

            long addpsi = limp;
            const double p2_double = p_double * p_double;
            double addpsi1 =
                    p_double * (std::pow(p_double, limp) - 1) /
                    (p_double - 1);
            double addpsi2 =
                    p2_double * (std::pow(p2_double, limp) - 1) /
                    (p2_double - 1);
            for (slong degree = 1;
                 degree < static_cast<slong>(degree_counts.size());
                 ++degree) {
                const slong count =
                        degree_counts[static_cast<std::size_t>(degree)];
                if (count == 0 || degree > limp) {
                    continue;
                }

                const double NP = std::pow(p_double, degree);
                double inverse_residue_term = 1 / NP;
                const long kmax = limp / degree;
                double NPk_degree = NP;
                for (long k = 2; k <= kmax; ++k) {
                    NPk_degree *= NP;
                    inverse_residue_term += 1 / (k * NPk_degree);
                }

                const double NP2 = NP * NP;
                log_inverse_residue_sum -= count * inverse_residue_term;
                addpsi -= static_cast<long>(count * degree * kmax);
                addpsi1 -= count *
                           (degree * NP *
                            (std::pow(NP, static_cast<double>(kmax)) - 1) /
                            (NP - 1));
                addpsi2 -= count *
                           (degree * NP2 *
                            (std::pow(NP2, static_cast<double>(kmax)) - 1) /
                            (NP2 - 1));
            }

            log_inverse_residue_sum -= (addpsi * c0 - addpsi1 * c1 + addpsi2 * c2) * logp;
            if (!std::isfinite(log_inverse_residue_sum)) {
                return false;
            }
        }
    }

    out = log_inverse_residue_sum;
    return std::isfinite(out);
}

bool class_regulator_product_estimate_impl(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        detail::ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        flint::Fmpz* torsion_order = nullptr,
        Element* torsion_generator = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.inverse_residue.class_regulator_product_estimate");
    if (!supported_order(order, precision) ||
        ((torsion_order == nullptr) != (torsion_generator == nullptr))) {
        return false;
    }

    if (order.degree() == 1) {
        flint::arb_set_si(out, 1);
        return true;
    }

    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    Signature sig;
    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    flint::Fmpz roots;
    Element roots_generator(*field);
    // The source relation-completion setup computes nfrootsof1 once before
    // compute_inverse_residue and carries the result to final unit publication.
    if (!sig.compute(*field) ||
        !roots_generator.is_defined() ||
        !roots_of_unity(flint::FmpzRef(roots), roots_generator, *field) ||
        !order.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }

    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    const double absD =
            flint::fmpz_get_d(flint::FmpzConstRef(abs_discriminant));
    const double roots_count = flint::fmpz_get_d(flint::FmpzConstRef(roots));
    if (!std::isfinite(absD) || absD <= 0 || !std::isfinite(roots_count) ||
        roots_count <= 0) {
        return false;
    }

    const double logD = std::log(absD);
    const long limres = prime_needed(
            static_cast<long>(order.degree()), static_cast<long>(sig.r1()),
            static_cast<long>(sig.r2()), logD);
    if (limres < 3) {
        return false;
    }

    double log_inverse_residue_sum = 0;
    if (!log_inverse_residue(log_inverse_residue_sum, order,
                                  static_cast<ulong>(limres), diagnostics,
                                  residue_degree_cache)) {
        return false;
    }

    static constexpr double kPi = 3.141592653589793238462643383279502884;
    const double log_invhr =
            (sig.r1() + sig.r2()) * std::log(2.0) +
            sig.r2() * std::log(kPi) - 0.5 * logD -
            std::log(roots_count) + log_inverse_residue_sum;
    const double product = std::exp(-log_invhr);
    if (!std::isfinite(product) || product <= 0) {
        return false;
    }

    flint::arb_set_d(out, product);
    if (!flint::arb_is_finite(flint::ArbConstRef(out.raw())) ||
        !flint::arb_is_positive(flint::ArbConstRef(out.raw()))) {
        return false;
    }
    if (torsion_generator != nullptr) {
        if ((!torsion_generator->has_parent(*field) &&
             !torsion_generator->define(*field)) ||
            !torsion_generator->set(roots_generator)) {
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(*torsion_order),
                        flint::FmpzConstRef(roots));
    }
    return true;
}

void bf_reset_residue_degree_coefficients(
        std::vector<slong>& coefficients) noexcept {
    for (slong& coefficient : coefficients) {
        coefficient = 0;
    }
}

bool bf_accumulate_weighted_prime_powers(
        flint::Arb& sum,
        flint::Arb& term,
        flint::Arb& weighted_term,
        BfSummandScratch& summand_scratch,
        ulong rational_prime,
        const std::vector<slong>& degree_coefficients,
        slong coefficient_sign,
        ulong cutoff,
        const flint::Arb& sqrt_cutoff_log_cutoff,
        slong precision) noexcept {
    if (coefficient_sign != 1 && coefficient_sign != -1) {
        return false;
    }

    for (slong degree = 1;
         degree < static_cast<slong>(degree_coefficients.size()); ++degree) {
        const slong coefficient =
                coefficient_sign *
                degree_coefficients[static_cast<std::size_t>(degree)];
        if (coefficient == 0) {
            continue;
        }

        ulong norm = 0;
        if (!bf_power_ui_less_than(norm, rational_prime, degree, cutoff)) {
            continue;
        }

        flint::arb_log_ui(summand_scratch.log_norm, norm, precision);
        const slong abs_coefficient =
                coefficient < 0 ? -coefficient : coefficient;
        ulong norm_power = norm;
        const ulong limit = cutoff - 1;
        for (slong exponent = 1;; ++exponent) {
            bf_summand_ui(term, summand_scratch, norm_power, exponent,
                          sqrt_cutoff_log_cutoff, precision);

            const flint::Arb* addend = &term;
            if (abs_coefficient != 1) {
                flint::arb_mul_ui(weighted_term, term,
                                  static_cast<ulong>(abs_coefficient),
                                  precision);
                addend = &weighted_term;
            }

            if (coefficient > 0) {
                flint::arb_add(sum, sum, *addend, precision);
            } else {
                flint::arb_sub(sum, sum, *addend, precision);
            }

            if (norm_power > limit / norm) {
                break;
            }
            norm_power *= norm;
        }
    }
    return true;
}

bool bf_add_prime_terms(flint::Arb& sum,
                        const Order& order,
                        const BfResidueDegreeContext& residue_context,
                        std::vector<slong>& residue_degrees,
                        std::vector<slong>& degree_scratch,
                        std::vector<slong>& degree_coefficients,
                        BfPrimeTermScratch& scratch,
                        ulong rational_prime,
                        ulong cutoff,
                        ulong cutoff_div_9,
                        const flint::Arb& sqrt_cutoff_log_cutoff,
                        const flint::Arb& sqrt_cutoff9_log_cutoff9,
                        slong precision,
                        const DiagnosticsContext* diagnostics,
                        detail::ZetaBfResidueDegreeCache*
                                residue_degree_cache = nullptr) noexcept {
    slong linear_factor_count = -1;
    bool have_complete_residue_degrees = false;
    bool have_linear_factor_count = false;

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.zeta_bf.term.prime_decomposition_type");
        const bool only_degree_one_contributes =
                cutoff > 1 &&
                rational_prime > (cutoff - 1) / rational_prime;
        if (only_degree_one_contributes) {
            have_complete_residue_degrees =
                    bf_residue_degree_cache_lookup(
                            residue_degrees, residue_degree_cache,
                            rational_prime, diagnostics);
            if (!have_complete_residue_degrees) {
                have_linear_factor_count =
                        bf_linear_residue_degree_count_from_defining_polynomial_ui(
                                linear_factor_count, residue_context,
                                rational_prime, diagnostics);
            }
        }
        if (!have_complete_residue_degrees &&
            !have_linear_factor_count) {
            if (!bf_prime_residue_degrees_ui(
                        residue_degrees, residue_context, order,
                        degree_scratch, rational_prime, diagnostics,
                        residue_degree_cache)) {
                return false;
            }
            have_complete_residue_degrees = true;
        }
    }

    bf_reset_residue_degree_coefficients(degree_coefficients);
    if (degree_coefficients.size() > 1) {
        degree_coefficients[1] -= 1;
    }
    if (have_linear_factor_count) {
        if (degree_coefficients.size() <= 1 ||
            linear_factor_count < 0 ||
            linear_factor_count > order.degree()) {
            return false;
        }
        // Higher residue degrees have norms outside both BF sums.
        degree_coefficients[1] += linear_factor_count;
    } else {
        if (!have_complete_residue_degrees) {
            return false;
        }
        for (slong residue_degree : residue_degrees) {
            if (residue_degree <= 0 ||
                residue_degree >=
                        static_cast<slong>(degree_coefficients.size())) {
                return false;
            }
            ++degree_coefficients[static_cast<std::size_t>(residue_degree)];
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.zeta_bf.term.weighted_prime_powers");
        if (!bf_accumulate_weighted_prime_powers(
                    sum, scratch.term, scratch.weighted_term,
                    scratch.summand, rational_prime, degree_coefficients, 1,
                    cutoff,
                    sqrt_cutoff_log_cutoff, precision)) {
            return false;
        }
    }

    if (rational_prime >= cutoff_div_9) {
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.zeta_bf.term.weighted_quotient_prime_powers");
        if (!bf_accumulate_weighted_prime_powers(
                    sum, scratch.term, scratch.weighted_term,
                    scratch.summand, rational_prime, degree_coefficients, -1,
                    cutoff_div_9,
                    sqrt_cutoff9_log_cutoff9, precision)) {
            return false;
        }
    }

    return true;
}

bool bf_term(flint::Arb& out,
             const Order& order,
             ulong cutoff,
             slong precision,
             const DiagnosticsContext* diagnostics,
             const FactorBase* residue_degree_base = nullptr,
             detail::ZetaBfResidueDegreeCache* residue_degree_cache = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.term");
    if (residue_degree_cache != nullptr) {
        residue_degree_cache->lookup_hint = 0;
    }
    cutoff = bf_adjust_cutoff(cutoff);
    const ulong cutoff_div_9 = cutoff / 9;
    if (cutoff_div_9 < 2) {
        return false;
    }

    flint::Arb sum;
    flint::Arb sqrt_cutoff;
    flint::Arb log_cutoff;
    flint::Arb sqrt_cutoff_log_cutoff;
    flint::Arb sqrt_cutoff9;
    flint::Arb log_cutoff9;
    flint::Arb sqrt_cutoff9_log_cutoff9;
    flint::Arb log_3cutoff;
    flint::Arb factor;
    BfResidueDegreeContext residue_context;
    std::vector<slong> residue_degrees;
    residue_degrees.reserve(static_cast<std::size_t>(order.degree()));
    std::vector<slong> residue_degree_coefficients(
            static_cast<std::size_t>(order.degree() + 1));
    BfPrimeTermScratch prime_term_scratch;
    bf_init_residue_degree_context(residue_context, order,
                                   residue_degree_base);
    std::vector<slong> degree_scratch(static_cast<std::size_t>(
            FLINT_MAX(order.degree(), WORD(1))));

    flint::arb_zero(sum);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.term.setup");
        flint::arb_set_ui(sqrt_cutoff, cutoff);
        flint::arb_sqrt(sqrt_cutoff, sqrt_cutoff, precision);
        flint::arb_log_ui(log_cutoff, cutoff, precision);
        flint::arb_mul(sqrt_cutoff_log_cutoff, sqrt_cutoff, log_cutoff,
                       precision);

        flint::arb_set_ui(sqrt_cutoff9, cutoff_div_9);
        flint::arb_sqrt(sqrt_cutoff9, sqrt_cutoff9, precision);
        flint::arb_log_ui(log_cutoff9, cutoff_div_9, precision);
        flint::arb_mul(sqrt_cutoff9_log_cutoff9, sqrt_cutoff9, log_cutoff9,
                       precision);
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.term.prime_loop");
        for (ulong p = 2; p < cutoff; p = n_nextprime(p, 1)) {
            if (!bf_add_prime_terms(sum, order, residue_context,
                                    residue_degrees, degree_scratch,
                                    residue_degree_coefficients,
                                    prime_term_scratch, p, cutoff, cutoff_div_9,
                                    sqrt_cutoff_log_cutoff,
                                    sqrt_cutoff9_log_cutoff9, precision,
                                    diagnostics, residue_degree_cache)) {
                return false;
            }
        }
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.term.final_scale");
        flint::arb_log_ui(log_3cutoff, 3 * cutoff, precision);
        flint::arb_mul(factor, sqrt_cutoff, log_3cutoff, precision);
        flint::arb_mul_2exp_si(factor, factor, 1);
        flint::arb_inv(factor, factor, precision);
        flint::arb_mul_ui(factor, factor, 3, precision);
        flint::arb_mul(out, factor, sum, precision);
    }
    return flint::arb_is_finite(out);
}

void bf_const_decimal(flint::Arb& out,
                      ulong numerator,
                      ulong denominator,
                      slong precision) noexcept {
    flint::arb_set_ui(out, numerator);
    flint::arb_div_ui(out, out, denominator, precision);
}

bool bf_error_bound(flint::Arb& out,
                    const Order& order,
                    ulong cutoff,
                    slong precision,
                    const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.error_bound");
    cutoff = bf_adjust_cutoff(cutoff);
    if (cutoff < kBfMinCutoff || cutoff / 9 <= 1) {
        return false;
    }

    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(abs_discriminant), 1) <= 0) {
        return false;
    }

    flint::Arb log_discriminant;
    flint::Arb sqrt_log_discriminant;
    flint::Arb sqrt_cutoff;
    flint::Arb log_3cutoff;
    flint::Arb log_cutoff9;
    flint::Arb c1;
    flint::Arb c2;
    flint::Arb c4;
    flint::Arb a1;
    flint::Arb a2;
    flint::Arb a3;
    flint::Arb a4;
    flint::Arb temp;

    flint::arb_log_fmpz(log_discriminant,
                        flint::FmpzConstRef(abs_discriminant), precision);
    flint::arb_sqrt(sqrt_log_discriminant, log_discriminant, precision);
    flint::arb_set_ui(sqrt_cutoff, cutoff);
    flint::arb_sqrt(sqrt_cutoff, sqrt_cutoff, precision);
    flint::arb_log_ui(log_3cutoff, 3 * cutoff, precision);
    flint::arb_log_ui(log_cutoff9, cutoff / 9, precision);

    bf_const_decimal(c1, 2324, 1000, precision);
    bf_const_decimal(c2, 388, 100, precision);
    bf_const_decimal(c4, 426, 100, precision);

    flint::arb_mul(a1, c1, log_discriminant, precision);
    flint::arb_mul(temp, sqrt_cutoff, log_3cutoff, precision);
    flint::arb_div(a1, a1, temp, precision);

    flint::arb_div(a2, c2, flint::ArbConstRef(log_cutoff9), precision);
    flint::arb_add_ui(a2, a2, 1, precision);

    flint::arb_inv(a3, sqrt_log_discriminant, precision);
    flint::arb_mul_2exp_si(a3, a3, 1);
    flint::arb_add_ui(a3, a3, 1, precision);
    flint::arb_mul(a3, a3, a3, precision);

    flint::arb_mul(a4, sqrt_cutoff, log_discriminant, precision);
    flint::arb_inv(a4, a4, precision);
    flint::arb_mul(a4, a4, c4, precision);
    flint::arb_mul_ui(a4, a4, static_cast<ulong>(order.degree() - 1),
                      precision);

    flint::arb_mul(temp, a2, a3, precision);
    flint::arb_add(temp, temp, a4, precision);
    flint::arb_mul(out, a1, temp, precision);
    return flint::arb_is_finite(out) && flint::arb_is_positive(out);
}

bool bf_select_cutoff(ulong& out,
                      const Order& order,
                      const flint::Arb& target,
                      ulong max_cutoff,
                      bool require_target,
                      slong precision,
                      const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.select_cutoff");
    ulong cutoff = bf_adjust_cutoff(kBfMinCutoff);
    max_cutoff = bf_adjust_cutoff(max_cutoff);
    if (max_cutoff < cutoff) {
        return false;
    }

    flint::Arb error_bound;
    for (;;) {
        if (!bf_error_bound(error_bound, order, cutoff, precision,
                            diagnostics)) {
            return false;
        }

        if (flint::arb_lt(error_bound, target)) {
            out = cutoff;
            return true;
        }

        if (cutoff >= max_cutoff) {
            if (!require_target) {
                out = cutoff;
                return true;
            }
            return false;
        }

        ulong next = 0;
        if (!bf_next_cutoff(next, cutoff, max_cutoff)) {
            return false;
        }
        cutoff = next;
    }
}

void bf_validation_targets(flint::Arb& validation_error,
                                 flint::Arb& tail_target,
                                 flint::Arb& numerical_radius_target,
                                 slong precision) noexcept {
    // reference src/NumFieldOrd/NfOrd/Clgp.jl validates with 0.6931 / 2,
    // and Zeta.jl:_residue_approx_bf reserves 2^-20 of that radius for
    // numerical evaluation of the BF term.
    flint::Fmpq error;
    flint::fmpq_set_si(error, WORD(6931), UWORD(20000));
    flint::arb_set_fmpq(validation_error, error, precision);
    flint::arb_one(numerical_radius_target);
    flint::arb_mul_2exp_si(numerical_radius_target,
                           numerical_radius_target,
                           -kValidationReserveBits);
    flint::arb_sub(tail_target, validation_error,
                   numerical_radius_target, precision);
}

bool bf_select_validation_cutoff(
        ulong& out,
        const Order& order,
        const flint::Arb& tail_target,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.zeta_validation.select_cutoff");
    ulong lower = bf_adjust_cutoff(kBfMinCutoff);
    max_cutoff = bf_adjust_cutoff(max_cutoff);
    if (max_cutoff < lower) {
        return false;
    }

    flint::Arb error_bound;
    if (!bf_error_bound(error_bound, order, lower, precision,
                        diagnostics)) {
        return false;
    }
    if (flint::arb_lt(error_bound, tail_target)) {
        out = lower;
        return true;
    }

    ulong upper = lower;
    for (;;) {
        if (!bf_next_cutoff(upper, lower, max_cutoff)) {
            return false;
        }
        if (!bf_error_bound(error_bound, order, upper, precision,
                            diagnostics)) {
            return false;
        }
        if (flint::arb_lt(error_bound, tail_target)) {
            break;
        }
        lower = upper;
    }

    // The doubling search above establishes an invalid/valid bracket.
    // Refine on the only cutoffs accepted by the BF term, multiples of 9,
    // retaining rigorous Arb comparisons at every candidate.
    ulong lower_index = lower / UWORD(9);
    ulong upper_index = upper / UWORD(9);
    while (upper_index - lower_index > 1) {
        const ulong middle_index =
                lower_index + (upper_index - lower_index) / 2;
        const ulong middle = UWORD(9) * middle_index;
        if (!bf_error_bound(error_bound, order, middle, precision,
                            diagnostics)) {
            return false;
        }
        if (flint::arb_lt(error_bound, tail_target)) {
            upper_index = middle_index;
        } else {
            lower_index = middle_index;
        }
    }

    out = UWORD(9) * upper_index;
    return true;
}

bool zeta_class_regulator_product_validation_impl(
        flint::ArbRef out,
        flint::ArbRef error_bound_out,
        ulong& cutoff_out,
        slong& work_precision_out,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        detail::ZetaBfResidueDegreeCache* residue_degree_cache = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.zeta_bf.zeta_validation.class_regulator_product");
    if (!supported_order(order, precision) || order.degree() <= 2 ||
        max_cutoff == 0) {
        return false;
    }

    slong work_precision = kValidationInitialWorkPrecision;
    ulong cutoff = 0;
    flint::Arb validation_error;
    flint::Arb tail_target;
    flint::Arb numerical_radius_target;
    bf_validation_targets(validation_error, tail_target,
                                numerical_radius_target, work_precision);
    if (!bf_select_validation_cutoff(
                cutoff, order, tail_target, max_cutoff, work_precision,
                diagnostics)) {
        return false;
    }

    flint::Arb term;
    flint::Arb log_residue;
    flint::Arb error_bound;
    for (;;) {
        bf_validation_targets(validation_error, tail_target,
                                    numerical_radius_target,
                                    work_precision);
        if (!bf_term(term, order, cutoff, work_precision, diagnostics,
                     residue_degree_base, residue_degree_cache) ||
            !bf_error_bound(error_bound, order, cutoff, work_precision,
                            diagnostics)) {
            return false;
        }

        flint::arb_set(flint::ArbRef(log_residue),
                       flint::ArbConstRef(term));
        flint::arb_add_error(log_residue, error_bound);
        if (bf_radius_lt(term, numerical_radius_target) &&
            bf_radius_lt(log_residue, validation_error)) {
            break;
        }

        if (work_precision >= kBfMaxWorkPrecision) {
            return false;
        }
        work_precision =
                FLINT_MIN(2 * work_precision, kBfMaxWorkPrecision);
    }

    flint::Arb residue;
    flint::Arb product;
    flint::Arb log_product;
    flint::arb_exp(residue, log_residue, precision);
    if (!flint::arb_is_finite(residue) ||
        !flint::arb_is_positive(residue) ||
        !class_regulator_product_from_residue(
                flint::ArbRef(product), order,
                flint::ArbConstRef(residue), precision, diagnostics)) {
        return false;
    }

    // The consumer proves an integer index from log(hR), so validate the
    // radius after the residue-to-hR conversion as well as before it.
    bf_validation_targets(validation_error, tail_target,
                                numerical_radius_target, precision);
    flint::arb_log(log_product, product, precision);
    if (!flint::arb_is_finite(log_product) ||
        !bf_radius_lt(log_product, validation_error)) {
        return false;
    }

    flint::arb_set(out, flint::ArbConstRef(product));
    flint::arb_set(error_bound_out,
                   flint::ArbConstRef(error_bound));
    cutoff_out = cutoff;
    work_precision_out = work_precision;
    return true;
}

bool bf_log_residue_cutoff(flint::ArbRef out,
                           const Order& order,
                           ulong max_cutoff,
                           bool require_target,
                           flint::ArbRef* audit_error_bound,
                           ulong* audit_cutoff,
                           slong* audit_work_precision,
                           slong precision,
                           const DiagnosticsContext* diagnostics,
                           const FactorBase* residue_degree_base = nullptr,
                           detail::ZetaBfResidueDegreeCache*
                                   residue_degree_cache = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.log_residue_cutoff");
    if (!supported_order(order, precision) || order.degree() <= 1) {
        return false;
    }

    flint::Arb value;
    flint::Arb error_bound;
    flint::Arb target;
    bf_target(target, precision);

    ulong cutoff = 0;
    slong work_precision = FLINT_MAX(precision + 64, WORD(192));
    if (!bf_select_cutoff(cutoff, order, target, max_cutoff, require_target,
                          work_precision, diagnostics)) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.refine_loop");
        for (;;) {
            if (!bf_term(value, order, cutoff, work_precision, diagnostics,
                         residue_degree_base, residue_degree_cache) ||
                !bf_error_bound(error_bound, order, cutoff, work_precision,
                                diagnostics)) {
                return false;
            }

            if (bf_radius_lt(value, target) ||
                bf_radius_lt(value, error_bound)) {
                break;
            }

            if (work_precision >= kBfMaxWorkPrecision) {
                return false;
            }
            work_precision =
                    FLINT_MIN(2 * work_precision, kBfMaxWorkPrecision);
        }
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.add_error");
        flint::arb_add_error(value, error_bound);
    }
    if (!flint::arb_is_finite(value)) {
        return false;
    }

    if (audit_error_bound != nullptr) {
        flint::arb_set(*audit_error_bound, flint::ArbConstRef(error_bound));
    }
    if (audit_cutoff != nullptr) {
        *audit_cutoff = cutoff;
    }
    if (audit_work_precision != nullptr) {
        *audit_work_precision = work_precision;
    }

    flint::arb_set(out, flint::ArbConstRef(value));
    return true;
}

bool bf_log_residue_default(flint::ArbRef out,
                            const Order& order,
                            slong precision,
                            const DiagnosticsContext* diagnostics,
                            const FactorBase* residue_degree_base = nullptr,
                            detail::ZetaBfResidueDegreeCache*
                                    residue_degree_cache = nullptr)
        noexcept {
    return bf_log_residue_cutoff(out, order, kBfDefaultMaxCutoff, false,
                                 nullptr, nullptr, nullptr, precision,
                                 diagnostics, residue_degree_base,
                                 residue_degree_cache);
}

}  // namespace

bool zeta_residue_impl(flint::ArbRef out,
                       const Order& order,
                       slong precision,
                       const DiagnosticsContext* diagnostics,
                       const FactorBase* residue_degree_base = nullptr,
                       detail::ZetaBfResidueDegreeCache*
                               residue_degree_cache = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.residue");
    if (!supported_order(order, precision)) {
        return false;
    }

    flint::Arb result;
    if (order.degree() == 1) {
        flint::arb_one(result);
    } else if (!quadratic_residue(result, order, precision)) {
        flint::Arb log_residue;
        if (!bf_log_residue_default(flint::ArbRef(log_residue),
                                    order, precision, diagnostics,
                                    residue_degree_base,
                                    residue_degree_cache)) {
            return false;
        }
        flint::arb_exp(result, log_residue, precision);
        if (!flint::arb_is_finite(result) ||
            !flint::arb_is_positive(result)) {
            return false;
        }
    }

    flint::arb_set(out, flint::ArbConstRef(result));
    return true;
}

bool zeta_residue(flint::ArbRef out,
                  const Order& order,
                  slong precision) noexcept {
    return zeta_residue_impl(out, order, precision, nullptr);
}

bool zeta_log_residue(flint::ArbRef out,
                      const Order& order,
                      slong precision) noexcept {
    if (!supported_order(order, precision)) {
        return false;
    }

    flint::Arb result;
    if (order.degree() == 1) {
        flint::arb_zero(result);
    } else {
        flint::Arb residue;
        if (quadratic_residue(residue, order, precision)) {
            if (!flint::arb_is_finite(residue) ||
                !flint::arb_is_positive(residue)) {
                return false;
            }
            flint::arb_log(result, residue, precision);
            if (!flint::arb_is_finite(result)) {
                return false;
            }
        } else if (!bf_log_residue_default(flint::ArbRef(result),
                                           order, precision, nullptr)) {
            return false;
        }
    }

    flint::arb_set(out, flint::ArbConstRef(result));
    return true;
}

bool zeta_log_residue_bf(flint::ArbRef out,
                         const Order& order,
                         ulong max_cutoff,
                         slong precision) noexcept {
    return bf_log_residue_cutoff(out, order, max_cutoff, true, nullptr,
                                 nullptr, nullptr, precision, nullptr);
}

bool zeta_log_residue_bf_audit(flint::ArbRef out,
                               flint::ArbRef error_bound,
                               ulong& cutoff,
                               slong& work_precision,
                               const Order& order,
                               ulong max_cutoff,
                               slong precision) noexcept {
    flint::ArbRef audit_error_bound(error_bound);
    return bf_log_residue_cutoff(out, order, max_cutoff, true,
                                 &audit_error_bound, &cutoff,
                                 &work_precision, precision, nullptr);
}

std::optional<ZetaBfAuditResult> zeta_log_residue_bf_audit(
        const Order& order,
        ulong max_cutoff,
        slong precision) noexcept {
    ZetaBfAuditResult result;
    if (!zeta_log_residue_bf_audit(flint::ArbRef(result.value),
                                   flint::ArbRef(result.error_bound),
                                   result.cutoff, result.work_precision,
                                   order, max_cutoff, precision)) {
        return std::nullopt;
    }
    return result;
}

bool zeta_residue_bf(flint::ArbRef out,
                     const Order& order,
                     ulong max_cutoff,
                     slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }

    flint::Arb log_residue;
    if (!zeta_log_residue_bf(flint::ArbRef(log_residue),
                             order, max_cutoff, precision)) {
        return false;
    }

    flint::arb_exp(log_residue, log_residue, precision);
    if (!flint::arb_is_finite(log_residue) ||
        !flint::arb_is_positive(log_residue)) {
        return false;
    }

    flint::arb_set(out, flint::ArbConstRef(log_residue));
    return true;
}

bool zeta_residue_bf_audit(flint::ArbRef out,
                           flint::ArbRef error_bound,
                           ulong& cutoff,
                           slong& work_precision,
                           const Order& order,
                           ulong max_cutoff,
                           slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }

    flint::Arb log_residue;
    flint::ArbRef audit_error_bound(error_bound);
    if (!bf_log_residue_cutoff(flint::ArbRef(log_residue), order,
                               max_cutoff, true, &audit_error_bound,
                               &cutoff, &work_precision, precision,
                               nullptr)) {
        return false;
    }

    flint::arb_exp(log_residue, log_residue, precision);
    if (!flint::arb_is_finite(log_residue) ||
        !flint::arb_is_positive(log_residue)) {
        return false;
    }

    flint::arb_set(out, flint::ArbConstRef(log_residue));
    return true;
}

std::optional<ZetaBfAuditResult> zeta_residue_bf_audit(
        const Order& order,
        ulong max_cutoff,
        slong precision) noexcept {
    ZetaBfAuditResult result;
    if (!zeta_residue_bf_audit(flint::ArbRef(result.value),
                               flint::ArbRef(result.error_bound),
                               result.cutoff, result.work_precision,
                               order, max_cutoff, precision)) {
        return std::nullopt;
    }
    return result;
}

bool zeta_class_regulator_product_impl(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        detail::ZetaBfResidueDegreeCache* residue_degree_cache = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.class_regulator_product");
    if (!supported_order(order, precision)) {
        return false;
    }

    flint::Arb result;
    if (order.degree() == 1) {
        flint::arb_one(result);
    } else {
        flint::Arb residue;
        if (!zeta_residue_impl(flint::ArbRef(residue), order, precision,
                               diagnostics, residue_degree_base,
                               residue_degree_cache) ||
            !class_regulator_product_from_residue(
                    flint::ArbRef(result), order, flint::ArbConstRef(residue),
                    precision, diagnostics)) {
            return false;
        }
    }

    flint::arb_set(out, flint::ArbConstRef(result));
    return true;
}

bool zeta_class_regulator_product(flint::ArbRef out,
                                  const Order& order,
                                  slong precision) noexcept {
    return zeta_class_regulator_product_impl(out, order, precision, nullptr);
}

bool zeta_class_regulator_product_bf(flint::ArbRef out,
                                     const Order& order,
                                     ulong max_cutoff,
                                     slong precision) noexcept {
    if (!supported_order(order, precision)) {
        return false;
    }

    if (order.degree() == 1) {
        flint::arb_set_si(out, 1);
        return true;
    }

    flint::Arb residue;
    return zeta_residue_bf(flint::ArbRef(residue), order, max_cutoff,
                           precision) &&
           class_regulator_product_from_residue(
                   out, order, flint::ArbConstRef(residue), precision,
                   nullptr);
}

bool zeta_class_regulator_product_bf_audit_impl(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        detail::ZetaBfResidueDegreeCache* residue_degree_cache = nullptr)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.zeta_bf.class_regulator_product_audit");
    if (!supported_order(order, precision)) {
        return false;
    }

    if (order.degree() == 1) {
        flint::arb_zero(error_bound);
        cutoff = 0;
        work_precision = precision;
        flint::arb_set_si(out, 1);
        return true;
    }

    flint::Arb log_residue;
    flint::ArbRef audit_error_bound(error_bound);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.audit_log_residue");
        if (!bf_log_residue_cutoff(flint::ArbRef(log_residue), order,
                                   max_cutoff, true, &audit_error_bound,
                                   &cutoff, &work_precision, precision,
                                   diagnostics, residue_degree_base,
                                   residue_degree_cache)) {
            return false;
        }
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.zeta_bf.exp_residue");
        flint::arb_exp(log_residue, log_residue, precision);
    }
    if (!flint::arb_is_finite(log_residue) ||
        !flint::arb_is_positive(log_residue)) {
        return false;
    }

    return class_regulator_product_from_residue(
            out, order, flint::ArbConstRef(log_residue), precision,
            diagnostics);
}

bool zeta_class_regulator_product_bf_audit(flint::ArbRef out,
                                           flint::ArbRef error_bound,
                                           ulong& cutoff,
                                           slong& work_precision,
                                           const Order& order,
                                           ulong max_cutoff,
                                           slong precision) noexcept {
    return zeta_class_regulator_product_bf_audit_impl(
            out, error_bound, cutoff, work_precision, order, max_cutoff,
            precision, nullptr);
}

std::optional<ZetaBfAuditResult> zeta_class_regulator_product_bf_audit(
        const Order& order,
        ulong max_cutoff,
        slong precision) noexcept {
    ZetaBfAuditResult result;
    if (!zeta_class_regulator_product_bf_audit(
                flint::ArbRef(result.value),
                flint::ArbRef(result.error_bound),
                result.cutoff, result.work_precision,
                order, max_cutoff, precision)) {
        return std::nullopt;
    }
    return result;
}

}  // namespace silex

namespace silex::detail {

bool grh_factor_base_bound_with_diagnostics(
        flint::FmpzRef out,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    return grh_factor_base_bound_impl(out, order, diagnostics);
}

bool zeta_class_regulator_product_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept {
    return zeta_class_regulator_product_impl(out, order, precision,
                                            diagnostics,
                                            residue_degree_base,
                                            residue_degree_cache);
}

bool zeta_class_regulator_product_bf_audit_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept {
    return zeta_class_regulator_product_bf_audit_impl(
            out, error_bound, cutoff, work_precision, order, max_cutoff,
            precision, diagnostics, residue_degree_base,
            residue_degree_cache);
}

bool zeta_class_regulator_product_validation_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base,
        ZetaBfResidueDegreeCache* residue_degree_cache) noexcept {
    return zeta_class_regulator_product_validation_impl(
            out, error_bound, cutoff, work_precision, order, max_cutoff,
            precision, diagnostics, residue_degree_base,
            residue_degree_cache);
}

bool zeta_bf_audit_cutoff_available(const Order& order,
                                    ulong max_cutoff,
                                    slong precision) noexcept {
    if (!supported_order(order, precision) || max_cutoff == 0) {
        return false;
    }
    if (order.degree() == 1) {
        return true;
    }

    flint::Arb target;
    bf_target(target, precision);
    ulong cutoff = 0;
    const slong work_precision = FLINT_MAX(precision + 64, WORD(192));
    return bf_select_cutoff(cutoff, order, target, max_cutoff, true,
                            work_precision, nullptr);
}

bool class_regulator_product_estimate_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        ZetaBfResidueDegreeCache* residue_degree_cache,
        flint::Fmpz* torsion_order,
        Element* torsion_generator) noexcept {
    return class_regulator_product_estimate_impl(out, order, precision,
                                                     diagnostics,
                                                     residue_degree_cache,
                                                     torsion_order,
                                                     torsion_generator);
}

}  // namespace silex::detail
