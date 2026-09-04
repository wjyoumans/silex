#include "relation_factor_base_plan_internal.hpp"

#include <silex/factor_base.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/signature.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace silex::detail {

ZetaBfResidueDegreeCache* relation_factor_base_plan_residue_degrees(
        RelationFactorBasePlan& plan,
        const Order& order) noexcept {
    return plan.valid && plan.order.has_same_data(order)
            ? &plan.residue_degrees
            : nullptr;
}

const ZetaBfResidueDegreeCache* relation_factor_base_plan_residue_degrees(
        const RelationFactorBasePlan& plan,
        const Order& order) noexcept {
    return plan.valid && plan.order.has_same_data(order)
            ? &plan.residue_degrees
            : nullptr;
}

namespace relation_search {
namespace {

struct ResidueDegreeGroup {
    slong residue_degree = 0;
    slong count = 0;
};

bool defining_polynomial_is_monic_integral(
        const fmpq_poly_t polynomial) noexcept {
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

void reduce_defining_polynomial(
        flint::FmpzModPoly& out,
        const fmpq_poly_t polynomial,
        const flint::FmpzModCtx& ctx) noexcept {
    flint::Fmpz coeff;
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpz(coeff.raw(), polynomial, i);
        fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, coeff.raw(), ctx.raw());
    }
}

bool push_residue_degree_group(std::vector<ResidueDegreeGroup>& groups,
                               slong residue_degree) noexcept {
    if (residue_degree <= 0) {
        return false;
    }
    auto found = std::find_if(
            groups.begin(), groups.end(),
            [residue_degree](const ResidueDegreeGroup& group) {
                return group.residue_degree == residue_degree;
            });
    if (found == groups.end()) {
        groups.push_back(ResidueDegreeGroup{residue_degree, 1});
    } else {
        ++found->count;
    }
    return true;
}

void sort_residue_degree_groups(
        std::vector<ResidueDegreeGroup>& groups) noexcept {
    std::sort(groups.begin(), groups.end(),
              [](const ResidueDegreeGroup& left,
                 const ResidueDegreeGroup& right) {
                  return left.residue_degree < right.residue_degree;
              });
}

bool prime_degree_groups_from_nonindex_polynomial(
        std::vector<ResidueDegreeGroup>& groups,
        const Order& order,
        flint::FmpzConstRef rational_prime,
        flint::FmpzConstRef equation_order_index) noexcept {
    groups.clear();
    const NumberField* field = order.parent();
    const nf_struct* raw_field =
            field == nullptr ? nullptr : field->raw_flint_field();
    if (raw_field == nullptr || !order.is_maximal() ||
        !defining_polynomial_is_monic_integral(raw_field->pol) ||
        fmpz_divisible(equation_order_index.raw(),
                       rational_prime.raw()) != 0) {
        return false;
    }

    // This fast path needs only residue degrees and multiplicities, so it
    // factors the defining polynomial modulo the non-index prime directly.
    flint::FmpzModCtx ctx(rational_prime.raw());
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPolyFactor factorization(ctx);
    reduce_defining_polynomial(reduced, raw_field->pol, ctx);
    fmpz_mod_poly_factor(factorization.raw(), reduced.raw(), ctx.raw());
    const slong num_factors = factorization.raw()->num;
    if (num_factors <= 0) {
        return false;
    }

    groups.reserve(static_cast<std::size_t>(num_factors));
    for (slong i = 0; i < num_factors; ++i) {
        const slong residue_degree =
                fmpz_mod_poly_degree(factorization.raw()->poly + i,
                                     ctx.raw());
        if (!push_residue_degree_group(groups, residue_degree)) {
            return false;
        }
    }
    sort_residue_degree_groups(groups);
    return true;
}

bool prime_degree_groups(
        std::vector<ResidueDegreeGroup>& groups,
        const Order& order,
        flint::FmpzConstRef rational_prime,
        flint::FmpzConstRef equation_order_index) noexcept {
    if (prime_degree_groups_from_nonindex_polynomial(
                groups, order, rational_prime, equation_order_index)) {
        return true;
    }

    groups.clear();
    PrimeIdealList decomposed;
    if (!decompose_prime(decomposed, order, rational_prime)) {
        return false;
    }

    for (slong i = 0; i < decomposed.size(); ++i) {
        const PrimeIdeal* prime = decomposed.at(i);
        if (prime == nullptr) {
            return false;
        }
        const slong residue_degree = prime->residue_degree();
        if (residue_degree <= 0) {
            return false;
        }
        if (!push_residue_degree_group(groups, residue_degree)) {
            return false;
        }
    }

    sort_residue_degree_groups(groups);
    return true;
}

struct PrimeDegreeCacheEntry {
    slong prime = 0;
    std::vector<ResidueDegreeGroup> groups;
};

struct PrimeDegreeCache {
    std::vector<PrimeDegreeCacheEntry> entries;
    flint::Fmpz last_cached_prime;
    flint::Fmpz equation_order_index;
};

bool export_prime_degree_cache(
        ZetaBfResidueDegreeCache& out,
        const PrimeDegreeCache& source) noexcept {
    ZetaBfResidueDegreeCache next;
    next.entries.reserve(source.entries.size());
    for (const PrimeDegreeCacheEntry& source_entry : source.entries) {
        if (source_entry.prime < 2 || source_entry.groups.empty()) {
            return false;
        }
        const std::size_t offset = next.residue_degrees.size();
        for (const ResidueDegreeGroup& group : source_entry.groups) {
            if (group.residue_degree <= 0 || group.count <= 0 ||
                static_cast<std::size_t>(group.count) >
                        next.residue_degrees.max_size() -
                                next.residue_degrees.size()) {
                return false;
            }
            for (slong i = 0; i < group.count; ++i) {
                next.residue_degrees.push_back(group.residue_degree);
            }
        }
        const std::size_t length = next.residue_degrees.size() - offset;
        if (length == 0) {
            return false;
        }
        next.entries.push_back(ZetaBfResidueDegreeCacheEntry{
                static_cast<ulong>(source_entry.prime), offset, length});
    }
    next.lookup_hint = 0;
    out = std::move(next);
    return !out.entries.empty();
}

bool append_next_prime_degree_cache(
        PrimeDegreeCache& cache,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    flint::Fmpz next_prime;
    flint::fmpz_nextprime(flint::FmpzRef(next_prime),
                          flint::FmpzConstRef(cache.last_cached_prime));
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(next_prime))) {
        return false;
    }

    PrimeDegreeCacheEntry entry;
    entry.prime = flint::fmpz_get_si(flint::FmpzConstRef(next_prime));
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_factor_base_plan.prime_decomposition");
        if (!prime_degree_groups(
                    entry.groups, order, flint::FmpzConstRef(next_prime),
                    flint::FmpzConstRef(cache.equation_order_index))) {
            return false;
        }
    }
    cache.entries.push_back(std::move(entry));
    flint::fmpz_set(flint::FmpzRef(cache.last_cached_prime),
                    flint::FmpzConstRef(next_prime));
    return true;
}

bool extend_prime_degree_cache(PrimeDegreeCache& cache,
                               const Order& order,
                               slong bound,
                               const DiagnosticsContext* diagnostics)
        noexcept {
    if (bound <= 1) {
        return true;
    }

    while (cache.entries.empty() || cache.entries.back().prime <= bound) {
        if (!append_next_prime_degree_cache(cache, order, diagnostics)) {
            return false;
        }
    }
    return true;
}

bool append_next_quadratic_prime_degree_cache(
        PrimeDegreeCache& cache,
        flint::FmpzConstRef discriminant) noexcept {
    flint::Fmpz next_prime;
    flint::fmpz_nextprime(flint::FmpzRef(next_prime),
                          flint::FmpzConstRef(cache.last_cached_prime));
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(next_prime))) {
        return false;
    }

    PrimeDegreeCacheEntry entry;
    entry.prime = flint::fmpz_get_si(flint::FmpzConstRef(next_prime));
    const int splitting = flint::fmpz_kronecker(
            discriminant, flint::FmpzConstRef(next_prime));
    if (splitting < 0) {
        entry.groups.push_back(ResidueDegreeGroup{2, 1});
    } else {
        entry.groups.push_back(
                ResidueDegreeGroup{1, splitting > 0 ? 2 : 1});
    }
    cache.entries.emplace_back(std::move(entry));
    flint::fmpz_set(flint::FmpzRef(cache.last_cached_prime),
                    flint::FmpzConstRef(next_prime));
    return true;
}

bool extend_quadratic_prime_degree_cache(
        PrimeDegreeCache& cache,
        flint::FmpzConstRef discriminant,
        slong bound) noexcept {
    if (bound <= 1) {
        return true;
    }
    while (cache.entries.empty() || cache.entries.back().prime <= bound) {
        if (!append_next_quadratic_prime_degree_cache(cache, discriminant)) {
            return false;
        }
    }
    return true;
}

// Preserve grouped residue-degree order and the in-place update of the
// retained lower-bound values. Sorted prime-ideal norms are not equivalent.
bool insert_nth_ideal_group(std::vector<flint::Fmpz>& values,
                            flint::FmpzConstRef norm,
                            slong count) noexcept {
    if (count <= 0) {
        return false;
    }
    const slong n = static_cast<slong>(values.size());
    slong k = 0;
    while (k < n &&
           flint::fmpz_cmp(
                   flint::FmpzConstRef(values[static_cast<std::size_t>(k)]),
                   norm) <= 0) {
        ++k;
    }
    if (k >= n) {
        return true;
    }

    for (slong l = k + count; l < n; ++l) {
        flint::fmpz_set(
                flint::FmpzRef(values[static_cast<std::size_t>(l)]),
                flint::FmpzConstRef(
                        values[static_cast<std::size_t>(l - count)]));
    }
    slong l = 0;
    for (; l < count && k + l < n; ++l) {
        flint::fmpz_set(
                flint::FmpzRef(values[static_cast<std::size_t>(k + l)]),
                norm);
    }
    while (l <= k) {
        flint::fmpz_set(
                flint::FmpzRef(values[static_cast<std::size_t>(l - 1)]),
                norm);
        ++l;
    }
    return true;
}

bool nth_ideal_lower_bound(flint::Fmpz& out,
                           const Order& order,
                           PrimeDegreeCache& prime_cache,
                           const DiagnosticsContext* diagnostics) noexcept {
    const slong degree = order.degree();
    if (degree <= 0) {
        return false;
    }
    if (degree == 1) {
        flint::fmpz_one(flint::FmpzRef(out));
        return true;
    }

    std::vector<flint::Fmpz> values(static_cast<std::size_t>(degree));
    for (flint::Fmpz& value : values) {
        flint::fmpz_set_si(flint::FmpzRef(value),
                           std::numeric_limits<slong>::max());
    }
    std::size_t prime_index = 0;
    while (true) {
        if (prime_index >= prime_cache.entries.size() &&
            !append_next_prime_degree_cache(prime_cache, order,
                                            diagnostics)) {
            return false;
        }
        const PrimeDegreeCacheEntry& entry =
                prime_cache.entries[prime_index];
        ++prime_index;

        if (!entry.groups.empty() &&
            entry.groups.front().residue_degree != degree) {
            for (std::size_t group_index = entry.groups.size();
                 group_index > 0;
                 --group_index) {
                const ResidueDegreeGroup& group =
                        entry.groups[group_index - 1];
                flint::Fmpz norm;
                flint::fmpz_set_ui(flint::FmpzRef(norm),
                                    static_cast<ulong>(entry.prime));
                flint::fmpz_pow_ui(flint::FmpzRef(norm),
                                    flint::FmpzConstRef(norm),
                                    static_cast<ulong>(
                                            group.residue_degree));
                if (!flint::fmpz_abs_fits_ui(flint::FmpzConstRef(norm))) {
                    continue;
                }
                if (!insert_nth_ideal_group(
                            values, flint::FmpzConstRef(norm),
                            group.count)) {
                    return false;
                }
            }
        }

        if (flint::fmpz_cmp_ui(
                    flint::FmpzConstRef(
                            values[static_cast<std::size_t>(degree - 1)]),
                    static_cast<ulong>(entry.prime)) < 0) {
            break;
        }
    }

    flint::fmpz_set(
            flint::FmpzRef(out),
            flint::FmpzConstRef(
                    values[static_cast<std::size_t>(degree - 1)]));
    return true;
}

bool analytic_bound_satisfied(PrimeDegreeCache& cache,
                              const Order& order,
                              slong degree,
                              slong real_embeddings,
                              double log_discriminant,
                              slong bound,
                              const DiagnosticsContext* diagnostics) noexcept {
    if (bound <= 1) {
        return false;
    }
    if (!extend_prime_degree_cache(cache, order, bound, diagnostics)) {
        return false;
    }

    constexpr double kPi = 3.141592653589793238462643383279502884;
    constexpr double kPiSquaredOverTwo = kPi * kPi / 2.0;
    constexpr double kC2 = 3.663862376709;
    constexpr double kC3 = 3.801387092431;
    const double c_n =
            static_cast<double>(real_embeddings) * kC2 +
            static_cast<double>(degree) * kPiSquaredOverTwo;
    const double c_d =
            log_discriminant - static_cast<double>(degree) * kC3 -
            static_cast<double>(real_embeddings) * kPi / 2.0;

    const double log_c = std::log(static_cast<double>(bound));
    double sa = 0.0;
    double sb = 0.0;

    for (const PrimeDegreeCacheEntry& entry : cache.entries) {
        const slong prime_si = entry.prime;
        if (prime_si > bound) {
            break;
        }

        const double log_p = std::log(static_cast<double>(prime_si));
        const double log_c_over_log_p = log_c / log_p;
        for (const ResidueDegreeGroup& group : entry.groups) {
            const slong residue_degree = group.residue_degree;
            if (static_cast<double>(residue_degree) > log_c_over_log_p) {
                break;
            }
            const double log_norm =
                    static_cast<double>(residue_degree) * log_p;
            const double norm =
                    std::pow(static_cast<double>(prime_si),
                             static_cast<double>(residue_degree));
            if (!std::isfinite(norm) || norm <= 0.0) {
                return false;
            }
            const double q = 1.0 / std::sqrt(norm);
            double a = log_norm * q;
            double b = log_norm * a;
            const slong m = static_cast<slong>(
                    log_c_over_log_p /
                    static_cast<double>(residue_degree));
            if (m > 1) {
                const double inv1_q = 1.0 / (1.0 - q);
                const double q_to_m =
                        std::pow(q, static_cast<double>(m));
                a *= (1.0 - q_to_m) * inv1_q;
                b *= (1.0 -
                      q_to_m *
                              (static_cast<double>(m + 1) -
                               static_cast<double>(m) * q)) *
                     inv1_q * inv1_q;
            }
            sa += static_cast<double>(group.count) * a;
            sb += static_cast<double>(group.count) * b;
        }
    }

    return c_d + (c_n + 2.0 * sb) / log_c - 2.0 * sa < -1e-8;
}

bool compute_analytic_bound(flint::Fmpz& out,
                            const Order& order,
                            slong degree,
                            slong real_embeddings,
                            double log_discriminant,
                            double log_discriminant_squared,
                            PrimeDegreeCache& prime_cache,
                            const DiagnosticsContext* diagnostics) noexcept {
    const double raw_max = 4.0 * log_discriminant_squared;
    if (!std::isfinite(raw_max) || raw_max < 1.0 ||
        raw_max >
                static_cast<double>(std::numeric_limits<slong>::max() / 2)) {
        return false;
    }
    const slong max_bound = static_cast<slong>(raw_max);
    slong low = 1;
    slong high = 1;
    while (!analytic_bound_satisfied(prime_cache, order, degree,
                                     real_embeddings, log_discriminant,
                                     high, diagnostics)) {
        low = high;
        if (high > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        high *= 2;
    }
    while (high - low > 1) {
        const slong test = low + (high - low) / 2;
        if (analytic_bound_satisfied(prime_cache, order, degree,
                                     real_embeddings, log_discriminant,
                                     test, diagnostics)) {
            high = test;
        } else {
            low = test;
        }
    }

    slong bound =
            (high == 2 &&
             analytic_bound_satisfied(prime_cache, order, degree,
                                      real_embeddings, log_discriminant,
                                      1, diagnostics))
                    ? 1
                    : high;
    if (bound > max_bound) {
        bound = max_bound;
    }
    if (bound < 1) {
        bound = 1;
    }
    flint::fmpz_set_si(flint::FmpzRef(out), bound);
    return true;
}

bool compute_quadratic_analytic_bound(
        flint::Fmpz& out,
        const Order& order,
        flint::FmpzConstRef discriminant,
        double log_discriminant,
        double log_discriminant_squared,
        PrimeDegreeCache& prime_cache) noexcept {
    const double raw_max = 4.0 * log_discriminant_squared;
    if (!std::isfinite(raw_max) || raw_max < 1.0 ||
        raw_max >
                static_cast<double>(std::numeric_limits<slong>::max() / 2)) {
        return false;
    }
    const slong max_bound = static_cast<slong>(raw_max);
    const auto bound_satisfied = [&](slong bound) noexcept {
        return extend_quadratic_prime_degree_cache(
                       prime_cache, discriminant, bound) &&
               analytic_bound_satisfied(
                       prime_cache, order, 2, 0, log_discriminant, bound,
                       nullptr);
    };

    slong low = 1;
    slong high = 1;
    while (!bound_satisfied(high)) {
        low = high;
        if (high > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        high *= 2;
    }
    while (high - low > 1) {
        const slong test = low + (high - low) / 2;
        if (bound_satisfied(test)) {
            high = test;
        } else {
            low = test;
        }
    }

    slong bound = high == 2 && bound_satisfied(1) ? 1 : high;
    if (bound > max_bound) {
        bound = max_bound;
    }
    if (bound < 1) {
        bound = 1;
    }
    flint::fmpz_set_si(flint::FmpzRef(out), bound);
    return true;
}

bool quadratic_prime_is_bad(flint::FmpzConstRef discriminant,
                            flint::FmpzConstRef rational_prime) noexcept {
    if (flint::fmpz_equal_si(rational_prime, 2)) {
        ulong residue = ::fmpz_fdiv_ui(discriminant.raw(), UWORD(16)) >> 1U;
        if (residue != 0 && flint::fmpz_sgn(discriminant) < 0) {
            residue = 8U - residue;
        }
        return residue < 4U;
    }

    flint::Fmpz square;
    flint::fmpz_mul(flint::FmpzRef(square), rational_prime, rational_prime);
    return flint::fmpz_divisible(discriminant,
                                flint::FmpzConstRef(square));
}

bool quadratic_nth_suitable_ideal_bound(
        flint::Fmpz& out,
        flint::FmpzConstRef discriminant,
        slong count) noexcept {
    if (count <= 0) {
        return false;
    }
    flint::Fmpz rational_prime;
    flint::fmpz_one(flint::FmpzRef(rational_prime));
    while (count > 0) {
        flint::fmpz_nextprime(flint::FmpzRef(rational_prime),
                              flint::FmpzConstRef(rational_prime));
        if (!quadratic_prime_is_bad(
                    discriminant, flint::FmpzConstRef(rational_prime)) &&
            flint::fmpz_kronecker(
                    discriminant,
                    flint::FmpzConstRef(rational_prime)) >= 0) {
            --count;
        }
    }
    flint::fmpz_set(flint::FmpzRef(out),
                    flint::FmpzConstRef(rational_prime));
    return true;
}

}  // namespace

bool build_relation_factor_base_plan(
        RelationFactorBasePlan& out,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_factor_base_plan");
    if (order.parent() == nullptr) {
        return false;
    }
    const slong degree = order.degree();
    if (degree <= 0) {
        return false;
    }

    RelationFactorBasePlan next;
    if (degree == 1) {
        flint::fmpz_one(flint::FmpzRef(next.working_bound));
        next.order = order;
        next.valid = true;
        out = std::move(next);
        return true;
    }

    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    const double discriminant_d =
            flint::fmpz_get_d(flint::FmpzConstRef(abs_discriminant));
    if (!std::isfinite(discriminant_d) || discriminant_d <= 1.0) {
        return false;
    }

    Signature signature_value;
    if (!signature(signature_value, *order.parent()) ||
        signature_value.degree() != degree) {
        return false;
    }
    const double log_discriminant = std::log(discriminant_d);
    const double log_discriminant_squared =
            log_discriminant * log_discriminant;

    PrimeDegreeCache prime_cache;
    if (order.is_maximal()) {
        Order equation = Order::equation_order(*order.parent());
        if (equation.is_defined()) {
            (void) order_index(
                    flint::FmpzRef(prime_cache.equation_order_index),
                    equation, order);
        }
    }
    flint::Fmpz analytic_bound;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_factor_base_plan.analytic_bound");
        if (!compute_analytic_bound(
                    analytic_bound, order, degree, signature_value.r1(),
                    log_discriminant, log_discriminant_squared, prime_cache,
                    diagnostics)) {
            return false;
        }
    }

    flint::Fmpz nth_ideal_bound;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_factor_base_plan.nth_ideal");
        if (!nth_ideal_lower_bound(nth_ideal_bound, order, prime_cache,
                                   diagnostics)) {
            return false;
        }
    }
    if (flint::fmpz_cmp(flint::FmpzConstRef(analytic_bound),
                       flint::FmpzConstRef(nth_ideal_bound)) >= 0) {
        flint::fmpz_set(flint::FmpzRef(next.working_bound),
                        flint::FmpzConstRef(analytic_bound));
    } else {
        flint::fmpz_set(flint::FmpzRef(next.working_bound),
                        flint::FmpzConstRef(nth_ideal_bound));
    }
    if (!export_prime_degree_cache(next.residue_degrees, prime_cache)) {
        return false;
    }
    next.order = order;
    next.valid = true;
    out = std::move(next);
    return true;
}

bool build_maximal_imaginary_quadratic_factor_base_plan(
        RelationFactorBasePlan& out,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.maximal_imaginary_quadratic_factor_base_plan");
    Signature signature_value;
    flint::Fmpz discriminant;
    flint::Fmpz conductor;
    if (order.parent() == nullptr || order.degree() != 2 ||
        !order.is_maximal() ||
        order.parent()->backend_kind() != NumberFieldBackendKind::quadratic ||
        !order.quadratic_conductor(flint::FmpzRef(conductor)) ||
        !flint::fmpz_is_one(flint::FmpzConstRef(conductor)) ||
        !signature(signature_value, *order.parent()) ||
        signature_value.r1() != 0 || signature_value.r2() != 1 ||
        !order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) >= 0 ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(discriminant))) {
        return false;
    }

    flint::Fmpz absolute_discriminant;
    flint::fmpz_abs(flint::FmpzRef(absolute_discriminant),
                    flint::FmpzConstRef(discriminant));
    const double discriminant_d =
            flint::fmpz_get_d(flint::FmpzConstRef(absolute_discriminant));
    if (!std::isfinite(discriminant_d) || discriminant_d <= 1.0) {
        return false;
    }
    const double log_discriminant = std::log(discriminant_d);
    const double log_discriminant_squared =
            log_discriminant * log_discriminant;

    PrimeDegreeCache prime_cache;
    flint::Fmpz analytic_bound;
    flint::Fmpz nth_ideal_bound;
    if (!compute_quadratic_analytic_bound(
                analytic_bound, order,
                flint::FmpzConstRef(discriminant), log_discriminant,
                log_discriminant_squared, prime_cache) ||
        !quadratic_nth_suitable_ideal_bound(
                nth_ideal_bound, flint::FmpzConstRef(discriminant), 2)) {
        return false;
    }

    RelationFactorBasePlan next;
    if (flint::fmpz_cmp(flint::FmpzConstRef(analytic_bound),
                        flint::FmpzConstRef(nth_ideal_bound)) >= 0) {
        flint::fmpz_set(flint::FmpzRef(next.working_bound),
                        flint::FmpzConstRef(analytic_bound));
    } else {
        flint::fmpz_set(flint::FmpzRef(next.working_bound),
                        flint::FmpzConstRef(nth_ideal_bound));
    }
    next.order = order;
    next.valid = true;
    if (!specialize_maximal_imaginary_quadratic_factor_base_plan(
                next, order, diagnostics) ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(next.working_bound)) ||
        !extend_quadratic_prime_degree_cache(
                prime_cache, flint::FmpzConstRef(discriminant),
                flint::fmpz_get_si(
                        flint::FmpzConstRef(next.working_bound))) ||
        !export_prime_degree_cache(next.residue_degrees, prime_cache)) {
        return false;
    }
    out = std::move(next);
    return true;
}

bool relation_factor_base_restart_limit(flint::Fmpz& out,
                                        const Order& order) noexcept {
    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant))) {
        return false;
    }
    flint::fmpz_abs(flint::FmpzRef(abs_discriminant),
                    flint::FmpzConstRef(discriminant));
    const double discriminant_d =
            flint::fmpz_get_d(flint::FmpzConstRef(abs_discriminant));
    if (!std::isfinite(discriminant_d) || discriminant_d <= 1.0) {
        return false;
    }

    const double log_discriminant = std::log(discriminant_d);
    const double raw_bound = 4.0 * log_discriminant * log_discriminant;
    if (!std::isfinite(raw_bound) || raw_bound < 2.0 ||
        raw_bound >
                static_cast<double>(std::numeric_limits<slong>::max())) {
        return false;
    }

    flint::fmpz_set_si(flint::FmpzRef(out),
                       static_cast<slong>(raw_bound));
    return true;
}

bool next_relation_factor_base_bound_to_limit(
        flint::Fmpz& out,
        flint::FmpzConstRef current,
        flint::FmpzConstRef limit) noexcept {
    if (!flint::fmpz_fits_si(current) ||
        !flint::fmpz_fits_si(limit)) {
        return false;
    }
    const slong current_si = flint::fmpz_get_si(current);
    const slong limit_si = flint::fmpz_get_si(limit);
    if (current_si < 1 || limit_si < 1 || current_si >= limit_si) {
        return false;
    }

    slong next = current_si;
    if (static_cast<double>(current_si) <=
        static_cast<double>(limit_si) / 13.333) {
        if (current_si > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        next = 2 * current_si;
    } else {
        const slong increment = limit_si / 20 > 1 ? limit_si / 20 : 1;
        if (current_si > std::numeric_limits<slong>::max() - increment) {
            return false;
        }
        next = current_si + increment;
    }
    if (next > limit_si) {
        next = limit_si;
    }
    if (next <= current_si) {
        return false;
    }

    flint::fmpz_set_si(flint::FmpzRef(out), next);
    return true;
}

bool next_relation_factor_base_bound(
        flint::Fmpz& out,
        flint::FmpzConstRef current,
        flint::FmpzConstRef limit) noexcept {
    if (!flint::fmpz_fits_si(current) ||
        !flint::fmpz_fits_si(limit)) {
        return false;
    }
    if (flint::fmpz_get_si(current) > flint::fmpz_get_si(limit) / 2) {
        return false;
    }
    return next_relation_factor_base_bound_to_limit(out, current, limit);
}

bool specialize_maximal_imaginary_quadratic_factor_base_plan(
        RelationFactorBasePlan& plan,
        const Order& order,
        const DiagnosticsContext* diagnostics) noexcept {
    Signature signature_value;
    flint::Fmpz discriminant;
    flint::Fmpz conductor;
    if (!plan.valid || !plan.order.has_same_data(order) ||
        order.parent() == nullptr || order.degree() != 2 ||
        !order.is_maximal() ||
        order.parent()->backend_kind() != NumberFieldBackendKind::quadratic ||
        !order.quadratic_conductor(flint::FmpzRef(conductor)) ||
        !flint::fmpz_is_one(flint::FmpzConstRef(conductor)) ||
        !signature(signature_value, *order.parent()) ||
        signature_value.r1() != 0 || signature_value.r2() != 1 ||
        !order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) >= 0 ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(discriminant))) {
        return false;
    }

    flint::Fmpz absolute_discriminant;
    flint::fmpz_abs(flint::FmpzRef(absolute_discriminant),
                    flint::FmpzConstRef(discriminant));
    const slong minimum_split_primes =
            ::fmpz_bits(absolute_discriminant.raw()) > 16 ? 3 : 2;

    flint::Fmpz restart_limit;
    if (!relation_factor_base_restart_limit(restart_limit, order)) {
        return false;
    }
    const bool exceptional_discriminant_minus_three =
            flint::fmpz_equal_si(
                    flint::FmpzConstRef(discriminant), -3);
    if (!exceptional_discriminant_minus_three) {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_factor_base_plan."
                "quadratic_split_restart");
        for (;;) {
            slong split_primes = 0;
            flint::Fmpz rational_prime;
            flint::fmpz_one(flint::FmpzRef(rational_prime));
            for (;;) {
                flint::fmpz_nextprime(
                        flint::FmpzRef(rational_prime),
                        flint::FmpzConstRef(rational_prime));
                if (flint::fmpz_cmp(
                            flint::FmpzConstRef(rational_prime),
                            flint::FmpzConstRef(plan.working_bound)) > 0) {
                    break;
                }
                if (flint::fmpz_kronecker(
                            flint::FmpzConstRef(discriminant),
                            flint::FmpzConstRef(rational_prime)) > 0) {
                    ++split_primes;
                }
            }
            if (split_primes >= minimum_split_primes) {
                break;
            }

            flint::Fmpz next_bound;
            if (!next_relation_factor_base_bound_to_limit(
                        next_bound,
                        flint::FmpzConstRef(plan.working_bound),
                        flint::FmpzConstRef(restart_limit))) {
                return false;
            }
            plan.working_bound.swap(next_bound);
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_factor_base_plan."
                "quadratic_generation_bound");
        flint::Fmpz generation_bound;
        if (!factor_base_class_group_bound(
                    flint::FmpzRef(generation_bound), order)) {
            return false;
        }
        if (flint::fmpz_cmp(
                    flint::FmpzConstRef(plan.working_bound),
                    flint::FmpzConstRef(generation_bound)) < 0) {
            plan.working_bound.swap(generation_bound);
        }
    }
    return true;
}

}  // namespace relation_search
}  // namespace silex::detail
