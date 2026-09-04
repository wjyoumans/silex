#include <silex/class_group.hpp>

#include <silex/archimedean.hpp>
#include <silex/embedding.hpp>
#include <silex/detail/class_relation_module_context.hpp>
#include <silex/fmpz_smat.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/lat.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>
#include <silex/residue_field.hpp>
#include <silex/signature.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include "class_group_internal.hpp"
#include "class_group_certification_internal.hpp"
#include "factor_base_proof_targets_internal.hpp"
#include "relation_saturation_internal.hpp"
#include "ideal_lattice_reduction_internal.hpp"
#include "../ideal_factorization/ideal_factorization_internal.hpp"
#include "../order_unit/class_unit_transaction_internal.hpp"
#include "../order_unit/relation_unit_internal.hpp"
#include "../order_unit/order_unit_internal.hpp"
#include "../order/order_internal.hpp"
#include "../residue_field/residue_field_internal.hpp"
#include "../zeta/zeta_internal.hpp"

#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>
#include <flint/qfb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace silex {

bool relation_saturation_torsion_needed(bool& out,
                                        const OrderUnitGroup& units,
                                        flint::FmpzConstRef ell) noexcept;

namespace detail {

bool imaginary_quadratic_class_number(flint::FmpzRef out,
                                      flint::FmpzConstRef discriminant)
        noexcept {
    if (flint::fmpz_sgn(discriminant) >= 0 ||
        !flint::fmpz_fits_si(discriminant)) {
        return false;
    }

    struct ReducedForms {
        ~ReducedForms() noexcept {
            if (forms != nullptr) {
                qfb_array_clear(&forms, count);
            }
        }

        qfb* forms = nullptr;
        slong count = 0;
    } reduced;

    reduced.count = qfb_reduced_forms(
            &reduced.forms, flint::fmpz_get_si(discriminant));
    if (reduced.count <= 0) {
        return false;
    }

    flint::fmpz_set_si(out, reduced.count);
    return true;
}

}  // namespace detail

namespace {

constexpr ulong kRowModuleModRankPrime = UWORD(27449);
constexpr slong kAnalyticProofAuxPrimeBound = 31;
constexpr slong kCompactMaxPrecision = WORD(1) << 18;
#if FLINT_BITS == 64
constexpr ulong kAuxiliaryPrimeStartPower = 60;
#else
constexpr ulong kAuxiliaryPrimeStartPower = 30;
#endif

struct RowModuleAddResult {
    bool rank_increased = false;
    bool index_refined = false;
};

slong max_slong_value(slong left, slong right) noexcept {
    return left < right ? right : left;
}

slong min_slong_value(slong left, slong right) noexcept {
    return left < right ? left : right;
}

bool same_factor_base(const FactorBase* left, const FactorBase* right) noexcept {
    return left != nullptr && right != nullptr && left->equal(*right);
}

bool valid_relation_source(ClassGroupRelationSource source) noexcept {
    switch (source) {
        case ClassGroupRelationSource::Unknown:
        case ClassGroupRelationSource::Search:
        case ClassGroupRelationSource::RandomProduct:
        case ClassGroupRelationSource::Supplied:
        case ClassGroupRelationSource::Saturation:
        case ClassGroupRelationSource::LargePrimeMatch:
        case ClassGroupRelationSource::ClassGenerator:
            return true;
    }
    return false;
}

bool valid_factor_base_generation_record_status(ProofState status) noexcept {
    return status == ProofState::unavailable ||
           status == ProofState::verified;
}

bool valid_relation_saturation_record_status(ProofState status) noexcept {
    return status == ProofState::unavailable ||
           status == ProofState::verified;
}

bool compact_multiply_power_fmpz(FactoredElement& accumulator,
                                 const FactoredElement& base,
                                 flint::FmpzConstRef exponent) noexcept;

bool multiply_element_power_fmpz(FactoredElement& accumulator,
                                 const Element& base,
                                 flint::FmpzConstRef exponent) noexcept;

bool push_relation_witnesses(FactoredElement& out,
                             const RelationMatrix& relations,
                             flint::FmpzMatConstRef coefficients,
                             slong row) noexcept {
    const Order* order = relations.parent();
    if (order == nullptr || order->parent() == nullptr ||
        out.parent() == nullptr ||
        !out.parent()->has_same_data(*order->parent()) ||
        row < 0 || row >= flint::fmpz_mat_nrows(coefficients) ||
        flint::fmpz_mat_ncols(coefficients) != relations.length()) {
        return false;
    }

    FactoredElement candidate(*order->parent());
    Element generator(*order->parent());
    if (!candidate.is_defined() || !generator.is_defined() ||
        !candidate.one()) {
        return false;
    }

    for (slong i = 0; i < relations.length(); ++i) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(coefficients, row, i);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!relations.generator(generator, i) ||
            !multiply_element_power_fmpz(candidate, generator, exponent)) {
            return false;
        }
    }

    candidate.normalize();
    out.swap(candidate);
    return true;
}

bool relation_exponents_are_zero(bool& is_zero,
                                 const Relation& relation) noexcept {
    flint::FmpzMat row(1, relation.length());
    if (!relation.exponents(flint::FmpzMatRef(row))) {
        return false;
    }

    is_zero = true;
    for (slong j = 0; j < relation.length(); ++j) {
        if (!flint::fmpz_is_zero(
                    flint::fmpz_mat_entry(flint::FmpzMatConstRef(row), 0, j))) {
            is_zero = false;
            return true;
        }
    }
    return true;
}

bool relation_matrix_contains_row(bool& contains,
                                  const RelationMatrix& relations,
                                  flint::FmpzMatConstRef row) noexcept {
    contains = false;
    if (!relations.is_defined() || flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != relations.ncols()) {
        return false;
    }

    flint::FmpzMat stored(1, relations.ncols());
    for (slong i = 0; i < relations.length(); ++i) {
        if (!relations.row(flint::FmpzMatRef(stored), i)) {
            return false;
        }
        if (flint::fmpz_mat_equal(flint::FmpzMatConstRef(stored), row)) {
            contains = true;
            return true;
        }
    }
    return true;
}

enum class DirectLargePrimeStatus {
    unsupported,
    found,
    no_candidate,
    out_of_bound,
};

bool defining_polynomial_is_nice(const NumberField& field) noexcept;

DirectLargePrimeStatus direct_large_prime_from_norm(
        flint::Fmpz& out,
        const Element& generator,
        const FactorBase& base,
        flint::FmpzConstRef factor_base_bound,
        const fmpq* known_norm,
        const fmpz_mat_struct* known_integral_coordinates) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !order->is_maximal() || !order->is_equation_order() ||
        !generator.has_parent(*order->parent()) ||
        flint::fmpz_cmp_ui(factor_base_bound, 2) < 0) {
        return DirectLargePrimeStatus::unsupported;
    }

    if (known_integral_coordinates == nullptr) {
        OrderElement integral_generator(*order);
        if (!integral_generator.is_defined() ||
            !integral_generator.set_element(generator)) {
            return DirectLargePrimeStatus::unsupported;
        }
    } else if (flint::fmpz_mat_nrows(
                       flint::FmpzMatConstRef(known_integral_coordinates)) != 1 ||
               flint::fmpz_mat_ncols(
                       flint::FmpzMatConstRef(known_integral_coordinates)) !=
                       order->degree()) {
        return DirectLargePrimeStatus::unsupported;
    }

    flint::Fmpq norm;
    const fmpq* norm_value = known_norm;
    if (norm_value == nullptr) {
        if (!generator.norm(flint::FmpqRef(norm))) {
            return DirectLargePrimeStatus::unsupported;
        }
        norm_value = norm.raw();
    }
    if (fmpz_is_one(fmpq_denref(norm_value)) == 0) {
        return DirectLargePrimeStatus::unsupported;
    }

    flint::Fmpz remainder;
    fmpz_abs(remainder.raw(), fmpq_numref(norm_value));
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(remainder), 1) <= 0) {
        return DirectLargePrimeStatus::no_candidate;
    }

    flint::Fmpz p;
    for (slong block = 0; block < base.rational_prime_block_count();
         ++block) {
        slong block_length = 0;
        if (!base.rational_prime_block_data(flint::FmpzRef(p),
                                            block_length, block)) {
            return DirectLargePrimeStatus::unsupported;
        }
        if (block_length <= 0) {
            continue;
        }
        while (::fmpz_divisible(remainder.raw(), p.raw()) != 0) {
            ::fmpz_divexact(remainder.raw(), remainder.raw(), p.raw());
        }
        if (flint::fmpz_is_one(flint::FmpzConstRef(remainder))) {
            return DirectLargePrimeStatus::no_candidate;
        }
    }

    flint::Fmpz norm_bound;
    flint::fmpz_mul(flint::FmpzRef(norm_bound), factor_base_bound,
                    factor_base_bound);
    if (flint::fmpz_cmp(flint::FmpzConstRef(remainder),
                        flint::FmpzConstRef(norm_bound)) > 0) {
        return DirectLargePrimeStatus::out_of_bound;
    }
    if (!flint::fmpz_is_prime(flint::FmpzConstRef(remainder))) {
        return DirectLargePrimeStatus::no_candidate;
    }

    flint::fmpz_set(flint::FmpzRef(out),
                    flint::FmpzConstRef(remainder));
    return DirectLargePrimeStatus::found;
}

bool fmpq_poly_set_fmpz_mod_poly_if_integral(
        flint::FmpzModPoly& out,
        const fmpq_poly_struct* polynomial,
        const flint::FmpzModCtx& ctx) noexcept {
    if (!out.is_initialized() || polynomial == nullptr) {
        return false;
    }

    ::fmpz_mod_poly_zero(out.raw(), ctx.raw());
    const slong degree = ::fmpq_poly_degree(polynomial);
    if (degree < 0) {
        return true;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        ::fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (!flint::fmpz_is_one(
                    flint::fmpq_den_ref(flint::FmpqConstRef(coeff)))) {
            return false;
        }
        ::fmpz_mod_poly_set_coeff_fmpz(
                out.raw(), i,
                flint::fmpq_num_ref(flint::FmpqConstRef(coeff)).raw(),
                ctx.raw());
    }
    return true;
}

DirectLargePrimeStatus direct_large_prime_residue_key(
        flint::FmpzPoly& out,
        const Element& generator,
        const FactorBase& base,
        flint::FmpzConstRef p,
        const fmpz_poly_struct* known_integral_polynomial) noexcept {
    const Order* order = base.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !order->is_maximal() ||
        !order->is_equation_order() ||
        !generator.has_parent(*field) || !flint::fmpz_is_prime(p) ||
        !defining_polynomial_is_nice(*field)) {
        return DirectLargePrimeStatus::unsupported;
    }

    const nf_struct* raw_field = field->raw_flint_field();
    if (raw_field == nullptr) {
        return DirectLargePrimeStatus::unsupported;
    }

    flint::FmpzModCtx ctx(p.raw());
    flint::FmpzModPoly element_mod(ctx);
    flint::FmpzModPoly defining_mod(ctx);
    flint::FmpzModPoly gcd_mod(ctx);
    if (!element_mod.is_initialized() || !defining_mod.is_initialized() ||
        !gcd_mod.is_initialized() ||
        !fmpq_poly_set_fmpz_mod_poly_if_integral(
                defining_mod, raw_field->pol, ctx)) {
        return DirectLargePrimeStatus::unsupported;
    }
    if (known_integral_polynomial != nullptr) {
        ::fmpz_mod_poly_set_fmpz_poly(element_mod.raw(),
                                      known_integral_polynomial, ctx.raw());
    } else {
        flint::FmpqPoly element_polynomial;
        if (!generator.get_fmpq_poly(flint::FmpqPolyRef(element_polynomial)) ||
            !fmpq_poly_set_fmpz_mod_poly_if_integral(
                    element_mod, element_polynomial.raw(), ctx)) {
            return DirectLargePrimeStatus::unsupported;
        }
    }

    ::fmpz_mod_poly_gcd(gcd_mod.raw(), element_mod.raw(),
                        defining_mod.raw(), ctx.raw());
    if (::fmpz_mod_poly_degree(gcd_mod.raw(), ctx.raw()) != 1) {
        return DirectLargePrimeStatus::unsupported;
    }

    ::fmpz_mod_poly_get_fmpz_poly(out.raw(), gcd_mod.raw(), ctx.raw());
    return DirectLargePrimeStatus::found;
}

struct DirectPartialResidueHit {
    slong index = -1;
    slong residue_degree = 0;
};

void direct_partial_coordinates_to_polynomial(
        flint::FmpzPoly& out,
        flint::FmpzMatConstRef coordinates) noexcept {
    ::fmpz_poly_zero(out.raw());
    const slong cols = flint::fmpz_mat_ncols(coordinates);
    for (slong j = 0; j < cols; ++j) {
        ::fmpz_poly_set_coeff_fmpz(
                out.raw(), j,
                flint::fmpz_mat_entry(coordinates, 0, j).raw());
    }
}

bool direct_partial_polynomial_is_zero_mod_residue(
        bool& out,
        flint::FmpzModPolyConstRef input,
        flint::FmpzModPolyConstRef modulus,
        flint::FmpzModPolyRef remainder,
        flint::FmpzModCtxConstRef ctx) noexcept {
    out = false;
    if (input.raw() == nullptr || modulus.raw() == nullptr ||
        remainder.raw() == nullptr || ctx.raw() == nullptr ||
        ::fmpz_mod_poly_is_zero(modulus.raw(), ctx.raw()) != 0) {
        return false;
    }

    ::fmpz_mod_poly_rem(remainder.raw(), input.raw(), modulus.raw(),
                        ctx.raw());
    out = ::fmpz_mod_poly_is_zero(remainder.raw(), ctx.raw()) != 0;
    return true;
}

detail::DirectPartialResidueBlockCacheEntry*
direct_partial_residue_block_cache(
        std::vector<detail::DirectPartialResidueBlockCacheEntry>& cache,
        const FactorBase& base,
        slong block_index,
        flint::FmpzConstRef rational_prime,
        slong length) noexcept {
    // reference FactorBaseSingleP owns one finite-ring polynomial context and its
    // residue factors for every rational-prime block in the factor base.
    for (auto& entry : cache) {
        if (entry.block_index == block_index &&
            flint::fmpz_equal(
                    flint::FmpzConstRef(entry.rational_prime),
                    rational_prime)) {
            return &entry;
        }
    }
    if (block_index < 0 || length <= 0) {
        return nullptr;
    }

    cache.emplace_back(block_index, rational_prime);
    auto& entry = cache.back();
    if (entry.context.raw() == nullptr) {
        cache.pop_back();
        return nullptr;
    }
    entry.residue_polynomials.reserve(static_cast<std::size_t>(length));
    entry.prime_indices.reserve(static_cast<std::size_t>(length));
    entry.residue_degrees.reserve(static_cast<std::size_t>(length));
    for (slong offset = 0; offset < length; ++offset) {
        slong index = -1;
        if (!base.rational_prime_block_index(index, block_index, offset)) {
            cache.pop_back();
            return nullptr;
        }
        const PrimeIdeal* prime = base.prime_at(index);
        flint::FmpzPoly residue_polynomial;
        if (prime == nullptr || prime->residue_degree() <= 0 ||
            !prime->residue_polynomial(
                    flint::FmpzPolyRef(residue_polynomial))) {
            cache.pop_back();
            return nullptr;
        }
        entry.residue_polynomials.emplace_back(entry.context);
        if (!entry.residue_polynomials.back().is_initialized()) {
            cache.pop_back();
            return nullptr;
        }
        ::fmpz_mod_poly_set_fmpz_poly(
                entry.residue_polynomials.back().raw(),
                residue_polynomial.raw(), entry.context.raw());
        if (::fmpz_mod_poly_is_zero(
                    entry.residue_polynomials.back().raw(),
                    entry.context.raw()) != 0) {
            cache.pop_back();
            return nullptr;
        }
        entry.prime_indices.push_back(index);
        entry.residue_degrees.push_back(prime->residue_degree());
    }
    return &entry;
}

DirectLargePrimeStatus direct_partial_relation_row(
        flint::FmpzMat& out_row,
        flint::Fmpz& out_rational_large_prime,
        flint::FmpzPoly& out_residue_key,
        const Element& generator,
        const FactorBase& base,
        flint::FmpzConstRef factor_base_bound,
        const fmpq* known_norm,
        const fmpz_mat_struct* known_integral_coordinates,
        const fmpz_poly_struct* known_integral_polynomial,
        std::vector<detail::DirectPartialResidueBlockCacheEntry>&
                residue_block_cache,
        const DiagnosticsContext* diagnostics) noexcept {
    // reference Rel_add.jl stores partial relations by the leftover rational
    // prime and `special_prime_ideal`, then subtracts the two `_factor!`
    // rows on a match.  This is the maximal-equation-order integral-element
    // residue-polynomial branch of that `_factor!` computation.
    auto unsupported = [&](const char* label) noexcept {
        (void)label;
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                            label);
        return DirectLargePrimeStatus::unsupported;
    };

    DirectLargePrimeStatus status = direct_large_prime_from_norm(
            out_rational_large_prime, generator, base, factor_base_bound,
            known_norm, known_integral_coordinates);
    if (status != DirectLargePrimeStatus::found) {
        if (status == DirectLargePrimeStatus::unsupported) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.partial_relation_direct_large_prime_unsupported_norm_screen");
        }
        return status;
    }

    const Order* order = base.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !order->is_maximal() ||
        !order->is_equation_order() || !generator.has_parent(*field) ||
        !defining_polynomial_is_nice(*field)) {
        return unsupported(
                "class_group.partial_relation_direct_large_prime_unsupported_order_shape");
    }

    OrderElement integral_generator(*order);
    if (!integral_generator.is_defined() ||
        (known_integral_coordinates == nullptr
                 ? !integral_generator.set_element(generator)
                 : !integral_generator.set_coordinates(
                           flint::FmpzMatConstRef(
                                   known_integral_coordinates)))) {
        return unsupported(
                "class_group.partial_relation_direct_large_prime_unsupported_integral_generator");
    }

    flint::Fmpq norm;
    const fmpq* norm_value = known_norm;
    if (norm_value == nullptr) {
        if (!generator.norm(flint::FmpqRef(norm))) {
            return unsupported(
                    "class_group.partial_relation_direct_large_prime_unsupported_norm");
        }
        norm_value = norm.raw();
    }
    if (fmpz_is_one(fmpq_denref(norm_value)) == 0) {
        return unsupported(
                "class_group.partial_relation_direct_large_prime_unsupported_norm");
    }

    flint::Fmpz remaining_norm;
    fmpz_abs(remaining_norm.raw(), fmpq_numref(norm_value));

    flint::FmpzMat coordinates(1, order->degree());
    const fmpz_mat_struct* coordinate_rows = known_integral_coordinates;
    if (coordinate_rows == nullptr) {
        if (!integral_generator.get_coordinates(
                    flint::FmpzMatRef(coordinates))) {
            return unsupported(
                    "class_group.partial_relation_direct_large_prime_unsupported_coordinates");
        }
        coordinate_rows = coordinates.raw();
    }

    flint::FmpzPoly input_polynomial_storage;
    const fmpz_poly_struct* input_polynomial = known_integral_polynomial;
    if (input_polynomial == nullptr) {
        direct_partial_coordinates_to_polynomial(
                input_polynomial_storage,
                flint::FmpzMatConstRef(coordinate_rows));
        input_polynomial = input_polynomial_storage.raw();
    }

    flint::FmpzMat candidate(1, base.length());
    flint::fmpz_mat_zero(flint::FmpzMatRef(candidate));

    flint::Fmpz p;
    std::vector<DirectPartialResidueHit> hits;
    for (slong block_index = 0;
         block_index < base.rational_prime_block_count() &&
         !flint::fmpz_is_one(flint::FmpzConstRef(remaining_norm));
         ++block_index) {
        slong length = 0;
        if (!base.rational_prime_block_data(flint::FmpzRef(p), length,
                                            block_index)) {
            return unsupported(
                    "class_group.partial_relation_direct_large_prime_unsupported_block_data");
        }

        slong remaining_vp = 0;
        while (::fmpz_divisible(remaining_norm.raw(), p.raw()) != 0) {
            ::fmpz_divexact(remaining_norm.raw(), remaining_norm.raw(),
                            p.raw());
            ++remaining_vp;
        }
        if (remaining_vp == 0) {
            continue;
        }

        detail::DirectPartialResidueBlockCacheEntry* residue_block =
                direct_partial_residue_block_cache(
                        residue_block_cache, base, block_index,
                        flint::FmpzConstRef(p), length);
        if (residue_block == nullptr ||
            !residue_block->input.is_initialized() ||
            !residue_block->remainder.is_initialized() ||
            residue_block->residue_polynomials.size() !=
                    static_cast<std::size_t>(length) ||
            residue_block->prime_indices.size() !=
                    static_cast<std::size_t>(length) ||
            residue_block->residue_degrees.size() !=
                    static_cast<std::size_t>(length)) {
            return unsupported(
                    "class_group.partial_relation_direct_large_prime_unsupported_residue_block_cache");
        }
        ::fmpz_mod_poly_set_fmpz_poly(residue_block->input.raw(),
                                      input_polynomial,
                                      residue_block->context.raw());

        hits.clear();
        bool block_supported = true;
        for (std::size_t offset = 0;
             offset < residue_block->residue_polynomials.size(); ++offset) {
            bool residue_zero = false;
            if (!direct_partial_polynomial_is_zero_mod_residue(
                        residue_zero,
                        flint::FmpzModPolyConstRef(residue_block->input),
                        flint::FmpzModPolyConstRef(
                                residue_block->residue_polynomials[offset]),
                        flint::FmpzModPolyRef(residue_block->remainder),
                        flint::FmpzModCtxConstRef(
                                residue_block->context))) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.partial_relation_direct_large_prime_unsupported_residue_reduction");
                block_supported = false;
                break;
            }
            if (residue_zero) {
                hits.push_back(DirectPartialResidueHit{
                        residue_block->prime_indices[offset],
                        residue_block->residue_degrees[offset]});
            }
        }

        if (!block_supported) {
            return DirectLargePrimeStatus::unsupported;
        }

        slong residue_hit_vp = 0;
        for (const DirectPartialResidueHit& hit : hits) {
            if (hit.index < 0 || hit.residue_degree <= 0 ||
                residue_hit_vp > WORD_MAX - hit.residue_degree) {
                return unsupported(
                        "class_group.partial_relation_direct_large_prime_unsupported_hit_data");
            }
            residue_hit_vp += hit.residue_degree;
        }

        if (residue_hit_vp == remaining_vp) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.partial_relation_direct_large_prime_path_residue_degree_exact");
            for (const DirectPartialResidueHit& hit : hits) {
                ::fmpz_set_si(flint::fmpz_mat_entry(candidate, 0, hit.index)
                                      .raw(),
                              1);
            }
            remaining_vp = 0;
            continue;
        }

        if (hits.empty()) {
            return unsupported(
                    "class_group.partial_relation_direct_large_prime_unsupported_exact_vp_leftover");
        }

        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.partial_relation_direct_large_prime_path_exact_valuation");
        slong exact_remaining_vp = remaining_vp;
        Ideal exact_valuation_principal(*order);
        bool have_exact_valuation_principal = false;
        for (const DirectPartialResidueHit& hit : hits) {
            const PrimeIdeal* prime = base.prime_at(hit.index);
            slong valuation = -1;
            if (prime == nullptr) {
                return unsupported(
                        "class_group.partial_relation_direct_large_prime_unsupported_exact_valuation");
            }
            flint::FmpzPoly prime_residue_polynomial;
            if (!prime->residue_polynomial(
                        flint::FmpzPolyRef(prime_residue_polynomial))) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.partial_relation_direct_large_prime_path_exact_valuation_no_residue_poly");
            }
            if (!prime->valuation(valuation, integral_generator)) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.partial_relation_direct_large_prime_path_exact_valuation_order_element_unavailable");
                // reference `fb_int_doit` and reference `divide_p_elt` compute exact
                // valuations only for the residue-hit primes over this
                // rational prime.  Until the broader anti-uniformizer element
                // valuation path is ported, reuse the existing exact ideal
                // valuation for just those hits instead of factoring the full
                // principal ideal.
                if (!have_exact_valuation_principal) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.partial_relation_direct_large_prime_path_exact_valuation_principal_ideal");
                    if (!exact_valuation_principal.is_defined() ||
                        !exact_valuation_principal.set_principal(
                                integral_generator)) {
                        return unsupported(
                                "class_group.partial_relation_direct_large_prime_unsupported_exact_valuation_principal_ideal");
                    }
                    have_exact_valuation_principal = true;
                }
                if (!prime->valuation(valuation, exact_valuation_principal)) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.partial_relation_direct_large_prime_unsupported_exact_valuation_kernel");
                    return unsupported(
                            "class_group.partial_relation_direct_large_prime_unsupported_exact_valuation");
                }
            }
            if (prime->ramification_index() != 1) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.partial_relation_direct_large_prime_path_exact_valuation_ramified");
            }
            if (valuation <= 0) {
                continue;
            }
            if (valuation > WORD_MAX / hit.residue_degree) {
                return unsupported(
                        "class_group.partial_relation_direct_large_prime_unsupported_exact_valuation_overflow");
            }
            ::fmpz_set_si(flint::fmpz_mat_entry(candidate, 0, hit.index)
                                  .raw(),
                          valuation);
            exact_remaining_vp -= valuation * hit.residue_degree;
        }
        if (exact_remaining_vp != 0) {
            return unsupported(
                    "class_group.partial_relation_direct_large_prime_unsupported_exact_vp_leftover");
        }
        remaining_vp = 0;
    }

    if (!flint::fmpz_equal(flint::FmpzConstRef(remaining_norm),
                           flint::FmpzConstRef(out_rational_large_prime))) {
        return unsupported(
                "class_group.partial_relation_direct_large_prime_unsupported_leftover_mismatch");
    }

    status = direct_large_prime_residue_key(
            out_residue_key, generator, base,
            flint::FmpzConstRef(out_rational_large_prime), input_polynomial);
    if (status != DirectLargePrimeStatus::found) {
        if (status == DirectLargePrimeStatus::unsupported) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.partial_relation_direct_large_prime_unsupported_residue_key");
        }
        return status;
    }

    out_row = flint::FmpzMat(1, base.length());
    flint::fmpz_mat_set(flint::FmpzMatRef(out_row),
                        flint::FmpzMatConstRef(candidate));
    return DirectLargePrimeStatus::found;
}

bool row_module_add_with_result(RowModuleAddResult& result,
                                fmpz_smat::HnfContext& module,
                                flint::FmpzMatConstRef row,
                                const DiagnosticsContext* diagnostics,
                                bool defer_dependent_refinement = false)
        noexcept {
    module.set_diagnostics(diagnostics);
    result = RowModuleAddResult{};
    if (flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != module.ambient_dim()) {
        return false;
    }

    bool independent = false;
    bool index_refined = false;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.row_module_add_row");
        if (defer_dependent_refinement) {
            if (!module.add_fmpz_mat_row_defer_dependent(&independent, row,
                                                         0)) {
                return false;
            }
        } else {
            if (!module.add_fmpz_mat_row(&independent, &index_refined, row,
                                         0)) {
                return false;
            }
        }
    }
    result.rank_increased = independent;
    result.index_refined = index_refined;

    return true;
}

bool should_keep_dependent_relation(const ClassGroupContext& context,
                                    const RowModuleAddResult& add_result,
                                    slong target_relation_kernel_units)
        noexcept {
    if (add_result.index_refined) {
        return true;
    }
    if (detail::ClassGroupRelationSearchAccess::
                defer_native_goal_publication(context)) {
        return detail::ClassGroupRelationSearchAccess::
                       relation_kernel_row_count(context) <
                target_relation_kernel_units;
    }
    return context.has_presentation() &&
           context.relation_kernel_unit_count() <
                   target_relation_kernel_units;
}

bool relation_kernel_target_open(const ClassGroupContext& context,
                                 slong target_relation_kernel_units) noexcept {
    if (detail::ClassGroupRelationSearchAccess::
                defer_native_goal_publication(context)) {
        return detail::ClassGroupRelationSearchAccess::
                       relation_kernel_row_count(context) <
                target_relation_kernel_units;
    }
    return context.relation_kernel_unit_count() <
            target_relation_kernel_units;
}

bool reset_row_module(fmpz_smat::HnfContext& module, slong ncols) noexcept {
    return module.reset(ncols, kRowModuleModRankPrime);
}

bool reset_relation_modules(
        fmpz_smat::HnfContext& row_module,
        detail::ClassRelationModuleContext& relation_module,
        slong ncols,
        const DiagnosticsContext* diagnostics) noexcept {
    row_module.set_diagnostics(diagnostics);
    relation_module.set_diagnostics(diagnostics);
    return reset_row_module(row_module, ncols) &&
           relation_module.reset(ncols, kRowModuleModRankPrime);
}

bool verify_dlog_full_rank_mod_ell(slong& rank,
                                   flint::FmpzMatConstRef dlog,
                                   flint::FmpzConstRef ell,
                                   slong target_rank) noexcept {
    rank = 0;
    if (!flint::fmpz_is_prime(ell) || target_rank <= 0 ||
        flint::fmpz_mat_nrows(dlog) < target_rank ||
        flint::fmpz_mat_ncols(dlog) <= 0) {
        return false;
    }

    flint::FmpzModCtx ctx(ell.raw());
    flint::FmpzModMat mod_dlog(flint::fmpz_mat_nrows(dlog),
                               flint::fmpz_mat_ncols(dlog), ctx);
    flint::fmpz_mod_mat_set_fmpz_mat(flint::FmpzModMatRef(mod_dlog),
                                     dlog, flint::FmpzModCtxConstRef(ctx));
    rank = flint::fmpz_mod_mat_rank(flint::FmpzModMatRef(mod_dlog),
                                    flint::FmpzModCtxConstRef(ctx));
    return rank >= target_rank;
}

bool dlog_kernel_from_matrix(flint::FmpzMat& out,
                             const flint::FmpzMat& matrix,
                             flint::FmpzConstRef ell) noexcept {
    if (!flint::fmpz_is_prime(ell)) {
        return false;
    }

    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong cols = flint::fmpz_mat_ncols(matrix);
    flint::FmpzModCtx ctx(ell.raw());
    flint::FmpzModMat mod_matrix(rows, cols, ctx);
    flint::FmpzModMat transpose(cols, rows, ctx);
    flint::FmpzModMat nullspace(rows, rows, ctx);
    flint::Fmpz entry;

    fmpz_mod_mat_set_fmpz_mat(mod_matrix.raw(), matrix.raw(), ctx.raw());
    fmpz_mod_mat_transpose(transpose.raw(), mod_matrix.raw(), ctx.raw());
    const slong nullity =
            fmpz_mod_mat_nullspace(nullspace.raw(), transpose.raw(),
                                   ctx.raw());

    flint::FmpzMat candidate(nullity, rows);
    for (slong j = 0; j < nullity; ++j) {
        for (slong i = 0; i < rows; ++i) {
            fmpz_mod_mat_get_entry(entry.raw(), nullspace.raw(), i, j,
                                   ctx.raw());
            flint::fmpz_set(flint::fmpz_mat_entry(candidate, j, i),
                            flint::FmpzConstRef(entry));
        }
    }

    out = std::move(candidate);
    return true;
}

void default_saturation_auxiliary_prime_start(
        flint::Fmpz& out) noexcept {
    fmpz_one(out.raw());
    fmpz_mul_2exp(out.raw(), out.raw(), kAuxiliaryPrimeStartPower);
    flint::fmpz_nextprime(flint::FmpzRef(out), flint::FmpzConstRef(out));
}

slong positive_mod_slong(slong value, slong modulus) noexcept {
    slong result = value % modulus;
    if (result < 0) {
        result += modulus;
    }
    return result;
}

bool mod_exponents(FactoredElement& out,
                         const FactoredElement& input,
                         slong modulus) noexcept {
    if (!out.is_defined() || input.parent() == nullptr ||
        out.parent() == nullptr || modulus <= 1 ||
        !out.parent()->has_same_data(*input.parent()) || !out.one()) {
        return false;
    }

    for (const auto& entry : input.factors()) {
        const slong exponent = positive_mod_slong(entry.exponent, modulus);
        if (exponent != 0 && !out.push(entry.factor, exponent)) {
            return false;
        }
    }
    out.normalize();
    return true;
}

bool defining_polynomial_is_nice(const NumberField& field) noexcept {
    const nf_struct* raw_field = field.raw_flint_field();
    if (raw_field == nullptr) {
        return false;
    }

    const fmpq_poly_struct* polynomial = raw_field->pol;
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree < 1) {
        return false;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (!flint::fmpz_is_one(flint::fmpq_den_ref(coeff))) {
            return false;
        }
    }
    return flint::fmpq_equal_si(flint::FmpqConstRef(coeff), 1);
}

bool order_index_divisor(bool& out,
                               const Order& order,
                               flint::FmpzConstRef q) noexcept {
    out = false;
    const NumberField* field = order.parent();
    if (field == nullptr || !flint::fmpz_is_prime(q) ||
        !defining_polynomial_is_nice(*field)) {
        return true;
    }

    Order equation = Order::equation_order(*field);
    flint::Fmpz index;
    if (!equation.is_defined() ||
        !order_index(flint::FmpzRef(index), equation, order)) {
        return true;
    }

    out = flint::fmpz_divisible(flint::FmpzConstRef(index), q);
    return true;
}

bool next_saturation_auxiliary_prime(
        flint::Fmpz& out,
        flint::FmpzConstRef lower_bound,
        ulong congruence_modulus) noexcept {
    if (congruence_modulus <= 1) {
        return false;
    }

    flint::fmpz_set(flint::FmpzRef(out), lower_bound);
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(out), 2) < 0) {
        flint::fmpz_set_ui(flint::FmpzRef(out), 2);
    }

    while (!flint::fmpz_is_prime(flint::FmpzConstRef(out)) ||
           flint::fmpz_fdiv_ui(flint::FmpzConstRef(out),
                               congruence_modulus) != 1) {
        flint::fmpz_nextprime(flint::FmpzRef(out),
                              flint::FmpzConstRef(out));
    }
    return true;
}

bool advance_saturation_auxiliary_prime(
        flint::Fmpz& out,
        ulong congruence_modulus) noexcept {
    if (congruence_modulus <= 1) {
        return false;
    }
    do {
        flint::fmpz_nextprime(flint::FmpzRef(out),
                              flint::FmpzConstRef(out));
    } while (flint::fmpz_fdiv_ui(flint::FmpzConstRef(out),
                                 congruence_modulus) != 1);
    return true;
}

bool auxiliary_prime_usable(const PrimeIdeal& prime,
                            flint::FmpzConstRef q,
                            flint::FmpzConstRef ell) noexcept;

bool degree_one_roots_for_auxiliary_prime(
        std::vector<flint::Fmpz>& out,
        const Order& order,
        flint::FmpzConstRef q) noexcept {
    out.clear();
    const NumberField* field = order.parent();
    if (field == nullptr || field->raw_flint_field() == nullptr ||
        !flint::fmpz_is_prime(q)) {
        return false;
    }

    flint::FmpzModCtx ctx(q.raw());
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPolyFactor roots(ctx);
    if (!reduced.is_initialized() ||
        !fmpq_poly_set_fmpz_mod_poly_if_integral(
                reduced, field->raw_flint_field()->pol, ctx)) {
        return false;
    }

    // The caller follows reference's saturation scan and has already rejected
    // index and discriminant primes, so the defining polynomial is squarefree
    // modulo q on this fast path.
    fmpz_mod_poly_roots(roots.raw(), reduced.raw(), 0, ctx.raw());
    out.reserve(static_cast<std::size_t>(roots.raw()->num));
    for (slong i = 0; i < roots.raw()->num; ++i) {
        if (fmpz_mod_poly_degree(roots.raw()->poly + i, ctx.raw()) != 1) {
            continue;
        }
        flint::Fmpz root;
        fmpz_mod_poly_get_coeff_fmpz(root.raw(), roots.raw()->poly + i, 0,
                                     ctx.raw());
        flint::fmpz_neg(flint::FmpzRef(root), flint::FmpzConstRef(root));
        ::fmpz_mod(root.raw(), root.raw(), q.raw());
        out.push_back(std::move(root));
    }
    return true;
}

bool degree_one_quotient_log_row_at_root(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& relations,
        flint::FmpzConstRef q,
        flint::FmpzConstRef root,
        slong n) noexcept {
    if (relations.empty() || n <= 1 ||
        !n_is_prime(static_cast<ulong>(n)) || !flint::fmpz_is_prime(q)) {
        return false;
    }

    flint::Fmpz ell;
    flint::fmpz_set_si(flint::FmpzRef(ell), n);
    flint::Fmpz cofactor;
    flint::Fmpz quotient_generator;
    if (!detail::quotient_log_mod_prime_setup(
                cofactor, quotient_generator, q, flint::FmpzConstRef(ell))) {
        return false;
    }

    flint::FmpzMat row(1, static_cast<slong>(relations.size()));
    flint::Fmpz value;
    for (slong j = 0; j < static_cast<slong>(relations.size()); ++j) {
        if (!detail::factored_value_at_degree_one_root(
                    value, relations[static_cast<std::size_t>(j)], q, root) ||
            !detail::quotient_log_mod_prime_apply(
                    flint::fmpz_mat_entry(row, 0, j),
                    flint::FmpzConstRef(value), flint::FmpzConstRef(cofactor),
                    flint::FmpzConstRef(quotient_generator), q,
                    flint::FmpzConstRef(ell))) {
            return false;
        }
    }

    out = std::move(row);
    return true;
}

bool reduce_saturation_candidate_space(flint::FmpzMat& candidates,
                                             flint::FmpzMatConstRef log_row,
                                             slong n) noexcept {
    const slong rows = flint::fmpz_mat_nrows(candidates);
    const slong cols = flint::fmpz_mat_ncols(candidates);
    if (n <= 1 || rows <= 0 || cols <= 0 ||
        flint::fmpz_mat_nrows(log_row) != 1 ||
        flint::fmpz_mat_ncols(log_row) != rows) {
        return false;
    }

    flint::Fmpz n_fmpz;
    flint::fmpz_set_si(flint::FmpzRef(n_fmpz), n);
    flint::FmpzModCtx ctx(n_fmpz.raw());
    flint::FmpzModMat mod_log(1, rows, ctx);
    flint::FmpzModMat mod_candidates(rows, cols, ctx);
    flint::FmpzModMat reduced_log(1, cols, ctx);
    flint::FmpzModMat kernel(cols, cols, ctx);

    fmpz_mod_mat_set_fmpz_mat(mod_log.raw(), log_row.raw(), ctx.raw());
    fmpz_mod_mat_set_fmpz_mat(mod_candidates.raw(), candidates.raw(),
                              ctx.raw());
    fmpz_mod_mat_mul(reduced_log.raw(), mod_log.raw(),
                     mod_candidates.raw(), ctx.raw());
    // reference `Saturate.jl:compute_candidates_for_saturate` computes
    // `kernel(z, side = :right)` for the row `z = z*A`.
    const slong nullity = fmpz_mod_mat_nullspace(
            kernel.raw(), reduced_log.raw(), ctx.raw());
    if (nullity == 0) {
        candidates = flint::FmpzMat(0, rows);
        return true;
    }

    flint::FmpzModMat kernel_basis(cols, nullity, ctx);
    flint::Fmpz entry;
    for (slong i = 0; i < cols; ++i) {
        for (slong j = 0; j < nullity; ++j) {
            fmpz_mod_mat_get_entry(entry.raw(), kernel.raw(), i, j,
                                   ctx.raw());
            fmpz_mod_mat_set_entry(kernel_basis.raw(), i, j, entry.raw(),
                                   ctx.raw());
        }
    }

    flint::FmpzModMat next(rows, nullity, ctx);
    fmpz_mod_mat_mul(next.raw(), mod_candidates.raw(),
                     kernel_basis.raw(), ctx.raw());
    flint::FmpzMat lifted(rows, nullity);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < nullity; ++j) {
            fmpz_mod_mat_get_entry(flint::fmpz_mat_entry(lifted, i, j).raw(),
                                   next.raw(), i, j, ctx.raw());
        }
    }

    candidates = std::move(lifted);
    return true;
}

bool saturation_candidates_from_relations(
        flint::FmpzMat& out,
        slong& auxiliary_prime_count,
        const std::vector<FactoredElement>& relations,
        const Order& order,
        slong n,
        double stable,
        flint::FmpzConstRef auxiliary_prime_start,
        slong max_auxiliary_primes,
        detail::SaturationCandidateSearchResult* search_result =
                nullptr,
        const DiagnosticsContext* diagnostics = nullptr) noexcept {
    auxiliary_prime_count = 0;
    if (search_result != nullptr) {
        *search_result = detail::SaturationCandidateSearchResult{};
    }
    if (relations.empty() || order.parent() == nullptr ||
        !order.is_maximal() || n <= 1 ||
        !n_is_prime(static_cast<ulong>(n)) || !std::isfinite(stable) ||
        stable <= 0.0 || max_auxiliary_primes < 0) {
        return false;
    }

    const NumberField& field = *order.parent();
    std::vector<FactoredElement> reduced_relations;
    reduced_relations.reserve(relations.size() + 1);
    for (const FactoredElement& relation : relations) {
        FactoredElement reduced(field);
        if (!reduced.is_defined() ||
            !mod_exponents(reduced, relation, n)) {
            return false;
        }
        reduced_relations.push_back(std::move(reduced));
    }

    flint::Fmpz torsion_order;
    if (!root_of_unity_order(flint::FmpzRef(torsion_order), field)) {
        return false;
    }
    if (flint::fmpz_fdiv_ui(flint::FmpzConstRef(torsion_order),
                            static_cast<ulong>(n)) == 0) {
        Element zeta(field);
        FactoredElement torsion(field);
        if (!zeta.is_defined() ||
            !root_of_unity_generator(zeta, field) ||
            !torsion.is_defined() || !torsion.set_element(zeta)) {
            return false;
        }
        reduced_relations.push_back(std::move(torsion));
    }

    const slong input_count =
            static_cast<slong>(reduced_relations.size());
    flint::FmpzMat candidates(input_count, input_count);
    flint::fmpz_mat_one(flint::FmpzMatRef(candidates));

    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }

    const double threshold = stable * static_cast<double>(input_count);
    slong stable_count = 1;
    slong current_cols = input_count;
    flint::Fmpz q;
    if (!next_saturation_auxiliary_prime(
                q, auxiliary_prime_start, static_cast<ulong>(n))) {
        return false;
    }

    while (auxiliary_prime_count < max_auxiliary_primes) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_saturation.auxiliary_prime_candidate");
        bool is_index_divisor = false;
        if (!order_index_divisor(is_index_divisor, order,
                                       flint::FmpzConstRef(q))) {
            return false;
        }
        if (is_index_divisor ||
            flint::fmpz_divisible(flint::FmpzConstRef(discriminant),
                                  flint::FmpzConstRef(q))) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_saturation.auxiliary_prime_rejected");
            if (!advance_saturation_auxiliary_prime(
                        q, static_cast<ulong>(n))) {
                return false;
            }
            continue;
        }

        std::vector<flint::Fmpz> degree_one_roots;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_saturation.auxiliary_prime_degree_one_"
                    "roots");
            // This is the degree-one specialization of reference
            // `prime_decomposition(OK, q, 1)` for the already-filtered
            // non-index, non-discriminant auxiliary prime.
            if (!degree_one_roots_for_auxiliary_prime(
                        degree_one_roots, order, flint::FmpzConstRef(q))) {
                return false;
            }
        }
        ++auxiliary_prime_count;
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_saturation.auxiliary_prime_used");
        if (search_result != nullptr) {
            search_result->auxiliary_prime_count = auxiliary_prime_count;
        }
        for (const flint::Fmpz& root : degree_one_roots) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_saturation.degree_one_prime");
            if (search_result != nullptr) {
                ++search_result->degree_one_prime_count;
            }
            flint::FmpzMat log_row(0, 0);
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_saturation.quotient_log_row_"
                        "degree_one");
                if (!degree_one_quotient_log_row_at_root(
                            log_row, reduced_relations,
                            flint::FmpzConstRef(q),
                            flint::FmpzConstRef(root), n)) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.relation_saturation.quotient_log_row_"
                            "skipped");
                    if (search_result != nullptr) {
                        ++search_result->skipped_quotient_log_row_count;
                    }
                    continue;
                }
            }
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_saturation.quotient_log_row");
            if (search_result != nullptr) {
                ++search_result->quotient_log_row_count;
            }
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_saturation.reduce_candidate_space");
                if (!reduce_saturation_candidate_space(
                            candidates, flint::FmpzMatConstRef(log_row), n)) {
                    return false;
                }
            }
            if (flint::fmpz_mat_nrows(candidates) == 0) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_saturation.candidate_space_empty");
                if (search_result != nullptr) {
                    search_result->collapsed_to_zero = true;
                    search_result->collapse_after_quotient_log_rows =
                            search_result->quotient_log_row_count;
                    search_result->final_candidate_rows = 0;
                    search_result->final_candidate_cols =
                            flint::fmpz_mat_ncols(
                                    flint::FmpzMatConstRef(candidates));
                    search_result->final_stable_count = stable_count;
                }
                out = std::move(candidates);
                return true;
            }

            const slong next_cols = flint::fmpz_mat_ncols(candidates);
            if (current_cols == next_cols) {
                ++stable_count;
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_saturation.candidate_space_stable");
            } else {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_saturation.candidate_space_reduced");
                stable_count = 0;
                current_cols = next_cols;
            }
            if (static_cast<double>(stable_count) > threshold) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_saturation.stable_threshold");
                if (search_result != nullptr) {
                    search_result->final_candidate_rows =
                            flint::fmpz_mat_nrows(
                                    flint::FmpzMatConstRef(candidates));
                    search_result->final_candidate_cols =
                            flint::fmpz_mat_ncols(
                                    flint::FmpzMatConstRef(candidates));
                    search_result->final_stable_count = stable_count;
                }
                out = std::move(candidates);
                return true;
            }
        }

        if (static_cast<double>(stable_count) > threshold) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_saturation.stable_threshold");
            if (search_result != nullptr) {
                search_result->final_candidate_rows =
                        flint::fmpz_mat_nrows(
                                flint::FmpzMatConstRef(candidates));
                search_result->final_candidate_cols =
                        flint::fmpz_mat_ncols(
                                flint::FmpzMatConstRef(candidates));
                search_result->final_stable_count = stable_count;
            }
            out = std::move(candidates);
            return true;
        }
        if (!advance_saturation_auxiliary_prime(
                    q, static_cast<ulong>(n))) {
            return false;
        }
    }

    if (search_result != nullptr) {
        search_result->final_candidate_rows =
                flint::fmpz_mat_nrows(flint::FmpzMatConstRef(candidates));
        search_result->final_candidate_cols =
                flint::fmpz_mat_ncols(flint::FmpzMatConstRef(candidates));
        search_result->final_stable_count = stable_count;
    }
    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.relation_saturation.auxiliary_prime_limit");
    out = std::move(candidates);
    return true;
}

bool auxiliary_prime_usable(const PrimeIdeal& prime,
                            flint::FmpzConstRef q,
                            flint::FmpzConstRef ell) noexcept {
    const slong f = prime.residue_degree();
    if (f <= 0 || !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::Fmpz order;
    flint::fmpz_pow_ui(flint::FmpzRef(order), q, static_cast<ulong>(f));
    flint::fmpz_sub_ui(flint::FmpzRef(order), flint::FmpzConstRef(order), 1);
    return flint::fmpz_divisible(flint::FmpzConstRef(order), ell);
}

bool compact_multiply_power_fmpz(FactoredElement& accumulator,
                                 const FactoredElement& base,
                                 flint::FmpzConstRef exponent) noexcept {
    if (flint::fmpz_is_zero(exponent)) {
        return true;
    }

    FactoredElement power(*accumulator.parent());
    FactoredElement product(*accumulator.parent());
    if (flint::fmpz_fits_si(exponent)) {
        return power.pow_si(base, flint::fmpz_get_si(exponent)) &&
               product.multiply(accumulator, power) &&
               accumulator.set(product);
    }

    Element base_value(*accumulator.parent());
    Element powered_value(*accumulator.parent());
    return base.evaluate(base_value) &&
           powered_value.pow_fmpz(base_value, exponent) &&
           power.set_element(powered_value) &&
           product.multiply(accumulator, power) &&
           accumulator.set(product);
}

bool multiply_element_power_fmpz(FactoredElement& accumulator,
                                 const Element& base,
                                 flint::FmpzConstRef exponent) noexcept {
    if (flint::fmpz_is_zero(exponent)) {
        return true;
    }
    if (flint::fmpz_fits_si(exponent)) {
        return accumulator.push(base, flint::fmpz_get_si(exponent));
    }

    FactoredElement factored_base(*accumulator.parent());
    return factored_base.set_element(base) &&
           compact_multiply_power_fmpz(accumulator, factored_base, exponent);
}

bool fmpz_mat_row_is_zero(bool& out,
                          flint::FmpzMatConstRef row) noexcept {
    if (flint::fmpz_mat_nrows(row) != 1) {
        return false;
    }

    out = true;
    for (slong j = 0; j < flint::fmpz_mat_ncols(row); ++j) {
        if (!flint::fmpz_is_zero(flint::fmpz_mat_entry(row, 0, j))) {
            out = false;
            return true;
        }
    }
    return true;
}

bool fmpz_mat_copy_row(flint::FmpzMatRef out,
                       slong out_row,
                       flint::FmpzMatConstRef in,
                       slong in_row) noexcept {
    if (out_row < 0 || out_row >= flint::fmpz_mat_nrows(out) ||
        in_row < 0 || in_row >= flint::fmpz_mat_nrows(in) ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(in)) {
        return false;
    }
    for (slong j = 0; j < flint::fmpz_mat_ncols(in); ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, out_row, j),
                        flint::fmpz_mat_entry(in, in_row, j));
    }
    return true;
}

bool hnf_row_first_nonzero_is_one(bool& out,
                                  flint::FmpzMatConstRef rows,
                                  slong row) noexcept {
    if (row < 0 || row >= flint::fmpz_mat_nrows(rows)) {
        return false;
    }
    out = false;
    for (slong j = 0; j < flint::fmpz_mat_ncols(rows); ++j) {
        flint::FmpzConstRef entry = flint::fmpz_mat_entry(rows, row, j);
        if (flint::fmpz_is_zero(entry)) {
            continue;
        }
        out = flint::fmpz_is_one(entry);
        return true;
    }
    return true;
}

slong fmpz_mat_row_first_nonzero_col(flint::FmpzMatConstRef rows,
                                     slong row,
                                     slong start = 0) noexcept {
    if (row < 0 || row >= flint::fmpz_mat_nrows(rows)) {
        return -1;
    }

    const slong cols = flint::fmpz_mat_ncols(rows);
    for (slong col = max_slong_value(WORD(0), start); col < cols; ++col) {
        if (!flint::fmpz_is_zero(flint::fmpz_mat_entry(rows, row, col))) {
            return col;
        }
    }
    return cols;
}

bool fmpz_mat_row_is_negative(flint::FmpzMatConstRef rows,
                              slong row) noexcept {
    const slong col = fmpz_mat_row_first_nonzero_col(rows, row);
    return col >= 0 && col < flint::fmpz_mat_ncols(rows) &&
           flint::fmpz_sgn(flint::fmpz_mat_entry(rows, row, col)) < 0;
}

void fmpz_mat_neg_row(flint::FmpzMatRef rows, slong row) noexcept {
    for (slong col = 0; col < flint::fmpz_mat_ncols(rows); ++col) {
        flint::FmpzRef entry = flint::fmpz_mat_entry(rows, row, col);
        ::fmpz_neg(entry.raw(), entry.raw());
    }
}

void fmpz_mat_zero_row(flint::FmpzMatRef rows, slong row) noexcept {
    for (slong col = 0; col < flint::fmpz_mat_ncols(rows); ++col) {
        flint::fmpz_zero(flint::fmpz_mat_entry(rows, row, col));
    }
}

bool fmpz_mat_copy_row_between(flint::FmpzMatRef out,
                               slong out_row,
                               flint::FmpzMatConstRef in,
                               slong in_row) noexcept {
    if (out_row < 0 || out_row >= flint::fmpz_mat_nrows(out) ||
        in_row < 0 || in_row >= flint::fmpz_mat_nrows(in) ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(in)) {
        return false;
    }
    for (slong col = 0; col < flint::fmpz_mat_ncols(in); ++col) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, out_row, col),
                        flint::fmpz_mat_entry(in, in_row, col));
    }
    return true;
}

bool fmpz_mat_add_scaled_row(flint::FmpzMatRef out,
                             slong out_row,
                             flint::FmpzMatConstRef in,
                             slong in_row,
                             flint::FmpzConstRef scalar) noexcept {
    if (out_row < 0 || out_row >= flint::fmpz_mat_nrows(out) ||
        in_row < 0 || in_row >= flint::fmpz_mat_nrows(in) ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(in)) {
        return false;
    }
    if (flint::fmpz_is_zero(scalar)) {
        return true;
    }

    for (slong col = 0; col < flint::fmpz_mat_ncols(in); ++col) {
        flint::fmpz_addmul(flint::fmpz_mat_entry(out, out_row, col),
                           scalar,
                           flint::fmpz_mat_entry(in, in_row, col));
    }
    return true;
}

bool fmpz_mat_transform_row_pair(flint::FmpzMatRef first,
                                 slong first_row,
                                 flint::FmpzMatRef second,
                                 slong second_row,
                                 flint::FmpzConstRef a,
                                 flint::FmpzConstRef b,
                                 flint::FmpzConstRef c,
                                 flint::FmpzConstRef d) noexcept {
    if (first_row < 0 || first_row >= flint::fmpz_mat_nrows(first) ||
        second_row < 0 || second_row >= flint::fmpz_mat_nrows(second) ||
        flint::fmpz_mat_ncols(first) != flint::fmpz_mat_ncols(second)) {
        return false;
    }

    flint::Fmpz first_next;
    flint::Fmpz second_next;
    for (slong col = 0; col < flint::fmpz_mat_ncols(first); ++col) {
        flint::FmpzConstRef first_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(first.raw()),
                        first_row, col);
        flint::FmpzConstRef second_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(second.raw()),
                        second_row, col);

        flint::fmpz_mul(flint::FmpzRef(first_next), a, first_entry);
        flint::fmpz_addmul(flint::FmpzRef(first_next), b, second_entry);
        flint::fmpz_mul(flint::FmpzRef(second_next), c, first_entry);
        flint::fmpz_addmul(flint::FmpzRef(second_next), d, second_entry);

        flint::fmpz_set(flint::fmpz_mat_entry(first, first_row, col),
                        flint::FmpzConstRef(first_next));
        flint::fmpz_set(flint::fmpz_mat_entry(second, second_row, col),
                        flint::FmpzConstRef(second_next));
    }
    return true;
}

slong dense_hnf_find_row_starting_with(flint::FmpzMatConstRef rows,
                                             slong active_rows,
                                             slong pivot_col) noexcept {
    for (slong row = 0; row < active_rows; ++row) {
        const slong pivot = fmpz_mat_row_first_nonzero_col(rows, row);
        if (pivot < 0) {
            return row;
        }
        if (pivot >= pivot_col) {
            return row;
        }
    }
    return active_rows;
}

bool dense_hnf_reduce_right_current(
        flint::FmpzMatRef basis,
        flint::FmpzMatRef coeffs,
        slong active_rows,
        flint::FmpzMatRef current,
        flint::FmpzMatRef current_coeff,
        slong start_col) noexcept {
    if (flint::fmpz_mat_nrows(current) != 1 ||
        flint::fmpz_mat_nrows(current_coeff) != 1 ||
        flint::fmpz_mat_ncols(basis) != flint::fmpz_mat_ncols(current) ||
        flint::fmpz_mat_ncols(coeffs) !=
                flint::fmpz_mat_ncols(current_coeff)) {
        return false;
    }

    flint::Fmpz quotient;
    flint::Fmpz remainder;
    for (slong col = max_slong_value(WORD(0), start_col);
         col < flint::fmpz_mat_ncols(current);) {
        const slong current_col =
                fmpz_mat_row_first_nonzero_col(
                        flint::FmpzMatConstRef(current.raw()), 0, col);
        if (current_col >= flint::fmpz_mat_ncols(current)) {
            break;
        }

        const slong basis_row = dense_hnf_find_row_starting_with(
                flint::FmpzMatConstRef(basis.raw()), active_rows,
                current_col);
        if (basis_row >= active_rows) {
            break;
        }

        const slong basis_col = fmpz_mat_row_first_nonzero_col(
                flint::FmpzMatConstRef(basis.raw()), basis_row);
        if (basis_col > current_col) {
            col = current_col + 1;
            continue;
        }
        if (basis_col != current_col) {
            return false;
        }

        flint::FmpzConstRef current_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(current.raw()), 0,
                        current_col);
        flint::FmpzConstRef basis_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(basis.raw()), basis_row,
                        basis_col);
        ::fmpz_tdiv_qr(quotient.raw(), remainder.raw(),
                       current_entry.raw(), basis_entry.raw());
        if (flint::fmpz_sgn(flint::FmpzConstRef(remainder)) < 0) {
            flint::fmpz_sub_ui(flint::FmpzRef(quotient),
                               flint::FmpzConstRef(quotient), 1);
            flint::fmpz_add(flint::FmpzRef(remainder),
                            flint::FmpzConstRef(remainder), basis_entry);
        }

        if (!flint::fmpz_is_zero(flint::FmpzConstRef(quotient))) {
            flint::fmpz_neg(flint::FmpzRef(quotient),
                            flint::FmpzConstRef(quotient));
            if (!fmpz_mat_add_scaled_row(
                        current, 0,
                        flint::FmpzMatConstRef(basis.raw()), basis_row,
                        flint::FmpzConstRef(quotient)) ||
                !fmpz_mat_add_scaled_row(
                        current_coeff, 0,
                        flint::FmpzMatConstRef(coeffs.raw()), basis_row,
                        flint::FmpzConstRef(quotient))) {
                return false;
            }
        }

        if (flint::fmpz_is_zero(flint::FmpzConstRef(remainder))) {
            col = current_col;
        } else {
            col = current_col + 1;
        }
    }

    if (fmpz_mat_row_is_negative(flint::FmpzMatConstRef(current.raw()), 0)) {
        fmpz_mat_neg_row(current, 0);
        fmpz_mat_neg_row(current_coeff, 0);
    }
    return true;
}

bool dense_hnf_reduce_right_basis_row(flint::FmpzMatRef basis,
                                            flint::FmpzMatRef coeffs,
                                            slong active_rows,
                                            slong row,
                                            slong start_col) noexcept {
    flint::Fmpz quotient;
    flint::Fmpz remainder;
    for (slong col = max_slong_value(WORD(0), start_col);
         col < flint::fmpz_mat_ncols(basis);) {
        const slong current_col =
                fmpz_mat_row_first_nonzero_col(
                        flint::FmpzMatConstRef(basis.raw()), row, col);
        if (current_col >= flint::fmpz_mat_ncols(basis)) {
            break;
        }

        const slong reducer = dense_hnf_find_row_starting_with(
                flint::FmpzMatConstRef(basis.raw()), active_rows,
                current_col);
        if (reducer >= active_rows) {
            break;
        }
        const slong reducer_col = fmpz_mat_row_first_nonzero_col(
                flint::FmpzMatConstRef(basis.raw()), reducer);
        if (reducer_col > current_col) {
            col = current_col + 1;
            continue;
        }
        if (reducer_col != current_col || reducer == row) {
            return false;
        }

        flint::FmpzConstRef current_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(basis.raw()), row,
                        current_col);
        flint::FmpzConstRef reducer_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(basis.raw()), reducer,
                        reducer_col);
        ::fmpz_tdiv_qr(quotient.raw(), remainder.raw(),
                       current_entry.raw(), reducer_entry.raw());
        if (flint::fmpz_sgn(flint::FmpzConstRef(remainder)) < 0) {
            flint::fmpz_sub_ui(flint::FmpzRef(quotient),
                               flint::FmpzConstRef(quotient), 1);
            flint::fmpz_add(flint::FmpzRef(remainder),
                            flint::FmpzConstRef(remainder), reducer_entry);
        }
        if (!flint::fmpz_is_zero(flint::FmpzConstRef(quotient))) {
            flint::fmpz_neg(flint::FmpzRef(quotient),
                            flint::FmpzConstRef(quotient));
            if (!fmpz_mat_add_scaled_row(
                        basis, row,
                        flint::FmpzMatConstRef(basis.raw()), reducer,
                        flint::FmpzConstRef(quotient)) ||
                !fmpz_mat_add_scaled_row(
                        coeffs, row,
                        flint::FmpzMatConstRef(coeffs.raw()), reducer,
                        flint::FmpzConstRef(quotient))) {
                return false;
            }
        }

        if (flint::fmpz_is_zero(flint::FmpzConstRef(remainder))) {
            col = current_col;
        } else {
            col = current_col + 1;
        }
    }

    if (fmpz_mat_row_is_negative(flint::FmpzMatConstRef(basis.raw()), row)) {
        fmpz_mat_neg_row(basis, row);
        fmpz_mat_neg_row(coeffs, row);
    }
    return true;
}

bool dense_hnf_reduce_full(flint::FmpzMatRef basis,
                                 flint::FmpzMatRef coeffs,
                                 slong active_rows,
                                 flint::FmpzMatRef current,
                                 flint::FmpzMatRef current_coeff,
                                 std::vector<slong>& changed_pivots) noexcept {
    changed_pivots.clear();
    if (flint::fmpz_mat_nrows(current) != 1 ||
        flint::fmpz_mat_nrows(current_coeff) != 1 ||
        flint::fmpz_mat_ncols(basis) != flint::fmpz_mat_ncols(current) ||
        flint::fmpz_mat_ncols(coeffs) !=
                flint::fmpz_mat_ncols(current_coeff)) {
        return false;
    }

    if (active_rows == 0) {
        if (fmpz_mat_row_is_negative(
                    flint::FmpzMatConstRef(current.raw()), 0)) {
            fmpz_mat_neg_row(current, 0);
            fmpz_mat_neg_row(current_coeff, 0);
        }
        return true;
    }

    flint::Fmpz quotient;
    flint::Fmpz remainder;
    flint::Fmpz gcd;
    flint::Fmpz a;
    flint::Fmpz b;
    flint::Fmpz c;
    flint::Fmpz d;

    while (fmpz_mat_row_first_nonzero_col(
                   flint::FmpzMatConstRef(current.raw()), 0) <
           flint::fmpz_mat_ncols(current)) {
        const slong pivot = fmpz_mat_row_first_nonzero_col(
                flint::FmpzMatConstRef(current.raw()), 0);
        const slong basis_row = dense_hnf_find_row_starting_with(
                flint::FmpzMatConstRef(basis.raw()), active_rows, pivot);
        if (basis_row >= active_rows) {
            if (fmpz_mat_row_is_negative(
                        flint::FmpzMatConstRef(current.raw()), 0)) {
                fmpz_mat_neg_row(current, 0);
                fmpz_mat_neg_row(current_coeff, 0);
            }
            return dense_hnf_reduce_right_current(
                    basis, coeffs, active_rows, current, current_coeff, 0);
        }

        const slong basis_pivot = fmpz_mat_row_first_nonzero_col(
                flint::FmpzMatConstRef(basis.raw()), basis_row);
        if (basis_pivot > pivot) {
            if (fmpz_mat_row_is_negative(
                        flint::FmpzMatConstRef(current.raw()), 0)) {
                fmpz_mat_neg_row(current, 0);
                fmpz_mat_neg_row(current_coeff, 0);
            }
            return dense_hnf_reduce_right_current(
                    basis, coeffs, active_rows, current, current_coeff, 0);
        }
        if (basis_pivot != pivot) {
            return false;
        }

        flint::FmpzConstRef current_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(current.raw()), 0, pivot);
        flint::FmpzConstRef basis_entry =
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(basis.raw()), basis_row,
                        basis_pivot);
        ::fmpz_tdiv_qr(quotient.raw(), remainder.raw(),
                       current_entry.raw(), basis_entry.raw());
        if (flint::fmpz_is_zero(flint::FmpzConstRef(remainder))) {
            if (flint::fmpz_is_zero(flint::FmpzConstRef(quotient))) {
                return false;
            }
            flint::fmpz_neg(flint::FmpzRef(quotient),
                            flint::FmpzConstRef(quotient));
            if (!fmpz_mat_add_scaled_row(
                        current, 0,
                        flint::FmpzMatConstRef(basis.raw()), basis_row,
                        flint::FmpzConstRef(quotient)) ||
                !fmpz_mat_add_scaled_row(
                        current_coeff, 0,
                        flint::FmpzMatConstRef(coeffs.raw()), basis_row,
                        flint::FmpzConstRef(quotient))) {
                return false;
            }
            continue;
        }

        ::fmpz_xgcd_canonical_bezout(gcd.raw(), a.raw(), b.raw(),
                                      basis_entry.raw(),
                                      current_entry.raw());
        if (flint::fmpz_sgn(flint::FmpzConstRef(gcd)) <= 0) {
            return false;
        }
        flint::fmpz_divexact(flint::FmpzRef(c), current_entry,
                             flint::FmpzConstRef(gcd));
        flint::fmpz_neg(flint::FmpzRef(c), flint::FmpzConstRef(c));
        flint::fmpz_divexact(flint::FmpzRef(d), basis_entry,
                             flint::FmpzConstRef(gcd));

        if (!fmpz_mat_transform_row_pair(
                    basis, basis_row, current, 0,
                    flint::FmpzConstRef(a), flint::FmpzConstRef(b),
                    flint::FmpzConstRef(c), flint::FmpzConstRef(d)) ||
            !fmpz_mat_transform_row_pair(
                    coeffs, basis_row, current_coeff, 0,
                    flint::FmpzConstRef(a), flint::FmpzConstRef(b),
                    flint::FmpzConstRef(c), flint::FmpzConstRef(d))) {
            return false;
        }

        changed_pivots.push_back(basis_pivot);
        if (!dense_hnf_reduce_right_basis_row(
                    basis, coeffs, active_rows, basis_row,
                    basis_pivot + 1)) {
            return false;
        }
    }

    if (!dense_hnf_reduce_right_current(
                basis, coeffs, active_rows, current, current_coeff, 0)) {
        return false;
    }
    if (fmpz_mat_row_is_negative(flint::FmpzMatConstRef(current.raw()), 0)) {
        fmpz_mat_neg_row(current, 0);
        fmpz_mat_neg_row(current_coeff, 0);
    }
    return true;
}

bool dense_hnf_reduce_up(flint::FmpzMatRef basis,
                               flint::FmpzMatRef coeffs,
                               slong active_rows,
                               std::vector<slong>& pivots) noexcept {
    if (pivots.empty()) {
        return true;
    }

    std::sort(pivots.begin(), pivots.end());
    const slong first_pivot = pivots.front();
    const slong last_pivot = pivots.back();
    const slong stop_row = dense_hnf_find_row_starting_with(
            flint::FmpzMatConstRef(basis.raw()), active_rows, last_pivot);
    for (slong row = stop_row - 1; row >= 0; --row) {
        const slong pivot = fmpz_mat_row_first_nonzero_col(
                flint::FmpzMatConstRef(basis.raw()), row);
        if (pivot < 0 || pivot >= flint::fmpz_mat_ncols(basis)) {
            return false;
        }
        if (!dense_hnf_reduce_right_basis_row(
                    basis, coeffs, active_rows, row,
                    max_slong_value(pivot + 1, first_pivot))) {
            return false;
        }
    }
    return true;
}

bool dense_hnf_insert_row(flint::FmpzMatRef basis,
                                flint::FmpzMatRef coeffs,
                                slong& active_rows,
                                slong position,
                                flint::FmpzMatConstRef current,
                                flint::FmpzMatConstRef current_coeff) noexcept {
    if (position < 0 || position > active_rows ||
        active_rows >= flint::fmpz_mat_nrows(basis)) {
        return false;
    }

    for (slong row = active_rows; row > position; --row) {
        if (!fmpz_mat_copy_row_between(
                    basis, row, flint::FmpzMatConstRef(basis.raw()),
                    row - 1) ||
            !fmpz_mat_copy_row_between(
                    coeffs, row, flint::FmpzMatConstRef(coeffs.raw()),
                    row - 1)) {
            return false;
        }
    }
    if (!fmpz_mat_copy_row_between(basis, position, current, 0) ||
        !fmpz_mat_copy_row_between(coeffs, position, current_coeff, 0)) {
        return false;
    }
    ++active_rows;
    return true;
}

bool dense_hnf_kannan_bachem_transform(
        flint::FmpzMat& hnf,
        flint::FmpzMat& transform,
        flint::FmpzMatConstRef input) noexcept {
    // Source trace: reference `Sparse/HNF.jl:hnf_kannan_bachem` with
    // `with_transform = Val(true)`, `truncate = true`, `full_hnf = true`.
    const slong input_rows = flint::fmpz_mat_nrows(input);
    const slong input_cols = flint::fmpz_mat_ncols(input);
    flint::FmpzMat basis(input_rows, input_cols);
    flint::FmpzMat coeffs(input_rows, input_rows);
    flint::FmpzMat current(1, input_cols);
    flint::FmpzMat current_coeff(1, input_rows);
    std::vector<slong> changed_pivots;
    slong active_rows = 0;

    for (slong row = 0; row < input_rows; ++row) {
        if (!fmpz_mat_copy_row_between(
                    flint::FmpzMatRef(current), 0, input, row)) {
            return false;
        }
        fmpz_mat_zero_row(flint::FmpzMatRef(current_coeff), 0);
        flint::fmpz_one(flint::fmpz_mat_entry(
                flint::FmpzMatRef(current_coeff), 0, row));

        if (!dense_hnf_reduce_full(
                    flint::FmpzMatRef(basis),
                    flint::FmpzMatRef(coeffs), active_rows,
                    flint::FmpzMatRef(current),
                    flint::FmpzMatRef(current_coeff), changed_pivots)) {
            return false;
        }

        const slong pivot = fmpz_mat_row_first_nonzero_col(
                flint::FmpzMatConstRef(current), 0);
        bool new_row = false;
        if (pivot < input_cols) {
            const slong position = dense_hnf_find_row_starting_with(
                    flint::FmpzMatConstRef(basis), active_rows, pivot);
            if (!dense_hnf_insert_row(
                        flint::FmpzMatRef(basis),
                        flint::FmpzMatRef(coeffs), active_rows, position,
                        flint::FmpzMatConstRef(current),
                        flint::FmpzMatConstRef(current_coeff))) {
                return false;
            }
            changed_pivots.push_back(pivot);
            new_row = true;
        }

        if (new_row && changed_pivots.size() == 1) {
            continue;
        }
        if (!dense_hnf_reduce_up(
                    flint::FmpzMatRef(basis),
                    flint::FmpzMatRef(coeffs), active_rows,
                    changed_pivots)) {
            return false;
        }
    }

    flint::FmpzMat out_hnf(active_rows, input_cols);
    flint::FmpzMat out_transform(active_rows, input_rows);
    for (slong row = 0; row < active_rows; ++row) {
        if (!fmpz_mat_copy_row_between(
                    flint::FmpzMatRef(out_hnf), row,
                    flint::FmpzMatConstRef(basis), row) ||
            !fmpz_mat_copy_row_between(
                    flint::FmpzMatRef(out_transform), row,
                    flint::FmpzMatConstRef(coeffs), row)) {
            return false;
        }
    }
    hnf = std::move(out_hnf);
    transform = std::move(out_transform);
    return true;
}

bool fmpz_mat_row_divexact_si(flint::FmpzMat& out,
                              bool& divisible,
                              flint::FmpzMatConstRef row,
                              slong divisor) noexcept {
    divisible = false;
    if (divisor <= 0 || flint::fmpz_mat_nrows(row) != 1) {
        return false;
    }

    flint::Fmpz divisor_fmpz;
    flint::fmpz_set_si(flint::FmpzRef(divisor_fmpz), divisor);
    flint::FmpzMat divided(1, flint::fmpz_mat_ncols(row));
    for (slong j = 0; j < flint::fmpz_mat_ncols(row); ++j) {
        flint::FmpzConstRef entry = flint::fmpz_mat_entry(row, 0, j);
        if (!flint::fmpz_divisible(entry,
                                   flint::FmpzConstRef(divisor_fmpz))) {
            return true;
        }
        flint::fmpz_divexact(
                flint::fmpz_mat_entry(divided, 0, j), entry,
                flint::FmpzConstRef(divisor_fmpz));
    }

    out = std::move(divided);
    divisible = true;
    return true;
}

bool saturation_reconstruct_candidate(
        FactoredElement& element,
        flint::FmpzMatRef factor_base_row,
        const std::vector<FactoredElement>& relations,
        flint::FmpzMatConstRef relation_rows,
        flint::FmpzMatConstRef candidates,
        slong column,
        const NumberField& field) noexcept {
    const slong relation_count = static_cast<slong>(relations.size());
    const slong candidate_rows = flint::fmpz_mat_nrows(candidates);
    const slong candidate_cols = flint::fmpz_mat_ncols(candidates);
    const slong generators = flint::fmpz_mat_ncols(relation_rows);
    if (!element.is_defined() ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(field) ||
        flint::fmpz_mat_nrows(factor_base_row) != 1 ||
        flint::fmpz_mat_ncols(factor_base_row) != generators ||
        flint::fmpz_mat_nrows(relation_rows) != relation_count ||
        column < 0 || column >= candidate_cols ||
        !(candidate_rows == relation_count ||
          candidate_rows == relation_count + 1) ||
        !element.one()) {
        return false;
    }

    flint::fmpz_mat_zero(factor_base_row);
    for (slong j = 0; j < relation_count; ++j) {
        flint::FmpzConstRef coefficient =
                flint::fmpz_mat_entry(candidates, j, column);
        if (flint::fmpz_is_zero(coefficient)) {
            continue;
        }
        if (!compact_multiply_power_fmpz(
                    element, relations[static_cast<std::size_t>(j)],
                    coefficient)) {
            return false;
        }
        for (slong k = 0; k < generators; ++k) {
            flint::fmpz_addmul(
                    flint::fmpz_mat_entry(factor_base_row, 0, k),
                    coefficient,
                    flint::fmpz_mat_entry(relation_rows, j, k));
        }
    }

    if (candidate_rows == relation_count + 1) {
        flint::FmpzConstRef torsion_coefficient =
                flint::fmpz_mat_entry(candidates, relation_count, column);
        if (!flint::fmpz_is_zero(torsion_coefficient)) {
            Element zeta(field);
            FactoredElement torsion(field);
            if (!zeta.is_defined() ||
                !root_of_unity_generator(zeta, field) ||
                !torsion.is_defined() || !torsion.set_element(zeta) ||
                !compact_multiply_power_fmpz(element, torsion,
                                             torsion_coefficient)) {
                return false;
            }
        }
    }

    element.normalize();
    return true;
}

bool relation_dlog_column(flint::FmpzMat& out,
                          const ClassGroupContext& context,
                          const PrimeIdeal& prime,
                          flint::FmpzConstRef ell) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_dlog_column");
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!context.has_factor_base() || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    const slong relation_len = context.relation_count();
    if (relation_len <= 0) {
        return false;
    }

    ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::FmpzMat candidate(relation_len, 1);
    ResidueFieldElement image(residue_field);
    ResidueFieldQuotientLog quotient_log(residue_field);
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.residue_quotient_log_setup");
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }
    for (slong i = 0; i < relation_len; ++i) {
        const Element* relation_generator = context.relation_generator_at(i);
        if (relation_generator == nullptr) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.relation_dlog_reduce_element");
            if (!image.set_element(*relation_generator)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.relation_dlog_apply");
            if (!quotient_log.apply(flint::fmpz_mat_entry(candidate, i, 0),
                                    image)) {
                return false;
            }
        }
    }

    out = std::move(candidate);
    return true;
}

bool relation_dlog_column_order_two(flint::FmpzMat& out,
                                    const ClassGroupContext& context,
                                    const PrimeIdeal& prime) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_dlog_column_order_two");
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!context.has_factor_base() || field == nullptr ||
        !same_order_parent(prime.parent(), order)) {
        return false;
    }

    const slong relation_len = context.relation_count();
    if (relation_len <= 0) {
        return false;
    }

    ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::Fmpz field_order;
    flint::Fmpz exponent;
    if (!residue_field.cardinality(flint::FmpzRef(field_order))) {
        return false;
    }
    flint::fmpz_sub_ui(flint::FmpzRef(field_order),
                       flint::FmpzConstRef(field_order), 1);
    if (::fmpz_divisible_ui(field_order.raw(), 2) == 0) {
        return false;
    }
    flint::fmpz_fdiv_q_2exp(flint::FmpzRef(exponent),
                            flint::FmpzConstRef(field_order), 1);
    flint::Fmpz characteristic;
    const bool degree_one_ui =
            residue_field.degree() == 1 &&
            residue_field.characteristic(flint::FmpzRef(characteristic)) &&
            flint::fmpz_abs_fits_ui(flint::FmpzConstRef(characteristic)) &&
            flint::fmpz_abs_fits_ui(flint::FmpzConstRef(exponent));
    const ulong characteristic_ui = degree_one_ui
            ? flint::fmpz_get_ui(flint::FmpzConstRef(characteristic))
            : UWORD(0);
    const ulong exponent_ui = degree_one_ui
            ? flint::fmpz_get_ui(flint::FmpzConstRef(exponent))
            : UWORD(0);
    const ulong characteristic_inverse =
            degree_one_ui ? n_preinvert_limb(characteristic_ui) : UWORD(0);

    flint::FmpzMat candidate(relation_len, 1);
    ResidueFieldElement image(residue_field);
    ResidueFieldElement powered(residue_field);
    ResidueFieldElement one(residue_field);
    ResidueFieldElement minus_one(residue_field);
    flint::Fmpz scalar;
    if (!image.is_defined() || !powered.is_defined() || !one.is_defined() ||
        !minus_one.is_defined() || !one.one() || !minus_one.negate(one)) {
        return false;
    }

    for (slong i = 0; i < relation_len; ++i) {
        const Element* relation_generator = context.relation_generator_at(i);
        if (relation_generator == nullptr) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_dlog_reduce_element_order_two");
            if (!image.set_element(*relation_generator)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_dlog_apply_order_two");
            if (degree_one_ui) {
                if (!image.degree_one_scalar(flint::FmpzRef(scalar)) ||
                    flint::fmpz_is_zero(flint::FmpzConstRef(scalar))) {
                    return false;
                }
                const ulong scalar_ui =
                        flint::fmpz_get_ui(flint::FmpzConstRef(scalar)) %
                        characteristic_ui;
                const ulong powered_ui = n_powmod2_ui_preinv(
                        scalar_ui, exponent_ui, characteristic_ui,
                        characteristic_inverse);
                if (powered_ui == UWORD(1) % characteristic_ui) {
                    flint::fmpz_set_ui(flint::fmpz_mat_entry(candidate, i, 0),
                                       0);
                } else if (powered_ui == characteristic_ui - UWORD(1)) {
                    flint::fmpz_set_ui(flint::fmpz_mat_entry(candidate, i, 0),
                                       1);
                } else {
                    return false;
                }
                continue;
            }
            if (!powered.pow_fmpz(image, flint::FmpzConstRef(exponent))) {
                return false;
            }
            if (powered.equal(one)) {
                flint::fmpz_set_ui(flint::fmpz_mat_entry(candidate, i, 0), 0);
            } else if (powered.equal(minus_one)) {
                flint::fmpz_set_ui(flint::fmpz_mat_entry(candidate, i, 0), 1);
            } else {
                return false;
            }
        }
    }

    out = std::move(candidate);
    return true;
}

bool residue_dlog_kernel(flint::FmpzMat& out,
                         const ClassGroupContext& context,
                         const PrimeIdeal& prime,
                         flint::FmpzConstRef ell) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.residue_dlog_kernel");
    if (!flint::fmpz_is_prime(ell) ||
        !same_order_parent(prime.parent(), context.parent())) {
        return false;
    }

    flint::FmpzMat column(0, 0);
    const bool order_two = ::fmpz_equal_ui(ell.raw(), 2) != 0;
    if (order_two) {
        if (!relation_dlog_column_order_two(column, context, prime)) {
            return false;
        }
    } else if (!relation_dlog_column(column, context, prime, ell)) {
        return false;
    }

    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.dlog_kernel_nullspace");
    return dlog_kernel_from_matrix(out, column, ell);
}

bool relation_dlog_column_with_units(flint::FmpzMat& out,
                                     const ClassGroupContext& context,
                                     const OrderUnitGroup& units,
                                     bool include_torsion,
                                     const PrimeIdeal& prime,
                                     flint::FmpzConstRef ell) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_dlog_column_with_units");
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!context.has_factor_base() || field == nullptr ||
        !same_order_parent(prime.parent(), order) || !units.is_set() ||
        !same_order_parent(units.parent(), order) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    bool torsion_needed = false;
    if (!relation_saturation_torsion_needed(torsion_needed, units, ell) ||
        include_torsion != torsion_needed) {
        return false;
    }

    const slong relation_len = context.relation_count();
    const slong free_rank = units.free_rank();
    const slong rows = relation_len + free_rank + (include_torsion ? 1 : 0);
    if (rows <= 0) {
        return false;
    }

    ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::FmpzMat candidate(rows, 1);
    ResidueFieldElement image(residue_field);
    ResidueFieldQuotientLog quotient_log(residue_field);
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.residue_quotient_log_setup_with_units");
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }
    for (slong i = 0; i < relation_len; ++i) {
        const Element* relation_generator = context.relation_generator_at(i);
        if (relation_generator == nullptr) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_dlog_reduce_element_with_units");
            if (!image.set_element(*relation_generator)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_dlog_apply_with_units");
            if (!quotient_log.apply(flint::fmpz_mat_entry(candidate, i, 0),
                                    image)) {
                return false;
            }
        }
    }

    for (slong i = 0; i < free_rank; ++i) {
        FactoredElement generator(*field);
        const slong row = relation_len + i;
        if (!units.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !quotient_log.apply(flint::fmpz_mat_entry(candidate, row, 0),
                                image)) {
            return false;
        }
    }

    if (include_torsion) {
        OrderElement torsion(*order);
        if (!units.torsion_generator(torsion) ||
            !image.set_order_element(torsion) ||
            !quotient_log.apply(flint::fmpz_mat_entry(
                                        candidate, relation_len + free_rank, 0),
                                image)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

bool residue_dlog_kernel_with_units(flint::FmpzMat& out,
                                    const ClassGroupContext& context,
                                    const OrderUnitGroup& units,
                                    bool include_torsion,
                                    const PrimeIdeal& prime,
                                    flint::FmpzConstRef ell) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.residue_dlog_kernel_with_units");
    if (!flint::fmpz_is_prime(ell) ||
        !same_order_parent(prime.parent(), context.parent())) {
        return false;
    }

    flint::FmpzMat column(0, 0);
    if (!relation_dlog_column_with_units(column, context, units,
                                         include_torsion, prime, ell)) {
        return false;
    }

    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.dlog_kernel_nullspace_with_units");
    return dlog_kernel_from_matrix(out, column, ell);
}

bool kernel_row_divisible(flint::FmpzMatConstRef kernel,
                          slong row,
                          flint::FmpzConstRef ell) noexcept {
    if (!flint::fmpz_is_prime(ell) || row < 0 ||
        row >= flint::fmpz_mat_nrows(kernel)) {
        return false;
    }

    for (slong j = 0; j < flint::fmpz_mat_ncols(kernel); ++j) {
        if (!flint::fmpz_divisible(flint::fmpz_mat_entry(kernel, row, j),
                                   ell)) {
            return false;
        }
    }
    return true;
}

bool saturation_relation_divided_row_from_rows(
        flint::FmpzMat& divided,
        bool& divisible,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef relation_rows,
        flint::FmpzMatConstRef kernel,
        slong row,
        flint::FmpzConstRef ell) noexcept {
    const slong relations = context.relation_count();
    const slong generators = context.generator_count();
    if (!context.has_factor_base() || !flint::fmpz_is_prime(ell) ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel) ||
        flint::fmpz_mat_ncols(kernel) < relations ||
        flint::fmpz_mat_nrows(relation_rows) != relations ||
        flint::fmpz_mat_ncols(relation_rows) != generators) {
        return false;
    }

    flint::FmpzMat candidate(1, generators);
    for (slong j = 0; j < relations; ++j) {
        flint::FmpzConstRef coefficient =
                flint::fmpz_mat_entry(kernel, row, j);
        if (flint::fmpz_is_zero(coefficient)) {
            continue;
        }
        for (slong i = 0; i < generators; ++i) {
            flint::fmpz_addmul(
                    flint::fmpz_mat_entry(candidate, 0, i), coefficient,
                    flint::fmpz_mat_entry(relation_rows, j, i));
        }
    }

    divisible = true;
    for (slong i = 0; i < generators; ++i) {
        if (!flint::fmpz_divisible(
                    flint::fmpz_mat_entry(flint::FmpzMatConstRef(candidate),
                                          0, i),
                    ell)) {
            divisible = false;
            return true;
        }
        ::fmpz_divexact(
                flint::fmpz_mat_entry(candidate, 0, i).raw(),
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(candidate), 0, i)
                        .raw(),
                ell.raw());
    }

    divided = std::move(candidate);
    return true;
}

bool saturation_relation_divided_row(flint::FmpzMat& divided,
                                     bool& divisible,
                                     const ClassGroupContext& context,
                                     flint::FmpzMatConstRef kernel,
                                     slong row,
                                     flint::FmpzConstRef ell) noexcept {
    const slong relations = context.relation_count();
    const slong generators = context.generator_count();
    flint::FmpzMat relation_rows(relations, generators);
    if (!context.relations(flint::FmpzMatRef(relation_rows))) {
        return false;
    }

    return saturation_relation_divided_row_from_rows(
            divided, divisible, context, flint::FmpzMatConstRef(relation_rows),
            kernel, row, ell);
}

bool saturation_relation_product_with_units(FactoredElement& product,
                                            const ClassGroupContext& context,
                                            flint::FmpzMatConstRef kernel,
                                            slong row,
                                            const OrderUnitGroup& units,
                                            bool include_torsion) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || product.parent() == nullptr ||
        !product.parent()->has_same_data(*field) || !product.one() ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel) ||
        flint::fmpz_mat_ncols(kernel) !=
                context.relation_count() + units.free_rank() +
                        (include_torsion ? 1 : 0)) {
        return false;
    }

    Element relation_generator(*field);
    for (slong j = 0; j < context.relation_count(); ++j) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(kernel, row, j);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!context.relation_generator(relation_generator, j) ||
            !multiply_element_power_fmpz(product, relation_generator,
                                         exponent)) {
            return false;
        }
    }

    const slong offset = context.relation_count();
    for (slong j = 0; j < units.free_rank(); ++j) {
        FactoredElement unit(*field);
        if (!units.free_generator(unit, j) ||
            !compact_multiply_power_fmpz(
                    product, unit, flint::fmpz_mat_entry(kernel, row,
                                                         offset + j))) {
            return false;
        }
    }

    if (include_torsion) {
        OrderElement torsion(*order);
        Element torsion_value(*field);
        FactoredElement torsion_factor(*field);
        if (!units.torsion_generator(torsion) ||
            !torsion.get_element(torsion_value) ||
            !torsion_factor.set_element(torsion_value) ||
            !compact_multiply_power_fmpz(
                    product, torsion_factor,
                    flint::fmpz_mat_entry(kernel, row,
                                          offset + units.free_rank()))) {
            return false;
        }
    }

    product.normalize();
    return true;
}

bool saturation_relation_product(FactoredElement& product,
                                 const ClassGroupContext& context,
                                 flint::FmpzMatConstRef kernel,
                                 slong row) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || product.parent() == nullptr ||
        !product.parent()->has_same_data(*field) || !product.one() ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel) ||
        flint::fmpz_mat_ncols(kernel) != context.relation_count()) {
        return false;
    }

    Element relation_generator(*field);
    for (slong j = 0; j < context.relation_count(); ++j) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(kernel, row, j);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!context.relation_generator(relation_generator, j) ||
            !multiply_element_power_fmpz(product, relation_generator,
                                         exponent)) {
            return false;
        }
    }

    product.normalize();
    return true;
}

bool saturation_relation(Relation& relation,
                         bool& is_relation,
                         const ClassGroupContext& context,
                         flint::FmpzMatConstRef kernel,
                         slong row,
                         flint::FmpzConstRef ell) noexcept {
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!same_factor_base(relation.factor_base(), base) ||
        field == nullptr || !flint::fmpz_is_prime(ell) ||
        !flint::fmpz_fits_si(ell) ||
        flint::fmpz_mat_ncols(kernel) != context.relation_count()) {
        return false;
    }

    flint::FmpzMat divided(1, base->length());
    bool divisible = false;
    if (!saturation_relation_divided_row(divided, divisible, context,
                                         kernel, row, ell)) {
        return false;
    }
    if (!divisible) {
        is_relation = false;
        return true;
    }

    FactoredElement product(*field);
    FactoredElement root(*field);
    if (!saturation_relation_product(product, context, kernel, row)) {
        return false;
    }

    bool is_power = false;
    if (!product.is_power_si(is_power, root, flint::fmpz_get_si(ell),
                             FactoredRootStrategy::reduced,
                             context.diagnostics())) {
        return false;
    }
    if (!is_power) {
        is_relation = false;
        return true;
    }

    Element beta(*field);
    Relation candidate(*base);
    flint::FmpzMat exponents(1, base->length());
    if (!root.evaluate(beta) ||
        !candidate.set_generator(beta, context.diagnostics()) ||
        !candidate.exponents(flint::FmpzMatRef(exponents))) {
        return false;
    }
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(exponents),
                               flint::FmpzMatConstRef(divided))) {
        is_relation = false;
        return true;
    }

    if (!relation.set(candidate)) {
        return false;
    }
    is_relation = true;
    return true;
}

bool saturation_relation_with_units(Relation& relation,
                                    bool& is_relation,
                                    const ClassGroupContext& context,
                                    flint::FmpzMatConstRef relation_rows,
                                    flint::FmpzMatConstRef kernel,
                                    slong row,
                                    const OrderUnitGroup& units,
                                    bool include_torsion,
                                    flint::FmpzConstRef ell) noexcept {
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!same_factor_base(relation.factor_base(), base) ||
        field == nullptr || !flint::fmpz_is_prime(ell) ||
        !flint::fmpz_fits_si(ell)) {
        return false;
    }

    flint::FmpzMat divided(1, base->length());
    bool divisible = false;
    {
        SILEX_PROFILE_SCOPE(
                context.diagnostics(), DiagnosticsModule::class_group,
                "class_group.saturation_relation_divided_row_with_units");
        if (!saturation_relation_divided_row_from_rows(
                    divided, divisible, context, relation_rows, kernel, row,
                    ell)) {
            return false;
        }
    }
    if (!divisible) {
        is_relation = false;
        return true;
    }

    FactoredElement product(*field);
    FactoredElement root(*field);
    {
        SILEX_PROFILE_SCOPE(
                context.diagnostics(), DiagnosticsModule::class_group,
                "class_group.saturation_relation_product_with_units");
        if (!saturation_relation_product_with_units(
                    product, context, kernel, row, units, include_torsion)) {
            return false;
        }
    }

    bool is_power = false;
    if (!product.is_power_si(is_power, root, flint::fmpz_get_si(ell),
                             FactoredRootStrategy::reduced,
                             context.diagnostics())) {
        return false;
    }
    if (!is_power) {
        is_relation = false;
        return true;
    }

    Element beta(*field);
    Relation candidate(*base);
    flint::FmpzMat exponents(1, base->length());
    if (!root.evaluate(beta) ||
        !candidate.set_generator(beta, context.diagnostics()) ||
        !candidate.exponents(flint::FmpzMatRef(exponents))) {
        return false;
    }
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(exponents),
                               flint::FmpzMatConstRef(divided))) {
        is_relation = false;
        return true;
    }

    if (!relation.set(candidate)) {
        return false;
    }
    is_relation = true;
    return true;
}

bool factor_base_covers_prime_decomposition(const FactorBase& base,
                                            flint::FmpzConstRef p,
                                            flint::FmpzConstRef bound)
        noexcept {
    const Order* order = base.parent();
    if (order == nullptr || !flint::fmpz_is_prime(p) ||
        flint::fmpz_sgn(bound) < 0) {
        return false;
    }

    PrimeIdealList primes;
    if (!decompose_prime(primes, *order, p)) {
        return false;
    }

    std::vector<slong> proof_targets;
    if (!detail::select_factor_base_proof_targets(
                proof_targets, primes, p, bound)) {
        return false;
    }

    for (slong index : proof_targets) {
        const PrimeIdeal* prime = primes.at(index);
        if (prime == nullptr || !base.contains(*prime)) {
            return false;
        }
    }

    return true;
}

bool relation_saturation_option_enabled(
        const detail::ClassGroupRelationOptions& options) noexcept {
    return options.relation_saturation_aux_prime_bound != 0 ||
           options.relation_saturation_max_appends_per_ell != 0;
}

struct IdealPowerEntry {
    explicit IdealPowerEntry(const Order& order) noexcept
        : ideal(order) {
    }

    IdealPowerEntry(const IdealPowerEntry&) = delete;
    IdealPowerEntry& operator=(const IdealPowerEntry&) = delete;
    IdealPowerEntry(IdealPowerEntry&&) noexcept = default;
    IdealPowerEntry& operator=(IdealPowerEntry&&) noexcept =
            default;

    Ideal ideal;
    PrimeIdeal prime;
    slong exponent = 0;
    bool has_prime = false;
};

bool checked_add_slong(slong& out, slong left, slong right) noexcept {
    if ((right > 0 && left > std::numeric_limits<slong>::max() - right) ||
        (right < 0 && left < std::numeric_limits<slong>::min() - right)) {
        return false;
    }

    out = left + right;
    return true;
}

bool checked_neg_slong(slong& out, slong value) noexcept {
    if (value == std::numeric_limits<slong>::min()) {
        return false;
    }

    out = -value;
    return true;
}

bool checked_mul_slong(slong& out, slong left, slong right) noexcept {
    const slong min = std::numeric_limits<slong>::min();
    const slong max = std::numeric_limits<slong>::max();
    if (left == 0 || right == 0) {
        out = 0;
        return true;
    }
    if ((left == -1 && right == min) ||
        (right == -1 && left == min)) {
        return false;
    }
    if (left > 0) {
        if (right > 0) {
            if (left > max / right) {
                return false;
            }
        } else if (right < min / left) {
            return false;
        }
    } else if (right > 0) {
        if (left < min / right) {
            return false;
        }
    } else if (left < max / right) {
        return false;
    }

    out = left * right;
    return true;
}

bool checked_pow_slong(slong& out, slong base, slong exponent) noexcept {
    if (base <= 0 || exponent < 0) {
        return false;
    }

    slong result = 1;
    for (slong i = 0; i < exponent; ++i) {
        if (!checked_mul_slong(result, result, base)) {
            return false;
        }
    }

    out = result;
    return true;
}

bool set_order_element_fmpz(OrderElement& out,
                            const Order& order,
                            flint::FmpzConstRef value) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr || !same_order_parent(out.parent(), &order)) {
        return false;
    }

    Element element(*field);
    return element.is_defined() && element.set_fmpz(value) &&
           out.set_element(element);
}

bool add_fmpz_to_order_element(OrderElement& out,
                               const OrderElement& input,
                               flint::FmpzConstRef value) noexcept {
    const Order* order = input.parent();
    if (order == nullptr || !same_order_parent(out.parent(), order)) {
        return false;
    }

    OrderElement scalar(*order);
    return scalar.is_defined() &&
           set_order_element_fmpz(scalar, *order, value) &&
           out.add(input, scalar);
}

bool order_element_norm_divisible_by(
        bool& out,
        const OrderElement& input,
        flint::FmpzConstRef divisor) noexcept {
    const Order* order = input.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || flint::fmpz_is_zero(divisor)) {
        return false;
    }

    Element element(*field);
    flint::Fmpq norm;
    flint::Fmpz numerator_abs;
    if (!element.is_defined() || !input.get_element(element) ||
        !element.norm(flint::FmpqRef(norm)) ||
        !flint::fmpz_is_one(flint::fmpq_den_ref(flint::FmpqConstRef(norm)))) {
        return false;
    }

    flint::fmpz_abs(flint::FmpzRef(numerator_abs),
                    flint::fmpq_num_ref(flint::FmpqConstRef(norm)));
    out = flint::fmpz_divisible(flint::FmpzConstRef(numerator_abs),
                                divisor);
    return true;
}

bool adjust_prime_two_generator(OrderElement& out,
                                      const OrderElement& input,
                                      const PrimeIdeal& prime) noexcept {
    const Order* order = input.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order) ||
        !same_order_parent(out.parent(), order)) {
        return false;
    }

    const slong ramification = prime.ramification_index();
    const slong residue_degree = prime.residue_degree();
    if (ramification <= 0 || residue_degree <= 0) {
        return false;
    }
    if (ramification != 1) {
        return out.set(input);
    }

    flint::Fmpz p;
    flint::Fmpz p_norm_times_p;
    bool divisible = false;
    if (!prime.rational_prime(flint::FmpzRef(p))) {
        return false;
    }

    flint::fmpz_pow_ui(flint::FmpzRef(p_norm_times_p),
                       flint::FmpzConstRef(p),
                       static_cast<ulong>(residue_degree + 1));
    if (!order_element_norm_divisible_by(
                divisible, input, flint::FmpzConstRef(p_norm_times_p))) {
        return false;
    }
    if (!divisible) {
        return out.set(input);
    }

    return add_fmpz_to_order_element(out, input, flint::FmpzConstRef(p));
}

bool order_element_positive_power_slong(OrderElement& out,
                                        const OrderElement& input,
                                        slong exponent) noexcept {
    const Order* order = input.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !same_order_parent(out.parent(), order) ||
        exponent < 0) {
        return false;
    }
    if (exponent == 0) {
        return out.one();
    }

    Element base(*field);
    Element power(*field);
    flint::Fmpz exponent_fmpz;
    if (!base.is_defined() || !power.is_defined() ||
        !input.get_element(base)) {
        return false;
    }
    flint::fmpz_set_si(flint::FmpzRef(exponent_fmpz), exponent);
    return power.pow_fmpz(base, flint::FmpzConstRef(exponent_fmpz)) &&
           out.set_element(power);
}

bool prime_ideal_power_positive_slong(
        Ideal& out,
        bool& handled,
        const PrimeIdeal& prime,
        slong exponent,
        const DiagnosticsContext* diagnostics = nullptr) noexcept {
    handled = false;
    const Order* order = out.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order) ||
        exponent < 0 || !prime.has_prime_data()) {
        return false;
    }
    if (exponent == 0) {
        handled = true;
        return out.one();
    }
    if (exponent == 1) {
        handled = true;
        return prime.get_ideal(out);
    }

    flint::FmpzMat generator_coordinates(1, order->degree());
    if (!prime.kummer_generator_coordinates(
                flint::FmpzMatRef(generator_coordinates))) {
        return true;
    }

    const slong ramification = prime.ramification_index();
    if (ramification <= 0) {
        return false;
    }

    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.prime_ideal_power");

    flint::Fmpz p;
    flint::Fmpz first_generator_integer;
    OrderElement raw_second_generator(*order);
    OrderElement second_generator(*order);
    OrderElement second_generator_power(*order);
    if (!prime.rational_prime(flint::FmpzRef(p)) ||
        !raw_second_generator.is_defined() ||
        !second_generator.is_defined() ||
        !second_generator_power.is_defined() ||
        !raw_second_generator.set_coordinates(
                flint::FmpzMatConstRef(generator_coordinates)) ||
        !adjust_prime_two_generator(
                second_generator, raw_second_generator, prime)) {
        return false;
    }

    const slong minimum_exponent = ((exponent - 1) / ramification) + 1;
    flint::fmpz_pow_ui(flint::FmpzRef(first_generator_integer),
                       flint::FmpzConstRef(p),
                       static_cast<ulong>(minimum_exponent));
    if (!order_element_positive_power_slong(
                second_generator_power, second_generator, exponent) ||
        !detail::set_known_two_generator_ideal(
                out, flint::FmpzConstRef(first_generator_integer),
                second_generator_power)) {
        return false;
    }

    handled = true;
    return true;
}

bool checked_abs_plus_one_slong(slong& out, slong value) noexcept {
    if (value == std::numeric_limits<slong>::min()) {
        return false;
    }

    const slong absolute = value < 0 ? -value : value;
    if (absolute == std::numeric_limits<slong>::max()) {
        return false;
    }

    out = absolute + 1;
    return true;
}

slong floor_div_positive_slong(slong value, slong divisor) noexcept {
    slong quotient = value / divisor;
    const slong remainder = value % divisor;
    if (remainder < 0) {
        --quotient;
    }
    return quotient;
}

bool floor_log_slong(slong& out, slong value, slong base) noexcept {
    if (value <= 0 || base <= 1) {
        return false;
    }

    slong current = 1;
    slong exponent = 0;
    while (current <= value / base) {
        current *= base;
        ++exponent;
    }

    out = exponent;
    return true;
}

bool ceil_log_slong(slong& out, slong value, slong base) noexcept {
    slong floor_value = 0;
    slong floor_power = 0;
    if (!floor_log_slong(floor_value, value, base) ||
        !checked_pow_slong(floor_power, base, floor_value)) {
        return false;
    }

    out = floor_value + (floor_power < value ? 1 : 0);
    return true;
}

bool decomposition_row_max_abs_plus_one(
        slong& out,
        flint::FmpzMatConstRef row) noexcept {
    if (flint::fmpz_mat_nrows(row) != 1) {
        return false;
    }

    slong maximum = 1;
    for (slong j = 0; j < flint::fmpz_mat_ncols(row); ++j) {
        flint::FmpzConstRef entry = flint::fmpz_mat_entry(row, 0, j);
        if (!flint::fmpz_fits_si(entry)) {
            return false;
        }
        slong candidate = 0;
        if (!checked_abs_plus_one_slong(candidate,
                                        flint::fmpz_get_si(entry))) {
            return false;
        }
        maximum = max_slong_value(maximum, candidate);
    }

    out = maximum;
    return true;
}

bool add_ideal_power(std::vector<IdealPowerEntry>& entries,
                     const Ideal& ideal,
                     slong exponent,
                     const PrimeIdeal* prime = nullptr) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf()) {
        return false;
    }
    if (prime != nullptr &&
        (!prime->has_prime_data() ||
         !same_order_parent(prime->parent(), order))) {
        return false;
    }
    if (exponent == 0) {
        return true;
    }

    for (IdealPowerEntry& entry : entries) {
        if (entry.ideal.equal(ideal)) {
            slong sum = 0;
            if (!checked_add_slong(sum, entry.exponent, exponent)) {
                return false;
            }
            entry.exponent = sum;
            if (prime != nullptr && !entry.has_prime) {
                if (!entry.prime.set(*prime)) {
                    return false;
                }
                entry.has_prime = true;
            }
            return true;
        }
    }

    entries.emplace_back(*order);
    IdealPowerEntry& entry = entries.back();
    if (!entry.ideal.is_defined() || !entry.ideal.set(ideal)) {
        entries.pop_back();
        return false;
    }
    if (prime != nullptr) {
        if (!entry.prime.set(*prime)) {
            entries.pop_back();
            return false;
        }
        entry.has_prime = true;
    }
    entry.exponent = exponent;
    return true;
}

bool factored_multiply_in_place(FactoredElement& target,
                                const FactoredElement& factor) noexcept {
    const NumberField* field = target.parent();
    if (field == nullptr || factor.parent() == nullptr ||
        !field->has_same_data(*factor.parent())) {
        return false;
    }

    FactoredElement product(*field);
    if (!product.multiply(target, factor)) {
        return false;
    }

    target.swap(product);
    return true;
}

bool factored_multiply_power_si(FactoredElement& target,
                                const FactoredElement& factor,
                                slong exponent) noexcept {
    const NumberField* field = target.parent();
    if (field == nullptr || factor.parent() == nullptr ||
        !field->has_same_data(*factor.parent())) {
        return false;
    }
    if (exponent == 0) {
        return true;
    }

    FactoredElement power(*field);
    return power.is_defined() && power.pow_si(factor, exponent) &&
           factored_multiply_in_place(target, power);
}

bool abs_order_discriminant(flint::Fmpz& out, const Order& order) noexcept {
    flint::Fmpz discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_is_zero(flint::FmpzConstRef(discriminant))) {
        return false;
    }

    flint::fmpz_abs(flint::FmpzRef(out),
                    flint::FmpzConstRef(discriminant));
    return !flint::fmpz_is_zero(flint::FmpzConstRef(out));
}

bool ideal_norm_product_gt(bool& out,
                           const Ideal& left,
                           const Ideal& right,
                           flint::FmpzConstRef bound) noexcept {
    flint::Fmpz left_norm;
    flint::Fmpz right_norm;
    flint::Fmpz product;
    if (!left.norm(flint::FmpzRef(left_norm)) ||
        !right.norm(flint::FmpzRef(right_norm))) {
        return false;
    }

    flint::fmpz_mul(flint::FmpzRef(product),
                    flint::FmpzConstRef(left_norm),
                    flint::FmpzConstRef(right_norm));
    out = flint::fmpz_cmp(flint::FmpzConstRef(product), bound) > 0;
    return true;
}

bool reduce_ideal_powers(
        Ideal& reduced,
        FactoredElement& multiplier,
        const std::vector<IdealPowerEntry>& powers,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = reduced.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) || precision <= 0) {
        return false;
    }

    flint::Fmpz abs_discriminant;
    if (!abs_order_discriminant(abs_discriminant, *order)) {
        return false;
    }

    Ideal current(*order);
    FactoredElement current_multiplier(*field);
    if (!current.is_defined() || !current_multiplier.is_defined() ||
        !current_multiplier.one()) {
        return false;
    }

    bool have_current = false;
    for (const IdealPowerEntry& entry : powers) {
        if (entry.exponent == 0) {
            continue;
        }
        if (entry.ideal.parent() == nullptr ||
            !entry.ideal.parent()->has_same_data(*order)) {
            return false;
        }

        Ideal part(*order);
        FactoredElement part_multiplier(*field);
        if (!part.is_defined() || !part_multiplier.is_defined() ||
            !detail::reduce_ideal_signed_power(
                    part, part_multiplier, entry.ideal, entry.exponent,
                    precision, diagnostics)) {
            return false;
        }

        if (!have_current) {
            if (!current.set(part) ||
                !current_multiplier.set(part_multiplier)) {
                return false;
            }
            have_current = true;
            continue;
        }

        if (!factored_multiply_in_place(current_multiplier,
                                        part_multiplier)) {
            return false;
        }

        bool reduce_product = false;
        if (!ideal_norm_product_gt(
                    reduce_product, current, part,
                    flint::FmpzConstRef(abs_discriminant))) {
            return false;
        }
        if (reduce_product) {
            Ideal next(*order);
            Element reduction_multiplier(*field);
            if (!next.is_defined() || !reduction_multiplier.is_defined() ||
                !detail::reduce_ideal_product(
                        next, reduction_multiplier, current, part,
                        precision, diagnostics) ||
                !current.set(next) ||
                !current_multiplier.push(reduction_multiplier, -1)) {
                return false;
            }
        } else {
            Ideal product(*order);
            if (!product.is_defined() ||
                !product.multiply(current, part) ||
                !current.set(product)) {
                return false;
            }
        }
    }

    if (!have_current && !current.one()) {
        return false;
    }

    current_multiplier.normalize();
    return reduced.set(current) && multiplier.set(current_multiplier);
}

bool ideal_pow_positive_slong(Ideal& out,
                              const Ideal& input,
                              slong exponent,
                              const DiagnosticsContext* diagnostics =
                                      nullptr) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.ideal_pow_positive_slong");
    const Order* order = input.parent();
    if (order == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*order) || exponent < 0 ||
        !input.has_hnf()) {
        return false;
    }
    if (exponent == 0) {
        return out.one();
    }
    if (exponent == 1) {
        return out.set(input);
    }
    if (exponent == 2) {
        return out.multiply(input, input);
    }

    Ideal result(*order);
    Ideal base(*order);
    Ideal product(*order);
    Ideal square(*order);
    if (!result.is_defined() || !base.is_defined() ||
        !product.is_defined() || !square.is_defined() || !result.one() ||
        !base.set(input)) {
        return false;
    }

    slong active_exponent = exponent;
    while (active_exponent > 0) {
        if ((active_exponent & 1) != 0) {
            if (!product.multiply(result, base)) {
                return false;
            }
            result.swap(product);
        }
        active_exponent >>= 1;
        if (active_exponent > 0) {
            if (!square.multiply(base, base)) {
                return false;
            }
            base.swap(square);
        }
    }

    return out.set(result);
}

bool compact_de_evaluate(
        Ideal& out,
        const std::vector<IdealPowerEntry>& decomposition,
        const DiagnosticsContext* diagnostics = nullptr) noexcept {
    const Order* order = out.parent();
    if (order == nullptr || !out.is_defined()) {
        return false;
    }

    const IdealPowerEntry* single_nonzero = nullptr;
    slong nonzero_entries = 0;
    for (const IdealPowerEntry& entry : decomposition) {
        if (entry.exponent == 0) {
            continue;
        }
        if (entry.exponent < 0 || entry.ideal.parent() == nullptr ||
            !entry.ideal.parent()->has_same_data(*order)) {
            return false;
        }
        ++nonzero_entries;
        if (nonzero_entries == 1) {
            single_nonzero = &entry;
        }
    }

    if (nonzero_entries == 0) {
        return out.one();
    }

    if (nonzero_entries == 1 && single_nonzero != nullptr) {
        if (single_nonzero->exponent == 1) {
            return out.set(single_nonzero->ideal);
        }

        Ideal power(*order);
        if (!power.is_defined()) {
            return false;
        }
        bool handled = false;
        if (single_nonzero->has_prime) {
            if (!prime_ideal_power_positive_slong(
                        power, handled, single_nonzero->prime,
                        single_nonzero->exponent, diagnostics)) {
                return false;
            }
            if (handled) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_final_evaluate_prime_power");
            }
        }
        if (!handled &&
            !ideal_pow_positive_slong(power, single_nonzero->ideal,
                                      single_nonzero->exponent,
                                      diagnostics)) {
            return false;
        }
        return out.set(power);
    }

    Ideal accumulator(*order);
    Ideal power(*order);
    Ideal product(*order);
    if (!accumulator.is_defined() || !power.is_defined() ||
        !product.is_defined() || !accumulator.one()) {
        return false;
    }

    for (const IdealPowerEntry& entry : decomposition) {
        if (entry.exponent == 0) {
            continue;
        }
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_final_evaluate.entry");
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_final_evaluate_power");
            bool handled = false;
            if (entry.has_prime) {
                if (!prime_ideal_power_positive_slong(
                            power, handled, entry.prime, entry.exponent,
                            diagnostics)) {
                    return false;
                }
                if (handled) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_final_evaluate_prime_power");
                }
            }
            if (!handled) {
                if (!ideal_pow_positive_slong(power, entry.ideal,
                                              entry.exponent, diagnostics)) {
                    return false;
                }
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_final_evaluate_multiply");
            if (!product.multiply(accumulator, power)) {
                return false;
            }
        }
        accumulator.swap(product);
    }

    return out.set(accumulator);
}

bool compact_de_insert_factorization(
        std::vector<IdealPowerEntry>& decomposition,
        const Ideal& ideal,
        slong exponent_multiplier) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || exponent_multiplier < 0) {
        return false;
    }
    if (exponent_multiplier == 0 || ideal.is_one()) {
        return true;
    }

    IdealFactorization factorization(*order);
    PrimeIdeal prime(*order);
    Ideal prime_ideal(*order);
    if (!factorization.is_defined() || !prime.is_defined() ||
        !prime_ideal.is_defined() || !factorization.factor(ideal)) {
        return false;
    }

    for (slong i = 0; i < factorization.length(); ++i) {
        slong factor_exponent = 0;
        slong exponent = 0;
        if (!factorization.prime(prime, i) ||
            !factorization.exponent(factor_exponent, i) ||
            !checked_mul_slong(exponent, factor_exponent,
                               exponent_multiplier) ||
            !prime.get_ideal(prime_ideal) ||
            !add_ideal_power(decomposition, prime_ideal, exponent, &prime)) {
            return false;
        }
    }
    return true;
}

bool compact_de_select(
        std::vector<IdealPowerEntry>& selected,
        const std::vector<IdealPowerEntry>& decomposition,
        slong power) noexcept {
    if (power <= 0) {
        return false;
    }

    selected.clear();
    for (const IdealPowerEntry& entry : decomposition) {
        if (entry.exponent >= power) {
            const slong quotient = entry.exponent / power;
            if (quotient > 0 &&
                !add_ideal_power(selected, entry.ideal, quotient,
                                 entry.has_prime ? &entry.prime : nullptr)) {
                return false;
            }
        } else if (entry.exponent < 0) {
            return false;
        }
    }
    return true;
}

bool compact_de_subtract_selected(
        std::vector<IdealPowerEntry>& decomposition,
        const std::vector<IdealPowerEntry>& selected,
        slong power) noexcept {
    if (power <= 0) {
        return false;
    }

    for (const IdealPowerEntry& entry : selected) {
        slong delta = 0;
        slong negative_delta = 0;
        if (!checked_mul_slong(delta, power, entry.exponent) ||
            !checked_neg_slong(negative_delta, delta) ||
            !add_ideal_power(decomposition, entry.ideal, negative_delta,
                             entry.has_prime ? &entry.prime : nullptr)) {
            return false;
        }
    }
    return true;
}

bool compact_factored_logs(
        flint::ArbVec& out,
        const FactoredElement& element,
        const FactoredElement& multiplier,
        EmbeddingContext& embeddings,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_factored_logs");
    const NumberField* field = element.parent();
    if (field == nullptr || multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) || precision <= 0) {
        return false;
    }

    Signature sig;
    if (!signature(sig, *field)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    flint::ArbVec element_logs(places);
    flint::ArbVec multiplier_logs(places);
    flint::ArbVec sum(places);
    if (!embeddings.is_defined() || !embeddings.refine(precision)) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_factored_logs_element");
        if (!element.logarithmic_embedding(
                    flint::ArbVecRef(element_logs), embeddings,
                    LogEmbeddingMode::product, precision, diagnostics)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_factored_logs_multiplier");
        if (!multiplier.logarithmic_embedding(
                    flint::ArbVecRef(multiplier_logs), embeddings,
                    LogEmbeddingMode::product, precision, diagnostics)) {
            return false;
        }
    }

    for (slong i = 0; i < places; ++i) {
        ::arb_add(sum.data() + i, element_logs.data() + i,
                  multiplier_logs.data() + i, precision);
    }
    out.swap(sum);
    return true;
}

bool compact_factored_log_abs_upper_bound(
        flint::Fmpz& out,
        const Order& order,
        const FactoredElement& element,
        slong precision,
        const DiagnosticsContext* diagnostics,
        DiagnosticsModule module,
        const char* profile_label) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, module, profile_label);
    (void)module;
    (void)profile_label;
    const NumberField* field = order.parent();
    if (field == nullptr || element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) || precision <= 0) {
        return false;
    }

    Signature sig;
    if (!signature(sig, *field)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    flint::ArbVec logs(places);
    EmbeddingContext embeddings(*field);
    if (!embeddings.is_defined() || !embeddings.refine(precision) ||
        !element.logarithmic_embedding(
                flint::ArbVecRef(logs), embeddings,
                LogEmbeddingMode::product, precision, diagnostics)) {
        return false;
    }

    flint::Fmpz maximum;
    flint::Fmpz bound;
    flint::Arf upper;
    flint::fmpz_zero(flint::FmpzRef(maximum));
    for (slong i = 0; i < places; ++i) {
        flint::arb_get_abs_ubound_arf(upper, logs.data() + i, precision);
        if (!flint::arf_is_finite(upper)) {
            return false;
        }
        flint::arf_get_fmpz(bound, upper, ARF_RND_CEIL);
        if (flint::fmpz_sgn(flint::FmpzConstRef(bound)) < 0) {
            return false;
        }
        if (flint::fmpz_cmp(flint::FmpzConstRef(bound),
                            flint::FmpzConstRef(maximum)) > 0) {
            flint::fmpz_set(flint::FmpzRef(maximum),
                            flint::FmpzConstRef(bound));
        }
    }

    flint::fmpz_set(flint::FmpzRef(out), flint::FmpzConstRef(maximum));
    return true;
}

bool compact_use_residual_power_test(
        bool& use_residual,
        const Order& order,
        const FactoredElement& element,
        const FactoredElement& leftover,
        const DiagnosticsContext* diagnostics,
        DiagnosticsModule module = DiagnosticsModule::unit_group,
        const char* split_profile_label =
                "unit_group.compact_unit_power_split_bound",
        const char* original_profile_label =
                "unit_group.compact_unit_original_log_bound",
        const char* residual_profile_label =
                "unit_group.compact_unit_residual_log_bound")
        noexcept {
    constexpr slong kPowerSplitLogPrecision = 64;
    SILEX_PROFILE_SCOPE(diagnostics, module, split_profile_label);
    (void)split_profile_label;
    flint::Fmpz original_bound;
    flint::Fmpz residual_bound;
    flint::Fmpz sqrt_original;
    if (!compact_factored_log_abs_upper_bound(
                original_bound, order, element, kPowerSplitLogPrecision,
                diagnostics, module, original_profile_label) ||
        !compact_factored_log_abs_upper_bound(
                residual_bound, order, leftover, kPowerSplitLogPrecision,
                diagnostics, module, residual_profile_label)) {
        return false;
    }

    flint::fmpz_sqrt(flint::FmpzRef(sqrt_original),
                     flint::FmpzConstRef(original_bound));
    use_residual = flint::fmpz_cmp(flint::FmpzConstRef(residual_bound),
                                   flint::FmpzConstRef(sqrt_original)) <= 0;
    return true;
}

bool compact_arb_radius_lt_2exp(const arb_struct* value,
                                      slong exponent,
                                      slong precision) noexcept {
    if (value == nullptr || precision <= 0 || ::arb_is_finite(value) == 0) {
        return false;
    }

    flint::Arb radius;
    flint::Arf upper;
    ::arb_get_rad_arb(radius.raw(), value);
    ::arb_get_ubound_arf(upper.raw(), radius.raw(), precision);
    return flint::arf_is_finite(upper) &&
           ::arf_cmpabs_2exp_si(upper.raw(), exponent) < 0;
}

bool compact_round_arb_to_fmpz(flint::FmpzRef out,
                                     const flint::Arb& value) noexcept {
    if (!flint::arb_is_finite(value)) {
        return false;
    }

    flint::Arf midpoint;
    flint::arf_set(midpoint, arb_midref(value.raw()));
    if (!flint::arf_is_finite(midpoint)) {
        return false;
    }
    flint::arf_get_fmpz(out, midpoint, ARF_RND_NEAR);
    return true;
}

bool compact_weight_vector(
        flint::FmpzMat& weights,
        const flint::ArbVec& logs,
        const Signature& signature,
        slong power,
        slong precision) noexcept {
    const slong r1 = signature.r1();
    const slong r2 = signature.r2();
    const slong places = r1 + r2;
    const slong degree = r1 + 2 * r2;
    if (power <= 0 || precision <= 0 || logs.length() != places ||
        flint::fmpz_mat_nrows(weights) != 1 ||
        flint::fmpz_mat_ncols(weights) != degree) {
        return false;
    }

    flint::Arb log_two;
    flint::Arb scaled_place;
    flint::Arb weight_value;
    flint::arb_log_ui(log_two, 2, precision);
    if (!flint::arb_is_finite(log_two) ||
        !flint::arb_is_positive(log_two)) {
        return false;
    }

    slong column = 0;
    for (slong place = 0; place < places; ++place) {
        ::arb_div_ui(scaled_place.raw(), logs.data() + place,
                     static_cast<ulong>(power), precision);
        if (!compact_arb_radius_lt_2exp(scaled_place.raw(), -5,
                                             precision)) {
            return false;
        }
        flint::arb_div(weight_value, scaled_place,
                       flint::ArbConstRef(log_two), precision);
        if (place >= r1) {
            flint::arb_div_ui(weight_value, weight_value, 2, precision);
        }
        if (!compact_round_arb_to_fmpz(
                    flint::fmpz_mat_entry(weights, 0, column),
                    weight_value)) {
            return false;
        }
        ++column;
        if (place >= r1) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(weights, 0, column),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(weights), 0,
                            column - 1));
            ++column;
        }
    }
    return column == degree;
}

bool compact_log_loop_bound(slong& out,
                                  const flint::ArbVec& logs,
                                  slong n,
                                  slong precision) noexcept {
    if (n <= 1 || precision <= 0) {
        return false;
    }

    flint::Arb maximum;
    flint::Arb absolute;
    flint::Arb one_plus;
    flint::Arb numerator;
    flint::Arb denominator;
    flint::Arb ratio;
    flint::arb_zero(maximum);
    for (slong i = 0; i < logs.length(); ++i) {
        ::arb_abs(absolute.raw(), logs.data() + i);
        flint::arb_max(maximum, maximum, absolute, precision);
    }

    flint::arb_add_ui(one_plus, maximum, 1, precision);
    flint::arb_log(numerator, one_plus, precision);
    flint::arb_log_ui(denominator, static_cast<ulong>(n), precision);
    flint::arb_div(ratio, numerator, denominator, precision);
    if (!flint::arb_is_finite(ratio)) {
        return false;
    }

    flint::Arf upper;
    flint::Fmpz bound;
    flint::arb_get_ubound_arf(upper, ratio, precision);
    if (!flint::arf_is_finite(upper)) {
        return false;
    }
    flint::arf_get_fmpz(bound, upper, ARF_RND_CEIL);
    if (flint::fmpz_sgn(flint::FmpzConstRef(bound)) < 0) {
        out = 0;
        return true;
    }
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(bound))) {
        return false;
    }
    out = flint::fmpz_get_si(flint::FmpzConstRef(bound));
    return true;
}

bool compact_short_norm_guard(
        const Element& short_element,
        const FractionalIdeal& inverse_ideal,
        flint::FmpzMatConstRef weights,
        const Order& order) noexcept {
    if (flint::fmpz_mat_nrows(weights) != 1 ||
        flint::fmpz_mat_ncols(weights) != order.degree()) {
        return false;
    }

    flint::Fmpq short_norm;
    flint::Fmpq inverse_norm;
    flint::Fmpq quotient;
    if (!short_element.norm(flint::FmpqRef(short_norm)) ||
        !inverse_ideal.norm(flint::FmpqRef(inverse_norm))) {
        return false;
    }
    flint::fmpq_div(quotient, short_norm, inverse_norm);
    flint::fmpq_abs(flint::FmpqRef(quotient),
                    flint::FmpqConstRef(quotient));

    flint::Fmpz sum_weights;
    for (slong column = 0; column < order.degree(); ++column) {
        flint::fmpz_add(flint::FmpzRef(sum_weights),
                        flint::FmpzConstRef(sum_weights),
                        flint::fmpz_mat_entry(weights, 0, column));
    }
    flint::fmpz_abs(flint::FmpzRef(sum_weights),
                    flint::FmpzConstRef(sum_weights));
    if (!flint::fmpz_abs_fits_ui(flint::FmpzConstRef(sum_weights))) {
        return false;
    }

    flint::Fmpz bound;
    if (!abs_order_discriminant(bound, order)) {
        return false;
    }
    const ulong shift =
            flint::fmpz_get_ui(flint::FmpzConstRef(sum_weights)) +
            static_cast<ulong>(order.degree());
    flint::fmpz_mul_2exp(flint::FmpzRef(bound),
                         flint::FmpzConstRef(bound), shift);

    flint::Fmpz rhs;
    flint::fmpz_mul(flint::FmpzRef(rhs),
                    flint::FmpzConstRef(bound),
                    flint::fmpq_den_ref(quotient));
    return flint::fmpz_cmp(flint::fmpq_num_ref(quotient),
                           flint::FmpzConstRef(rhs)) <= 0;
}

bool compact_multiply_integral_by_principal(
        Ideal& out,
        const Ideal& integral,
        const Element& multiplier,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.compact_multiply_integral_by_principal");
    const Order* order = integral.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*order) ||
        !multiplier.has_parent(*field)) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_multiply_principal_product");
        return detail::multiply_integral_ideal_by_element(
                out, integral, multiplier, diagnostics);
    }
}

bool compact_infinite_reduction(
        Ideal& residual_ideal,
        FactoredElement& multiplier,
        const Order& order,
        const FactoredElement& element,
        FactoredElement& accumulated_multiplier,
        std::vector<IdealPowerEntry>& decomposition,
        slong finite_bound,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_infinite_reduction");
    const NumberField* field = order.parent();
    if (field == nullptr || residual_ideal.parent() == nullptr ||
        !residual_ideal.parent()->has_same_data(order) ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) ||
        accumulated_multiplier.parent() == nullptr ||
        !accumulated_multiplier.parent()->has_same_data(*field) ||
        finite_bound < 0 || n <= 1 || arb_precision <= 0 ||
        short_precision <= 0) {
        return false;
    }

    Signature sig;
    if (!signature(sig, *field) || sig.degree() != order.degree()) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    EmbeddingContext embeddings(*field);
    flint::ArbVec logs(places);
    if (!compact_factored_logs(logs, element,
                                     accumulated_multiplier, embeddings,
                                     arb_precision, diagnostics)) {
        return false;
    }

    slong log_bound = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_log_loop_bound");
        while (!compact_log_loop_bound(log_bound, logs, n,
                                             arb_precision)) {
            if (arb_precision > kCompactMaxPrecision / 2) {
                return false;
            }
            arb_precision *= 2;
            if (!compact_factored_logs(logs, element,
                                             accumulated_multiplier,
                                             embeddings, arb_precision,
                                             diagnostics)) {
                return false;
            }
        }
    }

    const slong start_k = max_slong_value(finite_bound, log_bound);
    std::vector<slong> powers(static_cast<std::size_t>(start_k + 1));
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_power_table");
        for (slong k = 0; k <= start_k; ++k) {
            if (!checked_pow_slong(powers[static_cast<std::size_t>(k)], n, k)) {
                return false;
            }
        }
    }

    flint::Fmpz abs_discriminant;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_abs_discriminant");
        if (!abs_order_discriminant(abs_discriminant, order)) {
            return false;
        }
    }

    detail::OrderMinkowskiEmbeddingCache weighted_t2_cache;
    for (slong k = start_k; k >= 1; --k) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_infinite_step");
        const slong power = powers[static_cast<std::size_t>(k)];

        std::vector<IdealPowerEntry> selected;
        Ideal eA(order);
        FractionalIdeal eA_fractional(order);
        FractionalIdeal inverse_eA(order);
        Element short_element(*field);
        flint::FmpzMat weights(1, order.degree());
        if (!eA.is_defined() || !eA_fractional.is_defined() ||
            !inverse_eA.is_defined() || !short_element.is_defined()) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_select");
            if (!compact_de_select(selected, decomposition, power)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_evaluate");
            if (!compact_de_evaluate(eA, selected, diagnostics)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_step_integral_fractional");
            if (!detail::set_integral_ideal_known_hnf(eA_fractional, eA)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_invert");
            if (!inverse_eA.invert(eA_fractional)) {
                return false;
            }
        }

        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_weight_vector");
            while (!compact_weight_vector(weights, logs, sig, power,
                                                arb_precision)) {
                if (arb_precision > kCompactMaxPrecision / 2) {
                    return false;
                }
                arb_precision *= 2;
                if (!compact_factored_logs(logs, element,
                                                 accumulated_multiplier,
                                                 embeddings, arb_precision,
                                                 diagnostics)) {
                    return false;
                }
            }
        }

        while (true) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_step_short_element");
                if (!detail::weighted_ideal_lattice_short_element(
                            short_element, inverse_eA,
                            flint::FmpzMatConstRef(weights),
                            short_precision, diagnostics,
                            &weighted_t2_cache)) {
                    return false;
                }
            }
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_step_norm_guard");
                if (compact_short_norm_guard(
                            short_element, inverse_eA,
                            flint::FmpzMatConstRef(weights), order)) {
                    break;
                }
            }
            if (short_precision > kCompactMaxPrecision / 2) {
                return false;
            }
            short_precision *= 2;
        }

        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_subtract");
            if (!compact_de_subtract_selected(decomposition, selected,
                                                    power)) {
                return false;
            }
        }

        Ideal B1(order);
        flint::Fmpz B1_norm;
        if (!B1.is_defined()) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_multiply");
            if (!compact_multiply_integral_by_principal(
                        B1, eA, short_element, diagnostics)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_norm");
            if (!B1.norm(flint::FmpzRef(B1_norm)) ||
                flint::fmpz_cmp(flint::FmpzConstRef(B1_norm),
                                flint::FmpzConstRef(abs_discriminant)) > 0) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_step_insert_factorization");
            if (!compact_de_insert_factorization(decomposition, B1,
                                                       power)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_push");
            if (!accumulated_multiplier.push(short_element, power)) {
                return false;
            }
        }

        flint::ArbVec short_logs(places);
        FactoredElement short_factor(*field);
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_short_logs");
            if (!short_factor.is_defined() ||
                !short_factor.set_element(short_element) ||
                !short_factor.logarithmic_embedding(
                        flint::ArbVecRef(short_logs), embeddings,
                        LogEmbeddingMode::product, arb_precision,
                        diagnostics)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "class_group.compact_step_log_update");
            for (slong i = 0; i < places; ++i) {
                ::arb_addmul_si(logs.data() + i, short_logs.data() + i,
                                power, arb_precision);
            }
        }
    }

    Ideal residual(order);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_final_evaluate");
        if (!residual.is_defined() ||
            !compact_de_evaluate(residual, decomposition,
                                       diagnostics)) {
            return false;
        }
    }

    accumulated_multiplier.normalize();
    return residual_ideal.set(residual) &&
           multiplier.set(accumulated_multiplier);
}

}  // namespace

namespace detail {

class ClassGroupRelationAccess {
public:
    static bool relation_sources_consistent(
            const ClassGroupContext& context) noexcept {
        return context.relation_sources_.size() ==
               static_cast<std::size_t>(context.relation_count());
    }

    static const std::vector<char>* relation_basis_flags(
            const ClassGroupContext& context) noexcept {
        if (context.relation_basis_flags_.size() !=
            static_cast<std::size_t>(context.relation_count())) {
            return nullptr;
        }
        return &context.relation_basis_flags_;
    }

    static bool relation_row(flint::FmpzMatRef out,
                             const ClassGroupContext& context,
                             slong index) noexcept {
        return context.relations_.row(out, index);
    }

    static bool append_saturation_relation_keep(
            ClassGroupContext& context,
            const Relation& relation) noexcept {
        ClassGroupContext::RelationAppendOutcome outcome =
                ClassGroupContext::RelationAppendOutcome::none;
        return context.append_relation_with_outcome_(
                outcome, relation, ClassGroupRelationSource::Saturation,
                ClassGroupContext::DependentRelationPolicy::keep);
    }

    static bool push_relation_witnesses(
            FactoredElement& out,
            const ClassGroupContext& context,
            flint::FmpzMatConstRef coefficients,
            slong coefficient_row) noexcept {
        return silex::push_relation_witnesses(
                out, context.relations_, coefficients, coefficient_row);
    }
};

bool ClassGroupCertificationAccess::
        exact_imaginary_quadratic_class_order_for_run(
                flint::FmpzRef out,
                ClassGroupContext& context,
                flint::FmpzConstRef discriminant) noexcept {
    ClassUnitTransactionContext* const run_context =
            context.class_unit_transaction_context_;
    if (run_context != nullptr &&
        run_context->imaginary_quadratic_exact_order_valid &&
        flint::fmpz_equal(
                flint::FmpzConstRef(
                        run_context->
                                imaginary_quadratic_exact_order_discriminant),
                discriminant)) {
        flint::fmpz_set(
                out,
                flint::FmpzConstRef(
                        run_context->imaginary_quadratic_exact_order));
        return true;
    }

    if (!imaginary_quadratic_class_number(out, discriminant)) {
        return false;
    }
    if (run_context != nullptr) {
        flint::fmpz_set(
                flint::FmpzRef(
                        run_context->
                                imaginary_quadratic_exact_order_discriminant),
                discriminant);
        flint::fmpz_set(
                flint::FmpzRef(
                        run_context->imaginary_quadratic_exact_order),
                flint::FmpzConstRef(out.raw()));
        run_context->imaginary_quadratic_exact_order_valid = true;
    }
    return true;
}

bool ClassGroupCertificationAccess::
        try_certify_imaginary_quadratic_from_exact_order(
                ClassGroupContext& context,
                CertificationMode requested,
                flint::FmpzConstRef discriminant,
                flint::FmpzConstRef exact_order) noexcept {
    if (requested != CertificationMode::proven ||
        !context.has_presentation() || !context.parent_.is_maximal() ||
        context.factor_base_generation_status_ != ProofState::verified ||
        flint::fmpz_sgn(discriminant) >= 0 ||
        flint::fmpz_sgn(exact_order) <= 0 ||
        !detail::order_supports_exact_quadratic_class_certificate(
                context.parent_)) {
        return false;
    }

    flint::Fmpz context_discriminant;
    flint::Fmpz required_bound;
    flint::Fmpz class_order;
    if (!context.parent_.discriminant(
                flint::FmpzRef(context_discriminant)) ||
        !flint::fmpz_equal(
                flint::FmpzConstRef(context_discriminant), discriminant) ||
        !context.factor_base_generation_bound(
                flint::FmpzRef(required_bound)) ||
        !context.check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound)) ||
        !context.order(flint::FmpzRef(class_order)) ||
        !flint::fmpz_equal(flint::FmpzConstRef(class_order), exact_order)) {
        return false;
    }

    context.certification_ = CertificationMode::proven;
    context.relation_saturation_status_ = ProofState::verified;
    context.unit_proof_status_ = ProofState::verified;
    context.regulator_proof_status_ = ProofState::verified;
    return true;
}

bool multiply_factored_element_power_fmpz(
        FactoredElement& accumulator,
        const FactoredElement& base,
        flint::FmpzConstRef exponent) noexcept {
    return compact_multiply_power_fmpz(accumulator, base, exponent);
}

bool ClassGroupCertificationAccess::try_certify_class_unit_with_bf_audit(
        ClassGroupContext& context,
        OrderUnitGroup& units,
        flint::ArbConstRef analytic_class_regulator_product,
        flint::ArbConstRef error_bound,
        ulong cutoff,
        ulong max_cutoff,
        slong requested_precision,
        slong work_precision) noexcept {
    flint::Fmpz required_bound;
    if (!context.factor_base_generation_bound(flint::FmpzRef(required_bound)) ||
        !context.check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound))) {
        return false;
    }
    if (!context.record_analytic_class_unit_regulator_(
                units, analytic_class_regulator_product, requested_precision)) {
        return false;
    }
    return context.record_zeta_bf_audit_(error_bound, cutoff, max_cutoff,
                                         requested_precision,
                                         work_precision) &&
           context.try_promote_proven_certification_();
}

bool ClassGroupCertificationAccess::record_factor_base_honesty_proof(
        ClassGroupContext& context,
        flint::FmpzConstRef required_bound) noexcept {
    return context.record_factor_base_honesty_proof_(required_bound);
}

bool ClassGroupCertificationAccess::
        rank_zero_quadratic_class_index_bound(
                flint::FmpzRef out,
                ClassGroupContext& context,
                OrderUnitGroup& units,
                CertificationMode requested) noexcept {
    if ((requested != CertificationMode::grh &&
         requested != CertificationMode::proven) ||
        !units.is_set() ||
        units.free_rank() != 0 ||
        !same_order_parent(units.parent(), context.parent()) ||
        !context.has_presentation() || !context.parent_.is_maximal() ||
        (requested == CertificationMode::proven &&
         context.factor_base_generation_status_ != ProofState::verified)) {
        return false;
    }

    flint::Fmpz discriminant;
    if (!detail::order_supports_exact_quadratic_class_certificate(
                context.parent_) ||
        !context.parent_.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) >= 0) {
        return false;
    }

    flint::Fmpz required_bound;
    if (requested == CertificationMode::grh) {
        flint::Fmpz build_bound;
        if (!context.factor_base_build_bound(flint::FmpzRef(build_bound)) ||
            !detail::grh_factor_base_bound_with_diagnostics(
                    flint::FmpzRef(required_bound), context.parent_,
                    context.diagnostics()) ||
            flint::fmpz_cmp(flint::FmpzConstRef(build_bound),
                            flint::FmpzConstRef(required_bound)) < 0) {
            return false;
        }
    }
    flint::Fmpz candidate_order;
    flint::Fmpz exact_order;
    flint::Fmpz index;
    if ((requested == CertificationMode::proven &&
         (!context.factor_base_generation_bound(
                  flint::FmpzRef(required_bound)) ||
          !context.check_factor_base_generation_bound(
                  flint::FmpzConstRef(required_bound)))) ||
        !context.order(flint::FmpzRef(candidate_order)) ||
        !exact_imaginary_quadratic_class_order_for_run(
                flint::FmpzRef(exact_order), context,
                flint::FmpzConstRef(discriminant)) ||
        !flint::fmpz_divisible(flint::FmpzConstRef(candidate_order),
                               flint::FmpzConstRef(exact_order))) {
        return false;
    }
    flint::fmpz_divexact(flint::FmpzRef(index),
                         flint::FmpzConstRef(candidate_order),
                         flint::FmpzConstRef(exact_order));
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(index), 1) < 0) {
        return false;
    }

    if (requested == CertificationMode::proven &&
        flint::fmpz_is_one(flint::FmpzConstRef(index))) {
        context.certification_ = CertificationMode::proven;
        context.relation_saturation_status_ = ProofState::verified;
        context.unit_proof_status_ = ProofState::verified;
        context.regulator_proof_status_ = ProofState::verified;
        units.mark_certification_proven_();
    }
    flint::fmpz_set(out, flint::FmpzConstRef(index));
    return true;
}

bool ClassGroupCertificationAccess::
        saturate_relations_for_index_bound_with_units(
                bool& changed,
                bool& saturated,
                ClassGroupContext& context,
                const OrderUnitGroup& units,
                flint::FmpzConstRef index_bound,
                flint::FmpzConstRef aux_prime_bound,
                slong max_appends_per_ell,
                slong max_appends_total) noexcept {
    return context.saturate_relations_bounded_for_index_with_units_(
            changed, saturated, units, index_bound, aux_prime_bound,
            max_appends_per_ell, max_appends_total);
}

bool compact_first_reduction(
        Ideal& reduced,
        FactoredElement& multiplier,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_first_reduction");
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (base == nullptr || field == nullptr || reduced.parent() == nullptr ||
        !reduced.parent()->has_same_data(*order) ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) || n <= 1 ||
        precision <= 0 ||
        flint::fmpz_mat_nrows(decomposition_row) != 1 ||
        flint::fmpz_mat_ncols(decomposition_row) != base->length()) {
        return false;
    }

    slong max_abs_plus_one = 0;
    slong iterations = 0;
    if (!decomposition_row_max_abs_plus_one(max_abs_plus_one,
                                            decomposition_row) ||
        !floor_log_slong(iterations, max_abs_plus_one, n)) {
        return false;
    }

    std::vector<slong> powers(static_cast<std::size_t>(iterations + 1));
    for (slong k = 0; k <= iterations; ++k) {
        if (!checked_pow_slong(powers[static_cast<std::size_t>(k)], n, k)) {
            return false;
        }
    }

    Ideal active(*order);
    FactoredElement accumulated_multiplier(*field);
    if (!active.is_defined() || !accumulated_multiplier.is_defined() ||
        !active.one() || !accumulated_multiplier.one()) {
        return false;
    }

    for (slong k = iterations; k >= 0; --k) {
        const slong power = powers[static_cast<std::size_t>(k)];
        std::vector<IdealPowerEntry> support;
        support.reserve(static_cast<std::size_t>(base->length() + 1));

        for (slong j = 0; j < base->length(); ++j) {
            flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(decomposition_row, 0, j);
            if (!flint::fmpz_fits_si(entry)) {
                return false;
            }
            const slong valuation = flint::fmpz_get_si(entry);
            const slong digit = positive_mod_slong(
                    floor_div_positive_slong(valuation, power), n);
            if (digit == 0) {
                continue;
            }

            PrimeIdeal prime(*order);
            Ideal prime_ideal(*order);
            Ideal reduced_prime_power(*order);
            FactoredElement prime_multiplier(*field);
            if (!prime.is_defined() || !prime_ideal.is_defined() ||
                !reduced_prime_power.is_defined() ||
                !prime_multiplier.is_defined() ||
                !base->prime(prime, j) || !prime.get_ideal(prime_ideal) ||
                !reduce_ideal_signed_power(
                        reduced_prime_power, prime_multiplier, prime_ideal,
                        digit, precision, diagnostics) ||
                !add_ideal_power(support, reduced_prime_power, 1)) {
                return false;
            }

            slong negative_power = 0;
            if (!checked_neg_slong(negative_power, power) ||
                !factored_multiply_power_si(accumulated_multiplier,
                                            prime_multiplier,
                                            negative_power)) {
                return false;
            }
        }

        if (!add_ideal_power(support, active, n)) {
            return false;
        }

        Ideal next_active(*order);
        FactoredElement alpha(*field);
        if (!next_active.is_defined() || !alpha.is_defined() ||
            !reduce_ideal_powers(next_active, alpha, support,
                                       precision, diagnostics)) {
            return false;
        }

        slong negative_power = 0;
        if (!checked_neg_slong(negative_power, power) ||
            !factored_multiply_power_si(accumulated_multiplier, alpha,
                                        negative_power) ||
            !active.set(next_active)) {
            return false;
        }
    }

    accumulated_multiplier.normalize();
    return reduced.set(active) && multiplier.set(accumulated_multiplier);
}

bool compact_presentation_reduction(
        Ideal& residual_ideal,
        FactoredElement& multiplier,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_presentation_reduction");
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (base == nullptr || field == nullptr ||
        residual_ideal.parent() == nullptr ||
        !residual_ideal.parent()->has_same_data(*order) ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) || n <= 1 ||
        arb_precision <= 0 || short_precision <= 0 ||
        flint::fmpz_mat_nrows(decomposition_row) != 1 ||
        flint::fmpz_mat_ncols(decomposition_row) != base->length()) {
        return false;
    }

    Ideal active(*order);
    FactoredElement accumulated_multiplier(*field);
    if (!active.is_defined() || !accumulated_multiplier.is_defined() ||
        !compact_first_reduction(
                active, accumulated_multiplier, context, decomposition_row,
                element, n, short_precision, diagnostics)) {
        return false;
    }

    std::vector<IdealPowerEntry> decomposition;
    if (!compact_de_insert_factorization(decomposition, active, 1)) {
        return false;
    }

    slong max_abs_plus_one = 0;
    slong finite_bound = 0;
    if (!decomposition_row_max_abs_plus_one(max_abs_plus_one,
                                            decomposition_row) ||
        !ceil_log_slong(finite_bound,
                        max_slong_value(1, max_abs_plus_one - 1), n)) {
        return false;
    }

    return compact_infinite_reduction(
            residual_ideal, multiplier, *order, element,
            accumulated_multiplier, decomposition, finite_bound, n,
            arb_precision, short_precision, diagnostics);
}

bool compact_unit_presentation_reduction(
        Ideal& residual_ideal,
        FactoredElement& multiplier,
        const Order& order,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics,
        DiagnosticsModule module = DiagnosticsModule::unit_group,
        const char* profile_label =
                "unit_group.compact_unit_presentation_reduction")
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, module, profile_label);
    (void)module;
    (void)profile_label;
    const NumberField* field = order.parent();
    if (field == nullptr || residual_ideal.parent() == nullptr ||
        !residual_ideal.parent()->has_same_data(order) ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) ||
        element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) || n <= 1 ||
        arb_precision <= 0 || short_precision <= 0) {
        return false;
    }

    FactoredElement accumulated_multiplier(*field);
    std::vector<IdealPowerEntry> decomposition;
    if (!accumulated_multiplier.is_defined() ||
        !accumulated_multiplier.one()) {
        return false;
    }

    // reference `saturate!(U, n)` calls `is_power(a, n; decom = Dict())`.
    // For unit candidates this skips finite support reduction and runs the
    // compact-presentation infinite reduction from the unit ideal.
    return compact_infinite_reduction(
            residual_ideal, multiplier, order, element,
            accumulated_multiplier, decomposition, 0, n, arb_precision,
            short_precision, diagnostics);
}

void compact_floor_divrem_positive(
        slong& quotient,
        slong& remainder,
        slong value,
        slong divisor) noexcept {
    quotient = value / divisor;
    remainder = value % divisor;
    if (remainder < 0) {
        --quotient;
        remainder += divisor;
    }
}

bool compact_project_fmpq_poly_mod_p(
        flint::FmpzModPoly& out,
        flint::FmpqPolyConstRef input,
        const flint::FmpzModCtx& ctx) noexcept {
    if (!out.is_initialized()) {
        return false;
    }
    ::fmpz_mod_poly_zero(out.raw(), ctx.raw());
    const slong degree = ::fmpq_poly_degree(input.raw());
    if (degree < 0) {
        return true;
    }

    flint::Fmpq coeff;
    flint::Fmpz numerator;
    flint::Fmpz denominator;
    flint::Fmpz denominator_inv;
    flint::Fmpz value;
    for (slong i = 0; i <= degree; ++i) {
        ::fmpq_poly_get_coeff_fmpq(coeff.raw(), input.raw(), i);
        ::fmpz_mod_set_fmpz(numerator.raw(), fmpq_numref(coeff.raw()),
                            ctx.raw());
        ::fmpz_mod_set_fmpz(denominator.raw(), fmpq_denref(coeff.raw()),
                            ctx.raw());
        if (::fmpz_mod_is_invertible(denominator.raw(), ctx.raw()) == 0) {
            return false;
        }
        ::fmpz_mod_inv(denominator_inv.raw(), denominator.raw(), ctx.raw());
        ::fmpz_mod_mul(value.raw(), numerator.raw(),
                       denominator_inv.raw(), ctx.raw());
        ::fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, value.raw(),
                                       ctx.raw());
    }
    return true;
}

bool compact_project_element_mod_p(
        flint::FmpzModPoly& out,
        const Element& element,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx,
        const DiagnosticsContext* diagnostics) noexcept {
    flint::FmpqPoly polynomial;
    flint::FmpzModPoly projected(ctx);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_project_element_get_poly");
        if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
            return false;
        }
    }
    if (!projected.is_initialized()) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_project_element_fmpq_mod");
        if (!compact_project_fmpq_poly_mod_p(
                    projected, flint::FmpqPolyConstRef(polynomial), ctx)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_project_element_rem");
        ::fmpz_mod_poly_rem(out.raw(), projected.raw(), modulus.raw(),
                            ctx.raw());
    }
    return true;
}

bool compact_field_modulus_mod_p(
        flint::FmpzModPoly& modulus,
        const NumberField& field,
        const flint::FmpzModCtx& ctx,
        slong threshold) noexcept {
    flint::FmpzModPoly raw_modulus(ctx);
    if (!raw_modulus.is_initialized() ||
        !compact_project_fmpq_poly_mod_p(
                raw_modulus,
                flint::FmpqPolyConstRef(field.raw_flint_field()->pol),
                ctx)) {
        return false;
    }

    const slong degree = field.degree();
    if (::fmpz_mod_poly_degree(raw_modulus.raw(), ctx.raw()) != degree) {
        return false;
    }
    if (::fmpz_mod_poly_is_monic(raw_modulus.raw(), ctx.raw()) == 0) {
        ::fmpz_mod_poly_make_monic(raw_modulus.raw(), raw_modulus.raw(),
                                   ctx.raw());
    }
    if (::fmpz_mod_poly_is_squarefree(raw_modulus.raw(), ctx.raw()) == 0) {
        return false;
    }

    // reference `evaluate_mod` only uses the decomposition type to reject
    // factors with degree above `threshold`; when the whole field degree is
    // already within the threshold, that rejection is impossible.
    if (degree <= threshold) {
        ::fmpz_mod_poly_set(modulus.raw(), raw_modulus.raw(), ctx.raw());
        return true;
    }

    flint::FmpzModPolyFactor factorization(ctx);
    ::fmpz_mod_poly_factor(factorization.raw(), raw_modulus.raw(),
                           ctx.raw());
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        if (factorization.raw()->exp[i] != 1 ||
            ::fmpz_mod_poly_degree(factorization.raw()->poly + i,
                                   ctx.raw()) > threshold) {
            return false;
        }
    }

    ::fmpz_mod_poly_set(modulus.raw(), raw_modulus.raw(), ctx.raw());
    return true;
}

bool compact_field_modulus_mod_p_cached(
        bool& usable,
        flint::FmpzModPoly& modulus,
        const NumberField& field,
        const flint::FmpzModCtx& ctx,
        flint::FmpzConstRef prime,
        slong threshold,
        CompactFieldModulusCache* cache,
        const DiagnosticsContext* diagnostics) noexcept {
    usable = false;
    if (!modulus.is_initialized()) {
        return false;
    }

    if (cache == nullptr) {
        usable = compact_field_modulus_mod_p(
                modulus, field, ctx, threshold);
        return true;
    }

    if (cache->field == nullptr || !cache->field->has_same_data(field)) {
        cache->field = &field;
        cache->entries.clear();
    }

    for (const auto& entry : cache->entries) {
        if (entry.threshold != threshold ||
            ::fmpz_equal(entry.prime.raw(), prime.raw()) == 0) {
            continue;
        }
        if (!entry.usable) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_field_modulus_cache_reject");
            return true;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_field_modulus_cache_hit");
        ::fmpz_mod_poly_set_fmpz_poly(
                modulus.raw(), entry.modulus.raw(), ctx.raw());
        usable = true;
        return true;
    }

    usable = compact_field_modulus_mod_p(
            modulus, field, ctx, threshold);
    CompactFieldModulusCacheEntry entry;
    ::fmpz_set(entry.prime.raw(), prime.raw());
    entry.threshold = threshold;
    entry.usable = usable;
    if (usable) {
        ::fmpz_mod_poly_get_fmpz_poly(
                entry.modulus.raw(), modulus.raw(), ctx.raw());
    }
    cache->entries.push_back(std::move(entry));
    return true;
}

bool compact_modular_project_factored(
        flint::FmpzModPoly& out,
        const FactoredElement& input,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!out.is_initialized()) {
        return false;
    }
    ::fmpz_mod_poly_one(out.raw(), ctx.raw());

    struct ProjectedFactorPower {
        explicit ProjectedFactorPower(const flint::FmpzModCtx& context)
                noexcept
            : base(context) {
        }

        flint::FmpzModPoly base;
        flint::Fmpz exponent;
    };

    flint::FmpzModPoly base(ctx);
    flint::FmpzModPoly invertible_base(ctx);
    flint::FmpzModPoly power(ctx);
    flint::FmpzModPoly product(ctx);
    if (!base.is_initialized() || !invertible_base.is_initialized() ||
        !power.is_initialized() || !product.is_initialized()) {
        return false;
    }

    std::vector<ProjectedFactorPower> projected_factors;
    projected_factors.reserve(static_cast<std::size_t>(input.length()));
    for (const auto& entry : input.factors()) {
        if (entry.exponent == 0) {
            continue;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_project_factored_project");
            if (!compact_project_element_mod_p(
                        base, entry.factor, modulus, ctx, diagnostics)) {
                return false;
            }
        }
        bool combined = false;
        for (ProjectedFactorPower& projected : projected_factors) {
            if (::fmpz_mod_poly_equal(projected.base.raw(), base.raw(),
                                      ctx.raw()) != 0) {
                ::fmpz_add_si(projected.exponent.raw(),
                              projected.exponent.raw(), entry.exponent);
                combined = true;
                break;
            }
        }
        if (!combined) {
            projected_factors.emplace_back(ctx);
            ProjectedFactorPower& projected = projected_factors.back();
            if (!projected.base.is_initialized()) {
                return false;
            }
            ::fmpz_mod_poly_set(projected.base.raw(), base.raw(), ctx.raw());
            ::fmpz_set_si(projected.exponent.raw(), entry.exponent);
        }
    }

    flint::Fmpz absolute_exponent;
    for (const ProjectedFactorPower& projected : projected_factors) {
        if (::fmpz_is_zero(projected.exponent.raw()) != 0) {
            continue;
        }
        const flint::FmpzModPoly* exponent_base = &projected.base;
        if (::fmpz_sgn(projected.exponent.raw()) < 0) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_project_factored_invmod");
            if (::fmpz_mod_poly_invmod(
                        invertible_base.raw(), projected.base.raw(),
                        modulus.raw(), ctx.raw()) == 0) {
                return false;
            }
            exponent_base = &invertible_base;
        }
        ::fmpz_abs(absolute_exponent.raw(), projected.exponent.raw());

        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_project_factored_powmod");
            if (::fmpz_is_one(absolute_exponent.raw()) != 0) {
                ::fmpz_mod_poly_set(power.raw(), exponent_base->raw(),
                                    ctx.raw());
            } else if (::fmpz_equal_ui(absolute_exponent.raw(), 2) != 0) {
                ::fmpz_mod_poly_mulmod(
                        power.raw(), exponent_base->raw(),
                        exponent_base->raw(), modulus.raw(), ctx.raw());
            } else {
                ::fmpz_mod_poly_powmod_fmpz_binexp(
                        power.raw(), exponent_base->raw(),
                        absolute_exponent.raw(),
                        modulus.raw(), ctx.raw());
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_project_factored_mulmod");
            ::fmpz_mod_poly_mulmod(product.raw(), out.raw(), power.raw(),
                                   modulus.raw(), ctx.raw());
            ::fmpz_mod_poly_set(out.raw(), product.raw(), ctx.raw());
        }
    }

    return true;
}

bool compact_lift_modular_element(
        Element& out,
        const flint::FmpzModPoly& input,
        const NumberField& field,
        const flint::FmpzModCtx& ctx) noexcept {
    if (!out.has_parent(field)) {
        return false;
    }
    flint::FmpqPoly polynomial;
    flint::Fmpz coeff;
    const slong degree = field.degree();
    for (slong i = 0; i < degree; ++i) {
        ::fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), input.raw(), i,
                                       ctx.raw());
        ::fmpq_poly_set_coeff_fmpz(polynomial.raw(), i, coeff.raw());
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool compact_element_integer_coeff(
        flint::FmpzRef out,
        const Element& element,
        slong index) noexcept {
    flint::FmpqPoly polynomial;
    flint::Fmpq coeff;
    if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }
    ::fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial.raw(), index);
    if (::fmpz_is_one(fmpq_denref(coeff.raw())) == 0) {
        return false;
    }
    ::fmpz_set(out.raw(), fmpq_numref(coeff.raw()));
    return true;
}

void compact_inner_crt(
        flint::FmpzRef out,
        flint::FmpzConstRef left,
        flint::FmpzConstRef right,
        flint::FmpzConstRef multiplier,
        flint::FmpzConstRef modulus,
        flint::FmpzConstRef half_modulus) noexcept {
    flint::Fmpz value;
    ::fmpz_sub(value.raw(), right.raw(), left.raw());
    ::fmpz_mul(value.raw(), value.raw(), multiplier.raw());
    ::fmpz_add(value.raw(), value.raw(), left.raw());
    ::fmpz_fdiv_r(out.raw(), value.raw(), modulus.raw());
    if (::fmpz_cmp(out.raw(), half_modulus.raw()) > 0) {
        ::fmpz_sub(out.raw(), out.raw(), modulus.raw());
    }
}

bool compact_induce_inner_crt(
        Element& out,
        const Element& left,
        const Element& right,
        flint::FmpzConstRef multiplier,
        flint::FmpzConstRef modulus,
        flint::FmpzConstRef half_modulus) noexcept {
    const NumberField* field = left.parent();
    if (field == nullptr || !right.has_parent(*field) ||
        !out.has_parent(*field)) {
        return false;
    }

    flint::FmpqPoly polynomial;
    flint::Fmpz left_coeff;
    flint::Fmpz right_coeff;
    flint::Fmpz crt_coeff;
    const slong degree = field->degree();
    for (slong i = 0; i < degree; ++i) {
        if (!compact_element_integer_coeff(
                    flint::FmpzRef(left_coeff), left, i) ||
            !compact_element_integer_coeff(
                    flint::FmpzRef(right_coeff), right, i)) {
            return false;
        }
        compact_inner_crt(
                flint::FmpzRef(crt_coeff),
                flint::FmpzConstRef(left_coeff),
                flint::FmpzConstRef(right_coeff),
                multiplier, modulus, half_modulus);
        ::fmpq_poly_set_coeff_fmpz(polynomial.raw(), i, crt_coeff.raw());
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool compact_multiply_element_by_fmpz(
        Element& out,
        const Element& input,
        flint::FmpzConstRef scalar) noexcept {
    const NumberField* field = input.parent();
    if (field == nullptr || !out.has_parent(*field)) {
        return false;
    }
    Element scalar_element(*field);
    return scalar_element.is_defined() &&
           scalar_element.set_fmpz(scalar) &&
           out.multiply(input, scalar_element);
}

bool compact_element_scalar_div_fmpz(
        Element& out,
        const Element& input,
        flint::FmpzConstRef denominator) noexcept {
    const NumberField* field = input.parent();
    if (field == nullptr || !out.has_parent(*field) ||
        ::fmpz_is_zero(denominator.raw()) != 0) {
        return false;
    }
    flint::FmpqPoly polynomial;
    flint::FmpqPoly divided;
    if (!input.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }
    ::fmpq_poly_scalar_div_fmpz(divided.raw(), polynomial.raw(),
                                denominator.raw());
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(divided));
}

bool compact_element_denominator_in_order(
        flint::FmpzRef out,
        const Element& element,
        const Order& order) noexcept {
    if (order.parent() == nullptr || !element.has_parent(*order.parent())) {
        return false;
    }
    flint::FmpqMat coordinates(1, order.degree());
    if (!order.coordinates(flint::FmpqMatRef(coordinates), element)) {
        return false;
    }

    ::fmpz_one(out.raw());
    for (slong i = 0; i < order.degree(); ++i) {
        ::fmpz_lcm(out.raw(), out.raw(),
                   fmpq_denref(::fmpq_mat_entry(coordinates.raw(), 0, i)));
    }
    return true;
}

bool compact_order_denominator_adjustment(
        flint::FmpzRef out,
        const Order& order) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }
    Element generator(*field);
    flint::Fmpz generator_denominator;
    if (!generator.is_defined() || !generator.gen() ||
        !compact_element_denominator_in_order(
                flint::FmpzRef(generator_denominator), generator, order)) {
        return false;
    }
    ::fmpz_pow_ui(out.raw(), generator_denominator.raw(),
                  static_cast<ulong>(field->degree()));
    return true;
}

bool compact_mod_sym(
        Element& out,
        const Element& input,
        const Order& order,
        flint::FmpzConstRef modulus,
        flint::FmpzConstRef denominator_adjustment) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr || !input.has_parent(*field) ||
        !out.has_parent(*field) || ::fmpz_sgn(modulus.raw()) <= 0) {
        return false;
    }

    Element scaled(*field);
    OrderElement integral(order);
    flint::FmpzMat coordinates(1, order.degree());
    if (!scaled.is_defined() || !integral.is_defined() ||
        !compact_multiply_element_by_fmpz(
                scaled, input, denominator_adjustment) ||
        !integral.set_element(scaled) ||
        !integral.get_coordinates(flint::FmpzMatRef(coordinates))) {
        return false;
    }

    flint::Fmpz half_modulus;
    ::fmpz_fdiv_q_2exp(half_modulus.raw(), modulus.raw(), 1);
    for (slong i = 0; i < order.degree(); ++i) {
        fmpz* entry = ::fmpz_mat_entry(coordinates.raw(), 0, i);
        ::fmpz_fdiv_r(entry, entry, modulus.raw());
        if (::fmpz_cmp(entry, half_modulus.raw()) > 0) {
            ::fmpz_sub(entry, entry, modulus.raw());
        }
    }

    OrderElement reduced(order);
    return reduced.is_defined() &&
           reduced.set_coordinates(flint::FmpzMatConstRef(coordinates)) &&
           reduced.get_element(out);
}

bool compact_evaluate_mod_invariant_holds(
        const Element& evaluated,
        const Ideal& residual_ideal,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = residual_ideal.parent();
    if (order == nullptr) {
        return false;
    }
    FractionalIdeal principal(*order);
    FractionalIdeal residual_fractional(*order);
    if (!principal.is_defined() || !residual_fractional.is_defined()) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_evaluate_mod_principal");
        if (!principal.set_principal(evaluated, diagnostics)) {
            return false;
        }
    }
    return residual_fractional.set_integral(residual_ideal) &&
           principal.equal(residual_fractional);
}

bool compact_evaluate_mod(
        Element& out,
        const FactoredElement& compact_input,
        const Ideal& residual_ideal,
        const DiagnosticsContext* diagnostics,
        CompactFieldModulusCache* field_modulus_cache = nullptr) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_evaluate_mod");
    const Order* order = residual_ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !out.has_parent(*field) ||
        compact_input.parent() == nullptr ||
        !compact_input.parent()->has_same_data(*field)) {
        return false;
    }

    if (compact_input.length() == 0) {
        Element one(*field);
        if (!one.is_defined() || !one.one() ||
            !compact_evaluate_mod_invariant_holds(
                    one, residual_ideal, diagnostics)) {
            return false;
        }
        return out.set(one);
    }

    flint::Fmpz denominator_adjustment;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_evaluate_mod_denominator");
        if (!compact_order_denominator_adjustment(
                    flint::FmpzRef(denominator_adjustment), *order)) {
            return false;
        }
    }

    flint::Fmpz p_start;
    flint::Fmpz p;
    flint::Fmpz pp;
    ::fmpz_one(p_start.raw());
    ::fmpz_mul_2exp(p_start.raw(), p_start.raw(),
                    kAuxiliaryPrimeStartPower);
    ::fmpz_nextprime(p.raw(), p_start.raw(), 1);
    ::fmpz_one(pp.raw());

    const slong degree = field->degree();
    const slong threshold = degree > 30 ? degree / 10 : 3;
    Element reconstructed(*field);
    Element reduced(*field);
    if (!reconstructed.is_defined() || !reduced.is_defined()) {
        return false;
    }

    while (true) {
        if (::fmpz_bits(pp.raw()) >= 10000) {
            return false;
        }

        flint::FmpzModCtx ctx(p.raw());
        flint::FmpzModPoly modulus(ctx);
        flint::FmpzModPoly projected(ctx);
        if (ctx.raw() == nullptr || !modulus.is_initialized() ||
            !projected.is_initialized()) {
            return false;
        }
        bool projected_ok = true;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_evaluate_mod_field_modulus");
            if (!compact_field_modulus_mod_p_cached(
                        projected_ok, modulus, *field, ctx,
                        flint::FmpzConstRef(p), threshold,
                        field_modulus_cache, diagnostics)) {
                return false;
            }
        }
        if (projected_ok) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_evaluate_mod_project_factored");
            projected_ok = compact_modular_project_factored(
                    projected, compact_input, modulus, ctx, diagnostics);
        }
        if (!projected_ok) {
            ::fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        Element lifted(*field);
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.compact_evaluate_mod_lift");
            if (!lifted.is_defined() ||
                !compact_lift_modular_element(
                        lifted, projected, *field, ctx)) {
                return false;
            }
        }

        if (::fmpz_is_one(pp.raw()) != 0) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_evaluate_mod_initial");
                if (!reconstructed.set(lifted) ||
                    !compact_mod_sym(
                            reduced, reconstructed, *order, p,
                            flint::FmpzConstRef(denominator_adjustment))) {
                    return false;
                }
            }
            ::fmpz_set(pp.raw(), p.raw());
        } else {
            flint::Fmpz p2;
            flint::Fmpz pp_inverse_mod_p;
            flint::Fmpz multiplier;
            flint::Fmpz half_modulus;
            Element last(*field);
            if (!last.is_defined() || !last.set(reduced)) {
                return false;
            }
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_evaluate_mod_crt");
                ::fmpz_mul(p2.raw(), pp.raw(), p.raw());
                if (::fmpz_invmod(pp_inverse_mod_p.raw(),
                                  pp.raw(), p.raw()) == 0) {
                    ::fmpz_nextprime(p.raw(), p.raw(), 1);
                    continue;
                }
                ::fmpz_mul(multiplier.raw(), pp.raw(),
                           pp_inverse_mod_p.raw());
                ::fmpz_fdiv_q_2exp(half_modulus.raw(), p2.raw(), 1);
                if (!compact_induce_inner_crt(
                            reconstructed, reconstructed, lifted,
                            flint::FmpzConstRef(multiplier),
                            flint::FmpzConstRef(p2),
                            flint::FmpzConstRef(half_modulus))) {
                    return false;
                }
            }
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_evaluate_mod_mod_sym");
                if (!compact_mod_sym(
                            reduced, reconstructed, *order,
                            flint::FmpzConstRef(p2),
                            flint::FmpzConstRef(denominator_adjustment))) {
                    return false;
                }
            }
            bool stable = false;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.compact_evaluate_mod_stable");
                stable = reduced.equal(last);
            }
            if (stable) {
                Element divided(*field);
                {
                    SILEX_PROFILE_SCOPE(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.compact_evaluate_mod_finish");
                    if (!divided.is_defined() ||
                        !compact_element_scalar_div_fmpz(
                                divided, reduced,
                                flint::FmpzConstRef(
                                        denominator_adjustment)) ||
                        !compact_evaluate_mod_invariant_holds(
                                divided, residual_ideal, diagnostics)) {
                        return false;
                    }
                }
                return out.set(divided);
            }
            ::fmpz_set(pp.raw(), p2.raw());
        }

        ::fmpz_nextprime(p.raw(), p.raw(), 1);
    }
}

bool compact_split_power(
        FactoredElement& quotient,
        FactoredElement& remainder,
        const FactoredElement& compact,
        slong n) noexcept {
    if (quotient.parent() == nullptr || remainder.parent() == nullptr ||
        compact.parent() == nullptr ||
        !quotient.parent()->has_same_data(*compact.parent()) ||
        !remainder.parent()->has_same_data(*compact.parent()) || n <= 1 ||
        !quotient.one() || !remainder.one()) {
        return false;
    }

    for (const auto& entry : compact.factors()) {
        slong q = 0;
        slong r = 0;
        compact_floor_divrem_positive(q, r, entry.exponent, n);
        if (q != 0 && !quotient.push(entry.factor, q)) {
            return false;
        }
        if (r != 0 && !remainder.push(entry.factor, r)) {
            return false;
        }
    }

    quotient.normalize();
    remainder.normalize();
    return true;
}

bool compact_power_root_from_presentation(
        bool& is_power,
        FactoredElement& root,
        const Order& order,
        const FactoredElement& element,
        const Ideal& residual,
        const FactoredElement& multiplier,
        slong n,
        const DiagnosticsContext* diagnostics,
        CompactFieldModulusCache* field_modulus_cache = nullptr) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.compact_power_root_from_presentation");
    is_power = false;
    const NumberField* field = order.parent();
    if (field == nullptr || root.parent() == nullptr ||
        !root.parent()->has_same_data(*field) || element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) ||
        residual.parent() == nullptr ||
        !residual.parent()->has_same_data(order) ||
        multiplier.parent() == nullptr ||
        !multiplier.parent()->has_same_data(*field) || n <= 1) {
        return false;
    }

    FactoredElement compact_input(*field);
    Element evaluated(*field);
    FactoredElement inverse_multiplier(*field);
    FactoredElement compact(*field);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_prepare");
        if (!compact_input.is_defined() || !evaluated.is_defined() ||
            !inverse_multiplier.is_defined() || !compact.is_defined() ||
            !compact_input.multiply(element, multiplier)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_evaluate");
        if (!compact_evaluate_mod(
                    evaluated, compact_input, residual, diagnostics,
                    field_modulus_cache)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_compose");
        if (!inverse_multiplier.invert(multiplier) ||
            !compact.set(inverse_multiplier) ||
            !compact.push(evaluated, 1)) {
            return false;
        }
        compact.normalize();
    }

    FactoredElement quotient(*field);
    FactoredElement remainder(*field);
    Element remainder_value(*field);
    Element scaled_remainder(*field);
    Element remainder_root(*field);
    Element denominator_element(*field);
    flint::Fmpz exponent;
    flint::Fmpz denominator;
    flint::Fmpz adjusted_denominator;
    flint::Fmpz denominator_power;
    bool remainder_is_power = false;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_split");
        if (!quotient.is_defined() || !remainder.is_defined() ||
            !remainder_value.is_defined() || !scaled_remainder.is_defined() ||
            !remainder_root.is_defined() ||
            !denominator_element.is_defined() ||
            !compact_split_power(quotient, remainder, compact, n)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_remainder_evaluate");
        if (!remainder.evaluate(remainder_value)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_denominator");
        if (!compact_element_denominator_in_order(
                    flint::FmpzRef(denominator), remainder_value, order)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_scale");
        if (::fmpz_root(adjusted_denominator.raw(), denominator.raw(),
                        static_cast<ulong>(n)) == 0) {
            ::fmpz_set(adjusted_denominator.raw(), denominator.raw());
        }
        ::fmpz_pow_ui(denominator_power.raw(), adjusted_denominator.raw(),
                      static_cast<ulong>(n));
        if (!compact_multiply_element_by_fmpz(
                    scaled_remainder, remainder_value,
                    flint::FmpzConstRef(denominator_power))) {
            return false;
        }
    }
    flint::fmpz_set_si(flint::FmpzRef(exponent), n);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_element_is_power");
        if (!scaled_remainder.is_power(remainder_is_power, remainder_root,
                                       flint::FmpzConstRef(exponent))) {
            return false;
        }
    }
    if (!remainder_is_power) {
        return true;
    }

    FactoredElement remainder_root_factor(*field);
    FactoredElement denominator_factor(*field);
    FactoredElement candidate(*field);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.compact_power_root_build");
        if (!remainder_root_factor.is_defined() ||
            !denominator_factor.is_defined() || !candidate.is_defined() ||
            !remainder_root_factor.set_element(remainder_root) ||
            !candidate.multiply(quotient, remainder_root_factor)) {
            return false;
        }
        if (::fmpz_is_one(adjusted_denominator.raw()) == 0) {
            if (!denominator_element.set_fmpz(
                        flint::FmpzConstRef(adjusted_denominator)) ||
                !denominator_factor.one() ||
                !denominator_factor.push(denominator_element, -1) ||
                !candidate.multiply(candidate, denominator_factor)) {
                return false;
            }
        }
    }

    is_power = true;
    root.swap(candidate);
    return true;
}

bool compact_presentation_power_root(
        bool& is_power,
        FactoredElement& root,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "class_group.compact_presentation_power_root");
    is_power = false;
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (base == nullptr || field == nullptr || root.parent() == nullptr ||
        !root.parent()->has_same_data(*field) || element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) || n <= 1 ||
        arb_precision <= 0 || short_precision <= 0 ||
        flint::fmpz_mat_nrows(decomposition_row) != 1 ||
        flint::fmpz_mat_ncols(decomposition_row) != base->length()) {
        return false;
    }

    Ideal residual(*order);
    FactoredElement multiplier(*field);
    if (!residual.is_defined() || !multiplier.is_defined() ||
        !compact_presentation_reduction(
                residual, multiplier, context, decomposition_row, element, n,
                arb_precision, short_precision, diagnostics)) {
        return false;
    }

    return compact_power_root_from_presentation(
            is_power, root, *order, element, residual, multiplier, n,
            diagnostics);
}

bool compact_unit_presentation_power_root(
        bool& is_power,
        FactoredElement& root,
        const Order& order,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics,
        CompactFieldModulusCache* field_modulus_cache) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.compact_unit_presentation_power_root");
    is_power = false;
    const NumberField* field = order.parent();
    if (field == nullptr || root.parent() == nullptr ||
        !root.parent()->has_same_data(*field) || element.parent() == nullptr ||
        !element.parent()->has_same_data(*field) || n <= 1 ||
        arb_precision <= 0 || short_precision <= 0) {
        return false;
    }

    FactoredElement structural(*field);
    FactoredElement leftover(*field);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.compact_unit_formal_split");
        if (!structural.is_defined() || !leftover.is_defined() ||
            !compact_split_power(structural, leftover, element, n)) {
            return false;
        }
    }

    // reference `is_power(::FacElem, n)` first removes the formal exponent
    // quotient. It tests the residual factor only when its conjugate-log
    // bound is small relative to the original input; otherwise it computes
    // the compact presentation of the original factored element.
    if (leftover.length() == 0) {
        is_power = true;
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.compact_unit_structural_power");
        root.swap(structural);
        return true;
    }

    bool use_residual = false;
    if (!compact_use_residual_power_test(
                use_residual, order, element, leftover, diagnostics)) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.compact_unit_power_split_bound_failure");
        return false;
    }

    const FactoredElement& compact_input = use_residual ? leftover : element;
    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::unit_group,
            use_residual
                    ? "unit_group.compact_unit_residual_branch"
                    : "unit_group.compact_unit_original_branch");

    Ideal residual(order);
    FactoredElement multiplier(*field);
    if (!residual.is_defined() || !multiplier.is_defined() ||
        !compact_unit_presentation_reduction(
                residual, multiplier, order, compact_input, n, arb_precision,
                short_precision, diagnostics)) {
        return false;
    }

    FactoredElement compact_root(*field);
    if (!compact_root.is_defined() ||
        !compact_power_root_from_presentation(
                is_power, compact_root, order, compact_input, residual,
                multiplier, n, diagnostics, field_modulus_cache)) {
        return false;
    }
    if (!is_power) {
        return true;
    }

    if (!use_residual) {
        root.swap(compact_root);
        return true;
    }

    FactoredElement candidate(*field);
    if (!candidate.is_defined() ||
        !candidate.multiply(structural, compact_root)) {
        return false;
    }

    root.swap(candidate);
    return true;
}

ulong unit_candidate_random_next(ulong& state) noexcept {
    state = state * UWORD(6364136223846793005) +
            UWORD(1442695040888963407);
    return state;
}

slong unit_candidate_random_index(ulong& state,
                                        slong length) noexcept {
    if (length <= 0) {
        return -1;
    }
    return static_cast<slong>(
            unit_candidate_random_next(state) %
            static_cast<ulong>(length));
}

slong unit_candidate_used_count(
        const RelationUnitCandidateBatchState& state) noexcept {
    return static_cast<slong>(state.relations_used.size());
}

bool unit_candidate_relation_used(
        const RelationUnitCandidateBatchState& state,
        slong relation_index) noexcept {
    return std::find(state.relations_used.begin(),
                     state.relations_used.end(),
                     relation_index) != state.relations_used.end();
}

void unit_candidate_sync_relation_state(
        RelationUnitCandidateBatchState& state,
        const std::vector<slong>& extra_indices) noexcept {
    const slong rel_count = static_cast<slong>(extra_indices.size());
    if (state.extra_relation_count > rel_count ||
        static_cast<slong>(state.relations_used.size()) > rel_count) {
        state.relations_used.clear();
    } else {
        state.relations_used.erase(
                std::remove_if(
                        state.relations_used.begin(),
                        state.relations_used.end(),
                        [&extra_indices](slong relation_index) noexcept {
                            return !std::binary_search(
                                    extra_indices.begin(),
                                    extra_indices.end(), relation_index);
                        }),
                state.relations_used.end());
    }
    state.extra_relation_count = rel_count;
}

void try_select_unit_candidate_extra(
        std::vector<slong>& selected,
        std::vector<char>& selected_flags,
        RelationUnitCandidateBatchState& state,
        slong& used_count,
        const std::vector<slong>& extra_indices) noexcept {
    const slong rel_count = static_cast<slong>(extra_indices.size());
    const slong index =
            unit_candidate_random_index(state.random_state, rel_count);
    if (index < 0) {
        return;
    }

    const std::size_t pos = static_cast<std::size_t>(index);
    if (used_count != rel_count) {
        if (selected_flags[pos] != 0 ||
            unit_candidate_relation_used(state, extra_indices[pos])) {
            return;
        }
        state.relations_used.push_back(extra_indices[pos]);
        ++used_count;
    } else if (selected_flags[pos] != 0) {
        return;
    }

    selected_flags[pos] = 1;
    selected.push_back(extra_indices[pos]);
}

bool unit_candidate_extra_indices(
        std::vector<slong>& out,
        RelationUnitCandidateBatchState& state,
        const std::vector<slong>& extra_indices,
        slong unit_rank_value,
        slong add) noexcept {
    out.clear();
    if (add < 0) {
        return false;
    }

    const slong rel_count = static_cast<slong>(extra_indices.size());
    if (rel_count <= 0) {
        state.relations_used.clear();
        state.extra_relation_count = 0;
        return true;
    }
    // reference `find_candidates` stores selected `rel_gens` row numbers in
    // `u.relations_used`, not positions in a flag buffer.  Track native
    // relation-row identities for the same recycle-after-all-rows semantics.
    unit_candidate_sync_relation_state(state, extra_indices);

    slong nrel = max_slong_value(WORD(10), unit_rank_value);
    nrel = min_slong_value(nrel, rel_count);

    std::vector<char> selected_flags(extra_indices.size(), 0);
    slong used_count = unit_candidate_used_count(state);
    out.reserve(static_cast<std::size_t>(
            min_slong_value(rel_count, nrel + add)));

    while (static_cast<slong>(out.size()) < nrel) {
        try_select_unit_candidate_extra(out, selected_flags, state,
                                              used_count, extra_indices);
    }

    for (slong i = 0; i < add; ++i) {
        if (static_cast<slong>(out.size()) == rel_count) {
            break;
        }
        try_select_unit_candidate_extra(out, selected_flags, state,
                                              used_count, extra_indices);
    }

    return true;
}

bool saturate_full_row_rank_lattice(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef rows) noexcept {
    const slong row_count = flint::fmpz_mat_nrows(rows);
    const slong col_count = flint::fmpz_mat_ncols(rows);
    if (row_count <= 0 || col_count < row_count ||
        flint::fmpz_mat_rank(rows) != row_count) {
        return false;
    }

    // Source trace: reference `Sparse/HNF.jl:saturate` computes
    // `Hti = transpose(hnf(transpose(A)))`, keeps the leading square block,
    // solves against `transpose(A)`, then transposes the solution.  FLINT's
    // dense solve returns an integral numerator plus a common denominator, so
    // divide exactly when the reference solution is integral.
    flint::FmpzMat transpose(col_count, row_count);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(transpose), rows);

    flint::FmpzMat hnf(col_count, row_count);
    ::fmpz_mat_hnf(hnf.raw(), transpose.raw());

    flint::FmpzMat hti(row_count, col_count);
    fmpz_mat_transpose(hti.raw(), hnf.raw());

    flint::FmpzMat leading(row_count, row_count);
    for (slong i = 0; i < row_count; ++i) {
        for (slong j = 0; j < row_count; ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(leading, i, j),
                            flint::fmpz_mat_entry(
                                    flint::FmpzMatConstRef(hti), i, j));
        }
    }

    flint::Fmpz den;
    flint::FmpzMat saturated(row_count, col_count);
    if (::fmpz_mat_solve(saturated.raw(), den.raw(), leading.raw(),
                         rows.raw()) == 0 ||
        flint::fmpz_is_zero(flint::FmpzConstRef(den))) {
        return false;
    }
    if (!flint::fmpz_is_one(flint::FmpzConstRef(den))) {
        for (slong i = 0; i < row_count; ++i) {
            for (slong j = 0; j < col_count; ++j) {
                if (!flint::fmpz_divisible(
                            flint::fmpz_mat_entry(
                                    flint::FmpzMatConstRef(saturated), i, j),
                            flint::FmpzConstRef(den))) {
                    return false;
                }
            }
        }
        ::fmpz_mat_scalar_divexact_fmpz(saturated.raw(), saturated.raw(),
                                        den.raw());
    }
    out = std::move(saturated);
    return true;
}

bool relation_basis_extra_indices(std::vector<slong>& basis_indices,
                                        std::vector<slong>& extra_indices,
                                        const ClassGroupContext& context)
        noexcept {
    basis_indices.clear();
    extra_indices.clear();
    const std::vector<char>* relation_basis_flags =
            ClassGroupRelationAccess::relation_basis_flags(context);
    if (!context.has_factor_base() || relation_basis_flags == nullptr) {
        return false;
    }

    const slong generators = context.generator_count();
    if (generators <= 0 || context.relation_rank() < generators) {
        return true;
    }

    basis_indices.reserve(static_cast<std::size_t>(generators));
    extra_indices.reserve(static_cast<std::size_t>(
            max_slong_value(WORD(0), context.relation_count() - generators)));
    for (slong i = 0; i < context.relation_count(); ++i) {
        if ((*relation_basis_flags)[static_cast<std::size_t>(i)] != 0) {
            basis_indices.push_back(i);
        } else {
            extra_indices.push_back(i);
        }
    }

    if (static_cast<slong>(basis_indices.size()) != generators) {
        basis_indices.clear();
        extra_indices.clear();
    }
    return true;
}

bool saturation_relation_data(std::vector<FactoredElement>& relations,
                                    flint::FmpzMatRef relation_rows,
                                    const ClassGroupContext& context) noexcept {
    relations.clear();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!context.has_factor_base() || field == nullptr ||
        !ClassGroupRelationAccess::relation_sources_consistent(context)) {
        return false;
    }

    std::vector<slong> basis_indices;
    std::vector<slong> extra_indices;
    if (!relation_basis_extra_indices(basis_indices, extra_indices,
                                            context)) {
        return false;
    }

    const slong generators = context.generator_count();
    const bool use_storage_order =
            basis_indices.empty() && extra_indices.empty() &&
            context.relation_count() > 0 && context.relation_rank() < generators;
    const slong total_rows = use_storage_order
            ? context.relation_count()
            : static_cast<slong>(basis_indices.size() + extra_indices.size());
    if (flint::fmpz_mat_nrows(relation_rows) != total_rows ||
        flint::fmpz_mat_ncols(relation_rows) != generators) {
        return false;
    }
    relations.reserve(static_cast<std::size_t>(total_rows));

    Element generator(*field);
    flint::FmpzMat row(1, generators);
    auto append_index = [&](slong source_index, slong dest_index) noexcept {
        if (!ClassGroupRelationAccess::relation_row(
                    flint::FmpzMatRef(row), context, source_index) ||
            !context.relation_generator(generator, source_index)) {
            return false;
        }
        for (slong j = 0; j < generators; ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(relation_rows, dest_index, j),
                            flint::fmpz_mat_entry(
                                    flint::FmpzMatConstRef(row), 0, j));
        }

        relations.emplace_back(*field);
        if (!relations.back().is_defined() ||
            !relations.back().set_element(generator)) {
            return false;
        }
        return true;
    };

    slong dest = 0;
    if (use_storage_order) {
        for (slong index = 0; index < context.relation_count(); ++index) {
            if (!append_index(index, dest)) {
                return false;
            }
            ++dest;
        }
    } else {
        for (slong index : basis_indices) {
            if (!append_index(index, dest)) {
                return false;
            }
            ++dest;
        }
        for (slong index : extra_indices) {
            if (!append_index(index, dest)) {
                return false;
            }
            ++dest;
        }
    }

    return true;
}

bool simplified_saturation_relation_data(
        std::vector<FactoredElement>& relations,
        flint::FmpzMat& relation_rows,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        slong n,
        slong precision) noexcept {
    relations.clear();
    relation_rows = flint::FmpzMat(0, 0);
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (base == nullptr || field == nullptr || !units.is_set() ||
        !same_order_parent(units.parent(), order) ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*field) ||
        precision <= 0 || n <= 1 || !n_is_prime(static_cast<ulong>(n)) ||
        context.relation_rank() < context.generator_count() ||
        context.relation_count() < context.generator_count()) {
        return false;
    }

    const slong generators = context.generator_count();
    const slong relation_count = context.relation_count();
    flint::FmpzMat ordered_rows(relation_count, generators);
    std::vector<FactoredElement> ordered_relations;
    if (!saturation_relation_data(
                ordered_relations, flint::FmpzMatRef(ordered_rows),
                context)) {
        return false;
    }

    flint::FmpzMat hnf_full(0, 0);
    flint::FmpzMat transform(0, 0);
    if (!dense_hnf_kannan_bachem_transform(
                hnf_full, transform, flint::FmpzMatConstRef(ordered_rows))) {
        return false;
    }
    const slong hnf_rows = flint::fmpz_mat_nrows(
            flint::FmpzMatConstRef(hnf_full));
    if (hnf_rows != generators) {
        return false;
    }

    std::vector<slong> selected_rows;
    selected_rows.reserve(static_cast<std::size_t>(hnf_rows));
    for (slong i = 0; i < hnf_rows; ++i) {
        bool skip = false;
        if (!hnf_row_first_nonzero_is_one(
                    skip, flint::FmpzMatConstRef(hnf_full), i)) {
            return false;
        }
        if (skip) {
            continue;
        }
        selected_rows.push_back(i);
    }

    flint::FmpzMat selected(
            static_cast<slong>(selected_rows.size()), generators);
    relations.reserve(selected_rows.size());
    for (std::size_t out_row = 0; out_row < selected_rows.size(); ++out_row) {
        const slong hnf_row = selected_rows[out_row];
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(selected),
                    static_cast<slong>(out_row),
                    flint::FmpzMatConstRef(hnf_full), hnf_row)) {
            return false;
        }

        relations.emplace_back(*field);
        if (!relations.back().is_defined() || !relations.back().one()) {
            return false;
        }
        for (slong j = 0; j < relation_count; ++j) {
            flint::FmpzConstRef exponent =
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(transform), hnf_row, j);
            if (flint::fmpz_is_zero(exponent)) {
                continue;
            }
            if (!compact_multiply_power_fmpz(
                        relations.back(),
                        ordered_relations[static_cast<std::size_t>(j)],
                        exponent)) {
                return false;
            }
        }
        relations.back().normalize();
    }

    if (!reduce_relation_units_modulo(relations, units, embeddings, precision)) {
        return false;
    }

    relation_rows = std::move(selected);
    return true;
}

bool simplified_saturation_context(
        ClassGroupContext& simplified,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        const RelationUnitExtractionState& unit_state,
        slong n,
        slong precision) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.saturation_simplify");
    simplified.clear();
    const FactorBase* base = context.factor_base();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (base == nullptr || field == nullptr || !units.is_set() ||
        !same_order_parent(units.parent(), order)) {
        return false;
    }

    std::vector<FactoredElement> relations;
    flint::FmpzMat rows(0, 0);
    if (!simplified_saturation_relation_data(
                relations, rows, context, units, embeddings, n, precision)) {
        return false;
    }

    ClassGroupContext candidate(*order);
    candidate.set_diagnostics(context.diagnostics());
    if (!candidate.is_defined() || !candidate.set_factor_base(*base)) {
        return false;
    }

    const FactorBase* candidate_base = candidate.factor_base();
    if (candidate_base == nullptr) {
        return false;
    }

    Element generator(*field);
    for (slong i = 0; i < static_cast<slong>(relations.size()); ++i) {
        flint::FmpzMatConstWindow row(
                flint::FmpzMatConstRef(rows), i, 0, i + 1,
                context.generator_count());
        Relation relation(*candidate_base);
        if (!relation.is_defined() ||
            !relations[static_cast<std::size_t>(i)].evaluate(generator) ||
            !set_relation_from_known_row(relation, *candidate_base, generator,
                                         row.const_ref()) ||
            !ClassGroupRelationAccess::append_saturation_relation_keep(
                    candidate, relation)) {
            return false;
        }
    }

    flint::FmpzMat zero_row(1, context.generator_count());
    for (slong i = 0; i < units.free_rank(); ++i) {
        FactoredElement unit(*field);
        Relation relation(*candidate_base);
        // reference simplify appends its retained U.units directly.  Strict native
        // already owns the same basis in exact expanded form until a unit
        // mutation clears the dependent-unit cache.
        const Element* relation_generator =
                unit_state.expanded_reduce_mod_units_log_generator(
                        i, units.free_rank());
        if (relation_generator != nullptr &&
            !relation_generator->has_parent(*field)) {
            return false;
        }
        if (relation_generator == nullptr) {
            if (!unit.is_defined() || !units.free_generator(unit, i) ||
                !unit.evaluate(generator)) {
                return false;
            }
            relation_generator = &generator;
        }
        if (!relation.is_defined() ||
            !set_relation_from_known_row(
                    relation, *candidate_base, *relation_generator,
                    flint::FmpzMatConstRef(zero_row)) ||
            !ClassGroupRelationAccess::append_saturation_relation_keep(
                    candidate, relation)) {
            return false;
        }
    }

    simplified.swap(candidate);
    return true;
}

bool saturation_candidates(flint::FmpzMat& candidates,
                                 slong& auxiliary_prime_count,
                                 const ClassGroupContext& context,
                                 slong n,
                                 double stable,
                                 SaturationCandidateSearchResult*
                                         search_result) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.saturation_candidates");
    auxiliary_prime_count = 0;
    candidates = flint::FmpzMat(0, 0);
    const Order* order = context.parent();
    if (order == nullptr || context.generator_count() <= 0 ||
        context.relation_count() <= 0) {
        return false;
    }

    flint::FmpzMat relation_rows(context.relation_count(),
                                 context.generator_count());
    std::vector<FactoredElement> relations;
    if (!saturation_relation_data(
                relations, flint::FmpzMatRef(relation_rows), context)) {
        return false;
    }

    flint::Fmpz start;
    default_saturation_auxiliary_prime_start(start);
    return saturation_candidates_from_relations(
            candidates, auxiliary_prime_count, relations, *order, n, stable,
            flint::FmpzConstRef(start), WORD_MAX, search_result,
            diagnostics);
}

bool saturation_process_candidates_once(
        SaturationCandidateProcessingResult& result,
        ClassGroupContext& context,
        flint::FmpzMatConstRef candidates,
        slong n,
        bool append_class_relations,
        OrderUnitGroup* units,
        EmbeddingContext* embeddings,
        RelationUnitExtractionState* unit_state,
        slong precision,
        slong max_class_relation_appends) noexcept {
    return saturation_process_candidates_once(
            result, context, context, candidates, n, append_class_relations,
            units, embeddings, unit_state, precision,
            max_class_relation_appends);
}

bool saturation_process_candidates_once(
        SaturationCandidateProcessingResult& result,
        const ClassGroupContext& input_context,
        ClassGroupContext& append_context,
        flint::FmpzMatConstRef candidates,
        slong n,
        bool append_class_relations,
        OrderUnitGroup* units,
        EmbeddingContext* embeddings,
        RelationUnitExtractionState* unit_state,
        slong precision,
        slong max_class_relation_appends) noexcept {
    const DiagnosticsContext* diagnostics = input_context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.saturation_process_candidates");
    result = SaturationCandidateProcessingResult{};
    const FactorBase* base = input_context.factor_base();
    const Order* order = input_context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (base == nullptr || field == nullptr || n <= 1 ||
        !n_is_prime(static_cast<ulong>(n)) ||
        input_context.generator_count() <= 0 ||
        input_context.relation_count() <= 0 ||
        (append_class_relations &&
         (!same_factor_base(append_context.factor_base(), base) ||
          max_class_relation_appends < 0))) {
        return false;
    }

    const slong candidate_rows = flint::fmpz_mat_nrows(candidates);
    const slong candidate_cols = flint::fmpz_mat_ncols(candidates);
    if (candidate_rows == 0) {
        result.empty_candidate_space = true;
        return candidate_cols == input_context.relation_count() ||
               candidate_cols == input_context.relation_count() + 1;
    }

    flint::FmpzMat relation_rows(input_context.relation_count(),
                                 input_context.generator_count());
    std::vector<FactoredElement> relations;
    if (!saturation_relation_data(
                relations, flint::FmpzMatRef(relation_rows), input_context)) {
        return false;
    }
    const slong relation_count = static_cast<slong>(relations.size());
    if (!(candidate_rows == relation_count ||
          candidate_rows == relation_count + 1)) {
        return false;
    }

    for (slong column = candidate_cols - 1; column >= 0; --column) {
        if (append_class_relations &&
            result.appended_class_relations >= max_class_relation_appends) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.relation_saturation.class_relation_append_cap");
            break;
        }
        ++result.inspected_columns;
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.relation_saturation.candidate_inspected");

        FactoredElement element(*field);
        flint::FmpzMat fac_row(1, input_context.generator_count());
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturation_reconstruct_candidate");
            if (!saturation_reconstruct_candidate(
                        element, flint::FmpzMatRef(fac_row), relations,
                        flint::FmpzMatConstRef(relation_rows), candidates,
                        column, *field)) {
                return false;
            }
        }

        bool zero_row = false;
        if (!fmpz_mat_row_is_zero(zero_row, flint::FmpzMatConstRef(fac_row))) {
            return false;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                zero_row
                        ? "class_group.relation_saturation.candidate_zero_row"
                        : "class_group.relation_saturation.candidate_relation_row");

        flint::FmpzMat divided_row(1, input_context.generator_count());
        if (!zero_row) {
            bool divisible = false;
            if (!fmpz_mat_row_divexact_si(
                        divided_row, divisible,
                        flint::FmpzMatConstRef(fac_row), n)) {
                return false;
            }
            if (!divisible) {
                return false;
            }
        }

        FactoredElement root(*field);
        bool is_power = false;
        ++result.tested_power_candidates;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturation_candidate_power_test");
            // reference's `is_power(::FacElem, n; decom=fac_a)` first removes the
            // formal exponent quotient.  If the residual factor is small
            // enough in the conjugate-log bound, reference tests the residual
            // compact presentation without the ideal decomposition row.
            FactoredElement structural(*field);
            FactoredElement leftover(*field);
            if (!structural.is_defined() || !leftover.is_defined() ||
                !compact_split_power(structural, leftover, element, n)) {
                return false;
            }
            if (leftover.length() == 0) {
                is_power = true;
                ++result.structural_power_candidates;
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.saturation_structural_power");
                root.swap(structural);
            } else {
                ++result.compact_power_candidates;
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        "class_group.saturation_compact_power_test");

                bool use_residual = false;
                if (!compact_use_residual_power_test(
                            use_residual, *order, element, leftover,
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.saturation_power_split_bound",
                            "class_group.saturation_original_log_bound",
                            "class_group.saturation_residual_log_bound")) {
                    return false;
                }
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::class_group,
                        use_residual
                                ? "class_group.saturation_residual_power_branch"
                                : "class_group.saturation_original_power_branch");

                if (use_residual) {
                    Ideal residual(*order);
                    FactoredElement multiplier(*field);
                    FactoredElement compact_root(*field);
                    if (!residual.is_defined() || !multiplier.is_defined() ||
                        !compact_root.is_defined() ||
                        !compact_unit_presentation_reduction(
                                residual, multiplier, *order, leftover, n,
                                precision, precision, diagnostics,
                                DiagnosticsModule::class_group,
                                "class_group.saturation_residual_presentation_reduction") ||
                        !compact_power_root_from_presentation(
                                is_power, compact_root, *order, leftover,
                                residual, multiplier, n, diagnostics)) {
                        return false;
                    }
                    if (is_power) {
                        FactoredElement candidate(*field);
                        if (!candidate.is_defined() ||
                            !candidate.multiply(structural, compact_root)) {
                            return false;
                        }
                        root.swap(candidate);
                    }
                } else if (!compact_presentation_power_root(
                                   is_power, root, input_context,
                                   flint::FmpzMatConstRef(fac_row), element, n,
                                   precision, precision, diagnostics)) {
                    return false;
                }
            }
        }
        if (!is_power) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturation_candidate_not_power");
            result.saw_non_power = true;
            result.wasted = true;
            break;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.saturation_candidate_power");
        result.saw_power = true;

        if (zero_row) {
            result.saw_dependent_unit = true;
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturation_dependent_unit_candidate");
            if (units != nullptr || embeddings != nullptr ||
                unit_state != nullptr || precision != 0) {
                if (units == nullptr || embeddings == nullptr ||
                    unit_state == nullptr || precision <= 0 ||
                    !same_order_parent(units->parent(), order)) {
                    return false;
                }

                bool dependent_changed = false;
                if (!add_dependent_unit(
                            dependent_changed, *units, root, *embeddings,
                            *unit_state, precision)) {
                    return false;
                }
                if (dependent_changed) {
                    ++result.appended_dependent_units;
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.saturation_dependent_unit_appended");
                    unit_state->has_expected_regulator = false;
                } else {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::class_group,
                            "class_group.saturation_dependent_unit_unchanged");
                }
            }
            continue;
        }

        result.saw_class_relation = true;
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "class_group.saturation_class_relation_candidate");

        if (!append_class_relations) {
            continue;
        }

        Element generator(*field);
        const FactorBase* append_base = append_context.factor_base();
        if (append_base == nullptr) {
            return false;
        }
        Relation relation(*append_base);
        if (!generator.is_defined() || !relation.is_defined() ||
            !root.evaluate(generator) ||
            !set_relation_from_known_row(
                    relation, *append_base, generator,
                    flint::FmpzMatConstRef(divided_row))) {
            return false;
        }

        const slong relation_count_before = append_context.relation_count();
        if (!ClassGroupRelationAccess::append_saturation_relation_keep(
                    append_context, relation)) {
            return false;
        }
        if (append_context.relation_count() > relation_count_before) {
            ++result.appended_class_relations;
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturation_class_relation_appended");
            if (unit_state != nullptr) {
                unit_state->has_expected_regulator = false;
            }
        }
    }

    return true;
}

bool saturate_class_unit_context(
        bool& success,
        SaturationCandidateProcessingResult& last_result,
        ClassGroupContext& context,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        RelationUnitExtractionState& unit_state,
        slong n,
        double stable,
        slong precision) noexcept {
    const DiagnosticsContext* diagnostics = context.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "class_group.saturate_class_unit_context");
    success = false;
    last_result = SaturationCandidateProcessingResult{};
    if (n <= 1 || !n_is_prime(static_cast<ulong>(n)) ||
        !std::isfinite(stable) || stable <= 0.0 || precision <= 0 ||
        !units.is_set() || !same_order_parent(units.parent(), context.parent())) {
        return false;
    }

    double work_stable = stable;
    for (;;) {
        ClassGroupContext simplified;
        if (!simplified_saturation_context(
                    simplified, context, units, embeddings, unit_state, n,
                    precision)) {
            return false;
        }

        flint::FmpzMat candidates(0, 0);
        slong auxiliary_prime_count = 0;
        SaturationCandidateSearchResult search_result;
        if (!saturation_candidates(candidates, auxiliary_prime_count,
                                         simplified, n, work_stable,
                                         &search_result)) {
            return false;
        }
        (void)auxiliary_prime_count;

        auto record_candidate_generation =
                [&](SaturationCandidateProcessingResult& result)
                        noexcept {
                    result.simplified_relation_count =
                            simplified.relation_count();
                    result.simplified_generator_count =
                            simplified.generator_count();
                    result.candidate_rows = flint::fmpz_mat_nrows(
                            flint::FmpzMatConstRef(candidates));
                    result.candidate_cols = flint::fmpz_mat_ncols(
                            flint::FmpzMatConstRef(candidates));
                    result.auxiliary_prime_count = auxiliary_prime_count;
                    result.candidate_search = search_result;
                };


        if (flint::fmpz_mat_nrows(flint::FmpzMatConstRef(candidates)) == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturate_relation_lattice.empty_candidate_space");
            last_result.empty_candidate_space = true;
            record_candidate_generation(last_result);
            return true;
        }

        SaturationCandidateProcessingResult result;
        if (!saturation_process_candidates_once(
                    result, simplified, context,
                    flint::FmpzMatConstRef(candidates), n, true, &units,
                    &embeddings, &unit_state, precision)) {
            return false;
        }
        record_candidate_generation(result);
        last_result = result;
        if (result.appended_class_relations > 0 ||
            result.appended_dependent_units > 0) {
            success = true;
        }
        if (result.wasted) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturate_relation_lattice.wasted_candidate");
            if (work_stable > std::numeric_limits<double>::max() / 2.0) {
                return false;
            }
            work_stable *= 2.0;
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::class_group,
                    "class_group.saturate_relation_lattice.stable_doubled");
            continue;
        }

        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                success
                        ? "class_group.saturate_relation_lattice.success"
                        : "class_group.saturate_relation_lattice.no_new_power");
        return true;
    }
}


bool unit_candidate_witnesses(std::vector<FactoredElement>& out,
                                    RelationUnitCandidateBatchState& state,
                                    const ClassGroupContext& context,
                                    slong add)
        noexcept {
    out.clear();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!context.has_factor_base() || field == nullptr ||
        !ClassGroupRelationAccess::relation_sources_consistent(context) ||
        ClassGroupRelationAccess::relation_basis_flags(context) == nullptr) {
        return false;
    }

    const slong generators = context.generator_count();
    if (generators <= 0 || context.relation_rank() < generators) {
        return true;
    }

    std::vector<slong> basis_indices;
    std::vector<slong> extra_indices;
    if (!relation_basis_extra_indices(basis_indices, extra_indices,
                                            context)) {
        return false;
    }
    if (static_cast<slong>(basis_indices.size()) != generators) {
        state.relations_used.clear();
        state.extra_relation_count = 0;
        return true;
    }
    if (extra_indices.empty()) {
        state.relations_used.clear();
        state.extra_relation_count = 0;
        return true;
    }

    slong unit_rank_value = 0;
    if (!unit_rank(unit_rank_value, *field)) {
        return false;
    }
    std::vector<slong> selected_extra_indices;
    if (!unit_candidate_extra_indices(
                selected_extra_indices, state, extra_indices,
                unit_rank_value, add)) {
        return false;
    }
    const slong extra_count =
            static_cast<slong>(selected_extra_indices.size());
    if (extra_count <= 0) {
        return true;
    }
    const slong total_rows = generators + extra_count;
    flint::FmpzMat basis_rows(generators, generators);
    flint::FmpzMat extra_rows(extra_count, generators);
    flint::FmpzMat row(1, generators);
    for (slong i = 0; i < generators; ++i) {
        if (!ClassGroupRelationAccess::relation_row(
                    flint::FmpzMatRef(row), context, basis_indices[i])) {
            return false;
        }
        for (slong j = 0; j < generators; ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(basis_rows, i, j),
                            flint::fmpz_mat_entry(
                                    flint::FmpzMatConstRef(row), 0, j));
        }
    }
    for (slong i = 0; i < extra_count; ++i) {
        if (!ClassGroupRelationAccess::relation_row(
                    flint::FmpzMatRef(row), context,
                    selected_extra_indices[i])) {
            return false;
        }
        for (slong j = 0; j < generators; ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(extra_rows, i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(row), 0, j));
        }
    }

    if (flint::fmpz_mat_rank(flint::FmpzMatConstRef(basis_rows)) !=
        generators) {
        return true;
    }

    // Source trace: reference `FindUnits.jl:find_candidates` calls
    // `k, d = solve_dixon_sf(x.M.bas_gens, rel)` and then saturates
    // `hcat(k, (-d)*I)`.  Solve the equivalent transposed dense FLINT system
    // `basis_rows^T * k^T = d * extra_rows^T`.
    flint::FmpzMat basis_transpose(generators, generators);
    flint::FmpzMat extra_transpose(generators, extra_count);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(basis_transpose),
                              flint::FmpzMatConstRef(basis_rows));
    flint::fmpz_mat_transpose(flint::FmpzMatRef(extra_transpose),
                              flint::FmpzMatConstRef(extra_rows));

    flint::Fmpz denominator;
    flint::FmpzMat solution_transpose(generators, extra_count);
    if (::fmpz_mat_solve_dixon_den(
                solution_transpose.raw(), denominator.raw(),
                basis_transpose.raw(), extra_transpose.raw()) == 0 ||
        flint::fmpz_is_zero(flint::FmpzConstRef(denominator))) {
        return false;
    }

    if (flint::fmpz_sgn(flint::FmpzConstRef(denominator)) < 0) {
        flint::fmpz_neg(flint::FmpzRef(denominator),
                        flint::FmpzConstRef(denominator));
        ::fmpz_mat_neg(solution_transpose.raw(), solution_transpose.raw());
    }

    flint::FmpzMat solution(extra_count, generators);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(solution),
                              flint::FmpzMatConstRef(solution_transpose));

    flint::FmpzMat kernel_seed(extra_count, total_rows);
    flint::fmpz_mat_zero(flint::FmpzMatRef(kernel_seed));
    for (slong i = 0; i < extra_count; ++i) {
        for (slong j = 0; j < generators; ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(kernel_seed, i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(solution), i, j));
        }
        flint::fmpz_neg(
                flint::fmpz_mat_entry(kernel_seed, i, generators + i),
                flint::FmpzConstRef(denominator));
    }

    flint::FmpzMat kernel_rows(0, 0);
    if (!saturate_full_row_rank_lattice(
                kernel_rows, flint::FmpzMatConstRef(kernel_seed))) {
        return false;
    }
    const slong kernel_count =
            flint::fmpz_mat_nrows(flint::FmpzMatConstRef(kernel_rows));

    if (kernel_count > 1) {
        flint::FmpzMat transform(kernel_count, kernel_count);
        flint::FmpzLll lll;
        fmpz_lll(kernel_rows.raw(), transform.raw(), lll.raw());
    }

    flint::FmpzMat coefficients(1, context.relation_count());
    out.reserve(static_cast<std::size_t>(kernel_count));
    for (slong i = 0; i < kernel_count; ++i) {
        flint::fmpz_mat_zero(flint::FmpzMatRef(coefficients));
        bool nonzero = false;
        for (slong j = 0; j < generators; ++j) {
            flint::FmpzConstRef value =
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(kernel_rows), i, j);
            if (!flint::fmpz_is_zero(value)) {
                nonzero = true;
            }
            flint::fmpz_set(
                    flint::fmpz_mat_entry(coefficients, 0, basis_indices[j]),
                    value);
        }
        for (slong j = 0; j < extra_count; ++j) {
            flint::FmpzConstRef value =
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(kernel_rows), i,
                            generators + j);
            if (!flint::fmpz_is_zero(value)) {
                nonzero = true;
            }
            flint::fmpz_set(
                    flint::fmpz_mat_entry(
                            coefficients, 0, selected_extra_indices[j]),
                    value);
        }
        if (!nonzero) {
            continue;
        }

        out.emplace_back(*field);
        if (!ClassGroupRelationAccess::push_relation_witnesses(
                    out.back(), context,
                    flint::FmpzMatConstRef(coefficients), 0)) {
            return false;
        }
    }

    return true;
}

}  // namespace detail

bool class_regulator_index_is_one(const OrderUnitGroup& units,
                                  const ClassGroupContext& class_group,
                                  flint::ArbConstRef analytic_hR,
                                  slong precision) noexcept {
    if (precision <= 0 || !flint::arb_is_finite(analytic_hR) ||
        !flint::arb_is_positive(analytic_hR)) {
        return false;
    }

    flint::Arb candidate_hR;
    if (!units.class_regulator_product(
                flint::ArbRef(candidate_hR), class_group, precision)) {
        return false;
    }
    return detail::class_regulator_index_is_one_from_candidate_product(
            flint::ArbConstRef(candidate_hR), analytic_hR, precision,
            class_group.diagnostics());
}

void ClassGroupContext::reset_certification_metadata_() noexcept {
    certification_ = CertificationMode::unknown;
    relation_saturation_status_ = ProofState::not_checked;
    if (private_storage_ != nullptr) {
        private_storage_->relation_saturation_records.clear();
        private_storage_->relation_saturation_proof_records.clear();
    }
    analytic_class_regulator_status_ = ProofState::not_checked;
    zeta_bf_status_ = ProofState::not_checked;
    zeta_bf_cutoff_ = 0;
    zeta_bf_max_cutoff_ = 0;
    zeta_bf_requested_precision_ = 0;
    zeta_bf_work_precision_ = 0;
    flint::arb_zero(flint::ArbRef(zeta_bf_error_bound_));
    unit_proof_status_ = ProofState::not_checked;
    regulator_proof_status_ = ProofState::not_checked;
}

void ClassGroupContext::reset_factor_base_generation_check_() noexcept {
    if (private_storage_ != nullptr) {
        private_storage_->factor_base_generation_records.clear();
    }
    flint::fmpz_zero(
            flint::FmpzRef(factor_base_generation_checked_bound_));
    factor_base_generation_checked_status_ = ProofState::not_checked;
}

bool ClassGroupContext::mark_factor_base_generation_check_(
        flint::FmpzConstRef checked_bound,
        ProofState status) noexcept {
    if (flint::fmpz_sgn(checked_bound) < 0 ||
        !valid_factor_base_generation_record_status(status)) {
        return false;
    }

    if (status == ProofState::verified &&
        (factor_base_generation_status_ != ProofState::verified ||
         flint::fmpz_cmp(checked_bound,
                         flint::FmpzConstRef(factor_base_build_bound_)) >
                 0)) {
        return false;
    }

    if (status == ProofState::unavailable &&
        factor_base_generation_checked_status_ == ProofState::verified) {
        return true;
    }

    if (flint::fmpz_cmp(
                checked_bound,
                flint::FmpzConstRef(
                        factor_base_generation_checked_bound_)) >= 0 ||
        factor_base_generation_checked_status_ != ProofState::verified) {
        flint::fmpz_set(
                flint::FmpzRef(factor_base_generation_checked_bound_),
                checked_bound);
        factor_base_generation_checked_status_ = status;
    }

    return true;
}

bool ClassGroupContext::record_factor_base_generation_(
        flint::FmpzConstRef build_bound) noexcept {
    flint::fmpz_set(flint::FmpzRef(factor_base_build_bound_), build_bound);
    flint::fmpz_zero(flint::FmpzRef(factor_base_generation_bound_));
    factor_base_generation_status_ = ProofState::unavailable;
    reset_factor_base_generation_check_();

    flint::Fmpz required;
    if (parent_.is_defined() &&
        factor_base_class_group_bound(flint::FmpzRef(required), parent_)) {
        flint::fmpz_set(flint::FmpzRef(factor_base_generation_bound_),
                        flint::FmpzConstRef(required));
        if (flint::fmpz_cmp(build_bound, flint::FmpzConstRef(required)) >= 0) {
            factor_base_generation_status_ = ProofState::verified;
        }
    }

    return true;
}

bool ClassGroupContext::record_factor_base_honesty_proof_(
        flint::FmpzConstRef required_bound) noexcept {
    flint::Fmpz recorded_bound;
    if (!has_factor_base() || !has_presentation() || !parent_.is_maximal() ||
        flint::fmpz_sgn(required_bound) < 0 ||
        !factor_base_generation_bound(flint::FmpzRef(recorded_bound)) ||
        flint::fmpz_cmp(required_bound,
                        flint::FmpzConstRef(recorded_bound)) < 0) {
        return false;
    }

    reset_factor_base_generation_check_();
    factor_base_generation_status_ = ProofState::verified;
    flint::fmpz_set(
            flint::FmpzRef(factor_base_generation_checked_bound_),
            flint::FmpzConstRef(recorded_bound));
    factor_base_generation_checked_status_ = ProofState::verified;
    return true;
}

bool ClassGroupContext::relation_saturation_proof_prereqs_verified_()
        const noexcept {
    return quotient_.is_defined() &&
           factor_base_generation_checked_status_ == ProofState::verified &&
           unit_proof_status_ == ProofState::verified &&
           regulator_proof_status_ == ProofState::verified;
}

bool ClassGroupContext::record_relation_saturation_proof_(
        flint::FmpzConstRef ell,
        ProofState status,
        slong rank,
        slong target_rank,
        slong local_primes) noexcept {
    if (!ensure_private_storage_() || !flint::fmpz_is_prime(ell) ||
        !valid_relation_saturation_record_status(status) || rank < 0 ||
        target_rank <= 0 || local_primes < 0 ||
        (status == ProofState::verified && rank < target_rank)) {
        return false;
    }

    for (detail::RelationSaturationProofRecord& record :
         private_storage_->relation_saturation_proof_records) {
        if (!flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell)) {
            continue;
        }
        record.status = status;
        record.rank = rank;
        record.target_rank = target_rank;
        record.local_primes = local_primes;
        return true;
    }

    detail::RelationSaturationProofRecord record;
    flint::fmpz_set(flint::FmpzRef(record.ell), ell);
    record.status = status;
    record.rank = rank;
    record.target_rank = target_rank;
    record.local_primes = local_primes;
    private_storage_->relation_saturation_proof_records.push_back(
            std::move(record));
    return true;
}

bool ClassGroupContext::relation_saturation_proof_verified_(
        flint::FmpzConstRef ell) const noexcept {
    if (private_storage_ == nullptr) {
        return false;
    }
    for (const detail::RelationSaturationProofRecord& record :
         private_storage_->relation_saturation_proof_records) {
        if (flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell) &&
            record.status == ProofState::verified) {
            return true;
        }
    }
    return false;
}

bool ClassGroupContext::mark_relation_saturation_verified_(
        flint::FmpzConstRef ell) noexcept {
    return mark_relation_saturation_(ell, ProofState::verified);
}

bool ClassGroupContext::mark_relation_saturation_(
        flint::FmpzConstRef ell,
        ProofState status) noexcept {
    if (!ensure_private_storage_() || !flint::fmpz_is_prime(ell) ||
        !valid_relation_saturation_record_status(status)) {
        return false;
    }

    for (detail::RelationSaturationRecord& record :
         private_storage_->relation_saturation_records) {
        if (flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell)) {
            record.status = status;
            return true;
        }
    }

    detail::RelationSaturationRecord record;
    flint::fmpz_set(flint::FmpzRef(record.ell), ell);
    record.status = status;
    private_storage_->relation_saturation_records.push_back(std::move(record));
    return true;
}

bool ClassGroupContext::complete_relation_saturation_proof_(
        flint::FmpzConstRef ell) noexcept {
    if (!flint::fmpz_is_prime(ell)) {
        return false;
    }

    std::vector<flint::Fmpz> required_ells;
    required_ells.emplace_back();
    flint::fmpz_set(flint::FmpzRef(required_ells.back()), ell);
    return complete_relation_saturation_proof_(required_ells);
}

bool ClassGroupContext::complete_relation_saturation_proof_(
        const std::vector<flint::Fmpz>& required_ells) noexcept {
    if (!ensure_private_storage_() || required_ells.empty()) {
        return false;
    }

    private_storage_->relation_saturation_records.clear();
    for (const flint::Fmpz& ell : required_ells) {
        if (!flint::fmpz_is_prime(flint::FmpzConstRef(ell)) ||
            !relation_saturation_proof_verified_(flint::FmpzConstRef(ell))) {
            relation_saturation_status_ = ProofState::unavailable;
            return true;
        }
    }

    for (const flint::Fmpz& ell : required_ells) {
        if (!mark_relation_saturation_verified_(flint::FmpzConstRef(ell))) {
            relation_saturation_status_ = ProofState::unavailable;
            private_storage_->relation_saturation_records.clear();
            return true;
        }
    }

    relation_saturation_status_ = ProofState::verified;
    return true;
}

bool ClassGroupContext::relation_saturation_proof_complete_()
        const noexcept {
    if (private_storage_ == nullptr ||
        relation_saturation_status_ != ProofState::verified ||
        private_storage_->relation_saturation_records.empty()) {
        return false;
    }

    for (const detail::RelationSaturationRecord& record :
         private_storage_->relation_saturation_records) {
        if (record.status != ProofState::verified ||
            !relation_saturation_proof_verified_(
                    flint::FmpzConstRef(record.ell))) {
            return false;
        }
    }
    return true;
}

bool ClassGroupContext::quadratic_completeness_verified_() const noexcept {
    return parent_.degree() == 2 &&
           relation_saturation_status_ == ProofState::verified;
}

bool append_factor_base_generation_record(
        std::vector<detail::FactorBaseGenerationRecord>& records,
        flint::FmpzConstRef p,
        ProofState status) noexcept {
    if (!flint::fmpz_is_prime(p) ||
        !valid_factor_base_generation_record_status(status)) {
        return false;
    }

    detail::FactorBaseGenerationRecord record;
    flint::fmpz_set(flint::FmpzRef(record.p), p);
    record.status = status;
    records.push_back(std::move(record));
    return true;
}

bool relation_saturation_beta_rows(
        std::vector<FactoredElement>& beta_rows,
        const ClassGroupContext& context,
        flint::FmpzConstRef ell) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!context.has_presentation() || field == nullptr ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    std::vector<FactoredElement> candidate;
    candidate.reserve(static_cast<std::size_t>(context.invariant_count()));
    flint::Fmpz invariant;
    for (slong i = 0; i < context.invariant_count(); ++i) {
        if (!context.invariant(flint::FmpzRef(invariant), i)) {
            return false;
        }
        if (!flint::fmpz_divisible(flint::FmpzConstRef(invariant), ell)) {
            continue;
        }

        candidate.emplace_back(*field);
        if (!candidate.back().is_defined() ||
            !context.invariant_generator_power_witness(candidate.back(), i)) {
            return false;
        }
    }

    beta_rows.swap(candidate);
    return true;
}

bool relation_saturation_torsion_needed(bool& out,
                                        const OrderUnitGroup& units,
                                        flint::FmpzConstRef ell) noexcept {
    if (!units.is_set() || !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::Fmpz torsion_order;
    if (!units.torsion_order(flint::FmpzRef(torsion_order))) {
        return false;
    }

    out = flint::fmpz_divisible(flint::FmpzConstRef(torsion_order), ell);
    return true;
}

bool relation_saturation_beta_dlog_column(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& beta_rows,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        bool include_torsion,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (order == nullptr || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !units.is_set() || !same_order_parent(units.parent(), order) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    bool torsion_needed = false;
    if (!relation_saturation_torsion_needed(torsion_needed, units, ell) ||
        include_torsion != torsion_needed) {
        return false;
    }

    const slong beta_len = static_cast<slong>(beta_rows.size());
    const slong free_rank = units.free_rank();
    const slong rows = beta_len + free_rank + (include_torsion ? 1 : 0);
    if (rows <= 0) {
        return false;
    }

    ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::FmpzMat candidate(rows, 1);
    ResidueFieldElement image(residue_field);
    for (slong i = 0; i < beta_len; ++i) {
        if (!image.set_factored_element(
                    beta_rows[static_cast<std::size_t>(i)]) ||
            !image.quotient_log_mod_prime(
                    flint::fmpz_mat_entry(candidate, i, 0), ell)) {
            return false;
        }
    }

    for (slong i = 0; i < free_rank; ++i) {
        FactoredElement generator(*field);
        const slong row = beta_len + i;
        if (!units.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !image.quotient_log_mod_prime(
                    flint::fmpz_mat_entry(candidate, row, 0), ell)) {
            return false;
        }
    }

    if (include_torsion) {
        OrderElement torsion(*order);
        if (!units.torsion_generator(torsion) ||
            !image.set_order_element(torsion) ||
            !image.quotient_log_mod_prime(
                    flint::fmpz_mat_entry(candidate, beta_len + free_rank, 0),
                    ell)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

bool relation_saturation_proof_prime_dlog_column(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& beta_rows,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        bool include_torsion,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    const Order* order = context.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order) ||
        prime.residue_degree() != 1) {
        return false;
    }

    flint::Fmpz p;
    flint::Fmpz pminus;
    if (!prime.rational_prime(flint::FmpzRef(p))) {
        return false;
    }
    flint::fmpz_sub_ui(flint::FmpzRef(pminus), flint::FmpzConstRef(p), 1);
    if (!flint::fmpz_divisible(flint::FmpzConstRef(pminus), ell)) {
        return false;
    }

    return relation_saturation_beta_dlog_column(out, beta_rows, context, units,
                                                include_torsion, prime, ell);
}

bool append_dlog_column(flint::FmpzMat& out,
                        const std::vector<flint::FmpzMat>& columns,
                        const flint::FmpzMat& column,
                        slong rows) noexcept {
    if (rows <= 0 || flint::fmpz_mat_nrows(column) != rows ||
        flint::fmpz_mat_ncols(column) != 1) {
        return false;
    }

    const slong old_cols = static_cast<slong>(columns.size());
    flint::FmpzMat candidate(rows, old_cols + 1);
    for (slong j = 0; j < old_cols; ++j) {
        const flint::FmpzMat& old_column =
                columns[static_cast<std::size_t>(j)];
        if (flint::fmpz_mat_nrows(old_column) != rows ||
            flint::fmpz_mat_ncols(old_column) != 1) {
            return false;
        }
        for (slong i = 0; i < rows; ++i) {
            flint::fmpz_set(flint::fmpz_mat_entry(candidate, i, j),
                            flint::fmpz_mat_entry(old_column, i, 0));
        }
    }
    for (slong i = 0; i < rows; ++i) {
        flint::fmpz_set(flint::fmpz_mat_entry(candidate, i, old_cols),
                        flint::fmpz_mat_entry(column, i, 0));
    }

    out = std::move(candidate);
    return true;
}

static_assert(std::is_nothrow_default_constructible_v<
              detail::ClassGroupContextStorage>);

ClassGroupContext::ClassGroupContext() noexcept = default;

ClassGroupContext::ClassGroupContext(const Order& order) noexcept {
    define(order);
}

ClassGroupContext::~ClassGroupContext() noexcept = default;

ClassGroupContext::ClassGroupContext(ClassGroupContext&& other) noexcept {
    swap(other);
}

ClassGroupContext& ClassGroupContext::operator=(
        ClassGroupContext&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ClassGroupContext::swap(ClassGroupContext& other) noexcept {
    parent_.swap(other.parent_);
    base_.swap(other.base_);
    relations_.swap(other.relations_);
    relation_sources_.swap(other.relation_sources_);
    relation_basis_flags_.swap(other.relation_basis_flags_);
    generator_relation_scratch_.swap(other.generator_relation_scratch_);
    private_storage_.swap(other.private_storage_);
    quotient_.swap(other.quotient_);
    row_module_.swap(other.row_module_);
    relation_module_.swap(other.relation_module_);
    std::swap(relation_rank_, other.relation_rank_);
    std::swap(row_module_synced_, other.row_module_synced_);
    std::swap(skipped_dependent_relations_,
              other.skipped_dependent_relations_);
    std::swap(relation_kernel_units_target_,
              other.relation_kernel_units_target_);
    std::swap(generator_relation_policy_,
              other.generator_relation_policy_);
    std::swap(certification_, other.certification_);
    factor_base_build_bound_.swap(other.factor_base_build_bound_);
    factor_base_generation_bound_.swap(other.factor_base_generation_bound_);
    factor_base_generation_checked_bound_.swap(
            other.factor_base_generation_checked_bound_);
    std::swap(factor_base_generation_status_,
              other.factor_base_generation_status_);
    std::swap(factor_base_generation_checked_status_,
              other.factor_base_generation_checked_status_);
    std::swap(relation_saturation_status_,
              other.relation_saturation_status_);
    std::swap(analytic_class_regulator_status_,
              other.analytic_class_regulator_status_);
    std::swap(zeta_bf_status_, other.zeta_bf_status_);
    std::swap(zeta_bf_cutoff_, other.zeta_bf_cutoff_);
    std::swap(zeta_bf_max_cutoff_, other.zeta_bf_max_cutoff_);
    std::swap(zeta_bf_requested_precision_,
              other.zeta_bf_requested_precision_);
    std::swap(zeta_bf_work_precision_, other.zeta_bf_work_precision_);
    zeta_bf_error_bound_.swap(other.zeta_bf_error_bound_);
    std::swap(unit_proof_status_, other.unit_proof_status_);
    std::swap(regulator_proof_status_, other.regulator_proof_status_);
    std::swap(has_base_, other.has_base_);
    std::swap(diagnostics_, other.diagnostics_);
    row_module_.set_diagnostics(diagnostics_);
    relation_module_.set_diagnostics(diagnostics_);
    other.row_module_.set_diagnostics(other.diagnostics_);
    other.relation_module_.set_diagnostics(other.diagnostics_);
}

void ClassGroupContext::clear() noexcept {
    parent_.clear();
    base_.clear();
    relations_.clear();
    relation_sources_.clear();
    relation_basis_flags_.clear();
    generator_relation_scratch_.clear();
    if (private_storage_ != nullptr) {
        private_storage_->clear();
    }
    quotient_.clear();
    row_module_ = fmpz_smat::HnfContext();
    row_module_.set_diagnostics(diagnostics_);
    relation_module_.clear();
    relation_module_.set_diagnostics(diagnostics_);
    relation_rank_ = 0;
    row_module_synced_ = true;
    skipped_dependent_relations_ = 0;
    relation_kernel_units_target_ = 0;
    generator_relation_policy_ = DependentRelationPolicy::native;
    certification_ = CertificationMode::unknown;
    flint::fmpz_zero(flint::FmpzRef(factor_base_build_bound_));
    flint::fmpz_zero(flint::FmpzRef(factor_base_generation_bound_));
    flint::fmpz_zero(flint::FmpzRef(factor_base_generation_checked_bound_));
    factor_base_generation_status_ = ProofState::not_checked;
    factor_base_generation_checked_status_ = ProofState::not_checked;
    relation_saturation_status_ = ProofState::not_checked;
    analytic_class_regulator_status_ = ProofState::not_checked;
    zeta_bf_status_ = ProofState::not_checked;
    zeta_bf_cutoff_ = 0;
    zeta_bf_max_cutoff_ = 0;
    zeta_bf_requested_precision_ = 0;
    zeta_bf_work_precision_ = 0;
    flint::arb_zero(flint::ArbRef(zeta_bf_error_bound_));
    unit_proof_status_ = ProofState::not_checked;
    regulator_proof_status_ = ProofState::not_checked;
    has_base_ = false;
    diagnostics_ = nullptr;
}

bool ClassGroupContext::ensure_private_storage_() noexcept {
    if (private_storage_ == nullptr) {
        private_storage_.reset(
                new (std::nothrow) detail::ClassGroupContextStorage());
    }
    return private_storage_ != nullptr;
}

bool ClassGroupContext::define(const Order& order) noexcept {
    if (!detail::order_has_parented_basis(order)) {
        return false;
    }

    ClassGroupContext next;
    if (!next.ensure_private_storage_()) {
        return false;
    }
    next.diagnostics_ = diagnostics_;
    next.parent_ = order;
    next.base_ = FactorBase(order);
    if (!next.base_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool ClassGroupContext::is_defined() const noexcept {
    return private_storage_ != nullptr &&
           detail::order_has_parented_basis(parent_) && base_.is_defined();
}

const Order* ClassGroupContext::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

void ClassGroupContext::set_diagnostics(
        const DiagnosticsContext* diagnostics) noexcept {
    diagnostics_ = diagnostics;
    row_module_.set_diagnostics(diagnostics_);
    relation_module_.set_diagnostics(diagnostics_);
}

const DiagnosticsContext* ClassGroupContext::diagnostics() const noexcept {
    return diagnostics_;
}

slong ClassGroupContext::analytic_finish_precision() const noexcept {
    return private_storage_ == nullptr
            ? 0
            : private_storage_->analytic_finish_precision;
}

bool ClassGroupContext::analytic_finish_product(
        flint::ArbRef out,
        slong& precision) const noexcept {
    if (private_storage_ == nullptr ||
        !private_storage_->analytic_finish_product_valid ||
        private_storage_->analytic_finish_product_precision <= 0) {
        return false;
    }
    flint::arb_set(
            out,
            flint::ArbConstRef(
                    private_storage_->analytic_finish_product));
    precision = private_storage_->analytic_finish_product_precision;
    return true;
}

bool ClassGroupContext::sync_row_module_checkpoint_() noexcept {
    if (!has_factor_base()) {
        return false;
    }
    if (row_module_synced_) {
        return true;
    }

    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.relation_module_checkpoint");
    if (!relation_module_.checkpoint_context(row_module_)) {
        return false;
    }
    row_module_.set_diagnostics(diagnostics_);
    relation_rank_ = row_module_.rank();
    row_module_synced_ = true;
    return true;
}

bool ClassGroupContext::has_factor_base() const noexcept {
    return is_defined() && has_base_ && base_.is_defined();
}

const FactorBase* ClassGroupContext::factor_base() const noexcept {
    return has_factor_base() ? &base_ : nullptr;
}

slong ClassGroupContext::generator_count() const noexcept {
    return has_factor_base() ? base_.length() : 0;
}

bool ClassGroupContext::factor_base_prime(PrimeIdeal& out,
                                          slong index) const noexcept {
    return has_factor_base() && base_.prime(out, index);
}

bool ClassGroupContext::factor_base_prime_is_principal(
        bool& out,
        slong index) const noexcept {
    if (!has_factor_base() || index < 0 || index >= generator_count()) {
        return false;
    }

    bool is_principal = false;
    flint::FmpzMat row(1, generator_count());
    for (slong i = 0; i < relation_count(); ++i) {
        if (!relations_.row(flint::FmpzMatRef(row), i)) {
            return false;
        }

        slong nonzero_col = -1;
        slong nonzero_count = 0;
        for (slong j = 0; j < generator_count(); ++j) {
            flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(flint::FmpzMatConstRef(row), 0, j);
            if (flint::fmpz_is_zero(entry)) {
                continue;
            }
            nonzero_col = j;
            ++nonzero_count;
            if (nonzero_count > 1) {
                break;
            }
        }

        if (nonzero_count != 1 || nonzero_col != index) {
            continue;
        }

        flint::FmpzConstRef entry =
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(row),
                                      0, nonzero_col);
        if (flint::fmpz_equal_si(entry, 1) ||
            flint::fmpz_equal_si(entry, -1)) {
            is_principal = true;
            break;
        }
    }

    out = is_principal;
    return true;
}

bool ClassGroupContext::factor_base_prime_is_hnf_covered(
        bool& out,
        slong index) noexcept {
    if (!has_factor_base() || index < 0 || index >= generator_count()) {
        return false;
    }
    if (!sync_row_module_checkpoint_()) {
        return false;
    }

    bool is_covered = false;
    const slong rank = row_module_.rank();
    if (rank > 0) {
        flint::FmpzMat hnf(rank, generator_count());
        if (!row_module_.get_hnf_rows(flint::FmpzMatRef(hnf))) {
            return false;
        }

        for (slong i = 0; i < rank; ++i) {
            slong pivot_col = -1;
            for (slong j = 0; j < generator_count(); ++j) {
                if (!flint::fmpz_is_zero(
                            flint::fmpz_mat_entry(
                                    flint::FmpzMatConstRef(hnf), i, j))) {
                    pivot_col = j;
                    break;
                }
            }
            if (pivot_col == index) {
                is_covered = true;
                break;
            }
        }
    }

    out = is_covered;
    return true;
}

bool ClassGroupContext::build_factor_base(
        flint::FmpzConstRef bound) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.build_factor_base");
    if (!is_defined()) {
        return false;
    }

    ClassGroupContext candidate(parent_);
    candidate.set_diagnostics(diagnostics_);
    if (!candidate.is_defined() ||
        !candidate.base_.build(bound, candidate.diagnostics_)) {
        return false;
    }
    candidate.relations_ = RelationMatrix(candidate.base_);
    if (!candidate.relations_.is_defined() ||
        !reset_relation_modules(candidate.row_module_,
                                candidate.relation_module_,
                                candidate.base_.length(),
                                candidate.diagnostics_)) {
        return false;
    }

    candidate.has_base_ = true;
    if (!candidate.record_factor_base_generation_(bound)) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ClassGroupContext::build_relation_factor_base_(
        flint::FmpzConstRef bound) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.build_relation_factor_base");
    if (!is_defined()) {
        return false;
    }

    ClassGroupContext candidate(parent_);
    candidate.set_diagnostics(diagnostics_);
    const bool oriented_quadratic_base =
            class_unit_transaction_context_ != nullptr &&
            class_unit_transaction_context_->audit.policy.selected &&
            class_unit_transaction_context_->audit.policy.degree == 2 &&
            class_unit_transaction_context_->audit.policy.unit_rank == 0 &&
            class_unit_transaction_context_->audit.policy.relations ==
                    detail::NativeRelationStrategy::relation_completion_table;
    bool factor_base_built = false;
    if (candidate.is_defined()) {
        factor_base_built = oriented_quadratic_base
                ? candidate.base_.
                          build_maximal_imaginary_quadratic_relation_base(
                                  bound, candidate.diagnostics_)
                : candidate.base_.build_relation_completion_base(bound);
    }
    if (!factor_base_built) {
        return false;
    }
    candidate.relations_ = RelationMatrix(candidate.base_);
    if (!candidate.relations_.is_defined() ||
        !reset_relation_modules(candidate.row_module_,
                                candidate.relation_module_,
                                candidate.base_.length(),
                                candidate.diagnostics_)) {
        return false;
    }

    candidate.has_base_ = true;
    // reference's FBquad scan classifies each p once and retains the canonical
    // ramified/split representative.  This route has constructed that same
    // representative (with the generic decomposition path as its atomic
    // fallback), so its successful build is the canonical-bound coverage
    // receipt; do not decompose and compare every prime a second time.
    if (!candidate.record_factor_base_generation_(bound) ||
        (oriented_quadratic_base &&
         candidate.factor_base_generation_status_ == ProofState::verified &&
         !candidate.mark_factor_base_generation_check_(
                 flint::FmpzConstRef(
                         candidate.factor_base_generation_bound_),
                 ProofState::verified))) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ClassGroupContext::build_search_factor_base_(
        flint::FmpzConstRef bound,
        bool incomplete) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.build_search_factor_base");
    if (!is_defined()) {
        return false;
    }

    ClassGroupContext candidate(parent_);
    candidate.set_diagnostics(diagnostics_);
    if (!candidate.is_defined()) {
        return false;
    }
    const bool base_ok = incomplete
            ? candidate.base_.build_lll_relation_base_impl(
                      bound, true)
            : candidate.base_.build_lll_relation_base(bound);
    if (!base_ok) {
        return false;
    }
    candidate.relations_ = RelationMatrix(candidate.base_);
    if (!candidate.relations_.is_defined() ||
        !reset_relation_modules(candidate.row_module_,
                                candidate.relation_module_,
                                candidate.base_.length(),
                                candidate.diagnostics_)) {
        return false;
    }

    candidate.has_base_ = true;
    if (!candidate.record_factor_base_generation_(bound)) {
        return false;
    }
    swap(candidate);
    return true;
}

bool ClassGroupContext::set_factor_base(const FactorBase& base_in) noexcept {
    if (!base_in.is_defined() || base_in.parent() == nullptr) {
        return false;
    }
    if (is_defined() && !same_order_parent(&parent_, base_in.parent())) {
        return false;
    }

    ClassGroupContext candidate(*base_in.parent());
    candidate.set_diagnostics(diagnostics_);
    if (!candidate.is_defined() ||
        !candidate.base_.set(base_in)) {
        return false;
    }
    candidate.relations_ = RelationMatrix(candidate.base_);
    if (!candidate.relations_.is_defined() ||
        !reset_relation_modules(candidate.row_module_,
                                candidate.relation_module_,
                                candidate.base_.length(),
                                candidate.diagnostics_)) {
        return false;
    }

    candidate.has_base_ = true;
    swap(candidate);
    return true;
}

bool ClassGroupContext::append_relation(const Relation& relation) noexcept {
    return append_relation(relation, ClassGroupRelationSource::Supplied);
}

bool ClassGroupContext::append_relation(
        const Relation& relation,
        ClassGroupRelationSource source) noexcept {
    RelationAppendOutcome outcome = RelationAppendOutcome::none;
    return append_relation_with_outcome_(outcome, relation, source);
}

bool ClassGroupContext::append_relation_with_outcome_(
        RelationAppendOutcome& outcome,
        const Relation& relation,
        ClassGroupRelationSource source) noexcept {
    return append_relation_with_outcome_(
            outcome, relation, source, DependentRelationPolicy::native);
}

bool ClassGroupContext::append_relation_with_outcome_(
        RelationAppendOutcome& outcome,
        const Relation& relation,
        ClassGroupRelationSource source,
        DependentRelationPolicy dependent_policy) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.append_relation");
    outcome = RelationAppendOutcome::none;
    if (!has_factor_base() || !same_factor_base(relation.factor_base(), &base_)) {
        return false;
    }
    if (!valid_relation_source(source)) {
        return false;
    }

    flint::FmpzMat row(1, generator_count());
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.append_relation_row_extract");
        if (!relation.exponents(flint::FmpzMatRef(row))) {
            return false;
        }
    }

    if (dependent_policy == DependentRelationPolicy::keep_nonduplicate) {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.append_relation_duplicate_check");
        bool duplicate = false;
        if (!relation_matrix_contains_row(duplicate, relations_,
                                          flint::FmpzMatConstRef(row))) {
            return false;
        }
        if (duplicate) {
            ++skipped_dependent_relations_;
            outcome = RelationAppendOutcome::skipped_dependent;
            SILEX_LOG(diagnostics_, DiagnosticsModule::class_group,
                      LogLevel::detail, "skipped duplicate relation");
            return true;
        }
    }

    detail::ClassRelationModuleContext module_candidate;
    module_candidate.set_diagnostics(diagnostics_);
    detail::ClassRelationModuleAddResult add_result;
    RowModuleAddResult exact_add_result;
    fmpz_smat::HnfContext exact_module_candidate;
    exact_module_candidate.set_diagnostics(diagnostics_);
    const slong stored_relation_index = relation_count();
    bool relation_module_changed = false;
    bool exact_module_updated = false;
    bool retain_modular_dependent = true;
    bool exact_native_decision = false;
    switch (dependent_policy) {
        case DependentRelationPolicy::keep:
        case DependentRelationPolicy::keep_nonduplicate:
            retain_modular_dependent = true;
            break;
        case DependentRelationPolicy::native:
            if (detail::ClassGroupRelationSearchAccess::
                        defer_native_goal_publication(*this)) {
                exact_native_decision =
                        relation_rank_ < generator_count() ||
                        detail::ClassGroupRelationSearchAccess::
                                        relation_kernel_row_count(*this) >=
                                relation_kernel_units_target_;
            } else {
                exact_native_decision =
                        !has_presentation() ||
                        relation_kernel_unit_count() >=
                                relation_kernel_units_target_;
            }
            retain_modular_dependent = !exact_native_decision;
            break;
        case DependentRelationPolicy::skip:
            retain_modular_dependent = false;
            break;
    }
    if (exact_native_decision) {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.append_relation_module_classify");
        if (!relation_module_.classify_fmpz_mat_row(
                    add_result, flint::FmpzMatConstRef(row), 0,
                    retain_modular_dependent)) {
            return false;
        }
    } else {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.append_relation_module_update");
        {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.append_relation_module_copy");
            if (!module_candidate.set(relation_module_)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.append_relation_module_add");
            if (!module_candidate.add_fmpz_mat_row(
                        add_result, flint::FmpzMatConstRef(row), 0,
                        stored_relation_index, retain_modular_dependent)) {
                return false;
            }
        }
    }

    if (add_result.duplicate) {
        const bool keep_duplicate_relation_witness =
                dependent_policy == DependentRelationPolicy::keep ||
                (dependent_policy == DependentRelationPolicy::native &&
                 retain_modular_dependent);
        if (keep_duplicate_relation_witness) {
            // Explicit matrix replay may retain duplicate rows even though
            // the reference-style module stores each distinct row once.  The
            // native post-presentation defer path also keeps duplicate row
            // witnesses while searching for relation-kernel units.
        } else {
            ++skipped_dependent_relations_;
            outcome = RelationAppendOutcome::skipped_dependent;
            SILEX_LOG(diagnostics_, DiagnosticsModule::class_group,
                      LogLevel::detail, "skipped duplicate relation");
            return true;
        }
    }

    if (!add_result.duplicate && exact_native_decision) {
        bool keep_relation = false;
        bool exact_remainder_checked = false;
        bool exact_checkpoint_synced = false;
        if (!add_result.modular_independent) {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.append_relation_exact_remainder");
            bool exact_zero = false;
            if (!sync_row_module_checkpoint_()) {
                return false;
            }
            exact_checkpoint_synced = true;
            if (!row_module_.fmpz_mat_row_reduces_to_zero(
                        exact_zero, flint::FmpzMatConstRef(row), 0)) {
                return false;
            }
            if (exact_zero) {
                exact_remainder_checked = true;
                keep_relation = false;
            }
        }
        if (!exact_remainder_checked) {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.append_relation_exact_decision");
            if ((!exact_checkpoint_synced && !sync_row_module_checkpoint_()) ||
                !exact_module_candidate.set(row_module_) ||
                !row_module_add_with_result(
                        exact_add_result, exact_module_candidate,
                        flint::FmpzMatConstRef(row), diagnostics_)) {
                return false;
            }
            keep_relation = should_keep_dependent_relation(
                    *this, exact_add_result, relation_kernel_units_target_);
            if (exact_add_result.rank_increased) {
                keep_relation = true;
            }
        }

        if (!keep_relation) {
            // The rejected row did not change the exact relation lattice, but
            // the decision warmed the HNF/index checkpoint.  Keep that cache.
            if (!exact_remainder_checked) {
                row_module_.swap(exact_module_candidate);
                row_module_.set_diagnostics(diagnostics_);
            }
            relation_rank_ = row_module_.rank();
            row_module_synced_ = true;
            ++skipped_dependent_relations_;
            outcome = RelationAppendOutcome::skipped_dependent;
            SILEX_LOG(diagnostics_, DiagnosticsModule::class_group,
                      LogLevel::detail, "skipped dependent relation");
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.append_relation_exact_decision_cached");
            return true;
        }

        {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.append_relation_module_update");
            {
                SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                    "class_group.append_relation_module_copy");
                if (!module_candidate.set(relation_module_)) {
                    return false;
                }
            }
            {
                SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                    "class_group.append_relation_module_add");
                if (!module_candidate.add_fmpz_mat_row(
                            add_result, flint::FmpzMatConstRef(row), 0,
                            stored_relation_index, retain_modular_dependent)) {
                    return false;
                }
            }
        }

        if (!add_result.retained) {
            if (!module_candidate.add_fmpz_mat_row(
                        add_result, flint::FmpzMatConstRef(row), 0,
                        stored_relation_index, true) ||
                !add_result.retained || add_result.duplicate) {
                return false;
            }
        }
        exact_module_updated = true;
    }

    if (!add_result.duplicate && !add_result.retained) {
        ++skipped_dependent_relations_;
        outcome = RelationAppendOutcome::skipped_dependent;
        SILEX_LOG(diagnostics_, DiagnosticsModule::class_group,
                  LogLevel::detail, "skipped dependent relation");
        return true;
    }

    relation_module_changed = add_result.retained && !add_result.duplicate;

    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.append_relation_store");
        if (!relations_.append(relation)) {
            return false;
        }
        relation_sources_.push_back(source);
        relation_basis_flags_.push_back(
                add_result.modular_independent ? 1 : 0);
    }

    if (relation_module_changed) {
        relation_module_.swap(module_candidate);
        relation_rank_ = relation_module_.rank();
        row_module_synced_ = false;
    }
    if (exact_module_updated) {
        row_module_.swap(exact_module_candidate);
        row_module_.set_diagnostics(diagnostics_);
        relation_rank_ = row_module_.rank();
        row_module_synced_ = true;
    }
    quotient_.clear();
    reset_certification_metadata_();
    if (add_result.modular_independent || exact_add_result.rank_increased) {
        private_storage_->partial_relations_nonproductive_streak = 0;
        outcome = RelationAppendOutcome::rank;
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.relation_rank_increased");
    } else if (exact_add_result.index_refined) {
        private_storage_->partial_relations_nonproductive_streak = 0;
        outcome = RelationAppendOutcome::index;
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.relation_index_refined");
    } else {
        outcome = RelationAppendOutcome::kernel;
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.relation_kernel_row");
    }
    SILEX_DEBUG_CHECK(
            diagnostics_, DiagnosticsModule::class_group, DebugLevel::cheap,
            "relation storage/source counts match",
            relations_.length() ==
                    static_cast<slong>(relation_sources_.size()));
    SILEX_DEBUG_CHECK(
            diagnostics_, DiagnosticsModule::class_group, DebugLevel::cheap,
            "relation storage/basis-flag counts match",
            relations_.length() ==
                    static_cast<slong>(relation_basis_flags_.size()));
    return true;
}

bool ClassGroupContext::try_append_generator_relation(
        const Element& generator,
        ClassGroupRelationSource source) noexcept {
    bool partial_throttle_exit = false;
    return try_append_generator_relation(partial_throttle_exit, generator,
                                         source);
}

bool ClassGroupContext::try_append_generator_relation(
        bool& partial_throttle_exit,
        const Element& generator,
        ClassGroupRelationSource source) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_append_generator_relation");
    partial_throttle_exit = false;
    const FactorBase* base = factor_base();
    if (base == nullptr || parent_.parent() == nullptr ||
        !generator.has_parent(*parent_.parent()) ||
        !valid_relation_source(source)) {
        return false;
    }

    if (!generator_relation_scratch_.is_defined() &&
        !generator_relation_scratch_.define(*base)) {
        return false;
    }
    Relation& relation = generator_relation_scratch_;

    if (!relation.set_generator(generator, diagnostics_)) {
        if (!private_storage_->use_partial_relations &&
            !private_storage_->partial_relations_configured &&
            source != ClassGroupRelationSource::Supplied) {
            detail::ClassGroupRelationOptions options;
            configure_partial_relations_(options);
        }
        return try_partial_relation_(partial_throttle_exit, generator);
    }

    bool zero_row = true;
    if (!relation_exponents_are_zero(zero_row, relation)) {
        return false;
    }
    // reference `class_group_add_relation(...; always=true)` keeps a smooth
    // zero-row relation as a relation generator.  These rows are unit
    // witnesses for the later class/unit relation-kernel extraction.
    const bool keep_zero_unit_relation =
            zero_row &&
            (generator_relation_policy_ ==
                     DependentRelationPolicy::keep_nonduplicate ||
             generator_relation_policy_ == DependentRelationPolicy::keep ||
             relation_kernel_target_open(*this,
                                         relation_kernel_units_target_));
    if (zero_row && !keep_zero_unit_relation) {
        return true;
    }

    RelationAppendOutcome outcome = RelationAppendOutcome::none;
    return append_relation_with_outcome_(outcome, relation, source,
                                         generator_relation_policy_);
}

bool ClassGroupContext::try_append_generator_relation_with_norm(
        bool& partial_throttle_exit,
        const Element& generator,
        flint::FmpqConstRef norm,
        ClassGroupRelationSource source) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_append_generator_relation");
    partial_throttle_exit = false;
    const FactorBase* base = factor_base();
    if (base == nullptr || parent_.parent() == nullptr ||
        !generator.has_parent(*parent_.parent()) ||
        !valid_relation_source(source)) {
        return false;
    }

    if (!generator_relation_scratch_.is_defined() &&
        !generator_relation_scratch_.define(*base)) {
        return false;
    }
    Relation& relation = generator_relation_scratch_;

    if (!relation.set_generator_with_norm(generator, norm, diagnostics_)) {
        if (!private_storage_->use_partial_relations &&
            !private_storage_->partial_relations_configured &&
            source != ClassGroupRelationSource::Supplied) {
            detail::ClassGroupRelationOptions options;
            configure_partial_relations_(options);
        }
        return try_partial_relation_(partial_throttle_exit, generator);
    }

    bool zero_row = true;
    if (!relation_exponents_are_zero(zero_row, relation)) {
        return false;
    }
    const bool keep_zero_unit_relation =
            zero_row &&
            (generator_relation_policy_ ==
                     DependentRelationPolicy::keep_nonduplicate ||
             generator_relation_policy_ == DependentRelationPolicy::keep ||
             relation_kernel_target_open(*this,
                                         relation_kernel_units_target_));
    if (zero_row && !keep_zero_unit_relation) {
        return true;
    }

    RelationAppendOutcome outcome = RelationAppendOutcome::none;
    return append_relation_with_outcome_(outcome, relation, source,
                                         generator_relation_policy_);
}

bool ClassGroupContext::try_append_integral_generator_relation_(
        bool& partial_throttle_exit,
        const Element& generator,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial,
        ClassGroupRelationSource source) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_append_generator_relation");
    partial_throttle_exit = false;
    const FactorBase* base = factor_base();
    if (base == nullptr || parent_.parent() == nullptr ||
        !generator.has_parent(*parent_.parent()) ||
        flint::fmpz_mat_nrows(integral_coordinates) != 1 ||
        flint::fmpz_mat_ncols(integral_coordinates) != parent_.degree() ||
        !valid_relation_source(source)) {
        return false;
    }

    if (!generator_relation_scratch_.is_defined() &&
        !generator_relation_scratch_.define(*base)) {
        return false;
    }
    Relation& relation = generator_relation_scratch_;

    if (!detail::set_relation_from_integral_coordinates_and_norm(
                relation, generator, integral_coordinates, norm,
                integral_coordinate_polynomial, diagnostics_)) {
        if (!private_storage_->use_partial_relations &&
            !private_storage_->partial_relations_configured &&
            source != ClassGroupRelationSource::Supplied) {
            detail::ClassGroupRelationOptions options;
            configure_partial_relations_(options);
        }
        return try_partial_relation_(
                partial_throttle_exit, generator, norm.raw(),
                integral_coordinates.raw(),
                integral_coordinate_polynomial == nullptr
                        ? nullptr
                        : integral_coordinate_polynomial->raw());
    }

    bool zero_row = true;
    if (!relation_exponents_are_zero(zero_row, relation)) {
        return false;
    }
    const bool keep_zero_unit_relation =
            zero_row &&
            (generator_relation_policy_ ==
                     DependentRelationPolicy::keep_nonduplicate ||
             generator_relation_policy_ == DependentRelationPolicy::keep ||
             relation_kernel_target_open(*this,
                                         relation_kernel_units_target_));
    if (zero_row && !keep_zero_unit_relation) {
        return true;
    }

    RelationAppendOutcome outcome = RelationAppendOutcome::none;
    return append_relation_with_outcome_(outcome, relation, source,
                                         generator_relation_policy_);
}

void ClassGroupContext::reset_partial_relations_() noexcept {
    if (!is_defined() || private_storage_ == nullptr) {
        return;
    }
    private_storage_->partial_relations.clear();
    private_storage_->direct_partial_residue_blocks.clear();
    private_storage_->partial_relations_max = 0;
    private_storage_->partial_relations_nonproductive_streak = 0;
    private_storage_->use_partial_relations = false;
    private_storage_->partial_relations_configured = false;
}

void ClassGroupContext::configure_partial_relations_(
        const detail::ClassGroupRelationOptions& options) noexcept {
    reset_partial_relations_();
    if (private_storage_ == nullptr) {
        return;
    }
    private_storage_->partial_relations_configured = true;
    if (!has_factor_base() || parent_.degree() <= 2 ||
        options.max_candidates <= 0 || base_.length() <= 0) {
        return;
    }

    slong cap = max_slong_value(WORD(16), 4 * base_.length());
    cap = min_slong_value(cap, options.max_candidates);
    if (cap <= 0) {
        return;
    }

    private_storage_->partial_relations_max = cap;
    private_storage_->use_partial_relations = true;
}

bool ClassGroupContext::try_partial_relation_(
        bool& partial_throttle_exit,
        const Element& generator,
        const fmpq* known_norm,
        const fmpz_mat_struct* known_integral_coordinates,
        const fmpz_poly_struct* known_integral_polynomial) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation");
    partial_throttle_exit = false;
    if (!has_factor_base() || private_storage_ == nullptr ||
        !private_storage_->use_partial_relations ||
        parent_.parent() == nullptr ||
        !generator.has_parent(*parent_.parent())) {
        return true;
    }
    detail::ClassGroupContextStorage& partial_state = *private_storage_;

    PrimeIdeal large_prime(parent_);
    flint::Fmpz rational_large_prime;
    flint::FmpzPoly residue_key;
    flint::FmpzMat partial_row(1, base_.length());
    bool have_large_prime = false;
    bool have_residue_key = false;
    bool have_partial_row = false;
    if (!large_prime.is_defined()) {
        return false;
    }

    auto match_or_store_large_prime = [&]() noexcept -> bool {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.partial_relation_match_or_store");
        if (!have_large_prime && !have_residue_key) {
            return true;
        }

        std::size_t match_index = partial_state.partial_relations.size();
        bool duplicate_match = false;
        {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.partial_relation_match_scan");
            for (std::size_t i = 0;
                 i < partial_state.partial_relations.size(); ++i) {
                detail::PartialRelationEntry& stored =
                        partial_state.partial_relations[i];
                const bool residue_match =
                        have_residue_key && stored.has_residue_key &&
                        ::fmpz_equal(stored.rational_prime.raw(),
                                     rational_large_prime.raw()) != 0 &&
                        ::fmpz_poly_equal(stored.residue_key.raw(),
                                          residue_key.raw()) != 0;
                const bool prime_match =
                        have_large_prime && stored.prime.equal(large_prime);
                if (!residue_match && !prime_match) {
                    continue;
                }
                if (stored.generator.equal(generator)) {
                    duplicate_match = true;
                    break;
                }
                match_index = i;
                break;
            }
        }
        if (duplicate_match) {
            return true;
        }

        if (match_index != partial_state.partial_relations.size()) {
            detail::PartialRelationEntry& stored =
                    partial_state.partial_relations[match_index];

            flint::FmpzMat matched_row(0, 0);
            bool have_matched_row = false;
            if (have_partial_row && stored.has_row) {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_match_known_row");
                matched_row = flint::FmpzMat(1, base_.length());
                ::fmpz_mat_sub(matched_row.raw(), partial_row.raw(),
                               stored.row.raw());
                bool matched_row_zero = false;
                if (!fmpz_mat_row_is_zero(
                            matched_row_zero,
                            flint::FmpzMatConstRef(matched_row))) {
                    return false;
                }
                if (matched_row_zero) {
                    return true;
                }
                have_matched_row = true;
            }

            Element inverse(*generator.parent());
            Element quotient(*generator.parent());
            // The caller's scratch already owns this factor base.  Reuse it
            // just as the C context reuses C->R for a large-prime match.
            Relation& relation = generator_relation_scratch_;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_match_quotient");
                if (!inverse.is_defined() || !quotient.is_defined() ||
                    !relation.is_defined() ||
                    !inverse.invert(stored.generator) ||
                    !quotient.multiply(generator, inverse)) {
                    return true;
                }
            }

            if (have_matched_row) {
                if (!detail::set_relation_from_known_row(
                            relation, base_, quotient,
                            flint::FmpzMatConstRef(matched_row))) {
                    return false;
                }
            } else {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_match_factor");
                if (!relation.set_generator(quotient, diagnostics_)) {
                    return true;
                }
            }

            if (!have_matched_row) {
                bool zero_row = true;
                if (!relation_exponents_are_zero(zero_row, relation)) {
                    return false;
                }
                if (zero_row) {
                    return true;
                }
            }

            RelationAppendOutcome outcome = RelationAppendOutcome::none;
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_match_append");
                if (!append_relation_with_outcome_(
                            outcome, relation,
                            ClassGroupRelationSource::LargePrimeMatch)) {
                    return false;
                }
            }

            {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_match_erase");
                partial_state.partial_relations.erase(
                        partial_state.partial_relations.begin() +
                        static_cast<std::ptrdiff_t>(match_index));
            }

            if (outcome == RelationAppendOutcome::skipped_dependent) {
                ++partial_state.partial_relations_nonproductive_streak;
                const slong threshold =
                        max_slong_value(WORD(8), 2 * base_.length());
                if (partial_state.partial_relations_nonproductive_streak >=
                    threshold) {
                    partial_state.use_partial_relations = false;
                    partial_state.partial_relations.clear();
                    partial_throttle_exit = true;
                }
            }
            return true;
        }

        {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.partial_relation_store");
            if (partial_state.partial_relations_max > 0 &&
                static_cast<slong>(partial_state.partial_relations.size()) >=
                        partial_state.partial_relations_max) {
                return true;
            }

            detail::PartialRelationEntry entry(parent_);
            if (!entry.prime.is_defined() || !entry.generator.is_defined() ||
                !entry.generator.set(generator)) {
                return false;
            }
            if (have_large_prime && !entry.prime.set(large_prime)) {
                return false;
            }
            if (have_residue_key) {
                ::fmpz_set(entry.rational_prime.raw(),
                           rational_large_prime.raw());
                ::fmpz_poly_set(entry.residue_key.raw(), residue_key.raw());
                entry.has_residue_key = true;
            }
            if (have_partial_row) {
                entry.row = flint::FmpzMat(1, base_.length());
                flint::fmpz_mat_set(flint::FmpzMatRef(entry.row),
                                    flint::FmpzMatConstRef(partial_row));
                entry.has_row = true;
            }
            partial_state.partial_relations.emplace_back(std::move(entry));
        }
        return true;
    };

    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.partial_relation_direct_large_prime");
        // reference's large-prime route keys partial relations by the leftover
        // rational prime and residue polynomial for nice defining polynomials.
        // If this source-shaped screen cannot identify a unique degree-one
        // prime, the exact C-compatible factorization path below remains
        // authoritative.
        DirectLargePrimeStatus status = direct_partial_relation_row(
                partial_row, rational_large_prime, residue_key, generator,
                base_, flint::FmpzConstRef(factor_base_build_bound_),
                known_norm, known_integral_coordinates,
                known_integral_polynomial,
                partial_state.direct_partial_residue_blocks, diagnostics_);
        if (status == DirectLargePrimeStatus::found) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.partial_relation_direct_large_prime_found");
            have_residue_key = true;
            have_partial_row = true;
            return match_or_store_large_prime();
        }
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::class_group,
                status == DirectLargePrimeStatus::no_candidate
                        ? "class_group.partial_relation_direct_large_prime_no_candidate"
                : status == DirectLargePrimeStatus::out_of_bound
                        ? "class_group.partial_relation_direct_large_prime_out_of_bound"
                        : "class_group.partial_relation_direct_large_prime_unsupported");
        if (status == DirectLargePrimeStatus::out_of_bound) {
            return true;
        }
    }

    if (!parent_.is_equation_order() && known_norm != nullptr &&
        known_integral_coordinates != nullptr) {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "class_group.partial_relation_retained_integral_factor");
        OrderElement order_generator(parent_);
        detail::OneLargePrimeFactorStatus factor_status =
                detail::OneLargePrimeFactorStatus::no_candidate;
        if (order_generator.is_defined() &&
            order_generator.set_coordinates(flint::FmpzMatConstRef(
                    known_integral_coordinates)) &&
            detail::order_element_factor_over_base_with_one_large_prime(
                    factor_status, flint::FmpzMatRef(partial_row), large_prime,
                    order_generator, flint::FmpqConstRef(known_norm), base_,
                    diagnostics_)) {
            if (factor_status ==
                detail::OneLargePrimeFactorStatus::no_candidate) {
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_retained_integral_no_candidate");
                return true;
            }

            flint::Fmpz large_norm;
            flint::Fmpz norm_bound;
            flint::fmpz_mul(flint::FmpzRef(norm_bound),
                            flint::FmpzConstRef(factor_base_build_bound_),
                            flint::FmpzConstRef(factor_base_build_bound_));
            if (!large_prime.norm(flint::FmpzRef(large_norm)) ||
                flint::fmpz_cmp(flint::FmpzConstRef(large_norm),
                                flint::FmpzConstRef(norm_bound)) > 0) {
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::class_group,
                        "class_group.partial_relation_retained_integral_out_of_bound");
                return true;
            }

            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.partial_relation_retained_integral_found");
            have_large_prime = true;
            have_partial_row = true;
            return match_or_store_large_prime();
        }
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::class_group,
                "class_group.partial_relation_retained_integral_fallback");
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.partial_relation_full_factor");
        OrderElement order_generator(parent_);
        Ideal ideal(parent_);
        IdealFactorization factorization(parent_);
        PrimeIdeal factor_prime(parent_);
        if (!order_generator.is_defined() || !ideal.is_defined() ||
            !factorization.is_defined() || !factor_prime.is_defined() ||
            !order_generator.set_element(generator)) {
            return true;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.partial_relation_full_set_principal");
            if (!ideal.set_principal(order_generator)) {
                return true;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.partial_relation_full_factor_ideal");
            if (!factorization.factor(ideal)) {
                return true;
            }
        }

        have_large_prime = false;
        have_partial_row = false;
        slong large_exp = 0;
        flint::fmpz_mat_zero(flint::FmpzMatRef(partial_row));
        for (slong i = 0; i < factorization.length(); ++i) {
            slong exponent = 0;
            if (!factorization.prime(factor_prime, i) ||
                !factorization.exponent(exponent, i)) {
                return false;
            }
            const slong base_index = base_.index(factor_prime);
            if (base_index >= 0) {
                ::fmpz_add_si(flint::fmpz_mat_entry(partial_row, 0, base_index)
                                      .raw(),
                              flint::fmpz_mat_entry(partial_row, 0, base_index)
                                      .raw(),
                              exponent);
                continue;
            }
            if (!have_large_prime) {
                if (!large_prime.set(factor_prime)) {
                    return false;
                }
                large_exp = exponent;
                have_large_prime = true;
            } else if (large_prime.equal(factor_prime)) {
                large_exp += exponent;
            } else {
                return true;
            }
        }

        if (!have_large_prime || large_exp != 1) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.partial_relation_full_factor_no_candidate");
            return true;
        }

        flint::Fmpz large_norm;
        flint::Fmpz norm_bound;
        flint::fmpz_mul(flint::FmpzRef(norm_bound),
                        flint::FmpzConstRef(factor_base_build_bound_),
                        flint::FmpzConstRef(factor_base_build_bound_));
        if (!large_prime.norm(flint::FmpzRef(large_norm)) ||
            flint::fmpz_cmp(flint::FmpzConstRef(large_norm),
                            flint::FmpzConstRef(norm_bound)) > 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.partial_relation_full_factor_large_prime_out_of_bound");
            return true;
        }
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::class_group,
                "class_group.partial_relation_full_factor_found");
        have_partial_row = true;
    }

    return match_or_store_large_prime();
}

slong ClassGroupContext::relation_count() const noexcept {
    return has_factor_base() ? relations_.length() : 0;
}

slong ClassGroupContext::relation_rank() const noexcept {
    return has_factor_base() ? relation_rank_ : 0;
}

slong ClassGroupContext::skipped_dependent_relation_count() const noexcept {
    return has_factor_base() ? skipped_dependent_relations_ : 0;
}

bool ClassGroupContext::relation_source(
        ClassGroupRelationSource& out,
        slong index) const noexcept {
    if (!has_factor_base() || index < 0 || index >= relation_count() ||
        relation_sources_.size() !=
                static_cast<std::size_t>(relation_count())) {
        return false;
    }

    out = relation_sources_[static_cast<std::size_t>(index)];
    return true;
}

slong ClassGroupContext::relation_source_count(
        ClassGroupRelationSource source) const noexcept {
    if (!has_factor_base() || !valid_relation_source(source) ||
        relation_sources_.size() !=
                static_cast<std::size_t>(relation_count())) {
        return 0;
    }

    slong count = 0;
    for (ClassGroupRelationSource stored : relation_sources_) {
        if (stored == source) {
            ++count;
        }
    }
    return count;
}

bool ClassGroupContext::relations(flint::FmpzMatRef out) const noexcept {
    return has_factor_base() && relations_.rows(out);
}

std::optional<flint::FmpzMat> ClassGroupContext::relations() const noexcept {
    if (!has_factor_base()) {
        return std::nullopt;
    }
    flint::FmpzMat out(relation_count(), generator_count());
    if (!relations(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ClassGroupContext::relation_generator(Element& out,
                                           slong index) const noexcept {
    return has_factor_base() && relations_.generator(out, index);
}

const Element* ClassGroupContext::relation_generator_at(
        slong index) const noexcept {
    return has_factor_base() ? relations_.generator_at(index) : nullptr;
}

bool ClassGroupContext::set_relation_matrix(
        const RelationMatrix& matrix) noexcept {
    const FactorBase* matrix_base = matrix.factor_base();
    const Order* matrix_parent = matrix.parent();
    if (!matrix.is_defined() || matrix_base == nullptr ||
        matrix_parent == nullptr ||
        (is_defined() && !same_order_parent(&parent_, matrix_parent))) {
        return false;
    }

    ClassGroupContext candidate(*matrix_parent);
    candidate.set_diagnostics(diagnostics_);
    if (!candidate.is_defined() ||
        !candidate.base_.set(*matrix_base)) {
        return false;
    }
    candidate.relations_ = RelationMatrix(candidate.base_);
    if (!candidate.relations_.is_defined() ||
        !reset_relation_modules(candidate.row_module_,
                                candidate.relation_module_,
                                candidate.base_.length(),
                                candidate.diagnostics_)) {
        return false;
    }
    candidate.has_base_ = true;

    Element generator(*matrix_parent->parent());
    for (slong i = 0; i < matrix.length(); ++i) {
        Relation relation(candidate.base_);
        RelationAppendOutcome outcome = RelationAppendOutcome::none;
        if (!generator.is_defined() ||
            !matrix.generator(generator, i) ||
            !relation.set_generator(generator, diagnostics_) ||
            !candidate.append_relation_with_outcome_(
                    outcome, relation, ClassGroupRelationSource::Supplied,
                    DependentRelationPolicy::keep)) {
            return false;
        }
    }

    if (!candidate.relations_.to_abelian_group(
                candidate.quotient_)) {
        return false;
    }

    swap(candidate);
    return true;
}

bool ClassGroupContext::publish_presentation() noexcept {
    if (!has_factor_base()) {
        return false;
    }
    if (!sync_row_module_checkpoint_()) {
        return false;
    }

    flint::FmpzMat relation_rows(relation_count(), generator_count());
    flint::FmpzMat hnf_basis(generator_count(), generator_count());
    FiniteAbelianGroup candidate;
    if (!relations_.rows(flint::FmpzMatRef(relation_rows)) ||
        !row_module_.get_hnf_rows(flint::FmpzMatRef(hnf_basis)) ||
        !candidate.set_relation_matrix_with_hnf_basis(
                flint::FmpzMatConstRef(relation_rows),
                flint::FmpzMatConstRef(hnf_basis))) {
        return false;
    }

    quotient_.swap(candidate);
    return true;
}

bool ClassGroupContext::has_presentation() const noexcept {
    return is_defined() && quotient_.is_defined();
}

CertificationMode ClassGroupContext::certification_status() const noexcept {
    return has_presentation() ? certification_
                              : CertificationMode::unknown;
}

ProofState ClassGroupContext::factor_base_generation_status() const noexcept {
    return has_presentation() ? factor_base_generation_status_
                              : ProofState::not_checked;
}

ProofState ClassGroupContext::factor_base_generation_checked_status()
        const noexcept {
    return has_presentation() ? factor_base_generation_checked_status_
                              : ProofState::not_checked;
}

bool ClassGroupContext::factor_base_build_bound(
        flint::FmpzRef out) const noexcept {
    if (!has_presentation() ||
        factor_base_generation_status_ == ProofState::not_checked ||
        flint::fmpz_is_zero(
                flint::FmpzConstRef(factor_base_build_bound_))) {
        return false;
    }

    flint::fmpz_set(out,
                    flint::FmpzConstRef(factor_base_build_bound_));
    return true;
}

bool ClassGroupContext::factor_base_generation_bound(
        flint::FmpzRef out) const noexcept {
    if (!has_presentation() ||
        flint::fmpz_is_zero(
                flint::FmpzConstRef(factor_base_generation_bound_))) {
        return false;
    }

    flint::fmpz_set(
            out, flint::FmpzConstRef(factor_base_generation_bound_));
    return true;
}

bool ClassGroupContext::factor_base_generation_checked_bound(
        flint::FmpzRef out) const noexcept {
    if (!has_presentation() ||
        factor_base_generation_checked_status_ ==
                ProofState::not_checked) {
        return false;
    }

    flint::fmpz_set(
            out,
            flint::FmpzConstRef(factor_base_generation_checked_bound_));
    return true;
}

slong ClassGroupContext::factor_base_generation_record_count()
        const noexcept {
    return has_presentation() && private_storage_ != nullptr
            ? static_cast<slong>(
                      private_storage_->factor_base_generation_records.size())
            : 0;
}

bool ClassGroupContext::factor_base_generation_record(
        flint::FmpzRef p,
        ProofState& status,
        slong index) const noexcept {
    if (!has_presentation() || index < 0 ||
        index >= factor_base_generation_record_count()) {
        return false;
    }

    const detail::FactorBaseGenerationRecord& record =
            private_storage_->factor_base_generation_records[
                    static_cast<std::size_t>(index)];
    flint::fmpz_set(p, flint::FmpzConstRef(record.p));
    status = record.status;
    return true;
}

std::optional<ClassGroupFactorBaseGenerationRecord>
ClassGroupContext::factor_base_generation_record(slong index) const noexcept {
    ClassGroupFactorBaseGenerationRecord record;
    if (!factor_base_generation_record(flint::FmpzRef(record.p),
                                       record.status, index)) {
        return std::nullopt;
    }
    return record;
}

bool ClassGroupContext::check_factor_base_generation_bound(
        flint::FmpzConstRef required_bound) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::class_group,
            "class_group.factor_base_generation.check");
    if (!has_factor_base() || flint::fmpz_sgn(required_bound) < 0) {
        return false;
    }

    if (factor_base_generation_status_ == ProofState::verified &&
        factor_base_generation_checked_status_ == ProofState::verified &&
        flint::fmpz_cmp(
                flint::FmpzConstRef(factor_base_generation_checked_bound_),
                required_bound) >= 0) {
        return true;
    }

    reset_factor_base_generation_check_();
    if (private_storage_ == nullptr) {
        return false;
    }
    if (factor_base_generation_status_ != ProofState::verified ||
        flint::fmpz_cmp(required_bound,
                        flint::FmpzConstRef(
                                factor_base_build_bound_)) > 0) {
        mark_factor_base_generation_check_(
                required_bound, ProofState::unavailable);
        return false;
    }

    flint::Fmpz p;
    flint::fmpz_set_si(flint::FmpzRef(p), 2);
    while (flint::fmpz_cmp(flint::FmpzConstRef(p), required_bound) <= 0) {
        bool ok = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.factor_base_generation.decompose_and_match");
            ok = factor_base_covers_prime_decomposition(
                    base_, flint::FmpzConstRef(p), required_bound);
        }
        bool recorded = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.factor_base_generation.record");
            recorded = append_factor_base_generation_record(
                    private_storage_->factor_base_generation_records,
                    flint::FmpzConstRef(p),
                    ok ? ProofState::verified
                       : ProofState::unavailable);
        }
        if (!recorded) {
            mark_factor_base_generation_check_(
                    required_bound, ProofState::unavailable);
            return false;
        }

        if (!ok) {
            mark_factor_base_generation_check_(
                    required_bound, ProofState::unavailable);
            return false;
        }

        ::fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    return mark_factor_base_generation_check_(
            required_bound, ProofState::verified);
}

ProofState ClassGroupContext::relation_saturation_status() const noexcept {
    return has_presentation() ? relation_saturation_status_
                              : ProofState::not_checked;
}

ProofState ClassGroupContext::analytic_class_regulator_status()
        const noexcept {
    return has_presentation() ? analytic_class_regulator_status_
                              : ProofState::not_checked;
}

ProofState ClassGroupContext::zeta_bf_proof_status() const noexcept {
    return has_presentation() ? zeta_bf_status_
                              : ProofState::not_checked;
}

bool ClassGroupContext::zeta_bf_proof_record(
        ulong& cutoff,
        ulong& max_cutoff,
        slong& requested_precision,
        slong& work_precision,
        flint::ArbRef error_bound) const noexcept {
    if (!has_presentation() ||
        zeta_bf_status_ != ProofState::verified) {
        return false;
    }

    cutoff = zeta_bf_cutoff_;
    max_cutoff = zeta_bf_max_cutoff_;
    requested_precision = zeta_bf_requested_precision_;
    work_precision = zeta_bf_work_precision_;
    flint::arb_set(error_bound,
                   flint::ArbConstRef(zeta_bf_error_bound_));
    return true;
}

std::optional<ClassGroupZetaBfProofRecord>
ClassGroupContext::zeta_bf_proof_record() const noexcept {
    ClassGroupZetaBfProofRecord record;
    if (!zeta_bf_proof_record(record.cutoff, record.max_cutoff,
                              record.requested_precision,
                              record.work_precision,
                              flint::ArbRef(record.error_bound))) {
        return std::nullopt;
    }
    return record;
}

slong ClassGroupContext::relation_saturation_record_count()
        const noexcept {
    return has_presentation() && private_storage_ != nullptr
            ? static_cast<slong>(
                      private_storage_->relation_saturation_records.size())
            : 0;
}

bool ClassGroupContext::relation_saturation_record(
        flint::FmpzRef ell,
        ProofState& status,
        slong index) const noexcept {
    if (!has_presentation() || index < 0 ||
        index >= relation_saturation_record_count()) {
        return false;
    }

    const detail::RelationSaturationRecord& record =
            private_storage_->relation_saturation_records[
                    static_cast<std::size_t>(index)];
    flint::fmpz_set(ell, flint::FmpzConstRef(record.ell));
    status = record.status;
    return true;
}

std::optional<ClassGroupRelationSaturationRecord>
ClassGroupContext::relation_saturation_record(slong index) const noexcept {
    ClassGroupRelationSaturationRecord record;
    if (!relation_saturation_record(flint::FmpzRef(record.ell),
                                    record.status, index)) {
        return std::nullopt;
    }
    return record;
}

ProofState ClassGroupContext::unit_proof_status() const noexcept {
    return has_presentation() ? unit_proof_status_
                              : ProofState::not_checked;
}

ProofState ClassGroupContext::regulator_proof_status() const noexcept {
    return has_presentation() ? regulator_proof_status_
                              : ProofState::not_checked;
}

bool ClassGroupContext::try_certify_trivial_quotient(
        CertificationMode requested) noexcept {
    if (requested != CertificationMode::proven || !has_presentation() ||
        !parent_.is_maximal() ||
        unit_proof_status_ != ProofState::verified ||
        regulator_proof_status_ != ProofState::verified ||
        factor_base_generation_status_ != ProofState::verified ||
        flint::fmpz_is_zero(
                flint::FmpzConstRef(factor_base_generation_bound_))) {
        return false;
    }

    flint::Fmpz class_order;
    if (!order(flint::FmpzRef(class_order)) ||
        !flint::fmpz_is_one(flint::FmpzConstRef(class_order))) {
        return false;
    }

    flint::Fmpz required_bound;
    flint::fmpz_set(
            flint::FmpzRef(required_bound),
            flint::FmpzConstRef(factor_base_generation_bound_));
    if (!check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound))) {
        return false;
    }

    private_storage_->relation_saturation_records.clear();
    relation_saturation_status_ = ProofState::verified;
    if (factor_base_generation_checked_status_ ==
                ProofState::verified &&
        unit_proof_status_ == ProofState::verified &&
        regulator_proof_status_ == ProofState::verified) {
        certification_ = CertificationMode::proven;
        return true;
    }

    return false;
}

bool ClassGroupContext::try_certify_quadratic(
        CertificationMode requested) noexcept {
    if (requested != CertificationMode::proven || !has_presentation() ||
        !parent_.is_maximal() ||
        factor_base_generation_status_ != ProofState::verified) {
        return false;
    }

    flint::Fmpz discriminant;
    if (!detail::order_supports_exact_quadratic_class_certificate(parent_) ||
        !parent_.discriminant(flint::FmpzRef(discriminant))) {
        return false;
    }

    if (flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) < 0) {
        flint::Fmpz exact_order;
        return detail::ClassGroupCertificationAccess::
                exact_imaginary_quadratic_class_order_for_run(
                        flint::FmpzRef(exact_order), *this,
                        flint::FmpzConstRef(discriminant)) &&
               detail::ClassGroupCertificationAccess::
                       try_certify_imaginary_quadratic_from_exact_order(
                               *this, requested,
                               flint::FmpzConstRef(discriminant),
                               flint::FmpzConstRef(exact_order));
    }

    flint::Fmpz required_bound;
    flint::Fmpz class_order;
    if (!factor_base_generation_bound(flint::FmpzRef(required_bound)) ||
        !check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound)) ||
        !order(flint::FmpzRef(class_order)) ||
        !flint::fmpz_is_one(flint::FmpzConstRef(class_order))) {
        return false;
    }

    certification_ = CertificationMode::proven;
    relation_saturation_status_ = ProofState::verified;
    unit_proof_status_ = ProofState::verified;
    regulator_proof_status_ = ProofState::verified;
    return true;
}

bool ClassGroupContext::try_certify_with_units(
        const OrderUnitGroup& units,
        CertificationMode requested,
        slong precision) noexcept {
    return try_certify_with_units(units, requested, precision, 0);
}

bool ClassGroupContext::try_certify_with_units(
        const OrderUnitGroup& units,
        CertificationMode requested,
        slong precision,
        ulong zeta_bf_max_cutoff) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_certify_with_units");
    if (!detail::valid_certification_request(requested) || precision <= 0) {
        return false;
    }

    if (!has_presentation() || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_)) {
        return requested != CertificationMode::proven;
    }

    if (units.certification_status() != CertificationMode::proven) {
        return requested != CertificationMode::proven;
    }

    unit_proof_status_ = ProofState::verified;
    regulator_proof_status_ = ProofState::verified;

    if (requested != CertificationMode::proven) {
        return true;
    }

    if (try_promote_proven_certification_()) {
        return true;
    }

    if (try_certify_trivial_quotient(requested)) {
        return true;
    }

    flint::Arb analytic_hR;
    flint::Fmpz aux_bound;
    flint::fmpz_set_si(flint::FmpzRef(aux_bound), kAnalyticProofAuxPrimeBound);
    if (zeta_bf_max_cutoff != 0 && parent_.degree() > 2) {
        slong rank = -1;
        if (relation_saturation_status_ != ProofState::verified &&
            parent_.parent() != nullptr &&
            unit_rank(rank, *parent_.parent()) && rank == 1 &&
            zeta_class_regulator_product(flint::ArbRef(analytic_hR),
                                          parent_, precision)) {
            SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.pre_bf_analytic_index_bound");
            (void) try_analytic_index_bound_with_units(
                    units, flint::ArbConstRef(analytic_hR),
                    flint::FmpzConstRef(aux_bound), precision);
        }

        flint::Arb error_bound;
        ulong cutoff = 0;
        slong work_precision = 0;
        {
            SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                                "class_group.zeta_bf_audit");
            if (!zeta_class_regulator_product_bf_audit(
                        flint::ArbRef(analytic_hR),
                        flint::ArbRef(error_bound), cutoff, work_precision,
                        parent_, zeta_bf_max_cutoff, precision)) {
                return false;
            }
        }

        if (try_certify_analytic_class_regulator_(
                    units, flint::ArbConstRef(analytic_hR), precision)) {
            return record_zeta_bf_audit_(
                    flint::ArbConstRef(error_bound), cutoff,
                    zeta_bf_max_cutoff, precision, work_precision);
        }
        if (try_analytic_index_bound_with_units(
                    units, flint::ArbConstRef(analytic_hR),
                    flint::FmpzConstRef(aux_bound), precision)) {
            return record_zeta_bf_audit_(
                    flint::ArbConstRef(error_bound), cutoff,
                    zeta_bf_max_cutoff, precision, work_precision);
        }
        return false;
    }

    if (zeta_class_regulator_product(flint::ArbRef(analytic_hR),
                                      parent_, precision)) {
        if (try_certify_analytic_class_regulator_(
                    units, flint::ArbConstRef(analytic_hR), precision)) {
            return true;
        }
        if (try_analytic_index_bound_with_units(
                    units, flint::ArbConstRef(analytic_hR),
                    flint::FmpzConstRef(aux_bound), precision)) {
            return true;
        }
    }

    return false;
}

bool ClassGroupContext::prove_relation_saturation_dlog_ell_(
        const OrderUnitGroup& units,
        flint::FmpzConstRef ell,
        flint::FmpzConstRef aux_prime_bound) noexcept {
    if (!has_presentation() ||
        !flint::fmpz_is_prime(ell) ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 ||
        !units.is_set() || !same_order_parent(units.parent(), &parent_) ||
        units.certification_status() != CertificationMode::proven ||
        !relation_saturation_proof_prereqs_verified_()) {
        return false;
    }

    std::vector<FactoredElement> beta_rows;
    bool include_torsion = false;
    if (!relation_saturation_beta_rows(beta_rows, *this, ell) ||
        !relation_saturation_torsion_needed(include_torsion, units, ell)) {
        return false;
    }

    const slong target_rank =
            static_cast<slong>(beta_rows.size()) + units.free_rank() +
            (include_torsion ? 1 : 0);
    if (target_rank <= 0) {
        return false;
    }

    std::vector<flint::FmpzMat> columns;
    flint::FmpzMat current_dlog(target_rank, 0);
    slong best_rank = 0;

    flint::Fmpz p;
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (flint::fmpz_cmp(flint::FmpzConstRef(p), aux_prime_bound) <= 0) {
        PrimeIdealList local;
        if (decompose_prime(local, parent_, flint::FmpzConstRef(p))) {
            for (slong i = 0; i < local.size(); ++i) {
                const PrimeIdeal* prime = local.at(i);
                if (prime == nullptr) {
                    continue;
                }

                flint::FmpzMat column(target_rank, 1);
                if (!relation_saturation_proof_prime_dlog_column(
                            column, beta_rows, *this, units, include_torsion,
                            *prime, ell)) {
                    continue;
                }

                flint::FmpzMat candidate_dlog(target_rank, 0);
                if (!append_dlog_column(candidate_dlog, columns, column,
                                        target_rank)) {
                    return false;
                }

                slong rank = 0;
                (void) verify_dlog_full_rank_mod_ell(
                        rank, flint::FmpzMatConstRef(candidate_dlog), ell,
                        target_rank);
                best_rank = rank;
                current_dlog = std::move(candidate_dlog);
                columns.push_back(std::move(column));

                if (rank >= target_rank) {
                    if (!record_relation_saturation_proof_(
                                ell, ProofState::verified, rank, target_rank,
                                static_cast<slong>(columns.size()))) {
                        return false;
                    }
                    return true;
                }
            }
        }
        flint::fmpz_nextprime(flint::FmpzRef(p), flint::FmpzConstRef(p));
    }

    if (!columns.empty()) {
        (void) verify_dlog_full_rank_mod_ell(
                best_rank, flint::FmpzMatConstRef(current_dlog), ell,
                target_rank);
    }

    if (!record_relation_saturation_proof_(
                ell, ProofState::unavailable, best_rank, target_rank,
                static_cast<slong>(columns.size()))) {
        return false;
    }
    return true;
}

bool ClassGroupContext::try_promote_proven_certification_() noexcept {
    const bool unit_regulator_verified =
            unit_proof_status_ == ProofState::verified &&
            regulator_proof_status_ == ProofState::verified;
    const bool factor_base_verified =
            factor_base_generation_status_ == ProofState::verified &&
            factor_base_generation_checked_status_ == ProofState::verified;
    const bool factor_base_backed_proof =
            factor_base_verified && unit_regulator_verified &&
            (analytic_class_regulator_status_ == ProofState::verified ||
             quadratic_completeness_verified_() ||
             relation_saturation_proof_complete_());
    if (!has_presentation() || !factor_base_backed_proof) {
        return false;
    }

    certification_ = CertificationMode::proven;
    return true;
}

bool ClassGroupContext::try_certify_analytic_class_regulator_(
        const OrderUnitGroup& units,
        flint::ArbConstRef analytic_hR,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_certify_analytic_class_regulator");
    if (!has_presentation() ||
        !parent_.is_maximal() || precision <= 0 || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_) ||
        units.certification_status() != CertificationMode::proven) {
        return false;
    }

    flint::Fmpz required_bound;
    if (!factor_base_generation_bound(flint::FmpzRef(required_bound)) ||
        !check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound))) {
        return false;
    }

    slong rank = -1;
    if (parent_.parent() == nullptr ||
        !unit_rank(rank, *parent_.parent()) ||
        units.free_rank() != rank ||
        !class_regulator_index_is_one(units, *this, analytic_hR,
                                      precision)) {
        return false;
    }

    analytic_class_regulator_status_ = ProofState::verified;
    unit_proof_status_ = ProofState::verified;
    regulator_proof_status_ = ProofState::verified;
    return try_promote_proven_certification_();
}

bool ClassGroupContext::try_certify_analytic_class_unit_regulator_(
        OrderUnitGroup& units,
        flint::ArbConstRef analytic_hR,
        slong precision) noexcept {
    if (!record_analytic_class_unit_regulator_(units, analytic_hR,
                                               precision)) {
        return false;
    }

    return try_promote_proven_certification_();
}

bool ClassGroupContext::record_analytic_class_unit_regulator_(
        OrderUnitGroup& units,
        flint::ArbConstRef analytic_hR,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::class_group,
            "class_group.record_analytic_class_unit_regulator");
    if (!has_presentation() ||
        !parent_.is_maximal() || precision <= 0 || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_)) {
        return false;
    }

    flint::Fmpz required_bound;
    slong rank = -1;
    if (!factor_base_generation_bound(flint::FmpzRef(required_bound)) ||
        !check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound)) ||
        parent_.parent() == nullptr ||
        !unit_rank(rank, *parent_.parent()) ||
        units.free_rank() != rank ||
        !class_regulator_index_is_one(units, *this, analytic_hR,
                                      precision)) {
        return false;
    }

    units.mark_certification_proven_();
    analytic_class_regulator_status_ = ProofState::verified;
    unit_proof_status_ = ProofState::verified;
    regulator_proof_status_ = ProofState::verified;
    return true;
}

bool ClassGroupContext::try_certify_class_unit_with_units(
        OrderUnitGroup& units,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision) noexcept {
    if (!flint::arb_is_finite(analytic_class_regulator_product) ||
        !flint::arb_is_positive(analytic_class_regulator_product)) {
        return false;
    }
    return try_certify_analytic_class_unit_regulator_(
            units, analytic_class_regulator_product, precision);
}

bool ClassGroupContext::try_certify_class_unit_with_zeta(
        OrderUnitGroup& units,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_certify_class_unit_with_zeta");
    if (!has_presentation() || precision <= 0) {
        return false;
    }

    flint::Arb analytic_hR;
    if (!zeta_class_regulator_product(flint::ArbRef(analytic_hR),
                                      parent_, precision)) {
        return false;
    }
    return try_certify_analytic_class_unit_regulator_(
            units, flint::ArbConstRef(analytic_hR), precision);
}

bool ClassGroupContext::try_certify_class_unit_with_zeta_bf(
        OrderUnitGroup& units,
        ulong max_cutoff,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_certify_class_unit_with_zeta_bf");
    if (!has_presentation() || precision <= 0 ||
        max_cutoff == 0) {
        return false;
    }

    flint::Arb analytic_hR;
    flint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.zeta_bf_audit");
        if (!zeta_class_regulator_product_bf_audit(
                    flint::ArbRef(analytic_hR), flint::ArbRef(error_bound),
                    cutoff, work_precision, parent_, max_cutoff,
                    precision)) {
            return false;
        }
    }

    flint::Fmpz required_bound;
    if (!factor_base_generation_bound(flint::FmpzRef(required_bound)) ||
        !check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound))) {
        return false;
    }

    if (!record_analytic_class_unit_regulator_(
                units, flint::ArbConstRef(analytic_hR), precision)) {
        return false;
    }

    return record_zeta_bf_audit_(flint::ArbConstRef(error_bound), cutoff,
                                 max_cutoff, precision, work_precision) &&
           try_promote_proven_certification_();
}

bool ClassGroupContext::record_zeta_bf_audit_(
        flint::ArbConstRef error_bound,
        ulong cutoff,
        ulong max_cutoff,
        slong requested_precision,
        slong work_precision) noexcept {
    if (!has_presentation() || requested_precision <= 0 ||
        work_precision <= 0 || !flint::arb_is_finite(error_bound)) {
        return false;
    }

    analytic_class_regulator_status_ = ProofState::verified;
    zeta_bf_status_ = ProofState::verified;
    zeta_bf_cutoff_ = cutoff;
    zeta_bf_max_cutoff_ = max_cutoff;
    zeta_bf_requested_precision_ = requested_precision;
    zeta_bf_work_precision_ = work_precision;
    flint::arb_set(flint::ArbRef(zeta_bf_error_bound_), error_bound);
    return true;
}

bool ClassGroupContext::try_prove_relation_saturation_with_units(
        const OrderUnitGroup& units,
        flint::FmpzConstRef ell,
        flint::FmpzConstRef aux_prime_bound) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::class_group,
            "class_group.try_prove_relation_saturation_with_units");
    if (!prove_relation_saturation_dlog_ell_(units, ell, aux_prime_bound) ||
        !complete_relation_saturation_proof_(ell)) {
        return false;
    }
    return relation_saturation_status_ == ProofState::verified;
}

bool ClassGroupContext::try_prove_relation_saturation_index_bound_with_units(
        const OrderUnitGroup& units,
        flint::FmpzConstRef index_bound,
        flint::FmpzConstRef aux_prime_bound) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::class_group,
            "class_group.try_prove_relation_saturation_index_bound");
    if (!has_presentation() ||
        flint::fmpz_cmp_ui(index_bound, 1) <= 0 ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 ||
        !units.is_set() || !same_order_parent(units.parent(), &parent_) ||
        units.certification_status() != CertificationMode::proven ||
        !relation_saturation_proof_prereqs_verified_()) {
        return false;
    }

    std::vector<flint::Fmpz> required_ells;
    flint::Fmpz ell;
    flint::fmpz_one(flint::FmpzRef(ell));
    flint::fmpz_nextprime(flint::FmpzRef(ell), flint::FmpzConstRef(ell),
                          true);
    while (flint::fmpz_cmp(flint::FmpzConstRef(ell), index_bound) <= 0) {
        required_ells.emplace_back();
        flint::fmpz_set(flint::FmpzRef(required_ells.back()),
                        flint::FmpzConstRef(ell));
        (void) prove_relation_saturation_dlog_ell_(
                units, flint::FmpzConstRef(required_ells.back()),
                aux_prime_bound);
        flint::fmpz_nextprime(flint::FmpzRef(ell),
                              flint::FmpzConstRef(ell), true);
    }

    if (!complete_relation_saturation_proof_(required_ells) ||
        relation_saturation_status_ != ProofState::verified) {
        return false;
    }

    return try_promote_proven_certification_();
}

bool ClassGroupContext::try_analytic_index_bound_with_units(
        const OrderUnitGroup& units,
        flint::ArbConstRef analytic_class_regulator_product,
        flint::FmpzConstRef aux_prime_bound,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_analytic_index_bound_with_units");
    if (!has_presentation() ||
        !parent_.is_maximal() || precision <= 0 ||
        !flint::arb_is_finite(analytic_class_regulator_product) ||
        !flint::arb_is_positive(analytic_class_regulator_product) ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 ||
        factor_base_generation_status_ != ProofState::verified ||
        flint::fmpz_is_zero(
                flint::FmpzConstRef(factor_base_generation_bound_)) ||
        !units.is_set() || !same_order_parent(units.parent(), &parent_) ||
        units.certification_status() != CertificationMode::proven) {
        return false;
    }

    flint::Fmpz required_bound;
    if (!factor_base_generation_bound(flint::FmpzRef(required_bound)) ||
        !check_factor_base_generation_bound(
                flint::FmpzConstRef(required_bound))) {
        return false;
    }

    flint::Arb candidate_hR;
    if (!units.class_regulator_product(
                flint::ArbRef(candidate_hR), *this, precision)) {
        return false;
    }
    const bool index_is_one =
            detail::class_regulator_index_is_one_from_candidate_product(
                    flint::ArbConstRef(candidate_hR),
                    analytic_class_regulator_product, precision,
                    diagnostics_);

    unit_proof_status_ = ProofState::verified;
    regulator_proof_status_ = ProofState::verified;
    if (!index_is_one) {
        return false;
    }

    private_storage_->relation_saturation_records.clear();
    relation_saturation_status_ = ProofState::verified;
    analytic_class_regulator_status_ = ProofState::verified;
    return try_promote_proven_certification_();
}

bool ClassGroupContext::relation_row_refines_(
        const Relation& relation) noexcept {
    if (!has_factor_base() || !same_factor_base(relation.factor_base(), &base_)) {
        return false;
    }
    if (!sync_row_module_checkpoint_()) {
        return false;
    }

    flint::Fmpz index;
    if (!row_module_.full_rank_index(flint::FmpzRef(index))) {
        return true;
    }

    flint::FmpzMat row(1, generator_count());
    if (!relation.exponents(flint::FmpzMatRef(row))) {
        return false;
    }

    fmpz_smat::HnfContext module_candidate;
    module_candidate.set_diagnostics(diagnostics_);
    RowModuleAddResult add_result;
    if (!module_candidate.set(row_module_) ||
        !row_module_add_with_result(add_result, module_candidate,
                                    flint::FmpzMatConstRef(row),
                                    diagnostics_)) {
        return false;
    }
    return !add_result.rank_increased && add_result.index_refined;
}

bool ClassGroupContext::saturate_local_once_(
        bool& changed,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturate_local_once");
    changed = false;
    if (!has_factor_base() || !same_order_parent(prime.parent(), &parent_) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::FmpzMat kernel(0, 0);
    SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturation_residue_dlog");
    if (!residue_dlog_kernel(kernel, *this, prime, ell)) {
        return false;
    }

    bool saw_failure = false;
    const FactorBase* base = factor_base();
    if (base == nullptr) {
        return false;
    }

    for (slong i = 0; i < flint::fmpz_mat_nrows(kernel); ++i) {
        if (kernel_row_divisible(flint::FmpzMatConstRef(kernel), i, ell)) {
            continue;
        }

        Relation relation(*base);
        bool is_relation = false;
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                            "class_group.saturation_relation_candidate");
        if (!relation.is_defined() ||
            !saturation_relation(relation, is_relation, *this,
                                 flint::FmpzMatConstRef(kernel), i, ell)) {
            saw_failure = true;
            continue;
        }
        if (!is_relation) {
            continue;
        }
        if (!relation_row_refines_(relation)) {
            continue;
        }
        if (!append_relation(relation, ClassGroupRelationSource::Saturation)) {
            return false;
        }
        SILEX_LOG(diagnostics_, DiagnosticsModule::class_group,
                  LogLevel::detail, "saturation appended relation");
        changed = true;
        return true;
    }

    return !saw_failure;
}

bool ClassGroupContext::saturate_local_(
        bool& changed,
        bool& index_cleared,
        flint::FmpzConstRef ell,
        flint::FmpzConstRef aux_prime_bound,
        slong max_appends) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturate_local");
    changed = false;
    index_cleared = false;
    if (!has_factor_base() || !flint::fmpz_is_prime(ell) ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 || max_appends <= 0) {
        return false;
    }

    flint::Fmpz index;
    if (!sync_row_module_checkpoint_() ||
        !row_module_.full_rank_index(flint::FmpzRef(index))) {
        return false;
    }

    bool cleared = !flint::fmpz_divisible(flint::FmpzConstRef(index), ell);
    bool any_changed = false;
    slong appends = 0;
    const bool factor_base_covers_aux_primes =
            has_factor_base() &&
            flint::fmpz_cmp(flint::FmpzConstRef(factor_base_build_bound_),
                            aux_prime_bound) >= 0;
    while (!cleared && appends < max_appends) {
        bool appended = false;
        if (factor_base_covers_aux_primes) {
            flint::Fmpz q;
            for (slong block = 0;
                 block < base_.rational_prime_block_count(); ++block) {
                slong block_length = 0;
                if (!base_.rational_prime_block_data(flint::FmpzRef(q),
                                                      block_length, block)) {
                    return false;
                }
                if (flint::fmpz_cmp(flint::FmpzConstRef(q),
                                    aux_prime_bound) > 0) {
                    break;
                }
                for (slong offset = 0; offset < block_length; ++offset) {
                    slong prime_index = -1;
                    if (!base_.rational_prime_block_index(
                                prime_index, block, offset)) {
                        return false;
                    }
                    const PrimeIdeal* prime = base_.prime_at(prime_index);
                    if (prime == nullptr ||
                        !auxiliary_prime_usable(*prime,
                                                flint::FmpzConstRef(q), ell)) {
                        continue;
                    }

                    bool local_changed = false;
                    if (saturate_local_once_(local_changed, *prime, ell) &&
                        local_changed) {
                        appended = true;
                        break;
                    }
                }
                if (appended) {
                    break;
                }
            }
        } else {
            flint::Fmpz q;
            flint::fmpz_set_ui(flint::FmpzRef(q), 2);
            while (flint::fmpz_cmp(flint::FmpzConstRef(q),
                                   aux_prime_bound) <= 0) {
                PrimeIdealList local;
                if (decompose_prime(local, parent_,
                                    flint::FmpzConstRef(q))) {
                    for (slong i = 0; i < local.size(); ++i) {
                        const PrimeIdeal* prime = local.at(i);
                        if (prime == nullptr ||
                            !auxiliary_prime_usable(
                                    *prime, flint::FmpzConstRef(q), ell)) {
                            continue;
                        }

                        bool local_changed = false;
                        if (saturate_local_once_(local_changed, *prime, ell) &&
                            local_changed) {
                            appended = true;
                            break;
                        }
                    }
                }
                if (appended) {
                    break;
                }
                flint::fmpz_nextprime(flint::FmpzRef(q),
                                      flint::FmpzConstRef(q));
            }
        }

        if (!appended) {
            break;
        }

        any_changed = true;
        ++appends;
        if (!sync_row_module_checkpoint_() ||
            !row_module_.full_rank_index(flint::FmpzRef(index))) {
            return false;
        }
        cleared = !flint::fmpz_divisible(flint::FmpzConstRef(index), ell);
    }

    changed = any_changed;
    index_cleared = cleared;
    return true;
}

bool ClassGroupContext::saturate_relations_bounded_(
        bool& changed,
        bool& saturated,
        flint::FmpzConstRef aux_prime_bound,
        slong max_appends_per_ell,
        slong max_appends_total) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturate_relations_bounded");
    changed = false;
    saturated = false;
    if (!has_factor_base() || !has_presentation() ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 ||
        max_appends_per_ell <= 0 || max_appends_total <= 0) {
        return false;
    }

    flint::Fmpz index;
    if (!sync_row_module_checkpoint_() ||
        !row_module_.full_rank_index(flint::FmpzRef(index))) {
        return false;
    }

    flint::FmpzFactor factorization;
    flint::fmpz_factor(flint::FmpzFactorRef(factorization),
                       flint::FmpzConstRef(index));
    if (flint::fmpz_factor_num(flint::FmpzFactorConstRef(factorization)) ==
        0) {
        saturated = true;
        return true;
    }

    private_storage_->relation_saturation_records.clear();
    bool any_changed = false;
    bool all_cleared = true;
    slong remaining = max_appends_total;
    for (slong i = 0;
         i < flint::fmpz_factor_num(flint::FmpzFactorConstRef(factorization));
         ++i) {
        flint::Fmpz ell;
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(ell),
                flint::FmpzFactorConstRef(factorization), i);

        bool local_changed = false;
        bool index_cleared = false;
        if (remaining > 0) {
            const slong attempt_limit =
                    min_slong_value(max_appends_per_ell, remaining);
            const slong before = relation_count();
            if (!saturate_local_(local_changed, index_cleared,
                                 flint::FmpzConstRef(ell), aux_prime_bound,
                                 attempt_limit)) {
                return false;
            }
            remaining -= relation_count() - before;
        }

        any_changed = any_changed || local_changed;
        all_cleared = all_cleared && index_cleared;
        if (!mark_relation_saturation_(
                    flint::FmpzConstRef(ell),
                    index_cleared ? ProofState::verified
                                  : ProofState::unavailable)) {
            return false;
        }
    }

    if (any_changed && !publish_presentation()) {
        return false;
    }
    changed = any_changed;
    saturated = all_cleared;
    return true;
}

bool ClassGroupContext::try_auto_relation_saturation_(
        const Order& order,
        flint::FmpzConstRef factor_base_bound,
        const detail::ClassGroupRelationOptions& options,
        bool* changed_out) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.try_auto_relation_saturation");
    if (changed_out != nullptr) {
        *changed_out = false;
    }
    if (options.requested_certification != CertificationMode::unknown ||
        relation_saturation_option_enabled(options) || order.degree() <= 2 ||
        !has_presentation() || relation_count() >= options.max_relations) {
        return true;
    }

    if (class_unit_transaction_context_ != nullptr &&
        class_unit_transaction_context_
                ->defer_relation_saturation_until_units &&
        options.target_relation_kernel_units > 0) {
        // In the combined class/unit route, reference validates and saturates the
        // class/unit pair after relation-kernel units are available.  Avoid
        // doing class-only saturation before that validation step.
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::class_group,
                "class_group.auto_relation_saturation.skip_target_units");
        return true;
    }

    if (options.target_relation_kernel_units == 0 &&
        relation_kernel_unit_count() > 2) {
        return true;
    }

    flint::Fmpz class_order;
    if (!this->order(flint::FmpzRef(class_order)) ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(class_order), 1) <= 0) {
        return true;
    }

    const slong remaining = options.max_relations - relation_count();
    const slong total_cap = min_slong_value(remaining, WORD(4));
    if (total_cap <= 0) {
        return true;
    }

    flint::Fmpz aux_bound;
    flint::fmpz_set(flint::FmpzRef(aux_bound), factor_base_bound);
    if (flint::fmpz_cmp_ui(flint::FmpzConstRef(aux_bound), 32) < 0) {
        flint::fmpz_set_ui(flint::FmpzRef(aux_bound), 32);
    }

    bool changed = false;
    bool saturated = false;
    const bool ok = saturate_relations_bounded_(
            changed, saturated, flint::FmpzConstRef(aux_bound), 2,
            total_cap);
    if (ok && changed_out != nullptr) {
        *changed_out = changed;
    }
    return ok;
}

bool ClassGroupContext::saturate_local_once_with_units_(
        bool& changed,
        const OrderUnitGroup& units,
        bool include_torsion,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturate_local_once_with_units");
    changed = false;
    if (!has_factor_base() || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_) ||
        !same_order_parent(prime.parent(), &parent_) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    bool torsion_needed = false;
    if (!relation_saturation_torsion_needed(torsion_needed, units, ell) ||
        include_torsion != torsion_needed) {
        return false;
    }

    flint::FmpzMat kernel(0, 0);
    SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturation_residue_dlog_with_units");
    if (!residue_dlog_kernel_with_units(kernel, *this, units,
                                        include_torsion, prime, ell)) {
        return false;
    }

    bool saw_failure = false;
    const FactorBase* base = factor_base();
    if (base == nullptr) {
        return false;
    }

    flint::FmpzMat relation_rows(0, 0);
    bool have_relation_rows = false;
    Relation relation(*base);
    if (!relation.is_defined()) {
        return false;
    }
    for (slong i = 0; i < flint::fmpz_mat_nrows(kernel); ++i) {
        if (kernel_row_divisible(flint::FmpzMatConstRef(kernel), i, ell)) {
            continue;
        }

        if (!have_relation_rows) {
            relation_rows = flint::FmpzMat(relation_count(),
                                           generator_count());
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "class_group.saturation_relation_rows_with_units");
            if (!relations(flint::FmpzMatRef(relation_rows))) {
                return false;
            }
            have_relation_rows = true;
        }

        bool is_relation = false;
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::class_group,
                "class_group.saturation_relation_candidate_with_units");
        if (!saturation_relation_with_units(
                    relation, is_relation, *this,
                    flint::FmpzMatConstRef(relation_rows),
                    flint::FmpzMatConstRef(kernel), i, units,
                    include_torsion, ell)) {
            saw_failure = true;
            continue;
        }
        if (!is_relation) {
            continue;
        }
        if (!relation_row_refines_(relation)) {
            continue;
        }
        if (!append_relation(relation, ClassGroupRelationSource::Saturation)) {
            return false;
        }
        SILEX_LOG(diagnostics_, DiagnosticsModule::class_group,
                  LogLevel::detail,
                  "saturation with units appended relation");
        changed = true;
        return true;
    }

    return !saw_failure;
}

bool ClassGroupContext::saturate_local_with_units_(
        bool& changed,
        bool& index_cleared,
        const OrderUnitGroup& units,
        flint::FmpzConstRef ell,
        flint::FmpzConstRef aux_prime_bound,
        slong max_appends) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturate_local_with_units");
    changed = false;
    index_cleared = false;
    if (!has_factor_base() || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_) ||
        !flint::fmpz_is_prime(ell) ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 || max_appends <= 0) {
        return false;
    }

    flint::Fmpz index;
    if (!sync_row_module_checkpoint_() ||
        !row_module_.full_rank_index(flint::FmpzRef(index))) {
        return false;
    }

    bool include_torsion = false;
    if (!relation_saturation_torsion_needed(include_torsion, units, ell)) {
        return false;
    }

    bool cleared = !flint::fmpz_divisible(flint::FmpzConstRef(index), ell);
    bool any_changed = false;
    slong appends = 0;
    while (!cleared && appends < max_appends) {
        bool appended = false;
        flint::Fmpz q;
        flint::fmpz_set_ui(flint::FmpzRef(q), 2);
        while (flint::fmpz_cmp(flint::FmpzConstRef(q), aux_prime_bound) <= 0) {
            PrimeIdealList local;
            if (decompose_prime(local, parent_,
                                flint::FmpzConstRef(q))) {
                for (slong i = 0; i < local.size(); ++i) {
                    const PrimeIdeal* prime = local.at(i);
                    if (prime == nullptr ||
                        !auxiliary_prime_usable(*prime,
                                                flint::FmpzConstRef(q), ell)) {
                        continue;
                    }

                    bool local_changed = false;
                    if (saturate_local_once_with_units_(
                                local_changed, units, include_torsion, *prime,
                                ell) &&
                        local_changed) {
                        appended = true;
                        break;
                    }
                }
            }
            if (appended) {
                break;
            }
            flint::fmpz_nextprime(flint::FmpzRef(q), flint::FmpzConstRef(q));
        }

        if (!appended) {
            break;
        }

        any_changed = true;
        ++appends;
        if (!sync_row_module_checkpoint_() ||
            !row_module_.full_rank_index(flint::FmpzRef(index))) {
            return false;
        }
        cleared = !flint::fmpz_divisible(flint::FmpzConstRef(index), ell);
    }

    changed = any_changed;
    index_cleared = cleared;
    return true;
}

bool ClassGroupContext::saturate_relations_bounded_with_units(
        bool& changed,
        bool& saturated,
        const OrderUnitGroup& units,
        flint::FmpzConstRef aux_prime_bound,
        slong max_appends_per_ell,
        slong max_appends_total) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "class_group.saturate_relations_bounded_with_units");
    changed = false;
    saturated = false;
    if (!has_factor_base() || !has_presentation() || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_) ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 ||
        max_appends_per_ell <= 0 || max_appends_total <= 0) {
        return false;
    }

    flint::Fmpz index;
    if (!sync_row_module_checkpoint_() ||
        !row_module_.full_rank_index(flint::FmpzRef(index))) {
        return false;
    }

    return saturate_relations_bounded_for_index_with_units_(
            changed, saturated, units, flint::FmpzConstRef(index),
            aux_prime_bound, max_appends_per_ell, max_appends_total);
}

bool ClassGroupContext::saturate_relations_bounded_for_index_with_units_(
        bool& changed,
        bool& saturated,
        const OrderUnitGroup& units,
        flint::FmpzConstRef index_bound,
        flint::FmpzConstRef aux_prime_bound,
        slong max_appends_per_ell,
        slong max_appends_total) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::class_group,
            "class_group.saturate_relations_bounded_for_index_with_units");
    changed = false;
    saturated = false;
    if (!has_factor_base() || !has_presentation() || !units.is_set() ||
        !same_order_parent(units.parent(), &parent_) ||
        flint::fmpz_cmp_ui(index_bound, 1) < 0 ||
        flint::fmpz_cmp_ui(aux_prime_bound, 2) < 0 ||
        max_appends_per_ell <= 0 || max_appends_total <= 0) {
        return false;
    }

    flint::Fmpz current_index;
    if (!sync_row_module_checkpoint_() ||
        !row_module_.full_rank_index(flint::FmpzRef(current_index)) ||
        !flint::fmpz_divisible(flint::FmpzConstRef(current_index),
                               index_bound)) {
        return false;
    }

    flint::FmpzFactor factorization;
    flint::fmpz_factor(flint::FmpzFactorRef(factorization),
                       index_bound);
    if (flint::fmpz_factor_num(flint::FmpzFactorConstRef(factorization)) ==
        0) {
        saturated = true;
        return true;
    }

    private_storage_->relation_saturation_records.clear();
    bool any_changed = false;
    bool all_cleared = true;
    slong remaining = max_appends_total;
    const bool native_rank_zero_quadratic =
            units.free_rank() == 0 && class_unit_transaction_context_ != nullptr &&
            !detail::uses_class_unit_kernel(
                    class_unit_transaction_context_) &&
            parent_.is_maximal() &&
            detail::order_supports_exact_quadratic_class_certificate(parent_);
    for (slong i = 0;
         i < flint::fmpz_factor_num(flint::FmpzFactorConstRef(factorization));
         ++i) {
        flint::Fmpz ell;
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(ell),
                flint::FmpzFactorConstRef(factorization), i);

        bool local_changed = false;
        bool index_cleared = false;
        if (remaining > 0) {
            const slong attempt_limit =
                    min_slong_value(max_appends_per_ell, remaining);
            const slong before = relation_count();
            if (native_rank_zero_quadratic) {
                // Silex C keeps class-only and unit-aware residue-root
                // saturation as separate exact paths.  Rank-zero quadratic
                // native computation first tries the smaller class-only
                // kernel; if it finds no refining relation, the existing
                // torsion-aware path remains authoritative.
                if (!saturate_local_(
                            local_changed, index_cleared,
                            flint::FmpzConstRef(ell), aux_prime_bound,
                            attempt_limit)) {
                    return false;
                }
                if (!local_changed && !index_cleared &&
                    !saturate_local_with_units_(
                            local_changed, index_cleared, units,
                            flint::FmpzConstRef(ell), aux_prime_bound,
                            attempt_limit)) {
                    return false;
                }
            } else if (!saturate_local_with_units_(
                               local_changed, index_cleared, units,
                               flint::FmpzConstRef(ell), aux_prime_bound,
                               attempt_limit)) {
                return false;
            }
            remaining -= relation_count() - before;
        }

        any_changed = any_changed || local_changed;
        all_cleared = all_cleared && index_cleared;
        if (!mark_relation_saturation_(
                    flint::FmpzConstRef(ell),
                    index_cleared ? ProofState::verified
                                  : ProofState::unavailable)) {
            return false;
        }
    }

    if (any_changed && !publish_presentation()) {
        return false;
    }
    relation_saturation_status_ =
            all_cleared ? ProofState::verified : ProofState::unavailable;
    if (units.certification_status() == CertificationMode::proven) {
        unit_proof_status_ = ProofState::verified;
        regulator_proof_status_ = ProofState::verified;
    }
    changed = any_changed;
    saturated = all_cleared;
    return true;
}

bool ClassGroupContext::presentation(FiniteAbelianGroup& out) const noexcept {
    return has_presentation() && out.set(quotient_);
}

std::optional<FiniteAbelianGroup> ClassGroupContext::presentation()
        const noexcept {
    FiniteAbelianGroup out;
    if (!presentation(out)) {
        return std::nullopt;
    }
    return out;
}

slong ClassGroupContext::invariant_count() const noexcept {
    return has_presentation() ? quotient_.invariant_count() : 0;
}

bool ClassGroupContext::invariant(flint::FmpzRef out,
                                  slong index) const noexcept {
    return has_presentation() && quotient_.invariant(out, index);
}

std::optional<flint::Fmpz> ClassGroupContext::invariant(
        slong index) const noexcept {
    flint::Fmpz out;
    if (!invariant(flint::FmpzRef(out), index)) {
        return std::nullopt;
    }
    return out;
}

bool ClassGroupContext::invariants(flint::FmpzVecRef out) const noexcept {
    return has_presentation() && quotient_.invariants(out);
}

bool ClassGroupContext::order(flint::FmpzRef out) const noexcept {
    return has_presentation() && quotient_.order(out);
}

std::optional<flint::Fmpz> ClassGroupContext::order() const noexcept {
    flint::Fmpz out;
    if (!order(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ClassGroupContext::invariant_generator_matrix(
        flint::FmpzMatRef out) const noexcept {
    return has_presentation() &&
           quotient_.invariant_generator_matrix(out);
}

std::optional<flint::FmpzMat> ClassGroupContext::invariant_generator_matrix()
        const noexcept {
    if (!has_presentation()) {
        return std::nullopt;
    }
    flint::FmpzMat out(invariant_count(), generator_count());
    if (!invariant_generator_matrix(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool ClassGroupContext::invariant_generator(
        FractionalIdeal& out,
        slong index) const noexcept {
    if (!has_presentation() || !same_order_parent(out.parent(), &parent_) ||
        index < 0 || index >= invariant_count()) {
        return false;
    }

    FractionalIdeal accumulator(parent_);
    FractionalIdeal prime_ideal(parent_);
    FractionalIdeal power(parent_);
    PrimeIdeal prime(parent_);
    flint::FmpzMat coordinates(invariant_count(), generator_count());
    if (!accumulator.is_defined() || !prime_ideal.is_defined() ||
        !power.is_defined() || !prime.is_defined() ||
        !invariant_generator_matrix(flint::FmpzMatRef(coordinates)) ||
        !accumulator.one()) {
        return false;
    }

    for (slong i = 0; i < generator_count(); ++i) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(coordinates),
                                      index, i);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }

        if (!base_.prime(prime, i) ||
            !detail::prime_to_fractional_ideal(prime_ideal, prime) ||
            !power.pow_fmpz(prime_ideal, exponent) ||
            !accumulator.multiply(accumulator, power)) {
            return false;
        }
    }

    out.swap(accumulator);
    return true;
}

bool ClassGroupContext::ideal_class_coordinates(
        flint::FmpzMatRef out,
        const FractionalIdeal& ideal) const noexcept {
    if (!has_presentation() || !same_order_parent(ideal.parent(), &parent_) ||
        flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != invariant_count()) {
        return false;
    }

    flint::FmpzMat row(1, generator_count());
    flint::FmpzMat candidate(1, invariant_count());
    if (!ideal_factor_over_base(flint::FmpzMatRef(row), ideal, base_) ||
        !quotient_.invariant_coordinates(
                flint::FmpzMatRef(candidate), flint::FmpzMatConstRef(row))) {
        return false;
    }

    flint::fmpz_mat_set(out, flint::FmpzMatConstRef(candidate));
    return true;
}

bool ClassGroupContext::invariant_generator_power_witness(
        FactoredElement& out,
        slong index) const noexcept {
    if (!has_presentation() || parent_.parent() == nullptr ||
        out.parent() == nullptr ||
        !out.parent()->has_same_data(*parent_.parent()) ||
        index < 0 || index >= invariant_count()) {
        return false;
    }

    flint::FmpzMat coefficients(invariant_count(), relation_count());
    if (!quotient_.invariant_generator_relation_matrix(
                flint::FmpzMatRef(coefficients))) {
        return false;
    }
    return push_relation_witnesses(out, relations_,
                                   flint::FmpzMatConstRef(coefficients),
                                   index);
}

slong ClassGroupContext::relation_kernel_unit_count() const noexcept {
    return has_presentation() ? quotient_.relation_kernel_count() : 0;
}

bool ClassGroupContext::relation_kernel_unit(
        FactoredElement& out,
        slong index) const noexcept {
    if (!has_presentation() || parent_.parent() == nullptr ||
        out.parent() == nullptr ||
        !out.parent()->has_same_data(*parent_.parent()) ||
        index < 0 || index >= relation_kernel_unit_count()) {
        return false;
    }

    flint::FmpzMat coefficients(1, relation_count());
    if (!quotient_.relation_kernel_row(
                flint::FmpzMatRef(coefficients), index)) {
        return false;
    }
    return push_relation_witnesses(out, relations_,
                                   flint::FmpzMatConstRef(coefficients), 0);
}

}  // namespace silex
