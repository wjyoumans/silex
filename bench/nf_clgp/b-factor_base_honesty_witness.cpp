#include <benchmark/benchmark.h>

#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/lat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>

#include "benchmark_contract.hpp"
#include "class_group/factor_base_proof_targets_internal.hpp"
#include "class_group/relation_candidate_internal.hpp"
#include "class_group/relation_factor_base_plan_internal.hpp"
#include "ideal_factorization/ideal_factorization_internal.hpp"

#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

constexpr slong kHonestyScalarRadius = 8;

enum class HonestyWorkload {
    quartic,
    quintic,
};

constexpr std::size_t kQuarticScalarCandidates = 570;
constexpr std::size_t kQuarticEnumeratedCandidates = 102;
constexpr std::size_t kQuarticCandidateCount =
        kQuarticScalarCandidates + kQuarticEnumeratedCandidates;
constexpr std::size_t kQuinticScalarCandidates = 3240;
constexpr std::size_t kQuinticEnumeratedCandidates = 697;
constexpr std::size_t kQuinticCandidateCount =
        kQuinticScalarCandidates + kQuinticEnumeratedCandidates;

// Frozen output of exhaustive factorization over each deterministic fixture.
// Bit i (least-significant bit first) is the expected classification of
// candidate i.  Comparing every classification to these bitmaps prevents a
// count-preserving permutation of true and false results from passing.
constexpr std::array<std::uint64_t, 11> kQuarticExpectedPattern = {
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xfc00000000000000),
        UINT64_C(0xc081fe1840c18101),
        UINT64_C(0x00000000e010381f),
};

constexpr std::array<std::uint64_t, 62> kQuinticExpectedPattern = {
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x3007e10000000000),
        UINT64_C(0x38180200f00601c0),
        UINT64_C(0x000f81807e040204),
        UINT64_C(0x201c08101c000000),
        UINT64_C(0x87c07f00c3e0381c),
        UINT64_C(0x01ff00803e06020f),
        UINT64_C(0xe000000000181f82),
        UINT64_C(0x07fe01f003804080),
        UINT64_C(0x00f000fffc087808),
        UINT64_C(0x070000000000fffc),
        UINT64_C(0xffe01fc00fc0f008),
        UINT64_C(0x000000010381fff0),
};

static_assert(kQuarticExpectedPattern.size() ==
              (kQuarticCandidateCount + 63U) / 64U);
static_assert(kQuinticExpectedPattern.size() ==
              (kQuinticCandidateCount + 63U) / 64U);

struct WitnessPatternContract {
    const std::uint64_t* words = nullptr;
    std::size_t word_count = 0;
    std::size_t scalar_candidate_count = 0;
    std::size_t enumerated_candidate_count = 0;
    std::size_t searched_prime_count = 0;
    slong scalar_match_count = 0;
    slong total_match_count = 0;
};

WitnessPatternContract witness_pattern_contract(
        HonestyWorkload workload) noexcept {
    if (workload == HonestyWorkload::quartic) {
        return {kQuarticExpectedPattern.data(),
                kQuarticExpectedPattern.size(),
                kQuarticScalarCandidates,
                kQuarticEnumeratedCandidates,
                38,
                0,
                38};
    }
    return {kQuinticExpectedPattern.data(),
            kQuinticExpectedPattern.size(),
            kQuinticScalarCandidates,
            kQuinticEnumeratedCandidates,
            216,
            0,
            216};
}

bool expected_candidate_match(const WitnessPatternContract& contract,
                              std::size_t index) noexcept {
    return index / 64U < contract.word_count &&
           ((contract.words[index / 64U] >> (index % 64U)) & UINT64_C(1)) !=
                   UINT64_C(0);
}

struct FieldFixture {
    silex::NumberField field;
    silex::Order maximal_order;
};

enum class WitnessCandidateKind {
    scalar,
    enumerated,
};

struct WitnessCandidate {
    WitnessCandidate(const silex::Order& order,
                     WitnessCandidateKind candidate_kind) noexcept
        : scalar_element(*order.parent()),
          order_element(order),
          ideal(order),
          required_prime(order),
          kind(candidate_kind) {}

    silex::Element scalar_element;
    silex::OrderElement order_element;
    silex::Ideal ideal;
    silex::PrimeIdeal required_prime;
    WitnessCandidateKind kind = WitnessCandidateKind::scalar;
};

struct WitnessBatch {
    FieldFixture field;
    silex::FactorBase base;
    std::vector<WitnessCandidate> candidates;
    std::vector<silex::PrimeIdeal> searched_primes;
    std::size_t scalar_candidate_count = 0;
    slong scalar_match_count = 0;
    slong total_match_count = 0;
};

void set_workload_polynomial(sflint::FmpqPoly& polynomial,
                             HonestyWorkload workload) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    if (workload == HonestyWorkload::quartic) {
        sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 3, -8);
        sflint::fmpq_poly_set_coeff_si(polynomial, 2, 3);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, 4);
        return;
    }

    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, -7);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -6);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 3);
}

bool define_workload_field(FieldFixture& fixture,
                           HonestyWorkload workload) noexcept {
    sflint::FmpqPoly polynomial;
    set_workload_polynomial(polynomial, workload);

    silex::Order equation_order;
    return fixture.field.define_by_polynomial(
                   sflint::FmpqPolyConstRef(polynomial)) &&
           equation_order.define_equation_order(fixture.field) &&
           fixture.maximal_order.define(fixture.field) &&
           fixture.maximal_order.maximal_order(equation_order) &&
           fixture.maximal_order.is_maximal();
}

bool materializing_factor_over_base_with_required_prime(
        bool& matches,
        const silex::FactorBase& base,
        const WitnessCandidate& candidate) noexcept {
    matches = false;
    const silex::Order* order = base.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !silex::same_order_parent(candidate.order_element.parent(), order) ||
        !silex::same_order_parent(candidate.required_prime.parent(), order)) {
        return false;
    }

    silex::Element callback_element(*order->parent());
    const silex::Element* source = &candidate.scalar_element;
    if (candidate.kind == WitnessCandidateKind::enumerated) {
        if (!callback_element.is_defined() ||
            !candidate.order_element.get_element(callback_element)) {
            return false;
        }
        source = &callback_element;
    }

    silex::OrderElement generator(*order);
    silex::Ideal principal(*order);
    return generator.is_defined() && principal.is_defined() &&
           generator.set_element(*source) && principal.set_principal(generator) &&
           silex::detail::ideal_factor_over_base_with_required_prime(
                   matches, principal, base, candidate.required_prime);
}

bool candidate_order_element_factor_over_base_with_required_prime(
        bool& matches,
        const silex::FactorBase& base,
    const WitnessCandidate& candidate) noexcept {
    if (candidate.kind == WitnessCandidateKind::enumerated) {
        return silex::detail::
                order_element_factor_over_base_with_required_prime(
                matches, candidate.order_element, base,
                candidate.required_prime);
    }

    const silex::Order* order = base.parent();
    if (order == nullptr) {
        return false;
    }
    silex::OrderElement order_element(*order);
    return order_element.is_defined() &&
           order_element.set_element(candidate.scalar_element) &&
           silex::detail::order_element_factor_over_base_with_required_prime(
                   matches, order_element, base, candidate.required_prime);
}

bool append_scalar_candidate(std::vector<WitnessCandidate>& candidates,
                             const silex::Order& order,
                             const silex::PrimeIdeal& required_prime,
                             slong scalar) noexcept {
    WitnessCandidate candidate(order, WitnessCandidateKind::scalar);
    if (!candidate.ideal.is_defined() ||
        !candidate.required_prime.is_defined() ||
        !candidate.scalar_element.is_defined() ||
        !candidate.order_element.is_defined() ||
        !candidate.scalar_element.set_si(scalar) ||
        !candidate.order_element.set_element(candidate.scalar_element) ||
        !candidate.ideal.set_principal(candidate.order_element) ||
        !candidate.required_prime.set(required_prime)) {
        return false;
    }
    candidates.emplace_back(std::move(candidate));
    return true;
}

bool append_prime_scalar_candidates(
        std::vector<WitnessCandidate>& candidates,
        const silex::Order& order,
        const silex::PrimeIdeal& required_prime,
        slong rational_prime) noexcept {
    if (!append_scalar_candidate(candidates, order, required_prime,
                                 rational_prime)) {
        return false;
    }

    for (slong k = -kHonestyScalarRadius;
         k <= kHonestyScalarRadius; ++k) {
        if (k == -1 || k == 0 || k == 1 || k == rational_prime) {
            continue;
        }
        if (!append_scalar_candidate(candidates, order, required_prime, k)) {
            return false;
        }
    }
    return true;
}

bool append_rational_prime_candidates(
        WitnessBatch& batch,
        sflint::FmpzConstRef rational_prime,
        sflint::FmpzConstRef required_bound) noexcept {
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, batch.field.maximal_order,
                                rational_prime)) {
        return false;
    }

    std::vector<slong> proof_targets;
    if (!silex::detail::select_factor_base_proof_targets(
                proof_targets, primes, rational_prime, required_bound)) {
        return false;
    }

    if (!sflint::fmpz_fits_si(rational_prime)) {
        return false;
    }
    const slong rational_prime_si = sflint::fmpz_get_si(rational_prime);
    for (slong index : proof_targets) {
        const silex::PrimeIdeal* prime = primes.at(index);
        if (prime == nullptr) {
            return false;
        }
        if (!batch.base.contains(*prime)) {
            if (!append_prime_scalar_candidates(
                        batch.candidates, batch.field.maximal_order, *prime,
                        rational_prime_si)) {
                return false;
            }
            silex::PrimeIdeal searched(batch.field.maximal_order);
            if (!searched.is_defined() || !searched.set(*prime)) {
                return false;
            }
            batch.searched_primes.emplace_back(std::move(searched));
        }
    }
    return true;
}

struct EnumeratedCandidateContext {
    WitnessBatch* batch = nullptr;
    const silex::PrimeIdeal* required_prime = nullptr;
    const silex::Order* order = nullptr;
    sflint::FmpzMat* basis = nullptr;
    bool found = false;
    bool failed = false;
};

int append_enumerated_candidate(const fmpz_mat_t coefficients, void* user) {
    auto* context = static_cast<EnumeratedCandidateContext*>(user);
    if (context == nullptr || context->batch == nullptr ||
        context->required_prime == nullptr || context->order == nullptr ||
        context->basis == nullptr || context->found || context->failed) {
        return 0;
    }

    sflint::FmpzMat coordinates(1, context->order->degree());
    silex::detail::relation_search::coordinates_from_lattice_combination(
            coordinates, sflint::FmpzMatConstRef(coefficients),
            sflint::FmpzMatConstRef(*context->basis));

    WitnessCandidate candidate(*context->order,
                               WitnessCandidateKind::enumerated);
    if (!candidate.order_element.is_defined() ||
        !candidate.ideal.is_defined() ||
        !candidate.required_prime.is_defined() ||
        !candidate.order_element.set_coordinates(
                sflint::FmpzMatConstRef(coordinates)) ||
        !candidate.ideal.set_principal(candidate.order_element) ||
        !candidate.required_prime.set(*context->required_prime)) {
        context->failed = true;
        return 0;
    }
    context->batch->candidates.emplace_back(std::move(candidate));

    bool matches = false;
    if (!materializing_factor_over_base_with_required_prime(
                matches, context->batch->base,
                context->batch->candidates.back())) {
        context->failed = true;
        return 0;
    }
    if (matches) {
        context->found = true;
        return 0;
    }
    return 1;
}

bool append_prime_enumerated_candidates(
        WitnessBatch& batch,
        const silex::PrimeIdeal& required_prime) noexcept {
    const silex::Order* order = batch.base.parent();
    if (order == nullptr ||
        !silex::same_order_parent(required_prime.parent(), order) ||
        batch.base.contains(required_prime)) {
        return false;
    }

    silex::Ideal ideal(*order);
    sflint::FmpzMat hnf(order->degree(), order->degree());
    silex::lat::Lat lattice(order->degree());
    silex::lat::Lat reduced(order->degree());
    if (!ideal.is_defined() || !required_prime.get_ideal(ideal) ||
        !ideal.get_hnf(sflint::FmpzMatRef(hnf)) ||
        !lattice.set_basis(sflint::FmpzMatConstRef(hnf)) ||
        !lattice.lll_reduce(reduced) || reduced.nrows() <= 0) {
        return false;
    }

    sflint::FmpzMat basis(reduced.nrows(), order->degree());
    sflint::Fmpz max_diag;
    sflint::Fmpz gram_ii;
    sflint::Arb bound_sq;
    if (!reduced.get_basis(sflint::FmpzMatRef(basis))) {
        return false;
    }

    sflint::fmpz_zero(sflint::FmpzRef(max_diag));
    for (slong i = 0; i < sflint::fmpz_mat_nrows(basis); ++i) {
        sflint::fmpz_zero(sflint::FmpzRef(gram_ii));
        for (slong k = 0; k < sflint::fmpz_mat_ncols(basis); ++k) {
            const sflint::FmpzConstRef entry = sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(basis), i, k);
            sflint::fmpz_addmul(sflint::FmpzRef(gram_ii), entry, entry);
        }
        if (sflint::fmpz_cmp(sflint::FmpzConstRef(gram_ii),
                             sflint::FmpzConstRef(max_diag)) > 0) {
            sflint::fmpz_set(sflint::FmpzRef(max_diag),
                             sflint::FmpzConstRef(gram_ii));
        }
    }
    sflint::fmpz_mul_ui(sflint::FmpzRef(max_diag),
                        sflint::FmpzConstRef(max_diag),
                        static_cast<ulong>(reduced.nrows()));

    EnumeratedCandidateContext context;
    context.batch = &batch;
    context.required_prime = &required_prime;
    context.order = order;
    context.basis = &basis;
    for (slong tries = 0; tries < 4 && !context.found && !context.failed;
         ++tries) {
        sflint::arb_set_fmpz(bound_sq, max_diag);
        const bool enum_ok = reduced.enum_short_vectors_arb(
                sflint::ArbConstRef(bound_sq), kHonestyScalarRadius, 64,
                append_enumerated_candidate, &context);
        if (!enum_ok && !context.found && !context.failed) {
            (void) reduced.enum_short_vectors_arb(
                    sflint::ArbConstRef(bound_sq),
                    kHonestyScalarRadius, 256,
                    append_enumerated_candidate, &context);
        }
        sflint::fmpz_mul_ui(sflint::FmpzRef(max_diag),
                            sflint::FmpzConstRef(max_diag), UWORD(2));
    }

    return context.found && !context.failed;
}

bool build_witness_batch(WitnessBatch& batch,
                         HonestyWorkload workload) noexcept {
    const WitnessPatternContract contract =
            witness_pattern_contract(workload);
    if (!define_workload_field(batch.field, workload) ||
        !batch.base.define(batch.field.maximal_order)) {
        return false;
    }

    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    silex::detail::RelationFactorBasePlan factor_base_plan;
    if (!silex::detail::relation_search::build_relation_factor_base_plan(
                factor_base_plan, batch.field.maximal_order, nullptr)) {
        return false;
    }
    sflint::fmpz_set(
            sflint::FmpzRef(active_bound),
            sflint::FmpzConstRef(factor_base_plan.working_bound));
    if (!batch.base.build_relation_completion_base(
                sflint::FmpzConstRef(active_bound)) ||
        !silex::factor_base_class_group_bound(
                sflint::FmpzRef(required_bound),
                batch.field.maximal_order)) {
        return false;
    }

    batch.candidates.reserve(contract.scalar_candidate_count +
                             contract.enumerated_candidate_count);

    sflint::Fmpz rational_prime;
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(active_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(rational_prime), 2);
    } else {
        sflint::fmpz_nextprime(sflint::FmpzRef(rational_prime),
                               sflint::FmpzConstRef(active_bound), true);
    }
    while (sflint::fmpz_cmp(sflint::FmpzConstRef(rational_prime),
                            sflint::FmpzConstRef(required_bound)) <= 0) {
        if (!append_rational_prime_candidates(
                    batch, sflint::FmpzConstRef(rational_prime),
                    sflint::FmpzConstRef(required_bound))) {
            return false;
        }
        sflint::fmpz_nextprime(sflint::FmpzRef(rational_prime),
                               sflint::FmpzConstRef(rational_prime), true);
    }
    batch.scalar_candidate_count = batch.candidates.size();
    if (batch.scalar_candidate_count != contract.scalar_candidate_count) {
        return false;
    }

    for (const silex::PrimeIdeal& required_prime : batch.searched_primes) {
        if (!append_prime_enumerated_candidates(batch, required_prime)) {
            return false;
        }
    }
    return batch.candidates.size() ==
           contract.scalar_candidate_count +
                   contract.enumerated_candidate_count;
}

bool full_factorization_has_principal_witness(
        bool& matches,
        const silex::FactorBase& base,
        const WitnessCandidate& candidate) noexcept {
    matches = false;
    const silex::Order* order = base.parent();
    if (order == nullptr ||
        !silex::same_order_parent(candidate.ideal.parent(), order) ||
        !silex::same_order_parent(candidate.required_prime.parent(), order) ||
        base.contains(candidate.required_prime)) {
        return false;
    }

    silex::IdealFactorization factorization(*order);
    silex::PrimeIdeal factor(*order);
    if (!factorization.factor(candidate.ideal)) {
        return false;
    }

    slong required_exponent = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = 0;
        if (!factorization.prime(factor, i) ||
            !factorization.exponent(exponent, i)) {
            return false;
        }
        if (factor.equal(candidate.required_prime)) {
            required_exponent += exponent;
        } else if (!base.contains(factor)) {
            return true;
        }
    }
    matches = required_exponent == 1;
    return true;
}

bool validate_batch(WitnessBatch& batch,
                    HonestyWorkload workload) noexcept {
    const WitnessPatternContract contract =
            witness_pattern_contract(workload);
    const std::size_t expected_candidate_count =
            contract.scalar_candidate_count +
            contract.enumerated_candidate_count;
    if (batch.candidates.size() != expected_candidate_count ||
        batch.scalar_candidate_count != contract.scalar_candidate_count ||
        batch.candidates.size() - batch.scalar_candidate_count !=
                contract.enumerated_candidate_count ||
        batch.searched_primes.size() !=
                contract.searched_prime_count) {
        return false;
    }

    slong scalar_match_count = 0;
    slong match_count = 0;
    for (std::size_t i = 0; i < batch.candidates.size(); ++i) {
        const WitnessCandidate& candidate = batch.candidates[i];
        const WitnessCandidateKind expected_kind =
                i < batch.scalar_candidate_count
                        ? WitnessCandidateKind::scalar
                        : WitnessCandidateKind::enumerated;
        if (candidate.kind != expected_kind) {
            return false;
        }

        bool ideal_matches = false;
        bool materializing_matches = false;
        bool element_matches = false;
        bool full_matches = false;
        if (!silex::detail::ideal_factor_over_base_with_required_prime(
                    ideal_matches, candidate.ideal, batch.base,
                    candidate.required_prime) ||
            !materializing_factor_over_base_with_required_prime(
                    materializing_matches, batch.base, candidate) ||
            !candidate_order_element_factor_over_base_with_required_prime(
                    element_matches, batch.base, candidate) ||
            !full_factorization_has_principal_witness(
                    full_matches, batch.base, candidate)) {
            return false;
        }
        const bool expected_matches =
                expected_candidate_match(contract, i);
        if (ideal_matches != expected_matches ||
            materializing_matches != expected_matches ||
            element_matches != expected_matches ||
            full_matches != expected_matches) {
            return false;
        }
        match_count += ideal_matches ? 1 : 0;
        if (i < batch.scalar_candidate_count) {
            scalar_match_count += ideal_matches ? 1 : 0;
        }
    }
    if (match_count != contract.total_match_count ||
        scalar_match_count != contract.scalar_match_count) {
        return false;
    }
    batch.scalar_match_count = scalar_match_count;
    batch.total_match_count = match_count;
    return true;
}

enum class WitnessPredicate {
    full_factorization,
    prebuilt_required_prime,
    current_materializing,
    order_element_direct,
};

void benchmark_witness_batch(benchmark::State& state,
                             HonestyWorkload workload,
                             WitnessPredicate predicate,
                             bool scalar_only) {
    silex::bench_contract::initialize(state);
    WitnessBatch batch;
    bool setup_complete = false;
    slong matches = 0;
    for (auto _ : state) {
        if (!setup_complete) {
            state.PauseTiming();
            const bool setup_ok = build_witness_batch(batch, workload) &&
                                  validate_batch(batch, workload);
            state.ResumeTiming();
            if (!setup_ok) {
                silex::bench_contract::fail(
                        state, "honesty witness batch setup failed",
                        silex::bench_contract::FailureReason::setup);
                return;
            }
            setup_complete = true;
        }

        matches = 0;
        const std::size_t candidate_count =
                scalar_only ? batch.scalar_candidate_count
                            : batch.candidates.size();
        for (std::size_t i = 0; i < candidate_count; ++i) {
            const WitnessCandidate& candidate = batch.candidates[i];
            bool candidate_matches = false;
            switch (predicate) {
                case WitnessPredicate::full_factorization:
                    if (!full_factorization_has_principal_witness(
                                candidate_matches, batch.base, candidate)) {
                        silex::bench_contract::fail(
                                state,
                                "full factorization witness predicate failed",
                                silex::bench_contract::FailureReason::operation);
                        return;
                    }
                    break;
                case WitnessPredicate::prebuilt_required_prime:
                    if (!silex::detail::
                                ideal_factor_over_base_with_required_prime(
                                        candidate_matches, candidate.ideal,
                                        batch.base,
                                        candidate.required_prime)) {
                        silex::bench_contract::fail(
                                state,
                                "required-prime factor-over-base predicate failed",
                                silex::bench_contract::FailureReason::operation);
                        return;
                    }
                    break;
                case WitnessPredicate::current_materializing:
                    if (!materializing_factor_over_base_with_required_prime(
                                candidate_matches, batch.base, candidate)) {
                        silex::bench_contract::fail(
                                state,
                                "materializing witness predicate failed",
                                silex::bench_contract::FailureReason::operation);
                        return;
                    }
                    break;
                case WitnessPredicate::order_element_direct:
                    if (!candidate_order_element_factor_over_base_with_required_prime(
                                candidate_matches, batch.base, candidate)) {
                        silex::bench_contract::fail(
                                state,
                                "order-element witness predicate failed",
                                silex::bench_contract::FailureReason::operation);
                        return;
                    }
                    break;
            }
            matches += candidate_matches ? 1 : 0;
        }
        benchmark::DoNotOptimize(matches);
    }

    const slong expected_matches = scalar_only ? batch.scalar_match_count
                                               : batch.total_match_count;
    if (!setup_complete || matches != expected_matches) {
        silex::bench_contract::fail(
                state,
                "timed honesty witness pattern disagrees with fixture contract",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }

    state.SetItemsProcessed(
            state.iterations() *
            static_cast<std::int64_t>(
                    scalar_only ? batch.scalar_candidate_count
                                : batch.candidates.size()));
    state.counters["candidates"] =
            static_cast<double>(scalar_only ? batch.scalar_candidate_count
                                            : batch.candidates.size());
    state.counters["enumerated_candidates"] =
            static_cast<double>(scalar_only
                                        ? 0
                                        : batch.candidates.size() -
                                                  batch.scalar_candidate_count);
    state.counters["matches"] = static_cast<double>(matches);
    state.counters["reference_matches"] =
            static_cast<double>(expected_matches);
    state.counters["classification_mismatches"] = 0.0;
    state.counters["classification_pattern_candidates"] =
            static_cast<double>(scalar_only ? batch.scalar_candidate_count
                                            : batch.candidates.size());
    benchmark::ClobberMemory();
    silex::bench_contract::succeed(state);
}

void BM_honesty_scalar_full_factor_quartic(benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quartic,
                            WitnessPredicate::full_factorization, true);
}

void BM_honesty_scalar_required_prime_quartic(
        benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quartic,
                            WitnessPredicate::prebuilt_required_prime, true);
}

void BM_honesty_scalar_full_factor_quintic(benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quintic,
                            WitnessPredicate::full_factorization, true);
}

void BM_honesty_scalar_required_prime_quintic(
        benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quintic,
                            WitnessPredicate::prebuilt_required_prime, true);
}

void BM_honesty_current_materializing_quartic(benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quartic,
                            WitnessPredicate::current_materializing, false);
}

void BM_honesty_order_element_direct_quartic(benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quartic,
                            WitnessPredicate::order_element_direct, false);
}

void BM_honesty_current_materializing_quintic(benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quintic,
                            WitnessPredicate::current_materializing, false);
}

void BM_honesty_order_element_direct_quintic(benchmark::State& state) {
    benchmark_witness_batch(state, HonestyWorkload::quintic,
                            WitnessPredicate::order_element_direct, false);
}

}  // namespace

BENCHMARK(BM_honesty_scalar_full_factor_quartic);
BENCHMARK(BM_honesty_scalar_required_prime_quartic);
BENCHMARK(BM_honesty_scalar_full_factor_quintic);
BENCHMARK(BM_honesty_scalar_required_prime_quintic);
BENCHMARK(BM_honesty_current_materializing_quartic);
BENCHMARK(BM_honesty_order_element_direct_quartic);
BENCHMARK(BM_honesty_current_materializing_quintic);
BENCHMARK(BM_honesty_order_element_direct_quintic);

BENCHMARK_MAIN();
