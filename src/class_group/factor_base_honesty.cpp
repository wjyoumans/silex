#include "factor_base_honesty_internal.hpp"

#include <silex/ideal_factorization.hpp>
#include <silex/lat.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>

#include "class_group_internal.hpp"
#include "factor_base_proof_targets_internal.hpp"
#include "ideal_lattice_reduction_internal.hpp"
#include "ideal_t2_enumeration_internal.hpp"
#include "relation_candidate_internal.hpp"
#include "relation_search_internal.hpp"
#include "../ideal_factorization/ideal_factorization_internal.hpp"

#include <cstdio>
#include <vector>

namespace silex::detail::relation_search {

// Stage-0 witness search effort.  The 50 random tries follow the reference
// buch2.c:be_honest (maxtry_HONEST); the T2 caps are the shared per-ideal
// enumeration caps of ideal_t2_enumeration_internal.hpp.
constexpr slong kFactorBaseHonestySearchRadius = 8;
constexpr slong kFactorBaseHonestyMaxTwists = 16;
constexpr slong kFactorBaseHonestyMaxRandomTries = 50;

namespace {

struct FactorBaseWitnessSearchEffort {
    slong stage = 0;
    slong radius = kFactorBaseHonestySearchRadius;
    slong max_twists = kFactorBaseHonestyMaxTwists;
    slong random_tries = kFactorBaseHonestyMaxRandomTries;
    slong factor_attempts = kMaxFactorAttempts;
    slong element_steps = kMaxElementSteps;
};

// Each escalation stage doubles every cap.  The element-step limit doubles
// with the factor-attempt cap rather than following kMaxElementSteps =
// 4 * kMaxFactorAttempts^2, which would quadruple it per stage; this keeps
// the growth of every cap linear in the stage multiplier and is deliberate.
FactorBaseWitnessSearchEffort factor_base_witness_search_effort(
        slong stage) noexcept {
    FactorBaseWitnessSearchEffort effort;
    effort.stage = stage;
    for (slong i = 0; i < stage; ++i) {
        effort.radius *= 2;
        effort.max_twists *= 2;
        effort.random_tries *= 2;
        effort.factor_attempts *= 2;
        effort.element_steps *= 2;
    }
    return effort;
}

// Which per-ideal T2 cap stopped an exhausted direct witness search.
struct RequiredPrimeWitnessExhaustion {
    slong ideals = 0;
    slong factor_attempt_cap_hits = 0;
    slong element_step_cap_hits = 0;
};

bool full_factorization_has_principal_witness(
        const FactorBase& base,
        const PrimeIdeal& prime,
        const Ideal& principal) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order) ||
        !same_order_parent(principal.parent(), order) ||
        base.contains(prime)) {
        return false;
    }

    IdealFactorization factorization(*order);
    PrimeIdeal factor(*order);
    if (!factorization.is_defined() || !factor.is_defined() ||
        !factorization.factor(principal)) {
        return false;
    }

    slong prime_exponent = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = 0;
        if (!factorization.prime(factor, i) ||
            !factorization.exponent(exponent, i)) {
            return false;
        }
        if (factor.equal(prime)) {
            prime_exponent += exponent;
            continue;
        }
        if (!base.contains(factor)) {
            return false;
        }
    }

    return prime_exponent == 1;
}

enum class RequiredPrimeWitnessSearchResult {
    found,
    exhausted,
    failed,
};

bool multiply_back_multiplier(
        Element& back_multiplier,
        const Element& factor) noexcept {
    const NumberField* field = back_multiplier.parent();
    if (field == nullptr || !factor.has_parent(*field)) {
        return false;
    }

    Element product(*field);
    if (!product.is_defined() ||
        !product.multiply(back_multiplier, factor)) {
        return false;
    }
    back_multiplier.swap(product);
    return true;
}

}  // namespace

bool factor_base_principal_witness(
        const FactorBase& base,
        const PrimeIdeal& prime,
        const OrderElement& generator,
        FactorBaseWitnessPredicate predicate,
        const DiagnosticsContext* diagnostics,
        FactorBaseWitnessAudit* audit) noexcept {
    if (audit != nullptr) *audit = FactorBaseWitnessAudit{};
    const Order* order = base.parent();
    if (order == nullptr || !same_order_parent(generator.parent(), order) ||
        !same_order_parent(prime.parent(), order) || !prime.has_prime_data() ||
        base.contains(prime) || generator.equal_si(0)) return false;
    if (predicate == FactorBaseWitnessPredicate::selected) {
        predicate = FactorBaseWitnessPredicate::order_element_direct;
    }
    if (predicate == FactorBaseWitnessPredicate::order_element_direct) {
        // Source: buch2.c:divide_p_elt/can_factor account for the complete
        // norm using residue-degree-weighted prime valuations. Only evaluation
        // failure falls back; a completed negative classification is final.
        bool matches = false;
        if (detail::order_element_factor_over_base_with_required_prime(
                    matches, generator, base, prime, diagnostics)) return matches;
    } else if (predicate != FactorBaseWitnessPredicate::full_factorization) {
        return false;
    }
    if (audit != nullptr) audit->used_reference = true;
    Ideal principal(*order);
    return principal.is_defined() && principal.set_principal(generator) &&
           full_factorization_has_principal_witness(base, prime, principal);
}

bool factor_base_scalar_witness(
        const FactorBase& base,
        const PrimeIdeal& prime,
        flint::FmpzConstRef scalar,
        FactorBaseWitnessPredicate predicate,
        const DiagnosticsContext* diagnostics,
        FactorBaseWitnessAudit* audit) noexcept {
    if (audit != nullptr) *audit = FactorBaseWitnessAudit{};
    const Order* order = base.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !prime.has_prime_data() || base.contains(prime) ||
        flint::fmpz_is_zero(scalar)) return false;

    if (predicate == FactorBaseWitnessPredicate::selected) {
        bool matches = false;
        if (detail::scalar_factor_over_base_with_required_prime(
                    matches, scalar, base, prime, diagnostics)) {
            if (audit != nullptr) audit->used_scalar_direct = true;
            return matches;
        }
        // A completed negative is final; only failed evaluations fall back.
        predicate = FactorBaseWitnessPredicate::order_element_direct;
    }
    // Keep the general representation boundary for reference and fallback.
    Element alpha(*order->parent());
    OrderElement generator(*order);
    return alpha.is_defined() && generator.is_defined() &&
           alpha.set_fmpz(scalar) && generator.set_element(alpha) &&
           factor_base_principal_witness(
                   base, prime, generator, predicate, diagnostics, audit);
}

bool factor_base_honesty_primitive_part(
        Ideal& ideal,
        Element& back_multiplier) noexcept {
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong degree = order == nullptr ? 0 : order->degree();
    if (field == nullptr || degree <= 0 || !ideal.has_hnf() ||
        !back_multiplier.has_parent(*field)) {
        return false;
    }

    flint::FmpzMat hnf(degree, degree);
    if (!ideal.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }
    // The source column-HNF final diagonal maps to (0,0) under Silex's
    // row-HNF convention.
    if (flint::fmpz_is_one(flint::fmpz_mat_entry(
                flint::FmpzMatConstRef(hnf), 0, 0))) {
        return true;
    }

    flint::Fmpz content;
    flint::fmpz_zero(flint::FmpzRef(content));
    for (slong row = 0; row < degree; ++row) {
        for (slong column = 0; column < degree; ++column) {
            flint::fmpz_gcd(
                    flint::FmpzRef(content),
                    flint::FmpzConstRef(content),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(hnf), row, column));
        }
    }
    if (flint::fmpz_is_zero(flint::FmpzConstRef(content))) {
        return false;
    }
    if (flint::fmpz_is_one(flint::FmpzConstRef(content))) {
        return true;
    }

    ::fmpz_mat_scalar_divexact_fmpz(hnf.raw(), hnf.raw(), content.raw());
    Ideal primitive(*order);
    Element scalar(*field);
    if (!primitive.is_defined() || !scalar.is_defined() ||
        !primitive.set_hnf(flint::FmpzMatConstRef(hnf)) ||
        !scalar.set_fmpz(flint::FmpzConstRef(content)) ||
        !multiply_back_multiplier(back_multiplier, scalar)) {
        return false;
    }
    ideal.swap(primitive);
    return true;
}

bool factor_base_honesty_reduce_large_ideal(
        Ideal& ideal,
        Element& back_multiplier,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = ideal.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong degree = order == nullptr ? 0 : order->degree();
    if (field == nullptr || degree <= 0 || precision <= 0 ||
        !ideal.has_hnf() || !back_multiplier.has_parent(*field)) {
        return false;
    }

    flint::FmpzMat hnf(degree, degree);
    if (!ideal.get_hnf(flint::FmpzMatRef(hnf))) {
        return false;
    }
    // The source column-HNF first diagonal maps to the final diagonal
    // entry under Silex's row-HNF convention.
    const flint::FmpzConstRef source_first_diagonal =
            flint::fmpz_mat_entry(flint::FmpzMatConstRef(hnf), degree - 1,
                                  degree - 1);
    if (::fmpz_bits(source_first_diagonal.raw()) <= 101) {
        return true;
    }

    Ideal reduced(*order);
    Element reduction_multiplier(*field);
    Element inverse_multiplier(*field);
    if (!reduced.is_defined() || !reduction_multiplier.is_defined() ||
        !inverse_multiplier.is_defined() ||
        !silex::detail::reduce_ideal_lattice(
                reduced, reduction_multiplier, ideal, precision,
                diagnostics) ||
        !inverse_multiplier.invert(reduction_multiplier) ||
        !multiply_back_multiplier(
                back_multiplier, inverse_multiplier)) {
        return false;
    }
    ideal.swap(reduced);
    return true;
}

namespace {

RequiredPrimeWitnessSearchResult enumerate_required_prime_witness(
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const Ideal& ideal,
        const Element& back_multiplier,
        const FactorBaseWitnessSearchEffort& effort,
        RequiredPrimeWitnessExhaustion& exhaustion,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !same_order_parent(required_prime.parent(), order) ||
        !same_order_parent(ideal.parent(), order) || !ideal.has_hnf() ||
        !back_multiplier.has_parent(*field) || base.contains(required_prime)) {
        return RequiredPrimeWitnessSearchResult::failed;
    }

    FiniteIdealT2EnumerationData t2_context;
    if (!build_finite_ideal_t2_enumeration_data_with_retry(
                t2_context, ideal, diagnostics, embedding_cache,
                effort.factor_attempts)) {
        return RequiredPrimeWitnessSearchResult::failed;
    }

    FiniteQuadraticFormEnumerationContext enumeration;
    bool skip_first_scalar = false;
    if (!enumeration.reset(t2_context.quadratic_form_data, order->degree()) ||
        !reduced_basis_first_row_is_scalar_rational(
                skip_first_scalar, *order,
                flint::FmpzMatConstRef(t2_context.basis)) ||
        !enumeration.start(t2_context.initial_bound_value,
                           effort.element_steps, skip_first_scalar)) {
        return RequiredPrimeWitnessSearchResult::failed;
    }

    flint::FmpzMat coefficients(1, order->degree());
    flint::FmpzMat coordinates(1, order->degree());
    OrderElement order_element(*order);
    Element current_element(*field);
    Element original_element(*field);
    if (!order_element.is_defined() || !current_element.is_defined() ||
        !original_element.is_defined()) {
        return RequiredPrimeWitnessSearchResult::failed;
    }

    ++exhaustion.ideals;
    slong factor_attempts = 0;
    while (enumeration.next()) {
        if (!enumeration.current_row(flint::FmpzMatRef(coefficients))) {
            return RequiredPrimeWitnessSearchResult::failed;
        }
        if (!fmpz_mat_single_row_is_primitive(
                    flint::FmpzMatConstRef(coefficients))) {
            continue;
        }

        coordinates_from_lattice_combination(
                coordinates, flint::FmpzMatConstRef(coefficients),
                flint::FmpzMatConstRef(t2_context.basis));
        bool scalar = false;
        if (!order_element.set_coordinates(
                    flint::FmpzMatConstRef(coordinates)) ||
            !order_element.get_element(current_element) ||
            !element_is_scalar_rational(scalar, current_element)) {
            return RequiredPrimeWitnessSearchResult::failed;
        }
        if (scalar) {
            continue;
        }
        if (factor_attempts >= effort.factor_attempts) {
            ++exhaustion.factor_attempt_cap_hits;
            return RequiredPrimeWitnessSearchResult::exhausted;
        }
        ++factor_attempts;

        if (!original_element.multiply(back_multiplier, current_element) ||
            !order_element.set_element(original_element)) {
            return RequiredPrimeWitnessSearchResult::failed;
        }
        bool matches = false;
        if (!detail::order_element_factor_over_base_with_required_prime(
                    matches, order_element, base, required_prime,
                    diagnostics)) {
            return RequiredPrimeWitnessSearchResult::failed;
        }
        if (matches) {
            return RequiredPrimeWitnessSearchResult::found;
        }
    }

    if (enumeration.element_step_limit_reached()) {
        ++exhaustion.element_step_cap_hits;
    }
    return RequiredPrimeWitnessSearchResult::exhausted;
}

RequiredPrimeWitnessSearchResult find_required_prime_witness(
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const SubfactorBaseSchedule* subfactor_base_schedule,
        ulong& random_state,
        slong ideal_reduction_precision,
        const FactorBaseWitnessSearchEffort& effort,
        RequiredPrimeWitnessExhaustion& exhaustion,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || subfactor_base_schedule == nullptr ||
        ideal_reduction_precision <= 0 ||
        !same_order_parent(required_prime.parent(), order) ||
        base.contains(required_prime)) {
        return RequiredPrimeWitnessSearchResult::failed;
    }

    Ideal required_ideal(*order);
    Ideal current(*order);
    Ideal twist_ideal(*order);
    Ideal product(*order);
    PrimeIdeal twist_prime(*order);
    Element back_multiplier(*field);
    if (!required_ideal.is_defined() || !current.is_defined() ||
        !twist_ideal.is_defined() || !product.is_defined() ||
        !twist_prime.is_defined() || !back_multiplier.is_defined() ||
        !required_prime.get_ideal(required_ideal) ||
        !current.set(required_ideal) || !back_multiplier.one()) {
        return RequiredPrimeWitnessSearchResult::failed;
    }

    RequiredPrimeWitnessSearchResult result = enumerate_required_prime_witness(
            base, required_prime, current, back_multiplier, effort,
            exhaustion, embedding_cache, diagnostics);
    if (result != RequiredPrimeWitnessSearchResult::exhausted) {
        return result;
    }

    // Reset to the required prime before each bounded random
    // subfactor-base twist.
    for (slong trial = 0; trial < effort.random_tries; ++trial) {
        if (!current.set(required_ideal) || !back_multiplier.one()) {
            return RequiredPrimeWitnessSearchResult::failed;
        }

        for (slong index : subfactor_base_schedule->subfactor_base) {
            if (subfactor_base_index_is_excluded(
                        *subfactor_base_schedule, index) ||
                index >= base.length()) {
                return RequiredPrimeWitnessSearchResult::failed;
            }
            const slong exponent = static_cast<slong>(
                    next_relation_random_exponent(random_state));
            if (exponent == 0) {
                continue;
            }
            if (!base.prime(twist_prime, index) ||
                !twist_prime.get_ideal(twist_ideal)) {
                return RequiredPrimeWitnessSearchResult::failed;
            }
            for (slong power = 0; power < exponent; ++power) {
                if (!product.multiply(current, twist_ideal)) {
                    return RequiredPrimeWitnessSearchResult::failed;
                }
                current.swap(product);
            }
        }

        if (!factor_base_honesty_primitive_part(current, back_multiplier) ||
            !factor_base_honesty_reduce_large_ideal(
                    current, back_multiplier, ideal_reduction_precision,
                    diagnostics)) {
            return RequiredPrimeWitnessSearchResult::failed;
        }
        result = enumerate_required_prime_witness(
                base, required_prime, current, back_multiplier, effort,
                exhaustion, embedding_cache, diagnostics);
        if (result != RequiredPrimeWitnessSearchResult::exhausted) {
            return result;
        }
    }

    return RequiredPrimeWitnessSearchResult::exhausted;
}

struct PrimeReductionVisitContext {
    const FactorBase* base = nullptr;
    const PrimeIdeal* prime = nullptr;
    const Order* order = nullptr;
    flint::FmpzMat* basis = nullptr;
    const DiagnosticsContext* diagnostics = nullptr;
    FactorBaseWitnessPredicate predicate = FactorBaseWitnessPredicate::selected;
    bool found = false;
    bool failed = false;
};

int visit_prime_reduction_candidate(const fmpz_mat_t coefficients, void* user) {
    PrimeReductionVisitContext* visit =
            static_cast<PrimeReductionVisitContext*>(user);
    if (visit == nullptr || visit->base == nullptr ||
        visit->prime == nullptr || visit->order == nullptr ||
        visit->basis == nullptr || visit->found || visit->failed) {
        return 0;
    }

    flint::FmpzMat coordinates(1, visit->order->degree());
    coordinates_from_lattice_combination(
            coordinates, flint::FmpzMatConstRef(coefficients),
            flint::FmpzMatConstRef(*visit->basis));

    OrderElement order_element(*visit->order);
    if (!order_element.is_defined() ||
        !order_element.set_coordinates(flint::FmpzMatConstRef(coordinates))) {
        visit->failed = true;
        return 0;
    }

    const bool found = factor_base_principal_witness(
            *visit->base, *visit->prime, order_element, visit->predicate,
            visit->diagnostics);

    if (found) {
        visit->found = true;
        return 0;
    }
    return 1;
}

bool factor_base_reduces_prime_by_ideal_lattice_search(
        const FactorBase& base,
        const PrimeIdeal& prime,
        const Ideal& ideal,
        slong radius,
        FactorBaseWitnessPredicate predicate,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || order->parent() == nullptr ||
        radius <= 0 || !same_order_parent(prime.parent(), order) ||
        !same_order_parent(ideal.parent(), order) || base.contains(prime) ||
        !ideal.has_hnf()) {
        return false;
    }

    flint::FmpzMat hnf(order->degree(), order->degree());
    lat::Lat lattice(order->degree());
    lat::Lat reduced(order->degree());
    if (!ideal.get_hnf(flint::FmpzMatRef(hnf)) ||
        !lattice.set_basis(flint::FmpzMatConstRef(hnf)) ||
        !lattice.lll_reduce(reduced) || reduced.nrows() <= 0) {
        return false;
    }

    flint::FmpzMat basis(reduced.nrows(), order->degree());
    flint::Fmpz max_diag;
    flint::Fmpz gram_ii;
    flint::Arb bound_sq;
    if (!reduced.get_basis(flint::FmpzMatRef(basis))) {
        return false;
    }

    flint::fmpz_zero(flint::FmpzRef(max_diag));
    for (slong i = 0; i < flint::fmpz_mat_nrows(basis); ++i) {
        flint::fmpz_zero(flint::FmpzRef(gram_ii));
        for (slong k = 0; k < flint::fmpz_mat_ncols(basis); ++k) {
            flint::FmpzConstRef entry =
                    flint::fmpz_mat_entry(flint::FmpzMatConstRef(basis), i, k);
            flint::fmpz_addmul(flint::FmpzRef(gram_ii), entry, entry);
        }
        if (flint::fmpz_cmp(flint::FmpzConstRef(gram_ii),
                            flint::FmpzConstRef(max_diag)) > 0) {
            flint::fmpz_set(flint::FmpzRef(max_diag),
                            flint::FmpzConstRef(gram_ii));
        }
    }
    flint::fmpz_mul_ui(flint::FmpzRef(max_diag),
                       flint::FmpzConstRef(max_diag),
                       static_cast<ulong>(reduced.nrows()));

    PrimeReductionVisitContext visit;
    visit.base = &base;
    visit.prime = &prime;
    visit.order = order;
    visit.basis = &basis;
    visit.diagnostics = diagnostics;
    visit.predicate = predicate;
    for (slong tries = 0; tries < 4 && !visit.found && !visit.failed;
         ++tries) {
        flint::arb_set_fmpz(bound_sq, max_diag);
        bool enum_ok = reduced.enum_short_vectors_arb(
                flint::ArbConstRef(bound_sq), radius, 64,
                visit_prime_reduction_candidate, &visit);
        if (!enum_ok && !visit.found && !visit.failed) {
            (void) reduced.enum_short_vectors_arb(
                    flint::ArbConstRef(bound_sq), radius, 256,
                    visit_prime_reduction_candidate, &visit);
        }
        flint::fmpz_mul_ui(flint::FmpzRef(max_diag),
                           flint::FmpzConstRef(max_diag), UWORD(2));
    }

    return visit.found;
}

bool factor_base_reduces_prime_by_lattice_search(const FactorBase& base,
                                                 const PrimeIdeal& prime,
                                                 slong radius,
                                                 FactorBaseWitnessPredicate predicate,
                                                 const DiagnosticsContext*
                                                         diagnostics)
        noexcept {
    const Order* order = base.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order)) {
        return false;
    }

    Ideal ideal(*order);
    return ideal.is_defined() && prime.get_ideal(ideal) &&
           factor_base_reduces_prime_by_ideal_lattice_search(
                   base, prime, ideal, radius,
                   predicate, diagnostics);
}

bool factor_base_reduces_prime_by_twisted_lattice_search(
        const FactorBase& base,
        const PrimeIdeal& prime,
        slong radius,
        slong twist_cap,
        FactorBaseWitnessPredicate predicate,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order) ||
        base.contains(prime) || radius <= 0) {
        return false;
    }

    const slong max_twists =
            base.length() < twist_cap ? base.length() : twist_cap;
    PrimeIdeal twist_prime(*order);
    Ideal prime_ideal(*order);
    Ideal twist_ideal(*order);
    Ideal product(*order);
    if (!twist_prime.is_defined() || !prime_ideal.is_defined() ||
        !twist_ideal.is_defined() || !product.is_defined() ||
        !prime.get_ideal(prime_ideal)) {
        return false;
    }

    for (slong i = 0; i < max_twists; ++i) {
        if (!base.prime(twist_prime, i) ||
            !twist_prime.get_ideal(twist_ideal) ||
            !product.multiply(prime_ideal, twist_ideal)) {
            return false;
        }
        if (factor_base_reduces_prime_by_ideal_lattice_search(
                    base, prime, product, radius,
                    predicate, diagnostics)) {
            return true;
        }
    }
    return false;
}

bool factor_base_reduces_prime_by_principal_search(
        const FactorBase& base,
        const PrimeIdeal& prime,
        slong radius,
        slong twist_cap,
        FactorBaseWitnessPredicate predicate,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !same_order_parent(prime.parent(), order) || base.contains(prime)) {
        return false;
    }
    if (radius < 0) {
        radius = 0;
    }

    flint::Fmpz scalar;
    flint::Fmpz rational_prime;
    if (!prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    if (flint::fmpz_fits_si(flint::FmpzConstRef(rational_prime)) &&
        factor_base_scalar_witness(
                base, prime, flint::FmpzConstRef(rational_prime), predicate,
                diagnostics)) {
        return true;
    }

    const slong rational_prime_si =
            flint::fmpz_fits_si(flint::FmpzConstRef(rational_prime))
                    ? flint::fmpz_get_si(
                              flint::FmpzConstRef(rational_prime))
                    : 0;
    for (slong k = -radius; k <= radius; ++k) {
        if (k == -1 || k == 0 || k == 1 ||
            (flint::fmpz_fits_si(flint::FmpzConstRef(rational_prime)) &&
             k == rational_prime_si)) {
            continue;
        }
        flint::fmpz_set_si(scalar, k);
        if (factor_base_scalar_witness(
                    base, prime, flint::FmpzConstRef(scalar),
                    predicate, diagnostics)) {
            return true;
        }
    }

    return factor_base_reduces_prime_by_lattice_search(
                   base, prime, radius,
                   predicate, diagnostics) ||
           factor_base_reduces_prime_by_twisted_lattice_search(
                   base, prime, radius, twist_cap,
                   predicate, diagnostics);
}

bool factor_base_reduces_prime_by_random_subfactor_base_search(
        const FactorBase& base,
        const PrimeIdeal& prime,
        const SubfactorBaseSchedule* subfactor_base_schedule,
        slong radius,
        slong random_tries,
        ulong& random_state,
        FactorBaseWitnessPredicate predicate,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = base.parent();
    if (order == nullptr || !same_order_parent(prime.parent(), order) ||
        subfactor_base_schedule == nullptr ||
        subfactor_base_schedule->subfactor_base.empty() ||
        radius <= 0) {
        return false;
    }

    PrimeIdeal twist_prime(*order);
    Ideal prime_ideal(*order);
    Ideal current(*order);
    Ideal twist_ideal(*order);
    Ideal product(*order);
    if (!twist_prime.is_defined() || !prime_ideal.is_defined() ||
        !current.is_defined() || !twist_ideal.is_defined() ||
        !product.is_defined() || !prime.get_ideal(prime_ideal)) {
        return false;
    }

    for (slong trial = 0; trial < random_tries; ++trial) {
        if (!current.set(prime_ideal)) {
            return false;
        }

        for (slong index : subfactor_base_schedule->subfactor_base) {
            if (subfactor_base_index_is_excluded(
                        *subfactor_base_schedule, index) ||
                index >= base.length()) {
                return false;
            }

            const slong exponent = static_cast<slong>(
                    next_relation_random_exponent(random_state));
            if (exponent == 0) {
                continue;
            }

            if (!base.prime(twist_prime, index) ||
                !twist_prime.get_ideal(twist_ideal)) {
                return false;
            }
            for (slong k = 0; k < exponent; ++k) {
                if (!product.multiply(current, twist_ideal)) {
                    return false;
                }
                current.swap(product);
            }
        }

        if (factor_base_reduces_prime_by_ideal_lattice_search(
                    base, prime, current, radius,
                    predicate, diagnostics)) {
            return true;
        }
    }
    return false;
}

void log_unwitnessed_required_prime(
        const DiagnosticsContext* diagnostics,
        flint::FmpzConstRef p,
        const PrimeIdeal& prime,
        bool direct_witness_search,
        slong max_stage,
        const FactorBaseWitnessSearchEffort& effort,
        const RequiredPrimeWitnessExhaustion& exhaustion) noexcept {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (!log_enabled(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::detail)) {
        return;
    }
    const long p_value = flint::fmpz_fits_si(p)
            ? static_cast<long>(flint::fmpz_get_si(p))
            : -1L;
    char detail[320];
    if (direct_witness_search) {
        std::snprintf(detail, sizeof(detail),
                      "p=%ld residue_degree=%ld stage=%ld/%ld "
                      "search=t2 random_tries=%ld factor_attempts=%ld "
                      "element_steps=%ld ideals=%ld "
                      "factor_attempt_cap_hits=%ld "
                      "element_step_cap_hits=%ld",
                      p_value, static_cast<long>(prime.residue_degree()),
                      static_cast<long>(effort.stage),
                      static_cast<long>(max_stage),
                      static_cast<long>(effort.random_tries),
                      static_cast<long>(effort.factor_attempts),
                      static_cast<long>(effort.element_steps),
                      static_cast<long>(exhaustion.ideals),
                      static_cast<long>(exhaustion.factor_attempt_cap_hits),
                      static_cast<long>(exhaustion.element_step_cap_hits));
    } else {
        std::snprintf(detail, sizeof(detail),
                      "p=%ld residue_degree=%ld stage=%ld/%ld "
                      "search=lattice radius=%ld twists=%ld "
                      "random_tries=%ld",
                      p_value, static_cast<long>(prime.residue_degree()),
                      static_cast<long>(effort.stage),
                      static_cast<long>(max_stage),
                      static_cast<long>(effort.radius),
                      static_cast<long>(effort.max_twists),
                      static_cast<long>(effort.random_tries));
    }
    log_emit(diagnostics, DiagnosticsModule::class_group, LogLevel::detail,
             __func__, "be_honest witness search exhausted at final caps",
             detail);
#else
    (void) diagnostics;
    (void) p;
    (void) prime;
    (void) direct_witness_search;
    (void) max_stage;
    (void) effort;
    (void) exhaustion;
#endif
}

bool factor_base_honest_for_rational_prime(bool& honest,
                                           const FactorBase& base,
                                           flint::FmpzConstRef p,
                                           flint::FmpzConstRef bound,
                                           const SubfactorBaseSchedule*
                                                   subfactor_base_schedule,
                                           ulong& random_state,
                                           bool use_direct_required_prime_witness,
                                           slong ideal_reduction_precision,
                                           FactorBaseWitnessEscalation
                                                   escalation,
                                           detail::OrderMinkowskiEmbeddingCache*
                                                   embedding_cache,
                                           const DiagnosticsContext*
                                                   diagnostics,
                                           FactorBaseWitnessPredicate predicate,
                                           FactorBaseHonestyScanAudit* audit) noexcept {
    honest = false;
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
        if (prime == nullptr) {
            return false;
        }
        if (base.contains(*prime)) {
            continue;
        }
        if (audit != nullptr) ++audit->witness_targets;
        const slong max_stage =
                escalation == FactorBaseWitnessEscalation::bounded
                        ? kFactorBaseHonestyEscalationStages
                        : 0;
        bool witnessed = false;
        FactorBaseWitnessSearchEffort effort;
        RequiredPrimeWitnessExhaustion exhaustion;
        for (slong stage = 0;
             stage <= max_stage && !witnessed;
             ++stage) {
            effort = factor_base_witness_search_effort(stage);
            if (use_direct_required_prime_witness) {
                exhaustion = RequiredPrimeWitnessExhaustion{};
                const RequiredPrimeWitnessSearchResult result =
                        find_required_prime_witness(
                                base, *prime, subfactor_base_schedule,
                                random_state, ideal_reduction_precision,
                                effort, exhaustion, embedding_cache,
                                diagnostics);
                if (result == RequiredPrimeWitnessSearchResult::failed) {
                    return false;
                }
                witnessed =
                        result == RequiredPrimeWitnessSearchResult::found;
            } else {
                witnessed =
                        factor_base_reduces_prime_by_principal_search(
                                base, *prime, effort.radius,
                                effort.max_twists, predicate,
                                diagnostics) ||
                        factor_base_reduces_prime_by_random_subfactor_base_search(
                                base, *prime, subfactor_base_schedule,
                                effort.radius, effort.random_tries,
                                random_state, predicate, diagnostics);
            }
        }
        if (audit != nullptr && effort.stage > audit->max_search_stage) {
            audit->max_search_stage = effort.stage;
        }
        if (witnessed) {
            if (audit != nullptr) {
                ++audit->witnessed_targets;
                if (effort.stage > 0) ++audit->escalated_witnessed_targets;
            }
            continue;
        }
        log_unwitnessed_required_prime(
                diagnostics, p, *prime, use_direct_required_prime_witness,
                max_stage, effort, exhaustion);
        if (audit != nullptr) {
            FactorBaseHonestyFailure& failure = audit->unwitnessed;
            failure.recorded = true;
            flint::fmpz_set(flint::FmpzRef(failure.p), p);
            failure.residue_degree = prime->residue_degree();
            failure.direct_witness_search = use_direct_required_prime_witness;
            failure.stage = effort.stage;
            failure.max_stage = max_stage;
            failure.radius = effort.radius;
            failure.max_twists = effort.max_twists;
            failure.random_tries = effort.random_tries;
            failure.factor_attempts = effort.factor_attempts;
            failure.element_steps = effort.element_steps;
        }
        return true;
    }

    honest = true;
    return true;
}

}  // namespace

bool factor_base_honesty_check(bool& honest,
                               const FactorBase& base,
                               flint::FmpzConstRef active_bound,
                               flint::FmpzConstRef required_bound,
                               const SubfactorBaseSchedule*
                                       subfactor_base_schedule,
                               ulong random_seed,
                               bool use_direct_required_prime_witness,
                               slong ideal_reduction_precision,
                               FactorBaseWitnessEscalation escalation,
                               const DiagnosticsContext* diagnostics,
                               FactorBaseHonestyScanAudit* audit,
                               FactorBaseWitnessPredicate predicate)
        noexcept {
    honest = false;
    if (audit != nullptr) {
        *audit = FactorBaseHonestyScanAudit{};
        audit->final_random_state = random_seed;
    }
    if (base.parent() == nullptr || flint::fmpz_sgn(active_bound) < 0 ||
        flint::fmpz_sgn(required_bound) < 0) {
        return false;
    }

    ulong random_state = random_seed;
    // Retain one precision-keyed order embedding copy for the scan; every
    // ideal-specific setup stays local.
    detail::OrderMinkowskiEmbeddingCache embedding_cache;
    flint::Fmpz p;
    // The compact relation base can omit higher-residue-degree ideals whose
    // norms still lie below the required bound.  Check the complete rational-
    // prime interval and let base.contains() skip ideals already retained.
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (flint::fmpz_cmp(flint::FmpzConstRef(p), required_bound) <= 0) {
        if (audit != nullptr) {
            ++audit->rational_prime_checks;
            if (flint::fmpz_cmp(flint::FmpzConstRef(p), active_bound) <= 0) {
                ++audit->checks_at_or_below_active_bound;
            }
        }
        bool prime_honest = false;
        const bool checked = factor_base_honest_for_rational_prime(
                    prime_honest, base, flint::FmpzConstRef(p),
                    required_bound, subfactor_base_schedule, random_state,
                    use_direct_required_prime_witness,
                    ideal_reduction_precision, escalation, &embedding_cache,
                    diagnostics, predicate, audit);
        if (audit != nullptr) audit->final_random_state = random_state;
        if (!checked) {
            return false;
        }
        if (!prime_honest) {
            return true;
        }
        flint::fmpz_nextprime(flint::FmpzRef(p), flint::FmpzConstRef(p),
                              true);
    }

    honest = true;
    return true;
}

}  // namespace silex::detail::relation_search
