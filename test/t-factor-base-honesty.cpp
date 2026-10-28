#include "class_group/class_group_internal.hpp"
#include "class_group/factor_base_honesty_internal.hpp"
#include "class_group/factor_base_proof_targets_internal.hpp"
#include "class_group/ideal_minkowski_embedding_internal.hpp"
#include "class_group/ideal_t2_enumeration_internal.hpp"
#include "class_group/relation_candidate_internal.hpp"
#include "class_group/relation_search_driver_internal.hpp"
#include "factor_base/factor_base_internal.hpp"
#include "ideal/ideal_access_internal.hpp"
#include "ideal_factorization/ideal_factorization_internal.hpp"
#include "order/order_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "prime_ideal/prime_ideal_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"
#include "../tools/class_unit_failure_detail.hpp"

#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/lat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct FieldSetup {
    silex::NumberField field;
    silex::Order maximal_order;
};

FieldSetup cubic_x3_minus_2() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

FieldSetup quintic_lower_interval_fixture() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 6);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

FieldSetup truncated_decomposition_counterexample_fixture() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 4);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

FieldSetup ramified_last_fixture() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, -2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -4);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

silex::NumberField quadratic_x2_minus_2() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool first_prime_above(silex::PrimeIdeal& out,
                       const silex::Order& order,
                       slong rational_prime) noexcept {
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    silex::PrimeIdealList primes;
    const silex::PrimeIdeal* first = nullptr;
    return silex::decompose_prime(
                   primes, order, sflint::FmpzConstRef(p)) &&
           primes.size() > 0 && (first = primes.at(0)) != nullptr &&
           out.set(*first);
}

// The witness rule's one exception to "every other factor is in the base":
// the single prime ideal above a rational prime q other than the required
// prime's that the base lacks, when it has e = 1 and the base holds every
// other prime ideal above q (at least one).  Checked here from the complete
// decomposition of q, independently of the library helper.
bool lacking_prime_is_admissible(const silex::PrimeIdeal& factor,
                                 const silex::FactorBase& base,
                                 const silex::PrimeIdeal& required_prime) noexcept {
    sflint::Fmpz q, p;
    if (!factor.rational_prime(sflint::FmpzRef(q)) ||
        !required_prime.rational_prime(sflint::FmpzRef(p)) ||
        sflint::fmpz_equal(sflint::FmpzConstRef(q), sflint::FmpzConstRef(p)) ||
        factor.ramification_index() != 1) {
        return false;
    }
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, *base.parent(), sflint::FmpzConstRef(q))) {
        return false;
    }
    slong inside = 0;
    slong outside = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr) return false;
        if (base.contains(*prime)) {
            ++inside;
        } else {
            ++outside;
        }
    }
    return inside >= 1 && outside == 1;
}

// With strict set, the rule without the lacking-prime exception.
bool full_factorization_required_prime_result(
        bool& matches,
        const silex::Ideal& ideal,
        const silex::FactorBase& base,
        const silex::PrimeIdeal& required_prime,
        bool strict = false) noexcept {
    matches = false;
    const silex::Order* order = ideal.parent();
    if (order == nullptr ||
        !silex::same_order_parent(base.parent(), order) ||
        !silex::same_order_parent(required_prime.parent(), order) ||
        base.contains(required_prime)) {
        return false;
    }

    silex::IdealFactorization factorization(*order);
    silex::PrimeIdeal factor(*order);
    if (!factorization.factor(ideal)) {
        return false;
    }

    bool other_factors_are_in_base = true;
    slong required_exponent = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = 0;
        if (!factorization.prime(factor, i) ||
            !factorization.exponent(exponent, i)) {
            return false;
        }
        // Matched as ideals: the factorization may record residue data
        // other than required_prime's for the same prime.
        if (silex::detail::PrimeIdealAccess::same_prime_ideal(
                    factor, required_prime)) {
            required_exponent += exponent;
        } else if (!base.contains(factor) &&
                   (strict ||
                    !lacking_prime_is_admissible(factor, base, required_prime))) {
            other_factors_are_in_base = false;
        }
    }
    matches = other_factors_are_in_base && required_exponent == 1;
    return true;
}

bool prime_has_ef(const silex::PrimeIdealList& primes,
                  slong index,
                  slong ramification_index,
                  slong residue_degree) noexcept {
    const silex::PrimeIdeal* prime = primes.at(index);
    return prime != nullptr &&
           prime->ramification_index() == ramification_index &&
           prime->residue_degree() == residue_degree;
}

int test_proof_targets_require_complete_decomposition() {
    FieldSetup setup = truncated_decomposition_counterexample_fixture();
    sflint::Fmpz two;
    sflint::Fmpz truncated_bound;
    sflint::Fmpz complete_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(two), UWORD(2));
    sflint::fmpz_set_ui(sflint::FmpzRef(truncated_bound), UWORD(2));
    sflint::fmpz_set_ui(sflint::FmpzRef(complete_bound), UWORD(4));

    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, setup.maximal_order,
                                sflint::FmpzConstRef(two)) ||
        primes.size() != 3) {
        return 1;
    }

    std::vector<slong> degree_one_indices;
    slong degree_two_index = -1;
    for (slong i = 0; i < primes.size(); ++i) {
        if (prime_has_ef(primes, i, 1, 1)) {
            degree_one_indices.push_back(i);
        } else if (prime_has_ef(primes, i, 1, 2) &&
                   degree_two_index < 0) {
            degree_two_index = i;
        } else {
            return 1;
        }
    }
    if (degree_one_indices.size() != 2 || degree_two_index < 0) {
        return 1;
    }

    std::vector<slong> targets{99};
    if (!silex::detail::select_factor_base_proof_targets(
                targets, primes, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(truncated_bound)) ||
        targets != degree_one_indices) {
        return 1;
    }

    // Once the cutoff contains the full decomposition, any one unramified
    // ideal may be discharged by the principal factorization of 2.  Put a
    // degree-one ideal last to show that this does not require choosing the
    // maximum-residue-degree ideal.
    silex::PrimeIdealList controlled(setup.maximal_order, 3);
    if (!controlled.is_defined() ||
        !controlled.at(0)->set(*primes.at(degree_one_indices[0])) ||
        !controlled.at(1)->set(*primes.at(degree_two_index)) ||
        !controlled.at(2)->set(*primes.at(degree_one_indices[1]))) {
        return 1;
    }
    if (!silex::detail::select_factor_base_proof_targets(
                targets, controlled, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(truncated_bound)) ||
        targets.size() != 2 || targets[0] != 0 || targets[1] != 2) {
        return 1;
    }
    if (!silex::detail::select_factor_base_proof_targets(
                targets, controlled, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(complete_bound)) ||
        targets.size() != 2 || targets[0] != 0 || targets[1] != 1) {
        return 1;
    }

    sflint::Fmpz no_target_bound;
    sflint::fmpz_one(sflint::FmpzRef(no_target_bound));
    if (!silex::detail::select_factor_base_proof_targets(
                targets, controlled, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(no_target_bound)) ||
        !targets.empty()) {
        return 1;
    }

    // An inert prime is discharged: its only factor is pO_K, which is
    // principal, so no factor-base membership is required.
    silex::NumberField singleton_field = quadratic_x2_minus_2();
    silex::Order singleton_order =
            silex::test::equation_order(singleton_field);
    sflint::Fmpz three;
    sflint::Fmpz nine;
    sflint::fmpz_set_ui(sflint::FmpzRef(three), UWORD(3));
    sflint::fmpz_set_ui(sflint::FmpzRef(nine), UWORD(9));
    silex::PrimeIdealList singleton;
    if (!singleton_order.is_maximal() ||
        !silex::decompose_prime(singleton, singleton_order,
                                sflint::FmpzConstRef(three)) ||
        singleton.size() != 1 || !prime_has_ef(singleton, 0, 1, 2) ||
        !silex::detail::select_factor_base_proof_targets(
                targets, singleton, sflint::FmpzConstRef(three),
                sflint::FmpzConstRef(nine)) ||
        !targets.empty()) {
        return 1;
    }

    // A totally ramified singleton, (2) = P^2 in Q(sqrt(2)), stays a target:
    // the principal relation only gives 2[P] = 0.
    sflint::Fmpz ramified_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(ramified_bound), UWORD(2));
    if (!silex::decompose_prime(singleton, singleton_order,
                                sflint::FmpzConstRef(two)) ||
        singleton.size() != 1 || !prime_has_ef(singleton, 0, 2, 1) ||
        !silex::detail::select_factor_base_proof_targets(
                targets, singleton, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(ramified_bound)) ||
        targets.size() != 1 || targets[0] != 0) {
        return 1;
    }

    FieldSetup ramified_setup = ramified_last_fixture();
    if (!silex::decompose_prime(primes, ramified_setup.maximal_order,
                                sflint::FmpzConstRef(two)) ||
        primes.size() != 2) {
        return 1;
    }
    slong unramified_index = -1;
    slong ramified_index = -1;
    for (slong i = 0; i < primes.size(); ++i) {
        if (prime_has_ef(primes, i, 1, 1)) {
            unramified_index = i;
        } else if (prime_has_ef(primes, i, 2, 2)) {
            ramified_index = i;
        }
    }
    silex::PrimeIdealList ramified_last(ramified_setup.maximal_order, 2);
    if (unramified_index < 0 || ramified_index < 0 ||
        !ramified_last.is_defined() ||
        !ramified_last.at(0)->set(*primes.at(unramified_index)) ||
        !ramified_last.at(1)->set(*primes.at(ramified_index))) {
        return 1;
    }
    if (!silex::detail::select_factor_base_proof_targets(
                targets, ramified_last, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(truncated_bound)) ||
        targets.size() != 1 || targets[0] != 0 ||
        !silex::detail::select_factor_base_proof_targets(
                targets, ramified_last, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(complete_bound)) ||
        targets.size() != 2 || targets[0] != 0 || targets[1] != 1) {
        return 1;
    }

    sflint::Fmpz invalid_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(invalid_bound), -1);
    targets = {7, 11};
    if (silex::detail::select_factor_base_proof_targets(
                targets, ramified_last, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(invalid_bound)) ||
        targets.size() != 2 || targets[0] != 7 || targets[1] != 11) {
        return 1;
    }

    sflint::Fmpz zero;
    sflint::fmpz_zero(sflint::FmpzRef(zero));
    if (silex::detail::select_factor_base_proof_targets(
                targets, ramified_last, sflint::FmpzConstRef(zero),
                sflint::FmpzConstRef(complete_bound)) ||
        targets.size() != 2 || targets[0] != 7 || targets[1] != 11) {
        return 1;
    }

    silex::PrimeIdealList undefined_prime(ramified_setup.maximal_order, 1);
    if (!undefined_prime.is_defined() ||
        silex::detail::select_factor_base_proof_targets(
                targets, undefined_prime, sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(complete_bound)) ||
        targets.size() != 2 || targets[0] != 7 || targets[1] != 11) {
        return 1;
    }
    return 0;
}

int test_required_prime_predicates_match_full_factorization() {
    silex::NumberField field = quadratic_x2_minus_2();
    silex::Order order = silex::test::equation_order(field);
    if (!order.is_maximal()) {
        return 1;
    }

    sflint::Fmpz bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), UWORD(7));
    silex::FactorBase base(order);
    silex::PrimeIdeal required_prime(order);
    silex::PrimeIdeal other_prime(order);
    if (!base.is_defined() ||
        !base.build(sflint::FmpzConstRef(bound)) ||
        !first_prime_above(required_prime, order, 11) ||
        !first_prime_above(other_prime, order, 13) ||
        base.contains(required_prime) || base.contains(other_prime)) {
        return 1;
    }

    struct Case {
        slong value;
        bool expected;
    };
    constexpr std::array cases{
            Case{11, true}, Case{22, true}, Case{121, false},
            Case{143, false}, Case{13, false}};

    silex::OrderElement generator(order);
    silex::Ideal principal(order);
    if (!generator.is_defined() || !principal.is_defined()) {
        return 1;
    }
    for (const Case& test_case : cases) {
        if (!generator.set_si(test_case.value) ||
            !principal.set_principal(generator)) {
            return 1;
        }

        bool ideal_matches = !test_case.expected;
        bool element_matches = !test_case.expected;
        if (!silex::detail::ideal_factor_over_base_with_required_prime(
                    ideal_matches, principal, base, required_prime) ||
            !silex::detail::order_element_factor_over_base_with_required_prime(
                    element_matches, generator, base, required_prime)) {
            return 1;
        }
        bool full_matches = !test_case.expected;
        if (!full_factorization_required_prime_result(
                    full_matches, principal, base, required_prime) ||
            ideal_matches != test_case.expected ||
            element_matches != test_case.expected ||
            ideal_matches != full_matches || element_matches != full_matches) {
            return 1;
        }
    }

    // In Z[sqrt(2)], 1 + sqrt(2) is a unit. This genuinely non-scalar
    // generator therefore has the same required-prime support as (11), which
    // the independent full ideal factorization verifies below.
    sflint::FmpzMat coordinates(1, order.degree());
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 0),
            11);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 1),
            11);
    if (!generator.set_coordinates(sflint::FmpzMatConstRef(coordinates)) ||
        !principal.set_principal(generator)) {
        return 1;
    }
    bool ideal_matches = false;
    bool element_matches = false;
    bool full_matches = false;
    if (!silex::detail::ideal_factor_over_base_with_required_prime(
                ideal_matches, principal, base, required_prime) ||
        !silex::detail::order_element_factor_over_base_with_required_prime(
                element_matches, generator, base, required_prime) ||
        !full_factorization_required_prime_result(
                full_matches, principal, base, required_prime) ||
        !ideal_matches || !element_matches || !full_matches) {
        return 1;
    }
    return 0;
}

// Source: buch2.c:can_factor rejects a norm whose part prime to the factor
// base primes is not 1 (Z_ppo against prodZ) before calling absZ_factor. The
// cofactor below is the product of the Mersenne primes 2^61 - 1 and 2^89 - 1.
// The test asserts that the prefilter rejected it, through the internal audit
// in every build and through the reject_cofactor event in profiling builds.
// A smooth norm with large exponents still reaches the full
// classification and matches.
struct RejectCofactorCounter {
    slong count = 0;
};

void reject_cofactor_callback(void* user,
                              silex::DiagnosticsModule,
                              silex::ProfileEvent event,
                              const char*,
                              const char* label) noexcept {
    if (label != nullptr && event == silex::ProfileEvent::event &&
        std::strcmp(label, "ideal.required_prime_direct.reject_cofactor") ==
                0) {
        ++static_cast<RejectCofactorCounter*>(user)->count;
    }
}

int test_required_prime_direct_rejects_nonsmooth_cofactor() {
    silex::NumberField field = quadratic_x2_minus_2();
    silex::Order order = silex::test::equation_order(field);
    sflint::Fmpz bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), UWORD(7));
    silex::FactorBase base(order);
    silex::PrimeIdeal required_prime(order);
    if (!order.is_maximal() || !base.is_defined() ||
        !base.build(sflint::FmpzConstRef(bound)) ||
        !first_prime_above(required_prime, order, 11) ||
        base.contains(required_prime)) {
        return 1;
    }

    sflint::Fmpz two, m61, m89, value;
    sflint::fmpz_set_ui(sflint::FmpzRef(two), UWORD(2));
    sflint::fmpz_pow_ui(sflint::FmpzRef(m61), sflint::FmpzConstRef(two), 61);
    sflint::fmpz_sub_ui(sflint::FmpzRef(m61), sflint::FmpzConstRef(m61), 1);
    sflint::fmpz_pow_ui(sflint::FmpzRef(m89), sflint::FmpzConstRef(two), 89);
    sflint::fmpz_sub_ui(sflint::FmpzRef(m89), sflint::FmpzConstRef(m89), 1);

    sflint::FmpzMat coordinates(1, order.degree());
    silex::OrderElement generator(order);
    silex::Ideal principal(order);
    if (!generator.is_defined() || !principal.is_defined()) {
        return 1;
    }

    // 11 * 2^5 * 7^3 * (2^61 - 1) * (2^89 - 1): required prime and base
    // primes present, but the cofactor is not 1, so it is not a witness. The
    // cofactor primes are small enough that the full factorization returns
    // quickly, so a missing prefilter shows up as an unset audit flag.
    sflint::fmpz_mul(sflint::FmpzRef(value), sflint::FmpzConstRef(m61),
                     sflint::FmpzConstRef(m89));
    sflint::fmpz_mul_ui(sflint::FmpzRef(value), sflint::FmpzConstRef(value),
                        UWORD(11) * UWORD(32) * UWORD(343));
    sflint::fmpz_set(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 0),
            sflint::FmpzConstRef(value));
    bool matches = true;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    RejectCofactorCounter reject_counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::ideal),
            reject_cofactor_callback, &reject_counter);
    silex::detail::RequiredPrimeDirectAudit audit;
    if (!generator.set_coordinates(sflint::FmpzMatConstRef(coordinates)) ||
        !silex::detail::order_element_factor_over_base_with_required_prime(
                matches, generator, base, required_prime, &diagnostics,
                &audit) ||
        matches || !audit.rejected_by_cofactor) {
        return 1;
    }
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    if (reject_counter.count != 1) {
        return 1;
    }
#endif

    // 11 * 2^5 * 7^3 is smooth over the base, so the prefilter passes it on
    // and the full classification accepts it, as the ideal path does.
    sflint::fmpz_set_ui(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 0),
            UWORD(11) * UWORD(32) * UWORD(343));
    matches = false;
    bool ideal_matches = false;
    audit = silex::detail::RequiredPrimeDirectAudit{};
    if (!generator.set_coordinates(sflint::FmpzMatConstRef(coordinates)) ||
        !principal.set_principal(generator) ||
        !silex::detail::order_element_factor_over_base_with_required_prime(
                matches, generator, base, required_prime, nullptr, &audit) ||
        audit.rejected_by_cofactor ||
        !silex::detail::ideal_factor_over_base_with_required_prime(
                ideal_matches, principal, base, required_prime) ||
        !matches || !ideal_matches) {
        return 1;
    }
    return 0;
}

FieldSetup maximal_order_setup(const std::vector<slong>& coefficients) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

// Every nonzero coordinate vector in [-radius, radius]^n: the prefilter norm
// equals Element::norm, and the prefilter rejects exactly the elements the
// direct order-element test rejects by their norm cofactor. The box must hold
// witnesses, which a prefilter that dropped the required prime would reject.
int check_norm_prefilter_against_direct(const FieldSetup& setup,
                                        slong base_bound,
                                        slong required_rational_prime,
                                        slong radius,
                                        bool shared_rational_prime = false) {
    const silex::Order& order = setup.maximal_order;
    const slong n = order.degree();
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), base_bound);
    silex::FactorBase base(order);
    silex::PrimeIdeal required(order);
    if (!order.is_maximal() || !base.is_defined()) {
        return 1;
    }
    if (shared_rational_prime) {
        // Norm-bounded base: it holds the degree-one prime above the required
        // rational prime but not the required residue-degree-two prime, so
        // that rational prime occurs twice in the smooth product.
        sflint::Fmpz p;
        sflint::fmpz_set_si(sflint::FmpzRef(p), required_rational_prime);
        silex::PrimeIdealList primes;
        const silex::PrimeIdeal* above = nullptr;
        bool found = false;
        if (!base.build_prime_ideal_norm_bounded(
                    sflint::FmpzConstRef(bound)) ||
            !silex::decompose_prime(primes, order, sflint::FmpzConstRef(p))) {
            return 1;
        }
        for (slong i = 0; i < primes.size(); ++i) {
            above = primes.at(i);
            if (above != nullptr && above->residue_degree() == 2) {
                found = required.set(*above);
            }
        }
        if (!found || !base.contains(*primes.at(0)) ||
            primes.at(0)->residue_degree() != 1) {
            return 1;
        }
    } else if (!base.build(sflint::FmpzConstRef(bound)) ||
               !first_prime_above(required, order, required_rational_prime)) {
        return 1;
    }
    if (base.contains(required)) {
        return 1;
    }
    silex::detail::RequiredPrimeNormPrefilter prefilter;
    if (!silex::detail::required_prime_norm_prefilter_init(prefilter, base,
                                                          required) ||
        !prefilter.ready || prefilter.degree != n) {
        return 1;
    }
    if (shared_rational_prime) {
        sflint::Fmpz twice;
        sflint::fmpz_set_si(sflint::FmpzRef(twice),
                            required_rational_prime * required_rational_prime);
        if (!fmpz_divisible(prefilter.smooth_product.raw(), twice.raw())) {
            return 1;
        }
    }

    sflint::FmpzMat coordinates(1, n);
    silex::OrderElement element(order);
    silex::Element ambient(*order.parent());
    sflint::Fmpq field_norm;
    sflint::Fmpz prefilter_norm;
    std::vector<slong> digits(static_cast<std::size_t>(n), -radius);
    slong rejected_count = 0, kept_count = 0, witness_count = 0;
    for (;;) {
        bool nonzero = false;
        for (slong i = 0; i < n; ++i) {
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0,
                                           i),
                    digits[static_cast<std::size_t>(i)]);
            nonzero = nonzero || digits[static_cast<std::size_t>(i)] != 0;
        }
        if (nonzero) {
            bool matches = true;
            bool rejected = false;
            silex::detail::RequiredPrimeDirectAudit audit;
            if (!element.set_coordinates(
                        sflint::FmpzMatConstRef(coordinates)) ||
                !element.get_element(ambient) ||
                !ambient.norm(sflint::FmpqRef(field_norm)) ||
                !silex::detail::required_prime_norm_prefilter_norm(
                        sflint::FmpzRef(prefilter_norm), prefilter,
                        sflint::FmpzMatConstRef(coordinates)) ||
                fmpz_is_one(fmpq_denref(field_norm.raw())) == 0 ||
                fmpz_equal(fmpq_numref(field_norm.raw()),
                           prefilter_norm.raw()) == 0 ||
                !silex::detail::required_prime_norm_prefilter_rejects(
                        rejected, prefilter,
                        sflint::FmpzMatConstRef(coordinates)) ||
                !silex::detail::order_element_factor_over_base_with_required_prime(
                        matches, element, base, required, nullptr, &audit) ||
                rejected != audit.rejected_by_cofactor ||
                (rejected && matches)) {
                return 1;
            }
            rejected ? ++rejected_count : ++kept_count;
            witness_count += matches ? 1 : 0;
        }
        slong i = 0;
        while (i < n && digits[static_cast<std::size_t>(i)] == radius) {
            digits[static_cast<std::size_t>(i)] = -radius;
            ++i;
        }
        if (i == n) {
            break;
        }
        ++digits[static_cast<std::size_t>(i)];
    }
    return rejected_count > 0 && kept_count > rejected_count / 1000 &&
                   witness_count > 0
           ? 0
           : 1;
}

int test_required_prime_norm_prefilter_matches_direct() {
    // Degree 2 with a non-equation maximal order basis, degree 3, and two
    // quartics, one of them Q(sqrt 2, sqrt 3).
    const FieldSetup quadratic = maximal_order_setup({7, 0, 1});
    const FieldSetup cubic = cubic_x3_minus_2();
    const FieldSetup biquadratic = maximal_order_setup({1, 0, -10, 0, 1});
    const FieldSetup quartic = truncated_decomposition_counterexample_fixture();
    if (check_norm_prefilter_against_direct(quadratic, 7, 11, 20) != 0 ||
        check_norm_prefilter_against_direct(cubic, 7, 11, 5) != 0 ||
        check_norm_prefilter_against_direct(cubic, 5, 5, 5, true) != 0 ||
        check_norm_prefilter_against_direct(biquadratic, 7, 23, 3) != 0 ||
        check_norm_prefilter_against_direct(quartic, 7, 11, 3) != 0) {
        return 1;
    }

    // Not ready outside degrees 2 to 4, or when the required prime is in the
    // base: those keep the direct path.
    const FieldSetup quintic = ramified_last_fixture();
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 7);
    silex::FactorBase quintic_base(quintic.maximal_order);
    silex::PrimeIdeal quintic_required(quintic.maximal_order);
    silex::FactorBase quadratic_base(quadratic.maximal_order);
    silex::PrimeIdeal in_base(quadratic.maximal_order);
    silex::detail::RequiredPrimeNormPrefilter prefilter;
    if (quintic.maximal_order.degree() != 5 ||
        !quintic_base.build(sflint::FmpzConstRef(bound)) ||
        !first_prime_above(quintic_required, quintic.maximal_order, 11) ||
        quintic_base.contains(quintic_required) ||
        silex::detail::required_prime_norm_prefilter_init(
                prefilter, quintic_base, quintic_required) ||
        prefilter.ready || !quadratic_base.build(sflint::FmpzConstRef(bound)) ||
        !first_prime_above(in_base, quadratic.maximal_order, 2) ||
        !quadratic_base.contains(in_base) ||
        silex::detail::required_prime_norm_prefilter_init(
                prefilter, quadratic_base, in_base) ||
        prefilter.ready) {
        return 1;
    }
    return 0;
}

// The degree-2 word path of the norm prefilter against its fmpz path on one
// quadratic field: the word path applies exactly when the table and the
// coordinates are within the word limit and |N| < 2^FLINT_BITS, and then
// decides as the fmpz path does. Extra rows are coordinate pairs to add to a
// small grid.
int check_norm_prefilter_word_path(
        const FieldSetup& setup,
        slong base_bound,
        slong required_rational_prime,
        bool expect_word_table,
        const std::vector<std::array<sflint::Fmpz, 2>>& extra,
        const std::vector<std::array<sflint::Fmpz, 2>>& field_extra,
        slong& word_rejected,
        slong& word_kept,
        slong& fallbacks) {
    const silex::Order& order = setup.maximal_order;
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), base_bound);
    silex::FactorBase base(order);
    silex::PrimeIdeal required(order);
    silex::detail::RequiredPrimeNormPrefilter prefilter;
    if (order.degree() != 2 || !order.is_maximal() ||
        !base.build(sflint::FmpzConstRef(bound)) ||
        !first_prime_above(required, order, required_rational_prime) ||
        base.contains(required) ||
        !silex::detail::required_prime_norm_prefilter_init(prefilter, base,
                                                          required) ||
        !prefilter.ready ||
        prefilter.word_table_ready != expect_word_table) {
        return 1;
    }
    sflint::Fmpz word_bound;  // 2^FLINT_BITS
    sflint::fmpz_one(sflint::FmpzRef(word_bound));
    fmpz_mul_2exp(word_bound.raw(), word_bound.raw(), FLINT_BITS);
    sflint::Fmpz limit;
    sflint::fmpz_set_si(sflint::FmpzRef(limit),
                        silex::detail::required_prime_norm_prefilter_word_limit);

    sflint::FmpzMat coordinates(1, 2);
    sflint::Fmpz norm;
    const auto check = [&](const fmpz* a0, const fmpz* a1) -> bool {
        fmpz_set(fmpz_mat_entry(coordinates.raw(), 0, 0), a0);
        fmpz_set(fmpz_mat_entry(coordinates.raw(), 0, 1), a1);
        bool fmpz_rejected = true;
        bool word_rejected_now = true;
        bool combined_rejected = true;
        if (!silex::detail::required_prime_norm_prefilter_norm(
                    sflint::FmpzRef(norm), prefilter,
                    sflint::FmpzMatConstRef(coordinates)) ||
            !silex::detail::required_prime_norm_prefilter_rejects_fmpz(
                    fmpz_rejected, prefilter,
                    sflint::FmpzMatConstRef(coordinates)) ||
            !silex::detail::required_prime_norm_prefilter_rejects(
                    combined_rejected, prefilter,
                    sflint::FmpzMatConstRef(coordinates)) ||
            combined_rejected != fmpz_rejected) {
            return false;
        }
        const bool expect_word =
                expect_word_table && fmpz_cmpabs(a0, limit.raw()) <= 0 &&
                fmpz_cmpabs(a1, limit.raw()) <= 0 &&
                fmpz_cmpabs(norm.raw(), word_bound.raw()) < 0;
        const bool applied = silex::detail::required_prime_norm_prefilter_rejects_word(
                word_rejected_now, prefilter,
                sflint::FmpzMatConstRef(coordinates));
        if (applied != expect_word || (!applied && word_rejected_now) ||
            (applied && word_rejected_now != fmpz_rejected)) {
            return false;
        }
        if (!applied) {
            ++fallbacks;
        } else if (word_rejected_now) {
            ++word_rejected;
        } else {
            ++word_kept;
        }
        return true;
    };

    sflint::Fmpz x, y;
    for (slong i = -40; i <= 40; ++i) {
        for (slong j = -40; j <= 40; ++j) {
            sflint::fmpz_set_si(sflint::FmpzRef(x), i);
            sflint::fmpz_set_si(sflint::FmpzRef(y), j);
            if (!check(x.raw(), y.raw())) {
                return 1;
            }
        }
    }
    for (const auto* rows : {&extra, &field_extra}) {
    for (const auto& row : *rows) {
        // Each pair also with signs flipped and coordinates swapped.
        for (int sign = 0; sign < 4; ++sign) {
            fmpz_set(x.raw(), row[0].raw());
            fmpz_set(y.raw(), row[1].raw());
            if ((sign & 1) != 0) {
                fmpz_neg(x.raw(), x.raw());
            }
            if ((sign & 2) != 0) {
                fmpz_neg(y.raw(), y.raw());
            }
            if (!check(x.raw(), y.raw()) || !check(y.raw(), x.raw())) {
                return 1;
            }
        }
    }
    }
    return 0;
}

int test_required_prime_norm_prefilter_word_path_matches_fmpz() {
    const slong limit = silex::detail::required_prime_norm_prefilter_word_limit;
    if (limit != (WORD(1) << (FLINT_BITS / 2 - 1)) - 1) {
        return 1;
    }
    std::vector<std::array<sflint::Fmpz, 2>> extra;
    const auto add = [&extra](const fmpz* a0, const fmpz* a1) {
        extra.emplace_back();
        fmpz_set(extra.back()[0].raw(), a0);
        fmpz_set(extra.back()[1].raw(), a1);
    };
    const auto add_si = [&add](slong a0, slong a1) {
        sflint::Fmpz x, y;
        sflint::fmpz_set_si(sflint::FmpzRef(x), a0);
        sflint::fmpz_set_si(sflint::FmpzRef(y), a1);
        add(x.raw(), y.raw());
    };
    // Coordinates at and just past the limit, and far past it.
    sflint::Fmpz past, big;
    sflint::fmpz_set_si(sflint::FmpzRef(past), limit);
    fmpz_add_ui(past.raw(), past.raw(), 1);
    sflint::fmpz_one(sflint::FmpzRef(big));
    fmpz_mul_2exp(big.raw(), big.raw(), 3 * FLINT_BITS);
    for (slong small : {slong(0), slong(1), slong(3), slong(17)}) {
        add_si(limit, small);
        add_si(limit - 1, small);
        sflint::Fmpz s;
        sflint::fmpz_set_si(sflint::FmpzRef(s), small);
        add(past.raw(), s.raw());
        add(big.raw(), s.raw());
    }
    add_si(limit, limit);
    add(past.raw(), past.raw());
    // Smooth coordinates with norms near the word bound.
    add_si(WORD(1) << 30, 0);
    add_si(WORD(1162261467), 0);  // 3^19
    add_si(WORD(1162261467), WORD(1) << 30);
    add_si(WORD(30) * 30 * 30 * 30 * 30 * 30, 0);

    slong word_rejected = 0, word_kept = 0, fallbacks = 0;
    // Table entries within the limit: x^2 + 1009, the reported row
    // x^2 - x - 25000000, and x^2 + x + (2^31 - 1), whose table holds
    // -(2^31 - 1), exactly the limit on 64-bit words.
    const FieldSetup imaginary = maximal_order_setup({1009, 0, 1});
    std::vector<std::array<sflint::Fmpz, 2>> imaginary_extra;
    const std::vector<std::array<sflint::Fmpz, 2>> no_extra;
    {
        // |N| at 2^FLINT_BITS - 1 or below and at 2^FLINT_BITS or above:
        // for a fixed second coordinate the norm grows with the first.
        sflint::Fmpz word_bound;
        sflint::fmpz_one(sflint::FmpzRef(word_bound));
        fmpz_mul_2exp(word_bound.raw(), word_bound.raw(), FLINT_BITS);
        sflint::Fmpz bound_value;
        sflint::fmpz_set_si(sflint::FmpzRef(bound_value), 50);
        silex::FactorBase base(imaginary.maximal_order);
        silex::PrimeIdeal required(imaginary.maximal_order);
        silex::detail::RequiredPrimeNormPrefilter prefilter;
        if (!base.build(sflint::FmpzConstRef(bound_value)) ||
            !first_prime_above(required, imaginary.maximal_order, 53) ||
            !silex::detail::required_prime_norm_prefilter_init(prefilter, base,
                                                              required)) {
            return 1;
        }
        sflint::FmpzMat coordinates(1, 2);
        sflint::Fmpz norm;
        const auto norm_below = [&](slong a0, slong a1) -> int {
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 0),
                    a0);
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 1),
                    a1);
            if (!silex::detail::required_prime_norm_prefilter_norm(
                        sflint::FmpzRef(norm), prefilter,
                        sflint::FmpzMatConstRef(coordinates))) {
                return -1;
            }
            return fmpz_cmpabs(norm.raw(), word_bound.raw()) < 0 ? 1 : 0;
        };
        for (slong second : {WORD(120000000), WORD(125000001)}) {
            slong lo = 0, hi = limit;
            if (norm_below(lo, second) != 1 || norm_below(hi, second) != 0) {
                return 1;
            }
            while (hi - lo > 1) {
                const slong mid = lo + (hi - lo) / 2;
                const int below = norm_below(mid, second);
                if (below < 0) {
                    return 1;
                }
                (below == 1 ? lo : hi) = mid;
            }
            sflint::Fmpz x, y;
            for (slong first : {lo - 1, lo, hi, hi + 1}) {
                sflint::fmpz_set_si(sflint::FmpzRef(x), first);
                sflint::fmpz_set_si(sflint::FmpzRef(y), second);
                imaginary_extra.emplace_back();
                fmpz_set(imaginary_extra.back()[0].raw(), x.raw());
                fmpz_set(imaginary_extra.back()[1].raw(), y.raw());
            }
        }
    }
    // Negative norms at the 2^FLINT_BITS bound: in x^2 - x - 4 the norm
    // a0^2 + a0 a1 - 4 a1^2 is negative, with |N| just below 2^64 at
    // (0, 2^31 - 1) and growing past it as a0 decreases toward -a1 / 2.
    const FieldSetup negative = maximal_order_setup({-4, -1, 1});
    std::vector<std::array<sflint::Fmpz, 2>> negative_extra;
    {
        sflint::Fmpz word_bound;
        sflint::fmpz_one(sflint::FmpzRef(word_bound));
        fmpz_mul_2exp(word_bound.raw(), word_bound.raw(), FLINT_BITS);
        sflint::Fmpz bound_value;
        sflint::fmpz_set_si(sflint::FmpzRef(bound_value), 50);
        silex::FactorBase base(negative.maximal_order);
        silex::PrimeIdeal required(negative.maximal_order);
        silex::detail::RequiredPrimeNormPrefilter prefilter;
        if (!base.build(sflint::FmpzConstRef(bound_value)) ||
            !first_prime_above(required, negative.maximal_order, 53) ||
            !silex::detail::required_prime_norm_prefilter_init(prefilter, base,
                                                              required)) {
            return 1;
        }
        sflint::FmpzMat coordinates(1, 2);
        sflint::Fmpz norm;
        const slong second = limit;
        // 1 when the norm is negative with |N| < 2^64, 0 when |N| >= 2^64.
        const auto below = [&](slong a0) -> int {
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 0),
                    a0);
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 1),
                    second);
            if (!silex::detail::required_prime_norm_prefilter_norm(
                        sflint::FmpzRef(norm), prefilter,
                        sflint::FmpzMatConstRef(coordinates)) ||
                fmpz_sgn(norm.raw()) >= 0) {
                return -1;
            }
            return fmpz_cmpabs(norm.raw(), word_bound.raw()) < 0 ? 1 : 0;
        };
        slong lo = -(limit / 2), hi = 0;  // |N| >= 2^64 at lo, below at hi
        if (below(lo) != 0 || below(hi) != 1) {
            return 1;
        }
        while (hi - lo > 1) {
            const slong mid = lo + (hi - lo) / 2;
            const int b = below(mid);
            if (b < 0) {
                return 1;
            }
            (b == 1 ? hi : lo) = mid;
        }
        for (slong first : {lo - 1, lo, hi, hi + 1}) {
            negative_extra.emplace_back();
            sflint::fmpz_set_si(sflint::FmpzRef(negative_extra.back()[0]),
                                first);
            sflint::fmpz_set_si(sflint::FmpzRef(negative_extra.back()[1]),
                                second);
        }
    }
    const FieldSetup row = maximal_order_setup({-25000000, -1, 1});
    const FieldSetup at_limit = maximal_order_setup({2147483647, 1, 1});
    if (check_norm_prefilter_word_path(imaginary, 50, 53, true,
                                       extra, imaginary_extra, word_rejected,
                                       word_kept, fallbacks) != 0 ||
        check_norm_prefilter_word_path(negative, 50, 53, true, extra,
                                       negative_extra, word_rejected,
                                       word_kept, fallbacks) != 0 ||
        check_norm_prefilter_word_path(row, 50, 53, true, extra, no_extra,
                                       word_rejected, word_kept,
                                       fallbacks) != 0 ||
        check_norm_prefilter_word_path(at_limit, 50, 53, true, extra, no_extra,
                                       word_rejected, word_kept,
                                       fallbacks) != 0) {
        return 1;
    }
    // A table entry of 2^31, just past the limit: every candidate falls back.
    const slong before = fallbacks;
    const FieldSetup past_limit = maximal_order_setup({2147483648, 1, 1});
    slong past_rejected = 0, past_kept = 0;
    if (check_norm_prefilter_word_path(past_limit, 50, 53, false, extra,
                                       no_extra,
                                       past_rejected, past_kept,
                                       fallbacks) != 0 ||
        past_rejected != 0 || past_kept != 0 || fallbacks == before) {
        return 1;
    }
    return word_rejected > 0 && word_kept > 0 && before > 0 ? 0 : 1;
}

int test_selected_witness_edge_cases() {
    namespace search = silex::detail::relation_search;
    using Mode = search::FactorBaseWitnessPredicate;
    silex::NumberField field = quadratic_x2_minus_2();
    silex::Order order = silex::test::equation_order(field);
    silex::FactorBase empty(order), base(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_ui(bound, 2);
    assert(base.build(sflint::FmpzConstRef(bound)));
    silex::PrimeIdeal ramified(order), inert(order), split(order);
    assert(first_prime_above(ramified, order, 2));
    assert(first_prime_above(inert, order, 11));
    assert(first_prime_above(split, order, 7));
    assert(ramified.ramification_index() == 2 && inert.residue_degree() == 2);
    silex::OrderElement alpha(order);
    silex::Ideal principal(order);
    silex::DiagnosticsContext diagnostics;
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    slong debug_failures = 0;
    diagnostics.debug_failure_user = &debug_failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    const auto check = [&](const silex::FactorBase& fb,
                           const silex::PrimeIdeal& required, bool expected) {
        assert(principal.set_principal(alpha));
        bool reference = !expected;
        assert(full_factorization_required_prime_result(reference, principal, fb, required));
        assert(reference == expected);
        for (Mode mode : {Mode::selected, Mode::full_factorization, Mode::order_element_direct}) {
            search::FactorBaseWitnessAudit audit;
            assert(search::factor_base_principal_witness(fb, required, alpha, mode,
                                                        &diagnostics, &audit) == expected);
            assert(debug_failures == 0);
            if (mode != Mode::full_factorization) assert(!audit.used_reference);
            if (mode == Mode::full_factorization) assert(audit.used_reference);
        }
    };
    for (slong scalar : {1, -1, 11, -11, 121, 13}) {
        assert(alpha.set_si(scalar));
        check(empty, inert, scalar == 11 || scalar == -11);
    }
    assert(alpha.set_si(7));
    check(empty, split, false); // The other prime above 7 is missing from the base.
    assert(alpha.set_si(2));
    check(empty, ramified, false); // v_P(2)=2, not one.
    sflint::FmpzMat coordinates(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 1), 1);
    assert(alpha.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
    check(empty, ramified, true); // (theta)=P in Z[sqrt(2)].
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0, 0), 1);
    assert(alpha.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
    check(empty, inert, false); // A non-scalar unit has no required-prime support.
    fmpz_set_ui(fmpz_mat_entry(coordinates.raw(), 0, 0), 11);
    fmpz_mul_2exp(fmpz_mat_entry(coordinates.raw(), 0, 0),
                  fmpz_mat_entry(coordinates.raw(), 0, 0), 80);
    fmpz_zero(fmpz_mat_entry(coordinates.raw(), 0, 1));
    assert(alpha.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
    check(base, inert, true);

    assert(alpha.set_si(0));
    assert(!search::factor_base_principal_witness(empty, inert, alpha));
    assert(alpha.set_si(11));
    search::FactorBaseWitnessAudit audit;
    silex::OrderElement undefined;
    assert(!search::factor_base_principal_witness(empty, inert, undefined));
    FieldSetup other = cubic_x3_minus_2();
    silex::OrderElement wrong_parent(other.maximal_order);
    assert(wrong_parent.set_si(11));
    assert(!search::factor_base_principal_witness(empty, inert, wrong_parent));
    assert(!search::factor_base_principal_witness(base, ramified, alpha));
    assert(!search::factor_base_principal_witness(empty, inert, alpha, static_cast<Mode>(99)));
    // Unsupported direct scope must use the reference, not grant a proof.
    // The public API cannot withdraw maximality, so the internal hook builds
    // that state.
    silex::detail::OrderAccess::set_maximality_unchecked(order, false);
    const bool fallback = search::factor_base_principal_witness(
            empty, inert, alpha, Mode::order_element_direct, nullptr, &audit);
    assert(audit.used_reference);
    assert(fallback == search::factor_base_principal_witness(
            empty, inert, alpha, Mode::full_factorization));
    return 0;
}

// Copies `prime`, a degree-one prime with residue polynomial x - r in the
// alpha convention, into `out` with the residue polynomial x - r - p: the
// same residue map onto F_p, so the same prime ideal, recorded with an
// unreduced coefficient.  PrimeIdeal::equal tells the two apart.
bool set_unreduced_residue_variant(silex::PrimeIdeal& out,
                                   const silex::PrimeIdeal& prime) noexcept {
    using silex::detail::PrimeIdealAccess;
    sflint::FmpzPoly residue;
    sflint::Fmpz p;
    sflint::Fmpz constant;
    if (prime.residue_degree() != 1 ||
        PrimeIdealAccess::residue_uses_integral_generator(prime) ||
        !prime.residue_polynomial(sflint::FmpzPolyRef(residue)) ||
        !prime.rational_prime(sflint::FmpzRef(p)) || !out.set(prime)) {
        return false;
    }
    fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
    fmpz_add(constant.raw(), constant.raw(), p.raw());
    fmpz_poly_set_coeff_fmpz(residue.raw(), 0, constant.raw());
    return silex::detail::TestHookAccess::set_alpha_residue_polynomial(
                   out, sflint::FmpzPolyConstRef(residue)) &&
           !out.equal(prime);
}

// The witness rule matches the required prime as an ideal.  In x^3 - 2 the
// element 1 + alpha^2 has norm 5, so (1 + alpha^2) is a degree-one prime P
// above 5 and witnesses P over the empty base.  A copy of P recorded with
// other residue data is the same prime and has the same witness, under the
// full-factorization predicate and under the debug reference oracle that
// checks the direct predicates (diagnostics preset).
int test_witness_matches_required_prime_as_ideal() {
    namespace search = silex::detail::relation_search;
    using Mode = search::FactorBaseWitnessPredicate;
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    assert(order.is_maximal());
    silex::FactorBase empty(order);

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    silex::Element ambient(setup.field);
    silex::OrderElement generator(order);
    silex::Ideal principal(order);
    assert(ambient.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    assert(generator.set_element(ambient));
    assert(principal.set_principal(generator));

    silex::IdealFactorization factorization(order);
    silex::PrimeIdeal prime(order);
    slong exponent = 0;
    assert(factorization.factor(principal));
    assert(factorization.length() == 1);
    assert(factorization.prime(prime, 0));
    assert(factorization.exponent(exponent, 0) && exponent == 1);
    sflint::Fmpz p;
    assert(prime.rational_prime(sflint::FmpzRef(p)));
    assert(sflint::fmpz_equal_si(p, 5) && prime.residue_degree() == 1);

    silex::PrimeIdeal variant(order);
    assert(set_unreduced_residue_variant(variant, prime));
    assert(!variant.equal(prime));

    silex::DiagnosticsContext diagnostics;
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    slong debug_failures = 0;
    diagnostics.debug_failure_user = &debug_failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    for (const silex::PrimeIdeal* required : {&prime, &variant}) {
        bool reference = false;
        assert(full_factorization_required_prime_result(
                reference, principal, empty, *required));
        assert(reference);
        // The direct predicates first: under debug checks each one is
        // compared with the reference oracle.
        for (Mode mode : {Mode::order_element_direct, Mode::selected,
                          Mode::full_factorization}) {
            search::FactorBaseWitnessAudit audit;
            assert(search::factor_base_principal_witness(
                    empty, *required, generator, mode, &diagnostics, &audit));
            assert(debug_failures == 0);
            assert(audit.used_reference == (mode == Mode::full_factorization));
        }
    }
    return 0;
}

void check_scalar_reference(const silex::FactorBase& base,
                            const silex::PrimeIdeal& prime,
                            sflint::FmpzConstRef scalar) {
    namespace search = silex::detail::relation_search;
    using Mode = search::FactorBaseWitnessPredicate;
    silex::Element ambient(*base.parent()->parent());
    silex::OrderElement element(*base.parent());
    silex::Ideal principal(*base.parent());
    assert(ambient.set_fmpz(scalar) && element.set_element(ambient));
    assert(principal.set_principal(element));
    bool expected = false;
    assert(full_factorization_required_prime_result(expected, principal, base, prime));
    silex::DiagnosticsContext diagnostics;
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    slong failures = 0;
    diagnostics.debug_failure_user = &failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    bool direct = !expected;
    assert(silex::detail::scalar_factor_over_base_with_required_prime(
            direct, scalar, base, prime, &diagnostics));
    assert(direct == expected && failures == 0);
    for (Mode mode : {Mode::selected, Mode::order_element_direct, Mode::full_factorization}) {
        search::FactorBaseWitnessAudit audit;
        assert(search::factor_base_scalar_witness(base, prime, scalar, mode,
                                                   &diagnostics, &audit) == expected);
        assert(failures == 0);
        assert(audit.used_scalar_direct == (mode == Mode::selected));
        if (mode == Mode::full_factorization) assert(audit.used_reference);
    }
}

int test_scalar_witness_reference_cases() {
    namespace search = silex::detail::relation_search;
    silex::NumberField field = quadratic_x2_minus_2();
    silex::Order order = silex::test::equation_order(field);
    silex::FactorBase empty(order), complete(order), reordered(order);
    sflint::Fmpz scalar, bound;
    sflint::fmpz_set_ui(bound, 7);
    assert(complete.build(sflint::FmpzConstRef(bound)));
    assert(reordered.build_lll_relation_base(sflint::FmpzConstRef(bound)));
    silex::PrimeIdeal inert(order), split(order), ramified(order);
    assert(first_prime_above(inert, order, 11));
    assert(first_prime_above(split, order, 7));
    assert(first_prime_above(ramified, order, 2));
    for (slong value : {-143, -121, -22, -11, -7, -2, -1, 1, 2, 7, 11, 22, 121, 143}) {
        sflint::fmpz_set_si(scalar, value);
        for (const auto* base : {&empty, &complete, &reordered})
            check_scalar_reference(*base, inert, sflint::FmpzConstRef(scalar));
        check_scalar_reference(empty, split, sflint::FmpzConstRef(scalar));
        check_scalar_reference(empty, ramified, sflint::FmpzConstRef(scalar));
    }
    fmpz_set_ui(scalar.raw(), 11);
    fmpz_mul_2exp(scalar.raw(), scalar.raw(), 160);
    check_scalar_reference(complete, inert, sflint::FmpzConstRef(scalar));
    fmpz_neg(scalar.raw(), scalar.raw());
    check_scalar_reference(reordered, inert, sflint::FmpzConstRef(scalar));
    sflint::fmpz_zero(scalar);
    bool invalid_matches = true;
    assert(!silex::detail::scalar_factor_over_base_with_required_prime(
            invalid_matches, sflint::FmpzConstRef(scalar), empty, inert));
    assert(!invalid_matches);
    assert(!search::factor_base_scalar_witness(empty, inert, sflint::FmpzConstRef(scalar)));
    sflint::fmpz_set_ui(scalar, 11);
    silex::FactorBase undefined;
    silex::PrimeIdeal invalid(order);
    FieldSetup other = cubic_x3_minus_2();
    silex::PrimeIdeal wrong(other.maximal_order);
    assert(first_prime_above(wrong, other.maximal_order, 11));
    assert(!silex::detail::scalar_factor_over_base_with_required_prime(
            invalid_matches, sflint::FmpzConstRef(scalar), empty, wrong));
    assert(!silex::detail::scalar_factor_over_base_with_required_prime(
            invalid_matches, sflint::FmpzConstRef(scalar), empty, invalid));
    assert(!silex::detail::scalar_factor_over_base_with_required_prime(
            invalid_matches, sflint::FmpzConstRef(scalar), undefined, inert));
    assert(!silex::detail::scalar_factor_over_base_with_required_prime(
            invalid_matches, sflint::FmpzConstRef(scalar), complete, split));
    assert(!search::factor_base_scalar_witness(undefined, inert, sflint::FmpzConstRef(scalar)));
    assert(!search::factor_base_scalar_witness(empty, invalid, sflint::FmpzConstRef(scalar)));
    assert(!search::factor_base_scalar_witness(empty, wrong, sflint::FmpzConstRef(scalar)));
    assert(!search::factor_base_scalar_witness(complete, split, sflint::FmpzConstRef(scalar)));
    // Withdrawn maximality (internal hook; unreachable publicly) must refuse.
    silex::detail::OrderAccess::set_maximality_unchecked(order, false);
    assert(!silex::detail::scalar_factor_over_base_with_required_prime(
            invalid_matches, sflint::FmpzConstRef(scalar), empty, inert));
    search::FactorBaseWitnessAudit audit;
    assert(!search::factor_base_scalar_witness(empty, inert, sflint::FmpzConstRef(scalar),
            search::FactorBaseWitnessPredicate::selected, nullptr, &audit));
    assert(audit.used_reference && !audit.used_scalar_direct);

    // Same maximal lattice in a basis whose first vector is not one.
    silex::Order permuted(field);
    sflint::FmpqMat basis(2, 2);
    fmpq_one(fmpq_mat_entry(basis.raw(), 0, 1));
    fmpq_one(fmpq_mat_entry(basis.raw(), 1, 0));
    assert(permuted.set_basis(sflint::FmpqMatConstRef(basis)));
    permuted = silex::test::verified_maximal_order(permuted);
    silex::FactorBase permuted_base(permuted);
    sflint::fmpz_set_ui(bound, 2);
    assert(permuted_base.build(sflint::FmpzConstRef(bound)));
    silex::PrimeIdeal permuted_inert(permuted);
    assert(first_prime_above(permuted_inert, permuted, 11));
    for (slong value : {-121, -22, 11, 22, 121, 143}) {
        sflint::fmpz_set_si(scalar, value);
        check_scalar_reference(permuted_base, permuted_inert, sflint::FmpzConstRef(scalar));
    }
    return 0;
}

int test_scalar_noncontiguous_blocks() {
    FieldSetup setup = truncated_decomposition_counterexample_fixture();
    silex::FactorBase base(setup.maximal_order);
    sflint::Fmpz bound, p, scalar;
    sflint::fmpz_set_ui(bound, 7);
    assert(base.build_lll_relation_base(sflint::FmpzConstRef(bound)));
    bool saw_noncontiguous = false;
    for (slong block = 0; block < base.rational_prime_block_count(); ++block) {
        slong start = 0, length = 0;
        if (!base.rational_prime_block(p, start, length, block))
            saw_noncontiguous = true;
    }
    assert(saw_noncontiguous);
    silex::PrimeIdeal required(setup.maximal_order);
    assert(first_prime_above(required, setup.maximal_order, 11));
    for (slong value : {-2310, -22, 22, 2310, 121}) {
        sflint::fmpz_set_si(scalar, value);
        check_scalar_reference(base, required, sflint::FmpzConstRef(scalar));
    }
    return 0;
}

int test_witness_partial_rational_prime_blocks() {
    namespace search = silex::detail::relation_search;
    using Mode = search::FactorBaseWitnessPredicate;
    for (bool multiple_omissions : {false, true}) {
        FieldSetup setup = multiple_omissions ? quintic_lower_interval_fixture()
                                              : truncated_decomposition_counterexample_fixture();
        const auto& order = setup.maximal_order;
        sflint::Fmpz p;
        sflint::fmpz_set_ui(p, multiple_omissions ? 3 : 2);
        silex::FactorBase partial(order);
        assert(partial.build_prime_ideal_norm_bounded(sflint::FmpzConstRef(p)));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
        const silex::PrimeIdeal* required = nullptr;
        slong omitted = 0;
        for (slong i = 0; i < primes.size(); ++i) {
            const auto* prime = primes.at(i);
            assert(prime != nullptr);
            if (!partial.contains(*prime)) {
                assert(prime->residue_degree() == 2 && prime->ramification_index() == 1);
                required = prime;
                ++omitted;
            }
        }
        assert(required != nullptr && omitted == (multiple_omissions ? 2 : 1));
        silex::OrderElement alpha(order);
        silex::Ideal principal(order);
        assert(alpha.set_si(multiple_omissions ? 3 : 2));
        assert(principal.set_principal(alpha));
        bool reference = false;
        assert(full_factorization_required_prime_result(reference, principal, partial, *required));
        assert(reference == !multiple_omissions);
        check_scalar_reference(partial, *required, sflint::FmpzConstRef(p));
        search::FactorBaseWitnessAudit audit;
        assert(search::factor_base_principal_witness(partial, *required, alpha,
                Mode::order_element_direct, nullptr, &audit) == reference);
        assert(!audit.used_reference);
    }
    return 0;
}

int test_ideal_transforms_preserve_principal_witness() {
    FieldSetup setup = cubic_x3_minus_2();
    const silex::Order& order = setup.maximal_order;
    const silex::NumberField& field = setup.field;

    silex::OrderElement six_generator(order);
    silex::Ideal six(order);
    silex::Ideal primitive(order);
    silex::Ideal reconstructed(order);
    silex::Element primitive_back_multiplier(field);
    if (!six_generator.is_defined() || !six.is_defined() ||
        !primitive.is_defined() || !reconstructed.is_defined() ||
        !primitive_back_multiplier.is_defined() ||
        !six_generator.set_si(6) || !six.set_principal(six_generator) ||
        !primitive.set(six) || !primitive_back_multiplier.one() ||
        !silex::detail::relation_search::
                 factor_base_honesty_primitive_part(
                         primitive, primitive_back_multiplier) ||
        !primitive.is_one() ||
        !silex::detail::multiply_integral_ideal_by_element(
                reconstructed, primitive, primitive_back_multiplier) ||
        !reconstructed.equal(six)) {
        return 1;
    }

    sflint::Fmpz five;
    sflint::fmpz_set_ui(sflint::FmpzRef(five), UWORD(5));
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, order,
                                sflint::FmpzConstRef(five))) {
        return 1;
    }
    const silex::PrimeIdeal* degree_one_prime = nullptr;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime != nullptr && prime->ramification_index() == 1 &&
            prime->residue_degree() == 1) {
            degree_one_prime = prime;
            break;
        }
    }
    if (degree_one_prime == nullptr) {
        return 1;
    }

    silex::Ideal factor(order);
    sflint::FmpzMat factor_hnf(order.degree(), order.degree());
    if (!factor.is_defined() || !degree_one_prime->get_ideal(factor) ||
        !factor.get_hnf(sflint::FmpzMatRef(factor_hnf)) ||
        !sflint::fmpz_is_one(sflint::fmpz_mat_entry(
                sflint::FmpzMatConstRef(factor_hnf), 0, 0)) ||
        sflint::fmpz_is_one(sflint::fmpz_mat_entry(
                sflint::FmpzMatConstRef(factor_hnf), order.degree() - 1,
                order.degree() - 1))) {
        return 1;
    }

    sflint::Fmpz large_scalar;
    sflint::fmpz_one(sflint::FmpzRef(large_scalar));
    ::fmpz_mul_2exp(large_scalar.raw(), large_scalar.raw(), 102);
    silex::Element large_element(field);
    silex::OrderElement large_generator(order);
    silex::Ideal large_scalar_ideal(order);
    silex::Ideal large(order);
    if (!large_element.is_defined() || !large_generator.is_defined() ||
        !large_scalar_ideal.is_defined() || !large.is_defined() ||
        !large_element.set_fmpz(sflint::FmpzConstRef(large_scalar)) ||
        !large_generator.set_element(large_element) ||
        !large_scalar_ideal.set_principal(large_generator) ||
        !large.multiply(large_scalar_ideal, factor)) {
        return 1;
    }

    silex::Ideal original(order);
    silex::Ideal reduced_reconstructed(order);
    silex::Element reduction_back_multiplier(field);
    if (!original.is_defined() || !reduced_reconstructed.is_defined() ||
        !reduction_back_multiplier.is_defined() || !original.set(large) ||
        !reduction_back_multiplier.one() ||
        !silex::detail::relation_search::
                 factor_base_honesty_reduce_large_ideal(
                         large, reduction_back_multiplier, 200, nullptr) ||
        large.equal(original) || reduction_back_multiplier.equal_si(1) ||
        !silex::detail::multiply_integral_ideal_by_element(
                reduced_reconstructed, large, reduction_back_multiplier) ||
        !reduced_reconstructed.equal(original)) {
        return 1;
    }
    return 0;
}

void norm_prefilter_reject_callback(void* user,
                                    silex::DiagnosticsModule,
                                    silex::ProfileEvent event,
                                    const char*,
                                    const char* label) noexcept {
    if (label != nullptr && event == silex::ProfileEvent::event &&
        std::strcmp(label,
                    "class_group.honesty_principal_witness.norm_prefilter_reject") ==
                0) {
        ++static_cast<RejectCofactorCounter*>(user)->count;
    }
}

int test_legacy_scan_predicate_equivalence() {
    namespace search = silex::detail::relation_search;
    using Mode = search::FactorBaseWitnessPredicate;
    // The direct predicates use the lattice-search norm prefilter; in a
    // debug-check build each of its rejections is checked against the direct
    // order-element test.
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    slong debug_failures = 0;
    diagnostics.debug_failure_user = &debug_failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    RejectCofactorCounter prefilter_rejects;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::class_group),
            norm_prefilter_reject_callback, &prefilter_rejects);
    for (slong fixture = 0; fixture < 4; ++fixture) {
        FieldSetup setup = fixture == 0 ? cubic_x3_minus_2()
                : fixture == 1 ? truncated_decomposition_counterexample_fixture()
                : fixture == 2 ? quintic_lower_interval_fixture()
                               : ramified_last_fixture();
        sflint::Fmpz active, required;
        sflint::fmpz_set_ui(active, fixture == 1 ? 1 : 7);
        sflint::fmpz_set_ui(required, fixture == 1 ? 2 : 9);
        silex::FactorBase base(setup.maximal_order);
        if (fixture != 1 && !base.build_relation_completion_base(
                sflint::FmpzConstRef(active))) return 1;
        bool expected_honest = false;
        search::FactorBaseHonestyScanAudit expected;
        const bool expected_ok = search::factor_base_honesty_check(
                expected_honest, base, sflint::FmpzConstRef(active),
                sflint::FmpzConstRef(required), nullptr, UWORD(42), false,
                128, search::FactorBaseWitnessEscalation::none, nullptr,
                &expected, Mode::full_factorization);
        for (Mode mode : {Mode::selected, Mode::order_element_direct}) {
            bool honest = !expected_honest;
            search::FactorBaseHonestyScanAudit audit;
            const bool ok = search::factor_base_honesty_check(
                    honest, base, sflint::FmpzConstRef(active),
                    sflint::FmpzConstRef(required), nullptr, UWORD(42), false,
                    128, search::FactorBaseWitnessEscalation::none,
                    &diagnostics, &audit, mode);
            if (ok != expected_ok || honest != expected_honest ||
                audit.rational_prime_checks != expected.rational_prime_checks ||
                audit.checks_at_or_below_active_bound != expected.checks_at_or_below_active_bound ||
                audit.witness_targets != expected.witness_targets ||
                audit.witnessed_targets != expected.witnessed_targets ||
                audit.final_random_state != expected.final_random_state) return 1;
        }
        if (fixture == 1 && (!expected_ok || expected_honest ||
                expected.witness_targets != 2 || expected.witnessed_targets != 1)) return 1;
    }
    if (debug_failures != 0) return 1;
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    if (prefilter_rejects.count == 0) return 1;
#endif
    return 0;
}

int test_truncated_decomposition_withholds_honesty() {
    FieldSetup setup = truncated_decomposition_counterexample_fixture();
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::fmpz_one(sflint::FmpzRef(active_bound));
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), UWORD(2));

    // This is the exact norm-at-most-1 factor base: it is defined but empty.
    // At p=2 the first degree-one prime is principal, while the second is the
    // nontrivial class.  The degree-two prime has norm 4 and is outside the
    // requested cutoff, so it cannot justify omitting that second target.
    silex::FactorBase base(setup.maximal_order);
    if (!base.is_defined() || base.length() != 0) {
        return 1;
    }

    bool honest = true;
    silex::detail::relation_search::FactorBaseHonestyScanAudit audit;
    const bool check_succeeded =
            silex::detail::relation_search::factor_base_honesty_check(
                    honest, base, sflint::FmpzConstRef(active_bound),
                    sflint::FmpzConstRef(required_bound), nullptr, UWORD(0),
                    false, 128,
                    silex::detail::relation_search::
                            FactorBaseWitnessEscalation::none,
                    nullptr, &audit);

    // The scan itself succeeds, but it must withhold the honesty proof after
    // testing both relevant degree-one primes.  Skipping the final relevant
    // entry incorrectly reports this empty base as honest.
    return !check_succeeded || honest || audit.rational_prime_checks != 1 ||
                           audit.checks_at_or_below_active_bound != 0
            ? 1
            : 0;
}

int test_generation_bound_checker_rejects_missing_truncated_target() {
    FieldSetup setup = truncated_decomposition_counterexample_fixture();
    sflint::Fmpz two;
    sflint::Fmpz canonical_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(two), UWORD(2));
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(canonical_bound), setup.maximal_order) ||
        sflint::fmpz_cmp_ui(sflint::FmpzConstRef(canonical_bound), 2) < 0) {
        return 1;
    }

    silex::ClassGroupContext context(setup.maximal_order);
    if (!context.is_defined() ||
        !context.build_factor_base(sflint::FmpzConstRef(canonical_bound))) {
        return 1;
    }

    silex::PrimeIdealList primes;
    const silex::FactorBase* observed_base = context.factor_base();
    if (!silex::decompose_prime(primes, setup.maximal_order,
                                sflint::FmpzConstRef(two)) ||
        primes.size() != 3 || observed_base == nullptr) {
        return 1;
    }

    std::vector<slong> degree_one_indices;
    slong degree_two_index = -1;
    for (slong i = 0; i < primes.size(); ++i) {
        if (prime_has_ef(primes, i, 1, 1)) {
            degree_one_indices.push_back(i);
        } else if (prime_has_ef(primes, i, 1, 2) &&
                   degree_two_index < 0) {
            degree_two_index = i;
        } else {
            return 1;
        }
    }
    if (degree_one_indices.size() != 2 || degree_two_index < 0) {
        return 1;
    }
    const silex::PrimeIdeal* first_prime =
            primes.at(degree_one_indices.front());
    const silex::PrimeIdeal* second_prime =
            primes.at(degree_one_indices.back());
    const silex::PrimeIdeal* excluded_prime = primes.at(degree_two_index);
    const slong first_index = observed_base->index(*first_prime);
    const slong second_index = observed_base->index(*second_prime);
    const slong excluded_index = observed_base->index(*excluded_prime);
    if (first_index < 0 || second_index < 0 || excluded_index < 0 ||
        first_index == second_index || first_index == excluded_index ||
        second_index == excluded_index) {
        return 1;
    }

    // Deliberately simulate a proof base that lacks the second norm-2 target
    // and the norm-excluded degree-two prime while retaining the valid build
    // receipt and rational-prime block metadata.  This reaches the public
    // generation-bound gate without adding a test-only API.  The mutations
    // are made through the observer of the underlying non-const context and
    // are restored before returning.
    auto* mutable_base = const_cast<silex::FactorBase*>(observed_base);
    auto* second_slot = const_cast<silex::PrimeIdeal*>(
            mutable_base->prime_at(second_index));
    auto* excluded_slot = const_cast<silex::PrimeIdeal*>(
            mutable_base->prime_at(excluded_index));
    silex::PrimeIdeal saved_second(setup.maximal_order);
    silex::PrimeIdeal saved_excluded(setup.maximal_order);
    if (second_slot == nullptr || excluded_slot == nullptr ||
        !saved_second.set(*second_slot) ||
        !saved_excluded.set(*excluded_slot)) {
        return 1;
    }
    const bool second_mutated = second_slot->set(*first_prime);
    const bool excluded_mutated =
            second_mutated && excluded_slot->set(*first_prime);
    if (!second_mutated || !excluded_mutated) {
        static_cast<void>(second_slot->set(saved_second));
        static_cast<void>(excluded_slot->set(saved_excluded));
        return 1;
    }

    const bool expected_membership =
            mutable_base->contains(*first_prime) &&
            !mutable_base->contains(*second_prime) &&
            !mutable_base->contains(*excluded_prime);
    const bool corrupted_check = context.check_factor_base_generation_bound(
            sflint::FmpzConstRef(two));

    const bool second_restored = second_slot->set(saved_second);
    const bool excluded_restored = excluded_slot->set(saved_excluded);
    const bool restored =
            second_restored && excluded_restored &&
            mutable_base->contains(*first_prime) &&
            mutable_base->contains(*second_prime) &&
            mutable_base->contains(*excluded_prime);
    const bool restored_check =
            restored && context.check_factor_base_generation_bound(
                                sflint::FmpzConstRef(two));
    return !expected_membership || corrupted_check || !restored_check ? 1 : 0;
}

int test_ramified_target_is_never_omitted() {
    FieldSetup setup = ramified_last_fixture();
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::Fmpz two;
    sflint::Fmpz three;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), UWORD(3));
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), UWORD(4));
    sflint::fmpz_set_ui(sflint::FmpzRef(two), UWORD(2));
    sflint::fmpz_set_ui(sflint::FmpzRef(three), UWORD(3));

    silex::FactorBase base(setup.maximal_order);
    if (!base.is_defined() ||
        !base.build_prime_ideal_norm_bounded(
                sflint::FmpzConstRef(active_bound))) {
        return 1;
    }

    silex::PrimeIdealList primes;
    for (sflint::FmpzConstRef p :
         {sflint::FmpzConstRef(two), sflint::FmpzConstRef(three)}) {
        if (!silex::decompose_prime(primes, setup.maximal_order, p) ||
            primes.size() != 2) {
            return 1;
        }
        slong retained_unramified = 0;
        slong omitted_ramified = 0;
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            if (prime == nullptr) {
                return 1;
            }
            if (prime_has_ef(primes, i, 1, 1) && base.contains(*prime)) {
                ++retained_unramified;
            } else if (prime_has_ef(primes, i, 2, 2) &&
                       !base.contains(*prime)) {
                ++omitted_ramified;
            } else {
                return 1;
            }
        }
        if (retained_unramified != 1 || omitted_ramified != 1) {
            return 1;
        }
    }

    // The only omitted target up to the required bound is the ramified
    // norm-4 prime above 2.  Both witness searches must target it,
    // independently of the order in which the decomposition producer
    // returned the two primes, and witness it.
    for (bool direct : {false, true}) {
        bool honest = false;
        silex::detail::relation_search::FactorBaseHonestyScanAudit audit;
        const bool check_succeeded =
                silex::detail::relation_search::factor_base_honesty_check(
                        honest, base, sflint::FmpzConstRef(active_bound),
                        sflint::FmpzConstRef(required_bound), nullptr,
                        UWORD(0), direct, 128,
                        silex::detail::relation_search::
                                FactorBaseWitnessEscalation::none,
                        nullptr, &audit);
        if (!check_succeeded || !honest || audit.witness_targets != 1 ||
            audit.witnessed_targets != 1 || audit.unwitnessed.recorded ||
            audit.rational_prime_checks != 2 ||
            audit.checks_at_or_below_active_bound != 2) {
            return 1;
        }
    }
    return 0;
}

int test_scan_detects_lower_interval_omissions() {
    FieldSetup setup = quintic_lower_interval_fixture();
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::Fmpz two;
    sflint::Fmpz three;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), UWORD(7));
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), UWORD(9));
    sflint::fmpz_set_ui(sflint::FmpzRef(two), UWORD(2));
    sflint::fmpz_set_ui(sflint::FmpzRef(three), UWORD(3));

    silex::FactorBase base(setup.maximal_order);
    if (!base.is_defined() ||
        !base.build_relation_completion_base(
                sflint::FmpzConstRef(active_bound))) {
        return 1;
    }

    // The polynomial is Eisenstein at 2, so the sole prime above 2 is
    // retained. Modulo 3 it has residue degrees 1, 2, 2: the norm-3 ideal is
    // retained while both norm-9 ideals lie in the lower omitted interval.
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, setup.maximal_order,
                                sflint::FmpzConstRef(two)) ||
        primes.size() != 1 || primes.at(0) == nullptr ||
        !base.contains(*primes.at(0)) ||
        !silex::decompose_prime(primes, setup.maximal_order,
                                sflint::FmpzConstRef(three))) {
        return 1;
    }

    slong retained_degree_one = 0;
    slong omitted_degree_two = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr) {
            return 1;
        }
        if (prime->residue_degree() == 1 && base.contains(*prime)) {
            ++retained_degree_one;
        } else if (prime->residue_degree() == 2 &&
                   !base.contains(*prime)) {
            ++omitted_degree_two;
        }
    }
    if (retained_degree_one != 1 || omitted_degree_two != 2) {
        return 1;
    }

    // The norm-9 ideals above 3 hold one proof target (the last unramified
    // factor of (3) is discharged by the others).  Both witness searches
    // must target it, so the scan reaches p=3 and goes on to the required
    // bound, and witness it.
    for (bool direct : {false, true}) {
        bool honest = false;
        silex::detail::relation_search::FactorBaseHonestyScanAudit audit;
        const bool check_succeeded =
                silex::detail::relation_search::factor_base_honesty_check(
                        honest, base, sflint::FmpzConstRef(active_bound),
                        sflint::FmpzConstRef(required_bound), nullptr,
                        UWORD(0), direct, 128,
                        silex::detail::relation_search::
                                FactorBaseWitnessEscalation::none,
                        nullptr, &audit);
        if (!check_succeeded || !honest || audit.witness_targets != 1 ||
            audit.witnessed_targets != 1 || audit.unwitnessed.recorded ||
            audit.rational_prime_checks != 4 ||
            audit.checks_at_or_below_active_bound != 4) {
            return 1;
        }
    }
    return 0;
}

FieldSetup imaginary_quadratic_fixture(slong d) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, d);
    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

struct HonestyLogCapture {
    slong exhausted_messages = 0;
    std::string last_detail;
};

void capture_honesty_log(void* user,
                         silex::DiagnosticsModule,
                         silex::LogLevel,
                         const char*,
                         const char* message,
                         const char* detail) noexcept {
    auto* capture = static_cast<HonestyLogCapture*>(user);
    if (message != nullptr &&
        std::strcmp(message,
                    "be_honest witness search exhausted at final caps") ==
                0) {
        ++capture->exhausted_messages;
        capture->last_detail = detail == nullptr ? "" : detail;
    }
}

namespace search = silex::detail::relation_search;

// Expected detail-log text for an unwitnessed prime on the lattice route.
constexpr const char* kStageZeroCaps =
        "stage=0/0 search=lattice radius=8 twists=16 random_tries=50";
constexpr const char* kFinalStageCaps =
        "stage=3/3 search=lattice radius=64 twists=128 random_tries=400";

void enable_honesty_log(silex::DiagnosticsContext& diagnostics,
                        HonestyLogCapture& capture) noexcept {
    silex::diagnostics_context_init(diagnostics);
    silex::diagnostics_set_logging(
            diagnostics, silex::LogLevel::detail,
            silex::diagnostics_module_bit(
                    silex::DiagnosticsModule::class_group),
            capture_honesty_log, &capture);
}

bool logged_caps(const HonestyLogCapture& capture,
                 const char* caps) noexcept {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    return capture.exhausted_messages == 1 &&
           capture.last_detail.find(caps) != std::string::npos;
#else
    (void) capture;
    (void) caps;
    return true;
#endif
}

// An active bound of 1 keeps the defined but empty norm-at-most-1 base.
bool run_legacy_honesty_scan(bool& honest,
                             search::FactorBaseHonestyScanAudit& audit,
                             HonestyLogCapture& capture,
                             slong d,
                             ulong active,
                             ulong required,
                             search::FactorBaseWitnessEscalation escalation)
        noexcept {
    FieldSetup setup = imaginary_quadratic_fixture(d);
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), active);
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), required);
    silex::FactorBase base(setup.maximal_order);
    if (!base.is_defined() ||
        (active > 1 && !base.build_relation_completion_base(
                               sflint::FmpzConstRef(active_bound)))) {
        return false;
    }
    silex::DiagnosticsContext diagnostics;
    enable_honesty_log(diagnostics, capture);
    return search::factor_base_honesty_check(
            honest, base, sflint::FmpzConstRef(active_bound),
            sflint::FmpzConstRef(required_bound), nullptr, UWORD(0), false,
            128, escalation, &diagnostics, &audit);
}

// Q(sqrt(-1559)) has class number 51.  With the norm <= 3 relation base,
// one required prime below 300 has no witness at the stage-0 caps (radius 8,
// 16 twists, 50 random tries).  Without escalation the scan fails closed
// there, as the reference be_honest does before enlarging the factor base;
// with bounded escalation the prime is witnessed at a later stage and the
// scan proves honesty.
int test_escalated_witness_search_proves_large_class_group_base() {
    bool honest = true;
    search::FactorBaseHonestyScanAudit audit;
    HonestyLogCapture capture;
    if (!run_legacy_honesty_scan(honest, audit, capture, 1559, 3, 300,
                                 search::FactorBaseWitnessEscalation::none) ||
        honest || audit.max_search_stage != 0 ||
        audit.escalated_witnessed_targets != 0 ||
        !logged_caps(capture, kStageZeroCaps) ||
        !audit.unwitnessed.recorded || audit.unwitnessed.stage != 0 ||
        audit.unwitnessed.max_stage != 0 ||
        audit.unwitnessed.radius != 8 ||
        audit.unwitnessed.max_twists != 16 ||
        audit.unwitnessed.random_tries != 50) {
        return 1;
    }

    honest = false;
    capture = HonestyLogCapture{};
    if (!run_legacy_honesty_scan(
                honest, audit, capture, 1559, 3, 300,
                search::FactorBaseWitnessEscalation::bounded)) {
        return 1;
    }
    return honest && audit.witness_targets > 0 &&
                   audit.witnessed_targets == audit.witness_targets &&
                   audit.escalated_witnessed_targets >= 1 &&
                   audit.max_search_stage >= 1 &&
                   audit.max_search_stage <=
                           search::kFactorBaseHonestyEscalationStages &&
                   capture.exhausted_messages == 0 &&
                   !audit.unwitnessed.recorded
            ? 0
            : 1;
}

// The final-stage lattice caps recorded for an unwitnessed degree-one prime
// above p.
bool final_stage_lattice_failure(
        const silex::detail::FactorBaseHonestyFailure& failure,
        slong p) noexcept {
    return failure.recorded && sflint::fmpz_equal_si(failure.p, p) &&
           failure.residue_degree == 1 && !failure.direct_witness_search &&
           failure.stage == search::kFactorBaseHonestyEscalationStages &&
           failure.max_stage == search::kFactorBaseHonestyEscalationStages &&
           failure.radius == 64 && failure.max_twists == 128 &&
           failure.random_tries == 400;
}

// Q(sqrt(-5)) has class number 2 and the prime above 2 is the nontrivial
// class.  The empty norm-at-most-1 base generates only the trivial class, so
// no search effort can witness that prime: the test checks that bounded
// escalation runs every stage, fails closed, and logs the final stage's caps,
// not that a cap happened to be too small.
int test_bounded_escalation_reports_final_caps() {
    bool honest = true;
    search::FactorBaseHonestyScanAudit audit;
    HonestyLogCapture capture;
    return run_legacy_honesty_scan(
                   honest, audit, capture, 5, 1, 2,
                   search::FactorBaseWitnessEscalation::bounded) &&
                   !honest && audit.witness_targets == 1 &&
                   audit.witnessed_targets == 0 &&
                   audit.max_search_stage ==
                           search::kFactorBaseHonestyEscalationStages &&
                   logged_caps(capture, kFinalStageCaps) &&
                   final_stage_lattice_failure(audit.unwitnessed, 2)
            ? 0
            : 1;
}

// The final-stage T2 caps recorded for an unwitnessed degree-one prime
// above p: the stage-0 caps (50 random tries, 500 factored candidates and
// 10^6 enumeration steps per ideal) doubled at each of the three stages.
bool final_stage_t2_failure(
        const silex::detail::FactorBaseHonestyFailure& failure,
        slong p) noexcept {
    return failure.recorded && sflint::fmpz_equal_si(failure.p, p) &&
           failure.residue_degree == 1 && failure.direct_witness_search &&
           failure.stage == search::kFactorBaseHonestyEscalationStages &&
           failure.max_stage == search::kFactorBaseHonestyEscalationStages &&
           failure.random_tries == 400 && failure.factor_attempts == 4000 &&
           failure.element_steps == 8000000;
}

namespace driver = silex::detail::relation_search;

struct HonestyCheckpointRun {
    driver::AnalyticFinishDecision decision =
            driver::AnalyticFinishDecision::incomplete;
    silex::detail::ClassUnitTransactionReport audit;
};

// Runs the relation-search honesty checkpoint for Q(sqrt(-21)) on the
// norm-at-most-2 base, optionally inside a class/unit transaction.  The class
// group is Z/2 x Z/2; the base holds only the ramified prime above 2, and the
// ramified prime above 3 lies in another nontrivial class (neither 3 nor 6 is
// a norm x^2 + 21 y^2), so no search effort can witness it and the scan always
// stops there.  The sub-factor base for the random twists is that prime.
bool run_honesty_checkpoint(HonestyCheckpointRun& run,
                            bool factor_base_restart_available,
                            bool attach_transaction,
                            const silex::DiagnosticsContext* diagnostics =
                                    nullptr) noexcept {
    FieldSetup setup = imaginary_quadratic_fixture(21);
    sflint::Fmpz active_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), 2);
    silex::ClassGroupContext context(setup.maximal_order);
    if (!context.is_defined() ||
        !context.build_factor_base(sflint::FmpzConstRef(active_bound)) ||
        context.factor_base() == nullptr) {
        return false;
    }
    // The checkpoint needs a published presentation: (2) = P2^2.
    silex::Relation relation(*context.factor_base());
    silex::Element two(setup.field);
    if (!two.set_si(2) || !relation.set_generator(two) ||
        !context.append_relation(relation) ||
        !context.publish_presentation()) {
        return false;
    }
    run.audit.reset();
    silex::detail::ClassUnitTransactionContext transaction{
            run.audit, nullptr, {}, {}, false};
    if (attach_transaction) {
        silex::detail::ClassUnitTransactionAccess::set_run_context(
                context, &transaction);
    }
    silex::detail::RelationSearchControlState route_state;
    silex::detail::SubfactorBaseSchedule schedule;
    schedule.subfactor_base.push_back(0);
    schedule.excluded_indices.assign(
            static_cast<std::size_t>(context.factor_base()->length()), 0);
    run.decision = driver::apply_honesty_check(
            context, route_state, true, setup.maximal_order, &schedule,
            sflint::FmpzConstRef(active_bound),
            factor_base_restart_available, diagnostics);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return true;
}

// When the honesty checkpoint ends the route (no factor-base restart left),
// the class/unit transaction report names the unwitnessed prime and the
// final stage's caps, and a failed transaction reports
// factor_base_honesty_unwitnessed.  With a restart available the checkpoint
// asks for it and records nothing; clearing the record at the next attempt
// restores class_unit_computation_failed.
//
// The checkpoint searches witnesses with the capped T2 enumeration of
// buch2.c:be_honest: each ideal stops after maxtry_FACT = 500 factored
// candidates (Fincke_Pohst_ideal), and the prime is retried with at most
// maxtry_HONEST = 50 random sub-factor-base twists.  The prime above 3 has
// no witness, so every ideal searched stops at the candidate cap: the 51 of
// stage 0 before the restart request, and the 401 of the final stage, at
// the doubled cap, before the terminal failure.  The uncapped
// scalar/lattice search would log search=lattice instead.
int test_terminal_honesty_failure_names_transaction_reason() {
    HonestyLogCapture capture;
    silex::DiagnosticsContext diagnostics;
    enable_honesty_log(diagnostics, capture);
    HonestyCheckpointRun terminal;
    if (!run_honesty_checkpoint(terminal, false, true, &diagnostics) ||
        terminal.decision != driver::AnalyticFinishDecision::failed ||
        !logged_caps(capture,
                     "stage=3/3 search=t2 random_tries=400 "
                     "factor_attempts=4000 element_steps=8000000 "
                     "ideals=401 factor_attempt_cap_hits=401 "
                     "element_step_cap_hits=0")) {
        return 1;
    }
    const silex::detail::FactorBaseHonestyFailure& failure =
            terminal.audit.factor_base_honesty_failure;
    if (!final_stage_t2_failure(failure, 3) ||
        std::strcmp(silex::detail::class_unit_computation_failure_reason(
                            terminal.audit),
                    "factor_base_honesty_unwitnessed") != 0) {
        return 1;
    }

    capture = HonestyLogCapture{};
    HonestyCheckpointRun restart;
    if (!run_honesty_checkpoint(restart, true, true, &diagnostics) ||
        restart.decision !=
                driver::AnalyticFinishDecision::restart_factor_base ||
        !logged_caps(capture,
                     "stage=0/0 search=t2 random_tries=50 "
                     "factor_attempts=500 element_steps=1000000 ideals=51 "
                     "factor_attempt_cap_hits=51 element_step_cap_hits=0") ||
        restart.audit.factor_base_honesty_failure.recorded ||
        std::strcmp(silex::detail::class_unit_computation_failure_reason(
                            restart.audit),
                    "class_unit_computation_failed") != 0) {
        return 1;
    }

    HonestyCheckpointRun detached;
    if (!run_honesty_checkpoint(detached, false, false) ||
        detached.decision != driver::AnalyticFinishDecision::failed ||
        detached.audit.factor_base_honesty_failure.recorded) {
        return 1;
    }

    FieldSetup setup = imaginary_quadratic_fixture(21);
    silex::ClassGroupContext context(setup.maximal_order);
    silex::detail::ClassUnitTransactionContext transaction{
            terminal.audit, nullptr, {}, {}, false};
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               &transaction);
    silex::detail::ClassGroupRelationSearchAccess::
            clear_factor_base_honesty_failure(context);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return !terminal.audit.factor_base_honesty_failure.recorded &&
                   std::strcmp(
                           silex::detail::class_unit_computation_failure_reason(
                                   terminal.audit),
                           "class_unit_computation_failed") == 0
            ? 0
            : 1;
}

// A later class-group attempt in the same class/unit transaction clears the
// record an earlier attempt's terminal honesty checkpoint left, so the
// transaction's reason belongs to its last attempt.  The clear runs before
// the attempt's guards: an extension of a context with no factor base and a
// candidate computation with an invalid search radius are both rejected
// before any search, and both still clear.
int test_later_attempt_clears_honesty_failure() {
    HonestyCheckpointRun earlier;
    if (!run_honesty_checkpoint(earlier, false, true) ||
        !final_stage_t2_failure(
                earlier.audit.factor_base_honesty_failure, 3)) {
        return 1;
    }
    using Access = silex::detail::ClassGroupRelationSearchAccess;

    FieldSetup setup = imaginary_quadratic_fixture(21);
    silex::ClassGroupContext context(setup.maximal_order);
    silex::detail::ClassUnitTransactionContext transaction{
            earlier.audit, nullptr, {}, {}, false};
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               &transaction);
    sflint::Fmpz bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), 20);
    silex::detail::ClassGroupRelationOptions options =
            silex::detail::class_group_relation_options(
                    setup.maximal_order, silex::ClassGroupComputeOptions{});

    const bool extended = Access::extend_relation_kernel_units(
            context, setup.maximal_order, sflint::FmpzConstRef(bound),
            options);
    const bool extension_cleared =
            !extended && !earlier.audit.factor_base_honesty_failure.recorded;

    silex::detail::FactorBaseHonestyFailure present;
    present.recorded = true;
    sflint::fmpz_set_si(sflint::FmpzRef(present.p), 3);
    Access::record_factor_base_honesty_failure(context, std::move(present));
    silex::detail::ClassGroupRelationOptions invalid_options = options;
    invalid_options.coordinate_search_radius = 0;
    const bool rejected = Access::compute_relation_candidate(
            context, setup.maximal_order, sflint::FmpzConstRef(bound),
            invalid_options);
    const bool rejection_cleared =
            !rejected && !earlier.audit.factor_base_honesty_failure.recorded;

    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return extension_cleared && rejection_cleared &&
                   std::strcmp(
                           silex::detail::class_unit_computation_failure_reason(
                                   earlier.audit),
                           "class_unit_computation_failed") == 0
            ? 0
            : 1;
}

// A Q(sqrt(-21)) context on the norm-at-most-`active` base (2 or 3) with
// (2) = P2^2 and, for 3, (3) = P3^2 published.  P2, P3 and P2 P3 are all
// nonprincipal (none of 2, 3, 6 is a norm x^2 + 21 y^2), so the {P2, P3}
// base generates the class group Z/2 x Z/2 and its honesty scan passes,
// while the {P2} base cannot witness P3 and its scan fails.
bool d21_checkpoint_context(silex::ClassGroupContext& context,
                            const FieldSetup& setup,
                            ulong active) noexcept {
    sflint::Fmpz active_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), active);
    if (!context.is_defined() ||
        !context.build_factor_base(sflint::FmpzConstRef(active_bound)) ||
        context.factor_base() == nullptr) {
        return false;
    }
    for (slong p = 2; p <= static_cast<slong>(active); ++p) {
        silex::Relation relation(*context.factor_base());
        silex::Element generator(setup.field);
        if (!generator.set_si(p) || !relation.set_generator(generator) ||
            !context.append_relation(relation)) {
            return false;
        }
    }
    return context.publish_presentation();
}

driver::AnalyticFinishDecision d21_checkpoint(
        silex::ClassGroupContext& context,
        const FieldSetup& setup,
        ulong active,
        silex::detail::ClassUnitTransactionContext& transaction) noexcept {
    sflint::Fmpz active_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), active);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               &transaction);
    silex::detail::RelationSearchControlState route_state;
    const driver::AnalyticFinishDecision decision = driver::apply_honesty_check(
            context, route_state, true, setup.maximal_order, nullptr,
            sflint::FmpzConstRef(active_bound), false, nullptr);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return decision;
}

// Once a route has passed the honesty checkpoint
// (RelationSearchControlState::honesty_checked), a later call finishes only
// while factor-base generation stays checked up to the recorded
// requirement.  Q(sqrt(-26)) has requirement floor(sqrt(104 / 3)) = 5.  The
// public check accepts any bound up to the build bound, so a check at 2
// marks the generation check verified without showing that the base
// generates the class group: the checkpoint must then fail.  The same
// candidate checked at 5 finishes.
int test_passed_honesty_checkpoint_requires_required_bound() {
    FieldSetup setup = imaginary_quadratic_fixture(26);
    if (!setup.maximal_order.is_maximal()) {
        return 1;
    }
    sflint::Fmpz bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), 5);
    for (const ulong checked : {UWORD(2), UWORD(5)}) {
        silex::ClassGroupContext context;
        sflint::Fmpz active_bound;
        sflint::Fmpz required;
        sflint::Fmpz checked_bound;
        sflint::Fmpz recorded_check;
        sflint::fmpz_set_ui(sflint::FmpzRef(checked_bound), checked);
        if (!context.compute_candidate(setup.maximal_order,
                                       sflint::FmpzConstRef(bound)) ||
            !context.factor_base_build_bound(
                    sflint::FmpzRef(active_bound)) ||
            !context.factor_base_generation_bound(sflint::FmpzRef(required)) ||
            sflint::fmpz_cmp_ui(sflint::FmpzConstRef(required), 5) != 0 ||
            !context.check_factor_base_generation_bound(
                    sflint::FmpzConstRef(checked_bound)) ||
            !context.factor_base_generation_checked_bound(
                    sflint::FmpzRef(recorded_check)) ||
            sflint::fmpz_cmp_ui(sflint::FmpzConstRef(recorded_check),
                                checked) != 0) {
            return 1;
        }
        silex::detail::RelationSearchControlState route_state;
        route_state.honesty_checked = true;
        const driver::AnalyticFinishDecision decision =
                driver::apply_honesty_check(
                        context, route_state, true, setup.maximal_order,
                        nullptr, sflint::FmpzConstRef(active_bound), false,
                        nullptr);
        const driver::AnalyticFinishDecision expected = checked >= 5
                ? driver::AnalyticFinishDecision::finished
                : driver::AnalyticFinishDecision::failed;
        if (decision != expected) {
            return 1;
        }
    }
    return 0;
}

// Relation saturation proofs do not use factor-base generation: after the
// public check at 2 on Q(sqrt(-26)) (requirement 5, h = 6) the ell = 3 and
// ell = 2 proofs still verify, but the class group is not labelled proven,
// because no route publishes `proven` from the saturation status alone.
int test_saturation_proof_without_bounded_generation_is_not_proven() {
    FieldSetup setup = imaginary_quadratic_fixture(26);
    if (!setup.maximal_order.is_maximal()) {
        return 1;
    }
    silex::OrderUnitGroup units;
    if (!units.compute(setup.maximal_order)) {
        return 1;
    }
    sflint::Fmpz bound;
    sflint::Fmpz checked_bound;
    sflint::Fmpz aux_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), 5);
    sflint::fmpz_set_ui(sflint::FmpzRef(checked_bound), 2);
    sflint::fmpz_set_ui(sflint::FmpzRef(aux_bound), 200);
    silex::ClassGroupContext context;
    if (!context.compute_candidate(setup.maximal_order,
                                   sflint::FmpzConstRef(bound)) ||
        !context.check_factor_base_generation_bound(
                sflint::FmpzConstRef(checked_bound)) ||
        !context.try_certify_with_units(
                units, silex::CertificationMode::unknown, 80)) {
        return 1;
    }
    for (const ulong ell_value : {UWORD(3), UWORD(2)}) {
        sflint::Fmpz ell;
        sflint::fmpz_set_ui(sflint::FmpzRef(ell), ell_value);
        if (!context.try_prove_relation_saturation_with_units(
                    units, sflint::FmpzConstRef(ell),
                    sflint::FmpzConstRef(aux_bound))) {
            return 1;
        }
    }
    if (context.relation_saturation_status() != silex::ProofState::verified ||
        context.relation_saturation_record_count() != 2) {
        return 1;
    }
    // Pinned: the saturation proofs alone do not publish `proven`.
    if (context.certification_status() == silex::CertificationMode::proven) {
        return 1;
    }
    return 0;
}

// A passed honesty scan is kept in the class/unit transaction with its
// factor base and bounds, and a later context of the same transaction takes
// the verdict only for an equal base at the same active and required bounds.
// A failed scan keeps nothing.  The reuse probe records a verdict for the
// {P2} base, whose own scan fails: the checkpoint finishes only when the
// record matches exactly, which shows that the scan was skipped then and run
// otherwise.
int test_transaction_reuses_passed_honesty_scan() {
    using Access = silex::detail::ClassGroupRelationSearchAccess;
    FieldSetup setup = imaginary_quadratic_fixture(21);
    silex::detail::ClassUnitTransactionReport audit;
    audit.reset();
    silex::detail::ClassUnitTransactionContext transaction{
            audit, nullptr, {}, {}, false};

    silex::ClassGroupContext failing(setup.maximal_order);
    if (!d21_checkpoint_context(failing, setup, 2) ||
        d21_checkpoint(failing, setup, 2, transaction) !=
                driver::AnalyticFinishDecision::failed ||
        transaction.factor_base_honesty.base.is_defined()) {
        return 1;
    }

    silex::ClassGroupContext passing(setup.maximal_order);
    sflint::Fmpz required;
    if (!d21_checkpoint_context(passing, setup, 3) ||
        !Access::recorded_generation_requirement(
                sflint::FmpzRef(required), passing, setup.maximal_order) ||
        sflint::fmpz_cmp_ui(sflint::FmpzConstRef(required), 3) <= 0 ||
        d21_checkpoint(passing, setup, 3, transaction) !=
                driver::AnalyticFinishDecision::finished ||
        !transaction.factor_base_honesty.base.equal(*passing.factor_base()) ||
        sflint::fmpz_cmp_ui(sflint::FmpzConstRef(
                                    transaction.factor_base_honesty
                                            .active_bound),
                            3) != 0 ||
        !sflint::fmpz_equal(
                sflint::FmpzConstRef(
                        transaction.factor_base_honesty.required_bound),
                sflint::FmpzConstRef(required))) {
        return 1;
    }

    // A later context with the {P2} base: the kept {P2, P3} verdict is for
    // another base, so the scan runs and fails.
    silex::ClassGroupContext other_base(setup.maximal_order);
    if (!d21_checkpoint_context(other_base, setup, 2) ||
        d21_checkpoint(other_base, setup, 2, transaction) !=
                driver::AnalyticFinishDecision::failed) {
        return 1;
    }

    sflint::Fmpz two;
    sflint::Fmpz required_plus_one;
    sflint::fmpz_set_ui(sflint::FmpzRef(two), 2);
    sflint::fmpz_add_ui(sflint::FmpzRef(required_plus_one),
                        sflint::FmpzConstRef(required), 1);
    struct Probe {
        sflint::FmpzConstRef active;
        sflint::FmpzConstRef required;
        driver::AnalyticFinishDecision expected;
    };
    const Probe probes[] = {
            {sflint::FmpzConstRef(two), sflint::FmpzConstRef(required),
             driver::AnalyticFinishDecision::finished},
            {sflint::FmpzConstRef(two),
             sflint::FmpzConstRef(required_plus_one),
             driver::AnalyticFinishDecision::failed},
            {sflint::FmpzConstRef(required), sflint::FmpzConstRef(required),
             driver::AnalyticFinishDecision::failed},
    };
    // Same bounds, different base: keep the {P2, P3} base at active bound 2
    // and the required bound, then a {P2} context at those same bounds must
    // still scan (and fail), not reuse.
    {
        silex::ClassGroupContext keeper(setup.maximal_order);
        if (!d21_checkpoint_context(keeper, setup, 3)) {
            return 1;
        }
        silex::detail::ClassUnitTransactionAccess::set_run_context(
                keeper, &transaction);
        Access::record_transaction_factor_base_honest(
                keeper, *keeper.factor_base(), sflint::FmpzConstRef(two),
                sflint::FmpzConstRef(required));
        silex::detail::ClassUnitTransactionAccess::set_run_context(keeper,
                                                                   nullptr);
        silex::ClassGroupContext narrow(setup.maximal_order);
        if (!d21_checkpoint_context(narrow, setup, 2) ||
            d21_checkpoint(narrow, setup, 2, transaction) !=
                    driver::AnalyticFinishDecision::failed) {
            return 1;
        }
    }
    sflint::Fmpz passing_checked;
    if (!passing.factor_base_generation_checked_bound(
                sflint::FmpzRef(passing_checked))) {
        return 1;
    }
    for (const Probe& probe : probes) {
        silex::ClassGroupContext probed(setup.maximal_order);
        if (!d21_checkpoint_context(probed, setup, 2)) {
            return 1;
        }
        silex::detail::ClassUnitTransactionAccess::set_run_context(
                probed, &transaction);
        Access::record_transaction_factor_base_honest(
                probed, *probed.factor_base(), probe.active, probe.required);
        silex::detail::ClassUnitTransactionAccess::set_run_context(probed,
                                                                   nullptr);
        const driver::AnalyticFinishDecision decision =
                d21_checkpoint(probed, setup, 2, transaction);
        if (decision != probe.expected ||
            (decision == driver::AnalyticFinishDecision::finished &&
             probed.factor_base_generation_checked_status() !=
                     silex::ProofState::verified)) {
            return 1;
        }
        if (decision == driver::AnalyticFinishDecision::finished) {
            // The reused receipt matches a fresh scan's.
            sflint::Fmpz checked;
            if (!probed.factor_base_generation_checked_bound(
                        sflint::FmpzRef(checked)) ||
                !sflint::fmpz_equal(sflint::FmpzConstRef(checked),
                                    sflint::FmpzConstRef(passing_checked)) ||
                probed.factor_base_generation_status() !=
                        passing.factor_base_generation_status()) {
                return 1;
            }
        }
    }

    // Without a transaction nothing is kept or reused.
    silex::ClassGroupContext detached(setup.maximal_order);
    return d21_checkpoint_context(detached, setup, 2) &&
                   !Access::transaction_factor_base_honest(
                           detached, *detached.factor_base(),
                           sflint::FmpzConstRef(two),
                           sflint::FmpzConstRef(required))
            ? 0
            : 1;
}

// The instance JSON's failure_detail value for a class/unit transaction
// report.
std::string failure_detail_json(
        const silex::detail::ClassUnitTransactionReport& report) {
    std::ostringstream out;
    silex_tools::write_class_unit_failure_detail_json(out, report);
    return out.str();
}

// failure_detail names the unwitnessed prime and the final stage's caps when
// the reason is factor_base_honesty_unwitnessed: the lattice caps the real
// checkpoint recorded, or the T2 caps of a T2 record.  The prime is a
// decimal string, small or beyond a machine word, never a number or null.  It
// is null for any other reason, even with a record present, and with no
// reason.
int test_failure_detail_json() {
    HonestyCheckpointRun terminal;
    if (!run_honesty_checkpoint(terminal, false, true)) {
        return 1;
    }
    silex::detail::ClassUnitTransactionReport& report = terminal.audit;
    report.failure_stage = silex::detail::ClassUnitStage::total;
    report.failure_reason =
            silex::detail::class_unit_computation_failure_reason(report);
    const bool checkpoint =
            failure_detail_json(report) ==
            "{\n"
            "    \"p\": \"3\",\n"
            "    \"residue_degree\": 1,\n"
            "    \"search\": \"t2\",\n"
            "    \"stage\": 3,\n"
            "    \"max_stage\": 3,\n"
            "    \"random_tries\": 400,\n"
            "    \"factor_attempts\": 4000,\n"
            "    \"element_steps\": 8000000\n"
            "  }";

    // The lattice route's caps, as a direct scan caller may still record.
    silex::detail::FactorBaseHonestyFailure& failure =
            report.factor_base_honesty_failure;
    failure.direct_witness_search = false;
    failure.radius = 64;
    failure.max_twists = 128;
    failure.random_tries = 400;
    const bool lattice =
            failure_detail_json(report) ==
            "{\n"
            "    \"p\": \"3\",\n"
            "    \"residue_degree\": 1,\n"
            "    \"search\": \"lattice\",\n"
            "    \"stage\": 3,\n"
            "    \"max_stage\": 3,\n"
            "    \"radius\": 64,\n"
            "    \"twists\": 128,\n"
            "    \"random_tries\": 400\n"
            "  }";

    // 2^89 - 1, a prime beyond any machine word, prints in full.
    if (!sflint::fmpz_set_str(sflint::FmpzRef(failure.p),
                              "618970019642690137449562111")) {
        return 1;
    }
    failure.residue_degree = 2;
    failure.direct_witness_search = true;
    failure.stage = 1;
    failure.max_stage = 3;
    failure.random_tries = 100;
    failure.factor_attempts = 1000;
    failure.element_steps = 2000;
    const bool t2 =
            failure_detail_json(report) ==
            "{\n"
            "    \"p\": \"618970019642690137449562111\",\n"
            "    \"residue_degree\": 2,\n"
            "    \"search\": \"t2\",\n"
            "    \"stage\": 1,\n"
            "    \"max_stage\": 3,\n"
            "    \"random_tries\": 100,\n"
            "    \"factor_attempts\": 1000,\n"
            "    \"element_steps\": 2000\n"
            "  }";

    report.failure_reason = "unsupported_certification_request";
    const bool other_reason = failure_detail_json(report) == "null";
    report.failure_reason = nullptr;
    const bool no_reason = failure_detail_json(report) == "null";

    HonestyCheckpointRun restart;
    if (!run_honesty_checkpoint(restart, true, true)) {
        return 1;
    }
    restart.audit.failure_stage = silex::detail::ClassUnitStage::total;
    restart.audit.failure_reason =
            silex::detail::class_unit_computation_failure_reason(
                    restart.audit);
    const bool generic = failure_detail_json(restart.audit) == "null";
    return checkpoint && lattice && t2 && other_reason && no_reason &&
                           generic
            ? 0
            : 1;
}

// A candidate or relation cap that stopped the relation search names the
// reason relation_search_budget_exhausted, and failure_detail names the cap
// by its option name, after the generic failure's cause when a validation
// ran.  Without a cap the reason stays generic; the size
// limit and an unwitnessed prime keep their own reasons.
int test_budget_failure_reason_and_detail_json() {
    using silex::detail::RelationSearchBudgetCap;
    silex::detail::ClassUnitTransactionReport report;
    const auto reason = [&report]() {
        report.failure_stage = silex::detail::ClassUnitStage::total;
        report.failure_reason =
                silex::detail::class_unit_computation_failure_reason(report);
        return std::string(report.failure_reason);
    };
    const bool generic = reason() == "class_unit_computation_failed" &&
                         failure_detail_json(report) == "null";

    report.relation_search_budget_cap = RelationSearchBudgetCap::max_candidates;
    const bool candidates =
            reason() == "relation_search_budget_exhausted" &&
            failure_detail_json(report) ==
                    "{\n"
                    "    \"cap\": \"max_candidates\"\n"
                    "  }";
    report.relation_search_budget_cap = RelationSearchBudgetCap::max_relations;
    const bool relations =
            reason() == "relation_search_budget_exhausted" &&
            failure_detail_json(report) ==
                    "{\n"
                    "    \"cap\": \"max_relations\"\n"
                    "  }";
    // The cap adds to the generic failure's cause; it does not replace it.
    report.last_validation_outcome = "index_one_publication_failed";
    const bool cause_kept =
            reason() == "relation_search_budget_exhausted" &&
            failure_detail_json(report) ==
                    "{\n"
                    "    \"cause\": \"last unit-group validation outcome: "
                    "index_one_publication_failed\",\n"
                    "    \"cap\": \"max_relations\"\n"
                    "  }";

    report.factor_base_honesty_failure.recorded = true;
    const bool honesty_first = reason() == "factor_base_honesty_unwitnessed";
    report.factor_base_size_limit_exceeded = true;
    const bool size_limit_first =
            reason() == "factor_base_size_limit_exceeded" &&
            failure_detail_json(report) == "null";

    report.reset();
    const bool cleared = report.relation_search_budget_cap ==
                                 RelationSearchBudgetCap::none &&
                         reason() == "class_unit_computation_failed";
    return generic && candidates && relations && cause_kept &&
                           honesty_first && size_limit_first && cleared
            ? 0
            : 1;
}

// The LLL relation slice of the proven continuation refuses to run once the
// relation count has reached the relation cap, and notes that cap in the
// run's transaction, so a proven computation that then gives up reports
// relation_search_budget_exhausted.  A refusal for an invalid argument is
// not a cap stop and notes nothing.
int test_lll_relation_slice_cap_refusal_is_noted() {
    using Access = silex::detail::ClassGroupRelationSearchAccess;
    using silex::detail::RelationSearchBudgetCap;
    FieldSetup setup = imaginary_quadratic_fixture(21);
    sflint::Fmpz bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    silex::ClassGroupContext context(setup.maximal_order);
    if (!context.is_defined() ||
        !context.build_factor_base(sflint::FmpzConstRef(bound)) ||
        context.factor_base() == nullptr) {
        return 1;
    }
    // The slice needs a published presentation: (2) = P2^2.
    silex::Relation relation(*context.factor_base());
    silex::Element two(setup.field);
    if (!two.set_si(2) || !relation.set_generator(two) ||
        !context.append_relation(relation) ||
        !context.publish_presentation() || context.relation_count() != 1) {
        return 1;
    }
    silex::detail::ClassUnitTransactionReport audit;
    silex::detail::ClassUnitTransactionContext transaction{
            audit, nullptr, {}, {}, false};
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               &transaction);
    silex::detail::ClassGroupRelationOptions options =
            silex::detail::class_group_relation_options(
                    setup.maximal_order, silex::ClassGroupComputeOptions{});
    options.requested_certification = silex::CertificationMode::proven;

    silex::detail::ClassGroupRelationOptions invalid = options;
    invalid.max_candidates = -1;
    const bool invalid_refused =
            !Access::extend_lll_relation_slice(
                    context, setup.maximal_order, sflint::FmpzConstRef(bound),
                    invalid) &&
            transaction.relation_search_budget_cap_reached ==
                    RelationSearchBudgetCap::none;

    options.max_relations = context.relation_count();
    const bool capped_refused =
            !Access::extend_lll_relation_slice(
                    context, setup.maximal_order, sflint::FmpzConstRef(bound),
                    options) &&
            context.relation_count() == 1 &&
            transaction.relation_search_budget_cap_reached ==
                    RelationSearchBudgetCap::max_relations;

    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return invalid_refused && capped_refused ? 0 : 1;
}

// Escalation runs only where an unwitnessed prime is final: a caller with a
// factor-base restart or a fallback route keeps the stage-0 caps of the
// reference be_honest.
int test_escalation_only_without_recovery() {
    return search::factor_base_witness_escalation(true) ==
                           search::FactorBaseWitnessEscalation::none &&
                   search::factor_base_witness_escalation(false) ==
                           search::FactorBaseWitnessEscalation::bounded
            ? 0
            : 1;
}

// Runs one honesty scan with the capped T2 search (direct) or the
// scalar/lattice search, with no escalation.
bool run_fallback_scan(bool& honest,
                       search::FactorBaseHonestyScanAudit& audit,
                       const silex::FactorBase& base,
                       ulong active,
                       ulong required,
                       const silex::detail::SubfactorBaseSchedule* schedule,
                       bool direct) noexcept {
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), active);
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), required);
    honest = !honest;
    return search::factor_base_honesty_check(
            honest, base, sflint::FmpzConstRef(active_bound),
            sflint::FmpzConstRef(required_bound), schedule, UWORD(42), direct,
            128, search::FactorBaseWitnessEscalation::none, nullptr, &audit);
}

// The capped T2 witness search keeps the scalar/lattice search's outcomes
// where the per-candidate test or the twist schedule cannot be used: a
// direct-evaluation error falls back to the full factorization, and a null
// schedule or a schedule entry outside the base means no further twists, so
// the scan still completes with the same honesty verdict instead of failing.
int test_direct_witness_search_keeps_lattice_fallbacks() {
    using silex::detail::SubfactorBaseSchedule;
    // Z[x]/(x^3 - x - 1) (discriminant -23) is the maximal order.  The hook
    // marks it maximal for the base build and then withdraws that, so the
    // direct test cannot evaluate any candidate while the equation-order
    // prime decomposition stays available.  The T2 search must still
    // complete the scan, through the full-factorization fallback, with the
    // scalar/lattice search's verdict.
    {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_zero(polynomial);
        sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
        const silex::NumberField field = silex::test::field_by_polynomial(
                sflint::FmpqPolyConstRef(polynomial));
        silex::Order order = silex::test::equation_order(field);
        if (!order.is_equation_order()) {
            return 1;
        }
        silex::detail::OrderAccess::set_maximality_unchecked(order, true);
        sflint::Fmpz active;
        sflint::fmpz_set_ui(active, 7);
        silex::FactorBase base(order);
        if (!base.build_relation_completion_base(sflint::FmpzConstRef(active)) ||
            base.length() == 0) {
            return 1;
        }
        SubfactorBaseSchedule schedule;
        schedule.subfactor_base.push_back(0);
        schedule.excluded_indices.assign(
                static_cast<std::size_t>(base.length()), 0);
        silex::detail::OrderAccess::set_maximality_unchecked(order, false);
        bool lattice_honest = false;
        bool direct_honest = false;
        search::FactorBaseHonestyScanAudit lattice_audit;
        search::FactorBaseHonestyScanAudit direct_audit;
        if (!run_fallback_scan(lattice_honest, lattice_audit, base, 7, 30,
                               &schedule, false) ||
            !run_fallback_scan(direct_honest, direct_audit, base, 7, 30,
                               &schedule, true) ||
            direct_honest != lattice_honest ||
            direct_audit.witness_targets == 0 ||
            direct_audit.witness_targets != lattice_audit.witness_targets ||
            direct_audit.witnessed_targets !=
                    lattice_audit.witnessed_targets ||
            direct_audit.rational_prime_checks !=
                    lattice_audit.rational_prime_checks ||
            direct_audit.unwitnessed.recorded !=
                    lattice_audit.unwitnessed.recorded ||
            !sflint::fmpz_equal(direct_audit.unwitnessed.p,
                                lattice_audit.unwitnessed.p)) {
            return 1;
        }
    }
    // Q(sqrt(-5)) on the empty norm-at-most-1 base: the prime above 2 has no
    // witness.  A null schedule, an excluded entry and an entry past the end
    // of the base each end the search with no witness.
    {
        FieldSetup setup = imaginary_quadratic_fixture(5);
        silex::FactorBase base(setup.maximal_order);
        if (!base.is_defined() || base.length() != 0) {
            return 1;
        }
        SubfactorBaseSchedule past_end;
        past_end.subfactor_base.push_back(0);
        past_end.excluded_indices.assign(1, 0);
        SubfactorBaseSchedule excluded;
        excluded.subfactor_base.push_back(0);
        excluded.excluded_indices.assign(1, 1);
        for (const SubfactorBaseSchedule* schedule :
             {static_cast<const SubfactorBaseSchedule*>(nullptr),
              static_cast<const SubfactorBaseSchedule*>(&past_end),
              static_cast<const SubfactorBaseSchedule*>(&excluded)}) {
            for (bool direct : {false, true}) {
                bool honest = true;
                search::FactorBaseHonestyScanAudit audit;
                if (!run_fallback_scan(honest, audit, base, 1, 2, schedule,
                                       direct) ||
                    honest || audit.witness_targets != 1 ||
                    audit.witnessed_targets != 0 ||
                    !audit.unwitnessed.recorded ||
                    audit.unwitnessed.direct_witness_search != direct) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

}  // namespace

FieldSetup maximal_setup(std::initializer_list<slong> coefficients) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    slong degree = 0;
    for (slong c : coefficients) {
        sflint::fmpq_poly_set_coeff_si(polynomial, degree++, c);
    }
    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(equation_order) ||
        !setup.maximal_order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

// The prime ideals above p, in decomposition order.
std::vector<silex::PrimeIdeal> primes_above(const silex::Order& order,
                                            slong p) {
    sflint::Fmpz rational_prime;
    sflint::fmpz_set_si(sflint::FmpzRef(rational_prime), p);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order,
                                  sflint::FmpzConstRef(rational_prime)));
    std::vector<silex::PrimeIdeal> out;
    for (slong i = 0; i < primes.size(); ++i) {
        out.emplace_back(order);
        assert(primes.at(i) != nullptr && out.back().set(*primes.at(i)));
    }
    return out;
}

silex::FactorBase base_of(const silex::Order& order,
                          std::initializer_list<const silex::PrimeIdeal*> primes) {
    silex::FactorBase base(order);
    for (const silex::PrimeIdeal* prime : primes) {
        assert(silex::detail::TestHookAccess::append_prime(
                base, *prime));
    }
    return base;
}

// Every witness predicate, direct and reference, against the expected
// classification, with the expensive debug differentials enabled, and the
// integral-scalar predicate too when alpha is the nonzero scalar_value.  When
// through_lacking is set, the witness is accepted only through the lacking
// ideal: the rule without that exception refuses it.
void check_lacking_prime_witness(const silex::FactorBase& base,
                                 const silex::PrimeIdeal& required,
                                 const silex::OrderElement& alpha,
                                 bool expected,
                                 bool through_lacking,
                                 slong scalar_value = 0) {
    namespace search = silex::detail::relation_search;
    using Mode = search::FactorBaseWitnessPredicate;
    const silex::Order& order = *base.parent();
    silex::Ideal principal(order);
    assert(principal.set_principal(alpha));
    bool reference = !expected;
    assert(full_factorization_required_prime_result(reference, principal, base,
                                                    required));
    assert(reference == expected);
    bool strict = true;
    assert(full_factorization_required_prime_result(strict, principal, base,
                                                    required, true));
    assert(strict == (expected && !through_lacking));

    silex::DiagnosticsContext diagnostics;
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    slong failures = 0;
    diagnostics.debug_failure_user = &failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    bool ideal_matches = !expected;
    bool element_matches = !expected;
    assert(silex::detail::ideal_factor_over_base_with_required_prime(
            ideal_matches, principal, base, required, &diagnostics));
    assert(silex::detail::order_element_factor_over_base_with_required_prime(
            element_matches, alpha, base, required, &diagnostics));
    assert(ideal_matches == expected && element_matches == expected);
    for (Mode mode : {Mode::selected, Mode::full_factorization,
                      Mode::order_element_direct}) {
        search::FactorBaseWitnessAudit audit;
        assert(search::factor_base_principal_witness(
                       base, required, alpha, mode, &diagnostics, &audit) ==
               expected);
        assert(audit.used_reference == (mode == Mode::full_factorization));
    }

    if (scalar_value != 0) {
        sflint::Fmpz integer;
        sflint::fmpz_set_si(sflint::FmpzRef(integer), scalar_value);
        bool scalar_matches = !expected;
        assert(silex::detail::scalar_factor_over_base_with_required_prime(
                scalar_matches, sflint::FmpzConstRef(integer), base, required,
                &diagnostics));
        assert(scalar_matches == expected);
        for (Mode mode : {Mode::selected, Mode::order_element_direct,
                          Mode::full_factorization}) {
            assert(search::factor_base_scalar_witness(
                           base, required, sflint::FmpzConstRef(integer), mode,
                           &diagnostics) == expected);
        }
    }
    assert(failures == 0);
}

// A witness may carry valuation at the one prime ideal above q that the
// base lacks when e = 1, the base holds every other ideal above q, and q is
// not the required prime's rational prime: qO_K = Q_1^e_1 ... Q_g^e_g with
// e_g = 1 puts [Q_g] in the subgroup of the others, the step
// select_factor_base_proof_targets uses (2.17.3 buch2.c:be_honest and
// bnftestprimes skip the last unramified ideal).  Each refusal below
// breaks exactly one condition.
int test_witness_through_lacking_conjugate() {
    // Z[i]: 5 and 13 split, class number 1.
    {
        FieldSetup setup = maximal_setup({1, 0, 1});
        const silex::Order& order = setup.maximal_order;
        assert(order.is_defined());
        const auto p13 = primes_above(order, 13);
        const auto p5 = primes_above(order, 5);
        assert(p13.size() == 2 && p5.size() == 2);
        // A non-scalar a + b i with (alpha) = P13 * Q5, Q5 the ideal above 5
        // outside the base {Q5'}.
        const silex::FactorBase base = base_of(order, {&p5[0]});
        silex::OrderElement alpha(order);
        sflint::FmpzMat coordinates(1, 2);
        bool found = false;
        for (slong a = -8; a <= 8 && !found; ++a) {
            for (slong b = 1; b <= 8 && !found; ++b) {
                if (a * a + b * b != 65) continue;
                sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(coordinates), 0, 0), a);
                sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(coordinates), 0, 1), b);
                assert(alpha.set_coordinates(
                        sflint::FmpzMatConstRef(coordinates)));
                silex::Ideal principal(order);
                slong v13 = -1, v5 = -1;
                assert(principal.set_principal(alpha));
                assert(p13[0].valuation(v13, principal));
                assert(p5[1].valuation(v5, principal));
                found = v13 == 1 && v5 == 1;
            }
        }
        assert(found);
        check_lacking_prime_witness(base, p13[0], alpha, true, true);
        // The base holds no ideal above 5: refused.
        check_lacking_prime_witness(base_of(order, {&p13[1]}), p13[0], alpha,
                                    false, false);
        // The scalar 65 = P13 P13' Q5 Q5' with P13' and Q5' in the base.
        assert(alpha.set_si(65));
        check_lacking_prime_witness(base_of(order, {&p13[1], &p5[0]}), p13[0],
                                    alpha, true, true, 65);
        // P13' lacking as well: the second lacking ideal lies above the
        // required prime's own rational prime 13, refused.
        check_lacking_prime_witness(base, p13[0], alpha, false, false, 65);
        // 3 is inert, so (3 beta) = P13 * (3) with (3) the only ideal above
        // 3 and outside the base: refused, as in buch2.c:can_factor.
        found = false;
        for (slong a = -4; a <= 4 && !found; ++a) {
            for (slong b = 1; b <= 4 && !found; ++b) {
                if (a * a + b * b != 13) continue;
                sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(coordinates), 0, 0), 3 * a);
                sflint::fmpz_set_si(sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(coordinates), 0, 1), 3 * b);
                assert(alpha.set_coordinates(
                        sflint::FmpzMatConstRef(coordinates)));
                silex::Ideal principal(order);
                slong v13 = -1;
                assert(principal.set_principal(alpha));
                assert(p13[0].valuation(v13, principal));
                found = v13 == 1;
            }
        }
        assert(found);
        const auto p3 = primes_above(order, 3);
        assert(p3.size() == 1 && p3[0].ramification_index() == 1 &&
               p3[0].residue_degree() == 2);
        check_lacking_prime_witness(base, p13[0], alpha, false, false);
    }
    // x^3 - x - 1 (discriminant -23, class number 1): 23 = P^2 Q and
    // 5 = R1 R2 with f(R1) = 1, f(R2) = 2.  Required R1, scalar 5 * 23.
    {
        FieldSetup setup = maximal_setup({-1, -1, 0, 1});
        const silex::Order& order = setup.maximal_order;
        assert(order.is_defined());
        const auto p5 = primes_above(order, 5);
        const auto p23 = primes_above(order, 23);
        assert(p5.size() == 2 && p23.size() == 2);
        const silex::PrimeIdeal& r1 = p5[0].residue_degree() == 1 ? p5[0] : p5[1];
        const silex::PrimeIdeal& r2 = p5[0].residue_degree() == 1 ? p5[1] : p5[0];
        const silex::PrimeIdeal& ramified =
                p23[0].ramification_index() == 2 ? p23[0] : p23[1];
        const silex::PrimeIdeal& unramified =
                p23[0].ramification_index() == 2 ? p23[1] : p23[0];
        assert(r1.residue_degree() == 1 && r2.residue_degree() == 2);
        assert(ramified.ramification_index() == 2 &&
               unramified.ramification_index() == 1);
        silex::OrderElement alpha(order);
        assert(alpha.set_si(5 * 23));
        // Lacking Q (e = 1): accepted.
        check_lacking_prime_witness(base_of(order, {&r2, &ramified}), r1, alpha,
                                    true, true, 5 * 23);
        // Lacking P (e = 2): refused.
        check_lacking_prime_witness(base_of(order, {&r2, &unramified}), r1,
                                    alpha, false, false, 5 * 23);
    }
    // x^3 - 2: 31 = P1 P2 P3 and 5 = R1 R2 with f(R1) = 1, f(R2) = 2.
    {
        FieldSetup setup = cubic_x3_minus_2();
        const silex::Order& order = setup.maximal_order;
        assert(order.is_defined());
        const auto p5 = primes_above(order, 5);
        const auto p31 = primes_above(order, 31);
        assert(p5.size() == 2 && p31.size() == 3);
        const silex::PrimeIdeal& r1 = p5[0].residue_degree() == 1 ? p5[0] : p5[1];
        const silex::PrimeIdeal& r2 = p5[0].residue_degree() == 1 ? p5[1] : p5[0];
        assert(r1.residue_degree() == 1 && r2.residue_degree() == 2);
        silex::OrderElement alpha(order);
        assert(alpha.set_si(5 * 31));
        // One lacking ideal above 31: accepted.
        check_lacking_prime_witness(base_of(order, {&r2, &p31[0], &p31[1]}), r1,
                                    alpha, true, true, 5 * 31);
        // Two lacking ideals above 31: refused.
        check_lacking_prime_witness(base_of(order, {&r2, &p31[1]}), r1, alpha,
                                    false, false, 5 * 31);
        // The lacking ideal above 5 has f = 2, e = 1: accepted.
        check_lacking_prime_witness(base_of(order, {&r1, &p31[1], &p31[2]}),
                                    p31[0], alpha, true, true, 5 * 31);
        // 31 = P1 P2 P3 with P2 in the base and P3 lacking: P3 lies above the
        // required prime's own rational prime, refused.
        assert(alpha.set_si(31));
        check_lacking_prime_witness(base_of(order, {&p31[1]}), p31[0], alpha,
                                    false, false, 31);
    }
    return 0;
}

// The honesty scan skips a rational prime p of a maximal quadratic order
// when the Kronecker symbol (D/p) is -1.  The symbol must classify every p
// exactly as decompose_prime does: -1 for one prime of e = 1, f = 2
// (inert), 0 for one prime of e = 2 (ramified), +1 for two primes of
// e = f = 1 (split).  The fields cover p = 2 split (D = 17, -7), inert
// (D = 5, -3), and ramified (D = -4, 40), and odd p dividing D (D = 40,
// -4036 = -4 * 1009, 100000005 = 3 * 5 * 7 * 952381).  Polynomials other
// than x^2 - d with d squarefree use the generic number-field backend, and
// x^2 - 20 has the inert p = 2 dividing the index of Z[sqrt(20)].
int test_honesty_quadratic_kronecker_matches_decomposition() {
    struct Case {
        std::vector<slong> coefficients;
        slong discriminant;
    };
    const std::vector<Case> cases = {
            {{1, 0, 1}, -4},    {{-17, 0, 1}, 17},   {{-5, 0, 1}, 5},
            {{7, 0, 1}, -7},    {{3, 0, 1}, -3},     {{-10, 0, 1}, 40},
            {{2, 1, 1}, -7},    {{-20, 0, 1}, 5},    {{1009, 0, 1}, -4036},
            {{-25000001, -1, 1}, 100000005},
    };
    for (const Case& test_case : cases) {
        const FieldSetup setup = maximal_order_setup(test_case.coefficients);
        if (!setup.maximal_order.is_defined()) {
            return 1;
        }
        sflint::Fmpz discriminant;
        if (!silex::detail::relation_search::
                    factor_base_honesty_quadratic_discriminant(
                            sflint::FmpzRef(discriminant),
                            setup.maximal_order) ||
            ::fmpz_cmp_si(discriminant.raw(),
                                test_case.discriminant) != 0) {
            return 2;
        }
        sflint::Fmpz p;
        sflint::fmpz_set_ui(sflint::FmpzRef(p), UWORD(2));
        while (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(p), UWORD(200)) <= 0) {
            silex::PrimeIdealList primes;
            if (!silex::decompose_prime(primes, setup.maximal_order,
                                        sflint::FmpzConstRef(p))) {
                return 3;
            }
            const int symbol = sflint::fmpz_kronecker(
                    sflint::FmpzConstRef(discriminant),
                    sflint::FmpzConstRef(p));
            const bool matches =
                    symbol < 0   ? primes.size() == 1 &&
                                         prime_has_ef(primes, 0, 1, 2)
                    : symbol == 0 ? primes.size() == 1 &&
                                         prime_has_ef(primes, 0, 2, 1)
                                  : primes.size() == 2 &&
                                         prime_has_ef(primes, 0, 1, 1) &&
                                         prime_has_ef(primes, 1, 1, 1);
            if (!matches) {
                return 4;
            }
            sflint::fmpz_nextprime(sflint::FmpzRef(p),
                                   sflint::FmpzConstRef(p), true);
        }
    }

    // No pre-check outside the maximal-quadratic route: a nonmaximal
    // quadratic order and a cubic maximal order.
    const silex::NumberField nonmaximal_field =
            maximal_order_setup({-5, 0, 1}).field;
    const silex::Order nonmaximal =
            silex::test::equation_order(nonmaximal_field);
    const FieldSetup cubic = maximal_order_setup({-2, 0, 0, 1});
    sflint::Fmpz unused;
    if (silex::detail::relation_search::
                factor_base_honesty_quadratic_discriminant(
                        sflint::FmpzRef(unused), nonmaximal) ||
        !cubic.maximal_order.is_defined() ||
        silex::detail::relation_search::
                factor_base_honesty_quadratic_discriminant(
                        sflint::FmpzRef(unused), cubic.maximal_order)) {
        return 5;
    }

    // D = 5: p = 2 and 3 are inert, so an empty factor base is honest up to
    // 3 with both primes counted and no witness target.  At bound 5 the
    // ramified prime above 5 is still decomposed and is the one target; it
    // is principal (class number 1), so the empty base witnesses it.
    const FieldSetup golden = maximal_order_setup({-1, -1, 1});
    silex::FactorBase base(golden.maximal_order);
    if (!base.is_defined() || base.length() != 0) {
        return 6;
    }
    for (const ulong bound : {UWORD(3), UWORD(5)}) {
        sflint::Fmpz required_bound;
        sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), bound);
        bool honest = false;
        silex::detail::relation_search::FactorBaseHonestyScanAudit audit;
        if (!silex::detail::relation_search::factor_base_honesty_check(
                    honest, base, sflint::FmpzConstRef(required_bound),
                    sflint::FmpzConstRef(required_bound), nullptr, UWORD(7),
                    false, 128,
                    silex::detail::relation_search::
                            FactorBaseWitnessEscalation::none,
                    nullptr, &audit) ||
            !honest ||
            audit.rational_prime_checks != (bound == UWORD(3) ? 2 : 3) ||
            audit.witness_targets != (bound == UWORD(3) ? 0 : 1) ||
            audit.witnessed_targets != audit.witness_targets) {
            return 7;
        }
    }
    return 0;
}

// Runs one honesty scan on the capped T2 search (the route no production
// scan with a null schedule uses on its own) over `coefficients`, whose
// maximal order has base norm bound `active` (1 keeps the empty base).
bool run_direct_scan(bool& honest,
                     search::FactorBaseHonestyScanAudit& audit,
                     HonestyLogCapture& capture,
                     const std::vector<slong>& coefficients,
                     ulong active,
                     ulong required,
                     search::FactorBaseWitnessEscalation escalation) noexcept {
    FieldSetup setup = maximal_order_setup(coefficients);
    if (!setup.maximal_order.is_defined()) {
        return false;
    }
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), active);
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), required);
    silex::FactorBase base(setup.maximal_order);
    if (!base.is_defined() ||
        (active > 1 && !base.build_relation_completion_base(
                               sflint::FmpzConstRef(active_bound)))) {
        return false;
    }
    silex::DiagnosticsContext diagnostics;
    enable_honesty_log(diagnostics, capture);
    honest = true;
    return search::factor_base_honesty_check(
            honest, base, sflint::FmpzConstRef(active_bound),
            sflint::FmpzConstRef(required_bound), nullptr, UWORD(0), true,
            128, escalation, &diagnostics, &audit);
}

// Escalation on the capped T2 route (use_direct_required_prime_witness).
// Z[x]/(x^3 - 15) has class number 2 and the prime above 2 (ramified,
// 2 O = P^3) is not principal, so with the empty base no search can witness
// it.  Without escalation the scan stops at the stage-0 caps (50 random
// tries, 500 factored candidates and 10^6 enumeration steps per ideal); with
// bounded escalation it runs three more stages, doubling each cap, and
// reports the final caps in the audit and the log.  The null schedule means
// no twists, so one ideal is searched, and it stops at the candidate cap:
// factor_attempt_cap_hits counts that stop and element_step_cap_hits stays 0.
int test_direct_route_escalation_caps_and_counters() {
    const std::vector<slong> cubic = {-15, 0, 0, 1};
    {
        bool honest = true;
        search::FactorBaseHonestyScanAudit audit;
        HonestyLogCapture capture;
        if (!run_direct_scan(honest, audit, capture, cubic, 1, 30,
                             search::FactorBaseWitnessEscalation::none) ||
            honest || audit.max_search_stage != 0 ||
            audit.witness_targets != 1 || audit.witnessed_targets != 0 ||
            !audit.unwitnessed.recorded ||
            !audit.unwitnessed.direct_witness_search ||
            audit.unwitnessed.stage != 0 || audit.unwitnessed.max_stage != 0 ||
            audit.unwitnessed.random_tries != 50 ||
            audit.unwitnessed.factor_attempts != 500 ||
            audit.unwitnessed.element_steps != 1000000 ||
            !logged_caps(capture,
                         "stage=0/0 search=t2 random_tries=50 "
                         "factor_attempts=500 element_steps=1000000 "
                         "ideals=1 factor_attempt_cap_hits=1 "
                         "element_step_cap_hits=0")) {
            return 1;
        }
    }
    {
        bool honest = true;
        search::FactorBaseHonestyScanAudit audit;
        HonestyLogCapture capture;
        if (!run_direct_scan(honest, audit, capture, cubic, 1, 30,
                             search::FactorBaseWitnessEscalation::bounded) ||
            honest || audit.max_search_stage !=
                              search::kFactorBaseHonestyEscalationStages ||
            audit.witnessed_targets != 0 ||
            !final_stage_t2_failure(audit.unwitnessed, 2) ||
            !logged_caps(capture,
                         "stage=3/3 search=t2 random_tries=400 "
                         "factor_attempts=4000 element_steps=8000000 "
                         "ideals=1 factor_attempt_cap_hits=1 "
                         "element_step_cap_hits=0")) {
            return 1;
        }
    }
    return 0;
}

// Which escalation mode each production caller passes to the honesty scan.
// The mode matters only when a prime has no witness, which no cheap field
// produces on these routes, so the test reads the per-scan call counts
// (TestHookAccess).  The relation-production prepass is the
// caller with no factor-base restart and no fallback: bounded escalation on
// the T2 route.  The relation-search driver, with a factor-base restart
// available, keeps the stage-0 reference flow: no escalation on the T2 route.
// Both pass the capped T2 search; no production caller passes the lattice
// search.
bool class_unit_honesty_call_counts(slong& direct_none,
                                    slong& direct_bounded,
                                    slong& lattice_total,
                                    const FieldSetup& setup,
                                    ulong bound) noexcept {
    if (!setup.maximal_order.is_defined()) {
        return false;
    }
    sflint::Fmpz factor_base_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), bound);
    silex::ClassGroupComputeOptions options;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::TestHookAccess::factor_base_honesty_call_counts_reset();
    if (!units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128)) {
        return false;
    }
    using E = search::FactorBaseWitnessEscalation;
    using Hook = silex::detail::TestHookAccess;
    direct_none = Hook::factor_base_honesty_call_count(true, E::none);
    direct_bounded =
            Hook::factor_base_honesty_call_count(true, E::bounded);
    lattice_total =
            Hook::factor_base_honesty_call_count(false, E::none) +
            Hook::factor_base_honesty_call_count(false, E::bounded);
    return true;
}

bool class_unit_honesty_call_counts(slong& direct_none,
                                    slong& direct_bounded,
                                    slong& lattice_total,
                                    const std::vector<slong>& coefficients,
                                    ulong bound) noexcept {
    return class_unit_honesty_call_counts(
            direct_none, direct_bounded, lattice_total,
            maximal_order_setup(coefficients), bound);
}

// The T2 witness enumeration's norm prefilter against the direct order-element
// test, candidate by candidate. Each candidate is back * x for x in an ideal
// the search enumerates (the required prime P, P times a rational scalar
// ideal after factor_base_honesty_primitive_part, and P times a large base
// prime power after factor_base_honesty_reduce_large_ideal), given by the
// order coordinates of x as in enumerate_required_prime_witness. The
// prefilter must decide every integral candidate, reject exactly the
// candidates the direct test rejects by their norm cofactor, and classify
// rational scalars as element_is_scalar_rational does.
struct T2PrefilterCounts {
    slong rejected = 0;
    slong kept = 0;
    slong witnesses = 0;
    slong scalars = 0;
    slong back_rejected = 0;
};

int check_t2_prefilter_ideal(const silex::FactorBase& base,
                             const silex::PrimeIdeal& required,
                             search::RequiredPrimeWitnessPrefilter& prefilter,
                             const silex::Ideal& ideal,
                             const silex::Element& back,
                             slong radius,
                             T2PrefilterCounts& counts) {
    const silex::Order& order = *base.parent();
    const slong n = order.degree();
    sflint::FmpzMat hnf(n, n);
    silex::lat::Lat lattice(n);
    silex::lat::Lat reduced(n);
    sflint::FmpzMat basis(n, n);
    if (!search::required_prime_witness_prefilter_set_back_multiplier(
                prefilter, back) ||
        !ideal.get_hnf(sflint::FmpzMatRef(hnf)) ||
        !lattice.set_basis(sflint::FmpzMatConstRef(hnf)) ||
        !lattice.lll_reduce(reduced) || reduced.nrows() != n ||
        !reduced.get_basis(sflint::FmpzMatRef(basis))) {
        return 1;
    }
    const bool back_is_one = back.equal_si(1);

    sflint::FmpzMat coordinates(1, n);
    silex::OrderElement x(order);
    silex::OrderElement candidate(order);
    silex::Element element(*order.parent());
    silex::Element product(*order.parent());
    const auto check = [&]() -> bool {
        bool scalar = true;
        bool reference_scalar = false;
        if (!x.set_coordinates(sflint::FmpzMatConstRef(coordinates)) ||
            !x.get_element(element) ||
            !silex::detail::relation_search::element_is_scalar_rational(
                    reference_scalar, element) ||
            !search::required_prime_witness_prefilter_is_scalar(
                    scalar, prefilter,
                    sflint::FmpzMatConstRef(coordinates)) ||
            scalar != reference_scalar) {
            return false;
        }
        if (scalar) {
            ++counts.scalars;
            return true;
        }
        bool matches = true;
        bool rejected = false;
        silex::detail::RequiredPrimeDirectAudit audit;
        if (!product.multiply(back, element) ||
            !candidate.set_element(product) ||
            !silex::detail::order_element_factor_over_base_with_required_prime(
                    matches, candidate, base, required, nullptr, &audit) ||
            !search::required_prime_witness_prefilter_rejects(
                    rejected, prefilter,
                    sflint::FmpzMatConstRef(coordinates)) ||
            rejected != audit.rejected_by_cofactor ||
            (rejected && matches)) {
            return false;
        }
        rejected ? ++counts.rejected : ++counts.kept;
        counts.back_rejected += rejected && !back_is_one ? 1 : 0;
        counts.witnesses += matches ? 1 : 0;
        return true;
    };

    // Every nonzero combination in [-radius, radius]^n of the reduced basis.
    std::vector<slong> digits(static_cast<std::size_t>(n), -radius);
    for (;;) {
        bool nonzero = false;
        for (slong j = 0; j < n; ++j) {
            sflint::fmpz_zero(
                    sflint::fmpz_mat_entry(sflint::FmpzMatRef(coordinates), 0,
                                           j));
        }
        for (slong i = 0; i < n; ++i) {
            const slong d = digits[static_cast<std::size_t>(i)];
            nonzero = nonzero || d != 0;
            for (slong j = 0; j < n; ++j) {
                fmpz_addmul_si(fmpz_mat_entry(coordinates.raw(), 0, j),
                               fmpz_mat_entry(basis.raw(), i, j), d);
            }
        }
        if (nonzero && !check()) {
            return 1;
        }
        slong i = 0;
        while (i < n && digits[static_cast<std::size_t>(i)] == radius) {
            digits[static_cast<std::size_t>(i)] = -radius;
            ++i;
        }
        if (i == n) {
            break;
        }
        ++digits[static_cast<std::size_t>(i)];
    }
    // Rational scalars k * 1, k in [-3, 3] (k = 0 skipped); is_scalar does not
    // depend on ideal membership.
    silex::OrderElement one(order);
    sflint::FmpzMat one_coordinates(1, n);
    if (!one.one() ||
        !one.get_coordinates(sflint::FmpzMatRef(one_coordinates))) {
        return 1;
    }
    for (slong k = -3; k <= 3; ++k) {
        if (k == 0) {
            continue;
        }
        for (slong j = 0; j < n; ++j) {
            fmpz_mul_si(fmpz_mat_entry(coordinates.raw(), 0, j),
                        fmpz_mat_entry(one_coordinates.raw(), 0, j), k);
        }
        bool scalar = false;
        if (!search::required_prime_witness_prefilter_is_scalar(
                    scalar, prefilter,
                    sflint::FmpzMatConstRef(coordinates)) ||
            !scalar) {
            return 1;
        }
    }
    return 0;
}

int check_t2_prefilter_field(const FieldSetup& setup,
                             slong base_bound,
                             std::initializer_list<slong> required_primes,
                             slong radius,
                             T2PrefilterCounts& counts) {
    const silex::Order& order = setup.maximal_order;
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), base_bound);
    silex::FactorBase base(order);
    search::RequiredPrimeWitnessPrefilter prefilter;
    silex::PrimeIdeal twist_prime(order);
    if (!order.is_maximal() || !base.build(sflint::FmpzConstRef(bound)) ||
        base.length() == 0 || !base.prime(twist_prime, base.length() - 1) ||
        !search::required_prime_witness_prefilter_init_base(
                prefilter, base, search::FactorBaseWitnessPredicate::selected)) {
        return 1;
    }
    for (slong p : required_primes) {
        silex::PrimeIdeal required(order);
        if (!first_prime_above(required, order, p) ||
            base.contains(required) ||
            !search::required_prime_witness_prefilter_set_required_prime(
                    prefilter, base, required)) {
            return 1;
        }
        silex::Ideal prime_ideal(order);
        silex::Ideal twist(order);
        silex::Ideal product(order);
        silex::Ideal scalar_ideal(order);
        silex::Element back(*order.parent());
        silex::Element scalar(*order.parent());
        silex::OrderElement scalar_generator(order);
        if (!required.get_ideal(prime_ideal) || !back.one() ||
            check_t2_prefilter_ideal(base, required, prefilter, prime_ideal,
                                     back, radius, counts) != 0) {
            return 1;
        }
        // P * 6O: the primitive part is P with back multiplier 6.
        if (!scalar.set_si(6) || !scalar_generator.set_element(scalar) ||
            !scalar_ideal.set_principal(scalar_generator) ||
            !product.multiply(prime_ideal, scalar_ideal) || !back.one() ||
            !search::factor_base_honesty_primitive_part(product, back) ||
            back.equal_si(1) ||
            check_t2_prefilter_ideal(base, required, prefilter, product, back,
                                     radius, counts) != 0) {
            return 1;
        }
        // P * Q^40 for a base prime Q, reduced by the large-ideal step.
        if (!twist_prime.get_ideal(twist) || !product.set(prime_ideal)) {
            return 1;
        }
        for (slong e = 0; e < 40; ++e) {
            silex::Ideal next(order);
            if (!next.multiply(product, twist)) {
                return 1;
            }
            product.swap(next);
        }
        if (!back.one() ||
            !search::factor_base_honesty_primitive_part(product, back) ||
            !search::factor_base_honesty_reduce_large_ideal(product, back,
                                                            128, nullptr) ||
            back.equal_si(1) ||
            check_t2_prefilter_ideal(base, required, prefilter, product, back,
                                     radius, counts) != 0) {
            return 1;
        }
    }
    return 0;
}

int test_t2_witness_prefilter_matches_direct() {
    // Corpus field real_quadratic_d1000000000009 with its proven factor-base
    // bound 293, an imaginary quadratic, a cubic, and two quartics.
    T2PrefilterCounts counts;
    const FieldSetup d1e12 = maximal_order_setup({-250000000002, -1, 1});
    if (check_t2_prefilter_field(d1e12, 293, {307, 311, 1009}, 25, counts) !=
                0 ||
        check_t2_prefilter_field(maximal_order_setup({7, 0, 1}), 7,
                                 {11, 23}, 20, counts) != 0 ||
        check_t2_prefilter_field(cubic_x3_minus_2(), 7, {11, 31}, 5,
                                 counts) != 0 ||
        check_t2_prefilter_field(maximal_order_setup({1, 0, -10, 0, 1}), 7,
                                 {23}, 3, counts) != 0 ||
        check_t2_prefilter_field(
                truncated_decomposition_counterexample_fixture(), 7, {11}, 3,
                counts) != 0) {
        return 1;
    }
    if (counts.rejected == 0 || counts.kept == 0 || counts.witnesses == 0 ||
        counts.back_rejected == 0) {
        return 2;
    }

    // Not ready for the full-factorization predicate or outside degrees 2 to
    // 4; such scans build every candidate.
    search::RequiredPrimeWitnessPrefilter prefilter;
    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 7);
    silex::FactorBase quadratic_base(d1e12.maximal_order);
    const FieldSetup quintic = ramified_last_fixture();
    silex::FactorBase quintic_base(quintic.maximal_order);
    if (!quadratic_base.build(sflint::FmpzConstRef(bound)) ||
        search::required_prime_witness_prefilter_init_base(
                prefilter, quadratic_base,
                search::FactorBaseWitnessPredicate::full_factorization) ||
        prefilter.norm.base_ready || prefilter.norm.ready ||
        !quintic_base.build(sflint::FmpzConstRef(bound)) ||
        search::required_prime_witness_prefilter_init_base(
                prefilter, quintic_base,
                search::FactorBaseWitnessPredicate::selected) ||
        prefilter.norm.base_ready) {
        return 3;
    }
    return 0;
}

// Whole T2-route honesty scans with the prefilter (the direct predicates)
// against the same scans without it (the full-factorization predicate): the
// same verdicts, audits and random-state stream, with the prefilter
// rejecting candidates only on the direct predicates. In a debug-check
// build every rejection is also checked against the direct test.
int test_t2_scan_prefilter_matches_full_factorization() {
    using Mode = search::FactorBaseWitnessPredicate;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    slong debug_failures = 0;
    diagnostics.debug_failure_user = &debug_failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    struct Row {
        std::vector<slong> coefficients;
        ulong active;
        ulong required;
    };
    const std::vector<Row> rows = {
            {{-250000000002, -1, 1}, 293, 700},
            {{-6236, 0, 1}, 3, 60},
            {{-15, 0, 0, 1}, 2, 30},
            {{-2, 0, 0, 1}, 7, 40},
            {{1, 0, -10, 0, 1}, 7, 30},
    };
    slong rejections = 0;
    for (const Row& row : rows) {
        FieldSetup setup = maximal_order_setup(row.coefficients);
        sflint::Fmpz active, required;
        sflint::fmpz_set_ui(sflint::FmpzRef(active), row.active);
        sflint::fmpz_set_ui(sflint::FmpzRef(required), row.required);
        silex::FactorBase base(setup.maximal_order);
        if (!setup.maximal_order.is_defined() ||
            !base.build_relation_completion_base(
                    sflint::FmpzConstRef(active)) ||
            base.length() == 0) {
            return 1;
        }
        // Twists over the whole base, so later ideals carry a back
        // multiplier.
        silex::detail::SubfactorBaseSchedule schedule;
        for (slong i = 0; i < base.length(); ++i) {
            schedule.subfactor_base.push_back(i);
        }
        schedule.excluded_indices.assign(
                static_cast<std::size_t>(base.length()), 0);
        // Bounded escalation only raises the caps; it is left out to keep
        // the full-factorization scans short.
        const auto escalation = search::FactorBaseWitnessEscalation::none;
        bool expected_honest = false;
        search::FactorBaseHonestyScanAudit expected;
        const bool expected_ok = search::factor_base_honesty_check(
                expected_honest, base, sflint::FmpzConstRef(active),
                sflint::FmpzConstRef(required), &schedule, UWORD(42),
                true, 128, escalation, nullptr, &expected,
                Mode::full_factorization);
        if (!expected_ok || expected.t2_norm_prefilter_rejections != 0) {
            return 2;
        }
        for (Mode mode : {Mode::selected, Mode::order_element_direct}) {
            bool honest = !expected_honest;
            search::FactorBaseHonestyScanAudit audit;
            const bool ok = search::factor_base_honesty_check(
                    honest, base, sflint::FmpzConstRef(active),
                    sflint::FmpzConstRef(required), &schedule, UWORD(42),
                    true, 128, escalation, &diagnostics, &audit, mode);
            if (!ok || honest != expected_honest ||
                audit.rational_prime_checks !=
                        expected.rational_prime_checks ||
                audit.witness_targets != expected.witness_targets ||
                audit.witnessed_targets != expected.witnessed_targets ||
                audit.escalated_witnessed_targets !=
                        expected.escalated_witnessed_targets ||
                audit.max_search_stage != expected.max_search_stage ||
                audit.final_random_state != expected.final_random_state ||
                audit.unwitnessed.recorded !=
                        expected.unwitnessed.recorded ||
                !sflint::fmpz_equal(audit.unwitnessed.p,
                                    expected.unwitnessed.p)) {
                return 3;
            }
            rejections += audit.t2_norm_prefilter_rejections;
        }
    }
    // The d1000000000009 and x^2 - 6236 rows reject candidates.
    return debug_failures == 0 && rejections > 0 ? 0 : 4;
}

int test_routes_pass_intended_escalation_mode() {
    slong direct_none = 0;
    slong direct_bounded = 0;
    slong lattice = 0;
    // Prepass: Z[x]/(x^3 - 15), base bound 2.
    if (!class_unit_honesty_call_counts(direct_none, direct_bounded, lattice,
                                        {-15, 0, 0, 1}, 2) ||
        direct_bounded != 1 || direct_none != 0 || lattice != 0) {
        return 1;
    }
    // Relation-search driver with a factor-base restart available:
    // Z[x]/(x^2 - 6236) (real quadratic, d = 4 * 1559), base bound 3.
    if (!class_unit_honesty_call_counts(direct_none, direct_bounded, lattice,
                                        {-6236, 0, 1}, 3) ||
        direct_none != 1 || direct_bounded != 0 || lattice != 0) {
        return 2;
    }
    return 0;
}

// The imaginary quadratic backend runs no honesty scan: a `proven` class
// group of the maximal order of Q(sqrt(-4027)) reaches the generation bound
// with the default plan, and with the column limit lowered the plan keeps
// the compact base on the structure-certificate route.  A compact base
// without that route fails closed in the backend, and no plan produces one.
int test_quadratic_backend_runs_no_honesty_scan() {
    using Hook = silex::detail::TestHookAccess;
    slong direct_none = 0;
    slong direct_bounded = 0;
    slong lattice = 0;
    FieldSetup setup;
    setup.field = silex::test::quadratic_field(-4027);
    setup.maximal_order = silex::Order(setup.field);
    if (!setup.maximal_order.maximal_order(
                silex::test::equation_order(setup.field))) {
        return 1;
    }
    if (!class_unit_honesty_call_counts(direct_none, direct_bounded, lattice,
                                        setup, 3) ||
        direct_none + direct_bounded + lattice != 0) {
        return 2;
    }
    Hook::set_quadratic_structure_certificate_column_limit(1);
    const bool ok = class_unit_honesty_call_counts(
            direct_none, direct_bounded, lattice, setup, 3);
    Hook::set_quadratic_structure_certificate_column_limit(0);
    return ok && direct_none + direct_bounded + lattice == 0 ? 0 : 3;
}

// Lowering the T2 element-step cap makes a direct witness enumeration stop at
// the cap: the hit counter is zero at the default cap and positive when the
// cap is one.
int test_element_step_cap_hits_observed() {
    using Hook = silex::detail::TestHookAccess;
    FieldSetup setup = maximal_order_setup({-15, 0, 0, 1});
    if (!setup.maximal_order.is_defined()) {
        return 1;
    }
    sflint::Fmpz active, required;
    sflint::fmpz_set_ui(sflint::FmpzRef(active), 2);
    sflint::fmpz_set_ui(sflint::FmpzRef(required), 30);
    silex::FactorBase base(setup.maximal_order);
    if (!base.build_relation_completion_base(sflint::FmpzConstRef(active))) {
        return 2;
    }
    slong hits[2] = {0, 0};
    for (int i = 0; i < 2; ++i) {
        Hook::set_factor_base_honesty_element_step_cap(i == 0 ? 0 : 1);
        Hook::factor_base_honesty_element_step_cap_hits_reset();
        bool honest = false;
        const bool ok = search::factor_base_honesty_check(
                honest, base, sflint::FmpzConstRef(active),
                sflint::FmpzConstRef(required), nullptr, UWORD(42), true, 128,
                search::FactorBaseWitnessEscalation::none, nullptr);
        hits[i] = Hook::factor_base_honesty_element_step_cap_hits();
        Hook::set_factor_base_honesty_element_step_cap(0);
        if (!ok) {
            return 3;
        }
    }
    return hits[0] == 0 && hits[1] > 0 ? 0 : 4;
}

// The T2 double-row setup (T2EnumerationGram::double_rows) against the Arb
// setup on one ideal, both without an embedding cache so both start at the
// first rounding scale.  When the double route keeps its data, its basis is
// the Arb route's (the same LLL at the same accepted scale), its Gram is
// finite, and its LLL-bound test rejects a sheared basis whose Gram diagonal
// is outside the product bound.  When it falls back, the data is the Arb route's.
struct T2DoubleRowsCounts {
    slong accepted = 0;
    slong fallback = 0;
    // Accepted ideals whose sheared basis the LLL-bound test rejects.
    slong sheared = 0;
};

int check_t2_double_rows_ideal(const silex::Ideal& ideal,
                               T2DoubleRowsCounts& counts) {
    const silex::Order& order = *ideal.parent();
    const slong n = order.degree();
    search::FiniteIdealT2EnumerationData arb_data;
    search::FiniteIdealT2EnumerationData double_data;
    if (!search::build_finite_ideal_t2_enumeration_data_with_retry(
                arb_data, ideal, nullptr, nullptr, search::kMaxFactorAttempts,
                search::kT2EnumerationProofInitialPrecision,
                search::T2EnumerationGram::arb) ||
        arb_data.double_rows ||
        !search::build_finite_ideal_t2_enumeration_data_with_retry(
                double_data, ideal, nullptr, nullptr,
                search::kMaxFactorAttempts,
                search::kT2EnumerationProofInitialPrecision,
                search::T2EnumerationGram::double_rows)) {
        return 1;
    }
    const bool same_basis =
            sflint::fmpz_mat_equal(
                    sflint::FmpzMatConstRef(double_data.basis),
                    sflint::FmpzMatConstRef(arb_data.basis)) &&
            double_data.rounding_scale == arb_data.rounding_scale;
    if (!double_data.double_rows) {
        ++counts.fallback;
        return same_basis &&
                               double_data.quadratic_form_data ==
                                       arb_data.quadratic_form_data &&
                               double_data.initial_bound_value ==
                                       arb_data.initial_bound_value
                       ? 0
                       : 2;
    }
    ++counts.accepted;
    if (!same_basis) {
        return 3;
    }

    // The order rows as the setup reads them: midpoints of the rows at its
    // embedding precision (start precision plus 64 guard bits).
    sflint::FmpzMat identity(n, n);
    sflint::fmpz_mat_one(sflint::FmpzMatRef(identity));
    sflint::ArbMat order_rows(n, n);
    std::vector<double> order_rows_double;
    std::vector<double> gram;
    if (!silex::detail::build_ideal_minkowski_embedding_rows(
                order_rows, sflint::FmpzMatConstRef(identity), order,
                search::kT2EnumerationProofInitialPrecision + 64) ||
        !search::arb_mat_midpoints_double(order_rows_double, order_rows) ||
        !search::reduced_basis_t2_gram_double(
                gram, sflint::FmpzMatConstRef(double_data.basis),
                order_rows_double)) {
        return 4;
    }
    // No fixed Gram tolerance is asserted: the d9 twisted ideal sits at
    // 1.24x of n 2^-50 |r_i| |r_j| from a 300-bit Arb Gram because of
    // cancellation in the order-basis rows, not a defect of the dot
    // products.  The check is that the double Gram is finite and that the
    // selected basis equals the Arb route's (above).
    for (const double entry : gram) {
        if (!std::isfinite(entry)) {
            return 6;
        }
    }

    // One scale of the double route on the accepted basis gives the data
    // the setup kept.
    const double log2_determinant = search::t2_lattice_log2_determinant(ideal);
    std::vector<double> data;
    if (search::t2_double_rows_basis_data(
                data, sflint::FmpzMatConstRef(double_data.basis),
                order_rows_double, log2_determinant) !=
                search::T2DoubleRowsBasis::accepted ||
        data != double_data.quadratic_form_data) {
        return 7;
    }

    // The sheared basis b_0, b_i + 16 b_0 (i >= 1) of the same lattice has
    // the same first vector and the same Gram-Schmidt pivots, whose product
    // is d.  When its Gram diagonal product exceeds 2^(n (n - 1) / 2) d, the
    // product half of the LLL-bound test, which reads the Gram diagonal
    // |b_i|^2 and not the pivots |b_i*|^2, must reject it.  On some ideals
    // the shear's cancellation moves the double pivot product by more than
    // the 2^-20 tolerance, and the data is rejected before that test.
    sflint::FmpzMat sheared(n, n);
    fmpz_mat_set(sheared.raw(), double_data.basis.raw());
    for (slong i = 1; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_addmul_ui(fmpz_mat_entry(sheared.raw(), i, j),
                           fmpz_mat_entry(sheared.raw(), 0, j), 16);
        }
    }
    std::vector<double> sheared_gram;
    if (!search::reduced_basis_t2_gram_double(
                sheared_gram, sflint::FmpzMatConstRef(sheared),
                order_rows_double)) {
        return 8;
    }
    double log2_diagonal_product = 0.0;
    for (slong i = 0; i < n; ++i) {
        log2_diagonal_product += std::log2(
                sheared_gram[static_cast<std::size_t>(i * n + i)]);
    }
    const double nd = static_cast<double>(n);
    if (log2_diagonal_product >
        nd * (nd - 1.0) / 2.0 + log2_determinant + 1.0) {
        const search::T2DoubleRowsBasis status =
                search::t2_double_rows_basis_data(
                        data, sflint::FmpzMatConstRef(sheared),
                        order_rows_double, log2_determinant);
        if (status == search::T2DoubleRowsBasis::outside_lll_bounds) {
            ++counts.sheared;
        } else if (status !=
                   search::T2DoubleRowsBasis::determinant_mismatch) {
            return 9;
        }
    }
    return 0;
}

// Every prime ideal above each rational prime up to 71, and each twisted by
// Q^20 for the first prime Q above 2 (reduced as the honesty search reduces
// its twisted ideals), as the honesty witness search enumerates them.
int check_t2_double_rows_field(const FieldSetup& setup,
                               T2DoubleRowsCounts& counts) {
    const silex::Order& order = setup.maximal_order;
    silex::PrimeIdeal twist_prime(order);
    silex::Ideal twist(order);
    if (!order.is_defined() || !first_prime_above(twist_prime, order, 2) ||
        !twist_prime.get_ideal(twist)) {
        return 1;
    }
    for (const slong p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                          47, 53, 59, 61, 67, 71}) {
        sflint::Fmpz rational_prime;
        sflint::fmpz_set_si(sflint::FmpzRef(rational_prime), p);
        silex::PrimeIdealList primes;
        if (!silex::decompose_prime(primes, order,
                                    sflint::FmpzConstRef(rational_prime))) {
            return 1;
        }
        for (slong k = 0; k < primes.size(); ++k) {
            const silex::PrimeIdeal* prime = primes.at(k);
            silex::Ideal prime_ideal(order);
            silex::Ideal product(order);
            silex::Element back(*order.parent());
            if (prime == nullptr || !prime->get_ideal(prime_ideal) ||
                !product.set(prime_ideal) || !back.one()) {
                return 1;
            }
            for (slong e = 0; e < 20; ++e) {
                silex::Ideal next(order);
                if (!next.multiply(product, twist)) {
                    return 1;
                }
                product.swap(next);
            }
            if (!search::factor_base_honesty_primitive_part(product, back) ||
                !search::factor_base_honesty_reduce_large_ideal(product, back,
                                                                128, nullptr)) {
                return 1;
            }
            const int prime_status =
                    check_t2_double_rows_ideal(prime_ideal, counts);
            if (prime_status != 0) {
                return 10 + prime_status;
            }
            const int twisted_status = check_t2_double_rows_ideal(product, counts);
            if (twisted_status != 0) {
                return 20 + twisted_status;
            }
        }
    }
    return 0;
}

// Corpus fields random_d9_h4_s8, cyclotomic_real_q10_c11_deg10 and
// lmfdb_8_8_73116160000_2.  The double route keeps its data on the first and
// last; on the second, whose maximal-order basis gives reduced bases with
// coordinates above 2^53 and ill-conditioned double rows, it may fall back,
// and then its data is the Arb route's.
int test_t2_double_rows_gram_matches_arb() {
    T2DoubleRowsCounts d9;
    T2DoubleRowsCounts q10;
    T2DoubleRowsCounts l8;
    if (check_t2_double_rows_field(
                maximal_order_setup({2, -2, -3, -2, 4, -2, -4, 1, 2, 1}),
                d9) != 0) {
        return 1;
    }
    if (check_t2_double_rows_field(
                maximal_order_setup({238831, -228566, -128567, 146666, -15113,
                                     -13352, 2746, 334, -97, -2, 1}),
                q10) != 0) {
        return 2;
    }
    if (check_t2_double_rows_field(
                maximal_order_setup({-1, 68, 60, -232, 57, 68, -18, -4, 1}),
                l8) != 0) {
        return 3;
    }
    return d9.accepted > 0 && l8.accepted > 0 && d9.sheared > 0 &&
                   l8.sheared > 0
           ? 0
           : 4;
}

// Corpus field multiquadratic_m23_m31_m47, Q(sqrt(-23), sqrt(-31),
// sqrt(-47)) from a defining polynomial with large coefficients.  Its order
// rows are large (up to about 5e7 on the honesty search's ideals) and the
// reduced bases have large coordinates (up to about 1.5e11), so the double
// rows cancel to entries many orders smaller than their terms and keep
// rounding errors comparable to those entries.  Such a Gram can pass the
// pivot and LLL-bound tests; the check of the pivot product against the
// exact d = |disc| N(I)^2 must send every setup here to the Arb route, and
// the data is then the Arb route's.
int test_t2_double_rows_cancellation_takes_arb_route() {
    T2DoubleRowsCounts counts;
    if (check_t2_double_rows_field(
                maximal_order_setup({38188288, 1822160, -1462612, -62792,
                                     24143, 602, -220, -2, 1}),
                counts) != 0) {
        return 1;
    }
    return counts.fallback > 0 && counts.accepted == 0 ? 0 : 2;
}

// Counts the T2 setup and fallback profile events of the double-row route.
struct T2DoubleRowsEvents {
    slong setups = 0;
    slong setup_fallbacks = 0;
    slong enumeration_retries = 0;
};

void t2_double_rows_profile_callback(void* user,
                                     silex::DiagnosticsModule,
                                     silex::ProfileEvent event,
                                     const char*,
                                     const char* label) noexcept {
    auto* events = static_cast<T2DoubleRowsEvents*>(user);
    if (label == nullptr) {
        return;
    }
    if (event == silex::ProfileEvent::begin_scope &&
        std::strcmp(label, "class_group.honesty_t2_enumeration_data") == 0) {
        ++events->setups;
    } else if (event == silex::ProfileEvent::event &&
               std::strcmp(label,
                           "class_group.finite_enumeration_setup."
                           "double_rows_arb_fallback") == 0) {
        ++events->setup_fallbacks;
    } else if (event == silex::ProfileEvent::event &&
               std::strcmp(label,
                           "class_group.honesty_enumerate_witness."
                           "double_rows_arb_retry") == 0) {
        ++events->enumeration_retries;
    }
}

struct T2DoubleRowsScan {
    bool ok = false;
    bool honest = false;
    search::FactorBaseHonestyScanAudit audit;
    HonestyLogCapture capture;
    slong step_cap_hits = 0;
    slong debug_failures = 0;
    T2DoubleRowsEvents events;
};

// One T2-route honesty scan of `setup` (active and required norm bounds)
// under a T2 double-row seam mode.  With `twists` the sub-factor base is
// the whole base; otherwise there are no twists.  Debug checks run at the
// expensive level, so in a debug-check build every witness the scan accepts
// is also classified by the full factorization.
void run_t2_double_rows_scan(T2DoubleRowsScan& scan,
                             const FieldSetup& setup,
                             ulong active,
                             ulong required,
                             bool twists,
                             search::FactorBaseWitnessEscalation escalation,
                             search::T2DoubleRowsTestSeams::Mode mode) {
    using Hook = silex::detail::TestHookAccess;
    sflint::Fmpz active_bound;
    sflint::Fmpz required_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(active_bound), active);
    sflint::fmpz_set_ui(sflint::FmpzRef(required_bound), required);
    silex::FactorBase base(setup.maximal_order);
    if (!setup.maximal_order.is_defined() || !base.is_defined() ||
        (active > 1 && !base.build_relation_completion_base(
                               sflint::FmpzConstRef(active_bound)))) {
        return;
    }
    silex::detail::SubfactorBaseSchedule schedule;
    for (slong i = 0; i < base.length(); ++i) {
        schedule.subfactor_base.push_back(i);
    }
    schedule.excluded_indices.assign(static_cast<std::size_t>(base.length()),
                                     0);
    silex::DiagnosticsContext diagnostics;
    enable_honesty_log(diagnostics, scan.capture);
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(
                    silex::DiagnosticsModule::class_group),
            t2_double_rows_profile_callback, &scan.events);
    diagnostics.debug_level = silex::DebugLevel::expensive;
    diagnostics.debug_modules = silex::diagnostics_all_modules;
    diagnostics.debug_failure_user = &scan.debug_failures;
    diagnostics.debug_failure_callback = [](
            void* user, silex::DiagnosticsModule, silex::DebugLevel,
            const char*, const char*, int, const char*, const char*) noexcept {
        ++*static_cast<slong*>(user);
    };
    Hook::set_t2_double_rows_mode(mode);
    Hook::factor_base_honesty_element_step_cap_hits_reset();
    scan.honest = false;
    scan.ok = search::factor_base_honesty_check(
            scan.honest, base, sflint::FmpzConstRef(active_bound),
            sflint::FmpzConstRef(required_bound),
            twists && base.length() > 0 ? &schedule : nullptr, UWORD(42),
            true, 128, escalation, &diagnostics, &scan.audit);
    scan.step_cap_hits = Hook::factor_base_honesty_element_step_cap_hits();
    Hook::set_t2_double_rows_mode(search::T2DoubleRowsTestSeams::Mode::off);
}

bool same_t2_double_rows_scan(const T2DoubleRowsScan& left,
                              const T2DoubleRowsScan& right) {
    const search::FactorBaseHonestyScanAudit& a = left.audit;
    const search::FactorBaseHonestyScanAudit& b = right.audit;
    return left.ok == right.ok && left.honest == right.honest &&
           a.rational_prime_checks == b.rational_prime_checks &&
           a.checks_at_or_below_active_bound ==
                   b.checks_at_or_below_active_bound &&
           a.witness_targets == b.witness_targets &&
           a.witnessed_targets == b.witnessed_targets &&
           a.escalated_witnessed_targets == b.escalated_witnessed_targets &&
           a.max_search_stage == b.max_search_stage &&
           a.final_random_state == b.final_random_state &&
           a.t2_norm_prefilter_rejections == b.t2_norm_prefilter_rejections &&
           a.unwitnessed.recorded == b.unwitnessed.recorded &&
           sflint::fmpz_equal(a.unwitnessed.p, b.unwitnessed.p) &&
           a.unwitnessed.stage == b.unwitnessed.stage &&
           a.unwitnessed.factor_attempts == b.unwitnessed.factor_attempts &&
           a.unwitnessed.element_steps == b.unwitnessed.element_steps &&
           left.capture.exhausted_messages ==
                   right.capture.exhausted_messages &&
           left.capture.last_detail == right.capture.last_detail &&
           left.step_cap_hits == right.step_cap_hits;
}

// Seam tests of the double-row fallbacks on whole T2-route honesty scans.
// When every double-row setup sees non-finite rows (fallback to the Arb
// setup) or every enumeration of double-row data is treated as exhausted
// (re-enumeration from the Arb setup), the verdict, audit, random-state
// stream, exhaustion log and step-cap count are those of the Arb-only scan.
// The production scan (double rows kept) reaches the same verdict, and every
// witness it accepts passes the full-factorization classification (checked
// in a debug-check build).  The x^3 - 15 row has an unwitnessable prime (its
// prime above 2 is not principal and the base is empty), so its ideals are
// exhausted at the caps.
int test_t2_double_rows_fallbacks_match_arb_scan() {
    using Mode = search::T2DoubleRowsTestSeams::Mode;
    struct Row {
        std::vector<slong> coefficients;
        ulong active;
        ulong required;
        bool twists;
        search::FactorBaseWitnessEscalation escalation;
    };
    const auto none = search::FactorBaseWitnessEscalation::none;
    const std::vector<Row> rows = {
            {{-250000000002, -1, 1}, 293, 700, true, none},
            {{-6236, 0, 1}, 3, 60, true, none},
            {{-15, 0, 0, 1}, 2, 30, true, none},
            {{-2, 0, 0, 1}, 7, 40, true, none},
            {{1, 0, -10, 0, 1}, 7, 30, true, none},
            {{-15, 0, 0, 1}, 1, 30, false, none},
    };
    slong double_setups = 0;
    slong retries = 0;
    for (const Row& row : rows) {
        const FieldSetup setup = maximal_order_setup(row.coefficients);
        T2DoubleRowsScan expected;
        run_t2_double_rows_scan(expected, setup, row.active, row.required,
                                row.twists, row.escalation, Mode::arb_only);
        if (!expected.ok || expected.debug_failures != 0 ||
            expected.events.setup_fallbacks != 0 ||
            expected.events.enumeration_retries != 0) {
            return 1;
        }
        for (const Mode mode : {Mode::nonfinite_rows, Mode::exhausted}) {
            T2DoubleRowsScan scan;
            run_t2_double_rows_scan(scan, setup, row.active, row.required,
                                    row.twists, row.escalation, mode);
            if (!same_t2_double_rows_scan(scan, expected) ||
                scan.debug_failures != 0) {
                return 2;
            }
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
            if (mode == Mode::nonfinite_rows &&
                (scan.events.setup_fallbacks != scan.events.setups ||
                 scan.events.enumeration_retries != 0)) {
                return 3;
            }
            // Each setup that kept double rows is set up again on the Arb
            // route after its enumeration.
            if (mode == Mode::exhausted &&
                scan.events.setups != 2 * scan.events.enumeration_retries +
                                              scan.events.setup_fallbacks) {
                return 4;
            }
#endif
        }
        T2DoubleRowsScan production;
        run_t2_double_rows_scan(production, setup, row.active, row.required,
                                row.twists, row.escalation, Mode::off);
        if (!production.ok || production.honest != expected.honest ||
            production.debug_failures != 0 ||
            (production.honest && production.audit.witnessed_targets !=
                                          production.audit.witness_targets)) {
            return 5;
        }
        double_setups +=
                production.events.setups - production.events.setup_fallbacks;
        retries += production.events.enumeration_retries;
    }
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The production scans keep double-row data, and the x^3 - 15 rows
    // re-enumerate their exhausted ideals from the Arb setup.
    if (double_setups == 0 || retries == 0) {
        return 6;
    }
#else
    (void) double_setups;
    (void) retries;
#endif
    return 0;
}

int main() {
    assert(test_quadratic_backend_runs_no_honesty_scan() == 0);
    assert(test_element_step_cap_hits_observed() == 0);
    assert(test_scalar_noncontiguous_blocks() == 0);
    assert(test_scalar_witness_reference_cases() == 0);
    assert(test_selected_witness_edge_cases() == 0);
    assert(test_witness_matches_required_prime_as_ideal() == 0);
    assert(test_witness_partial_rational_prime_blocks() == 0);
    assert(test_legacy_scan_predicate_equivalence() == 0);
    assert(test_witness_through_lacking_conjugate() == 0);
    assert(test_t2_double_rows_gram_matches_arb() == 0);
    assert(test_t2_double_rows_cancellation_takes_arb_route() == 0);
    assert(test_t2_double_rows_fallbacks_match_arb_scan() == 0);
    return test_proof_targets_require_complete_decomposition() != 0 ||
                   test_required_prime_predicates_match_full_factorization() !=
                           0 ||
                   test_required_prime_direct_rejects_nonsmooth_cofactor() !=
                           0 ||
                   test_required_prime_norm_prefilter_matches_direct() != 0 ||
                   test_required_prime_norm_prefilter_word_path_matches_fmpz() !=
                           0 ||
                   test_ideal_transforms_preserve_principal_witness() != 0 ||
                   test_truncated_decomposition_withholds_honesty() != 0 ||
                   test_generation_bound_checker_rejects_missing_truncated_target() !=
                           0 ||
                   test_ramified_target_is_never_omitted() != 0 ||
                   test_scan_detects_lower_interval_omissions() !=
                           0 ||
                   test_escalated_witness_search_proves_large_class_group_base() !=
                           0 ||
                   test_bounded_escalation_reports_final_caps() != 0 ||
                   test_escalation_only_without_recovery() != 0 ||
                   test_direct_route_escalation_caps_and_counters() != 0 ||
                   test_routes_pass_intended_escalation_mode() != 0 ||
                   test_terminal_honesty_failure_names_transaction_reason() !=
                           0 ||
                   test_later_attempt_clears_honesty_failure() != 0 ||
                   test_transaction_reuses_passed_honesty_scan() != 0 ||
                   test_passed_honesty_checkpoint_requires_required_bound() !=
                           0 ||
                   test_saturation_proof_without_bounded_generation_is_not_proven() !=
                           0 ||
                   test_failure_detail_json() != 0 ||
                   test_budget_failure_reason_and_detail_json() != 0 ||
                   test_lll_relation_slice_cap_refusal_is_noted() != 0 ||
                   test_direct_witness_search_keeps_lattice_fallbacks() != 0 ||
                   test_honesty_quadratic_kronecker_matches_decomposition() !=
                           0 ||
                   test_t2_witness_prefilter_matches_direct() != 0 ||
                   test_t2_scan_prefilter_matches_full_factorization() != 0
            ? 1
            : 0;
}
