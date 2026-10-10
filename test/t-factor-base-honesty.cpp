#include "class_group/class_group_internal.hpp"
#include "class_group/factor_base_honesty_internal.hpp"
#include "class_group/factor_base_proof_targets_internal.hpp"
#include "class_group/relation_search_driver_internal.hpp"
#include "ideal_factorization/ideal_factorization_internal.hpp"
#include "order/order_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "test_support.hpp"
#include "../tools/class_unit_failure_detail.hpp"

#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>

#include <array>
#include <cassert>
#include <cstring>
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

bool full_factorization_required_prime_result(
        bool& matches,
        const silex::Ideal& ideal,
        const silex::FactorBase& base,
        const silex::PrimeIdeal& required_prime) noexcept {
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
        if (factor.equal(required_prime)) {
            required_exponent += exponent;
        } else if (!base.contains(factor)) {
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

    bool honest = true;
    silex::detail::relation_search::FactorBaseHonestyScanAudit audit;
    const bool check_succeeded =
            silex::detail::relation_search::factor_base_honesty_check(
                    honest, base, sflint::FmpzConstRef(active_bound),
                    sflint::FmpzConstRef(required_bound), nullptr, UWORD(0),
                    true, 128,
                    silex::detail::relation_search::
                            FactorBaseWitnessEscalation::none,
                    nullptr, &audit);

    // The first omitted target is the ramified norm-4 prime above 2.  With no
    // witness schedule the strict check must fail there, independently of the
    // order in which the decomposition producer returned the two primes.
    return check_succeeded || honest || audit.rational_prime_checks != 1 ||
                           audit.checks_at_or_below_active_bound != 1
            ? 1
            : 0;
}

int test_scan_detects_lower_interval_omissions_and_fails_closed() {
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

    bool honest = true;
    silex::detail::relation_search::FactorBaseHonestyScanAudit audit;
    const bool check_succeeded =
            silex::detail::relation_search::factor_base_honesty_check(
                    honest, base, sflint::FmpzConstRef(active_bound),
                    sflint::FmpzConstRef(required_bound), nullptr, UWORD(0),
                    true, 128,
                    silex::detail::relation_search::
                            FactorBaseWitnessEscalation::none,
                    nullptr, &audit);

    // With no sub-factor-base witness state, the first omitted norm-9 ideal
    // fails closed. The audit proves that the scan reached it through p=3.
    if (check_succeeded || honest) {
        return 1;
    }
    return audit.rational_prime_checks != 2 ||
                   audit.checks_at_or_below_active_bound != 2
            ? 1
            : 0;
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
            sflint::FmpzConstRef(active_bound), false, false, nullptr);
    silex::detail::ClassUnitTransactionAccess::set_run_context(context,
                                                               nullptr);
    return decision;
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

}  // namespace

int main() {
    assert(test_scalar_noncontiguous_blocks() == 0);
    assert(test_scalar_witness_reference_cases() == 0);
    assert(test_selected_witness_edge_cases() == 0);
    assert(test_witness_partial_rational_prime_blocks() == 0);
    assert(test_legacy_scan_predicate_equivalence() == 0);
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
                   test_scan_detects_lower_interval_omissions_and_fails_closed() !=
                           0 ||
                   test_escalated_witness_search_proves_large_class_group_base() !=
                           0 ||
                   test_bounded_escalation_reports_final_caps() != 0 ||
                   test_escalation_only_without_recovery() != 0 ||
                   test_terminal_honesty_failure_names_transaction_reason() !=
                           0 ||
                   test_later_attempt_clears_honesty_failure() != 0 ||
                   test_transaction_reuses_passed_honesty_scan() != 0 ||
                   test_failure_detail_json() != 0
            ? 1
            : 0;
}
