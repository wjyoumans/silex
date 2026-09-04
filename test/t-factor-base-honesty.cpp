#include "class_group/factor_base_honesty_internal.hpp"
#include "class_group/factor_base_proof_targets_internal.hpp"
#include "ideal_factorization/ideal_factorization_internal.hpp"
#include "test_support.hpp"

#include <silex/class_group.hpp>
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

#include <array>
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

    // Keep the conservative singleton behavior on a complete inert
    // decomposition: even an unramified principal ideal remains explicit.
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
                    false, 128, nullptr, &audit);

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
                    true, 128, nullptr, &audit);

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
                    true, 128, nullptr, &audit);

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

}  // namespace

int main() {
    return test_proof_targets_require_complete_decomposition() != 0 ||
                   test_required_prime_predicates_match_full_factorization() !=
                           0 ||
                   test_ideal_transforms_preserve_principal_witness() != 0 ||
                   test_truncated_decomposition_withholds_honesty() != 0 ||
                   test_generation_bound_checker_rejects_missing_truncated_target() !=
                           0 ||
                   test_ramified_target_is_never_omitted() != 0 ||
                   test_scan_detects_lower_interval_omissions_and_fails_closed() !=
                           0
            ? 1
            : 0;
}
