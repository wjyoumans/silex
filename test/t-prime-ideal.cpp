#include <flint/fmpq_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/nmod_mat.hpp>
#include <silex/factored_element.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_ring.hpp>

#include "ideal/prime_inverse_internal.hpp"
#include "order/order_internal.hpp"
#include "prime_ideal/prime_ideal_internal.hpp"
#include "residue_field/residue_field_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"

#include <array>
#include <cassert>
#include <cstring>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;
using PrimeIdealAccess = silex::detail::PrimeIdealAccess;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

void poly_x3_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

void poly_cubic_disc1724(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -8);
}

void poly_quartic_disc835584(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -8);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -3);
}

void poly_quartic_disc41184(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 6);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 4);
}

void poly_quintic_disc48755152(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 4);
}

void poly_x4_minus_x_minus_1(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
}

void poly_index3_quartic(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -5);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -9);
}

void poly_degree6_benchmark_field(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 6, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

silex::Order order_by_polynomial(silex::NumberField& field,
                                  sflint::FmpqPoly& polynomial) noexcept {
    field = field_by_polynomial(polynomial);
    return silex::test::equation_order(field);
}

bool fmpz_poly_coeff_is_si(const sflint::FmpzPoly& polynomial,
                           slong index,
                           slong value) noexcept {
    sflint::Fmpz coeff;
    fmpz_poly_get_coeff_fmpz(coeff.raw(), polynomial.raw(), index);
    return sflint::fmpz_equal_si(coeff, value);
}

bool fmpz_mat_entry_is_si(const sflint::FmpzMat& matrix,
                          slong row,
                          slong col,
                          slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(matrix), row, col),
            value);
}

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

void qmat_entry_si(sflint::FmpqMat& matrix,
                   slong row,
                   slong col,
                   slong value) noexcept {
    fmpq_set_si(fmpq_mat_entry(matrix.raw(), row, col), value, 1);
}

void qmat_entry_frac_si(sflint::FmpqMat& matrix,
                        slong row,
                        slong col,
                        slong num,
                        slong den) noexcept {
    fmpq_set_si(fmpq_mat_entry(matrix.raw(), row, col), num, den);
}

bool set_rational_principal(silex::Ideal& ideal, slong value) noexcept {
    const silex::Order* order = ideal.parent();
    if (order == nullptr) {
        return false;
    }

    silex::OrderElement generator(*order);
    return generator.set_si(value) && ideal.set_principal(generator);
}

slong count_with_ef(const silex::PrimeIdealList& primes,
                    slong ramification_index,
                    slong residue_degree) noexcept {
    slong count = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime != nullptr &&
            prime->ramification_index() == ramification_index &&
            prime->residue_degree() == residue_degree) {
            ++count;
        }
    }
    return count;
}

bool check_ef_sum(const silex::PrimeIdealList& primes, slong degree) noexcept {
    slong sum = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr) {
            return false;
        }
        sum += prime->ramification_index() * prime->residue_degree();
    }
    return sum == degree;
}

bool check_norms(const silex::PrimeIdealList& primes) noexcept {
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr || prime->parent() == nullptr) {
            return false;
        }

        sflint::Fmpz prime_norm;
        sflint::Fmpz ideal_norm;
        silex::Ideal ideal(*prime->parent());
        if (!prime->norm(sflint::FmpzRef(prime_norm)) ||
            !prime->get_ideal(ideal) ||
            !ideal.norm(sflint::FmpzRef(ideal_norm)) ||
            !sflint::fmpz_equal(prime_norm, ideal_norm)) {
            return false;
        }
    }
    return true;
}

bool first_prime_has_data(const silex::PrimeIdealList& primes,
                          slong list_size,
                          slong rational_prime,
                          slong ramification_index,
                          slong residue_degree,
                          slong norm) noexcept {
    if (!primes.is_defined() || primes.size() != list_size) {
        return false;
    }

    const silex::PrimeIdeal* prime = primes.at(0);
    if (prime == nullptr ||
        prime->ramification_index() != ramification_index ||
        prime->residue_degree() != residue_degree) {
        return false;
    }

    sflint::Fmpz stored_p;
    sflint::Fmpz stored_norm;
    return prime->rational_prime(sflint::FmpzRef(stored_p)) &&
           sflint::fmpz_equal_si(stored_p, rational_prime) &&
           prime->norm(sflint::FmpzRef(stored_norm)) &&
           sflint::fmpz_equal_si(stored_norm, norm);
}

bool degree_one_root_from_residue_polynomial(
        sflint::Fmpz& root,
        const silex::PrimeIdeal& prime) noexcept {
    sflint::Fmpz p;
    sflint::FmpzPoly residue;
    sflint::Fmpz constant;
    sflint::Fmpz leading;
    sflint::Fmpz inverse_leading;
    if (!prime.rational_prime(sflint::FmpzRef(p)) ||
        !prime.residue_polynomial(sflint::FmpzPolyRef(residue)) ||
        fmpz_poly_degree(residue.raw()) != 1) {
        return false;
    }

    fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
    fmpz_poly_get_coeff_fmpz(leading.raw(), residue.raw(), 1);
    fmpz_mod(leading.raw(), leading.raw(), p.raw());
    if (fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) == 0) {
        return false;
    }

    fmpz_neg(root.raw(), constant.raw());
    fmpz_mod(root.raw(), root.raw(), p.raw());
    fmpz_mul(root.raw(), root.raw(), inverse_leading.raw());
    fmpz_mod(root.raw(), root.raw(), p.raw());
    return true;
}

int test_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(prime->ramification_index() == 1);
    assert(prime->residue_degree() == 1);

    sflint::Fmpz stored_p;
    sflint::Fmpz norm;
    assert(prime->rational_prime(sflint::FmpzRef(stored_p)));
    assert(sflint::fmpz_equal_si(stored_p, 7));
    assert(prime->norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 7));

    sflint::FmpzPoly residue_polynomial;
    assert(prime->residue_polynomial(sflint::FmpzPolyRef(
            residue_polynomial)));
    assert(fmpz_poly_degree(residue_polynomial.raw()) == 1);
    assert(fmpz_poly_coeff_is_si(residue_polynomial, 1, 1));

    silex::Element sixteen(field);
    sflint::FmpzPoly reduced;
    assert(sixteen.set_si(16));
    assert(prime->reduce(sflint::FmpzPolyRef(reduced), sixteen));
    assert(fmpz_poly_degree(reduced.raw()) == 0);
    assert(fmpz_poly_coeff_is_si(reduced, 0, 2));

    silex::PrimeIdeal copy(order);
    assert(copy.set(*prime));
    assert(copy.equal(*prime));
    silex::PrimeIdeal moved(std::move(copy));
    assert(moved.equal(*prime));

    return 0;
}

// The internal degree-one constructor must not trust its root: on
// x^2 - 2, (7, theta - 1) is the unit ideal because 1 - 2 = -1 is not 0 mod 7,
// and at p = 2 the reduction x^2 is not squarefree (2 is ramified), so
// (2, theta) is not a degree-one prime with e = 1.  Both fail and preserve the
// output; the true root 3 (9 - 2 = 7) still builds the prime of norm 7.
int test_degree_one_prime_from_root_rejects_non_roots() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    sflint::Fmpz root;
    silex::PrimeIdeal prime(order);

    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(root, 1));
    assert(!PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
            prime, order, sflint::FmpzConstRef(p),
            sflint::FmpzConstRef(root)));
    assert(!prime.has_prime_data());

    assert(set_fmpz_si(p, 2));
    assert(set_fmpz_si(root, 0));
    assert(!PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
            prime, order, sflint::FmpzConstRef(p),
            sflint::FmpzConstRef(root)));
    assert(!prime.has_prime_data());

    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(root, 3));
    assert(PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
            prime, order, sflint::FmpzConstRef(p),
            sflint::FmpzConstRef(root)));
    sflint::Fmpz norm;
    assert(prime.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 7));
    assert(prime.residue_degree() == 1);

    return 0;
}

// On a maximal order with [O : Z[theta]] = 2, the root 1 of
// f = x^3 - 3x^2 + 2x - 8 mod 2 does not give a Dedekind-Kummer prime
// (2 divides the index; f = x^2 (x + 1) mod 2 is not squarefree), so the
// internal constructor refuses it and leaves its output unchanged.  The
// squarefree condition is sufficient but not necessary: 1 is a simple root
// mod 2 (2 does not divide f'(1)) and that alone makes (2, theta - 1)O a
// prime with e = f = 1 (Hensel lifts 1 to a root of f in Z_2, so f has a
// linear factor over Q_2).  The constructor refuses it all the same.
int test_degree_one_prime_from_root_rejects_index_prime() {
    sflint::FmpqPoly polynomial;
    poly_cubic_disc1724(polynomial);
    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Fmpz p;
    sflint::Fmpz root;
    assert(set_fmpz_si(p, 2));
    assert(set_fmpz_si(root, 1));
    silex::PrimeIdeal prime(maximal);
    assert(!PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
            prime, maximal, sflint::FmpzConstRef(p),
            sflint::FmpzConstRef(root)));
    assert(!prime.has_prime_data());

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 3);
    for (slong i = 0; i < 3; ++i) {
        assert(primes.at(i)->ramification_index() == 1);
        assert(primes.at(i)->residue_degree() == 1);
    }

    silex::Element theta(field);
    assert(theta.gen());
    silex::OrderElement theta_minus_one(maximal);
    silex::OrderElement one(maximal);
    silex::OrderElement two(maximal);
    assert(theta_minus_one.set_element(theta));
    assert(one.set_si(1));
    assert(theta_minus_one.subtract(theta_minus_one, one));
    assert(two.set_si(2));
    silex::Ideal two_ideal(maximal);
    silex::Ideal theta_ideal(maximal);
    silex::Ideal generated(maximal);
    assert(two_ideal.set_principal(two));
    assert(theta_ideal.set_principal(theta_minus_one));
    assert(generated.add(two_ideal, theta_ideal));
    sflint::Fmpz norm;
    assert(generated.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 2));
    slong valuation_one = 0;
    for (slong i = 0; i < 3; ++i) {
        slong value = -1;
        assert(primes.at(i)->valuation(value, generated));
        assert(value == 0 || value == 1);
        valuation_one += value;
    }
    assert(valuation_one == 1);
    return 0;
}

// On the cubic_disc1724 maximal order ([O : Z[alpha]] = 2) the alpha-
// convention reduction needs p not dividing the index.  At p = 3 the inert
// prime reduces alpha.  No construction route records alpha residue data at
// p = 2, so the test forces such a record through the internal hook and
// checks that reduce() and residue_variable_numerator() refuse it.  Only
// reduce() is checked to leave its outputs unchanged; the outputs of
// residue_variable_numerator() are unspecified on failure.
int test_alpha_reduction_refuses_index_prime() {
    sflint::FmpqPoly polynomial;
    poly_cubic_disc1724(polynomial);
    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));

    silex::Element alpha(field);
    assert(alpha.gen());

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 3));
    silex::PrimeIdealList primes3;
    assert(silex::decompose_prime(primes3, maximal, sflint::FmpzConstRef(p)));
    assert(primes3.size() == 1);
    assert(!PrimeIdealAccess::residue_uses_integral_generator(*primes3.at(0)));
    assert(PrimeIdealAccess::order_basis_is_p_integral(*primes3.at(0)));
    sflint::FmpzPoly reduced;
    assert(primes3.at(0)->reduce(sflint::FmpzPolyRef(reduced), alpha));

    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes2;
    assert(silex::decompose_prime(primes2, maximal, sflint::FmpzConstRef(p)));
    assert(primes2.size() > 0);
    silex::PrimeIdeal prime(maximal);
    assert(prime.set(*primes2.at(0)));
    sflint::FmpzPoly residue_polynomial;
    assert(!prime.residue_polynomial(
            sflint::FmpzPolyRef(residue_polynomial)));

    // Record g = x in the alpha convention at p = 2.
    fmpz_poly_set_coeff_si(residue_polynomial.raw(), 1, 1);
    assert(silex::detail::TestHookAccess::set_alpha_residue_polynomial(
            prime, sflint::FmpzPolyConstRef(residue_polynomial)));
    assert(!PrimeIdealAccess::residue_uses_integral_generator(prime));
    assert(!PrimeIdealAccess::order_basis_is_p_integral(prime));

    sflint::FmpzPoly untouched;
    fmpz_poly_set_si(untouched.raw(), 77);
    assert(!prime.reduce(sflint::FmpzPolyRef(untouched), alpha));
    assert(fmpz_poly_degree(untouched.raw()) == 0);
    assert(fmpz_poly_coeff_is_si(untouched, 0, 77));

    sflint::FmpzPoly numerator;
    sflint::Fmpz denominator;
    assert(!silex::detail::residue_variable_numerator(numerator, denominator,
                                                      prime, alpha));
    return 0;
}

int test_degree_one_prime_from_root_matches_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_x4_minus_x_minus_1(polynomial);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p), 1));
    assert(primes.size() > 0);

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        assert(prime->residue_degree() == 1);

        sflint::Fmpz root;
        assert(degree_one_root_from_residue_polynomial(root, *prime));
        silex::PrimeIdeal reconstructed(order);
        assert(PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
                reconstructed, order, sflint::FmpzConstRef(p),
                sflint::FmpzConstRef(root)));
        assert(reconstructed.equal(*prime));
        ++checked;
    }

    assert(checked == primes.size());
    return 0;
}

int test_invalid_prime_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(first_prime_has_data(primes, 2, 11, 1, 1, 11));

    assert(set_fmpz_si(p, 1));
    assert(!silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(first_prime_has_data(primes, 2, 11, 1, 1, 11));

    assert(set_fmpz_si(p, -3));
    assert(!silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(first_prime_has_data(primes, 2, 11, 1, 1, 11));

    return 0;
}

int test_unsupported_order_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order equation_order;
    equation_order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, equation_order,
                                  sflint::FmpzConstRef(p)));
    assert(first_prime_has_data(primes, 2, 11, 1, 1, 11));

    silex::Order explicit_order(field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 3, 2);
    qmat_entry_frac_si(basis, 1, 1, 3, 2);
    assert(explicit_order.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(!explicit_order.is_equation_order());
    assert(!explicit_order.is_maximal());

    assert(set_fmpz_si(p, 3));
    assert(!silex::decompose_prime(primes, explicit_order,
                                   sflint::FmpzConstRef(p)));
    assert(first_prime_has_data(primes, 2, 11, 1, 1, 11));

    assert(set_fmpz_si(p, 1));
    assert(!silex::decompose_prime(primes, explicit_order,
                                   sflint::FmpzConstRef(p)));
    assert(first_prime_has_data(primes, 2, 11, 1, 1, 11));

    return 0;
}

int test_quadratic_splitting_types() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p3;
    assert(set_fmpz_si(p3, 3));
    silex::PrimeIdealList primes3;
    assert(silex::decompose_prime(primes3, order, sflint::FmpzConstRef(p3)));
    assert(primes3.size() == 1);
    const silex::PrimeIdeal* prime3 = primes3.at(0);
    assert(prime3 != nullptr);
    assert(prime3->ramification_index() == 1);
    assert(prime3->residue_degree() == 2);

    silex::Element theta(field);
    sflint::FmpzPoly reduced_theta;
    assert(theta.gen());
    assert(prime3->reduce(sflint::FmpzPolyRef(reduced_theta), theta));
    assert(fmpz_poly_degree(reduced_theta.raw()) == 1);
    assert(fmpz_poly_coeff_is_si(reduced_theta, 1, 1));

    sflint::Fmpz p11;
    assert(set_fmpz_si(p11, 11));
    silex::PrimeIdealList primes11;
    assert(silex::decompose_prime(primes11, order, sflint::FmpzConstRef(p11)));
    assert(primes11.size() == 2);

    for (slong i = 0; i < primes11.size(); ++i) {
        const silex::PrimeIdeal* split_prime = primes11.at(i);
        assert(split_prime != nullptr);
        sflint::Fmpz norm;
        sflint::FmpzPoly root;
        assert(split_prime->ramification_index() == 1);
        assert(split_prime->residue_degree() == 1);
        assert(split_prime->norm(sflint::FmpzRef(norm)));
        assert(sflint::fmpz_equal_si(norm, 11));
        assert(split_prime->reduce(sflint::FmpzPolyRef(root), theta));
        assert(fmpz_poly_degree(root.raw()) == 0);

        sflint::Fmpz coeff;
        fmpz_poly_get_coeff_fmpz(coeff.raw(), root.raw(), 0);
        const ulong r = fmpz_fdiv_ui(coeff.raw(), 11);
        assert((r * r) % 11 == 5);
    }

    return 0;
}

int test_explicit_maximal_order_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order equation_order;
    equation_order = order_by_polynomial(field, polynomial);

    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());
    assert(!maximal_order.is_equation_order());

    silex::Element theta(field);
    assert(theta.gen());

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes2;
    assert(silex::decompose_prime(primes2, maximal_order,
                                  sflint::FmpzConstRef(p)));
    assert(primes2.size() == 1);
    assert(count_with_ef(primes2, 1, 2) == 1);
    assert(check_ef_sum(primes2, 2));
    assert(check_norms(primes2));
    sflint::FmpzPoly residue_polynomial;
    sflint::FmpzPoly reduced;
    assert(primes2.at(0)->residue_polynomial(
            sflint::FmpzPolyRef(residue_polynomial)));
    assert(primes2.at(0)->reduce(sflint::FmpzPolyRef(reduced), theta));
    assert(fmpz_poly_degree(reduced.raw()) == 0);
    assert(fmpz_poly_coeff_is_si(reduced, 0, 1));

    assert(set_fmpz_si(p, 5));
    silex::PrimeIdealList primes5;
    assert(silex::decompose_prime(primes5, maximal_order,
                                  sflint::FmpzConstRef(p)));
    assert(primes5.size() == 1);
    assert(count_with_ef(primes5, 2, 1) == 1);
    assert(check_ef_sum(primes5, 2));
    assert(check_norms(primes5));

    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes11;
    assert(silex::decompose_prime(primes11, maximal_order,
                                  sflint::FmpzConstRef(p)));
    assert(primes11.size() == 2);
    assert(count_with_ef(primes11, 1, 1) == 2);
    assert(check_ef_sum(primes11, 2));
    assert(check_norms(primes11));

    silex::PrimeIdealList degree_one11;
    assert(silex::decompose_prime(degree_one11, maximal_order,
                                  sflint::FmpzConstRef(p), 1));
    assert(degree_one11.size() == 2);
    assert(count_with_ef(degree_one11, 1, 1) == 2);
    assert(check_ef_sum(degree_one11, 2));
    assert(check_norms(degree_one11));
    assert(degree_one11.at(0)->residue_polynomial(
            sflint::FmpzPolyRef(residue_polynomial)));
    assert(degree_one11.at(0)->reduce(sflint::FmpzPolyRef(reduced), theta));

    for (slong i = 0; i < primes11.size(); ++i) {
        const silex::PrimeIdeal* split_prime = primes11.at(i);
        assert(split_prime != nullptr);

        sflint::Fmpz omega_root;
        sflint::Fmpz theta_root;
        assert(degree_one_root_from_residue_polynomial(
                omega_root, *split_prime));
        fmpz_mul_ui(theta_root.raw(), omega_root.raw(), 2);
        fmpz_sub_ui(theta_root.raw(), theta_root.raw(), 1);
        fmpz_mod(theta_root.raw(), theta_root.raw(), p.raw());

        silex::PrimeIdeal reconstructed(maximal_order);
        assert(PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
                reconstructed, maximal_order, sflint::FmpzConstRef(p),
                sflint::FmpzConstRef(theta_root)));
        assert(reconstructed.equal(*split_prime));
    }

    // reference `idealprimedec_kummer` expresses T/u in the maximal-order basis
    // before caching its multiplication matrix.  Exercise that exact basis
    // conversion on both split primes and on repeated valuation calls.
    for (slong i = 0; i < primes11.size(); ++i) {
        const silex::PrimeIdeal* split_prime = primes11.at(i);
        assert(split_prime != nullptr);
        assert(split_prime->reduce(sflint::FmpzPolyRef(reduced), theta));
        assert(fmpz_poly_degree(reduced.raw()) == 0);

        sflint::Fmpz root;
        fmpz_poly_get_coeff_fmpz(root.raw(), reduced.raw(), 0);
        silex::Element root_element(field);
        silex::Element difference(field);
        assert(root_element.set_si(
                static_cast<slong>(fmpz_fdiv_ui(root.raw(), 11))));
        assert(difference.subtract(theta, root_element));

        silex::OrderElement alpha(maximal_order);
        silex::OrderElement alpha_squared(maximal_order);
        assert(alpha.set_element(difference));
        assert(alpha_squared.multiply(alpha, alpha));

        slong value = -1;
        assert(split_prime->valuation(value, alpha_squared));
        assert(value == 2);
        value = -1;
        assert(split_prime->valuation(value, alpha_squared));
        assert(value == 2);
    }

    return 0;
}

int test_split_index_prime_uses_integral_quadratic_generator() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 17);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(
            primes, maximal, sflint::FmpzConstRef(p), 1));
    assert(primes.size() == 2);

    silex::Element theta(field);
    assert(theta.gen());
    bool saw_root_zero = false;
    bool saw_root_one = false;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);

        sflint::Fmpz omega_root;
        sflint::FmpzPoly reduced_theta;
        assert(degree_one_root_from_residue_polynomial(
                omega_root, *prime));
        const ulong root = fmpz_fdiv_ui(omega_root.raw(), 2);
        saw_root_zero = saw_root_zero || root == 0;
        saw_root_one = saw_root_one || root == 1;
        assert(prime->reduce(sflint::FmpzPolyRef(reduced_theta), theta));
        assert(fmpz_poly_degree(reduced_theta.raw()) == 0);
        assert(fmpz_poly_coeff_is_si(reduced_theta, 0, 1));
    }
    assert(saw_root_zero && saw_root_one);

    sflint::Fmpz theta_root;
    assert(set_fmpz_si(theta_root, 1));
    silex::PrimeIdeal from_theta_root(maximal);
    assert(!PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
            from_theta_root, maximal, sflint::FmpzConstRef(p),
            sflint::FmpzConstRef(theta_root)));

    return 0;
}

int test_index3_quartic_maximal_nonindex_prime_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_index3_quartic(polynomial);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    assert(!maximal.is_equation_order());

    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 3));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    assert(count_with_ef(primes, 1, 4) == 1);
    assert(check_ef_sum(primes, 4));
    assert(check_norms(primes));

    sflint::Fmpz norm;
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(prime->norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 16));

    sflint::FmpzPoly residue_polynomial;
    assert(prime->residue_polynomial(sflint::FmpzPolyRef(
            residue_polynomial)));
    assert(fmpz_poly_degree(residue_polynomial.raw()) == 4);
    for (slong i = 0; i <= 4; ++i) {
        assert(fmpz_poly_coeff_is_si(residue_polynomial, i, 1));
    }

    assert(set_fmpz_si(p, 3));
    silex::PrimeIdealList primes3;
    assert(silex::decompose_prime(primes3, maximal, sflint::FmpzConstRef(p)));
    assert(primes3.size() == 2);
    assert(count_with_ef(primes3, 1, 2) == 2);
    assert(check_ef_sum(primes3, 4));
    assert(check_norms(primes3));

    for (slong i = 0; i < primes3.size(); ++i) {
        const silex::PrimeIdeal* prime3 = primes3.at(i);
        assert(prime3 != nullptr);
        assert(prime3->norm(sflint::FmpzRef(norm)));
        assert(sflint::fmpz_equal_si(norm, 9));
        assert(!prime3->residue_polynomial(
                sflint::FmpzPolyRef(residue_polynomial)));
    }

    return 0;
}

int test_index2_cubic_unramified_prime_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_cubic_disc1724(polynomial);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    assert(!maximal.is_equation_order());

    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 3);
    assert(count_with_ef(primes, 1, 1) == 3);
    assert(check_ef_sum(primes, 3));
    assert(check_norms(primes));

    return 0;
}

int test_index8_quartic_ramified_prime_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_quartic_disc835584(polynomial);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    assert(!maximal.is_equation_order());

    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 8));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    assert(count_with_ef(primes, 4, 1) == 1);
    assert(check_ef_sum(primes, 4));
    assert(check_norms(primes));

    return 0;
}

int test_index2_quartic_mixed_ramified_prime_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_quartic_disc41184(polynomial);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    assert(!maximal.is_equation_order());

    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);
    assert(count_with_ef(primes, 2, 1) == 1);
    assert(count_with_ef(primes, 1, 2) == 1);
    assert(check_ef_sum(primes, 4));
    assert(check_norms(primes));

    return 0;
}

int test_index2_quintic_mixed_ramified_prime_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_quintic_disc48755152(polynomial);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    assert(!maximal.is_equation_order());

    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 3);
    assert(count_with_ef(primes, 1, 1) == 2);
    assert(count_with_ef(primes, 3, 1) == 1);
    assert(check_ef_sum(primes, 5));
    assert(check_norms(primes));

    return 0;
}

void poly_from_coefficients(sflint::FmpqPoly& polynomial,
                            const std::vector<slong>& low_to_high) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t i = 0; i < low_to_high.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       low_to_high[i]);
    }
}

// Decomposes p in the maximal order of the field and checks the (e, f)
// multiset, distinctness, norms and prod P^e = pO, and that a second call
// returns the same ideals in the same order.
int check_index_divisor_decomposition(const std::vector<slong>& low_to_high,
                                      slong rational_prime,
                                      slong count,
                                      slong ramification_index,
                                      slong residue_degree) {
    sflint::FmpqPoly polynomial;
    poly_from_coefficients(polynomial, low_to_high);
    const slong degree = static_cast<slong>(low_to_high.size()) - 1;

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(!maximal.is_equation_order());

    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(fmpz_divisible_si(index.raw(), rational_prime) != 0);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, rational_prime));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == count);
    assert(count_with_ef(primes, ramification_index, residue_degree) ==
           count);
    assert(check_ef_sum(primes, degree));
    assert(check_norms(primes));

    std::vector<silex::Ideal> ideals;
    for (slong i = 0; i < primes.size(); ++i) {
        ideals.emplace_back(maximal);
        assert(primes.at(i)->get_ideal(ideals.back()));
        for (slong j = 0; j < i; ++j) {
            assert(!ideals[static_cast<std::size_t>(j)].equal(ideals.back()));
        }
    }

    silex::Ideal product(maximal);
    assert(set_rational_principal(product, 1));
    for (const silex::Ideal& ideal : ideals) {
        for (slong k = 0; k < ramification_index; ++k) {
            silex::Ideal next(maximal);
            assert(next.multiply(product, ideal));
            product = std::move(next);
        }
    }
    silex::Ideal p_ideal(maximal);
    assert(set_rational_principal(p_ideal, rational_prime));
    assert(product.equal(p_ideal));

    silex::PrimeIdealList again;
    assert(silex::decompose_prime(again, maximal, sflint::FmpzConstRef(p)));
    assert(again.size() == count);
    for (slong i = 0; i < again.size(); ++i) {
        silex::Ideal ideal(maximal);
        assert(again.at(i)->get_ideal(ideal));
        assert(ideal.equal(ideals[static_cast<std::size_t>(i)]));
    }

    return 0;
}

// p^(number of primes above p) exceeds the bounded enumeration in these
// fields, so the recursive etale split decomposes p.  The (e, f) data match
// reference 2.17.3 idealprimedec.
int test_index_divisor_recursive_split_unramified_degree_one() {
    // LMFDB 6.6.46411625.1, p = 47: six primes with e = f = 1.
    return check_index_divisor_decomposition(
            {-1, -16, 58, 17, -18, -1, 1}, 47, 6, 1, 1);
}

int test_index_divisor_recursive_split_unramified_degree_two() {
    // Q(sqrt(-3), sqrt(-5), sqrt(-7)), p = 11: four primes with f = 2.
    return check_index_divisor_decomposition(
            {1201, 598, -524, -146, 69, 2, 4, -4, 1}, 11, 4, 1, 2);
}

int test_index_divisor_recursive_split_ramified() {
    // Q(sqrt(5), sqrt(-7), sqrt(-11)), p = 11: four primes with e = 2.
    return check_index_divisor_decomposition(
            {144, -300, 259, -196, 162, -48, 24, -2, 1}, 11, 4, 2, 1);
}

// Decomposes p in the maximal order and checks each prime's e and full HNF,
// in order.  `hnfs` holds one row-major HNF per prime.
int check_index_divisor_decomposition_hnfs(
        const std::vector<slong>& low_to_high,
        slong rational_prime,
        const std::vector<slong>& ramification_indices,
        const std::vector<std::vector<slong>>& hnfs) {
    sflint::FmpqPoly polynomial;
    poly_from_coefficients(polynomial, low_to_high);
    const slong degree = static_cast<slong>(low_to_high.size()) - 1;

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(!maximal.is_equation_order());

    sflint::Fmpz p;
    assert(set_fmpz_si(p, rational_prime));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == static_cast<slong>(hnfs.size()));
    for (slong i = 0; i < primes.size(); ++i) {
        const std::vector<slong>& expected = hnfs[static_cast<std::size_t>(i)];
        assert(primes.at(i)->ramification_index() ==
               ramification_indices[static_cast<std::size_t>(i)]);
        assert(primes.at(i)->residue_degree() == 1);
        silex::Ideal ideal(maximal);
        assert(primes.at(i)->get_ideal(ideal));
        sflint::FmpzMat hnf(degree, degree);
        assert(ideal.get_hnf(sflint::FmpzMatRef(hnf)));
        for (slong row = 0; row < degree; ++row) {
            for (slong col = 0; col < degree; ++col) {
                assert(fmpz_mat_entry_is_si(
                        hnf, row, col,
                        expected[static_cast<std::size_t>(row * degree +
                                                          col)]));
            }
        }
    }
    return 0;
}

// p^(number of primes) exceeds the bounded enumeration here too, but the
// maximal order has a monogenic generator, so the Dedekind-Kummer path
// decomposes p as before the recursive split was added and keeps its prime
// order.  17^3 g(x/17) for g = y^3 - 6y^2 + 11y + 11, where 17 splits
// completely.
int test_index_divisor_monogenic_order_unchanged_unramified() {
    return check_index_divisor_decomposition_hnfs(
            {54043, 3179, -102, 1}, 17, {1, 1, 1},
            {{1, 0, 15, 0, 1, 11, 0, 0, 17},
             {1, 0, 16, 0, 1, 16, 0, 0, 17},
             {1, 0, 4, 0, 1, 8, 0, 0, 17}});
}

// Pins the prime order the separating fixed element gives.  The maximal
// order of Q[x]/(x^3 - 5x^2 - 2x - 3) (field discriminant -239) has 3 as an
// index divisor with three unramified degree-one primes: nullity 3 <= p = 3
// and p^nullity = 27 <= 4096, so find_separating_fixed_element splits it and
// the prime order is the order of the factors of the minimal polynomial of
// the first accepted candidate.  A change that skips the search or reorders
// the primes changes the list below.  The expected HNFs were taken from the
// current decomposition.  An independent computation of the splitting type
// (three primes with e = f = 1) agrees, but it uses another basis and prime
// order, so it cannot confirm the matrices themselves.
int test_index_divisor_separating_element_order() {
    return check_index_divisor_decomposition_hnfs(
            {-3, -2, -5, 1}, 3, {1, 1, 1},
            {{1, 0, 2, 0, 1, 0, 0, 0, 3},
             {1, 2, 0, 0, 3, 0, 0, 0, 1},
             {1, 0, 2, 0, 1, 1, 0, 0, 3}});
}

// 17^4 g(x/17) for g = y (y - 1) (y - 2)^2 + 17: 17 = P1 P2 P3^2.
int test_index_divisor_monogenic_order_unchanged_ramified() {
    return check_index_divisor_decomposition_hnfs(
            {1419857, -19652, 2312, -85, 1}, 17, {1, 1, 2},
            {{17, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
             {1, 0, 0, 16, 0, 1, 0, 16, 0, 0, 1, 16, 0, 0, 0, 17},
             {1, 0, 0, 2, 0, 1, 0, 4, 0, 0, 1, 8, 0, 0, 0, 17}});
}

// Compares two decompositions of one prime over two orders with the same
// basis, prime by prime: e, f, the HNF, the known principal generator, and
// the residue polynomial and Kummer generator when present.
bool same_decomposition(const silex::PrimeIdealList& left,
                        const silex::Order& left_order,
                        const silex::PrimeIdealList& right,
                        const silex::Order& right_order) {
    if (left.size() != right.size()) {
        return false;
    }
    const slong n = left_order.degree();
    for (slong i = 0; i < left.size(); ++i) {
        const silex::PrimeIdeal& a = *left.at(i);
        const silex::PrimeIdeal& b = *right.at(i);
        if (a.ramification_index() != b.ramification_index() ||
            a.residue_degree() != b.residue_degree()) {
            return false;
        }
        silex::Ideal ideal_a(left_order);
        silex::Ideal ideal_b(right_order);
        sflint::FmpzMat hnf_a(n, n);
        sflint::FmpzMat hnf_b(n, n);
        if (!a.get_ideal(ideal_a) || !b.get_ideal(ideal_b) ||
            !ideal_a.get_hnf(sflint::FmpzMatRef(hnf_a)) ||
            !ideal_b.get_hnf(sflint::FmpzMatRef(hnf_b)) ||
            fmpz_mat_equal(hnf_a.raw(), hnf_b.raw()) == 0) {
            return false;
        }
        const silex::OrderElement* generator_a =
                ideal_a.known_principal_generator();
        const silex::OrderElement* generator_b =
                ideal_b.known_principal_generator();
        if ((generator_a == nullptr) != (generator_b == nullptr)) {
            return false;
        }
        if (generator_a != nullptr) {
            sflint::FmpzMat coordinates_a(1, n);
            sflint::FmpzMat coordinates_b(1, n);
            if (!generator_a->get_coordinates(
                        sflint::FmpzMatRef(coordinates_a)) ||
                !generator_b->get_coordinates(
                        sflint::FmpzMatRef(coordinates_b)) ||
                fmpz_mat_equal(coordinates_a.raw(), coordinates_b.raw()) ==
                        0) {
                return false;
            }
        }
        sflint::FmpzPoly residue_a;
        sflint::FmpzPoly residue_b;
        const bool has_residue_a =
                a.residue_polynomial(sflint::FmpzPolyRef(residue_a));
        if (has_residue_a !=
                    b.residue_polynomial(sflint::FmpzPolyRef(residue_b)) ||
            (has_residue_a &&
             fmpz_poly_equal(residue_a.raw(), residue_b.raw()) == 0)) {
            return false;
        }
        sflint::FmpzMat kummer_a(1, n);
        sflint::FmpzMat kummer_b(1, n);
        const bool has_kummer_a =
                a.kummer_generator_coordinates(sflint::FmpzMatRef(kummer_a));
        if (has_kummer_a != b.kummer_generator_coordinates(
                                    sflint::FmpzMatRef(kummer_b)) ||
            (has_kummer_a &&
             fmpz_mat_equal(kummer_a.raw(), kummer_b.raw()) == 0)) {
            return false;
        }
    }
    return true;
}

// decompose_prime keeps the whole decomposition of an index-divisor prime in
// the order's shared data and applies the residue-degree bound on read.  For
// each p, calls on one order with the bounds in `bounds`, in turn (the first
// fills the cache, possibly under a bound that drops primes), must give what
// the first call on a fresh maximal order gives for the same bound, prime
// order included, and the cache must hold the unbounded decomposition.
// A stored prime without a principal generator is rebuilt from its HNF
// without the release-build ideal test, so each such HNF must also pass the
// checked Ideal::set_hnf, and an unbounded read must publish that ideal.
// When `hnf_records` is not null it receives the number of such records.
int check_index_divisor_cache_matches_fresh(
        const std::vector<slong>& low_to_high,
        const std::vector<slong>& rational_primes,
        const std::vector<slong>& bounds,
        slong* hnf_records = nullptr) {
    slong hnf_record_count = 0;
    sflint::FmpqPoly polynomial;
    poly_from_coefficients(polynomial, low_to_high);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order cached(field);
    assert(cached.maximal_order(equation));
    assert(!cached.is_equation_order());

    for (const slong rational_prime : rational_primes) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, rational_prime));
        assert(silex::detail::OrderAccess::index_divisor_decomposition(
                       cached, sflint::FmpzConstRef(p)) == nullptr);

        silex::Order unbounded_order(field);
        assert(unbounded_order.maximal_order(equation));
        silex::PrimeIdealList unbounded;
        assert(silex::decompose_prime(unbounded, unbounded_order,
                                      sflint::FmpzConstRef(p)));

        for (const slong bound : bounds) {
            silex::PrimeIdealList from_cache;
            assert(silex::decompose_prime(from_cache, cached,
                                          sflint::FmpzConstRef(p), bound));
            const std::vector<silex::detail::IndexDivisorPrimeRecord>* stored =
                    silex::detail::OrderAccess::index_divisor_decomposition(
                            cached, sflint::FmpzConstRef(p));
            assert(stored != nullptr);
            assert(static_cast<slong>(stored->size()) == unbounded.size());
            // Order::set starts with nothing stored (the claim at
            // OrderAccess::index_divisor_decomposition), although `cached`
            // now holds a decomposition for p.
            silex::Order copy(field);
            assert(copy.set(cached));
            assert(silex::detail::OrderAccess::index_divisor_decomposition(
                           copy, sflint::FmpzConstRef(p)) == nullptr);

            silex::Order fresh(field);
            assert(fresh.maximal_order(equation));
            assert(!fresh.has_same_data(cached));
            silex::PrimeIdealList direct;
            assert(silex::decompose_prime(direct, fresh,
                                          sflint::FmpzConstRef(p), bound));
            assert(same_decomposition(from_cache, cached, direct, fresh));
        }

        const std::vector<silex::detail::IndexDivisorPrimeRecord>* stored =
                silex::detail::OrderAccess::index_divisor_decomposition(
                        cached, sflint::FmpzConstRef(p));
        assert(stored != nullptr);
        silex::PrimeIdealList all;
        assert(silex::decompose_prime(all, cached, sflint::FmpzConstRef(p)));
        assert(all.size() == static_cast<slong>(stored->size()));
        for (slong i = 0; i < all.size(); ++i) {
            const silex::detail::IndexDivisorPrimeRecord& record =
                    (*stored)[static_cast<std::size_t>(i)];
            if (record.has_principal_generator) {
                continue;
            }
            ++hnf_record_count;
            silex::Ideal checked(cached);
            assert(checked.set_hnf(sflint::FmpzMatConstRef(record.hnf)));
            silex::Ideal published(cached);
            assert(all.at(i)->get_ideal(published));
            assert(published.equal(checked));
        }
    }
    if (hnf_records != nullptr) {
        *hnf_records = hnf_record_count;
    }
    return 0;
}

// Q(sqrt(-23), sqrt(-31), sqrt(-47)): 2 splits into eight primes and 3 into
// four primes of degree two.  Both have more primes than p, so no separating
// fixed element exists and find_separating_fixed_element returns at once;
// p^g <= 4096, so the bounded enumeration splits them.  Each cached p is
// keyed separately on the one order.
int test_index_divisor_cache_more_primes_than_p() {
    slong hnf_records = 0;
    assert(check_index_divisor_cache_matches_fresh(
                   {38188288, 1822160, -1462612, -62792, 24143, 602, -220,
                    -2, 1},
                   {2, 3}, {1, 0, 2}, &hnf_records) == 0);
    // 8 + 4 primes, none with a known principal generator.
    assert(hnf_records == 12);
    return 0;
}

// Q(sqrt(-3), sqrt(-5), sqrt(-7)), p = 2: 2 ramifies (through sqrt(-5)) and
// divides the index, so the bounded ramified split builds the primes from
// the p-radical quotient and the cache rebuilds them from their HNFs.
int test_index_divisor_cache_ramified_hnf() {
    slong hnf_records = 0;
    assert(check_index_divisor_cache_matches_fresh(
                   {1201, 598, -524, -146, 69, 2, 4, -4, 1}, {2}, {1, 0, 2},
                   &hnf_records) == 0);
    // Two primes over 2 are cached, none with a principal
    // generator.
    assert(hnf_records == 2);
    return 0;
}

// Q(sqrt(-3), sqrt(-5), sqrt(-7)), p = 11: 11^4 > 4096, so the bounded
// enumeration does not apply; the recursive split decomposes p only if the
// monogenic generator search also fails for this order.  No existing probe
// observes that branch from this test (the profile scopes need a diagnostics
// build), so the test does not assert it.
int test_index_divisor_cache_recursive_split() {
    return check_index_divisor_cache_matches_fresh(
            {1201, 598, -524, -146, 69, 2, 4, -4, 1}, {11}, {1, 0, 2});
}

// The monogenic Dedekind-Kummer branch: 17^3 > 4096 primes of degree one,
// with and without a ramified prime, and 17^4 g(x/17) for
// g = y (y - 1) (y^2 + 3) + 17, where 17 has primes of degree 1, 1 and 2.
// The first call under bound 1 drops the degree-two prime.
int test_index_divisor_cache_monogenic() {
    assert(check_index_divisor_cache_matches_fresh({54043, 3179, -102, 1},
                                                   {17}, {1, 0}) == 0);
    assert(check_index_divisor_cache_matches_fresh(
                   {1419857, -19652, 2312, -85, 1}, {17}, {1, 0}) == 0);
    return check_index_divisor_cache_matches_fresh(
            {1419857, -14739, 867, -17, 1}, {17}, {1, 0, 2});
}

// x^3 + 4x + 8 = 8 g(x/2) with g = y^3 + y + 1 irreducible modulo 2 and of
// odd discriminant: 2 divides the index and is inert, so its one prime is
// 2O, built with a known principal generator, which the cache keeps.
int test_index_divisor_cache_inert_principal() {
    return check_index_divisor_cache_matches_fresh({8, 4, 0, 1}, {2},
                                                   {1, 0, 3});
}

bool same_index_divisor_records(
        const std::vector<silex::detail::IndexDivisorPrimeRecord>& left,
        const std::vector<silex::detail::IndexDivisorPrimeRecord>& right) {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        const silex::detail::IndexDivisorPrimeRecord& a = left[i];
        const silex::detail::IndexDivisorPrimeRecord& b = right[i];
        if (a.ramification_index != b.ramification_index ||
            a.residue_degree != b.residue_degree ||
            a.has_principal_generator != b.has_principal_generator ||
            fmpz_mat_equal(a.hnf.raw(), b.hnf.raw()) == 0 ||
            (a.has_principal_generator &&
             fmpz_mat_equal(a.principal_generator.raw(),
                            b.principal_generator.raw()) == 0)) {
            return false;
        }
    }
    return true;
}

std::vector<slong> published_residue_degrees(
        const silex::PrimeIdealList& primes) {
    std::vector<slong> degrees;
    for (slong i = 0; i < primes.size(); ++i) {
        degrees.push_back(primes.at(i)->residue_degree());
    }
    return degrees;
}

// detail::decompose_prime_residue_degrees must give the residue degrees
// decompose_prime publishes, in its order, on a cache miss, on a hit after
// its own miss, and on a hit after decompose_prime filled the cache.  Its
// miss must store the records decompose_prime stores on a fresh order, so a
// later decompose_prime on that order publishes the same primes.  Primes in
// `index_primes` must be cached index divisors; the others may take any
// branch.
int check_residue_degrees_match_decompose_prime(
        const std::vector<slong>& low_to_high,
        const std::vector<slong>& rational_primes,
        const std::vector<slong>& index_primes) {
    sflint::FmpqPoly polynomial;
    poly_from_coefficients(polynomial, low_to_high);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    for (const slong rational_prime : rational_primes) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, rational_prime));
        const sflint::FmpzConstRef p_ref(p);

        silex::Order reference(field);
        assert(reference.maximal_order(equation));
        silex::PrimeIdealList expected_primes;
        assert(silex::decompose_prime(expected_primes, reference, p_ref));
        const std::vector<slong> expected =
                published_residue_degrees(expected_primes);
        const std::vector<silex::detail::IndexDivisorPrimeRecord>*
                reference_records =
                        silex::detail::OrderAccess::index_divisor_decomposition(
                                reference, p_ref);
        bool is_index_prime = false;
        for (const slong index_prime : index_primes) {
            is_index_prime = is_index_prime || index_prime == rational_prime;
        }
        assert(!is_index_prime || reference_records != nullptr);

        // Miss, then hit, on one order.
        silex::Order order(field);
        assert(order.maximal_order(equation));
        assert(silex::detail::OrderAccess::index_divisor_decomposition(
                       order, p_ref) == nullptr);
        std::vector<slong> degrees;
        assert(silex::detail::decompose_prime_residue_degrees(degrees, order,
                                                              p_ref));
        assert(degrees == expected);
        const std::vector<silex::detail::IndexDivisorPrimeRecord>* stored =
                silex::detail::OrderAccess::index_divisor_decomposition(
                        order, p_ref);
        assert((stored == nullptr) == (reference_records == nullptr));
        if (stored != nullptr) {
            assert(same_index_divisor_records(*stored, *reference_records));
        }
        degrees.clear();
        assert(silex::detail::decompose_prime_residue_degrees(degrees, order,
                                                              p_ref));
        assert(degrees == expected);
        silex::PrimeIdealList after;
        assert(silex::decompose_prime(after, order, p_ref));
        assert(same_decomposition(after, order, expected_primes, reference));

        // Hit after decompose_prime filled the cache.
        degrees.clear();
        assert(silex::detail::decompose_prime_residue_degrees(degrees,
                                                              reference,
                                                              p_ref));
        assert(degrees == expected);
    }
    return 0;
}

// The index divisors 43, 67 and 1319 of the field of
// cyclotomic_real_q10_c11_deg10 split into ten primes of degree one, past the
// bounded enumeration, as do the discriminant divisors 2, 5 and 11, which
// decompose_prime handles outside its index-divisor branch.
int test_residue_degrees_cyclotomic_real_q10_c11() {
    return check_residue_degrees_match_decompose_prime(
            {238831, -228566, -128567, 146666, -15113, -13352, 2746, 334,
             -97, -2, 1},
            {43, 67, 1319, 2, 5, 11}, {43, 67, 1319});
}

// Bounded enumeration (2 and 3 in Q(sqrt(-23), sqrt(-31), sqrt(-47))), the
// recursive split or the monogenic branch (11 in Q(sqrt(-3), sqrt(-5),
// sqrt(-7))), the monogenic branch with an unramified and a ramified 17, and
// an inert index divisor with a principal generator.
int test_residue_degrees_index_divisor_routes() {
    assert(check_residue_degrees_match_decompose_prime(
                   {38188288, 1822160, -1462612, -62792, 24143, 602, -220, -2,
                    1},
                   {2, 3}, {2, 3}) == 0);
    assert(check_residue_degrees_match_decompose_prime(
                   {1201, 598, -524, -146, 69, 2, 4, -4, 1}, {11}, {11}) ==
           0);
    assert(check_residue_degrees_match_decompose_prime({54043, 3179, -102, 1},
                                                       {17}, {17}) == 0);
    assert(check_residue_degrees_match_decompose_prime(
                   {1419857, -19652, 2312, -85, 1}, {17}, {17}) == 0);
    return check_residue_degrees_match_decompose_prime({8, 4, 0, 1}, {2},
                                                       {2});
}

// The same prime, built the same way: ideal HNF, e, f, residue polynomial
// and Kummer generator.
bool same_prime_data(const silex::PrimeIdeal& a,
                     const silex::PrimeIdeal& b,
                     const silex::Order& order) {
    const slong n = order.degree();
    sflint::FmpzMat hnf_a(n, n);
    sflint::FmpzMat hnf_b(n, n);
    sflint::FmpzMat kummer_a(1, n);
    sflint::FmpzMat kummer_b(1, n);
    sflint::FmpzPoly residue_a;
    sflint::FmpzPoly residue_b;
    const bool has_kummer_a =
            a.kummer_generator_coordinates(sflint::FmpzMatRef(kummer_a));
    const bool has_kummer_b =
            b.kummer_generator_coordinates(sflint::FmpzMatRef(kummer_b));
    const bool has_residue_a =
            a.residue_polynomial(sflint::FmpzPolyRef(residue_a));
    const bool has_residue_b =
            b.residue_polynomial(sflint::FmpzPolyRef(residue_b));
    return a.ramification_index() == b.ramification_index() &&
           a.residue_degree() == b.residue_degree() &&
           PrimeIdealAccess::ideal(a).get_hnf(sflint::FmpzMatRef(hnf_a)) &&
           PrimeIdealAccess::ideal(b).get_hnf(sflint::FmpzMatRef(hnf_b)) &&
           fmpz_mat_equal(hnf_a.raw(), hnf_b.raw()) != 0 &&
           has_kummer_a == has_kummer_b &&
           (!has_kummer_a ||
            fmpz_mat_equal(kummer_a.raw(), kummer_b.raw()) != 0) &&
           has_residue_a == has_residue_b &&
           (!has_residue_a ||
            fmpz_poly_equal(residue_a.raw(), residue_b.raw()) != 0);
}

// decompose_prime_containing on `element`: on a Dedekind-Kummer route
// (`filtered`) exactly the primes of decompose_prime(.., p, 0) that contain
// the element, in that order and built identically; otherwise all of them.
void check_decompose_prime_containing(const silex::Order& order,
                                      const silex::PrimeIdealList& all,
                                      const silex::OrderElement& element,
                                      sflint::FmpzConstRef p,
                                      bool filtered) {
    std::vector<silex::PrimeIdeal> containing;
    assert(PrimeIdealAccess::decompose_prime_containing(containing, element,
                                                        p, nullptr));
    std::size_t next = 0;
    for (slong i = 0; i < all.size(); ++i) {
        const silex::PrimeIdeal& prime = *all.at(i);
        slong valuation = -1;
        assert(prime.valuation(valuation, element));
        if (filtered && valuation == 0) {
            continue;
        }
        assert(next < containing.size());
        assert(same_prime_data(containing[next], prime, order));
        ++next;
    }
    assert(next == containing.size());
}

// Elements of the maximal order of the field of `low_to_high`: small
// coordinate vectors, p, 1, and each prime's Kummer generator, which lies
// in that prime.  `index_primes` divide [O : Z[theta]] and keep the whole
// decomposition.
int check_decompose_prime_containing_routes(
        const std::vector<slong>& low_to_high,
        bool maximal_is_equation_order,
        const std::vector<slong>& rational_primes,
        const std::vector<slong>& index_primes) {
    sflint::FmpqPoly polynomial;
    poly_from_coefficients(polynomial, low_to_high);
    silex::NumberField field;
    silex::Order equation = order_by_polynomial(field, polynomial);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());
    assert(order.is_equation_order() == maximal_is_equation_order);
    const slong n = order.degree();

    for (const slong rational_prime : rational_primes) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, rational_prime));
        const sflint::FmpzConstRef p_ref(p);
        bool filtered = true;
        for (const slong index_prime : index_primes) {
            filtered = filtered && index_prime != rational_prime;
        }
        silex::PrimeIdealList all;
        assert(silex::decompose_prime(all, order, p_ref));
        assert(all.size() > 0);

        std::vector<silex::OrderElement> elements;
        sflint::FmpzMat coordinates(1, n);
        for (slong k = 0; k < 2 * n + 2; ++k) {
            for (slong j = 0; j < n; ++j) {
                fmpz_set_si(fmpz_mat_entry(coordinates.raw(), 0, j),
                            ((k + 1) * (j + 3) * (j + k + 1)) % 7 - 3);
            }
            if (fmpz_mat_is_zero(coordinates.raw()) != 0) {
                continue;
            }
            elements.emplace_back(order);
            assert(elements.back().set_coordinates(
                    sflint::FmpzMatConstRef(coordinates)));
        }
        elements.emplace_back(order);
        assert(elements.back().set_si(rational_prime));
        elements.emplace_back(order);
        assert(elements.back().set_si(1));
        for (slong i = 0; i < all.size(); ++i) {
            if (all.at(i)->kummer_generator_coordinates(
                        sflint::FmpzMatRef(coordinates)) &&
                fmpz_mat_is_zero(coordinates.raw()) == 0) {
                elements.emplace_back(order);
                assert(elements.back().set_coordinates(
                        sflint::FmpzMatConstRef(coordinates)));
            }
        }
        for (const silex::OrderElement& element : elements) {
            check_decompose_prime_containing(order, all, element, p_ref,
                                             filtered);
        }

    }
    return 0;
}

// The route decompose_prime takes, from the classifier that
// decompose_prime, decompose_prime_containing and
// decompose_prime_residue_degrees share, pinned to the routes each took
// before they shared it.
silex::detail::DecomposePrimeRoute route_at(const silex::Order& order,
                                            slong rational_prime,
                                            slong max_residue_degree) {
    sflint::Fmpz p;
    assert(set_fmpz_si(p, rational_prime));
    return silex::detail::TestHookAccess::decompose_prime_route(
            order, sflint::FmpzConstRef(p), max_residue_degree);
}

// decompose_prime_containing's use of the classifier, on the maximal order
// `order`: for the element 1 it builds no prime exactly on the
// Dedekind-Kummer routes, and the whole decomposition otherwise.
void check_containing_route(const silex::Order& order, slong rational_prime) {
    using Route = silex::detail::DecomposePrimeRoute;
    sflint::Fmpz p;
    assert(set_fmpz_si(p, rational_prime));
    const sflint::FmpzConstRef p_ref(p);
    silex::OrderElement one(order);
    assert(one.set_si(1));
    std::vector<silex::PrimeIdeal> containing;
    assert(PrimeIdealAccess::decompose_prime_containing(containing, one, p_ref,
                                                        nullptr));
    const Route route = route_at(order, rational_prime, 0);
    if (route == Route::nonindex_factor || route == Route::equation_factor) {
        assert(containing.empty());
    } else {
        silex::PrimeIdealList all;
        assert(silex::decompose_prime(all, order, p_ref));
        assert(static_cast<slong>(containing.size()) == all.size());
        assert(!containing.empty());
    }
}

int test_decompose_prime_route() {
    using Route = silex::detail::DecomposePrimeRoute;

    // Q: the degree-one order.
    sflint::FmpqPoly linear;
    poly_x(linear);
    silex::NumberField rationals;
    const silex::Order integers = order_by_polynomial(rationals, linear);
    assert(route_at(integers, 7, 0) == Route::degree_one);
    assert(route_at(integers, 7, 1) == Route::degree_one);
    assert(route_at(integers, 8, 0) == Route::refused);
    assert(route_at(integers, 7, -1) == Route::refused);

    // Q(sqrt(5)): the maximal order Z[(1 + sqrt(5)) / 2] takes the quadratic
    // integral-generator route for every bound; its equation order Z[sqrt(5)]
    // is not maximal and factors x^2 - 5.
    sflint::FmpqPoly quadratic;
    poly_from_coefficients(quadratic, {-5, 0, 1});
    silex::NumberField quadratic_field;
    const silex::Order quadratic_equation =
            order_by_polynomial(quadratic_field, quadratic);
    silex::Order quadratic_maximal(quadratic_field);
    assert(quadratic_maximal.maximal_order(quadratic_equation));
    assert(!quadratic_maximal.is_equation_order());
    for (const slong q : {2, 5, 11, 19}) {
        assert(route_at(quadratic_maximal, q, 0) ==
               Route::quadratic_integral_generator);
        assert(route_at(quadratic_maximal, q, 1) ==
               Route::quadratic_integral_generator);
        check_containing_route(quadratic_maximal, q);
    }
    assert(route_at(quadratic_maximal, 9, 0) == Route::refused);
    assert(route_at(quadratic_equation, 11, 0) == Route::equation_factor);
    assert(route_at(quadratic_equation, 11, 1) == Route::squarefree_roots);
    // 2 divides disc(x^2 - 5) = 20: no roots shortcut.
    assert(route_at(quadratic_equation, 2, 1) == Route::equation_factor);

    // x^3 - 2: the maximal order is the equation order.
    sflint::FmpqPoly cubic;
    poly_from_coefficients(cubic, {-2, 0, 0, 1});
    silex::NumberField cubic_field;
    const silex::Order cubic_equation =
            order_by_polynomial(cubic_field, cubic);
    silex::Order cubic_maximal(cubic_field);
    assert(cubic_maximal.maximal_order(cubic_equation));
    assert(cubic_maximal.is_equation_order());
    for (const slong q : {2, 3, 5, 31}) {
        assert(route_at(cubic_maximal, q, 0) == Route::equation_factor);
        check_containing_route(cubic_maximal, q);
    }
    assert(route_at(cubic_maximal, 31, 1) == Route::squarefree_roots);
    assert(route_at(cubic_maximal, 5, 1) == Route::squarefree_roots);
    // 2 and 3 divide disc(x^3 - 2) = -108.
    assert(route_at(cubic_maximal, 2, 1) == Route::equation_factor);
    assert(route_at(cubic_maximal, 3, 1) == Route::equation_factor);
    assert(route_at(cubic_maximal, 31, 2) == Route::equation_factor);

    // cyclotomic_real_q10_c11_deg10: a non-equation maximal order with
    // index divisor 43.
    sflint::FmpqPoly decic;
    poly_from_coefficients(decic, {238831, -228566, -128567, 146666, -15113,
                                   -13352, 2746, 334, -97, -2, 1});
    silex::NumberField decic_field;
    const silex::Order decic_equation =
            order_by_polynomial(decic_field, decic);
    silex::Order decic_maximal(decic_field);
    assert(decic_maximal.maximal_order(decic_equation));
    assert(!decic_maximal.is_equation_order());
    for (const slong q : {2, 3, 5, 11, 23, 89, 101}) {
        assert(route_at(decic_maximal, q, 0) == Route::nonindex_factor);
        check_containing_route(decic_maximal, q);
    }
    assert(route_at(decic_maximal, 43, 0) == Route::index_divisor);
    check_containing_route(decic_maximal, 43);
    assert(route_at(decic_maximal, 43, 1) == Route::index_divisor);
    assert(route_at(decic_maximal, 101, 2) == Route::nonindex_factor);
    assert(route_at(decic_equation, 43, 0) == Route::equation_factor);
    // The bound-1 squarefree test of a maximal non-equation order: 101 does
    // not divide the index or the discriminant.
    assert(route_at(decic_maximal, 101, 1) == Route::squarefree_roots);

    // A non-equation order not known to be maximal is refused, whatever the
    // prime (Q(sqrt(5)) with basis 1, 3 (1 + theta) / 2).
    silex::Order unknown_maximality(quadratic_field);
    sflint::FmpqMat basis(2, 2);
    qmat_entry_si(basis, 0, 0, 1);
    qmat_entry_frac_si(basis, 1, 0, 3, 2);
    qmat_entry_frac_si(basis, 1, 1, 3, 2);
    assert(unknown_maximality.set_basis(sflint::FmpqMatConstRef(basis)));
    assert(!unknown_maximality.is_equation_order());
    assert(!unknown_maximality.is_maximal());
    assert(route_at(unknown_maximality, 3, 0) == Route::refused);
    assert(route_at(unknown_maximality, 11, 0) == Route::refused);
    assert(route_at(unknown_maximality, 11, 1) == Route::refused);

    // decompose_prime_residue_degrees takes its own index-divisor path only
    // on the index_divisor route: elsewhere it must leave the order's
    // index-divisor cache empty.
    const auto check_no_index_path = [](const silex::Order& order,
                                        slong rational_prime) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, rational_prime));
        const sflint::FmpzConstRef p_ref(p);
        silex::PrimeIdealList all;
        assert(silex::decompose_prime(all, order, p_ref));
        const std::vector<slong> expected = published_residue_degrees(all);
        std::vector<slong> degrees;
        assert(silex::detail::decompose_prime_residue_degrees(degrees, order,
                                                              p_ref));
        assert(degrees == expected);
        assert(silex::detail::OrderAccess::index_divisor_decomposition(
                       order, p_ref) == nullptr);
    };
    for (const slong q : {2, 3, 5, 11, 23, 89, 101}) {
        silex::Order fresh(decic_field);
        assert(fresh.maximal_order(decic_equation));
        check_no_index_path(fresh, q);
    }
    for (const slong q : {2, 3, 5, 31}) {
        silex::Order fresh(cubic_field);
        assert(fresh.maximal_order(cubic_equation));
        check_no_index_path(fresh, q);
    }
    return 0;
}

// The non-equation maximal order of cyclotomic_real_q10_c11_deg10 at
// non-index primes (split, inert-like and ramified) and at the index divisor
// 43; x^3 - 2, whose maximal order is its equation order; the quadratic
// integral-generator route of Q(sqrt(5)).
int test_decompose_prime_containing() {
    assert(check_decompose_prime_containing_routes(
                   {238831, -228566, -128567, 146666, -15113, -13352, 2746,
                    334, -97, -2, 1},
                   false, {2, 3, 5, 11, 23, 43, 89, 101}, {43}) == 0);
    assert(check_decompose_prime_containing_routes({-2, 0, 0, 1}, true,
                                                   {2, 3, 5, 31}, {}) == 0);
    assert(check_decompose_prime_containing_routes({-5, 0, 1}, false,
                                                   {2, 5, 11, 19}, {2, 5, 11,
                                                                    19}) ==
           0);
    return 0;
}

// Splits O/(p, I_p) through the bounded idempotent enumeration and through
// the recursive split (both without the separating element) and checks that
// they give the same prime ideals in the same order, which is also the order
// decompose_prime publishes.  p is smaller than the number of primes above
// it, so a fixed element takes fewer values than there are fields and the
// recursive split needs more than one level.
int check_recursive_split_matches_enumeration(
        const std::vector<slong>& low_to_high,
        slong rational_prime,
        slong count) {
    sflint::FmpqPoly polynomial;
    poly_from_coefficients(polynomial, low_to_high);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(!maximal.is_equation_order());
    assert(rational_prime < count);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, rational_prime));
    std::vector<silex::Ideal> enumerated;
    std::vector<silex::Ideal> split;
    assert(silex::detail::TestHookAccess::index_divisor_prime_ideals(
            enumerated, maximal, sflint::FmpzConstRef(p), false));
    assert(silex::detail::TestHookAccess::index_divisor_prime_ideals(
            split, maximal, sflint::FmpzConstRef(p), true));
    assert(static_cast<slong>(enumerated.size()) == count);
    assert(split.size() == enumerated.size());

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == count);
    for (std::size_t i = 0; i < split.size(); ++i) {
        assert(split[i].equal(enumerated[i]));
        // Independent closure check that survives optimized library builds:
        // the checked set_hnf must accept the stored lattice unchanged.
        {
            sflint::FmpzMat stored(maximal.degree(), maximal.degree());
            assert(split[i].get_hnf(sflint::FmpzMatRef(stored)));
            silex::Ideal reloaded(maximal);
            assert(reloaded.set_hnf(sflint::FmpzMatConstRef(stored)));
            assert(reloaded.equal(split[i]));
        }
        silex::Ideal ideal(maximal);
        assert(primes.at(static_cast<slong>(i))->get_ideal(ideal));
        assert(ideal.equal(split[i]));
        for (std::size_t j = 0; j < i; ++j) {
            assert(!split[j].equal(split[i]));
        }
    }
    return 0;
}

int test_recursive_split_matches_enumeration_unramified() {
    // Q(sqrt(-7), sqrt(17)), p = 2: four primes with e = f = 1.
    assert(check_recursive_split_matches_enumeration({36, 0, -5, 0, 1}, 2,
                                                     4) == 0);
    // Q(sqrt(-2), sqrt(7)), p = 3: four primes with e = f = 1.
    return check_recursive_split_matches_enumeration({9, 0, 8, 0, 1}, 3, 4);
}

int test_recursive_split_matches_enumeration_ramified() {
    // Q(sqrt(-3), sqrt(7), sqrt(-2)), p = 3: four primes with e = 2, f = 1.
    return check_recursive_split_matches_enumeration(
            {81, 0, -72, 0, 55, 0, -8, 0, 1}, 3, 4);
}

// x^3 - c p^k with p > 2^64: the maximal order is not the equation order
// and p divides the index, so the index-divisor routes decompose p with the
// multiplication table modulo p held as fmpz residues (p does not fit a
// word).
silex::Order maximal_order_of_pure_cubic(silex::NumberField& field,
                                         slong c,
                                         const sflint::Fmpz& p,
                                         ulong k) noexcept {
    sflint::Fmpz constant;
    fmpz_pow_ui(constant.raw(), p.raw(), k);
    fmpz_mul_si(constant.raw(), constant.raw(), -c);
    sflint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, constant.raw());
    fmpq_poly_set_coeff_si(polynomial.raw(), 3, 1);
    silex::Order equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(!maximal.is_equation_order());
    return maximal;
}

// p = 18446744073709551667 > 2^64 is 1 mod 3 and 2 is a cube mod p, so p
// splits into three degree-one primes of Q(cbrt(2)), here given by
// x^3 - 2 p^3 (theta = p cbrt(2)).  The recursive etale split (Frobenius
// kernel and splitting products modulo p) must give, as a set, the primes
// decompose_prime finds through the monogenic generator cbrt(2).
int check_index_divisor_recursive_split_cubic(const char* prime_text,
                                              bool fits_word) {
    sflint::Fmpz p;
    fmpz_set_str(p.raw(), prime_text, 10);
    assert((fmpz_abs_fits_ui(p.raw()) != 0) == fits_word);
    silex::NumberField field;
    const silex::Order maximal = maximal_order_of_pure_cubic(field, 2, p, 3);

    std::vector<silex::Ideal> split;
    assert(silex::detail::TestHookAccess::index_divisor_prime_ideals(
            split, maximal, sflint::FmpzConstRef(p), true));
    assert(split.size() == 3);

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 3);
    assert(count_with_ef(primes, 1, 1) == 3);
    assert(check_norms(primes));
    for (const silex::Ideal& ideal : split) {
        slong matches = 0;
        for (slong i = 0; i < primes.size(); ++i) {
            silex::Ideal prime_ideal(maximal);
            assert(primes.at(i)->get_ideal(prime_ideal));
            matches += prime_ideal.equal(ideal) ? 1 : 0;
        }
        assert(matches == 1);
    }
    return 0;
}

int test_index_divisor_recursive_split_multiword_prime() {
    return check_index_divisor_recursive_split_cubic("18446744073709551667",
                                                     false);
}

// Same with p = 18446744073709551337, 2^63 < p < 2^64: p fits one word with
// nmod norm 0, so the word path of order_coords_mul_mod runs on a
// full-word modulus.
int test_index_divisor_recursive_split_full_word_prime() {
    return check_index_divisor_recursive_split_cubic("18446744073709551337",
                                                     true);
}

// x^3 - 3 p^4 with the same p: theta = p cbrt(3 p), and p is totally
// ramified in Q(cbrt(3 p)).  The ramified index-divisor route computes the
// p-radical as the kernel of x -> x^p modulo p; the prime it builds must
// be the one decompose_prime publishes (e = 3, f = 1, norm p), with cube
// p O.
int test_index_divisor_ramified_multiword_prime() {
    sflint::Fmpz p;
    fmpz_set_str(p.raw(), "18446744073709551667", 10);
    silex::NumberField field;
    const silex::Order maximal = maximal_order_of_pure_cubic(field, 3, p, 4);

    std::vector<silex::Ideal> radical;
    assert(silex::detail::TestHookAccess::index_divisor_prime_ideals(
            radical, maximal, sflint::FmpzConstRef(p), false));
    assert(radical.size() == 1);

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    assert(count_with_ef(primes, 3, 1) == 1);
    assert(check_norms(primes));
    sflint::Fmpz norm;
    assert(primes.at(0)->norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal(norm, p));

    silex::Ideal prime_ideal(maximal);
    assert(primes.at(0)->get_ideal(prime_ideal));
    assert(prime_ideal.equal(radical[0]));
    silex::Ideal cube(maximal);
    silex::Ideal square(maximal);
    assert(square.multiply(radical[0], radical[0]));
    assert(cube.multiply(square, radical[0]));
    sflint::FmpzMat p_hnf(3, 3);
    for (slong i = 0; i < 3; ++i) {
        fmpz_set(fmpz_mat_entry(p_hnf.raw(), i, i), p.raw());
    }
    silex::Ideal p_ideal(maximal);
    assert(p_ideal.set_hnf(sflint::FmpzMatConstRef(p_hnf)));
    assert(cube.equal(p_ideal));
    return 0;
}

int test_valuation_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::Ideal ideal(order);
    slong value = -1;
    assert(set_rational_principal(ideal, 1));
    assert(prime->valuation(value, ideal));
    assert(value == 0);

    value = -1;
    assert(set_rational_principal(ideal, 5));
    assert(prime->valuation(value, ideal));
    assert(value == 0);

    value = -1;
    assert(set_rational_principal(ideal, 7));
    assert(prime->valuation(value, ideal));
    assert(value == 1);

    silex::OrderElement generator(order);
    assert(generator.set_si(7));
    value = -1;
    assert(prime->valuation(value, generator));
    assert(value == 1);

    value = -1;
    assert(set_rational_principal(ideal, 343));
    assert(prime->valuation(value, ideal));
    assert(value == 3);

    assert(generator.set_si(343));
    value = -1;
    assert(prime->valuation(value, generator));
    assert(value == 3);

    return 0;
}

int test_valuation_order_element_split_repeated() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    silex::Element theta(field);
    silex::Element three(field);
    silex::Element theta_minus_three(field);
    silex::Element square(field);
    assert(theta.gen());
    assert(three.set_si(3));
    assert(theta_minus_three.subtract(theta, three));
    assert(square.multiply(theta_minus_three, theta_minus_three));

    silex::OrderElement alpha(order);
    assert(alpha.set_element(square));

    const silex::PrimeIdeal* root_three = nullptr;
    const silex::PrimeIdeal* root_four = nullptr;
    sflint::FmpzPoly reduced;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        assert(prime->ramification_index() == 1);
        assert(prime->residue_degree() == 1);
        assert(prime->reduce(sflint::FmpzPolyRef(reduced), theta));
        if (fmpz_poly_degree(reduced.raw()) == 0 &&
            fmpz_poly_coeff_is_si(reduced, 0, 3)) {
            root_three = prime;
        } else if (fmpz_poly_degree(reduced.raw()) == 0 &&
                   fmpz_poly_coeff_is_si(reduced, 0, 4)) {
            root_four = prime;
        }
    }
    assert(root_three != nullptr);
    assert(root_four != nullptr);

    slong value = -1;
    assert(root_three->valuation(value, alpha));
    assert(value == 2);

    slong repeated_value = -1;
    assert(root_three->valuation(repeated_value, alpha));
    assert(repeated_value == value);

    silex::Ideal ideal(order);
    slong ideal_value = -1;
    assert(ideal.set_principal(alpha));
    assert(root_three->valuation(ideal_value, ideal));
    assert(ideal_value == value);

    value = -1;
    assert(root_four->valuation(value, alpha));
    assert(value == 0);

    return 0;
}

int test_valuation_order_element_ramified_unique_prime_norm() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(prime->ramification_index() == 2);
    assert(prime->residue_degree() == 1);

    silex::Element theta(field);
    assert(theta.gen());
    silex::OrderElement alpha(order);
    assert(alpha.set_element(theta));

    slong value = 99;
    assert(prime->valuation(value, alpha));
    assert(value == 1);

    silex::Ideal ideal(order);
    assert(ideal.set_principal(alpha));
    slong ideal_value = -1;
    assert(prime->valuation(ideal_value, ideal));
    assert(ideal_value == value);

    assert(alpha.set_si(2));
    value = -1;
    assert(prime->valuation(value, alpha));
    assert(value == 2);

    return 0;
}

int test_valuation_order_element_ramified_split_coordinate_matrix() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 23));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order,
                                  sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    const silex::PrimeIdeal* ramified = nullptr;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        if (prime->ramification_index() == 2) {
            ramified = prime;
        }
    }
    assert(ramified != nullptr);
    assert(ramified->residue_degree() == 1);

    silex::OrderElement rational_p(order);
    assert(rational_p.set_si(23));
    slong value = -1;
    assert(ramified->valuation(value, rational_p));
    assert(value == 2);

    slong repeated_value = -1;
    assert(ramified->valuation(repeated_value, rational_p));
    assert(repeated_value == value);

    return 0;
}

int test_valuation_order_element_ramified_nonunique_containment() {
    sflint::FmpqPoly polynomial;
    poly_degree6_benchmark_field(polynomial);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());

    bool tested = false;
    sflint::Fmpz p;
    for (slong rational_prime = 2; rational_prime <= 31 && !tested;
         ++rational_prime) {
        assert(set_fmpz_si(p, rational_prime));
        silex::PrimeIdealList primes;
        if (!silex::decompose_prime(primes, maximal,
                                    sflint::FmpzConstRef(p))) {
            continue;
        }
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            if (prime == nullptr || prime->ramification_index() <= 1 ||
                prime->ramification_index() * prime->residue_degree() >=
                        maximal.degree()) {
                continue;
            }

            silex::OrderElement rational_p(maximal);
            silex::Ideal principal(maximal);
            assert(rational_p.set_si(rational_prime));
            assert(principal.set_principal(rational_p));

            slong direct_value = -1;
            slong ideal_value = -1;
            assert(prime->valuation(direct_value, rational_p));
            assert(prime->valuation(ideal_value, principal));
            assert(direct_value == ideal_value);
            assert(direct_value == prime->ramification_index());
            tested = true;
            break;
        }
    }

    assert(tested);
    return 0;
}

// Definitional valuation for the containment cross-check: the largest
// i <= bound with element in P^i, built from explicit ideal powers and
// Ideal::contains.  It shares Ideal::multiply and Ideal::contains with the
// code under test, so it checks the cached power sequence;
// test_valuation_power_containment_matches_gp pins the values themselves.
slong definitional_power_containment_valuation(
        const silex::Order& order,
        const silex::PrimeIdeal& prime,
        const silex::OrderElement& element,
        slong bound) {
    silex::Ideal prime_ideal(order);
    silex::Ideal power(order);
    assert(prime.get_ideal(prime_ideal));
    assert(power.set(prime_ideal));
    slong value = 0;
    while (value < bound && power.contains(element)) {
        ++value;
        silex::Ideal next(order);
        assert(next.multiply(power, prime_ideal));
        power.swap(next);
    }
    return value;
}

int test_valuation_power_containment_matches_definition() {
    // The maximal orders of x^4 - 2x^3 + 5x^2 - 4x + 474 and x^6 + 7 differ
    // from their equation orders, so the primes below include index
    // divisors and ramified primes with e*f < n.
    const std::vector<std::vector<slong>> fields = {
            {474, -4, 5, -2, 1},
            {7, 0, 0, 0, 0, 0, 1},
    };
    slong checked_positive = 0;
    for (const std::vector<slong>& coefficients : fields) {
        sflint::FmpqPoly polynomial;
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial,
                                           static_cast<slong>(i),
                                           coefficients[i]);
        }
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        const slong n = maximal.degree();

        for (slong rational_prime : {2, 3, 5, 7, 11}) {
            sflint::Fmpz p;
            assert(set_fmpz_si(p, rational_prime));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(primes, maximal,
                                          sflint::FmpzConstRef(p)));
            for (slong index = 0; index < primes.size(); ++index) {
                const silex::PrimeIdeal* prime = primes.at(index);
                assert(prime != nullptr);
                for (slong t = 0; t < 24; ++t) {
                    sflint::FmpzMat coordinates(1, n);
                    bool nonzero = false;
                    for (slong j = 0; j < n; ++j) {
                        const slong c = ((5 * t + 3 * j + j * j * t) % 9) - 4;
                        fmpz_set_si(fmpz_mat_entry(coordinates.raw(), 0, j),
                                    c);
                        nonzero = nonzero || c != 0;
                    }
                    if (!nonzero) {
                        continue;
                    }
                    silex::OrderElement base(maximal);
                    assert(base.set_coordinates(
                            sflint::FmpzMatConstRef(coordinates)));
                    // Powers of the sample and of p push the valuation
                    // above one so the cached P^k bases are exercised.
                    silex::OrderElement element(maximal);
                    assert(element.set(base));
                    for (slong k = 0; k < t % 3; ++k) {
                        silex::OrderElement next(maximal);
                        assert(next.multiply(element, base));
                        element = std::move(next);
                    }
                    if (t % 4 == 1) {
                        silex::OrderElement scale(maximal);
                        silex::OrderElement next(maximal);
                        assert(scale.set_si(rational_prime));
                        assert(next.multiply(element, scale));
                        element = std::move(next);
                    }
                    for (slong bound : {0, 1, 3, 7}) {
                        slong fast = -1;
                        assert(silex::detail::TestHookAccess::valuation_by_power_containment(
                                               fast, *prime, element,
                                               bound));
                        const slong expected =
                                definitional_power_containment_valuation(
                                        maximal, *prime, element, bound);
                        assert(fast == expected);
                        if (expected > 1) {
                            ++checked_positive;
                        }
                    }
                }
            }
        }
    }
    assert(checked_positive > 0);
    return 0;
}

// Independent oracle for the power-containment valuation, pinned prime by
// prime.  For each field f the table lists the primes P above p in
// {2, 3, 5, 7, 11} as idealprimedec(nf, p) returns them for nf = nfinit(f),
// each with its p, e, f and second generator g = P.gen[2] written as
// (sum c_i x^i) / d, so that P = p O + g O.  Each row is an element
// a = sum c_i x^i of Q[x]/(f) with its (prime index, v_P(a)) entries for
// v_P(a) > 0 from idealval(nf, a, P); every other listed P has v_P(a) = 0.
// The data are from an independent computer algebra system (gp 2.17.4),
// recorded by hand.  Each p also satisfies sum_P f_P v_P(a) = v_p(N(a)).
// Silex's primes above p are matched one to one to the listed primes by the
// ideal equality P = p O + g O, so a valuation assigned to the wrong prime of
// the same (p, e, f) fails.
struct GpPrime {
    slong p;
    slong e;
    slong f;
    slong denominator;
    std::vector<slong> generator;
};

struct GpValuationRow {
    std::vector<slong> element;
    std::vector<std::array<slong, 2>> positive;
};

bool set_order_element_from_coefficients(silex::OrderElement& out,
                                         const silex::NumberField& field,
                                         const std::vector<slong>& numerators,
                                         slong denominator) {
    sflint::FmpqPoly polynomial;
    for (std::size_t i = 0; i < numerators.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       numerators[i]);
    }
    fmpq_poly_scalar_div_si(polynomial.raw(), polynomial.raw(), denominator);
    silex::Element value(field);
    return value.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)) &&
           out.set_element(value);
}

int test_valuation_power_containment_matches_gp() {
    struct Field {
        std::vector<slong> polynomial;
        std::vector<GpPrime> primes;
        std::vector<GpValuationRow> rows;
    };
    const std::vector<Field> fields = {
            {{474, -4, 5, -2, 1},
             {
                     {2, 2, 1, 1, {0, 1}},
                     {2, 2, 1, 1, {1, 1}},
                     {3, 1, 1, 87, {-9, -113, 45, -1}},
                     {3, 1, 1, 1, {-1, 1}},
                     {3, 1, 1, 1, {0, 1}},
                     {3, 1, 1, 87, {183, 26, 42, 1}},
                     {5, 2, 2, 1, {2, -1, 1}},
                     {7, 1, 2, 87, {-39, -103, 108, -14}},
                     {7, 1, 2, 1, {5, 2, 1}},
                     {11, 1, 2, 1, {-3, -1, 1}},
                     {11, 1, 2, 1, {7, -1, 1}},
             },
             {
                     {{0, 1}, {{0, 1}, {4, 1}}},
                     {{-1, 1}, {{1, 1}, {3, 1}}},
                     {{2, 1, 1}, {{0, 1}, {1, 1}}},
                     {{3, 6}, {{2, 1}, {3, 4}, {4, 1}, {5, 1}}},
                     {{9, 6, 1}, {{1, 2}, {4, 4}}},
                     {{-25, 0, 5}, {{1, 2}, {6, 2}}},
                     {{15, -1, 0, 1}, {{2, 1}, {3, 1}, {4, 4}, {5, 1}}},
                     {{30, 20, 10}, {{0, 2}, {1, 7}, {3, 2}, {4, 1}, {6, 2}}},
                     {{1, 0, 3, 0, 3, 0, 1}, {{1, 9}}},
                     {{-7, 14}, {{2, 1}, {5, 1}, {7, 1}, {8, 1}}},
             }},
            {{7, 0, 0, 0, 0, 0, 1},
             {
                     {2, 1, 1, 2, {1, 3, 3, 1, 1, 1}},
                     {2, 1, 1, 2, {3, 1, 1, 1, 1, 1}},
                     {2, 1, 2, 2, {-2, 3, 3, 0, 1, 1}},
                     {2, 1, 2, 2, {2, 1, 1, 0, 1, 1}},
                     {3, 3, 2, 1, {1, 0, 1}},
                     {5, 1, 2, 1, {-2, -2, 1}},
                     {5, 1, 2, 1, {-2, 0, 1}},
                     {5, 1, 2, 1, {-2, 2, 1}},
                     {7, 6, 1, 1, {0, 1}},
                     {11, 1, 1, 1, {-4, 1}},
                     {11, 1, 1, 1, {4, 1}},
                     {11, 1, 2, 1, {5, -4, 1}},
                     {11, 1, 2, 1, {5, 4, 1}},
             },
             {
                     {{0, 1}, {{8, 1}}},
                     {{1, 1}, {{0, 2}, {1, 1}}},
                     {{1, 3, 3, 1}, {{0, 6}, {1, 3}}},
                     {{2, 4, 2}, {{0, 5}, {1, 3}, {2, 1}, {3, 1}}},
                     {{3, 1, 1}, {{2, 1}, {3, 1}}},
                     {{0, 49}, {{8, 13}}},
                     {{0, 0, 0, 1}, {{8, 3}}},
                     {{9, -3, 0, 0, 1}, {{9, 1}}},
                     {{6, 0, 3}, {{4, 3}}},
                     {{1, -2, 3, -2, 1}, {{2, 4}, {3, 2}}},
                     {{33, 15}, {{0, 1}, {1, 2}, {4, 3}}},
                     {{11, 0, 0, 0, 0, 1}, {{0, 1}, {1, 3}}},
             }},
    };
    // Above every expected valuation, so the bounded search returns v_P(a).
    const slong bound = 16;
    slong checked_positive = 0;
    for (const Field& data : fields) {
        sflint::FmpqPoly polynomial;
        for (std::size_t i = 0; i < data.polynomial.size(); ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                           data.polynomial[i]);
        }
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));

        std::vector<silex::OrderElement> elements;
        for (const GpValuationRow& row : data.rows) {
            silex::OrderElement element(maximal);
            assert(set_order_element_from_coefficients(element, field,
                                                       row.element, 1));
            elements.push_back(std::move(element));
        }

        for (slong rational_prime : {2, 3, 5, 7, 11}) {
            sflint::Fmpz p;
            assert(set_fmpz_si(p, rational_prime));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(primes, maximal,
                                          sflint::FmpzConstRef(p)));
            silex::OrderElement p_element(maximal);
            silex::Ideal p_ideal(maximal);
            assert(p_element.set_si(rational_prime));
            assert(p_ideal.set_principal(p_element));

            // gp_index[i] is the table prime equal to Silex's prime i.
            std::vector<slong> gp_index;
            std::vector<bool> used(data.primes.size(), false);
            for (slong index = 0; index < primes.size(); ++index) {
                const silex::PrimeIdeal* prime = primes.at(index);
                assert(prime != nullptr);
                silex::Ideal prime_ideal(maximal);
                assert(prime->get_ideal(prime_ideal));
                slong match = -1;
                for (std::size_t k = 0; k < data.primes.size(); ++k) {
                    const GpPrime& gp_prime = data.primes[k];
                    if (gp_prime.p != rational_prime) {
                        continue;
                    }
                    silex::OrderElement generator(maximal);
                    silex::Ideal generator_ideal(maximal);
                    silex::Ideal candidate(maximal);
                    assert(set_order_element_from_coefficients(
                            generator, field, gp_prime.generator,
                            gp_prime.denominator));
                    assert(generator_ideal.set_principal(generator));
                    assert(candidate.add(p_ideal, generator_ideal));
                    if (candidate.equal(prime_ideal)) {
                        assert(match < 0);
                        match = static_cast<slong>(k);
                    }
                }
                assert(match >= 0 && !used[static_cast<std::size_t>(match)]);
                used[static_cast<std::size_t>(match)] = true;
                const GpPrime& matched = data.primes[static_cast<std::size_t>(match)];
                assert(prime->ramification_index() == matched.e &&
                       prime->residue_degree() == matched.f);
                gp_index.push_back(match);
            }
            for (std::size_t k = 0; k < data.primes.size(); ++k) {
                assert(used[k] == (data.primes[k].p == rational_prime));
            }

            for (std::size_t r = 0; r < data.rows.size(); ++r) {
                for (slong index = 0; index < primes.size(); ++index) {
                    slong expected = 0;
                    for (const std::array<slong, 2>& entry :
                         data.rows[r].positive) {
                        if (entry[0] == gp_index[static_cast<std::size_t>(index)]) {
                            expected = entry[1];
                        }
                    }
                    slong v = -1;
                    assert(silex::detail::TestHookAccess::valuation_by_power_containment(
                                           v, *primes.at(index), elements[r],
                                           bound));
                    assert(v == expected);
                    if (v > 0) {
                        ++checked_positive;
                    }
                }
            }
        }
    }
    assert(checked_positive == 50);
    return 0;
}

// Differential sweep for the anti-uniformizer valuation route.  For every
// prime above p in a set of fields whose maximal orders differ from their
// equation orders (so index divisors appear, ramified or not) and fields
// with ramified primes, the coordinate-divisibility valuation with the
// anti-uniformizer taken from the HNF of P must equal the exact
// power-containment valuation at every bound, and so must the public
// element valuation, which now takes that route for primes without a Kummer
// quotient map (when its bound is at least two, the element lies in P, the
// norm shortcut does not apply, p fits in a word and P != pO; otherwise it
// falls back to containment).  Elements include products of basis elements of
// P so that valuations above one occur.  In a build with profiling enabled
// the sweep also asserts that the public valuation reached the
// path_anti_uniformizer event at least once.
struct AntiUniformizerRouteCounter {
    slong count = 0;
};

void anti_uniformizer_route_callback(void* user,
                                     silex::DiagnosticsModule,
                                     silex::ProfileEvent event,
                                     const char*,
                                     const char* label) noexcept {
    AntiUniformizerRouteCounter* counter =
            static_cast<AntiUniformizerRouteCounter*>(user);
    if (event == silex::ProfileEvent::event && label != nullptr &&
        std::strcmp(label,
                    "prime_ideal.valuation_order_element."
                    "path_anti_uniformizer") == 0) {
        ++counter->count;
    }
}

int test_valuation_anti_uniformizer_matches_containment() {
    const std::vector<std::vector<slong>> fields = {
            // x^3 - x^2 - 2x - 8: 2 is a common index divisor, split.
            {-8, -2, -1, 1},
            // x^3 - 10: 3 = P^2 Q divides the index.
            {-10, 0, 0, 1},
            // Q(sqrt(-7), sqrt(17)): 2 splits completely, index 2^3 * 3.
            {36, 0, -5, 0, 1},
            // Q(sqrt(2), sqrt(3)): 2 is ramified and divides the index.
            {1, 0, -10, 0, 1},
            // Quartic with a maximal order that is not the equation order.
            {474, -4, 5, -2, 1},
            // Sextic x^6 + 7 with a maximal order that is not the equation
            // order.
            {7, 0, 0, 0, 0, 0, 1},
            // polredabs of Q(sqrt(-23), sqrt(-31), sqrt(-47)): index
            // 2^20 * 3^2 * 2803; 23, 31 and 47 ramify.
            {34784, -16936, 19212, -4654, 2583, -310, 108, -4, 1},
    };
    slong checked_index_positive = 0;
    slong checked_ramified_positive = 0;
    slong checked_above_one = 0;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    AntiUniformizerRouteCounter route_counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(
                    silex::DiagnosticsModule::prime_ideal),
            anti_uniformizer_route_callback, &route_counter);
    for (const std::vector<slong>& coefficients : fields) {
        sflint::FmpqPoly polynomial;
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial,
                                           static_cast<slong>(i),
                                           coefficients[i]);
        }
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        const slong n = maximal.degree();
        sflint::Fmpz index;
        assert(equation.index_in(sflint::FmpzRef(index), maximal));

        for (slong rational_prime : {2, 3, 5, 7, 23, 31, 47}) {
            sflint::Fmpz p;
            assert(set_fmpz_si(p, rational_prime));
            const bool index_prime =
                    fmpz_divisible(index.raw(), p.raw()) != 0;
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(primes, maximal,
                                          sflint::FmpzConstRef(p)));
            for (slong prime_index = 0; prime_index < primes.size();
                 ++prime_index) {
                const silex::PrimeIdeal* prime = primes.at(prime_index);
                assert(prime != nullptr);
                silex::Ideal prime_ideal(maximal);
                sflint::FmpzMat hnf(n, n);
                assert(prime->get_ideal(prime_ideal));
                assert(prime_ideal.get_hnf(sflint::FmpzMatRef(hnf)));

                std::vector<silex::OrderElement> elements;
                for (slong t = 0; t < 12; ++t) {
                    sflint::FmpzMat coordinates(1, n);
                    bool nonzero = false;
                    for (slong j = 0; j < n; ++j) {
                        const slong c = ((7 * t + 3 * j + j * j * t) % 11) - 5;
                        fmpz_set_si(fmpz_mat_entry(coordinates.raw(), 0, j),
                                    c);
                        nonzero = nonzero || c != 0;
                    }
                    if (!nonzero) {
                        continue;
                    }
                    silex::OrderElement element(maximal);
                    assert(element.set_coordinates(
                            sflint::FmpzMatConstRef(coordinates)));
                    elements.push_back(std::move(element));
                }
                // Products of up to three HNF basis elements of P (each in
                // P), times a sample element, reach valuations above one.
                std::vector<silex::OrderElement> basis;
                for (slong i = 0; i < n; ++i) {
                    sflint::FmpzMat row(1, n);
                    for (slong j = 0; j < n; ++j) {
                        fmpz_set(fmpz_mat_entry(row.raw(), 0, j),
                                 fmpz_mat_entry(hnf.raw(), i, j));
                    }
                    silex::OrderElement element(maximal);
                    assert(element.set_coordinates(
                            sflint::FmpzMatConstRef(row)));
                    basis.push_back(std::move(element));
                }
                const std::size_t samples = elements.size();
                for (std::size_t i = 0; i < basis.size(); ++i) {
                    silex::OrderElement product(maximal);
                    assert(product.set(basis[i]));
                    for (std::size_t k = 1; k <= 2; ++k) {
                        silex::OrderElement next(maximal);
                        assert(next.multiply(
                                product,
                                basis[(i + 3 * k) % basis.size()]));
                        product = std::move(next);
                        silex::OrderElement mixed(maximal);
                        assert(mixed.multiply(product,
                                              elements[(i + k) % samples]));
                        elements.push_back(std::move(mixed));
                    }
                }

                for (const silex::OrderElement& element : elements) {
                    if (element.equal_si(0)) {
                        continue;
                    }
                    for (slong bound : {0, 1, 2, 5, 40}) {
                        slong expected = -1;
                        slong value = -1;
                        assert(silex::detail::TestHookAccess::valuation_by_power_containment(
                                               expected, *prime, element,
                                               bound));
                        if (!silex::detail::TestHookAccess::anti_uniformizer_valuation(
                                            value, *prime, element, bound)) {
                            // Only P = pO has no anti-uniformizer condition
                            // (every HNF row vanishes mod p); the helper
                            // declines there and the valuation falls back to
                            // containment.
                            assert(prime->ramification_index() == 1);
                            assert(prime->residue_degree() == n);
                            value = expected;
                        }
                        assert(value == expected);
                        if (bound == 40) {
                            assert(expected < bound);
                            slong public_value = -1;
                            assert(prime->valuation(public_value, element,
                                                    &diagnostics));
                            assert(public_value == expected);
                            if (expected > 0 && index_prime) {
                                ++checked_index_positive;
                            }
                            if (expected > 0 &&
                                prime->ramification_index() > 1) {
                                ++checked_ramified_positive;
                            }
                            if (expected > 1) {
                                ++checked_above_one;
                            }
                        }
                    }
                }
            }
        }
    }
    assert(checked_index_positive > 0);
    assert(checked_ramified_positive > 0);
    assert(checked_above_one > 0);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The public valuation (not only the test hook) took the route.
    assert(route_counter.count > 0);
#endif
    return 0;
}

// Differential check for the principal-ideal valuation route.  For an ideal
// with a known generator a, PrimeIdeal::valuation(Ideal) reuses the p-adic
// valuation of the ideal norm [O : aO] = |N(a)| instead of recomputing N(a).
// Its value must agree with the element route (which computes N(a)), with the
// containment route on the same HNF stored without a generator, and with the
// bounded power-containment definition.  The orders are maximal orders of a
// quartic and a sextic field whose maximal order is not the equation order,
// and of x^3 - 2, whose equation order is maximal.
int test_principal_ideal_valuation_matches_element_valuation() {
    const std::vector<std::vector<slong>> polynomials = {
            {474, -4, 5, -2, 1},
            {7, 0, 0, 0, 0, 0, 1},
            {-2, 0, 0, 1},
    };
    const std::vector<std::vector<slong>> seeds = {
            {0, 1},       {-1, 1},     {2, 1, 1},  {3, 6},
            {9, 6, 1},    {-25, 0, 5}, {30, 20, 10}, {-7, 14},
            {1, 0, 3, 0, 3, 0, 1},    {33, 15},   {12}, {1, -2, 3, -2, 1},
    };
    const slong bound = 64;
    slong checked = 0;
    slong checked_positive = 0;
    for (const std::vector<slong>& coefficients : polynomials) {
        sflint::FmpqPoly polynomial;
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                           coefficients[i]);
        }
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        const std::size_t degree = coefficients.size() - 1;

        std::vector<silex::OrderElement> elements;
        for (const std::vector<slong>& seed : seeds) {
            if (seed.size() > degree) {
                continue;
            }
            silex::OrderElement element(maximal);
            assert(set_order_element_from_coefficients(element, field, seed,
                                                       1));
            if (!element.equal_si(0)) {
                elements.push_back(std::move(element));
            }
        }
        // Products give larger and mixed valuations.
        const std::size_t seed_count = elements.size();
        for (std::size_t i = 0; i < seed_count; ++i) {
            for (std::size_t j = i; j < seed_count; ++j) {
                silex::OrderElement product(maximal);
                assert(product.multiply(elements[i], elements[j]));
                elements.push_back(std::move(product));
            }
        }

        for (slong rational_prime : {2, 3, 5, 7, 11}) {
            sflint::Fmpz p;
            assert(set_fmpz_si(p, rational_prime));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(primes, maximal,
                                          sflint::FmpzConstRef(p)));
            for (const silex::OrderElement& element : elements) {
                silex::Ideal principal(maximal);
                assert(principal.set_principal(element));
                assert(principal.known_principal_generator() != nullptr);
                sflint::FmpzMat hnf(static_cast<slong>(degree),
                                    static_cast<slong>(degree));
                silex::Ideal plain(maximal);
                assert(principal.get_hnf(sflint::FmpzMatRef(hnf)));
                assert(plain.set_hnf(sflint::FmpzMatConstRef(hnf)));
                assert(plain.known_principal_generator() == nullptr);

                for (slong index = 0; index < primes.size(); ++index) {
                    const silex::PrimeIdeal* prime = primes.at(index);
                    assert(prime != nullptr);
                    slong by_element = -1;
                    slong by_principal = -1;
                    slong by_plain = -1;
                    slong by_definition = -1;
                    assert(prime->valuation(by_element, element));
                    assert(prime->valuation(by_principal, principal));
                    assert(prime->valuation(by_plain, plain));
                    assert(silex::detail::TestHookAccess::valuation_by_power_containment(
                                           by_definition, *prime, element,
                                           bound));
                    assert(by_definition < bound);
                    assert(by_principal == by_element);
                    assert(by_principal == by_plain);
                    assert(by_principal == by_definition);
                    ++checked;
                    if (by_principal > 0) {
                        ++checked_positive;
                    }
                }
            }
        }
    }
    assert(checked > 0 && checked_positive > 0);
    return 0;
}

int test_valuation_fractional_element_and_factored() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(prime->ramification_index() == 2);

    silex::Element theta(field);
    silex::Element two(field);
    silex::Element half(field);
    silex::Element theta_over_two(field);
    assert(theta.gen());
    assert(two.set_si(2));
    assert(half.set_si_over_si(1, 2));
    assert(theta_over_two.scalar_div_si(theta, 2));

    slong value = 99;
    assert(prime->valuation(value, theta));
    assert(value == 1);
    assert(prime->valuation(value, two));
    assert(value == 2);
    assert(prime->valuation(value, half));
    assert(value == -2);
    assert(prime->valuation(value, theta_over_two));
    assert(value == -1);

    silex::Ideal unit_ideal(order);
    silex::Ideal prime_ideal(order);
    assert(unit_ideal.one());
    assert(prime->get_ideal(prime_ideal));

    sflint::Fmpz denominator;
    assert(set_fmpz_si(denominator, 2));
    silex::FractionalIdeal half_ideal(order);
    silex::FractionalIdeal prime_over_two(order);
    silex::FractionalIdeal principal_theta_over_two(order);
    assert(half_ideal.set_integral_den(unit_ideal,
                                       sflint::FmpzConstRef(denominator)));
    assert(prime_over_two.set_integral_den(
            prime_ideal, sflint::FmpzConstRef(denominator)));
    assert(principal_theta_over_two.set_principal(theta_over_two));

    assert(prime->valuation(value, half_ideal));
    assert(value == -2);
    assert(prime->valuation(value, prime_over_two));
    assert(value == -1);
    assert(prime->valuation(value, principal_theta_over_two));
    assert(value == -1);

    silex::FactoredElement factored(field);
    assert(factored.one());
    assert(factored.push(theta, 3));
    assert(factored.push(two, -2));
    assert(prime->valuation(value, factored));
    assert(value == -1);

    silex::Element expanded(field);
    slong expanded_value = 99;
    assert(factored.evaluate(expanded));
    assert(prime->valuation(expanded_value, expanded));
    assert(expanded_value == value);

    silex::FactoredElement one(field);
    assert(one.one());
    assert(prime->valuation(value, one));
    assert(value == 0);

    return 0;
}

int test_extended_valuation_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::Element zero(field);
    assert(zero.zero());
    slong value = 77;
    assert(!prime->valuation(value, zero));
    assert(value == 77);

    silex::FractionalIdeal undefined_fractional;
    assert(!prime->valuation(value, undefined_fractional));
    assert(value == 77);

    silex::FactoredElement undefined_factored;
    assert(!prime->valuation(value, undefined_factored));
    assert(value == 77);

    sflint::FmpqPoly other_polynomial;
    poly_x2_minus(other_polynomial, 5);
    silex::NumberField other_field;
    silex::Order other_equation =
            order_by_polynomial(other_field, other_polynomial);
    silex::Order other_order(other_field);
    assert(other_order.maximal_order(other_equation));

    silex::Element other_element(other_field);
    assert(other_element.set_si(2));
    assert(!prime->valuation(value, other_element));
    assert(value == 77);

    silex::FractionalIdeal other_fractional(other_order);
    assert(other_fractional.one());
    assert(!prime->valuation(value, other_fractional));
    assert(value == 77);

    silex::FactoredElement other_factored(other_field);
    assert(other_factored.one());
    assert(other_factored.push(other_element, 2));
    assert(!prime->valuation(value, other_factored));
    assert(value == 77);

    silex::FactoredElement one(field);
    assert(one.one());
    assert(!one.push(zero, 1));
    assert(prime->valuation(value, one));
    assert(value == 0);

    silex::Element two(field);
    silex::FactoredElement overflow(field);
    assert(two.set_si(2));
    assert(overflow.one());
    assert(overflow.push(two, WORD_MAX));
    value = 77;
    assert(!prime->valuation(value, overflow));
    assert(value == 77);

    return 0;
}

int test_valuation_maximal_quadratic() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 3));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::Ideal ideal(order);
    slong value = -1;
    assert(set_rational_principal(ideal, 9));
    assert(prime->valuation(value, ideal));
    assert(value == 2);

    value = -1;
    assert(set_rational_principal(ideal, 2));
    assert(prime->valuation(value, ideal));
    assert(value == 0);

    return 0;
}

int test_valuation_unmarked_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    assert(!order.is_maximal());

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 3));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 3));
    slong value = 99;
    assert(!prime->valuation(value, ideal));
    assert(value == 99);

    return 0;
}

int test_cubic_mixed_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 5));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    slong ef_sum = 0;
    bool saw_degree_one = false;
    bool saw_degree_two = false;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        assert(prime->ramification_index() == 1);
        ef_sum += prime->ramification_index() * prime->residue_degree();
        saw_degree_one = saw_degree_one || prime->residue_degree() == 1;
        saw_degree_two = saw_degree_two || prime->residue_degree() == 2;
    }
    assert(ef_sum == 3);
    assert(saw_degree_one);
    assert(saw_degree_two);

    return 0;
}

int test_residue_degree_limited_decomposition() {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    order = silex::test::verified_maximal_order(order);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 5));
    silex::PrimeIdealList full;
    assert(silex::decompose_prime(full, order, sflint::FmpzConstRef(p)));
    assert(full.size() == 2);
    assert(count_with_ef(full, 1, 1) == 1);
    assert(count_with_ef(full, 1, 2) == 1);

    silex::PrimeIdealList degree_one;
    assert(silex::decompose_prime(degree_one, order, sflint::FmpzConstRef(p),
                                  1));
    assert(degree_one.is_defined());
    assert(degree_one.size() == 1);
    assert(count_with_ef(degree_one, 1, 1) == 1);
    assert(check_norms(degree_one));

    silex::PrimeIdealList full_via_zero_limit;
    assert(silex::decompose_prime(full_via_zero_limit, order,
                                  sflint::FmpzConstRef(p), 0));
    assert(full_via_zero_limit.size() == full.size());
    assert(check_ef_sum(full_via_zero_limit, 3));

    return 0;
}

int test_kummer_generator_coordinates_are_centered_and_copied() {
    sflint::FmpqPoly polynomial;
    poly_x4_minus_x_minus_1(polynomial);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() >= 1);

    const silex::PrimeIdeal* degree_one = nullptr;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        if (prime != nullptr && prime->residue_degree() == 1) {
            degree_one = prime;
            break;
        }
    }
    assert(degree_one != nullptr);

    sflint::FmpzMat generator(1, order.degree());
    assert(degree_one->kummer_generator_coordinates(sflint::FmpzMatRef(
            generator)));
    assert(fmpz_mat_entry_is_si(generator, 0, 0, -3));
    assert(fmpz_mat_entry_is_si(generator, 0, 1, 1));
    assert(fmpz_mat_entry_is_si(generator, 0, 2, 0));
    assert(fmpz_mat_entry_is_si(generator, 0, 3, 0));

    silex::PrimeIdeal copied(order);
    assert(copied.set(*degree_one));
    sflint::FmpzMat copied_generator(1, order.degree());
    assert(copied.kummer_generator_coordinates(sflint::FmpzMatRef(
            copied_generator)));
    assert(fmpz_mat_entry_is_si(copied_generator, 0, 0, -3));
    assert(fmpz_mat_entry_is_si(copied_generator, 0, 1, 1));

    silex::PrimeIdeal moved(std::move(copied));
    sflint::FmpzMat moved_generator(1, order.degree());
    assert(moved.kummer_generator_coordinates(sflint::FmpzMatRef(
            moved_generator)));
    assert(fmpz_mat_entry_is_si(moved_generator, 0, 0, -3));
    assert(fmpz_mat_entry_is_si(moved_generator, 0, 1, 1));
    assert(!copied.is_defined());

    return 0;
}

int test_uncertified_repeated_factor_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);
    assert(!order.is_maximal());

    silex::PrimeIdealList primes(order, 1);
    assert(primes.is_defined());
    assert(primes.size() == 1);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 2));
    assert(!silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(prime->is_defined());
    assert(!prime->has_prime_data());

    return 0;
}

int test_residue_ring_prime_relation() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    silex::Ideal modulus(order);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);
    assert(prime->get_ideal(modulus));
    silex::ResidueRing ring(modulus);
    silex::ResidueElement theta(ring);
    silex::ResidueElement theta_squared(ring);
    silex::ResidueElement five(ring);
    silex::Element x(field);

    sflint::Fmpz cardinality;
    assert(ring.cardinality(sflint::FmpzRef(cardinality)));
    assert(sflint::fmpz_equal_si(cardinality, 11));

    assert(x.gen());
    assert(theta.set_element(x));
    assert(theta_squared.multiply(theta, theta));
    assert(x.set_si(5));
    assert(five.set_element(x));
    assert(theta_squared.equal(five));

    return 0;
}

int test_move_clear_and_redefine() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.is_defined());
    assert(primes.size() == 1);
    const silex::PrimeIdeal* prime = primes.at(0);
    assert(prime != nullptr);

    silex::PrimeIdeal copied(order);
    assert(copied.set(*prime));
    silex::PrimeIdeal moved(std::move(copied));
    assert(moved.is_defined());
    assert(moved.has_prime_data());
    assert(silex::same_order_parent(moved.parent(), &order));
    assert(moved.degree() == 1);
    assert(!copied.is_defined());
    assert(copied.parent() == nullptr);
    assert(copied.degree() == 0);

    silex::PrimeIdeal assigned;
    assigned = std::move(moved);
    assert(assigned.is_defined());
    assert(assigned.has_prime_data());
    assert(silex::same_order_parent(assigned.parent(), &order));
    assert(!moved.is_defined());
    assigned.clear();
    assert(!assigned.is_defined());
    assert(assigned.parent() == nullptr);
    assert(!assigned.has_prime_data());
    assert(assigned.define(order));
    assert(assigned.is_defined());
    assert(!assigned.has_prime_data());

    silex::PrimeIdealList list(order, 2);
    assert(list.is_defined());
    assert(list.size() == 2);
    silex::PrimeIdealList list_moved(std::move(list));
    assert(list_moved.is_defined());
    assert(list_moved.size() == 2);
    assert(!list.is_defined());
    assert(list.size() == 0);

    silex::PrimeIdealList list_assigned;
    list_assigned = std::move(list_moved);
    assert(list_assigned.is_defined());
    assert(list_assigned.size() == 2);
    assert(!list_moved.is_defined());
    list_assigned.clear();
    assert(!list_assigned.is_defined());
    assert(list_assigned.size() == 0);
    assert(list_assigned.define(order, 0));
    assert(list_assigned.is_defined());
    assert(list_assigned.size() == 0);

    return 0;
}

silex::PrimeIdeal local_prime_ideal() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);
    silex::NumberField field;
    silex::Order order;
    order = order_by_polynomial(field, polynomial);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    silex::PrimeIdeal prime(order);
    assert(prime.set(*primes.at(0)));
    return prime;
}

int test_keeps_parent_order_alive() {
    silex::PrimeIdeal prime = local_prime_ideal();
    assert(prime.is_defined());
    assert(prime.has_prime_data());
    assert(prime.parent() != nullptr);
    assert(prime.parent()->parent() != nullptr);
    assert(prime.degree() == 1);

    sflint::Fmpz p;
    sflint::Fmpz norm;
    assert(prime.rational_prime(sflint::FmpzRef(p)));
    assert(sflint::fmpz_equal_si(p, 11));
    assert(prime.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, 11));

    silex::Ideal ideal(*prime.parent());
    assert(prime.get_ideal(ideal));
    assert(ideal.parent() != nullptr);
    assert(silex::same_order_parent(ideal.parent(), prime.parent()));
    return 0;
}

// Residue maps on maximal orders whose basis is not a power basis.  The
// generic backend stores residue polynomials in the field generator alpha,
// so order coordinates must be mapped back through the order basis before
// reducing modulo (p, residue polynomial).

void element_alpha_over(silex::Element& out, slong numerator_constant,
                        slong numerator_alpha, slong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, numerator_constant);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, numerator_alpha);
    fmpq_poly_scalar_div_si(polynomial.raw(), polynomial.raw(), denominator);
    assert(out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool reduced_is_constant(const sflint::FmpzPoly& polynomial,
                         slong value) noexcept {
    return fmpz_poly_degree(polynomial.raw()) <= 0 &&
           fmpz_poly_coeff_is_si(polynomial, 0, value);
}

int test_generic_nonpower_basis_reduce_quadratic() {
    // x^2 - 12 is not on the quadratic backend (12 is not squarefree).  Its
    // maximal order has basis [1, alpha/2], index 2 over Z[alpha].
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 12);

    silex::NumberField field;
    silex::Order equation;
    equation = order_by_polynomial(field, polynomial);
    assert(field.backend_kind() == silex::NumberFieldBackendKind::generic);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    assert(!maximal.is_equation_order());
    sflint::Fmpz index;
    assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
    assert(sflint::fmpz_equal_si(index, 2));

    silex::Element alpha(field);
    silex::Element half_alpha(field);
    silex::Element third_alpha(field);
    assert(alpha.gen());
    element_alpha_over(half_alpha, 0, 1, 2);
    element_alpha_over(third_alpha, 0, 1, 3);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    bool saw_plus_one = false;
    bool saw_minus_one = false;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        assert(prime->residue_degree() == 1);

        sflint::FmpzPoly residue_polynomial;
        assert(prime->residue_polynomial(
                sflint::FmpzPolyRef(residue_polynomial)));
        assert(fmpz_poly_degree(residue_polynomial.raw()) == 1);

        sflint::FmpzPoly reduced_alpha;
        sflint::FmpzPoly reduced_half;
        assert(prime->reduce(sflint::FmpzPolyRef(reduced_alpha), alpha));
        assert(prime->reduce(sflint::FmpzPolyRef(reduced_half), half_alpha));

        // alpha/2 is the second order basis vector: coordinates (0, 1).
        silex::OrderElement half_order(maximal);
        sflint::FmpzMat coordinates(1, 2);
        fmpz_one(fmpz_mat_entry(coordinates.raw(), 0, 1));
        assert(half_order.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
        sflint::FmpzPoly reduced_half_order;
        assert(prime->reduce(sflint::FmpzPolyRef(reduced_half_order),
                             half_order));
        assert(fmpz_poly_equal(reduced_half_order.raw(),
                               reduced_half.raw()) != 0);

        if (reduced_is_constant(reduced_alpha, 1)) {
            // P = (11, alpha - 1): alpha/2 -> 1/2 = 6 mod 11.
            assert(fmpz_poly_coeff_is_si(residue_polynomial, 1, 1));
            assert(fmpz_poly_coeff_is_si(residue_polynomial, 0, 10));
            assert(reduced_is_constant(reduced_half, 6));
            saw_plus_one = true;
        } else {
            // P = (11, alpha + 1): alpha -> 10, alpha/2 -> -1/2 = 5.
            assert(reduced_is_constant(reduced_alpha, 10));
            assert(reduced_is_constant(reduced_half, 5));
            saw_minus_one = true;
        }

        // Only integral elements reduce through PrimeIdeal::reduce.
        sflint::FmpzPoly untouched;
        fmpz_poly_set_si(untouched.raw(), 77);
        assert(!prime->reduce(sflint::FmpzPolyRef(untouched), third_alpha));
        assert(reduced_is_constant(untouched, 77));
    }
    assert(saw_plus_one && saw_minus_one);

    // Primes dividing the index carry no alpha residue polynomial, so the
    // alpha-polynomial reduction is unavailable there.
    assert(set_fmpz_si(p, 2));
    silex::PrimeIdealList primes2;
    assert(silex::decompose_prime(primes2, maximal, sflint::FmpzConstRef(p)));
    for (slong i = 0; i < primes2.size(); ++i) {
        sflint::FmpzPoly reduced;
        assert(!primes2.at(i)->reduce(sflint::FmpzPolyRef(reduced), alpha));
    }

    return 0;
}

void poly_mulmod_prime(sflint::FmpzPoly& out,
                       const sflint::FmpzPoly& left,
                       const sflint::FmpzPoly& right,
                       const sflint::FmpzPoly& modulus,
                       const sflint::Fmpz& p) noexcept {
    sflint::FmpzPoly product;
    fmpz_poly_mul(product.raw(), left.raw(), right.raw());
    fmpz_poly_rem(out.raw(), product.raw(), modulus.raw());
    fmpz_poly_scalar_mod_fmpz(out.raw(), out.raw(), p.raw());
}

void poly_addmod_prime(sflint::FmpzPoly& out,
                       const sflint::FmpzPoly& left,
                       const sflint::FmpzPoly& right,
                       const sflint::Fmpz& p) noexcept {
    fmpz_poly_add(out.raw(), left.raw(), right.raw());
    fmpz_poly_scalar_mod_fmpz(out.raw(), out.raw(), p.raw());
}

// Checks that PrimeIdeal::reduce is a ring homomorphism O -> O/P on a
// maximal order with index > 1, at every prime above each listed p.
// Returns the number of primes checked.
slong check_reduce_homomorphism(const silex::Order& maximal,
                                const slong* rational_primes,
                                slong prime_count) noexcept {
    const slong n = maximal.degree();
    static const slong coordinate_rows[][5] = {
            {1, 0, 0, 0, 0},  {0, 1, 0, 0, 0},  {0, 0, 1, 0, 0},
            {0, 0, 0, 1, 0},  {2, -1, 3, 1, 0}, {-3, 2, 0, 5, 0},
            {1, 4, -2, -1, 0}, {0, 3, 1, 2, 0},
    };
    std::vector<silex::OrderElement> elements;
    for (const auto& row : coordinate_rows) {
        sflint::FmpzMat coordinates(1, n);
        bool nonzero = false;
        for (slong j = 0; j < n; ++j) {
            fmpz_set_si(fmpz_mat_entry(coordinates.raw(), 0, j), row[j]);
            nonzero = nonzero || row[j] != 0;
        }
        if (!nonzero) {
            continue;
        }
        elements.emplace_back(maximal);
        assert(elements.back().set_coordinates(
                sflint::FmpzMatConstRef(coordinates)));
    }

    slong checked = 0;
    for (slong k = 0; k < prime_count; ++k) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, rational_primes[k]));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, maximal,
                                      sflint::FmpzConstRef(p)));
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(prime != nullptr);
            sflint::FmpzPoly modulus;
            assert(prime->residue_polynomial(sflint::FmpzPolyRef(modulus)));

            silex::OrderElement one(maximal);
            sflint::FmpzPoly reduced_one;
            assert(one.one());
            assert(prime->reduce(sflint::FmpzPolyRef(reduced_one), one));
            assert(reduced_is_constant(reduced_one, 1));

            // p is in P, and a generator of P reduces to zero.
            silex::OrderElement p_element(maximal);
            sflint::FmpzPoly reduced_p;
            assert(p_element.set_si(rational_primes[k]));
            assert(prime->reduce(sflint::FmpzPolyRef(reduced_p), p_element));
            assert(fmpz_poly_is_zero(reduced_p.raw()) != 0);

            for (const silex::OrderElement& left : elements) {
                sflint::FmpzPoly reduced_left;
                assert(prime->reduce(sflint::FmpzPolyRef(reduced_left), left));
                silex::Element left_element(*maximal.parent());
                sflint::FmpzPoly reduced_left_element;
                assert(left.get_element(left_element));
                assert(prime->reduce(sflint::FmpzPolyRef(reduced_left_element),
                                     left_element));
                assert(fmpz_poly_equal(reduced_left.raw(),
                                       reduced_left_element.raw()) != 0);
                for (const silex::OrderElement& right : elements) {
                    silex::OrderElement product(maximal);
                    silex::OrderElement sum(maximal);
                    assert(product.multiply(left, right));
                    assert(sum.add(left, right));

                    sflint::FmpzPoly reduced_right;
                    sflint::FmpzPoly reduced_product;
                    sflint::FmpzPoly reduced_sum;
                    sflint::FmpzPoly expected;
                    assert(prime->reduce(sflint::FmpzPolyRef(reduced_right),
                                         right));
                    assert(prime->reduce(
                            sflint::FmpzPolyRef(reduced_product), product));
                    assert(prime->reduce(sflint::FmpzPolyRef(reduced_sum),
                                         sum));

                    poly_mulmod_prime(expected, reduced_left, reduced_right,
                                      modulus, p);
                    assert(fmpz_poly_equal(expected.raw(),
                                           reduced_product.raw()) != 0);
                    poly_addmod_prime(expected, reduced_left, reduced_right,
                                      p);
                    assert(fmpz_poly_equal(expected.raw(),
                                           reduced_sum.raw()) != 0);
                }
            }
            ++checked;
        }
    }
    return checked;
}

int test_reduce_homomorphism_on_index_maximal_orders() {
    {
        // Cubic, [O : Z[alpha]] = 2.
        sflint::FmpqPoly polynomial;
        poly_cubic_disc1724(polynomial);
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        assert(field.backend_kind() ==
               silex::NumberFieldBackendKind::generic);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        sflint::Fmpz index;
        assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
        assert(sflint::fmpz_equal_si(index, 2));

        const slong primes[] = {3, 5, 7, 11, 13, 431};
        assert(check_reduce_homomorphism(maximal, primes, 6) >= 6);
    }
    {
        // Quartic, [O : Z[alpha]] = 3.
        sflint::FmpqPoly polynomial;
        poly_index3_quartic(polynomial);
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        sflint::Fmpz index;
        assert(silex::order_index(sflint::FmpzRef(index), equation, maximal));
        assert(sflint::fmpz_equal_si(index, 3));

        const slong primes[] = {2, 5, 7, 11, 13, 17};
        assert(check_reduce_homomorphism(maximal, primes, 6) >= 6);
    }
    return 0;
}

int test_quadratic_backend_omega_reduce_unchanged() {
    // The direct maximal-quadratic path keeps residue polynomials in the
    // integral generator omega = (1 + theta)/2, and reduce() returns
    // omega-polynomials; theta = 2 omega - 1.
    silex::NumberField field = silex::test::quadratic_field(5);
    assert(field.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::Element theta(field);
    silex::Element omega(field);
    assert(theta.gen());
    element_alpha_over(omega, 1, 1, 2);

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        sflint::Fmpz omega_root;
        assert(degree_one_root_from_residue_polynomial(omega_root, *prime));
        const slong root = static_cast<slong>(fmpz_get_si(omega_root.raw()));

        sflint::FmpzPoly reduced;
        assert(prime->reduce(sflint::FmpzPolyRef(reduced), omega));
        assert(reduced_is_constant(reduced, root));
        assert(prime->reduce(sflint::FmpzPolyRef(reduced), theta));
        assert(reduced_is_constant(reduced, ((2 * root - 1) % 11 + 11) % 11));
    }

    // A degree-two residue field: omega itself reduces to the variable x.
    assert(set_fmpz_si(p, 7));
    silex::PrimeIdealList inert;
    assert(silex::decompose_prime(inert, maximal, sflint::FmpzConstRef(p)));
    assert(inert.size() == 1);
    sflint::FmpzPoly reduced;
    assert(inert.at(0)->reduce(sflint::FmpzPolyRef(reduced), omega));
    assert(fmpz_poly_degree(reduced.raw()) == 1);
    assert(fmpz_poly_coeff_is_si(reduced, 1, 1));
    assert(fmpz_poly_coeff_is_si(reduced, 0, 0));
    return 0;
}

// The residue-polynomial variable (omega or alpha) is fixed when the prime is
// built.  Order handles share their maximality record, so withdrawing the
// parent order's maximality afterwards (only reachable through the internal
// OrderAccess hook) must not make reduce() read omega-polynomial data as
// alpha-polynomial data.
int test_residue_convention_fixed_at_decomposition() {
    silex::NumberField field = silex::test::quadratic_field(5);
    assert(field.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::Element theta(field);
    silex::Element omega(field);
    assert(theta.gen());
    element_alpha_over(omega, 1, 1, 2);
    silex::OrderElement omega_order(maximal);
    assert(omega_order.set_element(omega));

    for (const slong prime_value : {11L, 7L}) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, prime_value));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, maximal,
                                      sflint::FmpzConstRef(p)));
        std::vector<sflint::FmpzPoly> theta_images(
                static_cast<std::size_t>(primes.size()));
        std::vector<sflint::FmpzPoly> omega_images(
                static_cast<std::size_t>(primes.size()));
        for (slong i = 0; i < primes.size(); ++i) {
            const auto k = static_cast<std::size_t>(i);
            assert(primes.at(i)->reduce(
                    sflint::FmpzPolyRef(theta_images[k]), theta));
            assert(primes.at(i)->reduce(
                    sflint::FmpzPolyRef(omega_images[k]), omega));
        }

        {
            silex::detail::OrderAccess::set_maximality_unchecked(maximal,
                                                                 false);
            for (slong i = 0; i < primes.size(); ++i) {
                const auto k = static_cast<std::size_t>(i);
                sflint::FmpzPoly reduced;
                assert(primes.at(i)->reduce(sflint::FmpzPolyRef(reduced),
                                            theta));
                assert(fmpz_poly_equal(reduced.raw(),
                                       theta_images[k].raw()) != 0);
                assert(primes.at(i)->reduce(sflint::FmpzPolyRef(reduced),
                                            omega));
                assert(fmpz_poly_equal(reduced.raw(),
                                       omega_images[k].raw()) != 0);
                assert(primes.at(i)->reduce(sflint::FmpzPolyRef(reduced),
                                            omega_order));
                assert(fmpz_poly_equal(reduced.raw(),
                                       omega_images[k].raw()) != 0);

                silex::PrimeIdeal copied(maximal);
                assert(copied.set(*primes.at(i)));
                assert(copied.reduce(sflint::FmpzPolyRef(reduced), omega));
                assert(fmpz_poly_equal(reduced.raw(),
                                       omega_images[k].raw()) != 0);
            }
        }
        silex::detail::OrderAccess::set_maximality_unchecked(maximal, true);
    }
    return 0;
}

// The degree-one root helper converts an omega root to the theta root from
// the convention recorded in the prime, so withdrawing the parent order's
// maximality afterwards (internal OrderAccess hook only) keeps the theta
// root: theta = 2 omega - 1 must still satisfy theta^2 = 5 mod p.
int test_degree_one_root_follows_recorded_convention() {
    silex::NumberField field = silex::test::quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);
    std::vector<slong> roots;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(PrimeIdealAccess::residue_uses_integral_generator(*prime));
        sflint::Fmpz root;
        sflint::Fmpz root_p;
        assert(silex::detail::degree_one_prime_root_mod_p(root, root_p,
                                                          *prime));
        const slong value = static_cast<slong>(fmpz_get_si(root.raw()));
        assert((value * value - 5) % 11 == 0);
        roots.push_back(value);
    }

    silex::detail::OrderAccess::set_maximality_unchecked(maximal, false);
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(PrimeIdealAccess::residue_uses_integral_generator(*prime));
        sflint::Fmpz root;
        sflint::Fmpz root_p;
        assert(silex::detail::degree_one_prime_root_mod_p(root, root_p,
                                                          *prime));
        assert(fmpz_get_si(root.raw()) ==
               roots[static_cast<std::size_t>(i)]);
    }
    silex::detail::OrderAccess::set_maximality_unchecked(maximal, true);
    return 0;
}

// The word route factored_value_at_degree_one_root_nmod reads each factor's
// stored form num(x) / den in place and inverts once at the end.  It must
// fail exactly where the fmpz route factored_value_at_degree_one_root does
// (a factor vanishing at the root, or a coefficient denominator divisible by
// p) and otherwise give the same residue, on the linear, quadratic and
// generic FLINT representations, for positive, negative and zero exponents.
// The roots range over all of F_p: both routes only evaluate polynomials.
int test_factored_value_nmod_matches_fmpz_route() {
    std::vector<silex::NumberField> fields;
    {
        sflint::FmpqPoly linear;
        sflint::fmpq_poly_set_coeff_si(linear, 1, 1);
        sflint::fmpq_poly_set_coeff_si(linear, 0, -4);
        fields.push_back(silex::test::field_by_polynomial(
                sflint::FmpqPolyConstRef(linear)));
    }
    fields.push_back(silex::test::quadratic_field(-7));
    {
        sflint::FmpqPoly cubic;
        poly_x3_minus(cubic, 2);
        fields.push_back(silex::test::field_by_polynomial(
                sflint::FmpqPolyConstRef(cubic)));
    }

    // Factor coefficients (numerator, denominator), low to high, padded to
    // the field degree by dropping higher terms.
    const std::array<std::array<std::array<slong, 2>, 3>, 5> factor_data = {{
            {{{7, 6}, {-5, 6}, {3, 6}}},
            {{{1, 1}, {1, 1}, {0, 1}}},
            {{{2, 5}, {1, 3}, {-1, 7}}},
            {{{-9, 4}, {0, 1}, {11, 10}}},
            {{{1, 15}, {0, 1}, {0, 1}}},
    }};
    const std::array<std::array<slong, 3>, 5> exponent_data = {{
            {1, 0, 0},
            {-2, 3, 0},
            {1, -1, 2},
            {0, 5, -3},
            {-1, -1, 4},
    }};

    slong checked = 0;
    slong failures = 0;
    for (const silex::NumberField& field : fields) {
        const slong degree = field.degree();
        std::vector<silex::Element> factors;
        for (const auto& data : factor_data) {
            sflint::FmpqPoly polynomial;
            for (slong i = 0; i < degree && i < 3; ++i) {
                sflint::Fmpq coefficient;
                fmpq_set_si(coefficient.raw(), data[i][0],
                            static_cast<ulong>(data[i][1]));
                ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), i,
                                           coefficient.raw());
            }
            silex::Element element(field);
            assert(element.set_fmpq_poly(
                    sflint::FmpqPolyConstRef(polynomial)));
            factors.push_back(std::move(element));
        }

        for (const ulong p : {2UL, 3UL, 5UL, 7UL, 11UL, 13UL}) {
            sflint::Fmpz p_fmpz;
            sflint::fmpz_set_ui(sflint::FmpzRef(p_fmpz), p);
            for (ulong root = 0; root < p; ++root) {
                sflint::Fmpz root_fmpz;
                sflint::fmpz_set_ui(sflint::FmpzRef(root_fmpz), root);
                for (const auto& exponents : exponent_data) {
                    for (std::size_t shift = 0; shift < factors.size();
                         ++shift) {
                        silex::FactoredElement element(field);
                        for (std::size_t k = 0; k < exponents.size(); ++k) {
                            assert(element.push(
                                    factors[(shift + k) % factors.size()],
                                    exponents[k]));
                        }
                        sflint::Fmpz expected;
                        const bool expected_ok =
                                silex::detail::factored_value_at_degree_one_root(
                                        expected, element,
                                        sflint::FmpzConstRef(p_fmpz),
                                        sflint::FmpzConstRef(root_fmpz));
                        ulong value = 0;
                        const bool ok = silex::detail::
                                factored_value_at_degree_one_root_nmod(
                                        value, element, p, root);
                        assert(ok == expected_ok);
                        if (ok) {
                            assert(fmpz_cmp_ui(expected.raw(), value) == 0);
                        } else {
                            ++failures;
                        }
                        ++checked;
                    }
                }
            }
        }
    }
    // Both outcomes are exercised.
    assert(failures > 0 && failures < checked);
    return 0;
}

// Same comparison at a prime just below 2^64 and for factors whose numerator
// coefficients (positive and negative) and denominator exceed a word.
int test_factored_value_nmod_matches_fmpz_route_large() {
    const char* const prime_text = "18446744073709551557";
    const char* const numerators[3][3] = {
            {"123456789012345678901234567890", "-98765432109876543210987654321",
             "55555555555555555555555555"},
            {"-340282366920938463463374607431768211457", "7",
             "-170141183460469231731687303715884105727"},
            {"18446744073709551557", "1", "-3"},
    };
    const char* const denominators[3] = {"1000000000000000000000007",
                                         "36893488147419103232",
                                         "99999999999999999999"};
    sflint::FmpqPoly cubic;
    poly_x3_minus(cubic, 2);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(cubic));
    std::vector<silex::Element> factors;
    for (int f = 0; f < 3; ++f) {
        sflint::FmpqPoly polynomial;
        for (slong i = 0; i < 3; ++i) {
            sflint::Fmpq coefficient;
            fmpz_set_str(fmpq_numref(coefficient.raw()), numerators[f][i],
                         10);
            fmpz_set_str(fmpq_denref(coefficient.raw()), denominators[f],
                         10);
            fmpq_canonicalise(coefficient.raw());
            ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), i,
                                       coefficient.raw());
        }
        silex::Element element(field);
        assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        factors.push_back(std::move(element));
    }

    sflint::Fmpz p_fmpz;
    fmpz_set_str(p_fmpz.raw(), prime_text, 10);
    const ulong p = fmpz_get_ui(p_fmpz.raw());
    const ulong roots[] = {0UL, 1UL, 2UL, 123456789UL, p / 2, p - 2, p - 1};
    const slong exponent_data[3][3] = {{1, 1, 1}, {2, -1, 3}, {-3, 2, -1}};
    slong checked = 0;
    for (const ulong root : roots) {
        sflint::Fmpz root_fmpz;
        sflint::fmpz_set_ui(sflint::FmpzRef(root_fmpz), root);
        for (const auto& exponents : exponent_data) {
            silex::FactoredElement element(field);
            for (std::size_t k = 0; k < 3; ++k) {
                assert(element.push(factors[k], exponents[k]));
            }
            sflint::Fmpz expected;
            const bool expected_ok =
                    silex::detail::factored_value_at_degree_one_root(
                            expected, element, sflint::FmpzConstRef(p_fmpz),
                            sflint::FmpzConstRef(root_fmpz));
            ulong value = 0;
            const bool ok = silex::detail::
                    factored_value_at_degree_one_root_nmod(value, element, p,
                                                           root);
            assert(ok == expected_ok);
            if (ok) {
                assert(fmpz_cmp_ui(expected.raw(), value) == 0);
            }
            ++checked;
        }
    }
    assert(checked == 21);
    return 0;
}

// Reverse misuse direction: an alpha-convention prime keeps that convention
// after its order is marked maximal, for reduce() and valuation (which
// needs the maximal marking, so it is checked afterwards).  For d = 3
// the equation order Z[theta] is maximal with conductor 1 and is recorded
// as maximal, so correct use gives the integral-generator convention.  The
// internal OrderAccess hook withdraws maximality to build an alpha prime and
// then restores it.  A prime decomposed after the restore uses the integral
// generator convention; PrimeIdeal::equal tells the two apart even though
// their ideals and residue-polynomial coefficients agree.
int test_alpha_convention_kept_after_marking_maximal() {
    silex::NumberField field = silex::test::quadratic_field(3);
    assert(field.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    silex::Order equation = silex::test::equation_order(field);
    assert(equation.is_maximal());

    silex::Element theta(field);
    assert(theta.gen());
    silex::Element shifted(field);
    silex::OrderElement shifted_order(equation);

    for (const slong prime_value : {11L, 5L}) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, prime_value));
        silex::detail::OrderAccess::set_maximality_unchecked(equation, false);
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, equation,
                                      sflint::FmpzConstRef(p)));
        const auto count = static_cast<std::size_t>(primes.size());
        std::vector<sflint::FmpzPoly> theta_images(count);
        for (slong i = 0; i < primes.size(); ++i) {
            const auto k = static_cast<std::size_t>(i);
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(!PrimeIdealAccess::residue_uses_integral_generator(*prime));
            assert(prime->reduce(sflint::FmpzPolyRef(theta_images[k]),
                                 theta));
        }

        silex::detail::OrderAccess::set_maximality_unchecked(equation, true);
        silex::PrimeIdealList marked;
        assert(silex::decompose_prime(marked, equation,
                                      sflint::FmpzConstRef(p)));
        assert(marked.size() == primes.size());
        for (slong i = 0; i < primes.size(); ++i) {
            const auto k = static_cast<std::size_t>(i);
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(!PrimeIdealAccess::residue_uses_integral_generator(*prime));
            sflint::FmpzPoly reduced;
            assert(prime->reduce(sflint::FmpzPolyRef(reduced), theta));
            assert(fmpz_poly_equal(reduced.raw(),
                                   theta_images[k].raw()) != 0);
            // theta - r lies in P (valuation 1, norm r^2 - 3) for the
            // degree-one prime with root r; for the inert prime theta + 1
            // has norm -2, a unit at P.
            slong root = 1;
            if (prime->residue_degree() == 1) {
                sflint::Fmpz stored_root;
                assert(degree_one_root_from_residue_polynomial(stored_root,
                                                               *prime));
                root = static_cast<slong>(fmpz_get_si(stored_root.raw()));
            }
            element_alpha_over(shifted, -root, 1, 1);
            assert(shifted_order.set_element(shifted));
            slong value = -1;
            assert(prime->valuation(value, shifted_order));
            assert(value == (prime->residue_degree() == 1 ? 1 : 0));

            bool found_same_ideal = false;
            for (slong j = 0; j < marked.size(); ++j) {
                const silex::PrimeIdeal* fresh = marked.at(j);
                assert(PrimeIdealAccess::residue_uses_integral_generator(
                        *fresh));
                sflint::FmpzPoly prime_poly;
                sflint::FmpzPoly fresh_poly;
                assert(prime->residue_polynomial(
                        sflint::FmpzPolyRef(prime_poly)));
                assert(fresh->residue_polynomial(
                        sflint::FmpzPolyRef(fresh_poly)));
                if (fmpz_poly_equal(prime_poly.raw(), fresh_poly.raw()) !=
                    0) {
                    silex::Ideal prime_ideal(equation);
                    silex::Ideal fresh_ideal(equation);
                    assert(prime->get_ideal(prime_ideal));
                    assert(fresh->get_ideal(fresh_ideal));
                    assert(prime_ideal.equal(fresh_ideal));
                    found_same_ideal = true;
                    slong fresh_value = -1;
                    assert(fresh->valuation(fresh_value, shifted_order));
                    assert(fresh_value == value);
                    assert(!prime->equal(*fresh));
                    assert(!fresh->equal(*prime));
                }
                assert(fresh->equal(*fresh));
            }
            assert(found_same_ideal);
            assert(prime->equal(*prime));
        }
    }
    return 0;
}

// For every prime P above each listed p that carries residue data, the
// stored Kummer generator lies in P (reduces to 0) and the images of the
// order basis span O/P = F_p^f, so reduce() is onto the residue field.
// Returns the number of primes checked.
slong check_kummer_generator_and_residue_span(const silex::Order& order,
                                              const slong* rational_primes,
                                              slong prime_count) noexcept {
    const slong n = order.degree();
    slong checked = 0;
    for (slong k = 0; k < prime_count; ++k) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, rational_primes[k]));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, order,
                                      sflint::FmpzConstRef(p)));
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(prime != nullptr);
            sflint::FmpzPoly modulus;
            if (!prime->residue_polynomial(sflint::FmpzPolyRef(modulus))) {
                continue;
            }
            const slong f = prime->residue_degree();
            assert(fmpz_poly_degree(modulus.raw()) == f);

            sflint::FmpzMat kummer(1, n);
            assert(prime->kummer_generator_coordinates(
                    sflint::FmpzMatRef(kummer)));
            silex::OrderElement generator(order);
            assert(generator.set_coordinates(sflint::FmpzMatConstRef(kummer)));
            sflint::FmpzPoly reduced;
            assert(prime->reduce(sflint::FmpzPolyRef(reduced), generator));
            assert(fmpz_poly_is_zero(reduced.raw()) != 0);

            sflint::NmodMat images(n, f, static_cast<ulong>(
                                                 rational_primes[k]));
            for (slong j = 0; j < n; ++j) {
                sflint::FmpzMat unit(1, n);
                fmpz_one(fmpz_mat_entry(unit.raw(), 0, j));
                silex::OrderElement basis_element(order);
                assert(basis_element.set_coordinates(
                        sflint::FmpzMatConstRef(unit)));
                assert(prime->reduce(sflint::FmpzPolyRef(reduced),
                                     basis_element));
                assert(fmpz_poly_degree(reduced.raw()) < f);
                for (slong c = 0; c < f; ++c) {
                    sflint::Fmpz coefficient;
                    fmpz_poly_get_coeff_fmpz(coefficient.raw(),
                                             reduced.raw(), c);
                    nmod_mat_entry(images.raw(), j, c) =
                            fmpz_fdiv_ui(coefficient.raw(),
                                         static_cast<ulong>(
                                                 rational_primes[k]));
                }
            }
            assert(nmod_mat_rank(images.raw()) == f);
            ++checked;
        }
    }
    return checked;
}

int test_kummer_generator_reduces_to_zero_and_residues_span() {
    {
        // Direct maximal-quadratic omega convention: split, inert, ramified.
        silex::NumberField field = silex::test::quadratic_field(5);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        const slong primes[] = {2, 5, 7, 11};
        assert(check_kummer_generator_and_residue_span(maximal, primes, 4) >=
               5);
    }
    {
        // Equation order, alpha convention.
        sflint::FmpqPoly polynomial;
        poly_x4_minus_x_minus_1(polynomial);
        silex::NumberField field;
        silex::Order order;
        order = order_by_polynomial(field, polynomial);
        const slong primes[] = {3, 7, 11};
        assert(check_kummer_generator_and_residue_span(order, primes, 3) >=
               3);
    }
    {
        // Maximal orders with [O : Z[alpha]] = 2 and 3, at primes not
        // dividing the index.
        sflint::FmpqPoly cubic;
        poly_cubic_disc1724(cubic);
        silex::NumberField cubic_field;
        silex::Order cubic_equation;
        cubic_equation = order_by_polynomial(cubic_field, cubic);
        silex::Order cubic_maximal(cubic_field);
        assert(cubic_maximal.maximal_order(cubic_equation));
        const slong cubic_primes[] = {3, 5, 7, 13};
        assert(check_kummer_generator_and_residue_span(cubic_maximal,
                                                       cubic_primes, 4) >= 4);

        sflint::FmpqPoly quartic;
        poly_index3_quartic(quartic);
        silex::NumberField quartic_field;
        silex::Order quartic_equation;
        quartic_equation = order_by_polynomial(quartic_field, quartic);
        silex::Order quartic_maximal(quartic_field);
        assert(quartic_maximal.maximal_order(quartic_equation));
        const slong quartic_primes[] = {2, 5, 7};
        assert(check_kummer_generator_and_residue_span(quartic_maximal,
                                                       quartic_primes, 3) >=
               3);
    }
    return 0;
}

// The prime-ideal inverse route used by the relation search must return
// exactly the normalized numerator and denominator of
// FractionalIdeal::invert.  Returns the number of ideals compared and counts
// prime ideals without a stored Kummer generator.
slong check_prime_inverse_matches_invert(const silex::Order& maximal,
                                         slong prime_bound,
                                         slong& without_kummer) {
    const slong n = maximal.degree();
    slong compared = 0;
    sflint::FmpzMat generator(1, n);
    for (slong value = 2; value < prime_bound; ++value) {
        if (n_is_prime(static_cast<ulong>(value)) == 0) {
            continue;
        }
        sflint::Fmpz p;
        assert(set_fmpz_si(p, value));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, maximal,
                                      sflint::FmpzConstRef(p)));
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal& prime = *primes.at(i);
            if (!prime.kummer_generator_coordinates(
                        sflint::FmpzMatRef(generator))) {
                ++without_kummer;
            }
            silex::Ideal ideal(maximal);
            assert(prime.get_ideal(ideal));

            silex::FractionalIdeal fractional(maximal);
            silex::FractionalIdeal inverse(maximal);
            silex::Ideal old_numerator(maximal);
            sflint::Fmpz old_den;
            assert(fractional.set_integral(ideal));
            assert(inverse.invert(fractional));
            assert(inverse.get_integral_den(old_numerator,
                                            sflint::FmpzRef(old_den)));

            silex::Ideal new_numerator(maximal);
            sflint::Fmpz new_den;
            assert(silex::detail::prime_ideal_inverse_num_den(
                    new_numerator, sflint::FmpzRef(new_den), ideal));

            sflint::FmpzMat old_hnf(n, n);
            sflint::FmpzMat new_hnf(n, n);
            assert(old_numerator.get_hnf(sflint::FmpzMatRef(old_hnf)));
            assert(new_numerator.get_hnf(sflint::FmpzMatRef(new_hnf)));
            assert(fmpz_mat_equal(old_hnf.raw(), new_hnf.raw()) != 0);
            assert(sflint::fmpz_equal(sflint::FmpzConstRef(old_den),
                                      sflint::FmpzConstRef(new_den)));
            assert(sflint::fmpz_equal(sflint::FmpzConstRef(new_den),
                                      sflint::FmpzConstRef(p)));
            ++compared;
        }

        // A product of two distinct primes above p also contains p O.  The
        // route may decline it (one kernel vector need not generate), but
        // whatever it accepts is exact.
        if (primes.size() >= 2) {
            silex::Ideal first(maximal);
            silex::Ideal second(maximal);
            silex::Ideal product(maximal);
            assert(primes.at(0)->get_ideal(first));
            assert(primes.at(1)->get_ideal(second));
            assert(product.multiply(first, second));
            silex::FractionalIdeal fractional(maximal);
            silex::FractionalIdeal inverse(maximal);
            silex::Ideal old_numerator(maximal);
            sflint::Fmpz old_den;
            assert(fractional.set_integral(product));
            assert(inverse.invert(fractional));
            assert(inverse.get_integral_den(old_numerator,
                                            sflint::FmpzRef(old_den)));
            silex::Ideal new_numerator(maximal);
            sflint::Fmpz new_den;
            if (!silex::detail::prime_ideal_inverse_num_den(
                        new_numerator, sflint::FmpzRef(new_den), product)) {
                continue;
            }
            sflint::FmpzMat old_hnf(n, n);
            sflint::FmpzMat new_hnf(n, n);
            assert(old_numerator.get_hnf(sflint::FmpzMatRef(old_hnf)));
            assert(new_numerator.get_hnf(sflint::FmpzMatRef(new_hnf)));
            assert(fmpz_mat_equal(old_hnf.raw(), new_hnf.raw()) != 0);
            assert(sflint::fmpz_equal(sflint::FmpzConstRef(old_den),
                                      sflint::FmpzConstRef(new_den)));
            ++compared;
        }
    }
    return compared;
}

silex::Order cyclotomic_maximal_order(silex::NumberField& field,
                                      ulong conductor) {
    sflint::FmpzPoly cyclotomic;
    fmpz_poly_cyclotomic(cyclotomic.raw(), conductor);
    sflint::FmpqPoly polynomial;
    fmpq_poly_set_fmpz_poly(polynomial.raw(), cyclotomic.raw());
    silex::Order equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    return maximal;
}

int test_prime_ideal_inverse_matches_invert() {
    slong without_kummer = 0;
    {
        // Search-base rows: cyclotomic_64_phi32 and cyclotomic_37_phi36.
        silex::NumberField field;
        silex::Order maximal = cyclotomic_maximal_order(field, 64);
        assert(check_prime_inverse_matches_invert(maximal, 30,
                                                  without_kummer) >= 8);
    }
    {
        silex::NumberField field;
        silex::Order maximal = cyclotomic_maximal_order(field, 37);
        assert(check_prime_inverse_matches_invert(maximal, 40,
                                                  without_kummer) >= 8);
    }

    // Index-divisor fields: x^3 - x^2 - 2x - 8 (2 splits completely and is a
    // common index divisor), cubic_disc1724 and the index-3 quartic.  Their
    // primes above the index divisors carry no Kummer generator.
    const slong before_index_fields = without_kummer;
    {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_zero(polynomial);
        sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 2, -1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -2);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, -8);
        silex::NumberField field;
        silex::Order equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(check_prime_inverse_matches_invert(maximal, 200,
                                                  without_kummer) >= 40);
    }
    {
        sflint::FmpqPoly polynomial;
        poly_cubic_disc1724(polynomial);
        silex::NumberField field;
        silex::Order equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(check_prime_inverse_matches_invert(maximal, 100,
                                                  without_kummer) >= 20);
    }
    {
        sflint::FmpqPoly polynomial;
        poly_index3_quartic(polynomial);
        silex::NumberField field;
        silex::Order equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(check_prime_inverse_matches_invert(maximal, 100,
                                                  without_kummer) >= 20);
    }
    assert(without_kummer > before_index_fields);
    return 0;
}

// The route declines inputs outside its contract and leaves its outputs
// unchanged: an order not known to be maximal, the unit ideal, and an
// ideal not containing p O.
int test_prime_ideal_inverse_declines_outside_contract() {
    sflint::FmpqPoly polynomial;
    poly_cubic_disc1724(polynomial);
    silex::NumberField field;
    silex::Order equation = order_by_polynomial(field, polynomial);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(!equation.is_maximal());

    sflint::Fmpz p;
    assert(set_fmpz_si(p, 5));
    silex::PrimeIdealList equation_primes;
    assert(silex::decompose_prime(equation_primes, equation,
                                  sflint::FmpzConstRef(p)));
    silex::Ideal equation_prime(equation);
    assert(equation_primes.at(0)->get_ideal(equation_prime));

    silex::Ideal numerator(equation);
    sflint::Fmpz den;
    assert(set_fmpz_si(den, 7));
    assert(numerator.one());
    assert(!silex::detail::prime_ideal_inverse_num_den(
            numerator, sflint::FmpzRef(den), equation_prime));
    assert(numerator.is_one() && sflint::fmpz_equal_si(den, 7));

    silex::Ideal unit(maximal);
    silex::Ideal maximal_numerator(maximal);
    assert(unit.one());
    assert(maximal_numerator.one());
    assert(!silex::detail::prime_ideal_inverse_num_den(
            maximal_numerator, sflint::FmpzRef(den), unit));
    assert(maximal_numerator.is_one() && sflint::fmpz_equal_si(den, 7));

    // P^2 for an unramified prime P of degree 1 has norm p^2 and does not
    // contain p O (P^2 does not divide p O), so the route must decline.  An
    // entry p^2 on the HNF diagonal is declined by the diagonal check; at
    // p = 2 the diagonal is (2, 2, 1), in {1, p}, so only the containment
    // check declines it.  Both kinds are asserted to occur.  The containment
    // check is backed up by later checks (kernel dimension, index of
    // p O + t O), so no ideal isolates it; the case guards the combined
    // decline.
    slong squares_checked = 0;
    slong containment_declines = 0;
    for (slong prime = 2; prime < 60; ++prime) {
        if (n_is_prime(static_cast<ulong>(prime)) == 0) {
            continue;
        }
        sflint::Fmpz q;
        assert(set_fmpz_si(q, prime));
        silex::PrimeIdealList primes;
        if (!silex::decompose_prime(primes, maximal,
                                    sflint::FmpzConstRef(q))) {
            continue;
        }
        for (slong i = 0; i < primes.size(); ++i) {
            if (primes.at(i)->residue_degree() != 1 ||
                primes.at(i)->ramification_index() != 1) {
                continue;
            }
            silex::Ideal prime_ideal(maximal);
            silex::Ideal square(maximal);
            assert(primes.at(i)->get_ideal(prime_ideal));
            assert(square.multiply(prime_ideal, prime_ideal));
            assert(!silex::detail::prime_ideal_inverse_num_den(
                    maximal_numerator, sflint::FmpzRef(den), square));
            assert(maximal_numerator.is_one() &&
                   sflint::fmpz_equal_si(den, 7));
            // Record the squares whose HNF diagonal lies in {1, q}: they
            // pass the diagonal check and are declined only because q O is
            // not contained in them.
            sflint::FmpzMat hnf(3, 3);
            assert(square.get_hnf(sflint::FmpzMatRef(hnf)));
            bool diagonal_in_one_q = true;
            for (slong d = 0; d < 3; ++d) {
                const fmpz* entry = fmpz_mat_entry(hnf.raw(), d, d);
                diagonal_in_one_q = diagonal_in_one_q &&
                        (fmpz_is_one(entry) != 0 ||
                         fmpz_equal_si(entry, prime) != 0);
            }
            if (diagonal_in_one_q) {
                ++containment_declines;
            }
            ++squares_checked;
        }
    }
    // The loop must find unramified degree-one primes, so this case cannot
    // silently skip.
    assert(squares_checked >= 3);
    assert(containment_declines >= 1);
    return 0;
}

// Minimal polynomial of the integral generator omega of Q(sqrt(d)) for
// squarefree d: x^2 - x - (d - 1)/4 when d = 1 mod 4, else x^2 - d.
void omega_minimal_polynomial(sflint::FmpzPoly& out,
                              const sflint::Fmpz& radicand) noexcept {
    sflint::Fmpz constant;
    fmpz_poly_zero(out.raw());
    fmpz_poly_set_coeff_si(out.raw(), 2, 1);
    if (fmpz_fdiv_ui(radicand.raw(), 4) == 1) {
        fmpz_poly_set_coeff_si(out.raw(), 1, -1);
        fmpz_sub_ui(constant.raw(), radicand.raw(), 1);
        fmpz_divexact_ui(constant.raw(), constant.raw(), 4);
        fmpz_neg(constant.raw(), constant.raw());
    } else {
        fmpz_neg(constant.raw(), radicand.raw());
    }
    fmpz_poly_set_coeff_fmpz(out.raw(), 0, constant.raw());
}

// Same factors with the same exponents, in any order.
bool same_factor_set(const sflint::FmpzModPolyFactor& a,
                     const sflint::FmpzModPolyFactor& b,
                     const sflint::FmpzModCtx& ctx) noexcept {
    if (a.raw()->num != b.raw()->num) {
        return false;
    }
    for (slong i = 0; i < a.raw()->num; ++i) {
        bool found = false;
        for (slong j = 0; j < b.raw()->num; ++j) {
            if (fmpz_mod_poly_equal(a.raw()->poly + i, b.raw()->poly + j,
                                    ctx.raw()) != 0 &&
                a.raw()->exp[i] == b.raw()->exp[j]) {
                found = true;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

// The root r of a linear factor x - r, in [0, p).
void linear_factor_root(sflint::Fmpz& root, const fmpz_mod_poly_struct* factor,
                        const sflint::FmpzModCtx& ctx) noexcept {
    fmpz_mod_poly_get_coeff_fmpz(root.raw(), factor, 0, ctx.raw());
    fmpz_mod_neg(root.raw(), root.raw(), ctx.raw());
}

struct ClosedFormTally {
    slong split = 0;
    slong ramified = 0;
    slong inert = 0;
};

void check_quadratic_closed_form_factor(
        const sflint::FmpzPoly& minimal_polynomial,
        const sflint::Fmpz& p,
        ClosedFormTally& tally) noexcept {
    sflint::FmpzModCtx ctx(p.raw());
    sflint::FmpzModPoly reduced(ctx);
    sflint::FmpzModPolyFactor closed(ctx);
    sflint::FmpzModPolyFactor generic(ctx);
    fmpz_mod_poly_set_fmpz_poly(reduced.raw(), minimal_polynomial.raw(),
                                ctx.raw());
    fmpz_mod_poly_factor(generic.raw(), reduced.raw(), ctx.raw());
    assert(silex::detail::monic_quadratic_factor_word(
            closed, sflint::FmpzModPolyConstRef(reduced), ctx));
    assert(same_factor_set(closed, generic, ctx));
    const slong num = closed.raw()->num;
    if (num == 2) {
        // Canonical order: increasing root.
        sflint::Fmpz r0;
        sflint::Fmpz r1;
        linear_factor_root(r0, closed.raw()->poly + 0, ctx);
        linear_factor_root(r1, closed.raw()->poly + 1, ctx);
        assert(fmpz_cmp(r0.raw(), r1.raw()) < 0);
        ++tally.split;
    } else if (closed.raw()->exp[0] == 2) {
        ++tally.ramified;
    } else {
        assert(fmpz_mod_poly_degree(closed.raw()->poly, ctx.raw()) == 2);
        ++tally.inert;
    }
}

int test_quadratic_closed_form_factor_matches_flint() {
    // The closed-form factorization of omega's minimal polynomial gives the
    // factors of fmpz_mod_poly_factor (as sets) with split factors by
    // increasing root, for p = 2, small odd p, p dividing d, and word-size
    // p near 2^63 and 2^64.
    const slong small_radicands[] = {-1, -2, -3, -5, -7, -15, -23, -47,
                                     2, 3, 5, 6, 7, 13, 17, 21, 33, 105,
                                     -100000000000319, 10000000005};
    std::vector<sflint::Fmpz> radicands;
    for (const slong d : small_radicands) {
        radicands.emplace_back();
        fmpz_set_si(radicands.back().raw(), d);
    }
    // Radicands that are themselves large primes, so p = |d| is ramified.
    const ulong mersenne61 = (UWORD(1) << 61) - 1;
    const ulong below_2_64 = UWORD(0xFFFFFFFFFFFFFFC5);  // 2^64 - 59, prime
    for (const ulong q : {mersenne61, below_2_64}) {
        radicands.emplace_back();
        fmpz_set_ui(radicands.back().raw(), q);
        radicands.emplace_back();
        fmpz_set_ui(radicands.back().raw(), q);
        fmpz_neg(radicands.back().raw(), radicands.back().raw());
    }

    std::vector<ulong> primes;
    for (ulong q = 2; q < 2000; q = n_nextprime(q, 1)) {
        primes.push_back(q);
    }
    ulong q = UWORD(1) << 62;
    for (int i = 0; i < 24; ++i) {
        q = n_nextprime(q, 1);
        primes.push_back(q);
    }
    q = UWORD(1) << 63;
    for (int i = 0; i < 24; ++i) {
        q = n_nextprime(q, 1);
        primes.push_back(q);
    }
    q = below_2_64 - 4000;
    while (q < below_2_64) {
        q = n_nextprime(q, 1);
        primes.push_back(q);
    }
    primes.push_back(mersenne61);
    assert(n_is_prime(mersenne61) != 0 && n_is_prime(below_2_64) != 0);

    // omega's minimal polynomials, and monic quadratics x^2 + b x + c with
    // other b (defining polynomials on the Dedekind-Kummer routes).
    std::vector<sflint::FmpzPoly> polynomials;
    for (const sflint::Fmpz& radicand : radicands) {
        polynomials.emplace_back();
        omega_minimal_polynomial(polynomials.back(), radicand);
    }
    const slong general[][2] = {{3, 7}, {2, 5}, {3, -27}, {1, 9}, {4, 1},
                                {-7, 10}, {5, 0}, {-1, 25000000000080},
                                {101, -2500000001}};
    for (const auto& bc : general) {
        polynomials.emplace_back();
        fmpz_poly_set_coeff_si(polynomials.back().raw(), 2, 1);
        fmpz_poly_set_coeff_si(polynomials.back().raw(), 1, bc[0]);
        fmpz_poly_set_coeff_si(polynomials.back().raw(), 0, bc[1]);
    }

    ClosedFormTally large_tally;
    ClosedFormTally tally;
    for (const sflint::FmpzPoly& polynomial : polynomials) {
        for (const ulong prime : primes) {
            sflint::Fmpz p;
            fmpz_set_ui(p.raw(), prime);
            check_quadratic_closed_form_factor(
                    polynomial, p,
                    prime > (UWORD(1) << 62) ? large_tally : tally);
        }
    }
    assert(tally.split > 0 && tally.ramified > 0 && tally.inert > 0);
    assert(large_tally.split > 0 && large_tally.ramified > 0 &&
           large_tally.inert > 0);

    // A prime beyond a word is left to the generic factorization.
    sflint::Fmpz radicand;
    sflint::Fmpz p;
    fmpz_set_si(radicand.raw(), -5);
    fmpz_set_ui(p.raw(), UWORD(1));
    fmpz_mul_2exp(p.raw(), p.raw(), 64);
    fmpz_nextprime(p.raw(), p.raw(), 1);
    sflint::FmpzPoly minimal_polynomial;
    omega_minimal_polynomial(minimal_polynomial, radicand);
    sflint::FmpzModCtx ctx(p.raw());
    sflint::FmpzModPoly reduced(ctx);
    fmpz_mod_poly_set_fmpz_poly(reduced.raw(), minimal_polynomial.raw(),
                                ctx.raw());
    sflint::FmpzModPolyFactor closed(ctx);
    assert(!silex::detail::monic_quadratic_factor_word(
            closed, sflint::FmpzModPolyConstRef(reduced), ctx));
    // So is a polynomial of another degree, or one that is not monic.
    sflint::FmpzModCtx small_ctx(UWORD(101));
    sflint::FmpzModPoly other(small_ctx);
    sflint::FmpzModPolyFactor other_factors(small_ctx);
    fmpz_mod_poly_set_coeff_ui(other.raw(), 3, 1, small_ctx.raw());
    fmpz_mod_poly_set_coeff_ui(other.raw(), 0, 2, small_ctx.raw());
    assert(!silex::detail::monic_quadratic_factor_word(
            other_factors, sflint::FmpzModPolyConstRef(other), small_ctx));
    fmpz_mod_poly_zero(other.raw(), small_ctx.raw());
    fmpz_mod_poly_set_coeff_ui(other.raw(), 2, 3, small_ctx.raw());
    fmpz_mod_poly_set_coeff_ui(other.raw(), 0, 2, small_ctx.raw());
    assert(!silex::detail::monic_quadratic_factor_word(
            other_factors, sflint::FmpzModPolyConstRef(other), small_ctx));
    return 0;
}

int test_quadratic_decompose_prime_closed_form_matches_flint_primes() {
    // On maximal quadratic-backend orders, decompose_prime and
    // set_first_degree_one_prime build the same prime ideals (HNF, residue
    // polynomial, Kummer generator) as the ideals built from FLINT's
    // factors, as sets; split primes come by increasing root of omega and
    // set_first_degree_one_prime keeps the first.  Includes p beyond a word,
    // which keeps the generic factorization.
    using Access = silex::detail::MaximalQuadraticPrimeAccess;
    using Kind = silex::detail::RetainedQuadraticPrimeKind;
    const slong radicands[] = {-1, -2, -5, -23, 5, 13, 6, -100000000000319};
    std::vector<sflint::Fmpz> primes;
    for (ulong q = 2; q < 260; q = n_nextprime(q, 1)) {
        primes.emplace_back();
        fmpz_set_ui(primes.back().raw(), q);
    }
    for (const ulong start : {UWORD(1) << 63, UWORD(0xFFFFFFFFFFFFF000)}) {
        ulong q = start;
        for (int i = 0; i < 6; ++i) {
            q = n_nextprime(q, 1);
            primes.emplace_back();
            fmpz_set_ui(primes.back().raw(), q);
        }
    }
    {
        sflint::Fmpz q;
        fmpz_one(q.raw());
        fmpz_mul_2exp(q.raw(), q.raw(), 64);
        for (int i = 0; i < 4; ++i) {
            fmpz_nextprime(q.raw(), q.raw(), 1);
            primes.emplace_back();
            fmpz_set(primes.back().raw(), q.raw());
        }
    }

    slong checked_split = 0;
    for (const slong d : radicands) {
        silex::NumberField field = silex::test::quadratic_field(d);
        assert(field.backend_kind() ==
               silex::NumberFieldBackendKind::quadratic);
        silex::Order equation = silex::test::equation_order(field);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        silex::Element omega(field);
        sflint::Fmpz radicand;
        fmpz_set_si(radicand.raw(), d);
        if (fmpz_fdiv_ui(radicand.raw(), 4) == 1) {
            element_alpha_over(omega, 1, 1, 2);
        } else {
            assert(omega.gen());
        }
        sflint::FmpzPoly minimal_polynomial;
        omega_minimal_polynomial(minimal_polynomial, radicand);

        for (const sflint::Fmpz& p : primes) {
            silex::PrimeIdealList list;
            assert(silex::decompose_prime(list, maximal,
                                          sflint::FmpzConstRef(p)));
            sflint::FmpzModCtx ctx(p.raw());
            sflint::FmpzModPoly reduced(ctx);
            sflint::FmpzModPolyFactor generic(ctx);
            fmpz_mod_poly_set_fmpz_poly(reduced.raw(),
                                        minimal_polynomial.raw(), ctx.raw());
            fmpz_mod_poly_factor(generic.raw(), reduced.raw(), ctx.raw());
            assert(list.size() == generic.raw()->num);
            for (slong j = 0; j < generic.raw()->num; ++j) {
                silex::PrimeIdeal reference(maximal);
                assert(Access::set_from_integral_generator_factor(
                        reference, maximal, sflint::FmpzConstRef(p),
                        sflint::FmpzModPolyConstRef(generic.raw()->poly + j),
                        generic.raw()->exp[j], omega, ctx));
                slong matches = 0;
                for (slong i = 0; i < list.size(); ++i) {
                    matches += same_prime_data(reference, *list.at(i),
                                               maximal) ? 1 : 0;
                }
                assert(matches == 1);
            }

            silex::PrimeIdeal retained(maximal);
            Kind kind = Kind::split;
            assert(Access::set_first_degree_one_prime(
                    retained, kind, maximal, sflint::FmpzConstRef(p),
                    nullptr));
            if (list.size() == 1 && list.at(0)->residue_degree() == 2) {
                assert(kind == Kind::inert);
                continue;
            }
            assert(kind == (list.size() == 2 ? Kind::split : Kind::ramified));
            assert(same_prime_data(retained, *list.at(0), maximal));
            if (list.size() == 2 && fmpz_abs_fits_ui(p.raw()) != 0) {
                sflint::Fmpz root0;
                sflint::Fmpz root1;
                assert(degree_one_root_from_residue_polynomial(
                        root0, *list.at(0)));
                assert(degree_one_root_from_residue_polynomial(
                        root1, *list.at(1)));
                assert(fmpz_cmp(root0.raw(), root1.raw()) < 0);
                ++checked_split;
            }
        }
    }
    assert(checked_split > 0);
    return 0;
}

// decompose_prime on an order of a quadratic field defined by a monic
// quadratic f that is not x^2 - d (so not the quadratic backend), on the
// Dedekind-Kummer routes: the primes are those built from FLINT's factors
// of f modulo p, as sets, and for word-size p they follow the closed form's
// factor order, a split prime's two primes by increasing root of f modulo p.
void check_quadratic_dedekind_kummer_route(const silex::Order& order,
                                           const sflint::FmpzPoly& f,
                                           const sflint::Fmpz& p,
                                           slong& checked_split) {
    using Route = silex::detail::DecomposePrimeRoute;
    const Route route = silex::detail::TestHookAccess::decompose_prime_route(
            order, sflint::FmpzConstRef(p), 0);
    if (route == Route::index_divisor) {
        return;
    }
    assert(route == Route::nonindex_factor || route == Route::equation_factor);
    sflint::FmpzModCtx ctx(p.raw());
    sflint::FmpzModPoly reduced(ctx);
    sflint::FmpzModPolyFactor generic(ctx);
    fmpz_mod_poly_set_fmpz_poly(reduced.raw(), f.raw(), ctx.raw());
    fmpz_mod_poly_factor(generic.raw(), reduced.raw(), ctx.raw());

    silex::PrimeIdealList list;
    if (!silex::decompose_prime(list, order, sflint::FmpzConstRef(p))) {
        // A non-maximal equation order at a prime where f is not squarefree.
        assert(!order.is_maximal() && route == Route::equation_factor);
        slong total = 0;
        for (slong j = 0; j < generic.raw()->num; ++j) {
            total += generic.raw()->exp[j];
        }
        assert(total > generic.raw()->num);
        return;
    }
    assert(list.size() == generic.raw()->num);
    for (slong j = 0; j < generic.raw()->num; ++j) {
        silex::PrimeIdeal reference(order);
        assert(PrimeIdealAccess::set_defining_polynomial_factor_prime(
                reference, order, sflint::FmpzConstRef(p),
                generic.raw()->poly + j, generic.raw()->exp[j], ctx, nullptr));
        slong matches = 0;
        for (slong i = 0; i < list.size(); ++i) {
            matches += same_prime_data(reference, *list.at(i), order) ? 1 : 0;
        }
        assert(matches == 1);
    }
    if (fmpz_abs_fits_ui(p.raw()) == 0) {
        return;
    }
    sflint::FmpzModPolyFactor closed(ctx);
    assert(silex::detail::monic_quadratic_factor_word(
            closed, sflint::FmpzModPolyConstRef(reduced), ctx));
    assert(closed.raw()->num == list.size());
    for (slong i = 0; i < closed.raw()->num; ++i) {
        silex::PrimeIdeal expected(order);
        assert(PrimeIdealAccess::set_defining_polynomial_factor_prime(
                expected, order, sflint::FmpzConstRef(p),
                closed.raw()->poly + i, closed.raw()->exp[i], ctx, nullptr));
        assert(same_prime_data(expected, *list.at(i), order));
    }
    if (closed.raw()->num == 2) {
        sflint::Fmpz r0;
        sflint::Fmpz r1;
        linear_factor_root(r0, closed.raw()->poly + 0, ctx);
        linear_factor_root(r1, closed.raw()->poly + 1, ctx);
        assert(fmpz_cmp(r0.raw(), r1.raw()) < 0);
        ++checked_split;
    }
}

int test_quadratic_dedekind_kummer_closed_form_matches_flint_primes() {
    // x^2 + x + 9 (D = -35) and x^2 - x - 2500000001 (D = 10000000005):
    // maximal equation orders.  x^2 + 2x + 5 (D = -16) and x^2 + 3x - 27
    // (D = 117 = 9 * 13): Z[theta] has index 2 and 3 in the maximal order,
    // which takes the nonindex route at other primes; their equation orders
    // take the equation route.  x^2 - x + 25000000000080 is
    // imaginary_quadratic_d100000000000319.
    const std::vector<std::vector<slong>> fields = {
            {9, 1, 1}, {-2500000001, -1, 1}, {5, 2, 1}, {-27, 3, 1},
            {25000000000080, -1, 1}};
    std::vector<sflint::Fmpz> primes;
    for (ulong q = 2; q < 200; q = n_nextprime(q, 1)) {
        primes.emplace_back();
        fmpz_set_ui(primes.back().raw(), q);
    }
    for (const ulong start : {UWORD(1) << 62, UWORD(1) << 63,
                              UWORD(0xFFFFFFFFFFFFF000)}) {
        ulong q = start;
        for (int i = 0; i < 6; ++i) {
            q = n_nextprime(q, 1);
            primes.emplace_back();
            fmpz_set_ui(primes.back().raw(), q);
        }
    }
    {
        sflint::Fmpz q;
        fmpz_one(q.raw());
        fmpz_mul_2exp(q.raw(), q.raw(), 64);
        for (int i = 0; i < 3; ++i) {
            fmpz_nextprime(q.raw(), q.raw(), 1);
            primes.emplace_back();
            fmpz_set(primes.back().raw(), q.raw());
        }
    }

    slong checked_split = 0;
    slong checked_nonindex = 0;
    for (const std::vector<slong>& low_to_high : fields) {
        sflint::FmpqPoly polynomial;
        poly_from_coefficients(polynomial, low_to_high);
        sflint::FmpzPoly f;
        for (std::size_t i = 0; i < low_to_high.size(); ++i) {
            fmpz_poly_set_coeff_si(f.raw(), static_cast<slong>(i),
                                   low_to_high[i]);
        }
        silex::NumberField field;
        silex::Order equation = order_by_polynomial(field, polynomial);
        assert(field.backend_kind() !=
               silex::NumberFieldBackendKind::quadratic);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        for (const sflint::Fmpz& p : primes) {
            check_quadratic_dedekind_kummer_route(maximal, f, p,
                                                  checked_split);
            if (!maximal.is_equation_order()) {
                check_quadratic_dedekind_kummer_route(equation, f, p,
                                                      checked_split);
                ++checked_nonindex;
            }
        }
    }
    assert(checked_split > 0 && checked_nonindex > 0);
    // decompose_prime_containing factors f through the same closed form, so
    // its primes are decompose_prime's, in the same order.
    assert(check_decompose_prime_containing_routes(
                   {9, 1, 1}, true, {2, 3, 5, 7, 11, 13, 29, 37}, {}) == 0);
    assert(check_decompose_prime_containing_routes(
                   {-27, 3, 1}, false, {2, 3, 5, 7, 13, 17, 23, 43}, {3}) ==
           0);
    return 0;
}

}  // namespace

// PrimeIdealAccess::same_prime_ideal, the same-prime test shared by the
// factor-base lookup and the class-group partial-relation pairing.  Over each
// rational prime with at least two primes above it, every prime matches
// itself and a copy, and no other prime above p (different ideals, whatever
// their e and f).  A copy given a different, equally valid residue
// polynomial (constant term shifted by p) still matches although
// PrimeIdeal::equal rejects it.  x^3 - x^2 - 2x - 8: 2 is a common index
// divisor and splits into three degree-one primes.  x^3 - 10: 3 = P^2 Q
// divides the index.
int test_same_prime_ideal_compares_ideals() {
    const std::vector<std::vector<slong>> fields = {
            {-8, -2, -1, 1},
            {-10, 0, 0, 1},
    };
    slong checked_index_divisor = 0;
    slong checked_non_index = 0;
    slong checked_residue_change = 0;
    for (const std::vector<slong>& coefficients : fields) {
        sflint::FmpqPoly polynomial;
        poly_from_coefficients(polynomial, coefficients);
        silex::NumberField field;
        silex::Order equation;
        equation = order_by_polynomial(field, polynomial);
        silex::Order maximal(field);
        assert(maximal.maximal_order(equation));
        assert(!maximal.is_equation_order());
        for (ulong q = 2; q < 60; q = n_nextprime(q, 0)) {
            sflint::Fmpz p;
            assert(set_fmpz_si(p, static_cast<slong>(q)));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(primes, maximal,
                                          sflint::FmpzConstRef(p)));
            if (primes.size() < 2) {
                continue;
            }
            const bool index_divisor =
                    (coefficients[0] == -8 && q == 2) ||
                    (coefficients[0] == -10 && q == 3);
            if (index_divisor) {
                ++checked_index_divisor;
            } else {
                ++checked_non_index;
            }
            for (slong i = 0; i < primes.size(); ++i) {
                const silex::PrimeIdeal& left = *primes.at(i);
                silex::PrimeIdeal copy(maximal);
                assert(copy.set(left));
                assert(PrimeIdealAccess::same_prime_ideal(left, left));
                assert(PrimeIdealAccess::same_prime_ideal(left, copy));
                assert(PrimeIdealAccess::same_prime_ideal(copy, left));
                for (slong j = 0; j < primes.size(); ++j) {
                    assert(PrimeIdealAccess::same_prime_ideal(
                                   left, *primes.at(j)) == (i == j));
                }

                sflint::FmpzPoly residue;
                if (!left.residue_polynomial(sflint::FmpzPolyRef(residue))) {
                    continue;
                }
                sflint::Fmpz constant;
                fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
                fmpz_add(constant.raw(), constant.raw(), p.raw());
                fmpz_poly_set_coeff_fmpz(residue.raw(), 0, constant.raw());
                assert(silex::detail::TestHookAccess::
                               set_alpha_residue_polynomial(
                                       copy,
                                       sflint::FmpzPolyConstRef(residue)));
                assert(!copy.equal(left));
                assert(PrimeIdealAccess::same_prime_ideal(copy, left));
                for (slong j = 0; j < primes.size(); ++j) {
                    assert(PrimeIdealAccess::same_prime_ideal(
                                   copy, *primes.at(j)) == (i == j));
                }
                ++checked_residue_change;
            }
        }
    }
    assert(checked_index_divisor == 2);
    assert(checked_non_index > 0);
    assert(checked_residue_change > 0);
    return 0;
}

// The word-size primitive-root search returns the smallest primitive root
// found by the fmpz search, and quotient_log_mod_prime_setup (word route for
// p below 2^64) gives the cofactor and quotient generator built from it.
// Covers every prime below 3000 and the primes next to 2^32, 2^63 and 2^64.
int test_primitive_root_word_size_matches_fmpz() {
    const auto previous_prime = [](ulong n) noexcept {
        do {
            --n;
        } while (n_is_prime(n) == 0);
        return n;
    };
    std::vector<ulong> primes;
    for (ulong p = 2; p < 3000; p = n_nextprime(p, 1)) {
        primes.push_back(p);
    }
    for (const ulong centre : {UWORD(1) << 32, UWORD(1) << 63}) {
        ulong below = centre;
        ulong above = centre;
        for (int k = 0; k < 6; ++k) {
            below = previous_prime(below);
            above = n_nextprime(above, 1);
            primes.push_back(below);
            primes.push_back(above);
        }
    }
    ulong top = UWORD_MAX;
    for (int k = 0; k < 6; ++k) {
        top = previous_prime(top);
        primes.push_back(top);
    }

    for (const ulong p : primes) {
        sflint::Fmpz p_fmpz;
        sflint::Fmpz expected;
        sflint::fmpz_set_ui(sflint::FmpzRef(p_fmpz), p);
        assert(silex::detail::primitive_root_mod_prime(
                expected, sflint::FmpzConstRef(p_fmpz)));
        assert(sflint::fmpz_abs_fits_ui(sflint::FmpzConstRef(expected)));
        const ulong root = silex::detail::primitive_root_mod_prime_ui(p);
        assert(root == sflint::fmpz_get_ui(sflint::FmpzConstRef(expected)));

        n_factor_t factors;
        n_factor_init(&factors);
        if (p > 2) {
            n_factor(&factors, p - 1, 1);
        }
        for (int i = 0; i < factors.num; ++i) {
            sflint::Fmpz ell;
            sflint::Fmpz cofactor;
            sflint::Fmpz quotient_generator;
            sflint::Fmpz expected_cofactor;
            sflint::Fmpz expected_generator;
            sflint::fmpz_set_ui(sflint::FmpzRef(ell), factors.p[i]);
            sflint::fmpz_set_ui(sflint::FmpzRef(expected_cofactor),
                                (p - 1) / factors.p[i]);
            fmpz_powm(expected_generator.raw(), expected.raw(),
                      expected_cofactor.raw(), p_fmpz.raw());
            assert(silex::detail::quotient_log_mod_prime_setup(
                    cofactor, quotient_generator,
                    sflint::FmpzConstRef(p_fmpz), sflint::FmpzConstRef(ell)));
            assert(fmpz_equal(cofactor.raw(), expected_cofactor.raw()) != 0);
            assert(fmpz_equal(quotient_generator.raw(),
                              expected_generator.raw()) != 0);
        }
    }

    // A prime ell that does not divide p - 1, or lies above the word size,
    // is refused as before.
    sflint::Fmpz cofactor;
    sflint::Fmpz quotient_generator;
    sflint::Fmpz p_fmpz;
    sflint::Fmpz ell;
    sflint::fmpz_set_ui(sflint::FmpzRef(p_fmpz), 13);
    sflint::fmpz_set_ui(sflint::FmpzRef(ell), 5);
    assert(!silex::detail::quotient_log_mod_prime_setup(
            cofactor, quotient_generator, sflint::FmpzConstRef(p_fmpz),
            sflint::FmpzConstRef(ell)));
    fmpz_set_ui(ell.raw(), UWORD_MAX);
    fmpz_nextprime(ell.raw(), ell.raw(), 1);
    assert(!silex::detail::quotient_log_mod_prime_setup(
            cofactor, quotient_generator, sflint::FmpzConstRef(p_fmpz),
            sflint::FmpzConstRef(ell)));

    // Above the word size the fmpz route still answers: p is the next prime
    // after 2^64 - 1 and ell = 2 divides p - 1.
    fmpz_set_ui(p_fmpz.raw(), UWORD_MAX);
    fmpz_nextprime(p_fmpz.raw(), p_fmpz.raw(), 1);
    sflint::fmpz_set_ui(sflint::FmpzRef(ell), 2);
    sflint::Fmpz generator;
    assert(silex::detail::primitive_root_mod_prime(
            generator, sflint::FmpzConstRef(p_fmpz)));
    assert(silex::detail::quotient_log_mod_prime_setup(
            cofactor, quotient_generator, sflint::FmpzConstRef(p_fmpz),
            sflint::FmpzConstRef(ell)));
    sflint::Fmpz expected_cofactor;
    sflint::Fmpz expected_generator;
    fmpz_sub_ui(expected_cofactor.raw(), p_fmpz.raw(), 1);
    fmpz_fdiv_q_2exp(expected_cofactor.raw(), expected_cofactor.raw(), 1);
    fmpz_powm(expected_generator.raw(), generator.raw(),
              expected_cofactor.raw(), p_fmpz.raw());
    assert(fmpz_equal(cofactor.raw(), expected_cofactor.raw()) != 0);
    assert(fmpz_equal(quotient_generator.raw(), expected_generator.raw()) !=
           0);
    return 0;
}

int main() {
    assert(test_degree_one() == 0);
    assert(test_same_prime_ideal_compares_ideals() == 0);
    assert(test_primitive_root_word_size_matches_fmpz() == 0);
    assert(test_degree_one_prime_from_root_matches_decomposition() == 0);
    assert(test_degree_one_prime_from_root_rejects_non_roots() == 0);
    assert(test_degree_one_prime_from_root_rejects_index_prime() == 0);
    assert(test_alpha_reduction_refuses_index_prime() == 0);
    assert(test_invalid_prime_failure_preserves_output() == 0);
    assert(test_unsupported_order_failure_preserves_output() == 0);
    assert(test_quadratic_splitting_types() == 0);
    assert(test_explicit_maximal_order_decomposition() == 0);
    assert(test_split_index_prime_uses_integral_quadratic_generator() == 0);
    assert(test_index3_quartic_maximal_nonindex_prime_decomposition() == 0);
    assert(test_index2_cubic_unramified_prime_decomposition() == 0);
    assert(test_index8_quartic_ramified_prime_decomposition() == 0);
    assert(test_index2_quartic_mixed_ramified_prime_decomposition() == 0);
    assert(test_index2_quintic_mixed_ramified_prime_decomposition() == 0);
    assert(test_index_divisor_recursive_split_unramified_degree_one() == 0);
    assert(test_index_divisor_recursive_split_unramified_degree_two() == 0);
    assert(test_index_divisor_recursive_split_ramified() == 0);
    assert(test_index_divisor_monogenic_order_unchanged_unramified() == 0);
    assert(test_index_divisor_monogenic_order_unchanged_ramified() == 0);
    assert(test_index_divisor_cache_more_primes_than_p() == 0);
    assert(test_index_divisor_cache_ramified_hnf() == 0);
    assert(test_index_divisor_cache_recursive_split() == 0);
    assert(test_index_divisor_cache_monogenic() == 0);
    assert(test_index_divisor_cache_inert_principal() == 0);
    assert(test_residue_degrees_cyclotomic_real_q10_c11() == 0);
    assert(test_residue_degrees_index_divisor_routes() == 0);
    assert(test_decompose_prime_route() == 0);
    assert(test_decompose_prime_containing() == 0);
    assert(test_index_divisor_separating_element_order() == 0);
    assert(test_recursive_split_matches_enumeration_unramified() == 0);
    assert(test_recursive_split_matches_enumeration_ramified() == 0);
    assert(test_index_divisor_recursive_split_multiword_prime() == 0);
    assert(test_index_divisor_recursive_split_full_word_prime() == 0);
    assert(test_index_divisor_ramified_multiword_prime() == 0);
    assert(test_valuation_degree_one() == 0);
    assert(test_valuation_order_element_split_repeated() == 0);
    assert(test_valuation_order_element_ramified_unique_prime_norm() == 0);
    assert(test_valuation_order_element_ramified_split_coordinate_matrix() ==
           0);
    assert(test_valuation_order_element_ramified_nonunique_containment() == 0);
    assert(test_valuation_power_containment_matches_definition() == 0);
    assert(test_valuation_power_containment_matches_gp() == 0);
    assert(test_valuation_anti_uniformizer_matches_containment() == 0);
    assert(test_principal_ideal_valuation_matches_element_valuation() ==
           0);
    assert(test_valuation_fractional_element_and_factored() == 0);
    assert(test_extended_valuation_failure_preserves_output() == 0);
    assert(test_valuation_maximal_quadratic() == 0);
    assert(test_valuation_unmarked_failure_preserves_output() == 0);
    assert(test_cubic_mixed_decomposition() == 0);
    assert(test_residue_degree_limited_decomposition() == 0);
    assert(test_kummer_generator_coordinates_are_centered_and_copied() == 0);
    assert(test_uncertified_repeated_factor_failure_preserves_output() == 0);
    assert(test_residue_ring_prime_relation() == 0);
    assert(test_move_clear_and_redefine() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    assert(test_generic_nonpower_basis_reduce_quadratic() == 0);
    assert(test_reduce_homomorphism_on_index_maximal_orders() == 0);
    assert(test_quadratic_backend_omega_reduce_unchanged() == 0);
    assert(test_residue_convention_fixed_at_decomposition() == 0);
    assert(test_kummer_generator_reduces_to_zero_and_residues_span() == 0);
    assert(test_degree_one_root_follows_recorded_convention() == 0);
    assert(test_factored_value_nmod_matches_fmpz_route() == 0);
    assert(test_factored_value_nmod_matches_fmpz_route_large() == 0);
    assert(test_alpha_convention_kept_after_marking_maximal() == 0);
    assert(test_prime_ideal_inverse_matches_invert() == 0);
    assert(test_prime_ideal_inverse_declines_outside_contract() == 0);
    assert(test_quadratic_closed_form_factor_matches_flint() == 0);
    assert(test_quadratic_decompose_prime_closed_form_matches_flint_primes() ==
           0);
    assert(test_quadratic_dedekind_kummer_closed_form_matches_flint_primes() ==
           0);
    return 0;
}
