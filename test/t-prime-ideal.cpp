#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/nmod_mat.hpp>
#include <silex/factored_element.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_ring.hpp>

#include "order/order_internal.hpp"
#include "prime_ideal/prime_ideal_internal.hpp"
#include "residue_field/residue_field_internal.hpp"
#include "test_support.hpp"

#include <cassert>
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
// internal constructor refuses it and leaves its output unchanged.
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

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_degree_one_prime_from_root_matches_decomposition() == 0);
    assert(test_degree_one_prime_from_root_rejects_non_roots() == 0);
    assert(test_degree_one_prime_from_root_rejects_index_prime() == 0);
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
    assert(test_valuation_degree_one() == 0);
    assert(test_valuation_order_element_split_repeated() == 0);
    assert(test_valuation_order_element_ramified_unique_prime_norm() == 0);
    assert(test_valuation_order_element_ramified_split_coordinate_matrix() ==
           0);
    assert(test_valuation_order_element_ramified_nonunique_containment() == 0);
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
    assert(test_alpha_convention_kept_after_marking_maximal() == 0);
    return 0;
}
