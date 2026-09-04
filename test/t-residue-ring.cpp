#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_ring.hpp>

#include "test_support.hpp"

#include <cassert>

namespace {
namespace sflint = silex::flint;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, col), value);
}

bool set_rational_principal(silex::Ideal& ideal, slong value) noexcept {
    const silex::Order* order = ideal.parent();
    if (order == nullptr) {
        return false;
    }

    silex::OrderElement generator(*order);
    return generator.set_si(value) && ideal.set_principal(generator);
}

void element_from_coeff(silex::Element& element,
                        slong index,
                        slong num,
                        slong den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coeff;
    sflint::fmpq_set_si(coeff, num, den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, index, coeff);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool elem_coord_is_si(const silex::ResidueElement& element,
                      slong col,
                      slong value) noexcept {
    sflint::FmpzMat row(1, element.degree());
    return element.get_coordinates(sflint::FmpzMatRef(row)) &&
           mat_entry_is_si(row, 0, col, value);
}

bool set_elem_si(silex::ResidueElement& element, slong value) noexcept {
    const silex::ResidueRing* ring = element.parent();
    if (ring == nullptr || ring->parent_order() == nullptr ||
        ring->parent_order()->parent() == nullptr) {
        return false;
    }

    silex::Element ambient(*ring->parent_order()->parent());
    return ambient.set_si(value) && element.set_element(ambient);
}

silex::ResidueElement local_residue_element() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 7));
    silex::ResidueRing ring(ideal);
    silex::ResidueElement element(ring);
    assert(set_elem_si(element, 16));
    return element;
}

int test_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    silex::Ideal copy(order);
    assert(set_rational_principal(ideal, 6));

    silex::ResidueRing ring(ideal);
    assert(ring.is_defined());
    sflint::Fmpz cardinality;
    assert(ring.cardinality(sflint::FmpzRef(cardinality)));
    assert(sflint::fmpz_equal_si(cardinality, 6));
    auto owned_cardinality = ring.cardinality();
    assert(owned_cardinality.has_value());
    assert(sflint::fmpz_equal_si(*owned_cardinality, 6));
    assert(ring.get_modulus(copy));
    assert(copy.equal(ideal));

    silex::ResidueElement a(ring);
    silex::ResidueElement b(ring);
    silex::ResidueElement c(ring);
    silex::ResidueElement one(ring);

    assert(set_elem_si(a, 16));
    assert(elem_coord_is_si(a, 0, 4));

    assert(set_elem_si(b, 5));
    assert(elem_coord_is_si(b, 0, 5));

    assert(one.one());
    assert(elem_coord_is_si(one, 0, 1));

    assert(c.add(a, b));
    assert(elem_coord_is_si(c, 0, 3));

    assert(c.multiply(a, b));
    assert(elem_coord_is_si(c, 0, 2));

    assert(c.negate(a));
    assert(elem_coord_is_si(c, 0, 2));

    assert(c.subtract(a, b));
    assert(elem_coord_is_si(c, 0, 5));

    assert(a.add(a, b));
    assert(elem_coord_is_si(a, 0, 3));

    assert(c.zero());
    sflint::FmpzMat row(1, 1);
    assert(c.get_coordinates(sflint::FmpzMatRef(row)));
    assert(mat_entry_is_si(row, 0, 0, 0));
    auto owned_coords = c.coordinates();
    assert(owned_coords.has_value());
    assert(sflint::fmpz_mat_equal(*owned_coords, row));

    return 0;
}

int test_element_keeps_parent_ring_alive() {
    silex::ResidueElement element = local_residue_element();
    assert(element.is_defined());
    assert(element.parent() != nullptr);
    assert(element.parent()->parent_order() != nullptr);
    assert(element.parent()->parent_order()->parent() != nullptr);
    assert(element.degree() == 1);
    assert(elem_coord_is_si(element, 0, 2));

    silex::Ideal modulus(*element.parent()->parent_order());
    sflint::Fmpz cardinality;
    assert(element.parent()->get_modulus(modulus));
    assert(element.parent()->cardinality(sflint::FmpzRef(cardinality)));
    assert(sflint::fmpz_equal_si(cardinality, 7));

    silex::OrderElement lift(*element.parent()->parent_order());
    assert(element.lift(lift));
    return 0;
}

int test_copy_swap_aliasing() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 7));
    silex::ResidueRing ring(ideal);
    silex::ResidueElement a(ring);
    silex::ResidueElement b(ring);
    silex::ResidueElement c(ring);

    assert(set_elem_si(a, 3));
    assert(b.set(a));
    assert(a.equal(b));

    assert(set_elem_si(c, 5));
    b.swap(c);
    assert(elem_coord_is_si(b, 0, 5));
    assert(elem_coord_is_si(c, 0, 3));

    assert(b.multiply(b, b));
    assert(elem_coord_is_si(b, 0, 4));
    assert(b.subtract(b, b));
    assert(elem_coord_is_si(b, 0, 0));

    return 0;
}

int test_quadratic_mod_two() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 2));
    silex::ResidueRing ring(ideal);
    silex::ResidueElement theta(ring);
    silex::ResidueElement theta2(ring);
    silex::ResidueElement five(ring);

    silex::Element x(field);
    assert(x.gen());
    assert(theta.set_element(x));
    assert(theta2.multiply(theta, theta));

    assert(x.set_si(5));
    assert(five.set_element(x));
    assert(theta2.equal(five));
    assert(elem_coord_is_si(five, 0, 1));

    return 0;
}

int test_failure_preserves_output() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal ideal(order);
    assert(set_rational_principal(ideal, 2));
    silex::ResidueRing ring(ideal);
    silex::ResidueElement a(ring);

    assert(set_elem_si(a, 3));
    assert(elem_coord_is_si(a, 0, 1));

    silex::Element theta_over_two(field);
    element_from_coeff(theta_over_two, 1, 1, 2);
    assert(!a.set_element(theta_over_two));
    assert(elem_coord_is_si(a, 0, 1));

    return 0;
}

int test_crt_degree_one() {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal i(order);
    silex::Ideal j(order);
    silex::Ideal product(order);
    silex::Ideal intersection(order);
    assert(set_rational_principal(i, 3));
    assert(set_rational_principal(j, 5));
    assert(product.multiply(i, j));
    assert(intersection.intersect(i, j));
    assert(product.equal(intersection));

    silex::ResidueRing ri(i);
    silex::ResidueRing rj(j);
    silex::ResidueRing rp(product);
    silex::ResidueRing wrong(i);
    silex::ResidueElement a(ri);
    silex::ResidueElement b(rj);
    silex::ResidueElement c(rp);
    silex::ResidueElement d(wrong);

    assert(set_elem_si(a, 2));
    assert(set_elem_si(b, 4));
    assert(set_elem_si(c, 7));
    assert(c.crt(a, b));
    assert(elem_coord_is_si(c, 0, 14));

    assert(set_elem_si(d, 1));
    assert(!d.crt(a, b));
    assert(elem_coord_is_si(d, 0, 1));

    assert(set_rational_principal(i, 6));
    assert(set_rational_principal(j, 10));
    assert(product.multiply(i, j));
    silex::ResidueRing r6(i);
    silex::ResidueRing r10(j);
    silex::ResidueRing r60(product);
    silex::ResidueElement a6(r6);
    silex::ResidueElement b10(r10);
    silex::ResidueElement c60(r60);
    assert(set_elem_si(a6, 1));
    assert(set_elem_si(b10, 3));
    assert(set_elem_si(c60, 7));
    assert(!c60.crt(a6, b10));
    assert(elem_coord_is_si(c60, 0, 7));

    return 0;
}

int test_crt_quadratic_equation_order() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 2);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    silex::Ideal i(order);
    silex::Ideal j(order);
    silex::Ideal product(order);
    assert(set_rational_principal(i, 3));
    assert(set_rational_principal(j, 5));
    assert(product.multiply(i, j));

    silex::ResidueRing ri(i);
    silex::ResidueRing rj(j);
    silex::ResidueRing rp(product);
    silex::ResidueElement a(ri);
    silex::ResidueElement b(rj);
    silex::ResidueElement c(rp);

    silex::Element x(field);
    assert(x.gen());
    assert(a.set_element(x));
    assert(x.set_si(2));
    assert(b.set_element(x));
    assert(c.crt(a, b));
    assert(elem_coord_is_si(c, 0, 12));
    assert(elem_coord_is_si(c, 1, 10));

    return 0;
}

int test_crt_quadratic_split_primes() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);

    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order order = silex::test::equation_order(field);

    sflint::Fmpz p;
    sflint::fmpz_set_si(p, 11);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    const silex::PrimeIdeal* first_prime = primes.at(0);
    const silex::PrimeIdeal* second_prime = primes.at(1);
    assert(first_prime != nullptr);
    assert(second_prime != nullptr);

    silex::Ideal i(order);
    silex::Ideal j(order);
    silex::Ideal product(order);
    assert(first_prime->get_ideal(i));
    assert(second_prime->get_ideal(j));
    assert(product.multiply(i, j));

    silex::ResidueRing ri(i);
    silex::ResidueRing rj(j);
    silex::ResidueRing rp(product);
    silex::ResidueElement a(ri);
    silex::ResidueElement b(rj);
    silex::ResidueElement c(rp);

    assert(set_elem_si(a, 1));
    assert(set_elem_si(b, 2));
    assert(c.crt(a, b));

    silex::Element lift(field);
    assert(c.lift(lift));

    silex::ResidueElement check_i(ri);
    silex::ResidueElement check_j(rj);
    assert(check_i.set_element(lift));
    assert(check_j.set_element(lift));
    assert(check_i.equal(a));
    assert(check_j.equal(b));

    return 0;
}

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_element_keeps_parent_ring_alive() == 0);
    assert(test_copy_swap_aliasing() == 0);
    assert(test_quadratic_mod_two() == 0);
    assert(test_failure_preserves_output() == 0);
    assert(test_crt_degree_one() == 0);
    assert(test_crt_quadratic_equation_order() == 0);
    assert(test_crt_quadratic_split_primes() == 0);
    return 0;
}
