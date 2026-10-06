#include <silex/flint/fmpq_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_ring.hpp>

#include "residue_ring/unit_group_internal.hpp"
#include "test_support.hpp"

#include <flint/fmpz.h>

#include <cassert>
#include <initializer_list>
#include <vector>

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


// Unit group (O_K/m)^x.  Expected cyclic structures are from GP 2.17.4,
// which lists them in decreasing order; ResidueUnitGroup lists them as
// d_1 | d_2 | ... | d_r.  Rational moduli:
//   K = bnfinit(pol); idealstar(K, m).cyc
// Prime-power moduli:
//   K = bnfinit(pol); P = idealprimedec(K, p)[i];
//   idealstar(K, idealpow(K, P, k)).cyc
// with i chosen for the stated ramification index and residue degree.

silex::Order maximal_order_of(std::initializer_list<slong> coefficients) noexcept {
    sflint::FmpqPoly polynomial;
    slong index = 0;
    for (slong coefficient : coefficients) {
        sflint::fmpq_poly_set_coeff_si(polynomial, index, coefficient);
        ++index;
    }
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    return maximal;
}

bool is_unit_mod(const silex::OrderElement& element,
                 const silex::Ideal& modulus) noexcept {
    if (element.equal_si(0)) {
        return modulus.is_one();
    }
    silex::Ideal principal(*modulus.parent());
    silex::Ideal sum(*modulus.parent());
    assert(principal.set_principal(element));
    assert(sum.add(principal, modulus));
    return sum.is_one();
}

bool congruent_mod(const silex::OrderElement& left,
                   const silex::OrderElement& right,
                   const silex::ResidueRing& ring) noexcept {
    silex::ResidueElement a(ring);
    silex::ResidueElement b(ring);
    assert(a.set_order_element(left));
    assert(b.set_order_element(right));
    return a.equal(b);
}

void check_dlog_round_trip(const silex::detail::ResidueUnitGroup& group,
                           const silex::ResidueRing& ring,
                           const silex::OrderElement& element,
                           const std::vector<slong>& invariants) noexcept {
    const slong r = group.generator_count();
    sflint::FmpzMat log(1, r);
    assert(group.discrete_log(sflint::FmpzMatRef(log), element));
    for (slong i = 0; i < r; ++i) {
        const fmpz* entry = sflint::fmpz_mat_entry(log, 0, i).raw();
        assert(::fmpz_sgn(entry) >= 0);
        assert(::fmpz_cmp_si(entry, invariants[static_cast<std::size_t>(i)]) <
               0);
    }
    silex::OrderElement back(*ring.parent_order());
    assert(group.exp(back, sflint::FmpzMatConstRef(log)));
    assert(congruent_mod(back, element, ring));
}

// Checks the structure of `group` against `expected` (d_1 | ... | d_r) and
// runs the discrete log on every residue class of O/m: units round-trip
// through exp, non-units fail and leave the output unchanged.  Together with
// g_i^(d_i) = 1 and #units = prod d_i this shows that
// exp: Z/d_1 x ... x Z/d_r -> (O/m)^x is an isomorphism.
void check_unit_group(const silex::detail::ResidueUnitGroup& group,
                      const silex::Ideal& modulus,
                      const std::vector<slong>& expected) noexcept {
    const silex::Order& order = *modulus.parent();
    const slong n = order.degree();
    const slong r = static_cast<slong>(expected.size());
    assert(group.is_defined());
    assert(group.parent_order() != nullptr);
    assert(group.generator_count() == r);

    silex::Ideal stored(order);
    assert(group.get_modulus(stored));
    assert(stored.equal(modulus));

    sflint::FmpzVec invariants(r);
    assert(group.invariants(sflint::FmpzVecRef(invariants)));
    slong product = 1;
    for (slong i = 0; i < r; ++i) {
        assert(::fmpz_equal_si(invariants.data() + i,
                               expected[static_cast<std::size_t>(i)]));
        product *= expected[static_cast<std::size_t>(i)];
    }
    sflint::Fmpz group_order;
    assert(group.order(sflint::FmpzRef(group_order)));
    assert(sflint::fmpz_equal_si(group_order, product));

    silex::ResidueRing ring(modulus);
    silex::OrderElement one(order);
    assert(one.one());
    for (slong i = 0; i < r; ++i) {
        silex::OrderElement generator(order);
        assert(group.generator(generator, i));
        // discrete_log(generator_i) is the i-th unit vector.
        sflint::FmpzMat log(1, r);
        assert(group.discrete_log(sflint::FmpzMatRef(log), generator));
        for (slong j = 0; j < r; ++j) {
            assert(sflint::fmpz_equal_si(sflint::fmpz_mat_entry(log, 0, j),
                                         i == j ? 1 : 0));
        }
        // generator_i has order exactly d_i.
        silex::ResidueElement g(ring);
        silex::ResidueElement power(ring);
        silex::ResidueElement unit(ring);
        assert(g.set_order_element(generator));
        assert(power.one());
        assert(unit.one());
        const slong d = expected[static_cast<std::size_t>(i)];
        for (slong k = 1; k <= d; ++k) {
            assert(power.multiply(power, g));
            assert(power.equal(unit) == (k == d));
        }
    }
    silex::OrderElement wrong_index(order);
    assert(!group.generator(wrong_index, r));

    // exp accepts any integers and reduces them.
    if (r > 0) {
        sflint::FmpzMat exponents(1, r);
        sflint::FmpzMat reduced(1, r);
        for (slong i = 0; i < r; ++i) {
            const slong d = expected[static_cast<std::size_t>(i)];
            sflint::fmpz_set_si(sflint::fmpz_mat_entry(exponents, 0, i),
                                -3 * d - i - 1);
            sflint::fmpz_set_si(sflint::fmpz_mat_entry(reduced, 0, i),
                                ((-i - 1) % d + d) % d);
        }
        silex::OrderElement a(order);
        silex::OrderElement b(order);
        assert(group.exp(a, sflint::FmpzMatConstRef(exponents)));
        assert(group.exp(b, sflint::FmpzMatConstRef(reduced)));
        assert(congruent_mod(a, b, ring));
    }
    sflint::FmpzMat wrong_shape(1, r + 1);
    silex::OrderElement unchanged(order);
    assert(!group.exp(unchanged, sflint::FmpzMatConstRef(wrong_shape)));
    assert(!group.discrete_log(sflint::FmpzMatRef(wrong_shape), one));

    // Every residue class: O/m is represented by sum c_j b_j with
    // 0 <= c_j < H_jj for the triangular HNF H of m.
    sflint::FmpzMat hnf(n, n);
    assert(modulus.get_hnf(sflint::FmpzMatRef(hnf)));
    std::vector<slong> diagonal(static_cast<std::size_t>(n));
    slong classes = 1;
    for (slong j = 0; j < n; ++j) {
        diagonal[static_cast<std::size_t>(j)] =
                ::fmpz_get_si(sflint::fmpz_mat_entry(hnf, j, j).raw());
        classes *= diagonal[static_cast<std::size_t>(j)];
    }
    sflint::Fmpz norm;
    assert(modulus.norm(sflint::FmpzRef(norm)));
    assert(sflint::fmpz_equal_si(norm, classes));

    std::vector<slong> digits(static_cast<std::size_t>(n), 0);
    slong units = 0;
    for (slong index = 0; index < classes; ++index) {
        sflint::FmpzMat coordinates(1, n);
        for (slong j = 0; j < n; ++j) {
            sflint::fmpz_set_si(sflint::fmpz_mat_entry(coordinates, 0, j),
                                digits[static_cast<std::size_t>(j)]);
        }
        silex::OrderElement x(order);
        assert(x.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
        if (is_unit_mod(x, modulus)) {
            ++units;
            check_dlog_round_trip(group, ring, x, expected);
        } else {
            sflint::FmpzMat log(1, r);
            for (slong i = 0; i < r; ++i) {
                sflint::fmpz_set_si(sflint::fmpz_mat_entry(log, 0, i), 77);
            }
            assert(!group.discrete_log(sflint::FmpzMatRef(log), x));
            for (slong i = 0; i < r; ++i) {
                assert(sflint::fmpz_equal_si(sflint::fmpz_mat_entry(log, 0, i),
                                             77));
            }
        }
        for (slong j = 0; j < n; ++j) {
            slong& digit = digits[static_cast<std::size_t>(j)];
            if (++digit < diagonal[static_cast<std::size_t>(j)]) {
                break;
            }
            digit = 0;
        }
    }
    assert(units == product);
}

void check_rational_modulus(const silex::Order& order,
                            slong m,
                            const std::vector<slong>& expected) noexcept {
    silex::Ideal modulus(order);
    assert(set_rational_principal(modulus, m));
    silex::detail::ResidueUnitGroup group;
    assert(group.set_modulus(modulus));
    check_unit_group(group, modulus, expected);
}

int test_unit_group_rational_moduli() {
    // pol = x^2 - 5: 2, 3 and 7 are inert, 5 ramifies.
    // m = 2: [3]; 3: [8]; 4: [6, 2]; 8: [12, 2, 2]; 9: [24, 3];
    // 27: [72, 9]; 24: [24, 4, 2, 2]; 7: [48]; 49: [336, 7].
    {
        silex::Order order = maximal_order_of({-5, 0, 1});
        check_rational_modulus(order, 2, {3});
        check_rational_modulus(order, 3, {8});
        check_rational_modulus(order, 4, {2, 6});
        check_rational_modulus(order, 8, {2, 2, 12});
        check_rational_modulus(order, 9, {3, 24});
        check_rational_modulus(order, 27, {9, 72});
        check_rational_modulus(order, 24, {2, 2, 4, 24});
        check_rational_modulus(order, 7, {48});
        check_rational_modulus(order, 49, {7, 336});
    }
    // pol = x^2 + 1: 2 ramifies.
    // m = 2: [2]; 4: [4, 2]; 8: [4, 4, 2]; 16: [8, 4, 4]; 3: [8].
    {
        silex::Order order = maximal_order_of({1, 0, 1});
        check_rational_modulus(order, 2, {2});
        check_rational_modulus(order, 4, {2, 4});
        check_rational_modulus(order, 8, {2, 4, 4});
        check_rational_modulus(order, 16, {4, 4, 8});
        check_rational_modulus(order, 3, {8});
    }
    // pol = x^2 + 3: 3 ramifies.  m = 3: [6]; 9: [6, 3, 3].
    {
        silex::Order order = maximal_order_of({3, 0, 1});
        check_rational_modulus(order, 3, {6});
        check_rational_modulus(order, 9, {3, 3, 6});
    }
    // pol = x^2 + 7: 2 splits.  m = 2: []; 8: [2, 2, 2, 2].
    {
        silex::Order order = maximal_order_of({7, 0, 1});
        check_rational_modulus(order, 2, {});
        check_rational_modulus(order, 8, {2, 2, 2, 2});
    }
    // pol = x^2 - 2: 2 ramifies, 3 inert.  m = 3: [8]; 6: [8, 2].
    {
        silex::Order order = maximal_order_of({-2, 0, 1});
        check_rational_modulus(order, 3, {8});
        check_rational_modulus(order, 6, {2, 8});
    }
    // pol = x^3 - 2: 2 and 3 are totally ramified.
    // m = 2: [4]; 3: [6, 3]; 6: [12, 6]; 9: [18, 9, 3].
    {
        silex::Order order = maximal_order_of({-2, 0, 0, 1});
        check_rational_modulus(order, 2, {4});
        check_rational_modulus(order, 3, {3, 6});
        check_rational_modulus(order, 6, {6, 12});
        check_rational_modulus(order, 9, {3, 9, 18});
    }
    return 0;
}

bool ideal_power(silex::Ideal& out, const silex::Ideal& base, slong k) noexcept {
    silex::Ideal result(*base.parent());
    if (!result.one()) {
        return false;
    }
    for (slong i = 0; i < k; ++i) {
        silex::Ideal next(*base.parent());
        if (!next.multiply(result, base)) {
            return false;
        }
        result.swap(next);
    }
    return out.set(result);
}

void check_prime_power(std::initializer_list<slong> coefficients,
                       slong p,
                       slong ramification,
                       slong degree,
                       slong k,
                       const std::vector<slong>& expected) noexcept {
    silex::Order order = maximal_order_of(coefficients);
    sflint::Fmpz prime_number;
    sflint::fmpz_set_si(prime_number, p);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order,
                                  sflint::FmpzConstRef(prime_number)));
    const silex::PrimeIdeal* prime = nullptr;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* candidate = primes.at(i);
        if (candidate->ramification_index() == ramification &&
            candidate->residue_degree() == degree) {
            prime = candidate;
            break;
        }
    }
    assert(prime != nullptr);

    silex::Ideal prime_ideal(order);
    silex::Ideal modulus(order);
    assert(prime->get_ideal(prime_ideal));
    assert(ideal_power(modulus, prime_ideal, k));

    silex::detail::ResidueUnitGroup group;
    assert(group.set_prime_power(*prime, k));
    check_unit_group(group, modulus, expected);

    // The same group through set_modulus.
    silex::detail::ResidueUnitGroup joined;
    assert(joined.set_modulus(modulus));
    check_unit_group(joined, modulus, expected);

    // Elements 1 + b of each layer (1 + p^j)/(1 + p^(j + 1)), 1 <= j < k,
    // for b running over the HNF rows of p^j.
    silex::ResidueRing ring(modulus);
    const slong n = order.degree();
    for (slong j = 1; j < k; ++j) {
        silex::Ideal layer(order);
        assert(ideal_power(layer, prime_ideal, j));
        sflint::FmpzMat hnf(n, n);
        assert(layer.get_hnf(sflint::FmpzMatRef(hnf)));
        for (slong row = 0; row < n; ++row) {
            sflint::FmpzMat coordinates(1, n);
            for (slong col = 0; col < n; ++col) {
                sflint::fmpz_set(sflint::fmpz_mat_entry(coordinates, 0, col),
                                 sflint::FmpzConstRef(sflint::fmpz_mat_entry(hnf, row, col).raw()));
            }
            silex::OrderElement b(order);
            silex::OrderElement one(order);
            silex::OrderElement x(order);
            assert(b.set_coordinates(sflint::FmpzMatConstRef(coordinates)));
            assert(one.one());
            assert(x.add(one, b));
            assert(is_unit_mod(x, modulus));
            check_dlog_round_trip(group, ring, x, expected);
            check_dlog_round_trip(joined, ring, x, expected);
        }
    }
}

int test_unit_group_prime_powers() {
    // x^2 + 1, p = 2 (e = 2, f = 1): k = 3: [4]; k = 5: [4, 2, 2].
    check_prime_power({1, 0, 1}, 2, 2, 1, 3, {4});
    check_prime_power({1, 0, 1}, 2, 2, 1, 5, {2, 2, 4});
    // x^3 - 2, p = 3 (e = 3, f = 1), k = 4: [18, 3].
    check_prime_power({-2, 0, 0, 1}, 3, 3, 1, 4, {3, 18});
    // x^2 - 5, p = 2 (e = 1, f = 2), k = 3: [12, 2, 2].
    check_prime_power({-5, 0, 1}, 2, 1, 2, 3, {2, 2, 12});
    // x^2 + 7, p = 2 (e = 1, f = 1), k = 5: [8, 2].
    check_prime_power({7, 0, 1}, 2, 1, 1, 5, {2, 8});
    // x^3 - 2, p = 5 (e = 1, f = 1), k = 2: [20].
    check_prime_power({-2, 0, 0, 1}, 5, 1, 1, 2, {20});

    // Two primes over the same rational prime: x^2 + 7, 2 = P Q,
    // m = P^3 Q^2: idealstar(K, idealmul(K, idealpow(K, P, 3),
    // idealpow(K, Q, 2))).cyc = [2, 2, 2] for either labeling.
    {
        silex::Order order = maximal_order_of({7, 0, 1});
        sflint::Fmpz two;
        sflint::fmpz_set_si(two, 2);
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, order,
                                      sflint::FmpzConstRef(two)));
        assert(primes.size() == 2);
        silex::Ideal p(order);
        silex::Ideal q(order);
        silex::Ideal p3(order);
        silex::Ideal q2(order);
        silex::Ideal modulus(order);
        assert(primes.at(0)->get_ideal(p));
        assert(primes.at(1)->get_ideal(q));
        assert(ideal_power(p3, p, 3));
        assert(ideal_power(q2, q, 2));
        assert(modulus.multiply(p3, q2));
        silex::detail::ResidueUnitGroup group;
        assert(group.set_modulus(modulus));
        check_unit_group(group, modulus, {2, 2, 2});
    }
    return 0;
}

int test_unit_group_failures() {
    sflint::FmpqPoly polynomial;
    poly_x2_minus(polynomial, 5);
    silex::NumberField field = field_by_polynomial(polynomial);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::detail::ResidueUnitGroup group;
    assert(!group.is_defined());
    assert(group.generator_count() == 0);
    silex::OrderElement unused(maximal);
    sflint::FmpzMat empty(1, 0);
    assert(!group.discrete_log(sflint::FmpzMatRef(empty), unused));
    assert(!group.exp(unused, sflint::FmpzMatConstRef(empty)));

    silex::Ideal four(maximal);
    assert(set_rational_principal(four, 4));
    assert(group.set_modulus(four));

    // The equation order Z[sqrt 5] is not maximal: refused, group unchanged.
    silex::Ideal non_maximal(equation);
    assert(set_rational_principal(non_maximal, 3));
    assert(!group.set_modulus(non_maximal));
    check_unit_group(group, four, {2, 6});

    // A non-positive exponent is refused.
    sflint::Fmpz three;
    sflint::fmpz_set_si(three, 3);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal,
                                  sflint::FmpzConstRef(three)));
    assert(primes.size() == 1);
    assert(!group.set_prime_power(*primes.at(0), 0));
    check_unit_group(group, four, {2, 6});

    // An element of another order is refused.
    silex::OrderElement foreign(equation);
    assert(foreign.one());
    sflint::FmpzMat log(1, 2);
    assert(!group.discrete_log(sflint::FmpzMatRef(log), foreign));
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
    assert(test_unit_group_rational_moduli() == 0);
    assert(test_unit_group_prime_powers() == 0);
    assert(test_unit_group_failures() == 0);
    return 0;
}
