#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>
#include <silex/residue_field.hpp>

#include <cassert>
#include <iostream>

namespace sflint = silex::flint;

namespace {

silex::NumberField rational_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    return silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

slong matrix_entry_si(const sflint::FmpzMat& matrix,
                      slong row,
                      slong col) noexcept {
    return sflint::fmpz_get_si(sflint::fmpz_mat_entry(matrix, row, col));
}

}  // namespace

int main() {
    silex::NumberField field = rational_field();
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());

    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 5);

    silex::FactorBase base(order);
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 3);
    assert(base.rational_prime_block_count() == 3);

    silex::PrimeIdeal prime5(order);
    sflint::Fmpz rational_prime;
    sflint::Fmpz prime_norm;
    assert(base.prime(prime5, 2));
    assert(prime5.rational_prime(sflint::FmpzRef(rational_prime)));
    assert(sflint::fmpz_equal_si(rational_prime, 5));
    assert(prime5.norm(sflint::FmpzRef(prime_norm)));
    assert(sflint::fmpz_equal_si(prime_norm, 5));

    silex::Element ten(field);
    assert(ten.set_si(10));

    silex::Relation relation(base);
    sflint::FmpzMat row(1, base.length());
    assert(relation.set_generator(ten));
    assert(relation.exponents(sflint::FmpzMatRef(row)));
    assert(matrix_entry_si(row, 0, 0) == 1);
    assert(matrix_entry_si(row, 0, 1) == 0);
    assert(matrix_entry_si(row, 0, 2) == 1);

    silex::RelationMatrix relation_matrix(base);
    assert(relation_matrix.append(relation));
    assert(relation_matrix.length() == 1);

    silex::ResidueField residue_field(prime5);
    auto cardinality = residue_field.cardinality();
    assert(cardinality.has_value());
    assert(sflint::fmpz_equal_si(*cardinality, 5));

    silex::OrderElement two(order);
    silex::OrderElement three(order);
    assert(two.set_si(2));
    assert(three.set_si(3));

    silex::ResidueFieldElement two_mod_p(residue_field);
    silex::ResidueFieldElement three_mod_p(residue_field);
    silex::ResidueFieldElement product(residue_field);
    silex::ResidueFieldElement one(residue_field);
    assert(two_mod_p.set_order_element(two));
    assert(three_mod_p.set_order_element(three));
    assert(product.multiply(two_mod_p, three_mod_p));
    assert(one.one());
    assert(product.equal(one));

    std::cout << "K = Q, O = Z\n";
    std::cout << "factor base bound = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(bound)) << "\n";
    std::cout << "factor base length = " << base.length() << "\n";
    std::cout << "selected prime norm = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(prime_norm))
              << "\n";
    std::cout << "relation row for 10 = ["
              << matrix_entry_si(row, 0, 0) << ", "
              << matrix_entry_si(row, 0, 1) << ", "
              << matrix_entry_si(row, 0, 2) << "]\n";
    std::cout << "residue field cardinality = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(*cardinality))
              << "\n";
    std::cout << "2 * 3 == 1 mod 5: " << product.equal(one) << "\n";
    return 0;
}
