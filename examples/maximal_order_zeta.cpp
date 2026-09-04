#include <silex/order.hpp>
#include <silex/zeta.hpp>

#include <silex/flint/fmpq_poly.hpp>

#include <cassert>
#include <iostream>

namespace sflint = silex::flint;

int main() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);

    silex::NumberField field =
            silex::NumberField::by_polynomial(
                    sflint::FmpqPolyConstRef(polynomial));
    assert(field.is_defined());

    silex::Order maximal_order(field);
    silex::Order equation_order = silex::Order::equation_order(field);
    assert(equation_order.is_defined());
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());

    sflint::Fmpz index;
    sflint::Fmpz discriminant;
    assert(silex::order_index(sflint::FmpzRef(index), equation_order,
                              maximal_order));
    assert(maximal_order.discriminant(sflint::FmpzRef(discriminant)));

    auto audit = silex::zeta_class_regulator_product_bf_audit(
            maximal_order, 20000, 128);
    assert(audit.has_value());

    std::cout << "K = Q(theta), theta^3 - theta - 1 = 0\n";
    std::cout << "degree = " << maximal_order.degree() << "\n";
    std::cout << "[O_K : Z[theta]] = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(index)) << "\n";
    std::cout << "disc(O_K) = "
              << sflint::fmpz_get_si(sflint::FmpzConstRef(discriminant))
              << "\n";
    std::cout << "BF cutoff = " << audit->cutoff << "\n";
    std::cout << "BF work precision = " << audit->work_precision << "\n";
    std::cout << "class-regulator product positive = "
              << sflint::arb_is_positive(audit->value) << "\n";
    std::cout << "BF error bound positive = "
              << sflint::arb_is_positive(audit->error_bound) << "\n";
    return 0;
}
