#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/element.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

bool qcoeff_equal_si(sflint::FmpqPolyConstRef polynomial,
                     slong i,
                     slong expected) noexcept {
    sflint::Fmpq coeff;
    sflint::fmpq_poly_get_coeff_fmpq(
            sflint::FmpqRef(coeff), polynomial, i);
    return sflint::fmpq_equal_si(coeff, expected);
}

void set_quadratic_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -5);
}

bool set_one_plus_theta(silex::Element& element) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    return element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

}  // namespace

int main() {
    silex::NumberField empty;
    assert(!empty.is_defined());
    assert(empty.degree() == 0);
    assert(empty.backend_kind() == silex::NumberFieldBackendKind::generic);
    assert(empty.raw_flint_field() == nullptr);

    sflint::FmpqPoly polynomial;
    set_quadratic_polynomial(polynomial);
    sflint::Fmpz stored_radicand;

    silex::NumberField from_rational_polynomial;
    assert(from_rational_polynomial.define_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial)));
    assert(from_rational_polynomial.is_defined());
    assert(from_rational_polynomial.degree() == 2);
    assert(from_rational_polynomial.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    assert(from_rational_polynomial.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));
    assert(from_rational_polynomial.raw_flint_field() != nullptr);
    assert(sflint::fmpq_poly_degree(sflint::FmpqPolyConstRef(
            from_rational_polynomial.raw_flint_field()->pol)) == 2);
    assert(qcoeff_equal_si(
            sflint::FmpqPolyConstRef(
                    from_rational_polynomial.raw_flint_field()->pol),
            0,
            -5));
    assert(qcoeff_equal_si(
            sflint::FmpqPolyConstRef(
                    from_rational_polynomial.raw_flint_field()->pol),
            2,
            1));
    assert(fmpq_poly_equal(from_rational_polynomial.raw_flint_field()->pol,
                           polynomial.raw()) != 0);

    silex::Element polynomial_generator(from_rational_polynomial);
    sflint::FmpqPoly polynomial_generator_value;
    assert(polynomial_generator.gen());
    assert(polynomial_generator.get_fmpq_poly(
            sflint::FmpqPolyRef(polynomial_generator_value)));
    assert(sflint::fmpq_poly_degree(polynomial_generator_value) == 1);
    assert(qcoeff_equal_si(
            sflint::FmpqPolyConstRef(polynomial_generator_value), 0, 0));
    assert(qcoeff_equal_si(
            sflint::FmpqPolyConstRef(polynomial_generator_value), 1, 1));

    silex::NumberField factory_rational = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    assert(factory_rational.is_defined());
    assert(factory_rational.degree() == 2);
    assert(factory_rational.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);

    sflint::FmpqPoly constant_polynomial;
    sflint::fmpq_poly_set_coeff_si(constant_polynomial, 0, 7);
    assert(!from_rational_polynomial.define_by_polynomial(
        sflint::FmpqPolyConstRef(constant_polynomial)));
    assert(from_rational_polynomial.is_defined());
    assert(from_rational_polynomial.degree() == 2);
    assert(!silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(constant_polynomial)).is_defined());

    sflint::FmpzPoly integral_polynomial;
    sflint::fmpz_poly_set_coeff_si(integral_polynomial, 2, 1);
    sflint::fmpz_poly_set_coeff_si(integral_polynomial, 0, -5);

    silex::NumberField from_integral_polynomial;
    assert(from_integral_polynomial.define_by_polynomial(
        sflint::FmpzPolyConstRef(integral_polynomial)));
    assert(from_integral_polynomial.degree() == 2);
    assert(from_integral_polynomial.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    assert(from_integral_polynomial.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));
    silex::NumberField factory_integral = silex::NumberField::by_polynomial(
            sflint::FmpzPolyConstRef(integral_polynomial));
    assert(factory_integral.is_defined());
    assert(factory_integral.degree() == 2);
    assert(factory_integral.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);

    sflint::FmpqPoly shifted_polynomial;
    sflint::fmpq_poly_set_coeff_si(shifted_polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(shifted_polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(shifted_polynomial, 0, 1);
    silex::NumberField shifted = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(shifted_polynomial));
    assert(shifted.is_defined());
    assert(shifted.backend_kind() == silex::NumberFieldBackendKind::generic);
    assert(!shifted.quadratic_radicand(stored_radicand));

    sflint::FmpqPoly nonsquarefree_polynomial;
    sflint::fmpq_poly_set_coeff_si(nonsquarefree_polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(nonsquarefree_polynomial, 0, -12);
    silex::NumberField nonsquarefree = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(nonsquarefree_polynomial));
    assert(nonsquarefree.is_defined());
    assert(nonsquarefree.backend_kind() ==
            silex::NumberFieldBackendKind::generic);
    assert(!nonsquarefree.quadratic_radicand(stored_radicand));

    sflint::FmpqPoly nonmonic_polynomial;
    sflint::fmpq_poly_set_coeff_si(nonmonic_polynomial, 2, 2);
    sflint::fmpq_poly_set_coeff_si(nonmonic_polynomial, 0, -10);
    silex::NumberField nonmonic = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(nonmonic_polynomial));
    assert(nonmonic.is_defined());
    assert(nonmonic.backend_kind() == silex::NumberFieldBackendKind::generic);
    assert(!nonmonic.quadratic_radicand(stored_radicand));

    sflint::FmpqPoly nonintegral_polynomial;
    sflint::Fmpq nonintegral_constant;
    sflint::fmpq_poly_set_coeff_si(nonintegral_polynomial, 2, 1);
    sflint::fmpq_set_si(nonintegral_constant, -5, 2);
    sflint::fmpq_poly_set_coeff_fmpq(
            nonintegral_polynomial, 0, nonintegral_constant);
    silex::NumberField nonintegral = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(nonintegral_polynomial));
    assert(nonintegral.is_defined());
    assert(nonintegral.backend_kind() ==
            silex::NumberFieldBackendKind::generic);
    assert(!nonintegral.quadratic_radicand(stored_radicand));

    sflint::FmpqPoly imaginary_polynomial;
    sflint::fmpq_poly_set_coeff_si(imaginary_polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(imaginary_polynomial, 0, 47);
    silex::NumberField imaginary = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(imaginary_polynomial));
    assert(imaginary.is_defined());
    assert(imaginary.backend_kind() ==
            silex::NumberFieldBackendKind::quadratic);
    assert(imaginary.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, -47));

    silex::NumberField copied;
    assert(copied.set(from_integral_polynomial));
    from_integral_polynomial.clear();
    assert(copied.is_defined());
    assert(copied.degree() == 2);
    assert(qcoeff_equal_si(sflint::FmpqPolyConstRef(
            copied.raw_flint_field()->pol), 0, -5));

    sflint::Fmpz radicand;
    sflint::fmpz_set_si(radicand, 5);

    silex::NumberField quadratic;
    assert(quadratic.define_quadratic(radicand));
    assert(quadratic.degree() == 2);
    assert(quadratic.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    assert(qcoeff_equal_si(sflint::FmpqPolyConstRef(
            quadratic.raw_flint_field()->pol), 0, -5));
    assert(qcoeff_equal_si(sflint::FmpqPolyConstRef(
            quadratic.raw_flint_field()->pol), 2, 1));
    silex::NumberField factory_quadratic =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(factory_quadratic.is_defined());
    assert(factory_quadratic.degree() == 2);
    assert(factory_quadratic.backend_kind() ==
           silex::NumberFieldBackendKind::quadratic);

    assert(quadratic.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));

    silex::NumberField moved(std::move(quadratic));
    assert(moved.is_defined());
    assert(!quadratic.is_defined());
    assert(moved.backend_kind() == silex::NumberFieldBackendKind::quadratic);

    silex::NumberField assigned;
    assigned = std::move(moved);
    assert(assigned.is_defined());
    assert(!moved.is_defined());
    assert(assigned.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));
    assert(sflint::fmpq_poly_degree(sflint::FmpqPolyConstRef(
            assigned.flint_field_ref().raw()->pol)) == 2);

    silex::NumberField swapped;
    swap(swapped, assigned);
    assert(swapped.is_defined());
    assert(swapped.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    assert(swapped.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));
    assert(!assigned.is_defined());

    assert(assigned.set(swapped));
    assert(assigned.is_defined());
    assert(assigned.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    assert(assigned.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));

    silex::NumberField undefined;
    assert(assigned.set(undefined));
    assert(!assigned.is_defined());
    assert(assigned.raw_flint_field() == nullptr);

    assert(assigned.set(swapped));
    assert(assigned.is_defined());
    assert(assigned.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    assert(assigned.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));

    sflint::fmpz_set_si(radicand, 0);
    assert(!assigned.define_quadratic(radicand));
    assert(!silex::NumberField::quadratic(
            sflint::FmpzConstRef(radicand)).is_defined());
    assert(assigned.backend_kind() == silex::NumberFieldBackendKind::quadratic);
    assert(assigned.quadratic_radicand(stored_radicand));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));

    sflint::fmpz_set_si(radicand, 4);
    assert(!assigned.define_quadratic(radicand));
    assert(assigned.backend_kind() == silex::NumberFieldBackendKind::quadratic);

    sflint::fmpz_set_si(radicand, 12);
    assert(!assigned.define_quadratic(radicand));
    assert(assigned.backend_kind() == silex::NumberFieldBackendKind::quadratic);

    silex::NumberField quadratic_paths =
            silex::NumberField::quadratic(sflint::FmpzConstRef(stored_radicand));
    assert(quadratic_paths.is_defined());

    silex::Order equation_order =
            silex::Order::equation_order(quadratic_paths);
    assert(equation_order.is_defined());
    silex::Order maximal_order(quadratic_paths);
    assert(maximal_order.maximal_order(equation_order));
    assert(maximal_order.is_maximal());
    assert(maximal_order.discriminant(sflint::FmpzRef(stored_radicand)));
    assert(sflint::fmpz_equal_si(stored_radicand, 5));

    sflint::Fmpz p;
    sflint::fmpz_set_si(p, 11);
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, maximal_order,
                                  sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);

    const silex::PrimeIdeal* first_prime = primes.at(0);
    assert(first_prime != nullptr);
    silex::Ideal prime_ideal(maximal_order);
    assert(first_prime->get_ideal(prime_ideal));
    assert(prime_ideal.norm(sflint::FmpzRef(stored_radicand)));
    assert(sflint::fmpz_equal_si(stored_radicand, 11));

    silex::Element generator(quadratic_paths);
    assert(set_one_plus_theta(generator));
    silex::OrderElement order_generator(maximal_order);
    assert(order_generator.set_element(generator));

    silex::Ideal principal(maximal_order);
    assert(principal.set_principal(order_generator));
    assert(principal.norm(sflint::FmpzRef(stored_radicand)));
    assert(sflint::fmpz_equal_si(stored_radicand, 4));

    sflint::FmpzMat hnf(2, 2);
    assert(principal.get_hnf(sflint::FmpzMatRef(hnf)));
    sflint::Fmpz det;
    fmpz_mat_det(det.raw(), hnf.raw());
    sflint::fmpz_abs(sflint::FmpzRef(det), sflint::FmpzConstRef(det));
    assert(sflint::fmpz_equal(det, stored_radicand));

    return 0;
}
