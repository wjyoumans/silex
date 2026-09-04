#include <silex/element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

#include <cassert>
#include <iostream>

namespace {
namespace sflint = silex::flint;

slong integer_value(const sflint::Fmpq& value) {
    const sflint::FmpzConstRef denominator = sflint::fmpq_den_ref(value);
    assert(sflint::fmpz_equal_si(denominator, 1));
    const sflint::FmpzConstRef numerator = sflint::fmpq_num_ref(value);
    assert(sflint::fmpz_fits_si(numerator));
    return sflint::fmpz_get_si(numerator);
}

int run() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(radicand, 5);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Element theta(field);
    assert(theta.gen());

    silex::Element alpha(field);
    assert(alpha.add_si(theta, 3));

    sflint::Fmpq trace;
    sflint::Fmpq alpha_trace;
    sflint::Fmpq norm;
    assert(alpha.trace(sflint::FmpqRef(trace)));
    assert(alpha.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(trace, 6));
    assert(sflint::fmpq_equal_si(norm, 4));
    sflint::fmpq_set(alpha_trace, trace);

    silex::Element square(field);
    assert(square.multiply(alpha, alpha));
    assert(square.trace(sflint::FmpqRef(trace)));
    assert(square.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(trace, 28));
    assert(sflint::fmpq_equal_si(norm, 16));

    silex::Element inverse(field);
    silex::Element check(field);
    assert(inverse.invert(alpha));
    assert(check.multiply(alpha, inverse));
    assert(check.equal_si(1));

    sflint::Fmpz exponent;
    sflint::fmpz_set_si(exponent, 3);

    silex::Element cube(field);
    assert(cube.pow_fmpz(alpha, sflint::FmpzConstRef(exponent)));
    assert(cube.norm(sflint::FmpqRef(norm)));
    assert(sflint::fmpq_equal_si(norm, 64));

    std::cout << "K = Q(sqrt(5))\n";
    std::cout << "alpha = theta + 3\n";
    std::cout << "trace(alpha) = " << integer_value(alpha_trace) << "\n";
    std::cout << "norm(alpha^3) = " << integer_value(norm) << "\n";
    std::cout << "alpha * alpha^-1 = 1: " << check.equal_si(1) << "\n";
    return 0;
}

}  // namespace

int main() {
    return run();
}
