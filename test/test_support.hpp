#pragma once

#include <flint/flint.h>

#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include <cassert>

namespace silex::test {

inline NumberField field_by_polynomial(
        flint::FmpqPolyConstRef polynomial) noexcept {
    NumberField field = NumberField::by_polynomial(polynomial);
    assert(field.is_defined());
    return field;
}

inline NumberField field_by_polynomial(
        flint::FmpzPolyConstRef polynomial) noexcept {
    NumberField field = NumberField::by_polynomial(polynomial);
    assert(field.is_defined());
    return field;
}

inline NumberField quadratic_field(slong radicand) noexcept {
    flint::Fmpz d;
    flint::fmpz_set_si(flint::FmpzRef(d), radicand);
    NumberField field = NumberField::quadratic(flint::FmpzConstRef(d));
    assert(field.is_defined());
    return field;
}

inline Order equation_order(const NumberField& parent) noexcept {
    Order order = Order::equation_order(parent);
    assert(order.is_defined());
    return order;
}

}  // namespace silex::test
