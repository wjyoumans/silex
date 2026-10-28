#pragma once

#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order_element.hpp>

namespace silex::detail {

class OrderElementAccess {
public:
    // Sets `out` to `element` with the 1 x n integral `coordinates` it
    // already has in the basis of `out`'s order, copying both instead of
    // recomputing one from the other (set_coordinates multiplies by the
    // basis; set_element solves for the coordinates).  The caller
    // guarantees that `coordinates` are the coordinates of `element`; only
    // the shape and the parent field are checked, except that a build with
    // debug checks verifies the guarantee when `diagnostics` enables the
    // order module's expensive checks, reporting a violation to the
    // debug-failure callback and returning false.
    static bool set_element_with_coordinates(
            OrderElement& out,
            const Element& element,
            flint::FmpzMatConstRef coordinates,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

    // The stored 1 x n coordinates of `element`, or nullptr when it is not
    // defined.
    static const fmpz_mat_struct* coordinates(
            const OrderElement& element) noexcept {
        return element.is_defined() ? element.coords_.raw() : nullptr;
    }
};

}  // namespace silex::detail
