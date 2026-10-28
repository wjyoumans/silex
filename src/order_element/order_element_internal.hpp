#pragma once

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
    // the shape and the parent field are checked.
    static bool set_element_with_coordinates(
            OrderElement& out,
            const Element& element,
            flint::FmpzMatConstRef coordinates) noexcept;
};

}  // namespace silex::detail
