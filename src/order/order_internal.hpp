#pragma once

#include <silex/order.hpp>

namespace silex::detail {

inline bool order_has_parented_basis(const Order& order) noexcept {
    return order.is_defined() && order.has_basis() && order.parent() != nullptr;
}

// Noninstalled maximality hook.  Public code cannot assert or withdraw
// maximality: Order records it only where it is computed (maximal_order,
// quadratic conductor metadata, degree one).  This hook writes the shared flag
// without any check, so every handle sharing the order's data sees the change.
// It exists for tests that must build states the public API cannot reach, such
// as an order whose maximality is withdrawn after dependent objects were
// built; library code must not use it to label an order maximal.
class OrderAccess {
public:
    static void set_maximality_unchecked(Order& order,
                                         bool is_maximal) noexcept {
        order.record_maximality(is_maximal);
    }
};

}  // namespace silex::detail
