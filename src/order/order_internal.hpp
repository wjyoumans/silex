#pragma once

#include <silex/order.hpp>

namespace silex::detail {

inline bool order_has_parented_basis(const Order& order) noexcept {
    return order.is_defined() && order.has_basis() && order.parent() != nullptr;
}

}  // namespace silex::detail
