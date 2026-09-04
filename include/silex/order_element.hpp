#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order.hpp>

#include <optional>

namespace silex {

class OrderElement {
public:
    OrderElement() noexcept = default;
    explicit OrderElement(const Order& parent) noexcept;
    ~OrderElement() noexcept;

    OrderElement(const OrderElement&) = delete;
    OrderElement& operator=(const OrderElement&) = delete;

    OrderElement(OrderElement&& other) noexcept;
    OrderElement& operator=(OrderElement&& other) noexcept;

    void swap(OrderElement& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const OrderElement& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong degree() const noexcept;

    bool zero() noexcept;
    bool one() noexcept;
    bool set_si(slong value) noexcept;
    bool set_element(const Element& element) noexcept;
    bool set_coordinates(flint::FmpzMatConstRef coordinates) noexcept;

    bool get_element(Element& out) const noexcept;
    bool get_coordinates(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::FmpzMat> coordinates() const noexcept;

    bool negate(const OrderElement& input) noexcept;
    bool add(const OrderElement& left, const OrderElement& right) noexcept;
    bool subtract(const OrderElement& left, const OrderElement& right) noexcept;
    bool multiply(const OrderElement& left, const OrderElement& right) noexcept;
    bool equal(const OrderElement& other) const noexcept;
    bool equal_si(slong value) const noexcept;

private:
    Order parent_;
    Element value_;
    flint::FmpzMat coords_{0, 0};
};

inline void swap(OrderElement& left, OrderElement& right) noexcept {
    left.swap(right);
}

}  // namespace silex
