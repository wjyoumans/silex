#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

#include <optional>

namespace silex {

class ResidueElement;

class ResidueRing {
public:
    ResidueRing() noexcept = default;
    explicit ResidueRing(const Ideal& modulus) noexcept;
    ~ResidueRing() noexcept;

    ResidueRing(const ResidueRing&) = delete;
    ResidueRing& operator=(const ResidueRing&) = delete;

    ResidueRing(ResidueRing&& other) noexcept;
    ResidueRing& operator=(ResidueRing&& other) noexcept;

    void swap(ResidueRing& other) noexcept;
    void clear() noexcept;
    bool define(const Ideal& modulus) noexcept;
    bool set(const ResidueRing& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent_order() const noexcept;
    const Ideal* modulus() const noexcept;
    slong degree() const noexcept;

    bool get_modulus(Ideal& out) const noexcept;
    bool cardinality(flint::FmpzRef out) const noexcept;
    std::optional<flint::Fmpz> cardinality() const noexcept;
    bool equal(const ResidueRing& other) const noexcept;

private:
    bool reduce_row(flint::FmpzMatRef row) const noexcept;
    bool multiply_rows(flint::FmpzMatRef out,
                       flint::FmpzMatConstRef left,
                       flint::FmpzMatConstRef right) const noexcept;

    friend class ResidueElement;

    Ideal modulus_;
    flint::FmpzMat multiplication_table_{0, 0};
};

class ResidueElement {
public:
    ResidueElement() noexcept = default;
    explicit ResidueElement(const ResidueRing& parent) noexcept;
    ~ResidueElement() noexcept;

    ResidueElement(const ResidueElement&) = delete;
    ResidueElement& operator=(const ResidueElement&) = delete;

    ResidueElement(ResidueElement&& other) noexcept;
    ResidueElement& operator=(ResidueElement&& other) noexcept;

    void swap(ResidueElement& other) noexcept;
    void clear() noexcept;
    bool define(const ResidueRing& parent) noexcept;
    bool set(const ResidueElement& other) noexcept;

    bool is_defined() const noexcept;
    const ResidueRing* parent() const noexcept;
    slong degree() const noexcept;

    bool zero() noexcept;
    bool one() noexcept;
    bool set_order_element(const OrderElement& element) noexcept;
    bool set_element(const Element& element) noexcept;
    bool get_coordinates(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::FmpzMat> coordinates() const noexcept;
    bool lift(OrderElement& out) const noexcept;
    bool lift(Element& out) const noexcept;

    bool equal(const ResidueElement& other) const noexcept;
    bool add(const ResidueElement& left, const ResidueElement& right) noexcept;
    bool negate(const ResidueElement& input) noexcept;
    bool subtract(const ResidueElement& left, const ResidueElement& right) noexcept;
    bool multiply(const ResidueElement& left, const ResidueElement& right) noexcept;
    bool crt(const ResidueElement& left, const ResidueElement& right) noexcept;

private:
    ResidueRing parent_;
    flint::FmpzMat coordinates_{0, 0};
};

inline void swap(ResidueRing& left, ResidueRing& right) noexcept {
    left.swap(right);
}

inline void swap(ResidueElement& left, ResidueElement& right) noexcept {
    left.swap(right);
}

}  // namespace silex
