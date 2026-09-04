#pragma once

#include <silex/element.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

namespace silex {

class FieldHom {
public:
    FieldHom() noexcept = default;
    FieldHom(const NumberField& domain,
             const NumberField& codomain) noexcept;
    ~FieldHom() noexcept;

    FieldHom(const FieldHom&) = delete;
    FieldHom& operator=(const FieldHom&) = delete;

    FieldHom(FieldHom&& other) noexcept;
    FieldHom& operator=(FieldHom&& other) noexcept;

    void swap(FieldHom& other) noexcept;
    void clear() noexcept;
    bool define(const NumberField& domain,
                const NumberField& codomain) noexcept;
    bool set(const FieldHom& other) noexcept;

    bool is_defined() const noexcept;
    const NumberField* domain() const noexcept;
    const NumberField* codomain() const noexcept;
    bool has_generator_image() const noexcept;

    bool set_generator_image(const Element& image) noexcept;
    bool generator_image(Element& out) const noexcept;
    bool apply(Element& out, const Element& input) const noexcept;
    bool is_identity() const noexcept;
    bool is_isomorphism() const noexcept;

private:
    NumberField domain_;
    NumberField codomain_;
    Element generator_image_;
    bool has_generator_image_ = false;
};

inline void swap(FieldHom& left, FieldHom& right) noexcept {
    left.swap(right);
}

class OrderHom {
public:
    OrderHom() noexcept = default;
    OrderHom(const Order& source, const Order& target) noexcept;
    ~OrderHom() noexcept;

    OrderHom(const OrderHom&) = delete;
    OrderHom& operator=(const OrderHom&) = delete;

    OrderHom(OrderHom&& other) noexcept;
    OrderHom& operator=(OrderHom&& other) noexcept;

    void swap(OrderHom& other) noexcept;
    void clear() noexcept;
    bool define(const Order& source, const Order& target) noexcept;
    bool set(const OrderHom& other) noexcept;

    bool is_defined() const noexcept;
    const Order* source_order() const noexcept;
    const Order* target_order() const noexcept;
    bool has_field_homomorphism() const noexcept;

    bool set_field_homomorphism(const FieldHom& field_hom) noexcept;
    bool apply(Element& out, const Element& input) const noexcept;
    bool field_homomorphism(FieldHom& out) const noexcept;
    bool image_matrix(flint::FmpzMatRef out) const noexcept;
    bool image_matrix(flint::FmpzMat& out) const noexcept;

private:
    Order source_;
    Order target_;
    FieldHom field_hom_;
    flint::FmpzMat image_matrix_{0, 0};
    bool has_field_hom_ = false;
};

inline void swap(OrderHom& left, OrderHom& right) noexcept {
    left.swap(right);
}

}  // namespace silex
