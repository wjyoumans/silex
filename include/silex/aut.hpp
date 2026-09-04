#pragma once

#include <silex/hom.hpp>
#include <silex/number_field.hpp>

namespace silex {

class FieldAutomorphism {
public:
    FieldAutomorphism() noexcept = default;
    explicit FieldAutomorphism(const NumberField& parent) noexcept;
    ~FieldAutomorphism() noexcept = default;

    FieldAutomorphism(const FieldAutomorphism&) = delete;
    FieldAutomorphism& operator=(const FieldAutomorphism&) = delete;

    FieldAutomorphism(FieldAutomorphism&& other) noexcept;
    FieldAutomorphism& operator=(FieldAutomorphism&& other) noexcept;

    void swap(FieldAutomorphism& other) noexcept;
    void clear() noexcept;
    bool define(const NumberField& parent) noexcept;
    bool set(const FieldAutomorphism& other) noexcept;

    bool is_defined() const noexcept;
    const NumberField* parent() const noexcept;
    bool has_homomorphism() const noexcept;

    bool set_identity() noexcept;
    bool set_quadratic_conjugation() noexcept;
    bool apply(Element& out, const Element& input) const noexcept;
    bool homomorphism(FieldHom& out) const noexcept;
    bool is_identity() const noexcept;

private:
    NumberField parent_;
    FieldHom hom_;
    bool has_homomorphism_ = false;
};

inline void swap(FieldAutomorphism& left,
                 FieldAutomorphism& right) noexcept {
    left.swap(right);
}

}  // namespace silex
