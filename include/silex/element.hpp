#pragma once

#include <flint/flint.h>
#include <flint/nf_elem.h>

#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/nf_elem.hpp>
#include <silex/number_field.hpp>

#include <optional>

namespace silex {

struct DiagnosticsContext;

class Element {
public:
    Element() noexcept = default;
    explicit Element(const NumberField& parent) noexcept;
    ~Element() noexcept;

    Element(const Element&) = delete;
    Element& operator=(const Element&) = delete;

    Element(Element&& other) noexcept;
    Element& operator=(Element&& other) noexcept;

    void swap(Element& other) noexcept;
    void clear() noexcept;
    bool define(const NumberField& parent) noexcept;
    bool set(const Element& other) noexcept;

    bool is_defined() const noexcept;
    const NumberField* parent() const noexcept;
    bool has_parent(const NumberField& parent) const noexcept;
    bool has_same_parent(const Element& other) const noexcept;

    bool zero() noexcept;
    bool one() noexcept;
    bool gen() noexcept;
    bool set_si(slong value) noexcept;
    bool set_fmpz(flint::FmpzConstRef value) noexcept;
    bool set_si_over_si(slong numerator, slong denominator) noexcept;
    bool set_fmpq_poly(flint::FmpqPolyConstRef polynomial) noexcept;
    bool get_fmpq_poly(flint::FmpqPolyRef polynomial) const noexcept;
    std::optional<flint::FmpqPoly> to_fmpq_poly() const noexcept;
    bool equal(const Element& other) const noexcept;
    bool equal_si(slong value) const noexcept;
    bool negate(const Element& input) noexcept;
    bool add(const Element& left, const Element& right) noexcept;
    bool add_si(const Element& input, slong value) noexcept;
    bool subtract(const Element& left, const Element& right) noexcept;
    bool multiply(const Element& left, const Element& right) noexcept;
    bool scalar_div_si(const Element& input, slong denominator) noexcept;
    bool invert(const Element& input) noexcept;
    bool pow_fmpz(const Element& input, flint::FmpzConstRef exponent) noexcept;
    bool is_square(bool& is_square,
                   Element& root,
                   const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;
    bool is_power(bool& is_power,
                  Element& root,
                  flint::FmpzConstRef exponent,
                  const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;
    bool trace(flint::FmpqRef out) const noexcept;
    bool norm(flint::FmpqRef out) const noexcept;
    bool conjugate(Element& out) const noexcept;

    // Low-level FLINT interop for bridge code, parity tests, and direct FLINT
    // call sites. Ordinary users should prefer Element domain operations.
    flint::NfElemRef flint_element_ref() noexcept;
    flint::NfElemConstRef flint_element_ref() const noexcept;
    nf_elem_struct* raw_flint_element() noexcept;
    const nf_elem_struct* raw_flint_element() const noexcept;

private:
    NumberField parent_;
    flint::NfElem value_;
};

inline void swap(Element& left, Element& right) noexcept {
    left.swap(right);
}

}  // namespace silex
