#pragma once

#include <flint/flint.h>
#include <flint/nf.h>

#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/nf.hpp>

#include <memory>

namespace silex {

namespace detail {
struct NumberFieldData;
}

enum class NumberFieldBackendKind {
    generic = 0,
    quadratic = 1,
};

class NumberField {
public:
    NumberField() noexcept = default;
    ~NumberField() noexcept;

    NumberField(const NumberField& other) noexcept;
    NumberField& operator=(const NumberField& other) noexcept;

    NumberField(NumberField&& other) noexcept;
    NumberField& operator=(NumberField&& other) noexcept;

    void swap(NumberField& other) noexcept;
    void clear() noexcept;
    bool set(const NumberField& other) noexcept;

    // Preferred public construction. Invalid input returns an undefined field.
    // Exact monic integral x^2-d inputs with squarefree nonsquare d retain
    // their supplied polynomial/generator and use the quadratic backend.
    static NumberField by_polynomial(
            flint::FmpqPolyConstRef polynomial) noexcept;
    static NumberField by_polynomial(
            flint::FmpzPolyConstRef polynomial) noexcept;
    static NumberField quadratic(flint::FmpzConstRef radicand) noexcept;

    // Compatibility/scratch-object construction helpers. Prefer factories in
    // ordinary user code unless mutation/failure preservation is required.
    bool define_by_polynomial(flint::FmpqPolyConstRef polynomial) noexcept;
    bool define_by_polynomial(flint::FmpzPolyConstRef polynomial) noexcept;
    bool define_quadratic(flint::FmpzConstRef radicand) noexcept;

    bool is_defined() const noexcept;
    bool has_same_data(const NumberField& other) const noexcept;
    slong degree() const noexcept;
    NumberFieldBackendKind backend_kind() const noexcept;
    bool quadratic_radicand(flint::FmpzRef out) const noexcept;

    // Low-level FLINT interop for bridge code, parity tests, and direct FLINT
    // call sites. Ordinary users should prefer NumberField domain operations.
    flint::NfRef flint_field_ref() noexcept;
    flint::NfConstRef flint_field_ref() const noexcept;
    nf_struct* raw_flint_field() noexcept;
    const nf_struct* raw_flint_field() const noexcept;

private:
    std::shared_ptr<detail::NumberFieldData> data_;
};

inline void swap(NumberField& left, NumberField& right) noexcept {
    left.swap(right);
}

}  // namespace silex
