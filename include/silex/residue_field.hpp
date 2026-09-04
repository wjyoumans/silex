#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>

#include <optional>

namespace silex {

class ResidueFieldElement;
class ResidueFieldQuotientLog;

class ResidueField {
public:
    ResidueField() noexcept = default;
    explicit ResidueField(const PrimeIdeal& prime) noexcept;
    ~ResidueField() noexcept;

    ResidueField(const ResidueField&) = delete;
    ResidueField& operator=(const ResidueField&) = delete;

    ResidueField(ResidueField&& other) noexcept;
    ResidueField& operator=(ResidueField&& other) noexcept;

    void swap(ResidueField& other) noexcept;
    void clear() noexcept;
    bool set(const ResidueField& other) noexcept;
    bool set_prime(const PrimeIdeal& prime) noexcept;

    bool is_defined() const noexcept;
    const Order* parent_order() const noexcept;
    const PrimeIdeal* prime() const noexcept;
    slong degree() const noexcept;

    bool get_prime(PrimeIdeal& out) const noexcept;
    bool characteristic(flint::FmpzRef out) const noexcept;
    bool cardinality(flint::FmpzRef out) const noexcept;
    bool modulus(flint::FmpzPolyRef out) const noexcept;
    std::optional<flint::Fmpz> characteristic() const noexcept;
    std::optional<flint::Fmpz> cardinality() const noexcept;
    std::optional<flint::FmpzPoly> modulus() const noexcept;
    bool equal(const ResidueField& other) const noexcept;
    bool multiplicative_generator(ResidueFieldElement& out) const noexcept;

private:
    explicit ResidueField(flint::FmpzConstRef characteristic) noexcept;

    PrimeIdeal prime_;
    flint::Fmpz p_;
    flint::FmpzModCtx ctx_;
    flint::FmpzModPoly modulus_;
    bool defined_ = false;

    friend class ResidueFieldElement;
    friend class ResidueFieldQuotientLog;
};

class ResidueFieldElement {
public:
    ResidueFieldElement() noexcept = default;
    explicit ResidueFieldElement(const ResidueField& parent) noexcept;
    ~ResidueFieldElement() noexcept;

    ResidueFieldElement(const ResidueFieldElement&) = delete;
    ResidueFieldElement& operator=(const ResidueFieldElement&) = delete;

    ResidueFieldElement(ResidueFieldElement&& other) noexcept;
    ResidueFieldElement& operator=(ResidueFieldElement&& other) noexcept;

    void swap(ResidueFieldElement& other) noexcept;
    void clear() noexcept;
    bool define(const ResidueField& parent) noexcept;
    bool set(const ResidueFieldElement& other) noexcept;

    bool is_defined() const noexcept;
    const ResidueField* parent() const noexcept;

    bool zero() noexcept;
    bool one() noexcept;
    bool is_zero() const noexcept;
    bool set_polynomial(flint::FmpzPolyConstRef polynomial) noexcept;
    bool set_order_element(const OrderElement& element) noexcept;
    bool set_element(const Element& element) noexcept;
    bool set_factored_element(const FactoredElement& element) noexcept;
    bool get_polynomial(flint::FmpzPolyRef out) const noexcept;
    std::optional<flint::FmpzPoly> polynomial() const noexcept;
    bool degree_one_scalar(flint::FmpzRef out) const noexcept;

    bool equal(const ResidueFieldElement& other) const noexcept;
    bool add(const ResidueFieldElement& left,
             const ResidueFieldElement& right) noexcept;
    bool negate(const ResidueFieldElement& input) noexcept;
    bool subtract(const ResidueFieldElement& left,
                  const ResidueFieldElement& right) noexcept;
    bool multiply(const ResidueFieldElement& left,
                  const ResidueFieldElement& right) noexcept;
    bool invert(const ResidueFieldElement& input) noexcept;
    bool pow_fmpz(const ResidueFieldElement& input,
                  flint::FmpzConstRef exponent) noexcept;
    bool multiplicative_order(flint::FmpzRef out) const noexcept;
    bool discrete_log(flint::FmpzRef out,
                      const ResidueFieldElement& base) const noexcept;
    bool quotient_log_mod_prime(flint::FmpzRef out,
                                flint::FmpzConstRef ell) const noexcept;

private:
    ResidueField parent_;
    flint::FmpzModPoly rep_;

    friend class ResidueFieldQuotientLog;
};

class ResidueFieldQuotientLog {
public:
    ResidueFieldQuotientLog() noexcept = default;
    explicit ResidueFieldQuotientLog(const ResidueField& parent) noexcept;
    ~ResidueFieldQuotientLog() noexcept;

    ResidueFieldQuotientLog(const ResidueFieldQuotientLog&) = delete;
    ResidueFieldQuotientLog& operator=(const ResidueFieldQuotientLog&) = delete;

    ResidueFieldQuotientLog(ResidueFieldQuotientLog&& other) noexcept;
    ResidueFieldQuotientLog& operator=(
            ResidueFieldQuotientLog&& other) noexcept;

    void swap(ResidueFieldQuotientLog& other) noexcept;
    void clear() noexcept;
    bool define(const ResidueField& parent) noexcept;
    bool set(const ResidueFieldQuotientLog& other) noexcept;

    bool is_defined() const noexcept;
    bool is_set() const noexcept;
    const ResidueField* parent() const noexcept;

    bool set_ell(flint::FmpzConstRef ell) noexcept;
    bool ell(flint::FmpzRef out) const noexcept;
    bool generator(ResidueFieldElement& out) const noexcept;
    bool quotient_generator(ResidueFieldElement& out) const noexcept;
    bool apply(flint::FmpzRef out,
               const ResidueFieldElement& element) const noexcept;

private:
    static bool degree_one_scalar_(flint::Fmpz& out,
                                   const ResidueFieldElement& element) noexcept;

    ResidueField parent_;
    flint::Fmpz ell_;
    flint::Fmpz order_;
    flint::Fmpz cofactor_;
    ResidueFieldElement generator_;
    ResidueFieldElement quotient_generator_;
    bool is_set_ = false;
};

inline void swap(ResidueField& left, ResidueField& right) noexcept {
    left.swap(right);
}

inline void swap(ResidueFieldElement& left, ResidueFieldElement& right) noexcept {
    left.swap(right);
}

inline void swap(ResidueFieldQuotientLog& left,
                 ResidueFieldQuotientLog& right) noexcept {
    left.swap(right);
}

}  // namespace silex
