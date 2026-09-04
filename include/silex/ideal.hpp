#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

namespace silex {

struct DiagnosticsContext;
class Ideal;
class FractionalIdeal;

namespace detail {

bool set_known_two_generator_ideal(Ideal& out,
                                   flint::FmpzConstRef scalar_generator,
                                   const OrderElement& element_generator) noexcept;
bool multiply_integral_ideal_by_two_generator(
        Ideal& out,
        const Ideal& ideal,
        flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool multiply_integral_ideal_by_element(
        Ideal& out,
        const Ideal& ideal,
        const Element& multiplier,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool set_integral_ideal_known_hnf(FractionalIdeal& out,
                                  const Ideal& numerator) noexcept;

}  // namespace detail

class Ideal {
public:
    Ideal() noexcept = default;
    explicit Ideal(const Order& parent) noexcept;
    ~Ideal() noexcept;

    Ideal(const Ideal&) = delete;
    Ideal& operator=(const Ideal&) = delete;

    Ideal(Ideal&& other) noexcept;
    Ideal& operator=(Ideal&& other) noexcept;

    void swap(Ideal& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const Ideal& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong degree() const noexcept;
    bool has_hnf() const noexcept;

    bool one() noexcept;
    bool is_one() const noexcept;
    bool set_hnf(flint::FmpzMatConstRef hnf) noexcept;
    bool get_hnf(flint::FmpzMatRef out) const noexcept;
    bool set_principal(const OrderElement& generator) noexcept;
    const OrderElement* known_principal_generator() const noexcept;
    bool contains(const OrderElement& element) const noexcept;
    bool contains(const Ideal& ideal) const noexcept;
    bool norm(flint::FmpzRef out) const noexcept;
    bool equal(const Ideal& other) const noexcept;
    bool is_coprime(bool& result, const Ideal& other) const noexcept;
    bool add_to_one(OrderElement& left_witness,
                    OrderElement& right_witness,
                    const Ideal& other) const noexcept;

    bool add(const Ideal& left, const Ideal& right) noexcept;
    bool intersect(const Ideal& left, const Ideal& right) noexcept;
    bool multiply(const Ideal& left, const Ideal& right) noexcept;
    bool multiplier_ring(Order& out) const noexcept;

private:
    bool set_hnf_direct(flint::FmpzMatConstRef hnf) noexcept;
    bool set_rows(flint::FmpzMatConstRef rows) noexcept;
    bool set_known_ideal_rows(flint::FmpzMatConstRef rows) noexcept;

    friend class Order;
    friend bool detail::set_known_two_generator_ideal(
            Ideal& out,
            flint::FmpzConstRef scalar_generator,
            const OrderElement& element_generator) noexcept;
    friend bool detail::multiply_integral_ideal_by_two_generator(
            Ideal& out,
            const Ideal& ideal,
            flint::FmpzConstRef scalar_generator,
            const OrderElement& element_generator,
            const DiagnosticsContext* diagnostics) noexcept;
    friend bool detail::multiply_integral_ideal_by_element(
            Ideal& out,
            const Ideal& ideal,
            const Element& multiplier,
            const DiagnosticsContext* diagnostics) noexcept;

    Order parent_;
    flint::FmpzMat hnf_{0, 0};
    OrderElement principal_generator_;
    mutable flint::Fmpz norm_;
    bool has_hnf_ = false;
    bool has_principal_generator_ = false;
    mutable bool has_norm_ = false;
};

class FractionalIdeal {
public:
    FractionalIdeal() noexcept = default;
    explicit FractionalIdeal(const Order& parent) noexcept;
    ~FractionalIdeal() noexcept;

    FractionalIdeal(const FractionalIdeal&) = delete;
    FractionalIdeal& operator=(const FractionalIdeal&) = delete;

    FractionalIdeal(FractionalIdeal&& other) noexcept;
    FractionalIdeal& operator=(FractionalIdeal&& other) noexcept;

    void swap(FractionalIdeal& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const FractionalIdeal& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong degree() const noexcept;
    bool has_integral_denominator() const noexcept;

    bool one() noexcept;
    bool set_principal(const Element& generator,
                       const DiagnosticsContext* diagnostics = nullptr) noexcept;
    bool set_integral(const Ideal& numerator) noexcept;
    bool set_integral_den(const Ideal& numerator, flint::FmpzConstRef den) noexcept;
    bool get_integral_den(Ideal& numerator, flint::FmpzRef den) const noexcept;
    bool contains(const Element& element) const noexcept;
    bool norm(flint::FmpqRef out) const noexcept;
    bool equal(const FractionalIdeal& other) const noexcept;

    bool add(const FractionalIdeal& left,
             const FractionalIdeal& right) noexcept;
    bool intersect(const FractionalIdeal& left,
                   const FractionalIdeal& right) noexcept;
    bool multiply(const FractionalIdeal& left,
                  const FractionalIdeal& right) noexcept;
    bool pow_fmpz(const FractionalIdeal& input,
                  flint::FmpzConstRef exponent) noexcept;
    bool colon(const FractionalIdeal& numerator,
               const FractionalIdeal& denominator) noexcept;
    bool invert(const FractionalIdeal& input) noexcept;

private:
    friend bool detail::set_integral_ideal_known_hnf(
            FractionalIdeal& out,
            const Ideal& numerator) noexcept;

    Order parent_;
    Ideal numerator_;
    flint::Fmpz den_;
    mutable flint::Fmpq norm_;
    bool has_num_ = false;
    mutable bool has_norm_ = false;
};

inline void swap(Ideal& left, Ideal& right) noexcept {
    left.swap(right);
}

inline void swap(FractionalIdeal& left, FractionalIdeal& right) noexcept {
    left.swap(right);
}

}  // namespace silex
