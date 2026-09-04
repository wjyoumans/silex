#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

#include <memory>

namespace silex {

struct DiagnosticsContext;
class FactoredElement;
class PrimeIdeal;
class PrimeIdealList;

namespace detail {

class MaximalQuadraticPrimeAccess;
bool set_degree_one_prime_ideal_from_root(PrimeIdeal& out,
                                          const Order& order,
                                          flint::FmpzConstRef p,
                                          flint::FmpzConstRef root) noexcept;
const flint::FmpzPoly* residue_polynomial_ptr(
        const PrimeIdeal& prime) noexcept;
const flint::Fmpz* linear_residue_root_ptr(
        const PrimeIdeal& prime) noexcept;
bool prime_ideal_valuation_with_norm_vp(
        slong& out,
        const PrimeIdeal& prime,
        const OrderElement& element,
        slong norm_vp,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

}  // namespace detail

class PrimeIdeal {
public:
    PrimeIdeal() noexcept;
    explicit PrimeIdeal(const Order& parent) noexcept;
    ~PrimeIdeal() noexcept;

    PrimeIdeal(const PrimeIdeal&) = delete;
    PrimeIdeal& operator=(const PrimeIdeal&) = delete;

    PrimeIdeal(PrimeIdeal&& other) noexcept;
    PrimeIdeal& operator=(PrimeIdeal&& other) noexcept;

    void swap(PrimeIdeal& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const PrimeIdeal& other) noexcept;

    bool is_defined() const noexcept;
    bool has_prime_data() const noexcept;
    const Order* parent() const noexcept;
    slong degree() const noexcept;

    bool rational_prime(flint::FmpzRef out) const noexcept;
    bool get_ideal(Ideal& out) const noexcept;
    bool kummer_generator_coordinates(flint::FmpzMatRef out) const noexcept;
    bool residue_polynomial(flint::FmpzPolyRef out) const noexcept;
    bool norm(flint::FmpzRef out) const noexcept;
    slong ramification_index() const noexcept;
    slong residue_degree() const noexcept;
    bool equal(const PrimeIdeal& other) const noexcept;

    bool reduce(flint::FmpzPolyRef out, const OrderElement& element) const noexcept;
    bool reduce(flint::FmpzPolyRef out, const Element& element) const noexcept;
    bool valuation(slong& out, const OrderElement& element) const noexcept;
    bool valuation(slong& out,
                   const OrderElement& element,
                   const DiagnosticsContext* diagnostics) const noexcept;
    bool valuation(slong& out, const Element& element) const noexcept;
    bool valuation(slong& out,
                   const Element& element,
                   const DiagnosticsContext* diagnostics) const noexcept;
    bool valuation(slong& out, const Ideal& ideal) const noexcept;
    bool valuation(slong& out,
                   const Ideal& ideal,
                   const DiagnosticsContext* diagnostics) const noexcept;
    bool valuation(slong& out, const FractionalIdeal& ideal) const noexcept;
    bool valuation(slong& out,
                   const FractionalIdeal& ideal,
                   const DiagnosticsContext* diagnostics) const noexcept;
    bool valuation(slong& out, const FactoredElement& element) const noexcept;
    bool valuation(slong& out,
                   const FactoredElement& element,
                   const DiagnosticsContext* diagnostics) const noexcept;

private:
    struct HenselValuationCache;
    struct ContainmentPowerCache;
    struct CoordinateValuationCache;

    bool set_data(flint::FmpzConstRef p,
                  slong ramification_index,
                  slong residue_degree,
                  const Ideal& ideal,
                  flint::FmpzPolyConstRef residue_polynomial) noexcept;
    bool set_data_no_residue(flint::FmpzConstRef p,
                             slong ramification_index,
                             slong residue_degree,
                             const Ideal& ideal) noexcept;
    bool set_kummer_generator(flint::FmpzMatConstRef generator) noexcept;
    const Ideal* containment_power_cached(
            slong exponent,
            const DiagnosticsContext* diagnostics) const
            noexcept;
    const flint::FmpzMat* coordinate_valuation_matrix_cached()
            const noexcept;
    bool valuation_by_power_containment(
            slong& out,
            const OrderElement& element,
            slong bound,
            const DiagnosticsContext* diagnostics) const noexcept;
    bool valuation_with_norm_vp_impl(
            slong& out,
            const OrderElement& element,
            slong norm_vp,
            const DiagnosticsContext* diagnostics) const noexcept;

    friend bool decompose_prime(PrimeIdealList& out,
                                const Order& order,
                                flint::FmpzConstRef p) noexcept;
    friend bool decompose_prime(PrimeIdealList& out,
                                const Order& order,
                                flint::FmpzConstRef p,
                                slong max_residue_degree) noexcept;
    friend bool decompose_prime(PrimeIdealList& out,
                                const Order& order,
                                flint::FmpzConstRef p,
                                slong max_residue_degree,
                                const DiagnosticsContext* diagnostics) noexcept;
    friend class detail::MaximalQuadraticPrimeAccess;
    friend bool detail::set_degree_one_prime_ideal_from_root(
            PrimeIdeal& out,
            const Order& order,
            flint::FmpzConstRef p,
            flint::FmpzConstRef root) noexcept;
    friend const flint::FmpzPoly* detail::residue_polynomial_ptr(
            const PrimeIdeal& prime) noexcept;
    friend const flint::Fmpz* detail::linear_residue_root_ptr(
            const PrimeIdeal& prime) noexcept;
    friend bool detail::prime_ideal_valuation_with_norm_vp(
            slong& out,
            const PrimeIdeal& prime,
            const OrderElement& element,
            slong norm_vp,
            const DiagnosticsContext* diagnostics) noexcept;

    Order parent_;
    flint::Fmpz p_;
    Ideal ideal_;
    flint::FmpzMat kummer_generator_{0, 0};
    flint::FmpzPoly residue_poly_;
    flint::Fmpz linear_residue_root_;
    mutable std::unique_ptr<HenselValuationCache> hensel_valuation_cache_;
    mutable std::unique_ptr<ContainmentPowerCache> containment_power_cache_;
    mutable std::unique_ptr<CoordinateValuationCache>
            coordinate_valuation_cache_;
    slong e_ = 0;
    slong f_ = 0;
    bool has_prime_ = false;
    bool has_kummer_generator_ = false;
    bool has_residue_poly_ = false;
    bool has_linear_residue_root_ = false;
};

class PrimeIdealList {
public:
    PrimeIdealList() noexcept = default;
    PrimeIdealList(const Order& parent, slong length) noexcept;
    ~PrimeIdealList() noexcept;

    PrimeIdealList(const PrimeIdealList&) = delete;
    PrimeIdealList& operator=(const PrimeIdealList&) = delete;

    PrimeIdealList(PrimeIdealList&& other) noexcept;
    PrimeIdealList& operator=(PrimeIdealList&& other) noexcept;

    void swap(PrimeIdealList& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent, slong length) noexcept;

    bool is_defined() const noexcept;
    slong size() const noexcept;
    PrimeIdeal* at(slong index) noexcept;
    const PrimeIdeal* at(slong index) const noexcept;

private:
    bool append(const Order& parent) noexcept;
    PrimeIdeal* back() noexcept;
    bool reserve(slong capacity) noexcept;

    friend bool decompose_prime(PrimeIdealList& out,
                                const Order& order,
                                flint::FmpzConstRef p) noexcept;
    friend bool decompose_prime(PrimeIdealList& out,
                                const Order& order,
                                flint::FmpzConstRef p,
                                slong max_residue_degree) noexcept;
    friend bool decompose_prime(PrimeIdealList& out,
                                const Order& order,
                                flint::FmpzConstRef p,
                                slong max_residue_degree,
                                const DiagnosticsContext* diagnostics) noexcept;

    std::unique_ptr<PrimeIdeal[]> entries_;
    slong size_ = 0;
    slong capacity_ = 0;
    bool defined_ = false;
};

bool decompose_prime(PrimeIdealList& out,
                     const Order& order,
                     flint::FmpzConstRef p) noexcept;
bool decompose_prime(PrimeIdealList& out,
                     const Order& order,
                     flint::FmpzConstRef p,
                     slong max_residue_degree) noexcept;
bool decompose_prime(PrimeIdealList& out,
                     const Order& order,
                     flint::FmpzConstRef p,
                     slong max_residue_degree,
                     const DiagnosticsContext* diagnostics) noexcept;

inline void swap(PrimeIdeal& left, PrimeIdeal& right) noexcept {
    left.swap(right);
}

inline void swap(PrimeIdealList& left, PrimeIdealList& right) noexcept {
    left.swap(right);
}

}  // namespace silex
