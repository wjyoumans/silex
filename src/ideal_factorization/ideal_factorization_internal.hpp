#pragma once

#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex {

class FactorBase;
class Ideal;
class OrderElement;
class PrimeIdeal;
struct DiagnosticsContext;

namespace detail {

enum class OneLargePrimeFactorStatus {
    no_candidate,
    found,
};

bool order_element_factor_over_base_with_one_large_prime(
        OneLargePrimeFactorStatus& status,
        flint::FmpzMatRef exponents,
        PrimeIdeal& large_prime,
        const OrderElement& element,
        flint::FmpqConstRef norm,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool ideal_factor_over_base_with_required_prime(
        bool& matches,
        const Ideal& ideal,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool order_element_factor_over_base_with_required_prime(
        bool& matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

}  // namespace detail
}  // namespace silex
