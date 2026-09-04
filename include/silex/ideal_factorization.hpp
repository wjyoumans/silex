#pragma once

#include <flint/flint.h>

#include <silex/factor_base.hpp>
#include <silex/factored.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>

#include <vector>

namespace silex {

struct DiagnosticsContext;

class IdealFactorization {
public:
    IdealFactorization() noexcept = default;
    explicit IdealFactorization(const Order& parent) noexcept;
    ~IdealFactorization() noexcept;

    IdealFactorization(const IdealFactorization&) = delete;
    IdealFactorization& operator=(const IdealFactorization&) = delete;

    IdealFactorization(IdealFactorization&& other) noexcept;
    IdealFactorization& operator=(IdealFactorization&& other) noexcept;

    void swap(IdealFactorization& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const IdealFactorization& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong length() const noexcept;

    bool factor(const Ideal& ideal) noexcept;
    bool prime(PrimeIdeal& out, slong index) const noexcept;
    bool exponent(slong& out, slong index) const noexcept;
    bool reconstruct(Ideal& out) const noexcept;

private:
    Order parent_;
    std::vector<FactorPower<PrimeIdeal>> factors_;
};

bool ideal_is_smooth(bool& smooth,
                     const Ideal& ideal,
                     const FactorBase& base) noexcept;

bool ideal_factor_over_base(flint::FmpzMatRef exponents,
                            const Ideal& ideal,
                            const FactorBase& base,
                            const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool ideal_factor_over_base(flint::FmpzMatRef exponents,
                            const FractionalIdeal& ideal,
                            const FactorBase& base,
                            const DiagnosticsContext* diagnostics = nullptr) noexcept;

inline void swap(IdealFactorization& left,
                 IdealFactorization& right) noexcept {
    left.swap(right);
}

}  // namespace silex
