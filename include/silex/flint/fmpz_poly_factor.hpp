#pragma once

#include <silex/flint/fmpz_poly.hpp>

#include <flint/fmpz_poly_factor.h>

namespace silex::flint {

class FmpzPolyFactorRef;
class FmpzPolyFactorConstRef;

class FmpzPolyFactor {
public:
    FmpzPolyFactor() noexcept { fmpz_poly_factor_init(value_); }
    ~FmpzPolyFactor() noexcept { fmpz_poly_factor_clear(value_); }

    FmpzPolyFactor(const FmpzPolyFactor&) = delete;
    FmpzPolyFactor& operator=(const FmpzPolyFactor&) = delete;
    FmpzPolyFactor(FmpzPolyFactor&&) = delete;
    FmpzPolyFactor& operator=(FmpzPolyFactor&&) = delete;

    fmpz_poly_factor_t& raw() noexcept { return value_; }
    const fmpz_poly_factor_t& raw() const noexcept { return value_; }

private:
    fmpz_poly_factor_t value_;
};

class FmpzPolyFactorConstRef {
public:
    explicit FmpzPolyFactorConstRef(
            const fmpz_poly_factor_t& value) noexcept
        : value_(value) {
    }

    explicit FmpzPolyFactorConstRef(const FmpzPolyFactor& value) noexcept
        : value_(value.raw()) {
    }

    const fmpz_poly_factor_struct* raw() const noexcept { return value_; }

private:
    const fmpz_poly_factor_struct* value_;
};

class FmpzPolyFactorRef {
public:
    explicit FmpzPolyFactorRef(fmpz_poly_factor_t& value) noexcept
        : value_(value) {
    }

    explicit FmpzPolyFactorRef(FmpzPolyFactor& value) noexcept
        : value_(value.raw()) {
    }

    fmpz_poly_factor_struct* raw() noexcept { return value_; }
    const fmpz_poly_factor_struct* raw() const noexcept { return value_; }

private:
    fmpz_poly_factor_struct* value_;
};

inline void fmpz_poly_factor(FmpzPolyFactorRef out,
                             FmpzPolyConstRef polynomial) noexcept {
    ::fmpz_poly_factor(out.raw(), polynomial.raw());
}

inline slong fmpz_poly_factor_num(
        FmpzPolyFactorConstRef factorization) noexcept {
    return factorization.raw()->num;
}

inline slong fmpz_poly_factor_exp(FmpzPolyFactorConstRef factorization,
                                  slong index) noexcept {
    return factorization.raw()->exp[index];
}

inline slong fmpz_poly_factor_poly_degree(
        FmpzPolyFactorConstRef factorization,
        slong index) noexcept {
    return ::fmpz_poly_degree(factorization.raw()->p + index);
}

}  // namespace silex::flint
