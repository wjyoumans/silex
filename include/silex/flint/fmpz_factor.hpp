#pragma once

#include <silex/flint/fmpz.hpp>

#include <flint/fmpz_factor.h>

namespace silex::flint {

class FmpzFactorRef;
class FmpzFactorConstRef;

class FmpzFactor {
public:
    FmpzFactor() noexcept { fmpz_factor_init(value_); }
    ~FmpzFactor() noexcept { fmpz_factor_clear(value_); }

    FmpzFactor(const FmpzFactor&) = delete;
    FmpzFactor& operator=(const FmpzFactor&) = delete;
    FmpzFactor(FmpzFactor&&) = delete;
    FmpzFactor& operator=(FmpzFactor&&) = delete;

    fmpz_factor_t& raw() noexcept { return value_; }
    const fmpz_factor_t& raw() const noexcept { return value_; }

private:
    fmpz_factor_t value_;
};

class FmpzFactorConstRef {
public:
    explicit FmpzFactorConstRef(const fmpz_factor_t& value) noexcept : value_(value) {}
    explicit FmpzFactorConstRef(const FmpzFactor& value) noexcept : value_(value.raw()) {}
    const fmpz_factor_struct* raw() const noexcept { return value_; }

private:
    const fmpz_factor_struct* value_;
};

class FmpzFactorRef {
public:
    explicit FmpzFactorRef(fmpz_factor_t& value) noexcept : value_(value) {}
    explicit FmpzFactorRef(FmpzFactor& value) noexcept : value_(value.raw()) {}
    fmpz_factor_struct* raw() noexcept { return value_; }
    const fmpz_factor_struct* raw() const noexcept { return value_; }

private:
    fmpz_factor_struct* value_;
};

inline void fmpz_factor(FmpzFactorRef out, FmpzConstRef value) noexcept {
    ::fmpz_factor(out.raw(), value.raw());
}

inline slong fmpz_factor_num(FmpzFactorConstRef factorization) noexcept {
    return factorization.raw()->num;
}

inline void fmpz_factor_get_fmpz(FmpzRef out,
                                 FmpzFactorConstRef factorization,
                                 slong index) noexcept {
    ::fmpz_factor_get_fmpz(out.raw(), factorization.raw(), index);
}

inline ulong fmpz_factor_exp(FmpzFactorConstRef factorization,
                             slong index) noexcept {
    return factorization.raw()->exp[index];
}

}  // namespace silex::flint
