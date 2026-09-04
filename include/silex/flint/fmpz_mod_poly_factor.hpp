#pragma once

#include <silex/flint/fmpz_mod_ctx.hpp>

#include <flint/fmpz_mod_poly_factor.h>

namespace silex::flint {

class FmpzModPolyFactorRef;
class FmpzModPolyFactorConstRef;

class FmpzModPolyFactor {
public:
    explicit FmpzModPolyFactor(const fmpz_mod_ctx_t& ctx) noexcept
        : FmpzModPolyFactor(static_cast<const fmpz_mod_ctx_struct*>(ctx)) {
    }

    explicit FmpzModPolyFactor(const FmpzModCtx& ctx) noexcept
        : FmpzModPolyFactor(static_cast<const fmpz_mod_ctx_struct*>(ctx.raw())) {
    }

    explicit FmpzModPolyFactor(const fmpz_mod_ctx_struct* ctx) noexcept
        : ctx_(ctx) {
        fmpz_mod_poly_factor_init(value_, ctx_);
    }

    ~FmpzModPolyFactor() noexcept { fmpz_mod_poly_factor_clear(value_, ctx_); }

    FmpzModPolyFactor(const FmpzModPolyFactor&) = delete;
    FmpzModPolyFactor& operator=(const FmpzModPolyFactor&) = delete;

    FmpzModPolyFactor(FmpzModPolyFactor&& other) noexcept
        : FmpzModPolyFactor(other.ctx_) {
        fmpz_mod_poly_factor_swap(value_, other.value_, ctx_);
    }

    FmpzModPolyFactor& operator=(FmpzModPolyFactor&& other) noexcept {
        if (this != &other && ctx_ == other.ctx_) {
            fmpz_mod_poly_factor_swap(value_, other.value_, ctx_);
        }
        return *this;
    }

    void swap(FmpzModPolyFactor& other) noexcept {
        if (ctx_ == other.ctx_) {
            fmpz_mod_poly_factor_swap(value_, other.value_, ctx_);
        }
    }

    fmpz_mod_poly_factor_t& raw() noexcept { return value_; }
    const fmpz_mod_poly_factor_t& raw() const noexcept { return value_; }
    const fmpz_mod_ctx_struct* context() const noexcept { return ctx_; }

private:
    const fmpz_mod_ctx_struct* ctx_;
    fmpz_mod_poly_factor_t value_;
};

class FmpzModPolyFactorConstRef {
public:
    explicit FmpzModPolyFactorConstRef(const fmpz_mod_poly_factor_t& value) noexcept
        : value_(value) {
    }

    explicit FmpzModPolyFactorConstRef(const FmpzModPolyFactor& value) noexcept
        : value_(value.raw()) {
    }

    const fmpz_mod_poly_factor_struct* raw() const noexcept { return value_; }

private:
    const fmpz_mod_poly_factor_struct* value_;
};

class FmpzModPolyFactorRef {
public:
    explicit FmpzModPolyFactorRef(fmpz_mod_poly_factor_t& value) noexcept : value_(value) {}
    explicit FmpzModPolyFactorRef(FmpzModPolyFactor& value) noexcept : value_(value.raw()) {}
    fmpz_mod_poly_factor_struct* raw() noexcept { return value_; }
    const fmpz_mod_poly_factor_struct* raw() const noexcept { return value_; }

private:
    fmpz_mod_poly_factor_struct* value_;
};

inline void swap(FmpzModPolyFactor& left, FmpzModPolyFactor& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
