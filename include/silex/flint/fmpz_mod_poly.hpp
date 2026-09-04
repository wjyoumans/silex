#pragma once

#include <silex/flint/fmpz_mod_ctx.hpp>

#include <flint/fmpz_mod_poly.h>

#include <utility>

namespace silex::flint {

class FmpzModPolyRef;
class FmpzModPolyConstRef;

class FmpzModPoly {
public:
    FmpzModPoly() noexcept = default;

    explicit FmpzModPoly(const fmpz_mod_ctx_t& ctx) noexcept
        : FmpzModPoly(static_cast<const fmpz_mod_ctx_struct*>(ctx)) {
    }

    explicit FmpzModPoly(const fmpz_mod_ctx_struct* ctx) noexcept
        : ctx_(ctx),
          value_(static_cast<fmpz_mod_poly_struct*>(
                  flint_malloc(sizeof(fmpz_mod_poly_struct)))) {
        if (value_ != nullptr) {
            fmpz_mod_poly_init(value_, ctx_);
        }
    }

    explicit FmpzModPoly(const FmpzModCtx& ctx) noexcept
        : FmpzModPoly(static_cast<const fmpz_mod_ctx_struct*>(ctx.raw())) {
    }

    ~FmpzModPoly() noexcept {
        if (value_ != nullptr) {
            fmpz_mod_poly_clear(value_, ctx_);
            flint_free(value_);
        }
    }

    FmpzModPoly(const FmpzModPoly&) = delete;
    FmpzModPoly& operator=(const FmpzModPoly&) = delete;

    FmpzModPoly(FmpzModPoly&& other) noexcept
        : ctx_(other.ctx_),
          value_(other.value_) {
        other.ctx_ = nullptr;
        other.value_ = nullptr;
    }

    FmpzModPoly& operator=(FmpzModPoly&& other) noexcept {
        if (this != &other) {
            FmpzModPoly tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(FmpzModPoly& other) noexcept {
        std::swap(ctx_, other.ctx_);
        std::swap(value_, other.value_);
    }

    bool is_initialized() const noexcept { return value_ != nullptr; }
    fmpz_mod_poly_struct* raw() noexcept { return value_; }
    const fmpz_mod_poly_struct* raw() const noexcept { return value_; }
    const fmpz_mod_ctx_struct* context() const noexcept { return ctx_; }

private:
    const fmpz_mod_ctx_struct* ctx_ = nullptr;
    fmpz_mod_poly_struct* value_ = nullptr;
};

class FmpzModPolyConstRef {
public:
    explicit FmpzModPolyConstRef(const fmpz_mod_poly_t& value) noexcept : value_(value) {}
    explicit FmpzModPolyConstRef(const fmpz_mod_poly_struct* value) noexcept
        : value_(value) {
    }
    explicit FmpzModPolyConstRef(const FmpzModPoly& value) noexcept : value_(value.raw()) {}
    const fmpz_mod_poly_struct* raw() const noexcept { return value_; }

private:
    const fmpz_mod_poly_struct* value_;
};

class FmpzModPolyRef {
public:
    explicit FmpzModPolyRef(fmpz_mod_poly_t& value) noexcept : value_(value) {}
    explicit FmpzModPolyRef(fmpz_mod_poly_struct* value) noexcept : value_(value) {}
    explicit FmpzModPolyRef(FmpzModPoly& value) noexcept : value_(value.raw()) {}
    fmpz_mod_poly_struct* raw() noexcept { return value_; }
    const fmpz_mod_poly_struct* raw() const noexcept { return value_; }

private:
    fmpz_mod_poly_struct* value_;
};

}  // namespace silex::flint
