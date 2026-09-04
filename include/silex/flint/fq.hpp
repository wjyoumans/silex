#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>

#include <flint/fq.h>

#include <utility>

namespace silex::flint {

class FqCtx {
public:
    FqCtx() noexcept = default;

    FqCtx(const fmpz_mod_poly_t modulus,
          const fmpz_mod_ctx_t ctx,
          const char* variable) noexcept
        : value_(static_cast<fq_ctx_struct*>(
                  flint_malloc(sizeof(fq_ctx_struct)))) {
        if (value_ != nullptr) {
            fq_ctx_init_modulus(value_, modulus, ctx, variable);
        }
    }

    FqCtx(const FmpzModPoly& modulus,
          const FmpzModCtx& ctx,
          const char* variable) noexcept
        : FqCtx(modulus.raw(), ctx.raw(), variable) {
    }

    ~FqCtx() noexcept {
        if (value_ != nullptr) {
            fq_ctx_clear(value_);
            flint_free(value_);
        }
    }

    FqCtx(const FqCtx&) = delete;
    FqCtx& operator=(const FqCtx&) = delete;

    FqCtx(FqCtx&& other) noexcept
        : value_(other.value_) {
        other.value_ = nullptr;
    }

    FqCtx& operator=(FqCtx&& other) noexcept {
        if (this != &other) {
            FqCtx tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(FqCtx& other) noexcept {
        std::swap(value_, other.value_);
    }

    bool is_initialized() const noexcept { return value_ != nullptr; }
    fq_ctx_struct* raw() noexcept { return value_; }
    const fq_ctx_struct* raw() const noexcept { return value_; }

private:
    fq_ctx_struct* value_ = nullptr;
};

class Fq {
public:
    Fq() noexcept = default;

    explicit Fq(const FqCtx& ctx) noexcept
        : ctx_(ctx.raw()),
          value_(static_cast<fq_struct*>(flint_malloc(sizeof(fq_struct)))) {
        if (value_ != nullptr && ctx_ != nullptr) {
            fq_init(value_, ctx_);
        } else if (value_ != nullptr) {
            flint_free(value_);
            value_ = nullptr;
        }
    }

    ~Fq() noexcept {
        if (value_ != nullptr) {
            fq_clear(value_, ctx_);
            flint_free(value_);
        }
    }

    Fq(const Fq&) = delete;
    Fq& operator=(const Fq&) = delete;

    Fq(Fq&& other) noexcept
        : ctx_(other.ctx_),
          value_(other.value_) {
        other.ctx_ = nullptr;
        other.value_ = nullptr;
    }

    Fq& operator=(Fq&& other) noexcept {
        if (this != &other) {
            Fq tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(Fq& other) noexcept {
        std::swap(ctx_, other.ctx_);
        std::swap(value_, other.value_);
    }

    bool is_initialized() const noexcept { return value_ != nullptr; }
    fq_struct* raw() noexcept { return value_; }
    const fq_struct* raw() const noexcept { return value_; }

private:
    const fq_ctx_struct* ctx_ = nullptr;
    fq_struct* value_ = nullptr;
};

inline void swap(FqCtx& left, FqCtx& right) noexcept {
    left.swap(right);
}

inline void swap(Fq& left, Fq& right) noexcept {
    left.swap(right);
}

inline void fq_ctx_order(FmpzRef out, const FqCtx& ctx) noexcept {
    ::fq_ctx_order(out.raw(), ctx.raw());
}

inline void fq_pow(Fq& out,
                   const Fq& base,
                   FmpzConstRef exponent,
                   const FqCtx& ctx) noexcept {
    ::fq_pow(out.raw(), base.raw(), exponent.raw(), ctx.raw());
}

inline void fq_inv(Fq& out, const Fq& in, const FqCtx& ctx) noexcept {
    ::fq_inv(out.raw(), in.raw(), ctx.raw());
}

inline bool fq_equal(const Fq& left,
                     const Fq& right,
                     const FqCtx& ctx) noexcept {
    return ::fq_equal(left.raw(), right.raw(), ctx.raw()) != 0;
}

inline bool fq_is_zero(const Fq& value, const FqCtx& ctx) noexcept {
    return ::fq_is_zero(value.raw(), ctx.raw()) != 0;
}

inline bool fq_is_one(const Fq& value, const FqCtx& ctx) noexcept {
    return ::fq_is_one(value.raw(), ctx.raw()) != 0;
}

inline void fq_get_fmpz_mod_poly(FmpzModPoly& out,
                                 const Fq& value,
                                 const FqCtx& ctx) noexcept {
    ::fq_get_fmpz_mod_poly(out.raw(), value.raw(), ctx.raw());
}

}  // namespace silex::flint
