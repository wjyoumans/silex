#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>

#include <flint/fmpz_mod_mat.h>

namespace silex::flint {

class FmpzModMatRef;
class FmpzModMatConstRef;

class FmpzModMat {
public:
    FmpzModMat(slong rows, slong cols, const fmpz_mod_ctx_t& ctx) noexcept
        : FmpzModMat(rows, cols,
                     static_cast<const fmpz_mod_ctx_struct*>(ctx)) {
    }

    FmpzModMat(slong rows,
               slong cols,
               const fmpz_mod_ctx_struct* ctx) noexcept
        : ctx_(ctx) {
        fmpz_mod_mat_init(value_, rows, cols, ctx_);
    }

    FmpzModMat(slong rows, slong cols, const FmpzModCtx& ctx) noexcept
        : FmpzModMat(rows, cols,
                     static_cast<const fmpz_mod_ctx_struct*>(ctx.raw())) {
    }

    ~FmpzModMat() noexcept { fmpz_mod_mat_clear(value_, ctx_); }

    FmpzModMat(const FmpzModMat&) = delete;
    FmpzModMat& operator=(const FmpzModMat&) = delete;
    FmpzModMat(FmpzModMat&&) = delete;
    FmpzModMat& operator=(FmpzModMat&&) = delete;

    fmpz_mod_mat_t& raw() noexcept { return value_; }
    const fmpz_mod_mat_t& raw() const noexcept { return value_; }
    const fmpz_mod_ctx_struct* context() const noexcept { return ctx_; }

private:
    const fmpz_mod_ctx_struct* ctx_;
    fmpz_mod_mat_t value_;
};

class FmpzModMatConstRef {
public:
    explicit FmpzModMatConstRef(const fmpz_mod_mat_t& value) noexcept : value_(value) {}
    explicit FmpzModMatConstRef(const FmpzModMat& value) noexcept : value_(value.raw()) {}
    const fmpz_mod_mat_struct* raw() const noexcept { return value_; }

private:
    const fmpz_mod_mat_struct* value_;
};

class FmpzModMatRef {
public:
    explicit FmpzModMatRef(fmpz_mod_mat_t& value) noexcept : value_(value) {}
    explicit FmpzModMatRef(FmpzModMat& value) noexcept : value_(value.raw()) {}
    fmpz_mod_mat_struct* raw() noexcept { return value_; }
    const fmpz_mod_mat_struct* raw() const noexcept { return value_; }

private:
    fmpz_mod_mat_struct* value_;
};

inline void fmpz_mod_mat_set_fmpz_mat(FmpzModMatRef out,
                                      FmpzMatConstRef in,
                                      FmpzModCtxConstRef ctx) noexcept {
    ::fmpz_mod_mat_set_fmpz_mat(out.raw(), in.raw(), ctx.raw());
}

inline slong fmpz_mod_mat_rank(FmpzModMatRef matrix,
                               FmpzModCtxConstRef ctx) noexcept {
    return ::fmpz_mod_mat_rank(matrix.raw(), ctx.raw());
}

inline void fmpz_mod_mat_get_entry(FmpzRef out,
                                   FmpzModMatConstRef matrix,
                                   slong row,
                                   slong col,
                                   FmpzModCtxConstRef ctx) noexcept {
    ::fmpz_mod_mat_get_entry(out.raw(), matrix.raw(), row, col, ctx.raw());
}

}  // namespace silex::flint
