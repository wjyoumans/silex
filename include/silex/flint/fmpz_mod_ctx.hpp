#pragma once

#include <flint/fmpz_mod.h>

#include <utility>

namespace silex::flint {

class FmpzModCtxRef;
class FmpzModCtxConstRef;

class FmpzModCtx {
public:
    FmpzModCtx() noexcept
        : FmpzModCtx(2UL) {
    }

    explicit FmpzModCtx(const fmpz_t modulus) noexcept {
        value_ = static_cast<fmpz_mod_ctx_struct*>(
                flint_malloc(sizeof(fmpz_mod_ctx_struct)));
        if (value_ != nullptr) {
            fmpz_mod_ctx_init(value_, modulus);
        }
    }

    explicit FmpzModCtx(ulong modulus) noexcept {
        value_ = static_cast<fmpz_mod_ctx_struct*>(
                flint_malloc(sizeof(fmpz_mod_ctx_struct)));
        if (value_ != nullptr) {
            fmpz_mod_ctx_init_ui(value_, modulus);
        }
    }

    ~FmpzModCtx() noexcept {
        if (value_ != nullptr) {
            fmpz_mod_ctx_clear(value_);
            flint_free(value_);
        }
    }

    FmpzModCtx(const FmpzModCtx&) = delete;
    FmpzModCtx& operator=(const FmpzModCtx&) = delete;

    FmpzModCtx(FmpzModCtx&& other) noexcept
        : value_(other.value_) {
        other.value_ = nullptr;
    }

    FmpzModCtx& operator=(FmpzModCtx&& other) noexcept {
        if (this != &other) {
            FmpzModCtx tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(FmpzModCtx& other) noexcept {
        std::swap(value_, other.value_);
    }

    fmpz_mod_ctx_struct* raw() noexcept { return value_; }
    const fmpz_mod_ctx_struct* raw() const noexcept { return value_; }

private:
    fmpz_mod_ctx_struct* value_ = nullptr;
};

class FmpzModCtxConstRef {
public:
    explicit FmpzModCtxConstRef(const fmpz_mod_ctx_t& value) noexcept : value_(value) {}
    explicit FmpzModCtxConstRef(const fmpz_mod_ctx_struct* value) noexcept
        : value_(value) {
    }
    explicit FmpzModCtxConstRef(const FmpzModCtx& value) noexcept : value_(value.raw()) {}
    const fmpz_mod_ctx_struct* raw() const noexcept { return value_; }

private:
    const fmpz_mod_ctx_struct* value_;
};

class FmpzModCtxRef {
public:
    explicit FmpzModCtxRef(fmpz_mod_ctx_t& value) noexcept : value_(value) {}
    explicit FmpzModCtxRef(fmpz_mod_ctx_struct* value) noexcept : value_(value) {}
    explicit FmpzModCtxRef(FmpzModCtx& value) noexcept : value_(value.raw()) {}
    fmpz_mod_ctx_struct* raw() noexcept { return value_; }
    const fmpz_mod_ctx_struct* raw() const noexcept { return value_; }

private:
    fmpz_mod_ctx_struct* value_;
};

}  // namespace silex::flint
