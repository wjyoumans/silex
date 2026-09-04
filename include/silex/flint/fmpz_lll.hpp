#pragma once

#include <flint/fmpz_lll.h>

namespace silex::flint {

class FmpzLllRef;
class FmpzLllConstRef;

class FmpzLll {
public:
    FmpzLll() noexcept { fmpz_lll_context_init_default(value_); }

    FmpzLll(double delta, double eta, rep_type representation, gram_type gram) noexcept {
        fmpz_lll_context_init(value_, delta, eta, representation, gram);
    }

    FmpzLll(const FmpzLll&) = delete;
    FmpzLll& operator=(const FmpzLll&) = delete;

    FmpzLll(FmpzLll&& other) noexcept {
        fmpz_lll_context_init_default(value_);
        swap(other);
    }

    FmpzLll& operator=(FmpzLll&& other) noexcept {
        if (this != &other) {
            swap(other);
        }
        return *this;
    }

    void swap(FmpzLll& other) noexcept {
        const fmpz_lll_struct tmp = value_[0];
        value_[0] = other.value_[0];
        other.value_[0] = tmp;
    }

    fmpz_lll_t& raw() noexcept { return value_; }
    const fmpz_lll_t& raw() const noexcept { return value_; }

private:
    fmpz_lll_t value_;
};

class FmpzLllConstRef {
public:
    explicit FmpzLllConstRef(const fmpz_lll_t& value) noexcept : value_(value) {}
    explicit FmpzLllConstRef(const FmpzLll& value) noexcept : value_(value.raw()) {}
    const fmpz_lll_struct* raw() const noexcept { return value_; }

private:
    const fmpz_lll_struct* value_;
};

class FmpzLllRef {
public:
    explicit FmpzLllRef(fmpz_lll_t& value) noexcept : value_(value) {}
    explicit FmpzLllRef(FmpzLll& value) noexcept : value_(value.raw()) {}
    fmpz_lll_struct* raw() noexcept { return value_; }
    const fmpz_lll_struct* raw() const noexcept { return value_; }

private:
    fmpz_lll_struct* value_;
};

inline void swap(FmpzLll& left, FmpzLll& right) noexcept { left.swap(right); }

}  // namespace silex::flint
