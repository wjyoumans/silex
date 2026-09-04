#pragma once

#include <flint/fmpq_mat.h>

#include <silex/flint/fmpq.hpp>

namespace silex::flint {

class FmpqMatRef;
class FmpqMatConstRef;

class FmpqMat {
public:
    FmpqMat(slong rows, slong cols) noexcept {
        fmpq_mat_init(value_, rows, cols);
    }

    ~FmpqMat() noexcept {
        fmpq_mat_clear(value_);
    }

    FmpqMat(const FmpqMat&) = delete;
    FmpqMat& operator=(const FmpqMat&) = delete;

    FmpqMat(FmpqMat&& other) noexcept
        : FmpqMat(0, 0) {
        fmpq_mat_swap(value_, other.value_);
    }

    FmpqMat& operator=(FmpqMat&& other) noexcept {
        if (this != &other) {
            fmpq_mat_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(FmpqMat& other) noexcept {
        fmpq_mat_swap(value_, other.value_);
    }

    fmpq_mat_t& raw() noexcept {
        return value_;
    }

    const fmpq_mat_t& raw() const noexcept {
        return value_;
    }

private:
    fmpq_mat_t value_;
};

class FmpqMatConstRef {
public:
    explicit FmpqMatConstRef(const fmpq_mat_t& value) noexcept
        : value_(value) {
    }

    explicit FmpqMatConstRef(const FmpqMat& value) noexcept
        : value_(value.raw()) {
    }

    const fmpq_mat_struct* raw() const noexcept {
        return value_;
    }

private:
    const fmpq_mat_struct* value_;
};

class FmpqMatRef {
public:
    explicit FmpqMatRef(fmpq_mat_t& value) noexcept
        : value_(value) {
    }

    explicit FmpqMatRef(FmpqMat& value) noexcept
        : value_(value.raw()) {
    }

    fmpq_mat_struct* raw() noexcept {
        return value_;
    }

    const fmpq_mat_struct* raw() const noexcept {
        return value_;
    }

private:
    fmpq_mat_struct* value_;
};

inline void swap(FmpqMat& left, FmpqMat& right) noexcept {
    left.swap(right);
}

inline void fmpq_mat_zero(FmpqMatRef matrix) noexcept {
    ::fmpq_mat_zero(matrix.raw());
}

inline void fmpq_mat_zero(FmpqMat& matrix) noexcept {
    ::fmpq_mat_zero(matrix.raw());
}

inline FmpqConstRef fmpq_mat_entry(
        FmpqMatConstRef matrix, slong row, slong col) noexcept {
    return FmpqConstRef(::fmpq_mat_entry(matrix.raw(), row, col));
}

inline FmpqConstRef fmpq_mat_entry(
        const FmpqMat& matrix, slong row, slong col) noexcept {
    return FmpqConstRef(::fmpq_mat_entry(matrix.raw(), row, col));
}

inline FmpqRef fmpq_mat_entry(FmpqMatRef matrix, slong row, slong col) noexcept {
    return FmpqRef(::fmpq_mat_entry(matrix.raw(), row, col));
}

inline FmpqRef fmpq_mat_entry(FmpqMat& matrix, slong row, slong col) noexcept {
    return FmpqRef(::fmpq_mat_entry(matrix.raw(), row, col));
}

inline bool fmpq_mat_equal(FmpqMatConstRef left, FmpqMatConstRef right) noexcept {
    return ::fmpq_mat_equal(left.raw(), right.raw()) != 0;
}

inline bool fmpq_mat_equal(const FmpqMat& left, const FmpqMat& right) noexcept {
    return ::fmpq_mat_equal(left.raw(), right.raw()) != 0;
}

}  // namespace silex::flint
