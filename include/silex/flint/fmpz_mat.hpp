#pragma once

#include <silex/flint/fmpz.hpp>

#include <flint/fmpz_mat.h>

namespace silex::flint {

class FmpzMatRef;
class FmpzMatConstRef;

class FmpzMat {
public:
    FmpzMat(slong rows, slong cols) noexcept {
        fmpz_mat_init(value_, rows, cols);
    }

    ~FmpzMat() noexcept {
        fmpz_mat_clear(value_);
    }

    FmpzMat(const FmpzMat&) = delete;
    FmpzMat& operator=(const FmpzMat&) = delete;

    FmpzMat(FmpzMat&& other) noexcept
        : FmpzMat(0, 0) {
        fmpz_mat_swap(value_, other.value_);
    }

    FmpzMat& operator=(FmpzMat&& other) noexcept {
        if (this != &other) {
            fmpz_mat_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(FmpzMat& other) noexcept {
        fmpz_mat_swap(value_, other.value_);
    }

    fmpz_mat_t& raw() noexcept {
        return value_;
    }

    const fmpz_mat_t& raw() const noexcept {
        return value_;
    }

private:
    fmpz_mat_t value_;
};

class FmpzMatConstRef {
public:
    FmpzMatConstRef(const fmpz_mat_struct* value) noexcept
        : value_(value) {
    }

    FmpzMatConstRef(const FmpzMat& value) noexcept
        : value_(value.raw()) {
    }

    const fmpz_mat_struct* raw() const noexcept {
        return value_;
    }

private:
    const fmpz_mat_struct* value_;
};

class FmpzMatRef {
public:
    FmpzMatRef(fmpz_mat_struct* value) noexcept
        : value_(value) {
    }

    FmpzMatRef(FmpzMat& value) noexcept
        : value_(value.raw()) {
    }

    fmpz_mat_struct* raw() noexcept {
        return value_;
    }

    const fmpz_mat_struct* raw() const noexcept {
        return value_;
    }

private:
    fmpz_mat_struct* value_;
};

class FmpzMatWindow {
public:
    FmpzMatWindow(FmpzMatRef matrix, slong r1, slong c1, slong r2, slong c2) noexcept {
        fmpz_mat_window_init(value_, matrix.raw(), r1, c1, r2, c2);
    }

    ~FmpzMatWindow() noexcept {
        fmpz_mat_window_clear(value_);
    }

    FmpzMatWindow(const FmpzMatWindow&) = delete;
    FmpzMatWindow& operator=(const FmpzMatWindow&) = delete;
    FmpzMatWindow(FmpzMatWindow&&) = delete;
    FmpzMatWindow& operator=(FmpzMatWindow&&) = delete;

    fmpz_mat_t& raw() noexcept {
        return value_;
    }

    const fmpz_mat_t& raw() const noexcept {
        return value_;
    }

    FmpzMatRef ref() noexcept {
        return FmpzMatRef(value_);
    }

    FmpzMatConstRef const_ref() const noexcept {
        return FmpzMatConstRef(value_);
    }

private:
    fmpz_mat_t value_;
};

class FmpzMatConstWindow {
public:
    FmpzMatConstWindow(
            FmpzMatConstRef matrix, slong r1, slong c1, slong r2, slong c2) noexcept {
        fmpz_mat_window_init(value_, matrix.raw(), r1, c1, r2, c2);
    }

    ~FmpzMatConstWindow() noexcept {
        fmpz_mat_window_clear(value_);
    }

    FmpzMatConstWindow(const FmpzMatConstWindow&) = delete;
    FmpzMatConstWindow& operator=(const FmpzMatConstWindow&) = delete;
    FmpzMatConstWindow(FmpzMatConstWindow&&) = delete;
    FmpzMatConstWindow& operator=(FmpzMatConstWindow&&) = delete;

    const fmpz_mat_t& raw() const noexcept {
        return value_;
    }

    FmpzMatConstRef const_ref() const noexcept {
        return FmpzMatConstRef(value_);
    }

private:
    fmpz_mat_t value_;
};

inline void swap(FmpzMat& left, FmpzMat& right) noexcept {
    left.swap(right);
}

inline slong fmpz_mat_nrows(FmpzMatConstRef matrix) noexcept {
    return ::fmpz_mat_nrows(matrix.raw());
}

inline slong fmpz_mat_nrows(const FmpzMat& matrix) noexcept {
    return ::fmpz_mat_nrows(matrix.raw());
}

inline slong fmpz_mat_nrows(FmpzMatRef matrix) noexcept {
    return ::fmpz_mat_nrows(matrix.raw());
}

inline slong fmpz_mat_ncols(FmpzMatConstRef matrix) noexcept {
    return ::fmpz_mat_ncols(matrix.raw());
}

inline slong fmpz_mat_ncols(const FmpzMat& matrix) noexcept {
    return ::fmpz_mat_ncols(matrix.raw());
}

inline slong fmpz_mat_ncols(FmpzMatRef matrix) noexcept {
    return ::fmpz_mat_ncols(matrix.raw());
}

inline void fmpz_mat_zero(FmpzMatRef matrix) noexcept {
    ::fmpz_mat_zero(matrix.raw());
}

inline void fmpz_mat_one(FmpzMatRef matrix) noexcept {
    ::fmpz_mat_one(matrix.raw());
}

inline slong fmpz_mat_rank(FmpzMatConstRef matrix) noexcept {
    return ::fmpz_mat_rank(matrix.raw());
}

inline FmpzConstRef fmpz_mat_entry(
        FmpzMatConstRef matrix, slong row, slong col) noexcept {
    return FmpzConstRef(::fmpz_mat_entry(matrix.raw(), row, col));
}

inline FmpzConstRef fmpz_mat_entry(
        const FmpzMat& matrix, slong row, slong col) noexcept {
    return FmpzConstRef(::fmpz_mat_entry(matrix.raw(), row, col));
}

inline FmpzRef fmpz_mat_entry(FmpzMatRef matrix, slong row, slong col) noexcept {
    return FmpzRef(::fmpz_mat_entry(matrix.raw(), row, col));
}

inline FmpzRef fmpz_mat_entry(FmpzMat& matrix, slong row, slong col) noexcept {
    return FmpzRef(::fmpz_mat_entry(matrix.raw(), row, col));
}

inline void fmpz_mat_set(FmpzMatRef out, FmpzMatConstRef in) noexcept {
    ::fmpz_mat_set(out.raw(), in.raw());
}

inline void fmpz_mat_neg(FmpzMatRef out, FmpzMatConstRef in) noexcept {
    ::fmpz_mat_neg(out.raw(), in.raw());
}

inline void fmpz_mat_mul(FmpzMatRef out,
                         FmpzMatConstRef left,
                         FmpzMatConstRef right) noexcept {
    ::fmpz_mat_mul(out.raw(), left.raw(), right.raw());
}

inline void fmpz_mat_transpose(FmpzMatRef out,
                               FmpzMatConstRef in) noexcept {
    ::fmpz_mat_transpose(out.raw(), in.raw());
}

inline void fmpz_mat_gram(FmpzMatRef out,
                          FmpzMatConstRef in) noexcept {
    ::fmpz_mat_gram(out.raw(), in.raw());
}

inline void fmpz_mat_scalar_tdiv_q_2exp(FmpzMatRef out,
                                        FmpzMatConstRef in,
                                        ulong exponent) noexcept {
    ::fmpz_mat_scalar_tdiv_q_2exp(out.raw(), in.raw(), exponent);
}

inline void fmpz_mat_det(FmpzRef out, FmpzMatConstRef matrix) noexcept {
    ::fmpz_mat_det(out.raw(), matrix.raw());
}

inline bool fmpz_mat_equal(FmpzMatConstRef left, FmpzMatConstRef right) noexcept {
    return ::fmpz_mat_equal(left.raw(), right.raw()) != 0;
}

inline void fmpz_mat_swap(FmpzMatRef left, FmpzMatRef right) noexcept {
    ::fmpz_mat_swap(left.raw(), right.raw());
}

inline void fmpz_mat_snf_transform(FmpzMatRef snf,
                                   FmpzMatRef left_transform,
                                   FmpzMatRef right_transform,
                                   FmpzMatConstRef matrix) noexcept {
    ::fmpz_mat_snf_transform(snf.raw(), left_transform.raw(),
                             right_transform.raw(), matrix.raw());
}

inline bool fmpz_mat_inv(FmpzMatRef out,
                         FmpzRef den,
                         FmpzMatConstRef matrix) noexcept {
    return ::fmpz_mat_inv(out.raw(), den.raw(), matrix.raw()) != 0;
}

}  // namespace silex::flint
