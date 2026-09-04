#pragma once

#include <silex/flint/acb.hpp>

#include <flint/acb_mat.h>

namespace silex::flint {

class AcbMatRef;
class AcbMatConstRef;

class AcbMat {
public:
    AcbMat(slong rows, slong cols) noexcept { acb_mat_init(value_, rows, cols); }
    ~AcbMat() noexcept { acb_mat_clear(value_); }

    AcbMat(const AcbMat&) = delete;
    AcbMat& operator=(const AcbMat&) = delete;

    AcbMat(AcbMat&& other) noexcept : AcbMat(0, 0) {
        acb_mat_swap(value_, other.value_);
    }

    AcbMat& operator=(AcbMat&& other) noexcept {
        if (this != &other) {
            acb_mat_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(AcbMat& other) noexcept { acb_mat_swap(value_, other.value_); }
    acb_mat_t& raw() noexcept { return value_; }
    const acb_mat_t& raw() const noexcept { return value_; }

private:
    acb_mat_t value_;
};

class AcbMatConstRef {
public:
    explicit AcbMatConstRef(const acb_mat_t& value) noexcept : value_(value) {}
    explicit AcbMatConstRef(const AcbMat& value) noexcept : value_(value.raw()) {}
    const acb_mat_struct* raw() const noexcept { return value_; }

private:
    const acb_mat_struct* value_;
};

class AcbMatRef {
public:
    explicit AcbMatRef(acb_mat_t& value) noexcept : value_(value) {}
    explicit AcbMatRef(AcbMat& value) noexcept : value_(value.raw()) {}
    acb_mat_struct* raw() noexcept { return value_; }
    const acb_mat_struct* raw() const noexcept { return value_; }

private:
    acb_mat_struct* value_;
};

inline void swap(AcbMat& left, AcbMat& right) noexcept { left.swap(right); }

inline slong acb_mat_nrows_value(AcbMatConstRef matrix) noexcept {
    return acb_mat_nrows(matrix.raw());
}

inline slong acb_mat_ncols_value(AcbMatConstRef matrix) noexcept {
    return acb_mat_ncols(matrix.raw());
}

inline slong acb_mat_nrows_value(const AcbMat& matrix) noexcept {
    return acb_mat_nrows(matrix.raw());
}

inline slong acb_mat_ncols_value(const AcbMat& matrix) noexcept {
    return acb_mat_ncols(matrix.raw());
}

inline AcbRef acb_mat_entry_ref(AcbMatRef matrix,
                                slong row,
                                slong col) noexcept {
    return AcbRef(acb_mat_entry(matrix.raw(), row, col));
}

inline AcbRef acb_mat_entry_ref(AcbMat& matrix,
                                slong row,
                                slong col) noexcept {
    return AcbRef(acb_mat_entry(matrix.raw(), row, col));
}

inline AcbConstRef acb_mat_entry_ref(AcbMatConstRef matrix,
                                     slong row,
                                     slong col) noexcept {
    return AcbConstRef(acb_mat_entry(matrix.raw(), row, col));
}

inline AcbConstRef acb_mat_entry_ref(const AcbMat& matrix,
                                     slong row,
                                     slong col) noexcept {
    return AcbConstRef(acb_mat_entry(matrix.raw(), row, col));
}

inline int acb_mat_solve(AcbMat& out,
                         const AcbMat& matrix,
                         const AcbMat& rhs,
                         slong precision) noexcept {
    return ::acb_mat_solve(out.raw(), matrix.raw(), rhs.raw(), precision);
}

}  // namespace silex::flint
