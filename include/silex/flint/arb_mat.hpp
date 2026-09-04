#pragma once

#include <silex/flint/arb.hpp>

#include <flint/arb_mat.h>

namespace silex::flint {

class ArbMatRef;
class ArbMatConstRef;

class ArbMat {
public:
    ArbMat(slong rows, slong cols) noexcept { arb_mat_init(value_, rows, cols); }
    ~ArbMat() noexcept { arb_mat_clear(value_); }

    ArbMat(const ArbMat&) = delete;
    ArbMat& operator=(const ArbMat&) = delete;

    ArbMat(ArbMat&& other) noexcept : ArbMat(0, 0) {
        arb_mat_swap(value_, other.value_);
    }

    ArbMat& operator=(ArbMat&& other) noexcept {
        if (this != &other) {
            arb_mat_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(ArbMat& other) noexcept { arb_mat_swap(value_, other.value_); }
    arb_mat_t& raw() noexcept { return value_; }
    const arb_mat_t& raw() const noexcept { return value_; }

private:
    arb_mat_t value_;
};

class ArbMatConstRef {
public:
    explicit ArbMatConstRef(const arb_mat_t& value) noexcept : value_(value) {}
    explicit ArbMatConstRef(const ArbMat& value) noexcept : value_(value.raw()) {}
    const arb_mat_struct* raw() const noexcept { return value_; }

private:
    const arb_mat_struct* value_;
};

class ArbMatRef {
public:
    explicit ArbMatRef(arb_mat_t& value) noexcept : value_(value) {}
    explicit ArbMatRef(ArbMat& value) noexcept : value_(value.raw()) {}
    arb_mat_struct* raw() noexcept { return value_; }
    const arb_mat_struct* raw() const noexcept { return value_; }

private:
    arb_mat_struct* value_;
};

inline void swap(ArbMat& left, ArbMat& right) noexcept { left.swap(right); }

inline slong arb_mat_ncols_value(ArbMatConstRef matrix) noexcept {
    return arb_mat_ncols(matrix.raw());
}

inline slong arb_mat_nrows_value(ArbMatConstRef matrix) noexcept {
    return arb_mat_nrows(matrix.raw());
}

inline slong arb_mat_nrows_value(const ArbMat& matrix) noexcept {
    return arb_mat_nrows(matrix.raw());
}

inline slong arb_mat_ncols_value(const ArbMat& matrix) noexcept {
    return arb_mat_ncols(matrix.raw());
}

inline ArbConstRef arb_mat_entry_ref(
        ArbMatConstRef matrix, slong row, slong col) noexcept {
    return ArbConstRef(arb_mat_entry(matrix.raw(), row, col));
}

inline ArbConstRef arb_mat_entry_ref(
        const ArbMat& matrix, slong row, slong col) noexcept {
    return ArbConstRef(arb_mat_entry(matrix.raw(), row, col));
}

inline ArbRef arb_mat_entry_ref(ArbMatRef matrix, slong row, slong col) noexcept {
    return ArbRef(arb_mat_entry(matrix.raw(), row, col));
}

inline ArbRef arb_mat_entry_ref(ArbMat& matrix, slong row, slong col) noexcept {
    return ArbRef(arb_mat_entry(matrix.raw(), row, col));
}

inline void arb_mat_set(ArbMatRef out, ArbMatConstRef in) noexcept {
    ::arb_mat_set(out.raw(), in.raw());
}

inline slong arb_mat_nrows_value(ArbMatRef matrix) noexcept {
    return arb_mat_nrows(matrix.raw());
}

inline slong arb_mat_ncols_value(ArbMatRef matrix) noexcept {
    return arb_mat_ncols(matrix.raw());
}

inline int arb_mat_cho(ArbMat& out, const ArbMat& in, slong prec) noexcept {
    return ::arb_mat_cho(out.raw(), in.raw(), prec);
}

}  // namespace silex::flint
