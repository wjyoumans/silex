#pragma once

#include <flint/nmod_mat.h>

namespace silex::flint {

class NmodMatRef;
class NmodMatConstRef;

class NmodMat {
public:
    NmodMat(slong rows, slong cols, ulong modulus) noexcept {
        nmod_mat_init(value_, rows, cols, modulus);
    }

    ~NmodMat() noexcept { nmod_mat_clear(value_); }

    NmodMat(const NmodMat&) = delete;
    NmodMat& operator=(const NmodMat&) = delete;

    NmodMat(NmodMat&& other) noexcept : NmodMat(0, 0, 2) {
        nmod_mat_swap(value_, other.value_);
    }

    NmodMat& operator=(NmodMat&& other) noexcept {
        if (this != &other) {
            nmod_mat_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(NmodMat& other) noexcept { nmod_mat_swap(value_, other.value_); }
    nmod_mat_t& raw() noexcept { return value_; }
    const nmod_mat_t& raw() const noexcept { return value_; }

private:
    nmod_mat_t value_;
};

class NmodMatConstRef {
public:
    NmodMatConstRef(const nmod_mat_struct* value) noexcept : value_(value) {}
    NmodMatConstRef(const NmodMat& value) noexcept : value_(value.raw()) {}
    const nmod_mat_struct* raw() const noexcept { return value_; }

private:
    const nmod_mat_struct* value_;
};

class NmodMatRef {
public:
    NmodMatRef(nmod_mat_struct* value) noexcept : value_(value) {}
    NmodMatRef(NmodMat& value) noexcept : value_(value.raw()) {}
    nmod_mat_struct* raw() noexcept { return value_; }
    const nmod_mat_struct* raw() const noexcept { return value_; }

private:
    nmod_mat_struct* value_;
};

inline void swap(NmodMat& left, NmodMat& right) noexcept { left.swap(right); }

inline slong nmod_mat_nrows(NmodMatConstRef matrix) noexcept {
    return ::nmod_mat_nrows(matrix.raw());
}

inline slong nmod_mat_nrows(const NmodMat& matrix) noexcept {
    return ::nmod_mat_nrows(matrix.raw());
}

inline slong nmod_mat_nrows(NmodMatRef matrix) noexcept {
    return ::nmod_mat_nrows(matrix.raw());
}

inline slong nmod_mat_ncols(NmodMatConstRef matrix) noexcept {
    return ::nmod_mat_ncols(matrix.raw());
}

inline slong nmod_mat_ncols(const NmodMat& matrix) noexcept {
    return ::nmod_mat_ncols(matrix.raw());
}

inline slong nmod_mat_ncols(NmodMatRef matrix) noexcept {
    return ::nmod_mat_ncols(matrix.raw());
}

inline void nmod_mat_zero(NmodMatRef matrix) noexcept {
    ::nmod_mat_zero(matrix.raw());
}

inline ulong nmod_mat_get_entry(NmodMatConstRef matrix, slong row, slong col) noexcept {
    return nmod_mat_entry(matrix.raw(), row, col);
}

inline ulong nmod_mat_get_entry(const NmodMat& matrix, slong row, slong col) noexcept {
    return nmod_mat_entry(matrix.raw(), row, col);
}

inline ulong nmod_mat_get_entry(NmodMatRef matrix, slong row, slong col) noexcept {
    return nmod_mat_entry(matrix.raw(), row, col);
}

inline void nmod_mat_set_entry(
        NmodMatRef matrix, slong row, slong col, ulong value) noexcept {
    nmod_mat_entry(matrix.raw(), row, col) = value;
}

inline void nmod_mat_set_entry(NmodMat& matrix, slong row, slong col, ulong value) noexcept {
    nmod_mat_entry(matrix.raw(), row, col) = value;
}

}  // namespace silex::flint
