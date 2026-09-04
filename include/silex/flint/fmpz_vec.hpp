#pragma once

#include <flint/fmpz_vec.h>

#include <utility>

namespace silex::flint {

class FmpzVecRef;
class FmpzVecConstRef;

class FmpzVec {
public:
    explicit FmpzVec(slong length) noexcept
        : data_(length > 0 ? _fmpz_vec_init(length) : nullptr), length_(length) {
    }

    ~FmpzVec() noexcept {
        if (data_ != nullptr) {
            _fmpz_vec_clear(data_, length_);
        }
    }

    FmpzVec(const FmpzVec&) = delete;
    FmpzVec& operator=(const FmpzVec&) = delete;

    FmpzVec(FmpzVec&& other) noexcept
        : data_(other.data_), length_(other.length_) {
        other.data_ = nullptr;
        other.length_ = 0;
    }

    FmpzVec& operator=(FmpzVec&& other) noexcept {
        if (this != &other) {
            FmpzVec tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(FmpzVec& other) noexcept {
        fmpz* const data = data_;
        const slong length = length_;
        data_ = other.data_;
        length_ = other.length_;
        other.data_ = data;
        other.length_ = length;
    }

    fmpz* data() noexcept { return data_; }
    const fmpz* data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    fmpz* data_;
    slong length_;
};

class FmpzVecConstRef {
public:
    FmpzVecConstRef(const fmpz* data, slong length) noexcept
        : data_(data), length_(length) {
    }

    explicit FmpzVecConstRef(const FmpzVec& value) noexcept
        : data_(value.data()), length_(value.length()) {
    }

    const fmpz* data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    const fmpz* data_;
    slong length_;
};

class FmpzVecRef {
public:
    FmpzVecRef(fmpz* data, slong length) noexcept : data_(data), length_(length) {}

    explicit FmpzVecRef(FmpzVec& value) noexcept
        : data_(value.data()), length_(value.length()) {
    }

    fmpz* data() noexcept { return data_; }
    const fmpz* data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    fmpz* data_;
    slong length_;
};

inline void swap(FmpzVec& left, FmpzVec& right) noexcept { left.swap(right); }

}  // namespace silex::flint
