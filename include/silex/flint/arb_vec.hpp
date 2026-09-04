#pragma once

#include <flint/arb.h>

#include <utility>

namespace silex::flint {

class ArbVecRef;
class ArbVecConstRef;

class ArbVec {
public:
    explicit ArbVec(slong length) noexcept
        : data_(length > 0 ? _arb_vec_init(length) : nullptr),
          length_(length > 0 ? length : 0) {
    }

    ~ArbVec() noexcept {
        if (data_ != nullptr) {
            _arb_vec_clear(data_, length_);
        }
    }

    ArbVec(const ArbVec&) = delete;
    ArbVec& operator=(const ArbVec&) = delete;

    ArbVec(ArbVec&& other) noexcept
        : data_(other.data_),
          length_(other.length_) {
        other.data_ = nullptr;
        other.length_ = 0;
    }

    ArbVec& operator=(ArbVec&& other) noexcept {
        if (this != &other) {
            ArbVec tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(ArbVec& other) noexcept {
        arb_ptr const data = data_;
        const slong length = length_;
        data_ = other.data_;
        length_ = other.length_;
        other.data_ = data;
        other.length_ = length;
    }

    bool set_from(ArbVecConstRef source) noexcept;

    arb_ptr data() noexcept { return data_; }
    arb_srcptr data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    arb_ptr data_;
    slong length_;
};

class ArbVecConstRef {
public:
    ArbVecConstRef(arb_srcptr data, slong length) noexcept
        : data_(data),
          length_(length) {
    }

    explicit ArbVecConstRef(const ArbVec& value) noexcept
        : data_(value.data()),
          length_(value.length()) {
    }

    arb_srcptr data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    arb_srcptr data_;
    slong length_;
};

class ArbVecRef {
public:
    ArbVecRef(arb_ptr data, slong length) noexcept
        : data_(data),
          length_(length) {
    }

    explicit ArbVecRef(ArbVec& value) noexcept
        : data_(value.data()),
          length_(value.length()) {
    }

    bool set_from(ArbVecConstRef source) noexcept {
        if (length_ != source.length()) {
            return false;
        }
        _arb_vec_set(data_, source.data(), length_);
        return true;
    }

    arb_ptr data() noexcept { return data_; }
    arb_srcptr data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    arb_ptr data_;
    slong length_;
};

inline bool ArbVec::set_from(ArbVecConstRef source) noexcept {
    return ArbVecRef(*this).set_from(source);
}

inline void swap(ArbVec& left, ArbVec& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
