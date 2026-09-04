#pragma once

#include <flint/acb.h>

#include <utility>

namespace silex::flint {

class AcbVecRef;
class AcbVecConstRef;

class AcbVec {
public:
    explicit AcbVec(slong length) noexcept
        : data_(length > 0 ? _acb_vec_init(length) : nullptr),
          length_(length > 0 ? length : 0) {
    }

    ~AcbVec() noexcept {
        if (data_ != nullptr) {
            _acb_vec_clear(data_, length_);
        }
    }

    AcbVec(const AcbVec&) = delete;
    AcbVec& operator=(const AcbVec&) = delete;

    AcbVec(AcbVec&& other) noexcept
        : data_(other.data_),
          length_(other.length_) {
        other.data_ = nullptr;
        other.length_ = 0;
    }

    AcbVec& operator=(AcbVec&& other) noexcept {
        if (this != &other) {
            AcbVec tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(AcbVec& other) noexcept {
        acb_ptr const data = data_;
        const slong length = length_;
        data_ = other.data_;
        length_ = other.length_;
        other.data_ = data;
        other.length_ = length;
    }

    bool set_from(AcbVecConstRef source) noexcept;

    acb_ptr data() noexcept { return data_; }
    acb_srcptr data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    acb_ptr data_;
    slong length_;
};

class AcbVecConstRef {
public:
    AcbVecConstRef(acb_srcptr data, slong length) noexcept
        : data_(data),
          length_(length) {
    }

    explicit AcbVecConstRef(const AcbVec& value) noexcept
        : data_(value.data()),
          length_(value.length()) {
    }

    acb_srcptr data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    acb_srcptr data_;
    slong length_;
};

class AcbVecRef {
public:
    AcbVecRef(acb_ptr data, slong length) noexcept
        : data_(data),
          length_(length) {
    }

    explicit AcbVecRef(AcbVec& value) noexcept
        : data_(value.data()),
          length_(value.length()) {
    }

    bool set_from(AcbVecConstRef source) noexcept {
        if (length_ != source.length()) {
            return false;
        }
        _acb_vec_set(data_, source.data(), length_);
        return true;
    }

    acb_ptr data() noexcept { return data_; }
    acb_srcptr data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    acb_ptr data_;
    slong length_;
};

inline bool AcbVec::set_from(AcbVecConstRef source) noexcept {
    return AcbVecRef(*this).set_from(source);
}

inline void swap(AcbVec& left, AcbVec& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
