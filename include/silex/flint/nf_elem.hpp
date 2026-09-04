#pragma once

#include <silex/flint/nf.hpp>

#include <flint/nf_elem.h>

#include <cstdlib>
#include <utility>

namespace silex::flint {

class NfElemRef;
class NfElemConstRef;
class NfElemVecRef;
class NfElemVecConstRef;

class NfElem {
public:
    NfElem() noexcept = default;

    explicit NfElem(const nf_struct* field) noexcept
        : field_(field),
          value_((field == nullptr)
                         ? nullptr
                         : static_cast<nf_elem_struct*>(
                                   std::malloc(sizeof(nf_elem_struct)))) {
        if (value_ == nullptr) {
            field_ = nullptr;
            return;
        }
        nf_elem_init(value_, field_);
    }

    explicit NfElem(const nf_t& field) noexcept
        : NfElem(static_cast<const nf_struct*>(field)) {
    }

    explicit NfElem(const Nf& field) noexcept
        : NfElem(static_cast<const nf_struct*>(field.raw())) {
    }

    ~NfElem() noexcept { clear(); }

    NfElem(const NfElem&) = delete;
    NfElem& operator=(const NfElem&) = delete;
    NfElem(NfElem&& other) noexcept
        : field_(other.field_),
          value_(other.value_) {
        other.field_ = nullptr;
        other.value_ = nullptr;
    }

    NfElem& operator=(NfElem&& other) noexcept {
        if (this != &other) {
            clear();
            field_ = other.field_;
            value_ = other.value_;
            other.field_ = nullptr;
            other.value_ = nullptr;
        }
        return *this;
    }

    void clear() noexcept {
        if (value_ != nullptr) {
            nf_elem_clear(value_, field_);
            std::free(value_);
            value_ = nullptr;
            field_ = nullptr;
        }
    }

    void swap(NfElem& other) noexcept {
        std::swap(field_, other.field_);
        std::swap(value_, other.value_);
    }

    bool is_defined() const noexcept { return value_ != nullptr; }

    nf_elem_struct* raw() noexcept { return value_; }
    const nf_elem_struct* raw() const noexcept { return value_; }
    const nf_struct* field() const noexcept { return field_; }

private:
    const nf_struct* field_ = nullptr;
    nf_elem_struct* value_ = nullptr;
};

class NfElemConstRef {
public:
    explicit NfElemConstRef(const nf_elem_t& value) noexcept : value_(value) {}
    explicit NfElemConstRef(const NfElem& value) noexcept : value_(value.raw()) {}
    const nf_elem_struct* raw() const noexcept { return value_; }

private:
    const nf_elem_struct* value_;
};

class NfElemRef {
public:
    explicit NfElemRef(nf_elem_t& value) noexcept : value_(value) {}
    explicit NfElemRef(NfElem& value) noexcept : value_(value.raw()) {}
    nf_elem_struct* raw() noexcept { return value_; }
    const nf_elem_struct* raw() const noexcept { return value_; }

private:
    nf_elem_struct* value_;
};

inline void swap(NfElem& left, NfElem& right) noexcept {
    left.swap(right);
}

class NfElemVec {
public:
    NfElemVec(slong length, const nf_struct* field) noexcept
        : data_((length > 0 && field != nullptr)
                        ? static_cast<nf_elem_struct*>(
                                  flint_malloc(length * sizeof(nf_elem_struct)))
                        : nullptr),
          length_(data_ == nullptr ? 0 : length),
          field_(field) {
        for (slong i = 0; i < length_; ++i) {
            nf_elem_init(data_ + i, field_);
        }
    }

    explicit NfElemVec(slong length, const Nf& field) noexcept
        : NfElemVec(length, static_cast<const nf_struct*>(field.raw())) {
    }

    ~NfElemVec() noexcept {
        for (slong i = length_ - 1; i >= 0; --i) {
            nf_elem_clear(data_ + i, field_);
        }
        flint_free(data_);
    }

    NfElemVec(const NfElemVec&) = delete;
    NfElemVec& operator=(const NfElemVec&) = delete;
    NfElemVec(NfElemVec&&) = delete;
    NfElemVec& operator=(NfElemVec&&) = delete;

    nf_elem_struct* data() noexcept { return data_; }
    const nf_elem_struct* data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }
    const nf_struct* field() const noexcept { return field_; }

private:
    nf_elem_struct* data_;
    slong length_;
    const nf_struct* field_;
};

class NfElemVecConstRef {
public:
    NfElemVecConstRef(const nf_elem_struct* data, slong length) noexcept
        : data_(data), length_(length) {
    }

    explicit NfElemVecConstRef(const NfElemVec& value) noexcept
        : data_(value.data()), length_(value.length()) {
    }

    const nf_elem_struct* data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    const nf_elem_struct* data_;
    slong length_;
};

class NfElemVecRef {
public:
    NfElemVecRef(nf_elem_struct* data, slong length) noexcept
        : data_(data), length_(length) {
    }

    explicit NfElemVecRef(NfElemVec& value) noexcept
        : data_(value.data()), length_(value.length()) {
    }

    nf_elem_struct* data() noexcept { return data_; }
    const nf_elem_struct* data() const noexcept { return data_; }
    slong length() const noexcept { return length_; }

private:
    nf_elem_struct* data_;
    slong length_;
};

}  // namespace silex::flint
