#pragma once

#include <flint/nf.h>

#include <cstdlib>
#include <utility>

namespace silex::flint {

class NfRef;
class NfConstRef;

class Nf {
public:
    Nf() noexcept = default;

    explicit Nf(const fmpq_poly_t polynomial) noexcept {
        define(polynomial);
    }

    ~Nf() noexcept {
        clear();
    }

    Nf(const Nf&) = delete;
    Nf& operator=(const Nf&) = delete;
    Nf(Nf&& other) noexcept
        : value_(other.value_) {
        other.value_ = nullptr;
    }

    Nf& operator=(Nf&& other) noexcept {
        if (this != &other) {
            clear();
            value_ = other.value_;
            other.value_ = nullptr;
        }
        return *this;
    }

    bool define(const fmpq_poly_t polynomial) noexcept {
        nf_struct* next = static_cast<nf_struct*>(std::malloc(sizeof(nf_struct)));
        if (next == nullptr) {
            return false;
        }

        nf_init(next, polynomial);
        clear();
        value_ = next;
        return true;
    }

    void clear() noexcept {
        if (value_ != nullptr) {
            nf_clear(value_);
            std::free(value_);
            value_ = nullptr;
        }
    }

    void swap(Nf& other) noexcept {
        std::swap(value_, other.value_);
    }

    bool is_defined() const noexcept {
        return value_ != nullptr;
    }

    nf_struct* raw() noexcept {
        return value_;
    }

    const nf_struct* raw() const noexcept {
        return value_;
    }

private:
    nf_struct* value_ = nullptr;
};

class NfConstRef {
public:
    explicit NfConstRef(const nf_t& value) noexcept
        : value_(value) {
    }

    explicit NfConstRef(const Nf& value) noexcept
        : value_(value.raw()) {
    }

    const nf_struct* raw() const noexcept {
        return value_;
    }

private:
    const nf_struct* value_;
};

class NfRef {
public:
    explicit NfRef(nf_t& value) noexcept
        : value_(value) {
    }

    explicit NfRef(Nf& value) noexcept
        : value_(value.raw()) {
    }

    nf_struct* raw() noexcept {
        return value_;
    }

    const nf_struct* raw() const noexcept {
        return value_;
    }

private:
    nf_struct* value_;
};

inline void swap(Nf& left, Nf& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
