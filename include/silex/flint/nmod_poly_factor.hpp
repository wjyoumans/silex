#pragma once

#include <flint/nmod_poly_factor.h>

namespace silex::flint {

class NmodPolyFactorRef;
class NmodPolyFactorConstRef;

class NmodPolyFactor {
public:
    NmodPolyFactor() noexcept {
        nmod_poly_factor_init(value_);
    }

    ~NmodPolyFactor() noexcept {
        nmod_poly_factor_clear(value_);
    }

    NmodPolyFactor(const NmodPolyFactor&) = delete;
    NmodPolyFactor& operator=(const NmodPolyFactor&) = delete;

    NmodPolyFactor(NmodPolyFactor&& other) noexcept : NmodPolyFactor() {
        nmod_poly_factor_swap(value_, other.value_);
    }

    NmodPolyFactor& operator=(NmodPolyFactor&& other) noexcept {
        if (this != &other) {
            nmod_poly_factor_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(NmodPolyFactor& other) noexcept {
        nmod_poly_factor_swap(value_, other.value_);
    }

    nmod_poly_factor_t& raw() noexcept {
        return value_;
    }

    const nmod_poly_factor_t& raw() const noexcept {
        return value_;
    }

private:
    nmod_poly_factor_t value_;
};

class NmodPolyFactorConstRef {
public:
    explicit NmodPolyFactorConstRef(
            const nmod_poly_factor_t& value) noexcept
        : value_(value) {
    }

    explicit NmodPolyFactorConstRef(const NmodPolyFactor& value) noexcept
        : value_(value.raw()) {
    }

    const nmod_poly_factor_struct* raw() const noexcept {
        return value_;
    }

private:
    const nmod_poly_factor_struct* value_;
};

class NmodPolyFactorRef {
public:
    explicit NmodPolyFactorRef(nmod_poly_factor_t& value) noexcept
        : value_(value) {
    }

    explicit NmodPolyFactorRef(NmodPolyFactor& value) noexcept
        : value_(value.raw()) {
    }

    nmod_poly_factor_struct* raw() noexcept {
        return value_;
    }

    const nmod_poly_factor_struct* raw() const noexcept {
        return value_;
    }

private:
    nmod_poly_factor_struct* value_;
};

inline void swap(NmodPolyFactor& left, NmodPolyFactor& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
