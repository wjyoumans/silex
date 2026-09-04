#pragma once

#include <flint/nmod_poly.h>

namespace silex::flint {

class NmodPolyRef;
class NmodPolyConstRef;

class NmodPoly {
public:
    explicit NmodPoly(ulong modulus) noexcept {
        nmod_poly_init(value_, modulus);
    }

    ~NmodPoly() noexcept {
        nmod_poly_clear(value_);
    }

    NmodPoly(const NmodPoly&) = delete;
    NmodPoly& operator=(const NmodPoly&) = delete;

    NmodPoly(NmodPoly&& other) noexcept : NmodPoly(2) {
        nmod_poly_swap(value_, other.value_);
    }

    NmodPoly& operator=(NmodPoly&& other) noexcept {
        if (this != &other) {
            nmod_poly_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(NmodPoly& other) noexcept {
        nmod_poly_swap(value_, other.value_);
    }

    nmod_poly_t& raw() noexcept {
        return value_;
    }

    const nmod_poly_t& raw() const noexcept {
        return value_;
    }

private:
    nmod_poly_t value_;
};

class NmodPolyConstRef {
public:
    explicit NmodPolyConstRef(const nmod_poly_t& value) noexcept
        : value_(value) {
    }

    explicit NmodPolyConstRef(const NmodPoly& value) noexcept
        : value_(value.raw()) {
    }

    const nmod_poly_struct* raw() const noexcept {
        return value_;
    }

private:
    const nmod_poly_struct* value_;
};

class NmodPolyRef {
public:
    explicit NmodPolyRef(nmod_poly_t& value) noexcept
        : value_(value) {
    }

    explicit NmodPolyRef(NmodPoly& value) noexcept
        : value_(value.raw()) {
    }

    nmod_poly_struct* raw() noexcept {
        return value_;
    }

    const nmod_poly_struct* raw() const noexcept {
        return value_;
    }

private:
    nmod_poly_struct* value_;
};

inline void swap(NmodPoly& left, NmodPoly& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
