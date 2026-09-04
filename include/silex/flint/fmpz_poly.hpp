#pragma once

#include <flint/fmpz_poly.h>

namespace silex::flint {

class FmpzPolyRef;
class FmpzPolyConstRef;

class FmpzPoly {
public:
    FmpzPoly() noexcept {
        fmpz_poly_init(value_);
    }

    ~FmpzPoly() noexcept {
        fmpz_poly_clear(value_);
    }

    FmpzPoly(const FmpzPoly&) = delete;
    FmpzPoly& operator=(const FmpzPoly&) = delete;

    FmpzPoly(FmpzPoly&& other) noexcept {
        fmpz_poly_init(value_);
        fmpz_poly_swap(value_, other.value_);
    }

    FmpzPoly& operator=(FmpzPoly&& other) noexcept {
        if (this != &other) {
            fmpz_poly_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(FmpzPoly& other) noexcept {
        fmpz_poly_swap(value_, other.value_);
    }

    fmpz_poly_t& raw() noexcept {
        return value_;
    }

    const fmpz_poly_t& raw() const noexcept {
        return value_;
    }

private:
    fmpz_poly_t value_;
};

class FmpzPolyConstRef {
public:
    explicit FmpzPolyConstRef(const fmpz_poly_t& value) noexcept
        : value_(value) {
    }

    explicit FmpzPolyConstRef(const FmpzPoly& value) noexcept
        : value_(value.raw()) {
    }

    const fmpz_poly_struct* raw() const noexcept {
        return value_;
    }

private:
    const fmpz_poly_struct* value_;
};

class FmpzPolyRef {
public:
    explicit FmpzPolyRef(fmpz_poly_t& value) noexcept
        : value_(value) {
    }

    explicit FmpzPolyRef(FmpzPoly& value) noexcept
        : value_(value.raw()) {
    }

    fmpz_poly_struct* raw() noexcept {
        return value_;
    }

    const fmpz_poly_struct* raw() const noexcept {
        return value_;
    }

private:
    fmpz_poly_struct* value_;
};

inline void swap(FmpzPoly& left, FmpzPoly& right) noexcept {
    left.swap(right);
}

inline void fmpz_poly_set_coeff_si(
        FmpzPolyRef polynomial, slong index, slong value) noexcept {
    ::fmpz_poly_set_coeff_si(polynomial.raw(), index, value);
}

inline void fmpz_poly_set_coeff_si(
        FmpzPoly& polynomial, slong index, slong value) noexcept {
    ::fmpz_poly_set_coeff_si(polynomial.raw(), index, value);
}

inline slong fmpz_poly_degree(FmpzPolyConstRef polynomial) noexcept {
    return ::fmpz_poly_degree(polynomial.raw());
}

inline slong fmpz_poly_degree(const FmpzPoly& polynomial) noexcept {
    return ::fmpz_poly_degree(polynomial.raw());
}

}  // namespace silex::flint
