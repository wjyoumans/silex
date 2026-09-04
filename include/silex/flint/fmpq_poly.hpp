#pragma once

#include <flint/fmpq_poly.h>

#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpq.hpp>

namespace silex::flint {

class FmpqPolyRef;
class FmpqPolyConstRef;

class FmpqPoly {
public:
    FmpqPoly() noexcept {
        fmpq_poly_init(value_);
    }

    ~FmpqPoly() noexcept {
        fmpq_poly_clear(value_);
    }

    FmpqPoly(const FmpqPoly&) = delete;
    FmpqPoly& operator=(const FmpqPoly&) = delete;

    FmpqPoly(FmpqPoly&& other) noexcept {
        fmpq_poly_init(value_);
        fmpq_poly_swap(value_, other.value_);
    }

    FmpqPoly& operator=(FmpqPoly&& other) noexcept {
        if (this != &other) {
            fmpq_poly_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(FmpqPoly& other) noexcept {
        fmpq_poly_swap(value_, other.value_);
    }

    fmpq_poly_t& raw() noexcept {
        return value_;
    }

    const fmpq_poly_t& raw() const noexcept {
        return value_;
    }

private:
    fmpq_poly_t value_;
};

class FmpqPolyConstRef {
public:
    explicit FmpqPolyConstRef(const fmpq_poly_t& value) noexcept
        : value_(value) {
    }

    explicit FmpqPolyConstRef(const FmpqPoly& value) noexcept
        : value_(value.raw()) {
    }

    const fmpq_poly_struct* raw() const noexcept {
        return value_;
    }

private:
    const fmpq_poly_struct* value_;
};

class FmpqPolyRef {
public:
    explicit FmpqPolyRef(fmpq_poly_t& value) noexcept
        : value_(value) {
    }

    explicit FmpqPolyRef(FmpqPoly& value) noexcept
        : value_(value.raw()) {
    }

    fmpq_poly_struct* raw() noexcept {
        return value_;
    }

    const fmpq_poly_struct* raw() const noexcept {
        return value_;
    }

private:
    fmpq_poly_struct* value_;
};

inline void swap(FmpqPoly& left, FmpqPoly& right) noexcept {
    left.swap(right);
}

inline void fmpq_poly_zero(FmpqPolyRef polynomial) noexcept {
    ::fmpq_poly_zero(polynomial.raw());
}

inline void fmpq_poly_zero(FmpqPoly& polynomial) noexcept {
    ::fmpq_poly_zero(polynomial.raw());
}

inline void fmpq_poly_set_coeff_si(
        FmpqPolyRef polynomial, slong index, slong value) noexcept {
    ::fmpq_poly_set_coeff_si(polynomial.raw(), index, value);
}

inline void fmpq_poly_set_coeff_si(
        FmpqPoly& polynomial, slong index, slong value) noexcept {
    ::fmpq_poly_set_coeff_si(polynomial.raw(), index, value);
}

inline void fmpq_poly_set_coeff_fmpq(
        FmpqPolyRef polynomial, slong index, FmpqConstRef value) noexcept {
    ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), index, value.raw());
}

inline void fmpq_poly_set_coeff_fmpq(
        FmpqPoly& polynomial, slong index, FmpqConstRef value) noexcept {
    ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), index, value.raw());
}

inline void fmpq_poly_set_coeff_fmpq(
        FmpqPoly& polynomial, slong index, const Fmpq& value) noexcept {
    ::fmpq_poly_set_coeff_fmpq(polynomial.raw(), index, value.raw());
}

inline slong fmpq_poly_degree(FmpqPolyConstRef polynomial) noexcept {
    return ::fmpq_poly_degree(polynomial.raw());
}

inline slong fmpq_poly_degree(const FmpqPoly& polynomial) noexcept {
    return ::fmpq_poly_degree(polynomial.raw());
}

inline void fmpq_poly_get_coeff_fmpq(
        FmpqRef out, FmpqPolyConstRef polynomial, slong index) noexcept {
    ::fmpq_poly_get_coeff_fmpq(out.raw(), polynomial.raw(), index);
}

inline void fmpq_poly_get_coeff_fmpq(
        FmpqRef out, const FmpqPoly& polynomial, slong index) noexcept {
    ::fmpq_poly_get_coeff_fmpq(out.raw(), polynomial.raw(), index);
}

inline void fmpq_poly_get_coeff_fmpz(
        FmpzRef out, FmpqPolyConstRef polynomial, slong index) noexcept {
    ::fmpq_poly_get_coeff_fmpz(out.raw(), polynomial.raw(), index);
}

inline void fmpq_poly_get_coeff_fmpz(
        FmpzRef out, const FmpqPoly& polynomial, slong index) noexcept {
    ::fmpq_poly_get_coeff_fmpz(out.raw(), polynomial.raw(), index);
}

inline void fmpq_poly_get_numerator(
        FmpzPolyRef out, FmpqPolyConstRef polynomial) noexcept {
    ::fmpq_poly_get_numerator(out.raw(), polynomial.raw());
}

inline void fmpq_poly_get_numerator(
        FmpzPoly& out, const FmpqPoly& polynomial) noexcept {
    ::fmpq_poly_get_numerator(out.raw(), polynomial.raw());
}

inline void fmpq_poly_scalar_div_ui(
        FmpqPolyRef out, FmpqPolyConstRef in, ulong value) noexcept {
    ::fmpq_poly_scalar_div_ui(out.raw(), in.raw(), value);
}

inline void fmpq_poly_scalar_div_ui(
        FmpqPoly& out, const FmpqPoly& in, ulong value) noexcept {
    ::fmpq_poly_scalar_div_ui(out.raw(), in.raw(), value);
}

inline void fmpq_poly_scalar_div_fmpz(
        FmpqPolyRef out, FmpqPolyConstRef in, FmpzConstRef value) noexcept {
    ::fmpq_poly_scalar_div_fmpz(out.raw(), in.raw(), value.raw());
}

inline void fmpq_poly_scalar_div_fmpz(
        FmpqPoly& out, const FmpqPoly& in, FmpzConstRef value) noexcept {
    ::fmpq_poly_scalar_div_fmpz(out.raw(), in.raw(), value.raw());
}

}  // namespace silex::flint
