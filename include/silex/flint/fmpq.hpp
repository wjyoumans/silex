#pragma once

#include <flint/fmpq.h>

#include <silex/flint/fmpz.hpp>

namespace silex::flint {

class FmpqRef;
class FmpqConstRef;

class Fmpq {
public:
    Fmpq() noexcept {
        fmpq_init(value_);
    }

    ~Fmpq() noexcept {
        fmpq_clear(value_);
    }

    Fmpq(const Fmpq&) = delete;
    Fmpq& operator=(const Fmpq&) = delete;

    Fmpq(Fmpq&& other) noexcept {
        fmpq_init(value_);
        fmpq_swap(value_, other.value_);
    }

    Fmpq& operator=(Fmpq&& other) noexcept {
        if (this != &other) {
            fmpq_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(Fmpq& other) noexcept {
        fmpq_swap(value_, other.value_);
    }

    fmpq_t& raw() noexcept {
        return value_;
    }

    const fmpq_t& raw() const noexcept {
        return value_;
    }

private:
    fmpq_t value_;
};

class FmpqConstRef {
public:
    explicit FmpqConstRef(const fmpq* value) noexcept
        : value_(value) {
    }

    explicit FmpqConstRef(const Fmpq& value) noexcept
        : value_(value.raw()) {
    }

    const fmpq* raw() const noexcept {
        return value_;
    }

private:
    const fmpq* value_;
};

class FmpqRef {
public:
    explicit FmpqRef(fmpq* value) noexcept
        : value_(value) {
    }

    explicit FmpqRef(Fmpq& value) noexcept
        : value_(value.raw()) {
    }

    fmpq* raw() noexcept {
        return value_;
    }

    const fmpq* raw() const noexcept {
        return value_;
    }

private:
    fmpq* value_;
};

inline void swap(Fmpq& left, Fmpq& right) noexcept {
    left.swap(right);
}

inline void fmpq_zero(FmpqRef out) noexcept {
    ::fmpq_zero(out.raw());
}

inline void fmpq_zero(Fmpq& out) noexcept {
    ::fmpq_zero(out.raw());
}

inline void fmpq_one(FmpqRef out) noexcept {
    ::fmpq_one(out.raw());
}

inline void fmpq_set_si(FmpqRef out, slong numerator, ulong denominator) noexcept {
    ::fmpq_set_si(out.raw(), numerator, denominator);
}

inline void fmpq_set_si(Fmpq& out, slong numerator, ulong denominator) noexcept {
    ::fmpq_set_si(out.raw(), numerator, denominator);
}

inline void fmpq_one(Fmpq& out) noexcept {
    ::fmpq_one(out.raw());
}

inline void fmpq_set(FmpqRef out, FmpqConstRef in) noexcept {
    ::fmpq_set(out.raw(), in.raw());
}

inline void fmpq_set(Fmpq& out, const Fmpq& in) noexcept {
    ::fmpq_set(out.raw(), in.raw());
}

inline void fmpq_set_fmpz(FmpqRef out, FmpzConstRef value) noexcept {
    ::fmpq_set_fmpz(out.raw(), value.raw());
}

inline void fmpq_set_fmpz(Fmpq& out, FmpzConstRef value) noexcept {
    ::fmpq_set_fmpz(out.raw(), value.raw());
}

inline void fmpq_set_fmpz_frac(FmpqRef out,
                               FmpzConstRef numerator,
                               FmpzConstRef denominator) noexcept {
    ::fmpq_set_fmpz_frac(out.raw(), numerator.raw(), denominator.raw());
}

inline void fmpq_set_fmpz_frac(Fmpq& out,
                               FmpzConstRef numerator,
                               FmpzConstRef denominator) noexcept {
    ::fmpq_set_fmpz_frac(out.raw(), numerator.raw(), denominator.raw());
}

inline void fmpq_div_2exp(FmpqRef out, FmpqConstRef in, ulong exponent) noexcept {
    ::fmpq_div_2exp(out.raw(), in.raw(), exponent);
}

inline void fmpq_div_2exp(Fmpq& out, const Fmpq& in, ulong exponent) noexcept {
    ::fmpq_div_2exp(out.raw(), in.raw(), exponent);
}

inline bool fmpq_equal(FmpqConstRef left, FmpqConstRef right) noexcept {
    return ::fmpq_equal(left.raw(), right.raw()) != 0;
}

inline bool fmpq_equal(const Fmpq& left, const Fmpq& right) noexcept {
    return ::fmpq_equal(left.raw(), right.raw()) != 0;
}

inline bool fmpq_equal_si(FmpqConstRef value, slong expected) noexcept {
    return ::fmpq_cmp_si(value.raw(), expected) == 0;
}

inline bool fmpq_equal_si(const Fmpq& value, slong expected) noexcept {
    return ::fmpq_cmp_si(value.raw(), expected) == 0;
}

inline void fmpq_neg(FmpqRef out, FmpqConstRef in) noexcept {
    ::fmpq_neg(out.raw(), in.raw());
}

inline void fmpq_neg(Fmpq& out, const Fmpq& in) noexcept {
    ::fmpq_neg(out.raw(), in.raw());
}

inline void fmpq_abs(FmpqRef out, FmpqConstRef in) noexcept {
    ::fmpq_abs(out.raw(), in.raw());
}

inline void fmpq_abs(Fmpq& out, const Fmpq& in) noexcept {
    ::fmpq_abs(out.raw(), in.raw());
}

inline FmpzConstRef fmpq_num_ref(FmpqConstRef value) noexcept {
    return FmpzConstRef(fmpq_numref(value.raw()));
}

inline FmpzConstRef fmpq_num_ref(const Fmpq& value) noexcept {
    return FmpzConstRef(fmpq_numref(value.raw()));
}

inline FmpzConstRef fmpq_den_ref(FmpqConstRef value) noexcept {
    return FmpzConstRef(fmpq_denref(value.raw()));
}

inline FmpzConstRef fmpq_den_ref(const Fmpq& value) noexcept {
    return FmpzConstRef(fmpq_denref(value.raw()));
}

inline void fmpq_div(FmpqRef out,
                     FmpqConstRef numerator,
                     FmpqConstRef denominator) noexcept {
    ::fmpq_div(out.raw(), numerator.raw(), denominator.raw());
}

inline void fmpq_div(Fmpq& out,
                     const Fmpq& numerator,
                     const Fmpq& denominator) noexcept {
    ::fmpq_div(out.raw(), numerator.raw(), denominator.raw());
}

}  // namespace silex::flint
