#pragma once

#include <flint/fmpz.h>

namespace silex::flint {

class FmpzRef;
class FmpzConstRef;

class Fmpz {
public:
    Fmpz() noexcept {
        fmpz_init(value_);
    }

    ~Fmpz() noexcept {
        fmpz_clear(value_);
    }

    Fmpz(const Fmpz&) = delete;
    Fmpz& operator=(const Fmpz&) = delete;

    Fmpz(Fmpz&& other) noexcept {
        fmpz_init(value_);
        fmpz_swap(value_, other.value_);
    }

    Fmpz& operator=(Fmpz&& other) noexcept {
        if (this != &other) {
            fmpz_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(Fmpz& other) noexcept {
        fmpz_swap(value_, other.value_);
    }

    fmpz_t& raw() noexcept {
        return value_;
    }

    const fmpz_t& raw() const noexcept {
        return value_;
    }

private:
    fmpz_t value_;
};

class FmpzConstRef {
public:
    FmpzConstRef(const fmpz* value) noexcept
        : value_(value) {
    }

    FmpzConstRef(const Fmpz& value) noexcept
        : value_(value.raw()) {
    }

    const fmpz* raw() const noexcept {
        return value_;
    }

private:
    const fmpz* value_;
};

class FmpzRef {
public:
    FmpzRef(fmpz* value) noexcept
        : value_(value) {
    }

    FmpzRef(Fmpz& value) noexcept
        : value_(value.raw()) {
    }

    fmpz* raw() noexcept {
        return value_;
    }

    const fmpz* raw() const noexcept {
        return value_;
    }

private:
    fmpz* value_;
};

inline void swap(Fmpz& left, Fmpz& right) noexcept {
    left.swap(right);
}

inline void fmpz_zero(FmpzRef out) noexcept {
    ::fmpz_zero(out.raw());
}

inline void fmpz_one(FmpzRef out) noexcept {
    ::fmpz_one(out.raw());
}

inline void fmpz_set(FmpzRef out, FmpzConstRef in) noexcept {
    ::fmpz_set(out.raw(), in.raw());
}

inline void fmpz_set_si(FmpzRef out, slong value) noexcept {
    ::fmpz_set_si(out.raw(), value);
}

inline void fmpz_set_ui(FmpzRef out, ulong value) noexcept {
    ::fmpz_set_ui(out.raw(), value);
}

inline bool fmpz_set_str(FmpzRef out,
                         const char* value,
                         int base = 10) noexcept {
    return ::fmpz_set_str(out.raw(), value, base) == 0;
}

inline bool fmpz_fits_si(FmpzConstRef value) noexcept {
    return ::fmpz_fits_si(value.raw()) != 0;
}

inline slong fmpz_get_si(FmpzConstRef value) noexcept {
    return ::fmpz_get_si(value.raw());
}

inline ulong fmpz_get_ui(FmpzConstRef value) noexcept {
    return ::fmpz_get_ui(value.raw());
}

inline double fmpz_get_d(FmpzConstRef value) noexcept {
    return ::fmpz_get_d(value.raw());
}

inline double fmpz_get_d_2exp(slong& exponent, FmpzConstRef value) noexcept {
    return ::fmpz_get_d_2exp(&exponent, value.raw());
}

inline bool fmpz_abs_fits_ui(FmpzConstRef value) noexcept {
    return ::fmpz_abs_fits_ui(value.raw()) != 0;
}

inline void fmpz_neg(FmpzRef out, FmpzConstRef in) noexcept {
    ::fmpz_neg(out.raw(), in.raw());
}

inline void fmpz_mul(FmpzRef out, FmpzConstRef left, FmpzConstRef right) noexcept {
    ::fmpz_mul(out.raw(), left.raw(), right.raw());
}

inline void fmpz_mul_ui(FmpzRef out, FmpzConstRef in, ulong value) noexcept {
    ::fmpz_mul_ui(out.raw(), in.raw(), value);
}

inline void fmpz_pow_ui(FmpzRef out, FmpzConstRef in, ulong exponent) noexcept {
    ::fmpz_pow_ui(out.raw(), in.raw(), exponent);
}

inline void fmpz_add(FmpzRef out, FmpzConstRef left, FmpzConstRef right) noexcept {
    ::fmpz_add(out.raw(), left.raw(), right.raw());
}

inline void fmpz_add_ui(FmpzRef out, FmpzConstRef in, ulong value) noexcept {
    ::fmpz_add_ui(out.raw(), in.raw(), value);
}

inline void fmpz_sub(FmpzRef out, FmpzConstRef left, FmpzConstRef right) noexcept {
    ::fmpz_sub(out.raw(), left.raw(), right.raw());
}

inline void fmpz_sub_ui(FmpzRef out, FmpzConstRef in, ulong value) noexcept {
    ::fmpz_sub_ui(out.raw(), in.raw(), value);
}

inline void fmpz_addmul(FmpzRef out,
                        FmpzConstRef left,
                        FmpzConstRef right) noexcept {
    ::fmpz_addmul(out.raw(), left.raw(), right.raw());
}

inline void fmpz_addmul_ui(FmpzRef out,
                           FmpzConstRef left,
                           ulong right) noexcept {
    ::fmpz_addmul_ui(out.raw(), left.raw(), right);
}

inline void fmpz_submul(FmpzRef out,
                        FmpzConstRef left,
                        FmpzConstRef right) noexcept {
    ::fmpz_submul(out.raw(), left.raw(), right.raw());
}

inline void fmpz_submul_ui(FmpzRef out,
                           FmpzConstRef left,
                           ulong right) noexcept {
    ::fmpz_submul_ui(out.raw(), left.raw(), right);
}

inline void fmpz_lcm(FmpzRef out,
                     FmpzConstRef left,
                     FmpzConstRef right) noexcept {
    ::fmpz_lcm(out.raw(), left.raw(), right.raw());
}

inline int fmpz_cmp(FmpzConstRef left, FmpzConstRef right) noexcept {
    return ::fmpz_cmp(left.raw(), right.raw());
}

inline int fmpz_cmpabs(FmpzConstRef left, FmpzConstRef right) noexcept {
    return ::fmpz_cmpabs(left.raw(), right.raw());
}

inline int fmpz_cmp_ui(FmpzConstRef left, ulong right) noexcept {
    return ::fmpz_cmp_ui(left.raw(), right);
}

inline void fmpz_fdiv_r(
        FmpzRef out, FmpzConstRef value, FmpzConstRef modulus) noexcept {
    ::fmpz_fdiv_r(out.raw(), value.raw(), modulus.raw());
}

inline bool fmpz_is_zero(FmpzConstRef value) noexcept {
    return ::fmpz_is_zero(value.raw()) != 0;
}

inline bool fmpz_is_zero(const Fmpz& value) noexcept {
    return ::fmpz_is_zero(value.raw()) != 0;
}

inline bool fmpz_is_zero(FmpzRef value) noexcept {
    return ::fmpz_is_zero(value.raw()) != 0;
}

inline bool fmpz_is_prime(FmpzConstRef value) noexcept {
    return ::fmpz_is_prime(value.raw()) != 0;
}

inline void fmpz_nextprime(FmpzRef out,
                           FmpzConstRef value,
                           bool proved = true) noexcept {
    ::fmpz_nextprime(out.raw(), value.raw(), proved ? 1 : 0);
}

inline bool fmpz_is_one(FmpzConstRef value) noexcept {
    return ::fmpz_is_one(value.raw()) != 0;
}

inline bool fmpz_is_one(const Fmpz& value) noexcept {
    return ::fmpz_is_one(value.raw()) != 0;
}

inline bool fmpz_is_one(FmpzRef value) noexcept {
    return ::fmpz_is_one(value.raw()) != 0;
}

inline bool fmpz_is_square(FmpzConstRef value) noexcept {
    return ::fmpz_is_square(value.raw()) != 0;
}

inline bool fmpz_is_pm1(FmpzConstRef value) noexcept {
    return ::fmpz_is_pm1(value.raw()) != 0;
}

inline bool fmpz_is_pm1(const Fmpz& value) noexcept {
    return ::fmpz_is_pm1(value.raw()) != 0;
}

inline bool fmpz_divisible(FmpzConstRef value,
                           FmpzConstRef divisor) noexcept {
    return ::fmpz_divisible(value.raw(), divisor.raw()) != 0;
}

inline slong fmpz_sgn(FmpzConstRef value) noexcept {
    return ::fmpz_sgn(value.raw());
}

inline slong fmpz_sgn(FmpzRef value) noexcept {
    return ::fmpz_sgn(value.raw());
}

inline ulong fmpz_fdiv_ui(FmpzConstRef value, ulong modulus) noexcept {
    return ::fmpz_fdiv_ui(value.raw(), modulus);
}

inline void fmpz_fdiv_q_2exp(FmpzRef out,
                             FmpzConstRef value,
                             ulong exponent) noexcept {
    ::fmpz_fdiv_q_2exp(out.raw(), value.raw(), exponent);
}

inline void fmpz_fdiv_q(FmpzRef out,
                        FmpzConstRef numerator,
                        FmpzConstRef denominator) noexcept {
    ::fmpz_fdiv_q(out.raw(), numerator.raw(), denominator.raw());
}

inline void fmpz_tdiv_qr(FmpzRef quotient,
                         FmpzRef remainder,
                         FmpzConstRef numerator,
                         FmpzConstRef denominator) noexcept {
    ::fmpz_tdiv_qr(quotient.raw(), remainder.raw(), numerator.raw(),
                   denominator.raw());
}

inline bool fmpz_tstbit(FmpzConstRef value, ulong bit) noexcept {
    return ::fmpz_tstbit(value.raw(), bit) != 0;
}

inline void fmpz_sqrtrem(FmpzRef root,
                         FmpzRef remainder,
                         FmpzConstRef input) noexcept {
    ::fmpz_sqrtrem(root.raw(), remainder.raw(), input.raw());
}

inline void fmpz_sqrt(FmpzRef out, FmpzConstRef input) noexcept {
    ::fmpz_sqrt(out.raw(), input.raw());
}

inline bool fmpz_equal_si(FmpzConstRef value, slong expected) noexcept {
    return ::fmpz_equal_si(value.raw(), expected) != 0;
}

inline bool fmpz_equal_si(const Fmpz& value, slong expected) noexcept {
    return ::fmpz_equal_si(value.raw(), expected) != 0;
}

inline bool fmpz_equal_si(FmpzRef value, slong expected) noexcept {
    return ::fmpz_equal_si(value.raw(), expected) != 0;
}

inline bool fmpz_equal(FmpzConstRef left, FmpzConstRef right) noexcept {
    return ::fmpz_equal(left.raw(), right.raw()) != 0;
}

inline bool fmpz_equal(const Fmpz& left, const Fmpz& right) noexcept {
    return ::fmpz_equal(left.raw(), right.raw()) != 0;
}

inline void fmpz_abs(FmpzRef out, FmpzConstRef in) noexcept {
    ::fmpz_abs(out.raw(), in.raw());
}

inline void fmpz_gcd(FmpzRef out,
                     FmpzConstRef left,
                     FmpzConstRef right) noexcept {
    ::fmpz_gcd(out.raw(), left.raw(), right.raw());
}

inline void fmpz_divexact(FmpzRef out,
                          FmpzConstRef left,
                          FmpzConstRef right) noexcept {
    ::fmpz_divexact(out.raw(), left.raw(), right.raw());
}

inline bool fmpz_invmod(FmpzRef out,
                        FmpzConstRef value,
                        FmpzConstRef modulus) noexcept {
    return ::fmpz_invmod(out.raw(), value.raw(), modulus.raw()) != 0;
}

inline int fmpz_kronecker(FmpzConstRef left, FmpzConstRef right) noexcept {
    return ::fmpz_kronecker(left.raw(), right.raw());
}

inline void fmpz_mul_2exp(FmpzRef out, FmpzConstRef in, ulong exponent) noexcept {
    ::fmpz_mul_2exp(out.raw(), in.raw(), exponent);
}

}  // namespace silex::flint
