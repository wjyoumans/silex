#pragma once

#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq.hpp>

#include <flint/arb.h>

namespace silex::flint {

class ArbRef;
class ArbConstRef;

class Arb {
public:
    Arb() noexcept {
        arb_init(value_);
    }

    ~Arb() noexcept {
        arb_clear(value_);
    }

    Arb(const Arb&) = delete;
    Arb& operator=(const Arb&) = delete;

    Arb(Arb&& other) noexcept {
        arb_init(value_);
        arb_swap(value_, other.value_);
    }

    Arb& operator=(Arb&& other) noexcept {
        if (this != &other) {
            arb_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(Arb& other) noexcept {
        arb_swap(value_, other.value_);
    }

    arb_t& raw() noexcept {
        return value_;
    }

    const arb_t& raw() const noexcept {
        return value_;
    }

private:
    arb_t value_;
};

class ArbConstRef {
public:
    ArbConstRef(const arb_struct* value) noexcept
        : value_(value) {
    }

    ArbConstRef(const Arb& value) noexcept
        : value_(value.raw()) {
    }

    const arb_struct* raw() const noexcept {
        return value_;
    }

private:
    const arb_struct* value_;
};

class ArbRef {
public:
    ArbRef(arb_struct* value) noexcept
        : value_(value) {
    }

    ArbRef(Arb& value) noexcept
        : value_(value.raw()) {
    }

    arb_struct* raw() noexcept {
        return value_;
    }

    const arb_struct* raw() const noexcept {
        return value_;
    }

private:
    arb_struct* value_;
};

inline void swap(Arb& left, Arb& right) noexcept {
    left.swap(right);
}

inline void arb_neg(Arb& out, const Arb& in) noexcept {
    ::arb_neg(out.raw(), in.raw());
}

inline bool arb_is_negative(const Arb& value) noexcept {
    return ::arb_is_negative(value.raw()) != 0;
}

inline bool arb_is_negative(const arb_struct* value) noexcept {
    return ::arb_is_negative(value) != 0;
}

inline bool arb_is_nonnegative(const arb_struct* value) noexcept {
    return ::arb_is_nonnegative(value) != 0;
}

inline bool arb_is_positive(const Arb& value) noexcept {
    return ::arb_is_positive(value.raw()) != 0;
}

inline bool arb_is_positive(ArbConstRef value) noexcept {
    return ::arb_is_positive(value.raw()) != 0;
}

inline bool arb_is_positive(const arb_struct* value) noexcept {
    return ::arb_is_positive(value) != 0;
}

inline void arb_zero(Arb& out) noexcept {
    ::arb_zero(out.raw());
}

inline void arb_zero(ArbRef out) noexcept {
    ::arb_zero(out.raw());
}

inline void arb_one(Arb& out) noexcept {
    ::arb_one(out.raw());
}

inline void arb_set_si(ArbRef out, slong value) noexcept {
    ::arb_set_si(out.raw(), value);
}

inline void arb_set_si(Arb& out, slong value) noexcept {
    ::arb_set_si(out.raw(), value);
}

inline void arb_set_si(arb_struct* out, slong value) noexcept {
    ::arb_set_si(out, value);
}

inline void arb_set_ui(ArbRef out, ulong value) noexcept {
    ::arb_set_ui(out.raw(), value);
}

inline void arb_set_ui(Arb& out, ulong value) noexcept {
    ::arb_set_ui(out.raw(), value);
}

inline void arb_set_d(Arb& out, double value) noexcept {
    ::arb_set_d(out.raw(), value);
}

inline void arb_set_d(ArbRef out, double value) noexcept {
    ::arb_set_d(out.raw(), value);
}

inline void arb_set(Arb& out, const arb_struct* in) noexcept {
    ::arb_set(out.raw(), in);
}

inline void arb_set(ArbRef out, ArbConstRef in) noexcept {
    ::arb_set(out.raw(), in.raw());
}

inline void arb_set(ArbRef out, ArbRef in) noexcept {
    ::arb_set(out.raw(), in.raw());
}

inline void arb_set(Arb& out, const Arb& in) noexcept {
    ::arb_set(out.raw(), in.raw());
}

inline void arb_set(Arb& out, ArbRef in) noexcept {
    ::arb_set(out.raw(), in.raw());
}

inline void arb_set(arb_struct* out, const Arb& in) noexcept {
    ::arb_set(out, in.raw());
}

inline void arb_set_fmpz(Arb& out, const fmpz* in) noexcept {
    ::arb_set_fmpz(out.raw(), in);
}

inline void arb_set_fmpz(Arb& out, FmpzConstRef in) noexcept {
    ::arb_set_fmpz(out.raw(), in.raw());
}

inline void arb_set_fmpz(Arb& out, const Fmpz& in) noexcept {
    ::arb_set_fmpz(out.raw(), in.raw());
}

inline void arb_set_fmpq(Arb& out, const Fmpq& in, slong prec) noexcept {
    ::arb_set_fmpq(out.raw(), in.raw(), prec);
}

inline void arb_abs(Arb& out, const Arb& in) noexcept {
    ::arb_abs(out.raw(), in.raw());
}

inline void arb_abs(ArbRef out, ArbConstRef in) noexcept {
    ::arb_abs(out.raw(), in.raw());
}

inline void arb_add(Arb& out, const Arb& left, const Arb& right, slong prec) noexcept {
    ::arb_add(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_add(Arb& out,
        const arb_struct* left,
        const Arb& right,
        slong prec) noexcept {
    ::arb_add(out.raw(), left, right.raw(), prec);
}

inline void arb_add(Arb& out,
        const Arb& left,
        const arb_struct* right,
        slong prec) noexcept {
    ::arb_add(out.raw(), left.raw(), right, prec);
}

inline void arb_add_ui(Arb& out, const Arb& left, ulong right, slong prec) noexcept {
    ::arb_add_ui(out.raw(), left.raw(), right, prec);
}

inline void arb_sub(Arb& out,
        const Arb& left,
        const Arb& right,
        slong prec) noexcept {
    ::arb_sub(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_sub(Arb& out,
        const arb_struct* left,
        const arb_struct* right,
        slong prec) noexcept {
    ::arb_sub(out.raw(), left, right, prec);
}

inline void arb_sub(Arb& out,
        const arb_struct* left,
        const Arb& right,
        slong prec) noexcept {
    ::arb_sub(out.raw(), left, right.raw(), prec);
}

inline void arb_mul(Arb& out, const Arb& left, const Arb& right, slong prec) noexcept {
    ::arb_mul(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_mul(Arb& out, const Arb& left, ArbConstRef right, slong prec) noexcept {
    ::arb_mul(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_mul(
        Arb& out, const arb_struct* left, const arb_struct* right, slong prec) noexcept {
    ::arb_mul(out.raw(), left, right, prec);
}

inline void arb_mul(
        Arb& out, ArbConstRef left, ArbConstRef right, slong prec) noexcept {
    ::arb_mul(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_mul(
        Arb& out, const Arb& left, const arb_struct* right, slong prec) noexcept {
    ::arb_mul(out.raw(), left.raw(), right, prec);
}

inline void arb_mul_ui(Arb& out, const Arb& left, ulong right, slong prec) noexcept {
    ::arb_mul_ui(out.raw(), left.raw(), right, prec);
}

inline void arb_mul_fmpz(Arb& out,
        const arb_struct* left,
        const fmpz* right,
        slong prec) noexcept {
    ::arb_mul_fmpz(out.raw(), left, right, prec);
}

inline void arb_mul_fmpz(Arb& out,
        const Arb& left,
        const fmpz* right,
        slong prec) noexcept {
    ::arb_mul_fmpz(out.raw(), left.raw(), right, prec);
}

inline void arb_mul_fmpz(Arb& out,
        const Arb& left,
        const Fmpz& right,
        slong prec) noexcept {
    ::arb_mul_fmpz(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_mul_fmpz(Arb& out,
        const Arb& left,
        FmpzConstRef right,
        slong prec) noexcept {
    ::arb_mul_fmpz(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_addmul_fmpz(ArbRef out,
        ArbConstRef left,
        FmpzConstRef right,
        slong prec) noexcept {
    ::arb_addmul_fmpz(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_div(Arb& out, const Arb& left, const Arb& right, slong prec) noexcept {
    ::arb_div(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_div(Arb& out, const Arb& left, ArbConstRef right, slong prec) noexcept {
    ::arb_div(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_div_ui(Arb& out, const Arb& left, ulong right, slong prec) noexcept {
    ::arb_div_ui(out.raw(), left.raw(), right, prec);
}

inline void arb_div_fmpz(Arb& out,
        const Arb& left,
        FmpzConstRef right,
        slong prec) noexcept {
    ::arb_div_fmpz(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_inv(Arb& out, const Arb& in, slong prec) noexcept {
    ::arb_inv(out.raw(), in.raw(), prec);
}

inline void arb_sqrt(Arb& out, const Arb& in, slong prec) noexcept {
    ::arb_sqrt(out.raw(), in.raw(), prec);
}

inline void arb_sqrt_ui(Arb& out, ulong value, slong prec) noexcept {
    ::arb_sqrt_ui(out.raw(), value, prec);
}

inline void arb_root_ui(Arb& out, const Arb& in, ulong exponent, slong prec) noexcept {
    ::arb_root_ui(out.raw(), in.raw(), exponent, prec);
}

inline void arb_sqr(Arb& out, const Arb& in, slong prec) noexcept {
    ::arb_sqr(out.raw(), in.raw(), prec);
}

inline void arb_log(Arb& out, const Arb& in, slong prec) noexcept {
    ::arb_log(out.raw(), in.raw(), prec);
}

inline void arb_log(Arb& out, ArbConstRef in, slong prec) noexcept {
    ::arb_log(out.raw(), in.raw(), prec);
}

inline void arb_log_ui(Arb& out, ulong value, slong prec) noexcept {
    ::arb_log_ui(out.raw(), value, prec);
}

inline void arb_log_fmpz(Arb& out, FmpzConstRef value, slong prec) noexcept {
    ::arb_log_fmpz(out.raw(), value.raw(), prec);
}

inline void arb_exp(Arb& out, const Arb& in, slong prec) noexcept {
    ::arb_exp(out.raw(), in.raw(), prec);
}

inline void arb_max(Arb& out, const Arb& left, const Arb& right, slong prec) noexcept {
    ::arb_max(out.raw(), left.raw(), right.raw(), prec);
}

inline void arb_mul_2exp_si(Arb& out, const arb_struct* in, slong exponent) noexcept {
    ::arb_mul_2exp_si(out.raw(), in, exponent);
}

inline void arb_mul_2exp_si(Arb& out, const Arb& in, slong exponent) noexcept {
    ::arb_mul_2exp_si(out.raw(), in.raw(), exponent);
}

inline void arb_const_pi(Arb& out, slong prec) noexcept {
    ::arb_const_pi(out.raw(), prec);
}

inline bool arb_lt(const Arb& left, const Arb& right) noexcept {
    return ::arb_lt(left.raw(), right.raw()) != 0;
}

inline bool arb_gt(const Arb& left, const Arb& right) noexcept {
    return ::arb_gt(left.raw(), right.raw()) != 0;
}

inline void arb_get_rad_arb(Arb& out, const Arb& in) noexcept {
    ::arb_get_rad_arb(out.raw(), in.raw());
}

inline void arb_add_error(Arb& out, const Arb& error) noexcept {
    ::arb_add_error(out.raw(), error.raw());
}

inline bool arb_contains_si(ArbConstRef value, slong expected) noexcept {
    return ::arb_contains_si(value.raw(), expected) != 0;
}

inline bool arb_contains_si(ArbRef value, slong expected) noexcept {
    return ::arb_contains_si(value.raw(), expected) != 0;
}

inline bool arb_contains_si(const Arb& value, slong expected) noexcept {
    return ::arb_contains_si(value.raw(), expected) != 0;
}

inline bool arb_contains_si(const arb_struct* value, slong expected) noexcept {
    return ::arb_contains_si(value, expected) != 0;
}

inline bool arb_contains(ArbConstRef left, ArbConstRef right) noexcept {
    return ::arb_contains(left.raw(), right.raw()) != 0;
}

inline bool arb_contains(const Arb& left, const Arb& right) noexcept {
    return ::arb_contains(left.raw(), right.raw()) != 0;
}

inline bool arb_contains_fmpq(ArbConstRef value, FmpqConstRef expected) noexcept {
    return ::arb_contains_fmpq(value.raw(), expected.raw()) != 0;
}

inline bool arb_contains_fmpq(ArbConstRef value, const Fmpq& expected) noexcept {
    return ::arb_contains_fmpq(value.raw(), expected.raw()) != 0;
}

inline bool arb_contains_fmpq(ArbRef value, FmpqConstRef expected) noexcept {
    return ::arb_contains_fmpq(value.raw(), expected.raw()) != 0;
}

inline bool arb_contains_fmpq(ArbRef value, const Fmpq& expected) noexcept {
    return ::arb_contains_fmpq(value.raw(), expected.raw()) != 0;
}

inline bool arb_contains_fmpq(const arb_struct* value, const Fmpq& expected) noexcept {
    return ::arb_contains_fmpq(value, expected.raw()) != 0;
}

inline bool arb_contains_zero(const arb_struct* value) noexcept {
    return ::arb_contains_zero(value) != 0;
}

inline bool arb_contains_zero(const Arb& value) noexcept {
    return ::arb_contains_zero(value.raw()) != 0;
}

inline bool arb_contains_zero(ArbConstRef value) noexcept {
    return ::arb_contains_zero(value.raw()) != 0;
}

inline bool arb_is_finite(const Arb& value) noexcept {
    return ::arb_is_finite(value.raw()) != 0;
}

inline bool arb_is_finite(ArbConstRef value) noexcept {
    return ::arb_is_finite(value.raw()) != 0;
}

inline bool arb_is_zero(const Arb& value) noexcept {
    return ::arb_is_zero(value.raw()) != 0;
}

inline bool arb_is_one(const Arb& value) noexcept {
    return ::arb_is_one(value.raw()) != 0;
}

inline bool arb_is_one(ArbConstRef value) noexcept {
    return ::arb_is_one(value.raw()) != 0;
}

inline bool arb_overlaps(ArbConstRef left, ArbConstRef right) noexcept {
    return ::arb_overlaps(left.raw(), right.raw()) != 0;
}

inline bool arb_overlaps(const arb_struct* left, const Arb& right) noexcept {
    return ::arb_overlaps(left, right.raw()) != 0;
}

inline bool arb_overlaps(const Arb& left, const Arb& right) noexcept {
    return ::arb_overlaps(left.raw(), right.raw()) != 0;
}

inline void arb_get_lbound_arf(Arf& out, const Arb& in, slong prec) noexcept {
    ::arb_get_lbound_arf(out.raw(), in.raw(), prec);
}

inline void arb_get_ubound_arf(Arf& out, const Arb& in, slong prec) noexcept {
    ::arb_get_ubound_arf(out.raw(), in.raw(), prec);
}

inline void arb_get_ubound_arf(Arf& out, const arb_struct* in, slong prec) noexcept {
    ::arb_get_ubound_arf(out.raw(), in, prec);
}

inline void arb_get_abs_ubound_arf(Arf& out, ArbConstRef in, slong prec) noexcept {
    ::arb_get_abs_ubound_arf(out.raw(), in.raw(), prec);
}

inline void arb_get_abs_ubound_arf(Arf& out, const Arb& in, slong prec) noexcept {
    ::arb_get_abs_ubound_arf(out.raw(), in.raw(), prec);
}

inline void arb_get_abs_ubound_arf(Arf& out, const arb_struct* in, slong prec) noexcept {
    ::arb_get_abs_ubound_arf(out.raw(), in, prec);
}

}  // namespace silex::flint
