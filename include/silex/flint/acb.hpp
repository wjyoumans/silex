#pragma once

#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>

#include <flint/acb.h>
#include <flint/acb_dirichlet.h>

namespace silex::flint {

class AcbRef;
class AcbConstRef;

class Acb {
public:
    Acb() noexcept { acb_init(value_); }
    ~Acb() noexcept { acb_clear(value_); }

    Acb(const Acb&) = delete;
    Acb& operator=(const Acb&) = delete;

    Acb(Acb&& other) noexcept {
        acb_init(value_);
        acb_swap(value_, other.value_);
    }

    Acb& operator=(Acb&& other) noexcept {
        if (this != &other) {
            acb_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(Acb& other) noexcept { acb_swap(value_, other.value_); }
    acb_t& raw() noexcept { return value_; }
    const acb_t& raw() const noexcept { return value_; }

private:
    acb_t value_;
};

class AcbConstRef {
public:
    explicit AcbConstRef(const acb_struct* value) noexcept : value_(value) {}
    explicit AcbConstRef(const Acb& value) noexcept : value_(value.raw()) {}
    const acb_struct* raw() const noexcept { return value_; }

private:
    const acb_struct* value_;
};

class AcbRef {
public:
    explicit AcbRef(acb_struct* value) noexcept : value_(value) {}
    explicit AcbRef(Acb& value) noexcept : value_(value.raw()) {}
    acb_struct* raw() noexcept { return value_; }
    const acb_struct* raw() const noexcept { return value_; }

private:
    acb_struct* value_;
};

inline void swap(Acb& left, Acb& right) noexcept { left.swap(right); }

inline void acb_set_si(AcbRef out, slong value) noexcept {
    ::acb_set_si(out.raw(), value);
}

inline void acb_set_si(Acb& out, slong value) noexcept {
    ::acb_set_si(out.raw(), value);
}

inline void acb_set_si(acb_struct* out, slong value) noexcept {
    ::acb_set_si(out, value);
}

inline void acb_add(
        Acb& out, const acb_struct* left, const acb_struct* right, slong prec) noexcept {
    ::acb_add(out.raw(), left, right, prec);
}

inline void acb_add(Acb& out, const Acb& left, const acb_struct* right, slong prec) noexcept {
    ::acb_add(out.raw(), left.raw(), right, prec);
}

inline void acb_mul(
        Acb& out, const acb_struct* left, const acb_struct* right, slong prec) noexcept {
    ::acb_mul(out.raw(), left, right, prec);
}

inline void acb_mul(Acb& out, AcbConstRef left, AcbConstRef right, slong prec) noexcept {
    ::acb_mul(out.raw(), left.raw(), right.raw(), prec);
}

inline void acb_mul(Acb& out, const Acb& left, const acb_struct* right, slong prec) noexcept {
    ::acb_mul(out.raw(), left.raw(), right, prec);
}

inline void acb_pow_ui(Acb& out, const acb_struct* in, ulong exponent, slong prec) noexcept {
    ::acb_pow_ui(out.raw(), in, exponent, prec);
}

inline void acb_pow_ui(Acb& out, AcbConstRef in, ulong exponent, slong prec) noexcept {
    ::acb_pow_ui(out.raw(), in.raw(), exponent, prec);
}

inline void acb_abs(Arb& out, AcbConstRef in, slong prec) noexcept {
    ::acb_abs(out.raw(), in.raw(), prec);
}

inline void acb_abs(Arb& out, const Acb& in, slong prec) noexcept {
    ::acb_abs(out.raw(), in.raw(), prec);
}

inline void acb_abs(Arb& out, const acb_struct* in, slong prec) noexcept {
    ::acb_abs(out.raw(), in, prec);
}

inline void acb_sub_si(Acb& out, const Acb& in, slong value, slong prec) noexcept {
    ::acb_sub_si(out.raw(), in.raw(), value, prec);
}

inline bool acb_equal_si(AcbConstRef value, slong expected) noexcept {
    return ::acb_equal_si(value.raw(), expected) != 0;
}

inline bool acb_equal_si(const Acb& value, slong expected) noexcept {
    return ::acb_equal_si(value.raw(), expected) != 0;
}

inline bool acb_equal_si(const acb_struct* value, slong expected) noexcept {
    return ::acb_equal_si(value, expected) != 0;
}

inline bool acb_contains_zero(AcbConstRef value) noexcept {
    return ::acb_contains_zero(value.raw()) != 0;
}

inline bool acb_contains_zero(const Acb& value) noexcept {
    return ::acb_contains_zero(value.raw()) != 0;
}

inline bool acb_is_finite(const Acb& value) noexcept {
    return ::acb_is_finite(value.raw()) != 0;
}

inline bool acb_contains_fmpq(AcbConstRef value, FmpqConstRef expected) noexcept {
    return ::acb_contains_fmpq(value.raw(), expected.raw()) != 0;
}

inline bool acb_contains_fmpq(const Acb& value, const Fmpq& expected) noexcept {
    return ::acb_contains_fmpq(value.raw(), expected.raw()) != 0;
}

inline const arb_struct* acb_realref_ptr(AcbConstRef value) noexcept {
    return acb_realref(value.raw());
}

inline ArbConstRef acb_real_part(AcbConstRef value) noexcept {
    return ArbConstRef(&value.raw()->real);
}

inline ArbConstRef acb_real_part(const Acb& value) noexcept {
    return ArbConstRef(&value.raw()->real);
}

inline const arb_struct* acb_realref_ptr(const Acb& value) noexcept {
    return acb_realref(value.raw());
}

inline const arb_struct* acb_realref_ptr(const acb_struct* value) noexcept {
    return acb_realref(value);
}

inline const arb_struct* acb_imagref_ptr(AcbConstRef value) noexcept {
    return acb_imagref(value.raw());
}

inline ArbConstRef acb_imag_part(AcbConstRef value) noexcept {
    return ArbConstRef(&value.raw()->imag);
}

inline ArbConstRef acb_imag_part(const Acb& value) noexcept {
    return ArbConstRef(&value.raw()->imag);
}

inline const arb_struct* acb_imagref_ptr(const Acb& value) noexcept {
    return acb_imagref(value.raw());
}

inline const arb_struct* acb_imagref_ptr(const acb_struct* value) noexcept {
    return acb_imagref(value);
}

inline void acb_dirichlet_l_fmpq(Acb& out,
                                 FmpqConstRef input,
                                 const dirichlet_group_struct* group,
                                 const dirichlet_char_struct* character,
                                 slong precision) noexcept {
    ::acb_dirichlet_l_fmpq(out.raw(), input.raw(), group, character,
                           precision);
}

}  // namespace silex::flint
