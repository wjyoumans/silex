#pragma once

#include <silex/flint/fmpz.hpp>

#include <flint/qfb.h>

namespace silex::flint {

class QfbRef;
class QfbConstRef;

class Qfb {
public:
    Qfb() noexcept {
        ::qfb_init(value_);
    }

    ~Qfb() noexcept {
        ::qfb_clear(value_);
    }

    Qfb(const Qfb&) = delete;
    Qfb& operator=(const Qfb&) = delete;

    Qfb(Qfb&& other) noexcept
        : Qfb() {
        swap(other);
    }

    Qfb& operator=(Qfb&& other) noexcept {
        if (this != &other) {
            swap(other);
        }
        return *this;
    }

    void swap(Qfb& other) noexcept {
        ::fmpz_swap(value_->a, other.value_->a);
        ::fmpz_swap(value_->b, other.value_->b);
        ::fmpz_swap(value_->c, other.value_->c);
    }

    ::qfb* raw() noexcept {
        return value_;
    }

    ::qfb* raw() const noexcept {
        return value_;
    }

private:
    // FLINT's qfb API is not uniformly const-correct.  Keeping the owned
    // storage mutable confines that legacy signature mismatch to this RAII
    // bridge while QfbConstRef still exposes only read-oriented operations.
    mutable qfb_t value_;
};

class QfbConstRef {
public:
    explicit QfbConstRef(::qfb* value) noexcept
        : value_(value) {
    }

    explicit QfbConstRef(const Qfb& value) noexcept
        : value_(value.raw()) {
    }

    ::qfb* raw() const noexcept {
        return value_;
    }

private:
    ::qfb* value_;
};

class QfbRef {
public:
    explicit QfbRef(::qfb* value) noexcept
        : value_(value) {
    }

    explicit QfbRef(Qfb& value) noexcept
        : value_(value.raw()) {
    }

    ::qfb* raw() noexcept {
        return value_;
    }

    const ::qfb* raw() const noexcept {
        return value_;
    }

private:
    ::qfb* value_;
};

inline void swap(Qfb& left, Qfb& right) noexcept {
    left.swap(right);
}

inline void qfb_set(QfbRef out, QfbConstRef input) noexcept {
    ::qfb_set(out.raw(), input.raw());
}

inline void qfb_prime_form(QfbRef out,
                           FmpzConstRef discriminant,
                           FmpzConstRef prime) noexcept {
    Fmpz discriminant_copy;
    Fmpz prime_copy;
    fmpz_set(FmpzRef(discriminant_copy), discriminant);
    fmpz_set(FmpzRef(prime_copy), prime);
    ::qfb_prime_form(out.raw(), discriminant_copy.raw(), prime_copy.raw());
}

inline void qfb_principal_form(QfbRef out,
                               FmpzConstRef discriminant) noexcept {
    Fmpz discriminant_copy;
    fmpz_set(FmpzRef(discriminant_copy), discriminant);
    ::qfb_principal_form(out.raw(), discriminant_copy.raw());
}

inline void qfb_reduce(QfbRef out,
                       QfbConstRef input,
                       FmpzConstRef discriminant) noexcept {
    Fmpz discriminant_copy;
    fmpz_set(FmpzRef(discriminant_copy), discriminant);
    ::qfb_reduce(out.raw(), input.raw(), discriminant_copy.raw());
}

inline void qfb_nucomp(QfbRef out,
                       QfbConstRef left,
                       QfbConstRef right,
                       FmpzConstRef discriminant,
                       FmpzConstRef fourth_root) noexcept {
    Fmpz discriminant_copy;
    Fmpz fourth_root_copy;
    fmpz_set(FmpzRef(discriminant_copy), discriminant);
    fmpz_set(FmpzRef(fourth_root_copy), fourth_root);
    ::qfb_nucomp(out.raw(), left.raw(), right.raw(), discriminant_copy.raw(),
                 fourth_root_copy.raw());
}

inline void qfb_pow_ui(QfbRef out,
                       QfbConstRef input,
                       FmpzConstRef discriminant,
                       ulong exponent) noexcept {
    Fmpz discriminant_copy;
    fmpz_set(FmpzRef(discriminant_copy), discriminant);
    ::qfb_pow_ui(out.raw(), input.raw(), discriminant_copy.raw(), exponent);
}

inline void qfb_inverse(QfbRef out, QfbConstRef input) noexcept {
    ::qfb_inverse(out.raw(), input.raw());
}

inline bool qfb_equal(QfbConstRef left, QfbConstRef right) noexcept {
    return ::qfb_equal(left.raw(), right.raw()) != 0;
}

inline bool qfb_is_principal_form(QfbConstRef form,
                                  FmpzConstRef discriminant) noexcept {
    Fmpz discriminant_copy;
    fmpz_set(FmpzRef(discriminant_copy), discriminant);
    return ::qfb_is_principal_form(form.raw(), discriminant_copy.raw()) != 0;
}

inline FmpzConstRef qfb_a(QfbConstRef form) noexcept {
    return FmpzConstRef(form.raw()->a);
}

inline FmpzConstRef qfb_b(QfbConstRef form) noexcept {
    return FmpzConstRef(form.raw()->b);
}

inline FmpzConstRef qfb_c(QfbConstRef form) noexcept {
    return FmpzConstRef(form.raw()->c);
}

}  // namespace silex::flint
