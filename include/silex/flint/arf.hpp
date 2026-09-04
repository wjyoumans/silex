#pragma once

#include <silex/flint/fmpz.hpp>

#include <flint/arf.h>

namespace silex::flint {

class ArfRef;
class ArfConstRef;

class Arf {
public:
    Arf() noexcept { arf_init(value_); }
    ~Arf() noexcept { arf_clear(value_); }

    Arf(const Arf&) = delete;
    Arf& operator=(const Arf&) = delete;

    Arf(Arf&& other) noexcept {
        arf_init(value_);
        arf_swap(value_, other.value_);
    }

    Arf& operator=(Arf&& other) noexcept {
        if (this != &other) {
            arf_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(Arf& other) noexcept { arf_swap(value_, other.value_); }
    arf_t& raw() noexcept { return value_; }
    const arf_t& raw() const noexcept { return value_; }

private:
    arf_t value_;
};

class ArfConstRef {
public:
    explicit ArfConstRef(const arf_t& value) noexcept : value_(value) {}
    explicit ArfConstRef(const Arf& value) noexcept : value_(value.raw()) {}
    const arf_struct* raw() const noexcept { return value_; }

private:
    const arf_struct* value_;
};

class ArfRef {
public:
    explicit ArfRef(arf_t& value) noexcept : value_(value) {}
    explicit ArfRef(Arf& value) noexcept : value_(value.raw()) {}
    arf_struct* raw() noexcept { return value_; }
    const arf_struct* raw() const noexcept { return value_; }

private:
    arf_struct* value_;
};

inline void swap(Arf& left, Arf& right) noexcept { left.swap(right); }

inline void arf_set(Arf& out, const arf_struct* in) noexcept {
    ::arf_set(out.raw(), in);
}

inline bool arf_is_nan(const Arf& value) noexcept {
    return ::arf_is_nan(value.raw()) != 0;
}

inline bool arf_is_inf(const Arf& value) noexcept {
    return ::arf_is_inf(value.raw()) != 0;
}

inline bool arf_is_finite(const Arf& value) noexcept {
    return ::arf_is_finite(value.raw()) != 0;
}

inline int arf_cmp(const Arf& left, const Arf& right) noexcept {
    return ::arf_cmp(left.raw(), right.raw());
}

inline int arf_cmp_ui(const Arf& left, ulong right) noexcept {
    return ::arf_cmp_ui(left.raw(), right);
}

inline double arf_get_d(const Arf& value, arf_rnd_t rnd) noexcept {
    return ::arf_get_d(value.raw(), rnd);
}

inline void arf_get_fmpz(Fmpz& out, const Arf& value, arf_rnd_t rnd) noexcept {
    ::arf_get_fmpz(out.raw(), value.raw(), rnd);
}

inline void arf_get_fmpz(FmpzRef out, const Arf& value, arf_rnd_t rnd) noexcept {
    ::arf_get_fmpz(out.raw(), value.raw(), rnd);
}

}  // namespace silex::flint
