#pragma once

#include <flint/dirichlet.h>

namespace silex::flint {

class DirichletGroupRef;
class DirichletGroupConstRef;
class DirichletCharRef;
class DirichletCharConstRef;

class DirichletGroup {
public:
    explicit DirichletGroup(ulong modulus) noexcept {
        initialized_ = dirichlet_group_init(value_, modulus) != 0;
    }

    ~DirichletGroup() noexcept {
        if (initialized_) {
            dirichlet_group_clear(value_);
        }
    }

    DirichletGroup(const DirichletGroup&) = delete;
    DirichletGroup& operator=(const DirichletGroup&) = delete;
    DirichletGroup(DirichletGroup&&) = delete;
    DirichletGroup& operator=(DirichletGroup&&) = delete;

    dirichlet_group_t& raw() noexcept { return value_; }
    const dirichlet_group_t& raw() const noexcept { return value_; }
    bool is_initialized() const noexcept { return initialized_; }
    ulong modulus() const noexcept { return value_->q; }
    ulong character_count() const noexcept { return value_->phi_q; }
    ulong exponent() const noexcept { return value_->expo; }

private:
    dirichlet_group_t value_;
    bool initialized_ = false;
};

class DirichletGroupConstRef {
public:
    explicit DirichletGroupConstRef(const dirichlet_group_t& value) noexcept : value_(value) {}
    explicit DirichletGroupConstRef(const DirichletGroup& value) noexcept : value_(value.raw()) {}
    const dirichlet_group_struct* raw() const noexcept { return value_; }

private:
    const dirichlet_group_struct* value_;
};

class DirichletGroupRef {
public:
    explicit DirichletGroupRef(dirichlet_group_t& value) noexcept : value_(value) {}
    explicit DirichletGroupRef(DirichletGroup& value) noexcept : value_(value.raw()) {}
    dirichlet_group_struct* raw() noexcept { return value_; }
    const dirichlet_group_struct* raw() const noexcept { return value_; }

private:
    dirichlet_group_struct* value_;
};

class DirichletChar {
public:
    explicit DirichletChar(const dirichlet_group_t& group) noexcept
        : group_(group) {
        dirichlet_char_init(value_, group_);
    }

    explicit DirichletChar(const DirichletGroup& group) noexcept
        : DirichletChar(group.raw()) {
    }

    ~DirichletChar() noexcept { dirichlet_char_clear(value_); }

    DirichletChar(const DirichletChar&) = delete;
    DirichletChar& operator=(const DirichletChar&) = delete;
    DirichletChar(DirichletChar&&) = delete;
    DirichletChar& operator=(DirichletChar&&) = delete;

    dirichlet_char_t& raw() noexcept { return value_; }
    const dirichlet_char_t& raw() const noexcept { return value_; }
    const dirichlet_group_struct* group() const noexcept { return group_; }

private:
    const dirichlet_group_struct* group_;
    dirichlet_char_t value_;
};

class DirichletCharConstRef {
public:
    explicit DirichletCharConstRef(const dirichlet_char_t& value) noexcept : value_(value) {}
    explicit DirichletCharConstRef(const DirichletChar& value) noexcept : value_(value.raw()) {}
    const dirichlet_char_struct* raw() const noexcept { return value_; }

private:
    const dirichlet_char_struct* value_;
};

class DirichletCharRef {
public:
    explicit DirichletCharRef(dirichlet_char_t& value) noexcept : value_(value) {}
    explicit DirichletCharRef(DirichletChar& value) noexcept : value_(value.raw()) {}
    dirichlet_char_struct* raw() noexcept { return value_; }
    const dirichlet_char_struct* raw() const noexcept { return value_; }

private:
    dirichlet_char_struct* value_;
};

inline void dirichlet_char_index(DirichletChar& out,
                                 const DirichletGroup& group,
                                 ulong index) noexcept {
    ::dirichlet_char_index(out.raw(), group.raw(), index);
}

inline void dirichlet_char_set(DirichletChar& out,
                               const DirichletGroup& group,
                               const DirichletChar& in) noexcept {
    ::dirichlet_char_set(out.raw(), group.raw(), in.raw());
}

inline bool dirichlet_char_is_real(const DirichletGroup& group,
                                   const DirichletChar& character) noexcept {
    return ::dirichlet_char_is_real(group.raw(), character.raw()) != 0;
}

inline bool dirichlet_char_is_primitive(
        const DirichletGroup& group,
        const DirichletChar& character) noexcept {
    return ::dirichlet_char_is_primitive(group.raw(), character.raw()) != 0;
}

inline ulong dirichlet_chi(const DirichletGroup& group,
                           const DirichletChar& character,
                           ulong n) noexcept {
    return ::dirichlet_chi(group.raw(), character.raw(), n);
}

}  // namespace silex::flint
