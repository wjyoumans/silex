#pragma once

#include <flint/flint.h>

#include <silex/number_field.hpp>

namespace silex {

class Signature {
public:
    Signature() noexcept = default;
    Signature(slong real_embeddings, slong complex_pairs) noexcept
        : r1_(real_embeddings),
          r2_(complex_pairs) {
    }

    slong r1() const noexcept { return r1_; }
    slong r2() const noexcept { return r2_; }
    slong degree() const noexcept { return r1_ + 2 * r2_; }

    void set(slong real_embeddings, slong complex_pairs) noexcept;
    bool compute(const NumberField& field) noexcept;

private:
    slong r1_ = 0;
    slong r2_ = 0;
};

bool signature(Signature& out, const NumberField& field) noexcept;
bool is_totally_real(const NumberField& field) noexcept;
bool is_totally_complex(const NumberField& field) noexcept;

}  // namespace silex
