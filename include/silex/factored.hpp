#pragma once

#include <flint/flint.h>

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace silex {

template <class Factor>
struct FactorPower {
    Factor factor;
    slong exponent = 0;
};

template <class Factor>
using FactorSpan = std::span<const FactorPower<Factor>>;

template <class Factor>
class Factored {
public:
    using factor_type = Factor;

    FactorSpan<Factor> factors() const noexcept {
        return FactorSpan<Factor>(factors_.data(), factors_.size());
    }

    void clear() noexcept {
        factors_.clear();
    }

    void reserve(std::size_t capacity) {
        factors_.reserve(capacity);
    }

    void push(Factor factor, slong exponent) {
        if (exponent == 0) {
            return;
        }
        factors_.push_back(FactorPower<Factor>{std::move(factor), exponent});
    }

    void normalize() {
        std::size_t write = 0;
        for (std::size_t read = 0; read < factors_.size(); ++read) {
            if (factors_[read].exponent == 0) {
                continue;
            }
            if (write != read) {
                factors_[write] = std::move(factors_[read]);
            }
            ++write;
        }
        factors_.resize(write);
    }

private:
    std::vector<FactorPower<Factor>> factors_;
};

}  // namespace silex
