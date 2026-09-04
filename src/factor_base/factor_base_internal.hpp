#pragma once

#include <silex/factor_base.hpp>

namespace silex::detail {

class FactorBaseBlockAccess {
public:
    static bool rational_prime_block_is_complete(
            bool& complete,
            const FactorBase& base,
            slong block_index) noexcept;
};

}  // namespace silex::detail
