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

// Proven generation bound used by factor_base_class_group_bound for degree
// n = real_places + 2 * complex_pairs >= 3: every ideal class of a field
// with this absolute discriminant and signature contains an integral ideal
// of norm at most out.  It is the smallest of Zimmert's bound (n <= 20),
// Minkowski's bound with (4/pi)^r2, and the exact integer Minkowski form
// with 2^r2.  Fails, leaving out unchanged, if n < 3, the signature is
// negative, or abs_discriminant is not positive.
bool generation_bound(flint::FmpzRef out,
                      flint::FmpzConstRef abs_discriminant,
                      slong real_places,
                      slong complex_pairs) noexcept;

// Zimmert's bound alone (Zimmert 1981, Satz 2 with alpha from Bemerkung 1):
// the ceiling of a rigorous upper bound for sqrt|d| exp(-L(r1, r2, gamma))
// with the tabulated gamma.  Fails, leaving out unchanged, outside
// 3 <= n <= 20 or on the inputs generation_bound rejects.
bool zimmert_generation_bound(flint::FmpzRef out,
                              flint::FmpzConstRef abs_discriminant,
                              slong real_places,
                              slong complex_pairs) noexcept;

}  // namespace silex::detail
