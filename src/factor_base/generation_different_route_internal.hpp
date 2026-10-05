#pragma once

namespace silex::detail {

// How a proven generation requirement handles the class [D] of the
// different (Zimmert 1981, Korollar 1, needs [D] in the generated subgroup).
// - none: the requirement is factor_base_class_group_bound (every class
//   contains a small ideal), so the different plays no part.
// - index_one: the defining polynomial f is monic integral and
//   O_K = Z[theta]; the conductor of Z[theta] is then O_K (Narkiewicz,
//   Elementary and Analytic Theory of Algebraic Numbers, 3rd ed., 2004,
//   Corollary to Proposition 4.18, p. 153), so the different is
//   f'(theta) O_K, a principal ideal.
// - conductor_primes: f is monic integral and [O_K : Z[theta]] > 1; the
//   requirement also covers every prime ideal above each p dividing the
//   index, which includes every prime dividing the conductor F of Z[theta],
//   and f'(theta) O_K = F D (Narkiewicz 2004, Proposition 4.18(i), p. 152,
//   with Proposition 4.17, p. 151) makes [D] = [F]^(-1).
enum class GenerationDifferentRoute {
    none = 0,
    index_one = 1,
    conductor_primes = 2
};

}  // namespace silex::detail
