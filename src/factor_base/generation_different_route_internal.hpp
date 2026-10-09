#pragma once

namespace silex::detail {

// How a proven generation requirement treats the class [D] of the
// different.
// - none: the requirement is factor_base_class_group_bound (every class
//   contains a small ideal), so the different plays no part.
// - kor1_generation: the requirement is Zimmert's twin-class bound
//   (Zimmert 1981, Korollar 1).  It generates the class group with no
//   hypothesis on [D], by GB-kor1-generation (silex-research,
//   research/proofs/generation/GB-kor1-generation.md at commit 8926120),
//   so [D] needs no check.
enum class GenerationDifferentRoute {
    none = 0,
    kor1_generation = 1
};

}  // namespace silex::detail
