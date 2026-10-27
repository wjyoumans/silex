#pragma once

#include <silex/diagnostics.hpp>
#include <silex/factored_element.hpp>

#include <flint/flint.h>

namespace silex::detail {

// Result of factored_element_is_power_si.
//
// `decided`: the test answered; `is_power` says whether `input` is an
// exponent-th power in its field, and on `true` `root` holds a root.  A
// non-power answer is a residue disproof, CompactElement's `not_power`, or
// Element::is_power's decided answer, and is only as strong as that source.
//
// `error`: the arguments were rejected (an undefined input or root, a root
// in another field, an exponent below one, or an unknown strategy).  The
// same arguments fail again on any retry.
//
// `undecided`: the arguments were valid but no answer was reached.  This
// covers every failure inside the test (a failed split, evaluation, product
// or allocation, CompactElement's unsupported root, and an undecided
// Element::is_power or is_square), because those helpers' bool results do
// not separate an internal error from an inconclusive computation.
//
// On `error` and `undecided`, `is_power` and `root` are left unchanged.
enum class FactoredPowerResult : unsigned char {
    error,
    undecided,
    decided,
};

// FactoredElement::is_power_si(is_power, root, exponent, strategy,
// diagnostics) on `input`, with its two failures told apart.  The public
// member returns true exactly when this returns `decided`, with the same
// outputs.
FactoredPowerResult factored_element_is_power_si(
        bool& is_power,
        FactoredElement& root,
        const FactoredElement& input,
        slong exponent,
        FactoredRootStrategy strategy,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

}  // namespace silex::detail
