#pragma once

#include <silex/diagnostics.hpp>
#include <flint/flint.h>

namespace silex {
namespace detail {

// Shared cap for precision-doubling loops on the class/unit path (compact
// independence, compact regulator to a tolerance, factored logarithmic
// embedding, conjugate-log cutoff inverse).  The loops terminate in theory,
// but memory grows with the compact exponent size, so each stops at
// max(64 * requested precision, 2^20) bits and fails closed through its
// existing failure path.  Not part of the installed headers.
inline constexpr slong kPrecisionDoublingCapFloorBits = WORD(1) << 20;
inline constexpr slong kPrecisionDoublingCapFactor = 64;

// The cap, in bits, for a loop whose caller requested `requested` bits.
slong precision_doubling_cap(slong requested) noexcept;

// True when `work_precision` may be doubled without exceeding `cap`.  When
// it may not, logs `reason` at detail level and returns false.
bool precision_doubling_allowed(slong work_precision,
                                slong cap,
                                const DiagnosticsContext* diagnostics,
                                const char* reason) noexcept;

// Test hook: a positive value replaces the cap computed by
// precision_doubling_cap; zero restores the default.
void set_precision_doubling_cap_for_testing(slong cap) noexcept;

}  // namespace detail
}  // namespace silex
