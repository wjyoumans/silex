#pragma once

#include <silex/diagnostics.hpp>
#include <flint/flint.h>

#include <atomic>

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

class TestHookAccess;

// Test-only override of the cap computed by precision_doubling_cap; zero
// means off.  Only detail::TestHookAccess (test/test_hook_access.hpp) can set
// it; production code only reads it.  Relaxed atomic so a reader on another
// thread never races with a test setting it.
class PrecisionCapTestSeams {
public:
    static slong forced_cap() noexcept {
        return forced_cap_.load(std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<slong> forced_cap_{0};
};

}  // namespace detail
}  // namespace silex
