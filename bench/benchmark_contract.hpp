#pragma once

#include <benchmark/benchmark.h>

namespace silex::bench_contract {

enum class FailureReason : int {
    not_completed = 1,
    setup = 2,
    operation = 3,
    invariant = 4,
    reference_mismatch = 5,
};

inline void initialize(benchmark::State& state) {
    state.counters["success"] = 0.0;
    state.counters["failure_reason"] =
        static_cast<double>(FailureReason::not_completed);
}

inline void fail(
    benchmark::State& state,
    const char* message,
    FailureReason reason) {
    state.counters["success"] = 0.0;
    state.counters["failure_reason"] = static_cast<double>(reason);
    state.SkipWithError(message);
}

inline void succeed(benchmark::State& state) {
    state.counters["success"] = 1.0;
    state.counters["failure_reason"] = 0.0;
}

}  // namespace silex::bench_contract
