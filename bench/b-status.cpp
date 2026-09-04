#include <benchmark/benchmark.h>

#include <silex/status.hpp>
#include <silex/version.hpp>

#include "benchmark_contract.hpp"

#include <string_view>

#if !defined(SILEX_BENCHMARK_EXPECTED_VERSION_MAJOR) || \
        !defined(SILEX_BENCHMARK_EXPECTED_VERSION_MINOR) || \
        !defined(SILEX_BENCHMARK_EXPECTED_VERSION_PATCH) || \
        !defined(SILEX_BENCHMARK_EXPECTED_VERSION_STRING)
#error "The status benchmark requires the configured project version"
#endif

namespace {

void BM_status_ok(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    for (auto _ : state) {
        benchmark::DoNotOptimize(silex::ok(silex::Status::ok));
    }
    if (!silex::ok(silex::Status::ok)) {
        silex::bench_contract::fail(
                state, "Status::ok is not successful",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_version_string(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    for (auto _ : state) {
        benchmark::DoNotOptimize(silex::version_string());
    }
    const char* observed = silex::version_string();
    const bool exact = observed != nullptr &&
            silex::version_major() ==
                    SILEX_BENCHMARK_EXPECTED_VERSION_MAJOR &&
            silex::version_minor() ==
                    SILEX_BENCHMARK_EXPECTED_VERSION_MINOR &&
            silex::version_patch() ==
                    SILEX_BENCHMARK_EXPECTED_VERSION_PATCH &&
            std::string_view(observed) ==
                    SILEX_BENCHMARK_EXPECTED_VERSION_STRING;
    if (!exact) {
        silex::bench_contract::fail(
                state, "version does not match the configured project version",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    state.counters["version_major"] =
            static_cast<double>(silex::version_major());
    state.counters["version_minor"] =
            static_cast<double>(silex::version_minor());
    state.counters["version_patch"] =
            static_cast<double>(silex::version_patch());
    silex::bench_contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_status_ok);
BENCHMARK(BM_version_string);

BENCHMARK_MAIN();
