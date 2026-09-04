#include <benchmark/benchmark.h>

#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/lat.hpp>

#include "lat/flatter_backend_internal.hpp"
#include "lat/fplll_backend_internal.hpp"
#include "benchmark_contract.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace {

constexpr slong kDegree14Dimension = 14;

struct BasisQuality {
    double first_norm_squared = 0.0;
    double max_norm_squared = 0.0;
    double sum_norm_squared = 0.0;
    double log2_orthogonality_defect = 0.0;
};

void set_entry_si(fmpz_mat_t matrix, slong row, slong col, slong value) noexcept {
    fmpz_set_si(fmpz_mat_entry(matrix, row, col), value);
}

silex::flint::FmpzMat degree14_hnf_basis() {
    constexpr slong diagonal_base = WORD(1000000000000037);

    silex::flint::FmpzMat basis(kDegree14Dimension, kDegree14Dimension);
    for (slong row = 0; row < kDegree14Dimension; ++row) {
        const slong diagonal = diagonal_base + row * WORD(1000003);
        fmpz_set_si(fmpz_mat_entry(basis.raw(), row, row), diagonal);
        for (slong col = row + 1; col < kDegree14Dimension; ++col) {
            const slong column_diagonal = diagonal_base + col * WORD(1000003);
            const slong value = column_diagonal / (row + 2) +
                                (row + 1) * (col + 1) * WORD(7919);
            fmpz_set_si(fmpz_mat_entry(basis.raw(), row, col), value);
        }
    }
    return basis;
}

silex::flint::FmpzMat native_degree14_ideal_hnf_basis() {
    // Zero-based capture 13 from the full native degree-14 ideal trajectory.
    constexpr std::array<slong,
                         kDegree14Dimension * kDegree14Dimension>
            entries = {
                    1, 0, 0, 0, 0, 0, 0, 5, 3, 49, 10, 0, 45, 41,
                    0, 1, 0, 0, 0, 0, 0, 10, 1, 30, 49, 16, 25, 11,
                    0, 0, 1, 0, 0, 0, 0, 3, 1, 26, 17, 30, 30, 20,
                    0, 0, 0, 1, 0, 0, 0, 9, 2, 37, 13, 12, 33, 15,
                    0, 0, 0, 0, 1, 0, 0, 10, 5, 3, 37, 19, 11, 51,
                    0, 0, 0, 0, 0, 1, 0, 8, 12, 50, 29, 21, 0, 15,
                    0, 0, 0, 0, 0, 0, 1, 10, 4, 13, 50, 35, 20, 18,
                    0, 0, 0, 0, 0, 0, 0, 13, 0, 26, 13, 0, 13, 26,
                    0, 0, 0, 0, 0, 0, 0, 0, 13, 26, 26, 13, 26, 13,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52,
            };

    silex::flint::FmpzMat basis(kDegree14Dimension, kDegree14Dimension);
    for (slong row = 0; row < kDegree14Dimension; ++row) {
        for (slong col = 0; col < kDegree14Dimension; ++col) {
            fmpz_set_si(
                    fmpz_mat_entry(basis.raw(), row, col),
                    entries[static_cast<std::size_t>(
                            row * kDegree14Dimension + col)]);
        }
    }
    return basis;
}

double log2_abs_fmpz(const fmpz* value) noexcept {
    slong exponent = 0;
    const double mantissa = fmpz_get_d_2exp(&exponent, value);
    return std::log2(std::abs(mantissa)) + static_cast<double>(exponent);
}

BasisQuality basis_quality(silex::flint::FmpzMatConstRef basis) noexcept {
    BasisQuality quality;
    silex::flint::Fmpz norm_squared;
    silex::flint::Fmpz determinant;
    double sum_log2_norm_squared = 0.0;
    for (slong row = 0; row < kDegree14Dimension; ++row) {
        fmpz_zero(norm_squared.raw());
        for (slong col = 0; col < kDegree14Dimension; ++col) {
            const fmpz* entry = fmpz_mat_entry(basis.raw(), row, col);
            fmpz_addmul(norm_squared.raw(), entry, entry);
        }
        const double norm_squared_value = fmpz_get_d(norm_squared.raw());
        if (row == 0) {
            quality.first_norm_squared = norm_squared_value;
        }
        quality.max_norm_squared =
            std::max(quality.max_norm_squared, norm_squared_value);
        quality.sum_norm_squared += norm_squared_value;
        sum_log2_norm_squared += log2_abs_fmpz(norm_squared.raw());
    }
    fmpz_mat_det(determinant.raw(), basis.raw());
    quality.log2_orthogonality_defect =
        0.5 * sum_log2_norm_squared - log2_abs_fmpz(determinant.raw());
    return quality;
}

bool row_transform_is_valid(silex::flint::FmpzMatConstRef input,
                            silex::flint::FmpzMatConstRef reduced,
                            silex::flint::FmpzMatConstRef transform) noexcept {
    silex::flint::FmpzMat expected(kDegree14Dimension, kDegree14Dimension);
    silex::flint::Fmpz determinant;
    fmpz_mat_mul(expected.raw(), transform.raw(), input.raw());
    fmpz_mat_det(determinant.raw(), transform.raw());
    return fmpz_mat_equal(expected.raw(), reduced.raw()) != 0 &&
           fmpz_is_pm1(determinant.raw()) != 0;
}

bool lll_result_is_valid(silex::flint::FmpzMatConstRef input,
                         silex::flint::FmpzMatConstRef reduced) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong columns = fmpz_mat_ncols(input.raw());
    if (fmpz_mat_nrows(reduced.raw()) != rows ||
        fmpz_mat_ncols(reduced.raw()) != columns) {
        return false;
    }

    silex::flint::FmpzMat input_hnf(rows, columns);
    silex::flint::FmpzMat reduced_hnf(rows, columns);
    silex::flint::FmpzLll config;
    fmpz_mat_hnf(input_hnf.raw(), input.raw());
    fmpz_mat_hnf(reduced_hnf.raw(), reduced.raw());
    return fmpz_mat_equal(input_hnf.raw(), reduced_hnf.raw()) != 0 &&
           fmpz_lll_is_reduced(reduced.raw(), config.raw(), 0) != 0;
}

void publish_quality(benchmark::State& state,
                     silex::flint::FmpzMatConstRef basis) noexcept {
    const BasisQuality quality = basis_quality(basis);
    state.counters["first_norm_sq"] = quality.first_norm_squared;
    state.counters["max_norm_sq"] = quality.max_norm_squared;
    state.counters["sum_norm_sq"] = quality.sum_norm_squared;
    state.counters["log2_orthogonality_defect"] =
        quality.log2_orthogonality_defect;
}

void BM_lat_lll(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::lat::Lat lattice(4);
    silex::lat::Lat reduced(4);
    fmpz_mat_t basis;

    fmpz_mat_init(basis, 4, 4);

    set_entry_si(basis, 0, 0, 105);
    set_entry_si(basis, 0, 1, 821);
    set_entry_si(basis, 0, 2, 17);
    set_entry_si(basis, 0, 3, 9);
    set_entry_si(basis, 1, 0, 37);
    set_entry_si(basis, 1, 1, 19);
    set_entry_si(basis, 1, 2, 401);
    set_entry_si(basis, 1, 3, 11);
    set_entry_si(basis, 2, 0, 2);
    set_entry_si(basis, 2, 1, 3);
    set_entry_si(basis, 2, 2, 5);
    set_entry_si(basis, 2, 3, 7);
    set_entry_si(basis, 3, 0, 13);
    set_entry_si(basis, 3, 1, 29);
    set_entry_si(basis, 3, 2, 31);
    set_entry_si(basis, 3, 3, 37);
    if (!lattice.set_basis(basis)) {
        silex::bench_contract::fail(
                state, "lattice setup failed",
                silex::bench_contract::FailureReason::setup);
        fmpz_mat_clear(basis);
        return;
    }

    bool operation_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = lattice.lll_reduce(reduced);
        operation_ok = iteration_ok && operation_ok;
        benchmark::DoNotOptimize(reduced.nrows());
    }

    benchmark::ClobberMemory();
    const bool valid =
            operation_ok && reduced.nrows() == 4 &&
            lll_result_is_valid(lattice.basis_ref(), reduced.basis_ref());
    if (!valid) {
        silex::bench_contract::fail(
                state, "LLL result failed lattice and reduction checks",
                operation_ok
                        ? silex::bench_contract::FailureReason::invariant
                        : silex::bench_contract::FailureReason::operation);
        fmpz_mat_clear(basis);
        return;
    }
    state.counters["rows"] = static_cast<double>(reduced.nrows());
    silex::bench_contract::succeed(state);
    fmpz_mat_clear(basis);
}

void BM_lat_lll_hnf_degree14(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::flint::FmpzMat basis = degree14_hnf_basis();
    silex::lat::Lat lattice(kDegree14Dimension);
    silex::lat::Lat reduced(kDegree14Dimension);
    if (!lattice.set_basis(basis) || !lattice.is_hnf()) {
        silex::bench_contract::fail(
                state, "degree-14 input is not HNF",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool operation_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = lattice.lll_reduce(reduced);
        operation_ok = iteration_ok && operation_ok;
        benchmark::DoNotOptimize(reduced.nrows());
    }

    benchmark::ClobberMemory();
    const bool valid = operation_ok &&
                       reduced.nrows() == kDegree14Dimension &&
                       lll_result_is_valid(lattice.basis_ref(),
                                           reduced.basis_ref());
    if (!valid) {
        silex::bench_contract::fail(
                state, "degree-14 LLL result failed validation",
                operation_ok
                        ? silex::bench_contract::FailureReason::invariant
                        : silex::bench_contract::FailureReason::operation);
        return;
    }
    publish_quality(state, reduced.basis_ref());
    state.counters["rows"] = static_cast<double>(reduced.nrows());
    silex::bench_contract::succeed(state);
}

void BM_native_ideal_hnf_degree14_flint(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::flint::FmpzMat basis = native_degree14_ideal_hnf_basis();
    silex::flint::FmpzMat reduced(kDegree14Dimension, kDegree14Dimension);
    silex::flint::FmpzMat transform(kDegree14Dimension, kDegree14Dimension);
    silex::flint::FmpzLll config;

    for (auto _ : state) {
        fmpz_mat_set(reduced.raw(), basis.raw());
        fmpz_mat_one(transform.raw());
        fmpz_lll(reduced.raw(), transform.raw(), config.raw());
        benchmark::DoNotOptimize(fmpz_mat_entry(reduced.raw(), 0, 0));
    }
    const bool valid =
            row_transform_is_valid(silex::flint::FmpzMatConstRef(basis),
                                   silex::flint::FmpzMatConstRef(reduced),
                                   silex::flint::FmpzMatConstRef(transform)) &&
            lll_result_is_valid(silex::flint::FmpzMatConstRef(basis),
                                silex::flint::FmpzMatConstRef(reduced));
    if (!valid) {
        silex::bench_contract::fail(
                state,
                "FLINT LLL transform or reduced basis failed validation",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    benchmark::ClobberMemory();
    publish_quality(state, silex::flint::FmpzMatConstRef(reduced));
    silex::bench_contract::succeed(state);
}

void BM_native_ideal_hnf_degree14_fplll(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::flint::FmpzMat basis = native_degree14_ideal_hnf_basis();
    silex::flint::FmpzMat reduced(kDegree14Dimension, kDegree14Dimension);
    silex::flint::FmpzMat transform(kDegree14Dimension, kDegree14Dimension);

    silex::lat::detail::FplllLllResult result;
    for (auto _ : state) {
        result = silex::lat::detail::fplll_row_lll_transform(
            silex::flint::FmpzMatRef(reduced),
            silex::flint::FmpzMatRef(transform),
            silex::flint::FmpzMatConstRef(basis), 0.99);
        benchmark::DoNotOptimize(result.backend_status);
    }
    if (result.status == silex::lat::detail::FplllBackendStatus::unavailable) {
        silex::bench_contract::fail(
                state, "fplll backend unavailable",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    if (result.status != silex::lat::detail::FplllBackendStatus::success ||
        !row_transform_is_valid(silex::flint::FmpzMatConstRef(basis),
                                silex::flint::FmpzMatConstRef(reduced),
                                silex::flint::FmpzMatConstRef(transform))) {
        silex::bench_contract::fail(
                state, "fplll LLL transform failed validation",
                result.status ==
                                silex::lat::detail::FplllBackendStatus::success
                        ? silex::bench_contract::FailureReason::invariant
                        : silex::bench_contract::FailureReason::operation);
        return;
    }
    benchmark::ClobberMemory();
    publish_quality(state, silex::flint::FmpzMatConstRef(reduced));
    silex::bench_contract::succeed(state);
}

void BM_native_ideal_hnf_degree14_flatter_rhf_1_02(
        benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::flint::FmpzMat basis = native_degree14_ideal_hnf_basis();
    silex::flint::FmpzMat transposed(kDegree14Dimension, kDegree14Dimension);
    silex::flint::FmpzMat reduced_transposed(kDegree14Dimension,
                                             kDegree14Dimension);
    silex::flint::FmpzMat column_transform(kDegree14Dimension,
                                           kDegree14Dimension);
    silex::flint::FmpzMat reduced(kDegree14Dimension, kDegree14Dimension);
    silex::flint::FmpzMat row_transform(kDegree14Dimension, kDegree14Dimension);
    fmpz_mat_transpose(transposed.raw(), basis.raw());

    silex::lat::detail::FlatterLllResult result;
    for (auto _ : state) {
        result = silex::lat::detail::flatter_column_lll_transform(
            silex::flint::FmpzMatRef(reduced_transposed),
            silex::flint::FmpzMatRef(column_transform),
            silex::flint::FmpzMatConstRef(transposed), 1.02, 1);
        fmpz_mat_transpose(reduced.raw(), reduced_transposed.raw());
        fmpz_mat_transpose(row_transform.raw(), column_transform.raw());
        benchmark::DoNotOptimize(result.rank);
    }
    if (result.status ==
        silex::lat::detail::FlatterBackendStatus::unavailable) {
        silex::bench_contract::fail(
                state, "flatter backend unavailable",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    if (result.status != silex::lat::detail::FlatterBackendStatus::success ||
        !row_transform_is_valid(silex::flint::FmpzMatConstRef(basis),
                                silex::flint::FmpzMatConstRef(reduced),
                                silex::flint::FmpzMatConstRef(row_transform))) {
        silex::bench_contract::fail(
                state, "flatter transform failed validation",
                result.status ==
                                silex::lat::detail::FlatterBackendStatus::success
                        ? silex::bench_contract::FailureReason::invariant
                        : silex::bench_contract::FailureReason::operation);
        return;
    }
    benchmark::ClobberMemory();
    publish_quality(state, silex::flint::FmpzMatConstRef(reduced));
    silex::bench_contract::succeed(state);
}

void BM_native_ideal_hnf_degree14_fplll_bkz(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    const int block_size = static_cast<int>(state.range(0));
    silex::flint::FmpzMat basis = native_degree14_ideal_hnf_basis();
    silex::flint::FmpzMat reduced(kDegree14Dimension, kDegree14Dimension);
    silex::flint::FmpzMat transform(kDegree14Dimension, kDegree14Dimension);

    silex::lat::detail::FplllLllResult result;
    for (auto _ : state) {
        result = silex::lat::detail::fplll_row_bkz_transform(
            silex::flint::FmpzMatRef(reduced),
            silex::flint::FmpzMatRef(transform),
            silex::flint::FmpzMatConstRef(basis), block_size, 1);
        benchmark::DoNotOptimize(result.backend_status);
    }
    if (result.status == silex::lat::detail::FplllBackendStatus::unavailable) {
        silex::bench_contract::fail(
                state, "fplll backend unavailable",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    if (result.status != silex::lat::detail::FplllBackendStatus::success ||
        !row_transform_is_valid(silex::flint::FmpzMatConstRef(basis),
                                silex::flint::FmpzMatConstRef(reduced),
                                silex::flint::FmpzMatConstRef(transform))) {
        silex::bench_contract::fail(
                state, "fplll BKZ transform failed validation",
                result.status ==
                                silex::lat::detail::FplllBackendStatus::success
                        ? silex::bench_contract::FailureReason::invariant
                        : silex::bench_contract::FailureReason::operation);
        return;
    }
    benchmark::ClobberMemory();
    publish_quality(state, silex::flint::FmpzMatConstRef(reduced));
    silex::bench_contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_lat_lll);
BENCHMARK(BM_lat_lll_hnf_degree14);
BENCHMARK(BM_native_ideal_hnf_degree14_flint);
BENCHMARK(BM_native_ideal_hnf_degree14_fplll);
BENCHMARK(BM_native_ideal_hnf_degree14_flatter_rhf_1_02);
BENCHMARK(BM_native_ideal_hnf_degree14_fplll_bkz)
    ->ArgName("block_size")
    ->Arg(4)
    ->Arg(8)
    ->Arg(10)
    ->Arg(12)
    ->Arg(14);
BENCHMARK_MAIN();
