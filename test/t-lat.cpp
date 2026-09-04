#include <silex/lat.hpp>

#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>

#include "lat/fplll_backend_internal.hpp"
#include "lat/flatter_backend_internal.hpp"

#include <type_traits>
#include <utility>

namespace {

struct EnumCounter {
    slong count = 0;
    slong abort_after = 0;
    slong max_abs_coord = 0;
};

struct DebugCounter {
    int failures = 0;
    int verbose = 0;
    const char* last_label = nullptr;
    const char* last_expression = nullptr;
};

void set_entry_si(fmpz_mat_t matrix, slong row, slong col, slong value) noexcept {
    fmpz_set_si(fmpz_mat_entry(matrix, row, col), value);
}

bool matrix_equals_si(const fmpz_mat_t matrix,
        const slong* expected,
        slong rows,
        slong cols) noexcept {
    if (fmpz_mat_nrows(matrix) != rows || fmpz_mat_ncols(matrix) != cols) {
        return false;
    }
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < cols; ++col) {
            if (fmpz_cmp_si(fmpz_mat_entry(matrix, row, col),
                        expected[row * cols + col]) != 0) {
                return false;
            }
        }
    }
    return true;
}

void debug_failure_callback(void* user,
        silex::DiagnosticsModule module,
        silex::DebugLevel,
        const char*,
        const char*,
        int,
        const char* label,
        const char* expression) noexcept {
    auto* counter = static_cast<DebugCounter*>(user);
    if (module == silex::DiagnosticsModule::lattice) {
        ++counter->failures;
        counter->last_label = label;
        counter->last_expression = expression;
    }
}

void verbose_callback(void* user,
        silex::DiagnosticsModule module,
        silex::VerboseLevel,
        const char*,
        const char*,
        const char*) noexcept {
    auto* counter = static_cast<DebugCounter*>(user);
    if (module == silex::DiagnosticsModule::lattice) {
        ++counter->verbose;
    }
}

int enum_count_callback(const fmpz_mat_t coeffs, void* user) {
    auto* counter = static_cast<EnumCounter*>(user);
    const slong cols = fmpz_mat_ncols(coeffs);
    for (slong j = 0; j < cols; ++j) {
        slong value = fmpz_get_si(fmpz_mat_entry(coeffs, 0, j));
        if (value < 0) {
            value = -value;
        }
        if (value > counter->max_abs_coord) {
            counter->max_abs_coord = value;
        }
    }

    ++counter->count;
    return counter->abort_after > 0 && counter->count >= counter->abort_after ? 0
                                                                              : 1;
}

bool lat_equal_hnf(const silex::lat::Lat& left, const silex::lat::Lat& right) {
    silex::lat::Lat left_hnf(left.ambient_dim());
    silex::lat::Lat right_hnf(right.ambient_dim());
    if (!left.hnf(left_hnf) || !right.hnf(right_hnf) ||
        left_hnf.nrows() != right_hnf.nrows()) {
        return false;
    }
    fmpz_mat_t left_basis;
    fmpz_mat_t right_basis;
    fmpz_mat_init(left_basis, left_hnf.nrows(), left_hnf.ambient_dim());
    fmpz_mat_init(right_basis, right_hnf.nrows(), right_hnf.ambient_dim());
    left_hnf.get_basis(left_basis);
    right_hnf.get_basis(right_basis);
    const bool equal = fmpz_mat_equal(left_basis, right_basis) != 0;
    fmpz_mat_clear(right_basis);
    fmpz_mat_clear(left_basis);
    return equal;
}

int test_init_set_swap_basis() {
    silex::lat::Lat lattice(2);
    silex::lat::Lat other(3);
    fmpz_mat_t basis;
    fmpz_mat_t got;

    fmpz_mat_init(basis, 2, 2);
    fmpz_mat_init(got, 2, 2);

    if (lattice.ambient_dim() != 2 || lattice.nrows() != 0 ||
        !lattice.is_hnf()) {
        return 1;
    }

    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 1, 1, 3);
    if (!lattice.set_basis(basis) || lattice.ambient_dim() != 2 ||
        lattice.nrows() != 2 || !lattice.get_basis(got) ||
        !fmpz_mat_equal(basis, got)) {
        return 1;
    }
    silex::flint::FmpzMat owned_basis = lattice.basis();
    if (fmpz_mat_equal(basis, owned_basis.raw()) == 0) {
        return 1;
    }

    if (!lattice.set_basis(lattice.basis_ref()) || !lattice.get_basis(got) ||
        !fmpz_mat_equal(basis, got)) {
        return 1;
    }

    if (!other.set(lattice) || other.ambient_dim() != 2 ||
        !other.get_basis(got) || !fmpz_mat_equal(basis, got)) {
        return 1;
    }
    if (!lattice.set(lattice) || lattice.ambient_dim() != 2 ||
        !lattice.get_basis(got) || !fmpz_mat_equal(basis, got)) {
        return 1;
    }
    lattice.swap(other);
    if (lattice.ambient_dim() != 2 || other.ambient_dim() != 2 ||
        !lattice.get_basis(got) || !fmpz_mat_equal(basis, got)) {
        return 1;
    }

    fmpz_mat_clear(got);
    fmpz_mat_clear(basis);
    return 0;
}

int test_hnf_and_transform() {
    silex::lat::Lat lattice(2);
    silex::lat::Lat hnf(2);
    silex::lat::Lat hnf_twice(2);
    silex::lat::Lat transformed(2);
    fmpz_mat_t basis;
    fmpz_mat_t transform;
    fmpz_mat_t product;

    fmpz_mat_init(basis, 3, 2);
    fmpz_mat_init(transform, 3, 3);
    fmpz_mat_init(product, 3, 2);
    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 1, 1, 2);
    set_entry_si(basis, 2, 0, 1);
    set_entry_si(basis, 2, 1, 1);

    if (!lattice.set_basis(basis) || !lattice.hnf(hnf) || hnf.nrows() != 2 ||
        hnf.ambient_dim() != 2 || !hnf.is_hnf() || !hnf.hnf(hnf_twice) ||
        !lat_equal_hnf(hnf, hnf_twice)) {
        return 1;
    }

    if (!lattice.hnf(lattice) || !lat_equal_hnf(lattice, hnf)) {
        return 1;
    }

    lattice.set_basis(basis);
    if (!lattice.hnf_transform(transformed, transform)) {
        return 1;
    }
    fmpz_mat_mul(product, transform, basis);
    if (!fmpz_mat_equal(product, transformed.basis_ref().raw()) ||
        fmpz_mat_is_in_hnf(transformed.basis_ref().raw()) == 0) {
        return 1;
    }

    lattice.set_basis(basis);
    if (!lattice.hnf_transform(lattice, transform)) {
        return 1;
    }
    fmpz_mat_mul(product, transform, basis);
    if (!fmpz_mat_equal(product, lattice.basis_ref().raw()) ||
        fmpz_mat_is_in_hnf(lattice.basis_ref().raw()) == 0) {
        return 1;
    }

    fmpz_mat_clear(product);
    fmpz_mat_clear(transform);
    fmpz_mat_clear(basis);
    return 0;
}

int test_contains() {
    silex::lat::Lat lattice(2);
    silex::lat::Lat zero(2);
    fmpz_mat_t basis;
    fmpz_mat_t row_vectors;
    fmpz_mat_t wrong_dim;
    fmpz vec[2];
    fmpz zeros[2];

    fmpz_init(vec + 0);
    fmpz_init(vec + 1);
    fmpz_init(zeros + 0);
    fmpz_init(zeros + 1);
    fmpz_mat_init(basis, 2, 2);
    fmpz_mat_init(row_vectors, 2, 2);
    fmpz_mat_init(wrong_dim, 1, 3);

    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 1, 1, 3);
    lattice.set_basis(basis);

    fmpz_set_si(vec + 0, 4);
    fmpz_set_si(vec + 1, -6);
    if (!lattice.contains(vec)) {
        return 1;
    }
    set_entry_si(row_vectors, 0, 0, 4);
    set_entry_si(row_vectors, 0, 1, -6);
    set_entry_si(row_vectors, 1, 0, 1);
    if (!lattice.contains_row(row_vectors, 0) ||
        lattice.contains_row(row_vectors, 1) ||
        lattice.contains_row(row_vectors, -1) ||
        lattice.contains_row(row_vectors, 2) ||
        lattice.contains_row(wrong_dim, 0)) {
        return 1;
    }
    fmpz_set_si(vec + 0, 1);
    fmpz_zero(vec + 1);
    if (lattice.contains(vec) || !zero.contains(zeros) || zero.contains(vec)) {
        return 1;
    }

    fmpz_mat_clear(wrong_dim);
    fmpz_mat_clear(row_vectors);
    fmpz_mat_clear(basis);
    fmpz_clear(zeros + 1);
    fmpz_clear(zeros + 0);
    fmpz_clear(vec + 1);
    fmpz_clear(vec + 0);
    return 0;
}

int test_sum_intersection_index() {
    silex::lat::Lat z2(2);
    silex::lat::Lat lattice(2);
    silex::lat::Lat other(2);
    silex::lat::Lat sum(2);
    silex::lat::Lat intersection(2);
    silex::lat::Lat expected(2);
    fmpz_mat_t basis;
    fmpz_mat_t other_basis;
    fmpz_t index;
    fmpz vec[2];

    fmpz_init(index);
    fmpz_init(vec + 0);
    fmpz_init(vec + 1);
    fmpz_mat_init(basis, 2, 2);
    fmpz_mat_init(other_basis, 2, 2);

    set_entry_si(basis, 0, 0, 1);
    set_entry_si(basis, 1, 1, 1);
    z2.set_basis(basis);

    fmpz_mat_zero(basis);
    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 1, 1, 1);
    lattice.set_basis(basis);

    fmpz_mat_zero(other_basis);
    set_entry_si(other_basis, 0, 0, 1);
    set_entry_si(other_basis, 1, 1, 3);
    other.set_basis(other_basis);

    if (!lattice.sum(sum, other) || !lat_equal_hnf(sum, z2) ||
        !lattice.sum(lattice, other) || !lat_equal_hnf(lattice, z2)) {
        return 1;
    }

    fmpz_mat_zero(basis);
    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 1, 1, 1);
    lattice.set_basis(basis);

    if (!lattice.intersection(intersection, other)) {
        return 1;
    }
    fmpz_mat_zero(other_basis);
    set_entry_si(other_basis, 0, 0, 2);
    set_entry_si(other_basis, 1, 1, 3);
    expected.set_basis(other_basis);
    if (!lat_equal_hnf(intersection, expected)) {
        return 1;
    }

    fmpz_set_si(vec + 0, 2);
    fmpz_set_si(vec + 1, 3);
    if (!intersection.contains(vec)) {
        return 1;
    }
    fmpz_set_si(vec + 0, 1);
    if (intersection.contains(vec)) {
        return 1;
    }

    if (!lattice.intersection(lattice, other) || !lat_equal_hnf(lattice, expected)) {
        return 1;
    }

    if (!z2.index(index, expected) || !fmpz_equal_si(index, 6)) {
        return 1;
    }
    fmpz_set_si(index, 6);
    if (expected.index(index, z2) || !fmpz_equal_si(index, 6) ||
        !other.index(index, expected) || !fmpz_equal_si(index, 2)) {
        return 1;
    }

    fmpz_mat_zero(basis);
    set_entry_si(basis, 0, 0, 1);
    set_entry_si(basis, 0, 1, 1);
    set_entry_si(basis, 1, 0, 1);
    set_entry_si(basis, 1, 1, -1);
    lattice.set_basis(basis);
    fmpz_mat_zero(other_basis);
    set_entry_si(other_basis, 0, 0, 2);
    set_entry_si(other_basis, 1, 1, 2);
    other.set_basis(other_basis);
    if (!lattice.intersection(intersection, other) ||
        !lat_equal_hnf(intersection, other)) {
        return 1;
    }

    fmpz_clear(vec + 1);
    fmpz_clear(vec + 0);
    fmpz_clear(index);
    fmpz_mat_clear(other_basis);
    fmpz_mat_clear(basis);
    return 0;
}

int test_saturate() {
    silex::lat::Lat lattice(2);
    silex::lat::Lat saturated(2);
    silex::lat::Lat expected(2);
    fmpz_mat_t basis;
    fmpz_mat_t expected_basis;
    fmpz_t p;
    fmpz vec[2];

    fmpz_mat_init(basis, 2, 2);
    fmpz_mat_init(expected_basis, 2, 2);
    fmpz_init(p);
    fmpz_init(vec + 0);
    fmpz_init(vec + 1);

    set_entry_si(basis, 0, 0, 4);
    set_entry_si(basis, 1, 1, 6);
    lattice.set_basis(basis);

    fmpz_set_si(p, 2);
    if (!lattice.saturate(saturated, p)) {
        return 1;
    }
    set_entry_si(expected_basis, 0, 0, 1);
    set_entry_si(expected_basis, 1, 1, 3);
    expected.set_basis(expected_basis);
    if (!lat_equal_hnf(saturated, expected)) {
        return 1;
    }

    fmpz_set_si(p, 3);
    if (!lattice.saturate(saturated, p)) {
        return 1;
    }
    fmpz_mat_zero(expected_basis);
    set_entry_si(expected_basis, 0, 0, 4);
    set_entry_si(expected_basis, 1, 1, 2);
    expected.set_basis(expected_basis);
    if (!lat_equal_hnf(saturated, expected)) {
        return 1;
    }

    fmpz_mat_zero(basis);
    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 0, 1, 2);
    lattice.set_basis(basis);
    fmpz_set_si(p, 2);
    if (!lattice.saturate(lattice, p)) {
        return 1;
    }
    fmpz_mat_zero(expected_basis);
    set_entry_si(expected_basis, 0, 0, 1);
    set_entry_si(expected_basis, 0, 1, 1);
    expected.set_basis(expected_basis);
    if (!lat_equal_hnf(lattice, expected)) {
        return 1;
    }

    fmpz_set_si(vec + 0, 1);
    fmpz_set_si(vec + 1, 1);
    if (!lattice.contains(vec)) {
        return 1;
    }

    saturated.set_basis(basis);
    fmpz_set_si(p, 4);
    if (lattice.saturate(saturated, p)) {
        return 1;
    }
    expected.set_basis(basis);
    if (!lat_equal_hnf(saturated, expected)) {
        return 1;
    }

    fmpz_clear(vec + 1);
    fmpz_clear(vec + 0);
    fmpz_clear(p);
    fmpz_mat_clear(expected_basis);
    fmpz_mat_clear(basis);
    return 0;
}

int test_lll_reduce() {
    silex::lat::Lat lattice(3);
    silex::lat::Lat reduced(3);
    fmpz_mat_t basis;

    fmpz_mat_init(basis, 3, 3);
    set_entry_si(basis, 0, 0, 105);
    set_entry_si(basis, 0, 1, 821);
    set_entry_si(basis, 0, 2, 17);
    set_entry_si(basis, 1, 0, 37);
    set_entry_si(basis, 1, 1, 19);
    set_entry_si(basis, 1, 2, 401);
    set_entry_si(basis, 2, 0, 2);
    set_entry_si(basis, 2, 1, 3);
    set_entry_si(basis, 2, 2, 5);
    lattice.set_basis(basis);

    if (!lattice.lll_reduce(reduced) || !lat_equal_hnf(reduced, lattice) ||
        reduced.is_hnf()) {
        return 1;
    }

    silex::flint::FmpzLll config;
    if (fmpz_lll_is_reduced(reduced.raw_basis(), config.raw(), 0) == 0) {
        return 1;
    }

    if (!lattice.lll_reduce(lattice) || !lat_equal_hnf(lattice, reduced)) {
        return 1;
    }

    silex::lat::Lat empty(3);
    if (!empty.lll_reduce(reduced) || reduced.nrows() != 0 ||
        reduced.ambient_dim() != 3) {
        return 1;
    }

    constexpr slong hnf_dimension = 14;
    silex::flint::FmpzMat hnf_basis(hnf_dimension + 1, hnf_dimension);
    for (slong row = 0; row < hnf_dimension; ++row) {
        const slong diagonal = WORD(1000003) + row * WORD(1009);
        set_entry_si(hnf_basis.raw(), row, row, diagonal);
        for (slong col = row + 1; col < hnf_dimension; ++col) {
            const slong column_diagonal = WORD(1000003) + col * WORD(1009);
            set_entry_si(hnf_basis.raw(), row, col,
                    column_diagonal / (row + 2) +
                            (row + 1) * (col + 1));
        }
    }

    silex::lat::Lat hnf_lattice(hnf_dimension);
    silex::lat::Lat hnf_reduced(hnf_dimension);
    if (!hnf_lattice.set_basis(hnf_basis) || !hnf_lattice.is_hnf()) {
        return 1;
    }

    silex::flint::FmpzMat expected(hnf_dimension, hnf_dimension);
    silex::flint::FmpzMat transform(hnf_dimension, hnf_dimension);
    silex::flint::FmpzMatConstWindow hnf_window(
            hnf_basis, 0, 0, hnf_dimension, hnf_dimension);
    silex::flint::fmpz_mat_set(expected, hnf_window.const_ref());
    silex::flint::FmpzLll hnf_config;
    fmpz_lll(expected.raw(), transform.raw(), hnf_config.raw());

    if (!hnf_lattice.lll_reduce(hnf_reduced) ||
        hnf_reduced.nrows() != hnf_dimension ||
        fmpz_mat_equal(hnf_reduced.raw_basis(), expected.raw()) == 0 ||
        !hnf_lattice.lll_reduce(hnf_lattice) ||
        fmpz_mat_equal(hnf_lattice.raw_basis(), expected.raw()) == 0) {
        return 1;
    }

    fmpz_mat_clear(basis);
    return 0;
}

int test_fplll_row_transform_boundary() {
    constexpr slong rows = 3;
    constexpr slong cols = 9;
    const slong input_entries[rows * cols] = {
            32373176721380998, 0, 0, -26155458657878188,
            -27972425631556, 976620610054834, 41529514241169880,
            25748372788561934, 7006,
            0, 32373176721380998, 0, 7703955618578414,
            -1359153410186158, -3363101440953616,
            -29118641599698938, 81558150271772994, 9420,
            0, 0, 32373176721380998, 10949198571930016,
            1275236133291490, 6292963271118118,
            56587654159065584, -52872796988158689, 11598,
    };
    const slong expected_transform[rows * rows] = {
            0, 1, 1,
            1, 0, 0,
            0, 0, -1,
    };

    silex::flint::FmpzMat input(rows, cols);
    silex::flint::FmpzMat reduced(rows, cols);
    silex::flint::FmpzMat transform(rows, rows);
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < cols; ++col) {
            set_entry_si(input.raw(), row, col,
                    input_entries[row * cols + col]);
        }
    }

    const auto result =
            silex::lat::detail::fplll_row_lll_transform(
                    silex::flint::FmpzMatRef(reduced),
                    silex::flint::FmpzMatRef(transform),
                    silex::flint::FmpzMatConstRef(input), 0.99);
    if (result.status == silex::lat::detail::FplllBackendStatus::unavailable) {
        return 0;
    }
    if (result.status != silex::lat::detail::FplllBackendStatus::success ||
        result.backend_status != 0 ||
        !matrix_equals_si(transform.raw(), expected_transform, rows, rows)) {
        return 1;
    }

    silex::flint::FmpzMat expected_reduced(rows, cols);
    fmpz_mat_mul(expected_reduced.raw(), transform.raw(), input.raw());
    return fmpz_mat_equal(expected_reduced.raw(), reduced.raw()) == 0 ? 1 : 0;
}

int test_fplll_column_image_transform_boundary() {
    constexpr slong rows = 3;
    constexpr slong cols = 9;
    const slong input_entries[rows * cols] = {
            32373176721380998, 0, 0, -26155458657878188,
            -27972425631556, 976620610054834, 41529514241169880,
            25748372788561934, 7006,
            0, 32373176721380998, 0, 7703955618578414,
            -1359153410186158, -3363101440953616,
            -29118641599698938, 81558150271772994, 9420,
            0, 0, 32373176721380998, 10949198571930016,
            1275236133291490, 6292963271118118,
            56587654159065584, -52872796988158689, 11598,
    };

    silex::flint::FmpzMat input(rows, cols);
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < cols; ++col) {
            set_entry_si(input.raw(), row, col,
                    input_entries[row * cols + col]);
        }
    }

    const slong rank = fmpz_mat_rank(input.raw());
    silex::flint::FmpzMat reduced(rows, rank);
    silex::flint::FmpzMat transform(cols, rank);
    const auto result =
            silex::lat::detail::fplll_column_image_lll_transform(
                    silex::flint::FmpzMatRef(reduced),
                    silex::flint::FmpzMatRef(transform),
                    silex::flint::FmpzMatConstRef(input), 0.99);
    if (result.status == silex::lat::detail::FplllBackendStatus::unavailable) {
        return 0;
    }
    if (result.status != silex::lat::detail::FplllBackendStatus::success ||
        result.backend_status != 0 ||
        fmpz_mat_rank(reduced.raw()) != rank) {
        return 1;
    }

    silex::flint::FmpzMat expected_reduced(rows, rank);
    fmpz_mat_mul(expected_reduced.raw(), input.raw(), transform.raw());
    return fmpz_mat_equal(expected_reduced.raw(), reduced.raw()) == 0 ? 1 : 0;
}

int test_fplll_bounded_bkz_row_transform_boundary() {
    constexpr slong rows = 4;
    constexpr slong cols = 4;
    const slong input_entries[rows * cols] = {
            105, 821, 17, 9,
            37, 19, 401, 11,
            2, 3, 5, 7,
            13, 29, 31, 37,
    };

    silex::flint::FmpzMat input(rows, cols);
    silex::flint::FmpzMat reduced(rows, cols);
    silex::flint::FmpzMat transform(rows, rows);
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < cols; ++col) {
            set_entry_si(input.raw(), row, col,
                         input_entries[row * cols + col]);
        }
    }

    const auto result = silex::lat::detail::fplll_row_bkz_transform(
            silex::flint::FmpzMatRef(reduced),
            silex::flint::FmpzMatRef(transform),
            silex::flint::FmpzMatConstRef(input), 3, 1);
    if (result.status == silex::lat::detail::FplllBackendStatus::unavailable) {
        return 0;
    }
    if (result.status != silex::lat::detail::FplllBackendStatus::success ||
        fmpz_mat_rank(reduced.raw()) != rows) {
        return 1;
    }

    silex::flint::Fmpz determinant;
    fmpz_mat_det(determinant.raw(), transform.raw());
    silex::flint::FmpzMat expected_reduced(rows, cols);
    fmpz_mat_mul(expected_reduced.raw(), transform.raw(), input.raw());
    return fmpz_is_pm1(determinant.raw()) == 0 ||
                   fmpz_mat_equal(expected_reduced.raw(), reduced.raw()) == 0
           ? 1
           : 0;
}

int test_flatter_full_rank_column_transform_boundary() {
    constexpr slong rows = 2;
    constexpr slong cols = 2;
    const slong input_entries[rows * cols] = {
            105, 37,
            821, 19,
    };

    silex::flint::FmpzMat input(rows, cols);
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < cols; ++col) {
            set_entry_si(input.raw(), row, col,
                    input_entries[row * cols + col]);
        }
    }

    const slong rank = fmpz_mat_rank(input.raw());
    silex::flint::FmpzMat reduced(rows, rank);
    silex::flint::FmpzMat transform(cols, rank);
    const auto result =
            silex::lat::detail::flatter_column_lll_transform(
                    silex::flint::FmpzMatRef(reduced),
                    silex::flint::FmpzMatRef(transform),
                    silex::flint::FmpzMatConstRef(input), 1.02, 1);
    if (result.status == silex::lat::detail::FlatterBackendStatus::unavailable) {
        return 0;
    }
    if (result.status != silex::lat::detail::FlatterBackendStatus::success ||
        result.rank != rank || fmpz_mat_rank(reduced.raw()) != rank) {
        return 1;
    }

    silex::flint::FmpzMat expected_reduced(rows, rank);
    fmpz_mat_mul(expected_reduced.raw(), input.raw(), transform.raw());
    return fmpz_mat_equal(expected_reduced.raw(), reduced.raw()) == 0 ? 1 : 0;
}

int test_flatter_wide_transform_boundary() {
    constexpr slong rows = 3;
    constexpr slong cols = 9;
    const slong input_entries[rows * cols] = {
            32373176721380998, 0, 0, -26155458657878188,
            -27972425631556, 976620610054834, 41529514241169880,
            25748372788561934, 7006,
            0, 32373176721380998, 0, 7703955618578414,
            -1359153410186158, -3363101440953616,
            -29118641599698938, 81558150271772994, 9420,
            0, 0, 32373176721380998, 10949198571930016,
            1275236133291490, 6292963271118118,
            56587654159065584, -52872796988158689, 11598,
    };

    silex::flint::FmpzMat input(rows, cols);
    for (slong row = 0; row < rows; ++row) {
        for (slong col = 0; col < cols; ++col) {
            set_entry_si(input.raw(), row, col,
                    input_entries[row * cols + col]);
        }
    }

    const slong rank = fmpz_mat_rank(input.raw());
    silex::flint::FmpzMat reduced(rows, rank);
    silex::flint::FmpzMat transform(cols, rank);
    const auto result =
            silex::lat::detail::flatter_column_lll_transform(
                    silex::flint::FmpzMatRef(reduced),
                    silex::flint::FmpzMatRef(transform),
                    silex::flint::FmpzMatConstRef(input), 1.02, 1);
    if (result.status == silex::lat::detail::FlatterBackendStatus::unavailable) {
        return 0;
    }
    if (result.status ==
        silex::lat::detail::FlatterBackendStatus::transform_unavailable) {
        return 0;
    }
    if (result.status != silex::lat::detail::FlatterBackendStatus::success ||
        result.rank != rank || fmpz_mat_rank(reduced.raw()) != rank) {
        return 1;
    }

    silex::flint::FmpzMat expected_reduced(rows, rank);
    fmpz_mat_mul(expected_reduced.raw(), input.raw(), transform.raw());
    return fmpz_mat_equal(expected_reduced.raw(), reduced.raw()) == 0 ? 1 : 0;
}

bool enum_count_for_basis(fmpz_mat_t basis,
        slong ambient_dim,
        slong bound_si,
    slong max_coord,
    EnumCounter& counter) {
    silex::lat::Lat lattice(ambient_dim);
    silex::flint::Arb bound;
    arb_set_si(bound.raw(), bound_si);

    const bool ok = lattice.set_basis(basis) &&
                    lattice.enum_short_vectors_arb(bound, max_coord, 128,
                            enum_count_callback, &counter);

    return ok;
}

int test_short_vector_enum() {
    fmpz_mat_t basis;
    EnumCounter counter;

    fmpz_mat_init(basis, 2, 2);

    set_entry_si(basis, 0, 0, 1);
    set_entry_si(basis, 1, 1, 1);
    counter = {};
    if (!enum_count_for_basis(basis, 2, 4, -1, counter) ||
        counter.count != 12 || counter.max_abs_coord != 2) {
        return 1;
    }
    counter = {};
    if (!enum_count_for_basis(basis, 2, 6, -1, counter) ||
        counter.count != 20 || counter.max_abs_coord != 2) {
        return 1;
    }
    counter = {};
    if (!enum_count_for_basis(basis, 2, 0, -1, counter) || counter.count != 0) {
        return 1;
    }
    counter = {};
    counter.abort_after = 3;
    if (!enum_count_for_basis(basis, 2, 100, -1, counter) ||
        counter.count != 3) {
        return 1;
    }
    counter = {};
    if (!enum_count_for_basis(basis, 2, 100, 1, counter) ||
        counter.count != 8 || counter.max_abs_coord != 1) {
        return 1;
    }

    fmpz_mat_clear(basis);
    fmpz_mat_init(basis, 3, 3);
    set_entry_si(basis, 0, 0, 1);
    set_entry_si(basis, 1, 1, 1);
    set_entry_si(basis, 2, 2, 1);
    counter = {};
    if (!enum_count_for_basis(basis, 3, 2, -1, counter) ||
        counter.count != 18 || counter.max_abs_coord != 1) {
        return 1;
    }

    fmpz_mat_clear(basis);
    fmpz_mat_init(basis, 4, 4);
    set_entry_si(basis, 0, 0, 1);
    set_entry_si(basis, 1, 1, 1);
    set_entry_si(basis, 2, 2, 1);
    set_entry_si(basis, 3, 3, 1);
    counter = {};
    if (!enum_count_for_basis(basis, 4, 4, -1, counter) ||
        counter.count != 88 || counter.max_abs_coord != 2) {
        return 1;
    }

    fmpz_mat_clear(basis);
    fmpz_mat_init(basis, 2, 2);
    set_entry_si(basis, 0, 0, 3);
    set_entry_si(basis, 0, 1, 0);
    set_entry_si(basis, 1, 0, 1);
    set_entry_si(basis, 1, 1, 2);
    counter = {};
    if (!enum_count_for_basis(basis, 2, 5, -1, counter) ||
        counter.count != 2) {
        return 1;
    }
    counter = {};
    if (!enum_count_for_basis(basis, 2, 6, -1, counter) ||
        counter.count != 2) {
        return 1;
    }

    silex::lat::Lat empty(2);
    silex::flint::Arb bound;
    arb_set_si(bound.raw(), -1);
    counter = {};
    if (!empty.enum_short_vectors_arb(bound, -1, 128, enum_count_callback,
                &counter) ||
        counter.count != 0) {
        return 1;
    }
    arb_set_si(bound.raw(), 1);
    if (empty.enum_short_vectors_arb(bound, -1, 1, enum_count_callback,
                &counter) ||
        empty.enum_short_vectors_arb(nullptr, -1, 128, enum_count_callback,
                &counter) ||
        empty.enum_short_vectors_arb(bound, -1, 128, nullptr, &counter)) {
        return 1;
    }

    fmpz_mat_clear(basis);
    return 0;
}

int test_lat_check() {
    silex::lat::Lat lattice(2);
    fmpz_mat_t basis;
    fmpz_mat_init(basis, 2, 2);
    set_entry_si(basis, 0, 0, 1);
    set_entry_si(basis, 1, 1, 1);
    lattice.set_basis(basis);

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    DebugCounter counter;
    const auto lattice_mask =
            silex::diagnostics_module_bit(silex::DiagnosticsModule::lattice);

    if (!lattice.check(nullptr) || !lattice.check(&diagnostics)) {
        return 1;
    }

    silex::diagnostics_set_verbose(diagnostics, silex::VerboseLevel::detail,
            lattice_mask, verbose_callback, &counter);
    if (!lattice.check(&diagnostics) || counter.verbose != 0) {
        return 1;
    }

    silex::diagnostics_set_debug_checks(diagnostics, silex::DebugLevel::normal,
            lattice_mask, debug_failure_callback, &counter);
    if (!lattice.check(&diagnostics) || counter.failures != 0 ||
        counter.verbose != 0) {
        return 1;
    }

    fmpz_mat_clear(basis);
    return 0;
}

int test_native_cpp_raii_call_sites() {
    static_assert(!std::is_copy_constructible_v<silex::lat::Lat>);
    static_assert(!std::is_copy_assignable_v<silex::lat::Lat>);

    silex::flint::FmpzMat basis(3, 2);
    silex::flint::FmpzMat got(3, 2);
    silex::flint::FmpzMat transform(3, 3);
    silex::flint::FmpzMat product(3, 2);
    silex::flint::FmpzMat z2_basis(2, 2);
    silex::flint::FmpzMat expected_basis(2, 2);
    silex::flint::Fmpz index;
    silex::flint::Fmpz p;
    silex::flint::Arb bound;
    EnumCounter counter;

    set_entry_si(basis.raw(), 0, 0, 2);
    set_entry_si(basis.raw(), 1, 1, 2);
    set_entry_si(basis.raw(), 2, 0, 1);
    set_entry_si(basis.raw(), 2, 1, 1);

    silex::lat::Lat lattice(2);
    silex::lat::Lat expected_from_basis(2);
    if (!lattice.set_basis(basis)) {
        return 1;
    }
    expected_from_basis.set_basis(basis);

    silex::lat::Lat moved(std::move(lattice));
    silex::lat::Lat assigned(0);
    assigned = std::move(moved);
    if (assigned.ambient_dim() != 2 || assigned.nrows() != 3 ||
        !assigned.get_basis(got) ||
        !lat_equal_hnf(assigned, expected_from_basis)) {
        return 1;
    }

    silex::lat::Lat hnf(2);
    if (!assigned.hnf_transform(hnf, transform)) {
        return 1;
    }
    fmpz_mat_mul(product.raw(), transform.raw(), basis.raw());
    if (fmpz_mat_equal(product.raw(), hnf.raw_basis()) == 0) {
        return 1;
    }

    set_entry_si(z2_basis.raw(), 0, 0, 1);
    set_entry_si(z2_basis.raw(), 1, 1, 1);
    silex::lat::Lat z2(2);
    if (!z2.set_basis(z2_basis) || !z2.index(index, assigned) ||
        !fmpz_equal_si(index.raw(), 2)) {
        return 1;
    }

    fmpz_set_si(p.raw(), 2);
    silex::lat::Lat saturated(2);
    if (!assigned.saturate(saturated, p) || !lat_equal_hnf(saturated, z2)) {
        return 1;
    }

    silex::lat::Lat reduced(2);
    if (!assigned.lll_reduce(reduced) || !lat_equal_hnf(reduced, assigned)) {
        return 1;
    }

    arb_set_si(bound.raw(), 2);
    counter = {};
    if (!z2.enum_short_vectors_arb(bound, 1, 128, enum_count_callback,
                &counter) ||
        counter.count != 8 || counter.max_abs_coord != 1) {
        return 1;
    }

    fmpz_mat_zero(expected_basis.raw());
    set_entry_si(expected_basis.raw(), 0, 0, 1);
    set_entry_si(expected_basis.raw(), 1, 1, 1);
    if (!z2.get_basis(expected_basis) ||
        fmpz_mat_equal(expected_basis.raw(), z2_basis.raw()) == 0) {
        return 1;
    }

    return 0;
}

}  // namespace

int main() {
    return test_init_set_swap_basis() != 0 || test_hnf_and_transform() != 0 ||
                   test_contains() != 0 || test_sum_intersection_index() != 0 ||
                   test_saturate() != 0 || test_lll_reduce() != 0 ||
                   test_fplll_row_transform_boundary() != 0 ||
                   test_fplll_column_image_transform_boundary() != 0 ||
                   test_fplll_bounded_bkz_row_transform_boundary() != 0 ||
                   test_flatter_full_rank_column_transform_boundary() != 0 ||
                   test_flatter_wide_transform_boundary() != 0 ||
                   test_short_vector_enum() != 0 || test_lat_check() != 0 ||
                   test_native_cpp_raii_call_sites() != 0
               ? 1
               : 0;
}
