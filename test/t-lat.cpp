#include <silex/lat.hpp>

#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>

#include "lat/fplll_backend_internal.hpp"
#include "lat/flatter_backend_internal.hpp"
#include "lll_reference.hpp"
#include "lat/lll_internal.hpp"

#include <algorithm>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

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

bool lll_basis_only_contract(const silex::flint::FmpzMat& input) {
    const slong columns = fmpz_mat_ncols(input.raw());
    silex::lat::Lat lattice(columns), normalized(columns), reduced(columns);
    if (!lattice.set_basis(input) || !lattice.hnf(normalized)) {
        return false;
    }
    const slong rank = normalized.nrows();
    silex::flint::FmpzMat legacy(rank, columns), tracked(rank, columns),
            untracked(rank, columns), zero(rank, rank), transform(rank, rank),
            product(rank, columns), canonical(rank, columns);
    fmpz_mat_set(legacy.raw(), normalized.raw_basis());
    fmpz_mat_set(tracked.raw(), normalized.raw_basis());
    fmpz_mat_set(untracked.raw(), normalized.raw_basis());
    silex::flint::FmpzLll config;
    if (rank != 0) {
        fmpz_mat_one(transform.raw());
        fmpz_lll(legacy.raw(), zero.raw(), config.raw());
        fmpz_lll(tracked.raw(), transform.raw(), config.raw());
        fmpz_lll(untracked.raw(), nullptr, config.raw());
        silex::flint::Fmpz determinant;
        fmpz_mat_det(determinant.raw(), transform.raw());
        fmpz_mat_mul(product.raw(), transform.raw(), normalized.raw_basis());
        if (!fmpz_mat_equal(product.raw(), tracked.raw()) ||
            !fmpz_is_pm1(determinant.raw()) ||
            !fmpz_mat_is_zero(zero.raw()) ||
            !silex::test::rational_lll_reduced(untracked) ||
            !fmpz_mat_is_reduced(untracked.raw(),
                    config.raw()->delta, config.raw()->eta)) {
            return false;
        }
    }
    fmpz_mat_hnf(canonical.raw(), untracked.raw());
    if (!fmpz_mat_equal(legacy.raw(), tracked.raw()) ||
        !fmpz_mat_equal(legacy.raw(), untracked.raw()) ||
        !fmpz_mat_equal(canonical.raw(), normalized.raw_basis()) ||
        !lattice.lll_reduce(reduced) || reduced.nrows() != rank ||
        reduced.ambient_dim() != columns ||
        !fmpz_mat_equal(reduced.raw_basis(), legacy.raw()) ||
        !lattice.lll_reduce(lattice) || lattice.nrows() != rank ||
        lattice.ambient_dim() != columns ||
        !fmpz_mat_equal(lattice.raw_basis(), legacy.raw())) {
        return false;
    }
    return true;
}

int test_lll_basis_only_contract() {
    for (slong columns : {0, 3}) {
        silex::flint::FmpzMat empty(0, columns), zero(4, columns);
        if (!lll_basis_only_contract(empty) || !lll_basis_only_contract(zero)) {
            return 1;
        }
    }
    silex::flint::FmpzMat dependent(4, 5);
    for (slong col = 0; col < 5; ++col) {
        set_entry_si(dependent.raw(), 0, col, (col % 2 ? -1 : 1) * (col + 2));
        fmpz_mul_si(fmpz_mat_entry(dependent.raw(), 1, col),
                fmpz_mat_entry(dependent.raw(), 0, col), -3);
    }
    if (!lll_basis_only_contract(dependent)) {
        return 1;
    }
    set_entry_si(dependent.raw(), 2, 3, 7);
    if (!lll_basis_only_contract(dependent)) {
        return 1;
    }

    // Exercise both sides of the upstream truncation threshold, on full-rank
    // HNF input so normalization cannot remove the large coefficients.
    for (slong bits : {16, 249, 250, 251, 512}) {
        for (slong dimension : {3, 4, 5, 14, 15}) {
            silex::flint::FmpzMat basis(dimension + 1, dimension);
            for (slong row = 0; row < dimension; ++row) {
                fmpz* diagonal = fmpz_mat_entry(basis.raw(), row, row);
                fmpz_one(diagonal);
                fmpz_mul_2exp(diagonal, diagonal, bits - 1);
                fmpz_add_ui(diagonal, diagonal, 101 + 2 * row);
                for (slong col = row + 1; col < dimension; ++col) {
                    fmpz* entry = fmpz_mat_entry(basis.raw(), row, col);
                    fmpz_tdiv_q_ui(entry, diagonal, row + 2);
                    fmpz_add_ui(entry, entry, (row + 1) * (col + 1));
                }
            }
            if (!lll_basis_only_contract(basis)) {
                return 1;
            }
        }
    }
    return 0;
}

int test_lll_certifier_boundaries() {
    silex::flint::FmpzMat basis(2, 2);
    silex::flint::Fmpz scale, cutoff, value;
    silex::flint::Fmpq parameter;
    silex::flint::FmpzLll config;
    fmpz_one(scale.raw());
    fmpz_mul_2exp(scale.raw(), scale.raw(), 60);
    auto agrees = [&](bool expected) {
        return silex::test::rational_lll_reduced(basis) == expected &&
               (fmpz_mat_is_reduced(basis.raw(), 0.99, 0.51) != 0) == expected &&
               (fmpz_lll_is_reduced(basis.raw(), config.raw(), 120) != 0) == expected;
    };
    silex::test::exact_binary_parameter(parameter, 0.51);
    fmpz_mul(cutoff.raw(), fmpq_numref(parameter.raw()), scale.raw());
    fmpz_divexact(cutoff.raw(), cutoff.raw(), fmpq_denref(parameter.raw()));
    for (slong offset : {-1, 0, 1}) {
        fmpz_mat_zero(basis.raw());
        fmpz_set(fmpz_mat_entry(basis.raw(), 0, 0), scale.raw());
        fmpz_set(fmpz_mat_entry(basis.raw(), 1, 1), scale.raw());
        fmpz_add_si(value.raw(), cutoff.raw(), offset);
        for (int sign : {-1, 1}) {
            fmpz_mul_si(fmpz_mat_entry(basis.raw(), 1, 0), value.raw(), sign);
            if (!agrees(offset <= 0)) return 1;
        }
    }
    silex::test::exact_binary_parameter(parameter, 0.99);
    fmpz_mul(cutoff.raw(), scale.raw(), scale.raw());
    fmpz_mul(cutoff.raw(), cutoff.raw(), fmpq_numref(parameter.raw()));
    fmpz_divexact(cutoff.raw(), cutoff.raw(), fmpq_denref(parameter.raw()));
    fmpz_sqrt(cutoff.raw(), cutoff.raw());
    for (slong offset : {0, 1}) {
        fmpz_mat_zero(basis.raw());
        fmpz_set(fmpz_mat_entry(basis.raw(), 0, 0), scale.raw());
        fmpz_add_ui(fmpz_mat_entry(basis.raw(), 1, 1), cutoff.raw(), offset);
        if (!agrees(offset == 1)) return 1;
    }
    // A nearly parallel input with a large cancellation in Gram-Schmidt.
    fmpz_set(fmpz_mat_entry(basis.raw(), 1, 0), scale.raw());
    fmpz_one(fmpz_mat_entry(basis.raw(), 1, 1));
    if (!agrees(false)) return 1;
    fmpz_lll(basis.raw(), nullptr, config.raw());
    return agrees(true) ? 0 : 1;
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

// Route selection in Lat::enum_short_vectors_arb is internal; these tests pick
// the route through its documented inputs.  The double route requires at most
// 32 rows and max_coord <= 10000 (or, with no cap, a coordinate bound within
// 10000); max_coord = 10001 or 33 rows forces the Arb route.
constexpr slong enum_double_route_cap = 10000;
constexpr slong enum_arb_route_cap = 10001;

struct EnumCollector {
    std::vector<std::string> vectors;
    slong abort_after = 0;
};

int enum_collect_callback(const fmpz_mat_t coeffs, void* user) {
    auto* collector = static_cast<EnumCollector*>(user);
    std::string key;
    for (slong j = 0; j < fmpz_mat_ncols(coeffs); ++j) {
        char* text = fmpz_get_str(nullptr, 10, fmpz_mat_entry(coeffs, 0, j));
        key += text;
        key += ',';
        flint_free(text);
    }
    collector->vectors.push_back(key);
    return collector->abort_after > 0 &&
                    static_cast<slong>(collector->vectors.size()) >=
                            collector->abort_after
            ? 0
            : 1;
}

bool enum_collect(const silex::lat::Lat& lattice,
        slong bound_si,
        slong max_coord,
        slong prec,
        EnumCollector& collector) {
    silex::flint::Arb bound;
    arb_set_si(bound.raw(), bound_si);
    return lattice.enum_short_vectors_arb(bound, max_coord, prec,
            enum_collect_callback, &collector);
}

bool enum_all_distinct(std::vector<std::string> vectors) {
    std::sort(vectors.begin(), vectors.end());
    return std::adjacent_find(vectors.begin(), vectors.end()) == vectors.end();
}

std::vector<std::string> enum_sorted(std::vector<std::string> vectors) {
    std::sort(vectors.begin(), vectors.end());
    return vectors;
}

// Exact reference: every nonzero coefficient row in [-radius, radius]^r whose
// exact squared norm is at most bound.  The caller picks radius at least
// the Fincke--Pohst coordinate bound (or the enumeration's max_coord) so the
// box holds every solution.
std::vector<std::string> enum_brute_force(const fmpz_mat_t basis,
        const fmpz_t bound,
        slong radius) {
    const slong rows = fmpz_mat_nrows(basis);
    const slong cols = fmpz_mat_ncols(basis);
    std::vector<slong> coeffs(static_cast<std::size_t>(rows), -radius);
    std::vector<std::string> out;
    silex::flint::Fmpz norm;
    silex::flint::Fmpz entry;
    for (;;) {
        bool nonzero = false;
        fmpz_zero(norm.raw());
        for (slong j = 0; j < cols; ++j) {
            fmpz_zero(entry.raw());
            for (slong i = 0; i < rows; ++i) {
                fmpz_addmul_si(entry.raw(), fmpz_mat_entry(basis, i, j),
                        coeffs[static_cast<std::size_t>(i)]);
            }
            fmpz_addmul(norm.raw(), entry.raw(), entry.raw());
        }
        std::string key;
        for (slong i = 0; i < rows; ++i) {
            nonzero = nonzero || coeffs[static_cast<std::size_t>(i)] != 0;
            key += std::to_string(coeffs[static_cast<std::size_t>(i)]);
            key += ',';
        }
        if (nonzero && fmpz_cmp(norm.raw(), bound) <= 0) {
            out.push_back(key);
        }

        slong position = 0;
        while (position < rows &&
                coeffs[static_cast<std::size_t>(position)] == radius) {
            coeffs[static_cast<std::size_t>(position)] = -radius;
            ++position;
        }
        if (position == rows) {
            break;
        }
        ++coeffs[static_cast<std::size_t>(position)];
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::vector<std::string> enum_brute_force(const fmpz_mat_t basis,
        slong bound_si,
        slong radius) {
    silex::flint::Fmpz bound;
    fmpz_set_si(bound.raw(), bound_si);
    return enum_brute_force(basis, bound.raw(), radius);
}

// Regression: a double-route range failure used to restart on the Arb route
// and deliver the same vectors again.
int test_short_vector_enum_no_duplicate_restart() {
    silex::flint::FmpzMat basis(3, 3);
    fmpz_one(fmpz_mat_entry(basis.raw(), 0, 0));
    fmpz_one(fmpz_mat_entry(basis.raw(), 1, 1));
    fmpz_one(fmpz_mat_entry(basis.raw(), 2, 2));
    fmpz_set_str(fmpz_mat_entry(basis.raw(), 2, 1), "10000000000000000000", 10);
    silex::lat::Lat lattice(3);
    if (!lattice.set_basis(basis)) {
        return 1;
    }

    const std::vector<std::string> capped = {
            "-1,0,0,", "0,-1,0,", "0,1,0,", "1,0,0,"};
    for (const slong max_coord :
            {slong(1), enum_double_route_cap, enum_arb_route_cap}) {
        EnumCollector collector;
        if (!enum_collect(lattice, 1, max_coord, 128, collector) ||
            enum_sorted(collector.vectors) != capped) {
            return 1;
        }
    }

    // A capped call whose uncapped coordinate interval is far outside the
    // slong range: the interval center (about 1e19 here) must be clamped in
    // double before any integer conversion.  The double
    // route's 2^-40 relative bound slack admits extra rows at this size, so
    // require a duplicate-free superset of the exact set.
    {
        silex::flint::Fmpz huge;
        fmpz_set_str(huge.raw(), "100000000000000000000000000000000000000", 10);
        silex::flint::Arb huge_bound;
        arb_set_fmpz(huge_bound.raw(), huge.raw());
        EnumCollector collector;
        if (!lattice.enum_short_vectors_arb(huge_bound, 1, 256,
                    enum_collect_callback, &collector) ||
            !enum_all_distinct(collector.vectors)) {
            return 1;
        }
        const std::vector<std::string> found = enum_sorted(collector.vectors);
        const std::vector<std::string> exact =
                enum_brute_force(basis.raw(), huge.raw(), 1);
        if (exact.empty() ||
            !std::includes(found.begin(), found.end(), exact.begin(),
                    exact.end())) {
            return 1;
        }
    }

    const std::vector<std::string> uncapped = {"-1,0,0,", "0,-1,0,",
            "0,-10000000000000000000,1,", "0,1,0,",
            "0,10000000000000000000,-1,", "1,0,0,"};
    EnumCollector collector;
    if (!enum_collect(lattice, 1, -1, 128, collector) ||
        enum_sorted(collector.vectors) != uncapped) {
        return 1;
    }
    return 0;
}

int test_short_vector_enum_arb_route() {
    // Vector set, including the capped boundary, on a skewed basis.
    silex::flint::FmpzMat skew(2, 2);
    set_entry_si(skew.raw(), 0, 0, 3);
    set_entry_si(skew.raw(), 1, 0, 1);
    set_entry_si(skew.raw(), 1, 1, 2);
    silex::lat::Lat skew_lattice(2);
    if (!skew_lattice.set_basis(skew)) {
        return 1;
    }
    EnumCollector collector;
    if (!enum_collect(skew_lattice, 9, enum_arb_route_cap, 128, collector) ||
        !enum_all_distinct(collector.vectors) ||
        enum_sorted(collector.vectors) != enum_brute_force(skew.raw(), 9, 4)) {
        return 1;
    }

    // A callback return of 0 aborts and the call still succeeds.
    silex::flint::FmpzMat z2(2, 2);
    set_entry_si(z2.raw(), 0, 0, 1);
    set_entry_si(z2.raw(), 1, 1, 1);
    silex::lat::Lat z2_lattice(2);
    if (!z2_lattice.set_basis(z2)) {
        return 1;
    }
    collector = {};
    collector.abort_after = 3;
    if (!enum_collect(z2_lattice, 100, enum_arb_route_cap, 128, collector) ||
        collector.vectors.size() != 3) {
        return 1;
    }

    // max_coord clamps coefficients on the Arb route.
    silex::flint::FmpzMat z1(1, 1);
    set_entry_si(z1.raw(), 0, 0, 1);
    silex::lat::Lat z1_lattice(1);
    if (!z1_lattice.set_basis(z1)) {
        return 1;
    }
    EnumCounter counter;
    silex::flint::Arb bound;
    arb_set_si(bound.raw(), 10005 * 10005);
    if (!z1_lattice.enum_short_vectors_arb(bound, enum_arb_route_cap, 128,
                enum_count_callback, &counter) ||
        counter.count != 2 * enum_arb_route_cap ||
        counter.max_abs_coord != enum_arb_route_cap) {
        return 1;
    }

    // More than 32 rows takes the Arb route with or without a cap.
    const slong rows = 33;
    silex::flint::FmpzMat identity(rows, rows);
    for (slong i = 0; i < rows; ++i) {
        fmpz_one(fmpz_mat_entry(identity.raw(), i, i));
    }
    silex::lat::Lat wide(rows);
    if (!wide.set_basis(identity)) {
        return 1;
    }
    for (const slong max_coord : {slong(-1), slong(1)}) {
        collector = {};
        if (!enum_collect(wide, 1, max_coord, 128, collector) ||
            collector.vectors.size() != static_cast<std::size_t>(2 * rows) ||
            !enum_all_distinct(collector.vectors)) {
            return 1;
        }
    }

    // Precision failure: at 8 bits the Arb route cannot decide the vectors
    // whose squared norm equals the bound and reports failure.  The double
    // route ignores prec after the Cholesky step and succeeds on the same
    // lattice, which also confirms the two caps select different routes.
    silex::flint::FmpzMat shear(2, 2);
    set_entry_si(shear.raw(), 0, 0, 1);
    set_entry_si(shear.raw(), 0, 1, 1);
    set_entry_si(shear.raw(), 1, 1, 1);
    silex::lat::Lat shear_lattice(2);
    if (!shear_lattice.set_basis(shear)) {
        return 1;
    }
    collector = {};
    if (enum_collect(shear_lattice, 2, enum_arb_route_cap, 8, collector)) {
        return 1;
    }
    collector = {};
    if (!enum_collect(shear_lattice, 2, enum_double_route_cap, 8, collector) ||
        enum_sorted(collector.vectors) != enum_brute_force(shear.raw(), 2, 3)) {
        return 1;
    }

    // Dependent rows make the Gram matrix singular: failure, no callbacks.
    silex::flint::FmpzMat dependent(2, 2);
    set_entry_si(dependent.raw(), 0, 0, 1);
    set_entry_si(dependent.raw(), 1, 0, 2);
    silex::lat::Lat dependent_lattice(2);
    if (!dependent_lattice.set_basis(dependent)) {
        return 1;
    }
    for (const slong max_coord : {enum_double_route_cap, enum_arb_route_cap}) {
        collector = {};
        if (enum_collect(dependent_lattice, 4, max_coord, 128, collector) ||
            !collector.vectors.empty()) {
            return 1;
        }
    }
    return 0;
}

// The double and Arb routes must report the same coefficient rows on inputs
// that both accept, and both must match an exact brute-force reference.
int test_short_vector_enum_cross_route() {
    struct Case {
        slong rows;
        slong entries[16];
        slong bound;
        slong radius;
    };
    const Case cases[] = {
            {2, {3, 0, 1, 2}, 20, 6},
            {3, {2, 1, 0, 1, 3, 1, 0, 1, 4}, 30, 8},
            {3, {1, 1, 1, 0, 2, 1, 1, 0, 3}, 12, 8},
            {4, {2, 0, 1, 0, 1, 3, 0, 1, 0, 1, 2, 1, 1, 0, 0, 3}, 16, 6},
    };
    for (const Case& test_case : cases) {
        silex::flint::FmpzMat basis(test_case.rows, test_case.rows);
        for (slong i = 0; i < test_case.rows; ++i) {
            for (slong j = 0; j < test_case.rows; ++j) {
                set_entry_si(basis.raw(), i, j,
                        test_case.entries[i * test_case.rows + j]);
            }
        }
        silex::lat::Lat lattice(test_case.rows);
        if (!lattice.set_basis(basis)) {
            return 1;
        }
        const std::vector<std::string> expected =
                enum_brute_force(basis.raw(), test_case.bound, test_case.radius);
        if (expected.empty()) {
            return 1;
        }
        for (const slong max_coord :
                {slong(-1), enum_double_route_cap, enum_arb_route_cap}) {
            EnumCollector collector;
            if (!enum_collect(lattice, test_case.bound, max_coord, 128,
                        collector) ||
                !enum_all_distinct(collector.vectors) ||
                enum_sorted(collector.vectors) != expected) {
                return 1;
            }
        }
    }
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

int test_lll_certified_routes() {
    using namespace silex::lat::detail;
    for (const slong rows : {1, 4, 14, 15}) {
        for (const slong cols : {rows, rows + 1}) {
            for (const slong bits : {16, 250, 251, 512}) {
                silex::flint::FmpzMat input(rows, cols), expected(rows, cols);
                silex::flint::FmpzMat reduced(rows, cols), transform(rows, rows);
                silex::flint::FmpzLll config;
                for (slong i = 0; i < rows; ++i) {
                    fmpz_one(fmpz_mat_entry(input.raw(), i, i));
                    fmpz_mul_2exp(fmpz_mat_entry(input.raw(), i, i),
                            fmpz_mat_entry(input.raw(), i, i), bits - 1);
                    if (i + 1 < cols) {
                        set_entry_si(input.raw(), i, i + 1, 7);
                    }
                }
                fmpz_mat_set(expected.raw(), input.raw());
                fmpz_lll(expected.raw(), transform.raw(), config.raw());
                const bool eligible = rows >= 2 && rows <= 14 &&
                        cols <= 14 && bits <= 250;
                for (const auto failure : {LllTestFailure::none,
                             LllTestFailure::reducer, LllTestFailure::certification}) {
                    const auto route = reduce_normalized_basis(reduced, input, failure);
                    const auto expected_route = !eligible ? LllRoute::ineligible
                            : failure == LllTestFailure::none ? LllRoute::certified
                                                             : LllRoute::fallback;
                    if (route != expected_route ||
                        !fmpz_mat_equal(reduced.raw(), expected.raw()) ||
                        !silex::test::rational_lll_reduced(reduced)) {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

int main() {
    return test_init_set_swap_basis() != 0 || test_hnf_and_transform() != 0 ||
                   test_contains() != 0 || test_sum_intersection_index() != 0 ||
                   test_saturate() != 0 || test_lll_reduce() != 0 ||
                   test_lll_basis_only_contract() != 0 ||
                   test_lll_certifier_boundaries() != 0 ||
                   test_lll_certified_routes() != 0 ||
                   test_fplll_row_transform_boundary() != 0 ||
                   test_fplll_column_image_transform_boundary() != 0 ||
                   test_fplll_bounded_bkz_row_transform_boundary() != 0 ||
                   test_flatter_full_rank_column_transform_boundary() != 0 ||
                   test_flatter_wide_transform_boundary() != 0 ||
                   test_short_vector_enum() != 0 ||
                   test_short_vector_enum_no_duplicate_restart() != 0 ||
                   test_short_vector_enum_arb_route() != 0 ||
                   test_short_vector_enum_cross_route() != 0 ||
                   test_lat_check() != 0 ||
                   test_native_cpp_raii_call_sites() != 0
               ? 1
               : 0;
}
