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
#include "test_hook_access.hpp"

#include <algorithm>
#include <cstdint>
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

// Deterministic pseudorandom source for the differential contains tests.
struct ContainsRng {
    std::uint64_t state;

    slong next(slong lo, slong hi) noexcept {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return lo + static_cast<slong>((state >> 33) %
                                       static_cast<std::uint64_t>(hi - lo + 1));
    }
};

slong hnf_nonzero_rows(const fmpz_mat_t hnf) noexcept {
    slong rows = fmpz_mat_nrows(hnf);
    while (rows > 0 && fmpz_mat_is_zero_row(hnf, rows - 1) != 0) {
        --rows;
    }
    return rows;
}

// The containment rule Lat::contains used before its echelon back-substitution
// path: v lies in the row lattice of B iff HNF(B) and HNF(B; v) agree after
// trailing zero rows are dropped.
bool reference_hnf_contains(const fmpz_mat_t basis, const fmpz* vector) {
    const slong rows = fmpz_mat_nrows(basis);
    const slong cols = fmpz_mat_ncols(basis);
    silex::flint::FmpzMat hnf(rows, cols);
    silex::flint::FmpzMat extended(rows + 1, cols);
    silex::flint::FmpzMat extended_hnf(rows + 1, cols);
    fmpz_mat_hnf(hnf.raw(), basis);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            fmpz_set(fmpz_mat_entry(extended.raw(), i, j),
                     fmpz_mat_entry(basis, i, j));
        }
    }
    for (slong j = 0; j < cols; ++j) {
        fmpz_set(fmpz_mat_entry(extended.raw(), rows, j), vector + j);
    }
    fmpz_mat_hnf(extended_hnf.raw(), extended.raw());
    const slong rank = hnf_nonzero_rows(hnf.raw());
    if (hnf_nonzero_rows(extended_hnf.raw()) != rank) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        for (slong j = 0; j < cols; ++j) {
            if (fmpz_equal(fmpz_mat_entry(hnf.raw(), i, j),
                           fmpz_mat_entry(extended_hnf.raw(), i, j)) == 0) {
                return false;
            }
        }
    }
    return true;
}

void random_lattice_basis(fmpz_mat_t basis, ContainsRng& rng, slong kind) {
    const slong rows = fmpz_mat_nrows(basis);
    const slong cols = fmpz_mat_ncols(basis);
    fmpz_mat_zero(basis);
    if (kind == 1) {
        // Rank-deficient: rows are small combinations of fewer generators.
        const slong rank = rng.next(0, rows < cols ? rows : cols);
        silex::flint::FmpzMat gens(rank, cols);
        for (slong i = 0; i < rank; ++i) {
            for (slong j = 0; j < cols; ++j) {
                fmpz_set_si(fmpz_mat_entry(gens.raw(), i, j), rng.next(-6, 6));
            }
        }
        for (slong i = 0; i < rows; ++i) {
            for (slong k = 0; k < rank; ++k) {
                silex::flint::Fmpz c;
                fmpz_set_si(c.raw(), rng.next(-3, 3));
                for (slong j = 0; j < cols; ++j) {
                    fmpz_addmul(fmpz_mat_entry(basis, i, j),
                                fmpz_mat_entry(gens.raw(), k, j), c.raw());
                }
            }
        }
        return;
    }
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            fmpz* entry = fmpz_mat_entry(basis, i, j);
            if (kind == 3 && rng.next(0, 2) != 0) {
                continue;
            }
            fmpz_set_si(entry, rng.next(-9, 9));
            if (kind == 2 && j % 2 == 0) {
                // Multi-limb entries on alternate columns.
                fmpz_mul_2exp(entry, entry, 80);
                fmpz_add_si(entry, entry, rng.next(-9, 9));
            }
        }
    }
}

int check_contains_against_reference(const silex::lat::Lat& lattice,
                                      const fmpz_mat_t basis,
                                      const fmpz* vector,
                                      slong& members,
                                      slong& nonmembers) {
    const bool expected = reference_hnf_contains(basis, vector);
    if (lattice.contains(vector) != expected) {
        return 1;
    }
    ++(expected ? members : nonmembers);
    return 0;
}

// Differential test of Lat::contains against the HNF-comparison rule over
// random full-rank, rank-deficient, multi-limb, sparse, empty and zero-width
// lattices, each stored three ways: the raw basis (not flagged HNF), its
// square-padded HNF with trailing zero rows, and the trimmed HNF from
// Lat::hnf.
int test_contains_matches_hnf_reference() {
    ContainsRng rng{0x5113c0de};
    slong members = 0;
    slong nonmembers = 0;
    slong raw_non_hnf = 0;
    for (slong trial = 0; trial < 600; ++trial) {
        const slong cols = rng.next(0, 6);
        const slong rows = rng.next(0, 8);
        const slong kind = rng.next(0, 3);
        silex::flint::FmpzMat basis(rows, cols);
        random_lattice_basis(basis.raw(), rng, kind);
        silex::flint::FmpzMat hnf(rows, cols);
        fmpz_mat_hnf(hnf.raw(), basis.raw());

        silex::lat::Lat raw(cols);
        silex::lat::Lat padded(cols);
        silex::lat::Lat trimmed(cols);
        if (!raw.set_basis(silex::flint::FmpzMatConstRef(basis)) ||
            !padded.set_basis(silex::flint::FmpzMatConstRef(hnf)) ||
            !raw.hnf(trimmed) || !padded.is_hnf() || !trimmed.is_hnf()) {
            return 1;
        }
        raw_non_hnf += raw.is_hnf() ? 0 : 1;

        silex::flint::FmpzMat vector(1, cols);
        // A 1 x 0 FLINT row may be a null pointer, which contains() rejects,
        // so zero-width lattices are probed through a separate (empty) array.
        silex::flint::Fmpz zero_width;
        fmpz* v = cols > 0 ? fmpz_mat_row(vector.raw(), 0) : zero_width.raw();
        silex::flint::Fmpz c;
        for (slong sample = 0; sample < 12; ++sample) {
            // Shapes: 0 lattice combination, 1 combination plus a unit
            // perturbation, 2 arbitrary vector, 3 combination divided by 2 or
            // 3 when exact (non-integral coordinates), and the zero vector.
            const slong shape = sample % 5;
            fmpz_mat_zero(vector.raw());
            if (shape == 0 || shape == 1 || shape == 3) {
                for (slong i = 0; i < rows; ++i) {
                    fmpz_set_si(c.raw(), rng.next(-5, 5));
                    for (slong j = 0; j < cols; ++j) {
                        fmpz_addmul(v + j, fmpz_mat_entry(basis.raw(), i, j),
                                    c.raw());
                    }
                }
            }
            if (shape == 1 && cols > 0) {
                fmpz_add_si(v + rng.next(0, cols - 1), v + rng.next(0, cols - 1),
                            rng.next(1, 3));
            }
            if (shape == 2) {
                for (slong j = 0; j < cols; ++j) {
                    fmpz_set_si(v + j, rng.next(-12, 12));
                }
            }
            if (shape == 3) {
                const slong divisor = rng.next(2, 3);
                bool exact = true;
                for (slong j = 0; j < cols; ++j) {
                    exact = exact && fmpz_divisible_si(v + j, divisor) != 0;
                }
                if (exact) {
                    for (slong j = 0; j < cols; ++j) {
                        fmpz_divexact_si(v + j, v + j, divisor);
                    }
                }
            }
            if (check_contains_against_reference(raw, basis.raw(), v,
                                                 members, nonmembers) != 0 ||
                check_contains_against_reference(padded, basis.raw(), v,
                                                 members, nonmembers) != 0 ||
                check_contains_against_reference(trimmed, basis.raw(), v,
                                                 members, nonmembers) != 0 ||
                (cols > 0 &&
                 !raw.contains_row(silex::flint::FmpzMatConstRef(vector), 0) !=
                         !raw.contains(v))) {
                return 1;
            }
        }
    }

    // A non-square echelon basis whose first pivot is not in column 0.
    silex::flint::FmpzMat echelon(2, 3);
    set_entry_si(echelon.raw(), 0, 1, 2);
    set_entry_si(echelon.raw(), 0, 2, 1);
    set_entry_si(echelon.raw(), 1, 2, 3);
    silex::lat::Lat shifted(3);
    if (!shifted.set_basis(silex::flint::FmpzMatConstRef(echelon)) ||
        !shifted.is_hnf()) {
        return 1;
    }
    const slong probes[6][3] = {{0, 4, 5}, {1, 0, 0}, {0, 2, 4},
                                {0, 1, 0}, {0, 0, 1}, {0, -2, 2}};
    silex::flint::FmpzMat probe(1, 3);
    for (const auto& entries : probes) {
        for (slong j = 0; j < 3; ++j) {
            set_entry_si(probe.raw(), 0, j, entries[j]);
        }
        slong unused = 0;
        if (check_contains_against_reference(shifted, echelon.raw(),
                                             fmpz_mat_row(probe.raw(), 0),
                                             unused, unused) != 0) {
            return 1;
        }
    }
    return (members > 0 && nonmembers > 0 && raw_non_hnf > 0) ? 0 : 1;
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

// Reference intersection through the Smith form route Lat::intersection used
// before it switched to one HNF transform: the right SNF transform columns
// past the rank span the integer kernel of [L^T | -R^T].
bool smith_intersection_reference(silex::flint::FmpzMat& out,
                                  const silex::lat::Lat& left,
                                  const silex::lat::Lat& right) {
    using silex::flint::FmpzMat;
    const slong n = left.ambient_dim();
    silex::lat::Lat left_hnf(n);
    silex::lat::Lat right_hnf(n);
    if (!left.hnf(left_hnf) || !right.hnf(right_hnf)) {
        return false;
    }
    const slong r = left_hnf.nrows();
    const slong s = right_hnf.nrows();
    if (r == 0 || s == 0) {
        out = FmpzMat(0, n);
        return true;
    }
    FmpzMat columns(n, r + s);
    FmpzMat snf(n, r + s);
    FmpzMat transform(r + s, r + s);
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(columns.raw(), j, i),
                     fmpz_mat_entry(left_hnf.raw_basis(), i, j));
        }
    }
    for (slong i = 0; i < s; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_neg(fmpz_mat_entry(columns.raw(), j, r + i),
                     fmpz_mat_entry(right_hnf.raw_basis(), i, j));
        }
    }
    fmpz_mat_snf_transform(snf.raw(), nullptr, transform.raw(), columns.raw());
    slong rank = 0;
    while (rank < n && rank < r + s &&
           !fmpz_is_zero(fmpz_mat_entry(snf.raw(), rank, rank))) {
        ++rank;
    }
    const slong nullity = r + s - rank;
    FmpzMat basis(nullity, n);
    for (slong k = 0; k < nullity; ++k) {
        for (slong i = 0; i < r; ++i) {
            for (slong j = 0; j < n; ++j) {
                fmpz_addmul(fmpz_mat_entry(basis.raw(), k, j),
                            fmpz_mat_entry(transform.raw(), i, rank + k),
                            fmpz_mat_entry(left_hnf.raw_basis(), i, j));
            }
        }
    }
    silex::lat::Lat tmp(n);
    if (!tmp.set_basis(basis) || !tmp.hnf(tmp)) {
        return false;
    }
    out = tmp.basis();
    return true;
}

// Lat::intersection must return exactly the HNF the Smith form route gave,
// on random lattices and on rank-deficient, full-rank, equal, nested and
// zero cases (trivial intersections arise only from random kind 0),
// including the aliased call.
int check_intersection_matches_smith(const silex::lat::Lat& left,
                                     const silex::lat::Lat& right) {
    silex::flint::FmpzMat expected(0, 0);
    if (!smith_intersection_reference(expected, left, right)) {
        return 1;
    }
    silex::lat::Lat result(left.ambient_dim());
    if (!left.intersection(result, right) ||
        fmpz_mat_nrows(result.raw_basis()) != fmpz_mat_nrows(expected.raw()) ||
        fmpz_mat_ncols(result.raw_basis()) != fmpz_mat_ncols(expected.raw()) ||
        !fmpz_mat_equal(result.raw_basis(), expected.raw())) {
        return 1;
    }
    silex::lat::Lat aliased(left.ambient_dim());
    if (!aliased.set_basis(left.basis_ref()) ||
        !aliased.intersection(aliased, right) ||
        fmpz_mat_nrows(aliased.raw_basis()) != fmpz_mat_nrows(expected.raw()) ||
        !fmpz_mat_equal(aliased.raw_basis(), expected.raw())) {
        return 1;
    }
    return 0;
}

int test_intersection_matches_smith_route() {
    using silex::flint::FmpzMat;
    flint_rand_t state;
    flint_rand_init(state);
    int failed = 0;
    slong nonzero_results = 0;
    slong zero_results = 0;
    for (slong iter = 0; iter < 400 && failed == 0; ++iter) {
        const slong n = 1 + static_cast<slong>(n_randint(state, 6));
        const slong kind = static_cast<slong>(n_randint(state, 6));
        const slong r = static_cast<slong>(n_randint(state, n + 3));
        FmpzMat left_basis(r, n);
        fmpz_mat_randtest(left_basis.raw(), state, 1 + n_randint(state, 12));
        silex::lat::Lat left(n);
        silex::lat::Lat right(n);
        if (!left.set_basis(left_basis)) {
            failed = 1;
            break;
        }
        if (kind == 0) {
            // Random, usually rank-deficient when s < n.
            const slong s = static_cast<slong>(n_randint(state, n + 3));
            FmpzMat right_basis(s, n);
            fmpz_mat_randtest(right_basis.raw(), state,
                              1 + n_randint(state, 12));
            failed = !right.set_basis(right_basis);
        } else if (kind == 1) {
            // Full rank: a nonsingular diagonal-dominant random basis.
            FmpzMat right_basis(n, n);
            fmpz_mat_randtest(right_basis.raw(), state, 6);
            for (slong i = 0; i < n; ++i) {
                fmpz_add_ui(fmpz_mat_entry(right_basis.raw(), i, i),
                            fmpz_mat_entry(right_basis.raw(), i, i), 1);
                fmpz_mul_ui(fmpz_mat_entry(right_basis.raw(), i, i),
                            fmpz_mat_entry(right_basis.raw(), i, i), 1000);
                fmpz_add_ui(fmpz_mat_entry(right_basis.raw(), i, i),
                            fmpz_mat_entry(right_basis.raw(), i, i), 1);
            }
            failed = !right.set_basis(right_basis);
        } else if (kind == 2) {
            // Equal lattices, given by different bases.
            FmpzMat mixed(r, n);
            fmpz_mat_set(mixed.raw(), left_basis.raw());
            if (r > 1) {
                fmpz_mat_t row_op;
                fmpz_mat_init(row_op, r, r);
                fmpz_mat_one(row_op);
                fmpz_set_si(fmpz_mat_entry(row_op, 0, r - 1),
                            static_cast<slong>(n_randint(state, 7)) - 3);
                fmpz_mat_mul(mixed.raw(), row_op, left_basis.raw());
                fmpz_mat_clear(row_op);
            }
            failed = !right.set_basis(mixed);
        } else if (kind == 3) {
            // Nested: an integer combination of the left rows.
            const slong s = static_cast<slong>(n_randint(state, r + 2));
            FmpzMat combination(s, r);
            fmpz_mat_randtest(combination.raw(), state, 4);
            FmpzMat right_basis(s, n);
            if (r > 0) {
                fmpz_mat_mul(right_basis.raw(), combination.raw(),
                             left_basis.raw());
            }
            failed = !right.set_basis(right_basis);
        } else if (kind == 4) {
            // Zero lattice: no rows, or only zero rows.
            FmpzMat right_basis(static_cast<slong>(n_randint(state, 3)), n);
            failed = !right.set_basis(right_basis);
        } else {
            // A scaled copy of the left lattice, which it contains.
            FmpzMat right_basis(r, n);
            fmpz_mat_scalar_mul_si(right_basis.raw(), left_basis.raw(),
                                   2 + static_cast<slong>(n_randint(state, 5)));
            failed = !right.set_basis(right_basis);
        }
        if (failed != 0) {
            break;
        }
        failed = check_intersection_matches_smith(left, right) != 0 ||
                 check_intersection_matches_smith(right, left) != 0;
        silex::lat::Lat result(n);
        if (failed == 0 && left.intersection(result, right)) {
            if (result.nrows() > 0) {
                ++nonzero_results;
            } else {
                ++zero_results;
            }
        }
    }
    flint_rand_clear(state);
    return (failed == 0 && nonzero_results > 100 && zero_results > 10) ? 0 : 1;
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
    // The Gram matrix certificate on B B^T (entries of about 120 bits, so
    // the integer test runs) against FLINT's exact fmpz_mat_is_reduced_gram.
    silex::flint::FmpzMat gram(2, 2);
    silex::flint::FmpzLll gram_config(0.99, 0.51, GRAM, APPROX);
    auto agrees = [&](bool expected) {
        fmpz_mat_gram(gram.raw(), basis.raw());
        return silex::test::rational_lll_reduced(basis) == expected &&
               (fmpz_mat_is_reduced(basis.raw(), 0.99, 0.51) != 0) == expected &&
               silex::lat::detail::lll_certificate_holds(basis,
                       silex::flint::FmpzLllConstRef(config)) == expected &&
               (fmpz_lll_is_reduced(basis.raw(), config.raw(), 120) != 0) == expected &&
               (fmpz_mat_is_reduced_gram(gram.raw(), 0.99, 0.51) != 0) == expected &&
               silex::lat::detail::integer_gram_is_lll_reduced(gram, 0.99,
                       0.51) == expected &&
               silex::lat::detail::lll_certificate_holds(gram,
                       silex::flint::FmpzLllConstRef(gram_config)) == expected;
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
    // double before any integer conversion.  The double route's 2^-40
    // relative bound slack admits extra rows at this size, so require a
    // duplicate-free superset of the exact set.
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

// Lat::lll_reduce's certified first phase returns fmpz_lll's basis on
// every route.  The 16-bit rows take the integer certificate and the
// 250-bit rows fmpz_mat_is_reduced (certified_lll_integer_test_max_bits).
int test_lll_certified_routes() {
    using namespace silex::lat::detail;
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::TestHookAccess;
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
                for (const auto failure : {CertifiedLllTestFailure::none,
                             CertifiedLllTestFailure::first_phase,
                             CertifiedLllTestFailure::certificate}) {
                    const slong injected =
                            TestHookAccess::certified_lll_failures_injected();
                    TestHookAccess::set_certified_lll_failure(failure);
                    const auto route = reduce_normalized_basis(reduced, input);
                    TestHookAccess::set_certified_lll_failure(
                            CertifiedLllTestFailure::none);
                    const auto expected_route =
                            !eligible ? CertifiedLllRoute::ineligible
                            : failure == CertifiedLllTestFailure::none
                                    ? CertifiedLllRoute::certified
                                    : CertifiedLllRoute::fallback;
                    const slong expected_injected =
                            eligible && failure != CertifiedLllTestFailure::none
                                    ? 1 : 0;
                    if (route != expected_route ||
                        TestHookAccess::certified_lll_failures_injected() !=
                                injected + expected_injected ||
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

// certified_fmpz_lll and lll_first_phase_certified against plain fmpz_lll
// on random bases of 4 to 280 bits, with a null, zero and identity
// transform (the forms the call sites pass) and every injected failure.
// The matrix and transform must equal fmpz_lll's whenever the helper
// reports them as final.  The run must certify phase outputs on both sides
// of the 16-bit certificate cutoff.
int test_certified_fmpz_lll_matches_fmpz_lll() {
    using namespace silex::lat::detail;
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::TestHookAccess;
    using silex::flint::FmpzMat;
    silex::flint::FmpzLll config;
    const silex::flint::FmpzLllConstRef config_ref(config);
    flint_rand_t state;
    flint_rand_init(state);
    slong certified_small = 0;
    slong certified_large = 0;
    slong fallback = 0;
    slong ineligible = 0;
    const flint_bitcnt_t bit_choices[] = {4, 12, 20, 40, 120, 245, 280};
    for (slong iter = 0; iter < 700; ++iter) {
        const slong rows = 1 + static_cast<slong>(n_randint(state, 10));
        const slong cols = rows + static_cast<slong>(n_randint(state, 6));
        const flint_bitcnt_t bits = bit_choices[n_randint(state, 7)];
        FmpzMat input(rows, cols);
        fmpz_mat_randtest(input.raw(), state, bits);
        if (n_randint(state, 2) == 0) {
            FmpzMat scratch(rows, rows);
            fmpz_lll(input.raw(), scratch.raw(), config.raw());
            if (n_randint(state, 2) == 0 && rows > 1) {
                // Break the reduction so the phase has work to do.
                for (slong j = 0; j < cols; ++j) {
                    fmpz_addmul_ui(fmpz_mat_entry(input.raw(), rows - 1, j),
                            fmpz_mat_entry(input.raw(), 0, j), 3);
                }
            }
        }
        const auto failure = static_cast<CertifiedLllTestFailure>(
                n_randint(state, 3));
        for (int form = 0; form < 3; ++form) {
            FmpzMat start_transform(rows, rows);
            if (form == 2) {
                fmpz_mat_one(start_transform.raw());
            }
            fmpz_mat_struct* const no_transform = nullptr;
            FmpzMat expected(rows, cols), expected_transform(rows, rows);
            fmpz_mat_set(expected.raw(), input.raw());
            fmpz_mat_set(expected_transform.raw(), start_transform.raw());
            fmpz_lll(expected.raw(),
                     form == 0 ? no_transform : expected_transform.raw(),
                     config.raw());

            FmpzMat actual(rows, cols), actual_transform(rows, rows);
            fmpz_mat_set(actual.raw(), input.raw());
            fmpz_mat_set(actual_transform.raw(), start_transform.raw());
            TestHookAccess::set_certified_lll_failure(failure);
            const CertifiedLllRoute route = certified_fmpz_lll(actual.raw(),
                    form == 0 ? no_transform : actual_transform.raw(),
                    config_ref);
            TestHookAccess::set_certified_lll_failure(
                    CertifiedLllTestFailure::none);
            if (!fmpz_mat_equal(actual.raw(), expected.raw()) ||
                (form != 0 && !fmpz_mat_equal(actual_transform.raw(),
                                              expected_transform.raw()))) {
                flint_rand_clear(state);
                return 1;
            }

            // The phase-only entry point, as the saturation kernel calls it
            // (null transform): certified means final.
            FmpzMat phase(rows, cols);
            fmpz_mat_set(phase.raw(), input.raw());
            TestHookAccess::set_certified_lll_failure(failure);
            const CertifiedLllRoute phase_route = lll_first_phase_certified(
                    phase.raw(), no_transform, config_ref);
            TestHookAccess::set_certified_lll_failure(
                    CertifiedLllTestFailure::none);
            if (phase_route != route ||
                (phase_route == CertifiedLllRoute::certified &&
                 !fmpz_mat_equal(phase.raw(), expected.raw())) ||
                (phase_route == CertifiedLllRoute::ineligible &&
                 !fmpz_mat_equal(phase.raw(), input.raw()))) {
                flint_rand_clear(state);
                return 1;
            }
            if (route == CertifiedLllRoute::certified) {
                const slong out_bits = fmpz_mat_max_bits(actual.raw());
                ((out_bits < 0 ? -out_bits : out_bits) <=
                                 certified_lll_integer_test_max_bits
                         ? certified_small
                         : certified_large) += 1;
            } else if (route == CertifiedLllRoute::fallback) {
                ++fallback;
            } else {
                ++ineligible;
            }
        }
    }
    flint_rand_clear(state);
    return certified_small > 0 && certified_large > 0 && fallback > 0 &&
                   ineligible > 0
           ? 0
           : 1;
}

// integer_gram_is_lll_reduced returns false on a singular or an indefinite
// symmetric matrix, so it certifies only positive definite Gram matrices.
// [[0, 0], [0, 1]] (d_1 = 0) and, with delta = 0.25 < eta^2 (outside
// FLINT's range), [[4, 2], [2, 1]] (d_2 = 0) satisfy the size and Lovasz
// inequalities and are rejected only by the positivity test of the leading
// minors.
int test_integer_gram_rejects_non_positive_definite() {
    using silex::flint::FmpzMat;
    struct Case {
        slong entries[4];
        double delta;
    };
    const Case cases[] = {
            {{1, 0, 0, 0}, 0.99},  {{1, 0, 0, -1}, 0.99},
            {{0, 0, 0, 1}, 0.99},  {{-1, 0, 0, 1}, 0.99},
            {{4, 2, 2, 1}, 0.99},  {{4, 2, 2, 1}, 0.25},
            {{4, 2, 2, 0}, 0.25},
    };
    for (const Case& c : cases) {
        FmpzMat gram(2, 2);
        for (slong k = 0; k < 4; ++k) {
            fmpz_set_si(fmpz_mat_entry(gram.raw(), k / 2, k % 2), c.entries[k]);
        }
        if (silex::lat::detail::integer_gram_is_lll_reduced(gram, c.delta,
                                                            0.51)) {
            return 1;
        }
    }
    return 0;
}
// gram_interval_screen_proves_lll_reduced is a sufficient test: wherever
// it holds, integer_gram_is_lll_reduced and FLINT's exact
// fmpz_mat_is_reduced_gram hold.  [[12, 4], [4, 9]] at delta = 0.75 has
// mu = 1/3 and Lovasz equality (r_11 = 23/3 = (delta - 1/9) 12): the
// balls cannot prove it, and the exact test accepts it.  [[12, 4], [4, 8]]
// breaks the Lovasz test, [[4, 3], [3, 4]] breaks only the size condition,
// and the non positive definite matrices of the test above are not
// accepted.  On random fmpz_lll outputs of 2 to 40 rows (4 to 120 bits) the
// screen must decide.  On unreduced or perturbed matrices soundness is
// checked against FLINT's fmpz_mat_is_reduced_gram (integer_gram_is_lll_reduced
// calls the screen itself, so it is no independent check there).
int test_gram_interval_screen_is_sufficient() {
    using namespace silex::lat::detail;
    using silex::flint::FmpzMat;
    auto gram2 = [](slong a, slong b, slong c) {
        FmpzMat gram(2, 2);
        fmpz_set_si(fmpz_mat_entry(gram.raw(), 0, 0), a);
        fmpz_set_si(fmpz_mat_entry(gram.raw(), 1, 0), b);
        fmpz_set_si(fmpz_mat_entry(gram.raw(), 0, 1), b);
        fmpz_set_si(fmpz_mat_entry(gram.raw(), 1, 1), c);
        return gram;
    };
    {
        const FmpzMat equality = gram2(12, 4, 9);
        if (gram_interval_screen_proves_lll_reduced(equality, 0.75, 0.51) ||
            !integer_gram_is_lll_reduced(equality, 0.75, 0.51) ||
            fmpz_mat_is_reduced_gram(equality.raw(), 0.75, 0.51) == 0) {
            return 1;
        }
        const FmpzMat broken = gram2(12, 4, 8);
        if (gram_interval_screen_proves_lll_reduced(broken, 0.75, 0.51) ||
            integer_gram_is_lll_reduced(broken, 0.75, 0.51)) {
            return 1;
        }
        // mu = 3/4 > eta while the Lovasz test holds (1.75 >= 1.71): only
        // the size condition fails.
        const FmpzMat oversize = gram2(4, 3, 4);
        if (gram_interval_screen_proves_lll_reduced(oversize, 0.99, 0.51) ||
            integer_gram_is_lll_reduced(oversize, 0.99, 0.51)) {
            return 1;
        }
        const FmpzMat reduced = gram2(4, 1, 4);
        if (!gram_interval_screen_proves_lll_reduced(reduced, 0.99, 0.51) ||
            !integer_gram_is_lll_reduced(reduced, 0.99, 0.51)) {
            return 1;
        }
        const slong non_definite[][3] = {
                {1, 0, 0}, {1, 0, -1}, {0, 0, 1}, {-1, 0, 1}, {4, 2, 1}};
        for (const auto& c : non_definite) {
            if (gram_interval_screen_proves_lll_reduced(gram2(c[0], c[1], c[2]),
                                                        0.25, 0.51)) {
                return 1;
            }
        }
        FmpzMat single(1, 1);
        fmpz_one(fmpz_mat_entry(single.raw(), 0, 0));
        FmpzMat wide(2, 3);
        if (gram_interval_screen_proves_lll_reduced(single, 0.99, 0.51) ||
            gram_interval_screen_proves_lll_reduced(wide, 0.99, 0.51)) {
            return 1;
        }
    }
    silex::flint::FmpzLll basis_config;
    flint_rand_t state;
    flint_rand_init(state);
    const flint_bitcnt_t bit_choices[] = {4, 20, 60, 120};
    slong decided = 0;
    slong rejected = 0;
    for (slong iter = 0; iter < 200; ++iter) {
        const slong rows = 2 + static_cast<slong>(n_randint(state, 39));
        const slong cols = rows + static_cast<slong>(n_randint(state, 4));
        FmpzMat basis(rows, cols);
        do {
            fmpz_mat_randtest(basis.raw(), state, bit_choices[iter % 4]);
        } while (fmpz_mat_rank(basis.raw()) < rows);
        const bool lll = iter % 3 != 0;
        if (lll) {
            fmpz_lll(basis.raw(), nullptr, basis_config.raw());
            if (iter % 3 == 2) {
                // Perturb one entry of the reduced basis by +-1.
                fmpz* entry = fmpz_mat_entry(
                        basis.raw(),
                        static_cast<slong>(
                                n_randint(state, static_cast<ulong>(rows))),
                        static_cast<slong>(
                                n_randint(state, static_cast<ulong>(cols))));
                fmpz_add_si(entry, entry, n_randint(state, 2) == 0 ? -1 : 1);
            }
        }
        FmpzMat gram(rows, rows);
        fmpz_mat_gram(gram.raw(), basis.raw());
        const bool screen =
                gram_interval_screen_proves_lll_reduced(gram, 0.99, 0.51);
        const bool exact = integer_gram_is_lll_reduced(gram, 0.99, 0.51);
        if (screen &&
            fmpz_mat_is_reduced_gram(gram.raw(), 0.99, 0.51) == 0) {
            flint_rand_clear(state);
            return 1;
        }
        if (lll && iter % 3 == 1 && !screen) {
            flint_rand_clear(state);
            return 1;
        }
        decided += screen ? 1 : 0;
        rejected += exact ? 0 : 1;
    }
    flint_rand_clear(state);
    return decided > 0 && rejected > 0 ? 0 : 1;
}

// certified_fmpz_lll on GRAM input (the Gram matrices B B^T of random
// full-rank bases of 4 to 200 bits, 2 to 34 rows) against plain fmpz_lll,
// with a null and an identity transform and every injected failure.  The
// matrix and transform must equal fmpz_lll's on every route.  The run must
// certify phase outputs with the integer Gram test and, above
// certified_lll_gram_integer_test_max_size, with FLINT's own check, and
// must reach the fallback.  Every certified-or-not decision of the integer
// Gram test is compared with fmpz_mat_is_reduced_gram and with
// integer_basis_is_lll_reduced on B, and a single row, a non square matrix
// and a GRAM EXACT context are ineligible.
int test_certified_fmpz_lll_gram_matches_fmpz_lll() {
    using namespace silex::lat::detail;
    using silex::detail::CertifiedLllTestFailure;
    using silex::detail::TestHookAccess;
    using silex::flint::FmpzMat;
    silex::flint::FmpzLll basis_config;
    silex::flint::FmpzLll config(0.99, 0.51, GRAM, APPROX);
    const silex::flint::FmpzLllConstRef config_ref(config);
    flint_rand_t state;
    flint_rand_init(state);
    slong certified_integer = 0;
    slong certified_flint = 0;
    slong fallback = 0;
    slong large_rejected = 0;
    const flint_bitcnt_t bit_choices[] = {4, 20, 60, 120, 200};
    for (slong iter = 0; iter < 260; ++iter) {
        const bool large = iter % 20 == 0;
        const slong rows = large ? 30 + static_cast<slong>(n_randint(state, 5))
                                 : 2 + static_cast<slong>(n_randint(state, 11));
        const slong cols = rows + static_cast<slong>(n_randint(state, 4));
        const flint_bitcnt_t bits =
                large ? 120 : bit_choices[n_randint(state, 5)];
        FmpzMat basis(rows, cols);
        do {
            fmpz_mat_randtest(basis.raw(), state, bits);
        } while (fmpz_mat_rank(basis.raw()) < rows);
        if (large || n_randint(state, 2) == 0) {
            fmpz_lll(basis.raw(), nullptr, basis_config.raw());
            if (n_randint(state, 2) == 0) {
                // Break the reduction so the phase has work to do.
                for (slong j = 0; j < cols; ++j) {
                    fmpz_addmul_ui(fmpz_mat_entry(basis.raw(), rows - 1, j),
                            fmpz_mat_entry(basis.raw(), 0, j), 3);
                }
            }
        }
        FmpzMat input(rows, rows);
        fmpz_mat_gram(input.raw(), basis.raw());

        // The integer Gram test against FLINT's exact Gram test and the
        // basis test, on the input and on the reduced input.
        FmpzMat reduced_basis(rows, cols);
        fmpz_mat_set(reduced_basis.raw(), basis.raw());
        fmpz_lll(reduced_basis.raw(), nullptr, basis_config.raw());
        for (const FmpzMat* pair : {&basis, &reduced_basis}) {
            FmpzMat gram(rows, rows);
            fmpz_mat_gram(gram.raw(), pair->raw());
            const bool integer = integer_gram_is_lll_reduced(gram, 0.99, 0.51);
            // The certificate on the Gram matrix: the integer test up to
            // the cutoff, FLINT's own check above it; on these positive
            // definite inputs both decide the exact predicate.
            const bool certificate = lll_certificate_holds(gram, config_ref);
            if (integer !=
                        (fmpz_mat_is_reduced_gram(gram.raw(), 0.99, 0.51) != 0) ||
                integer != integer_basis_is_lll_reduced(*pair, 0.99, 0.51) ||
                certificate != integer ||
                certificate !=
                        (fmpz_lll_is_reduced(gram.raw(), config.raw(), 120) != 0) ||
                (pair == &reduced_basis && !certificate)) {
                flint_rand_clear(state);
                return 1;
            }
            const slong gram_bits = fmpz_mat_max_bits(gram.raw());
            if (large && pair == &basis &&
                rows * (gram_bits < 0 ? -gram_bits : gram_bits) >
                        certified_lll_gram_integer_test_max_size) {
                // Above the cutoff, FLINT's check must also be seen to
                // reject a matrix the phase has not reduced.
                large_rejected += certificate ? 0 : 1;
            }
        }

        const auto failure = static_cast<CertifiedLllTestFailure>(
                n_randint(state, 3));
        for (int form = 0; form < 2; ++form) {
            fmpz_mat_struct* const no_transform = nullptr;
            FmpzMat expected(rows, rows), expected_transform(rows, rows);
            fmpz_mat_set(expected.raw(), input.raw());
            fmpz_mat_one(expected_transform.raw());
            fmpz_lll(expected.raw(),
                     form == 0 ? no_transform : expected_transform.raw(),
                     config.raw());

            FmpzMat actual(rows, rows), actual_transform(rows, rows);
            fmpz_mat_set(actual.raw(), input.raw());
            fmpz_mat_one(actual_transform.raw());
            TestHookAccess::set_certified_lll_failure(failure);
            const CertifiedLllRoute route = certified_fmpz_lll(actual.raw(),
                    form == 0 ? no_transform : actual_transform.raw(),
                    config_ref);
            TestHookAccess::set_certified_lll_failure(
                    CertifiedLllTestFailure::none);
            if (!fmpz_mat_equal(actual.raw(), expected.raw()) ||
                (form != 0 && !fmpz_mat_equal(actual_transform.raw(),
                                              expected_transform.raw())) ||
                route == CertifiedLllRoute::ineligible ||
                (failure != CertifiedLllTestFailure::none &&
                 route != CertifiedLllRoute::fallback)) {
                flint_rand_clear(state);
                return 1;
            }
            if (route == CertifiedLllRoute::certified) {
                const slong out_bits = fmpz_mat_max_bits(actual.raw());
                (rows * (out_bits < 0 ? -out_bits : out_bits) <=
                                 certified_lll_gram_integer_test_max_size
                         ? certified_integer
                         : certified_flint) += 1;
            } else {
                ++fallback;
            }
        }
    }
    flint_rand_clear(state);
    if (certified_integer == 0 || certified_flint == 0 || fallback == 0 ||
        large_rejected == 0) {
        return 1;
    }

    // Ineligible inputs go to fmpz_lll unchanged.
    const auto ineligible = [](FmpzMat& matrix,
                               const silex::flint::FmpzLll& context) {
        const slong rows = fmpz_mat_nrows(matrix.raw());
        FmpzMat expected(rows, fmpz_mat_ncols(matrix.raw()));
        FmpzMat expected_transform(rows, rows);
        fmpz_mat_set(expected.raw(), matrix.raw());
        fmpz_mat_one(expected_transform.raw());
        fmpz_lll(expected.raw(), expected_transform.raw(), context.raw());
        FmpzMat transform(rows, rows);
        fmpz_mat_one(transform.raw());
        return certified_fmpz_lll(matrix.raw(), transform.raw(),
                       silex::flint::FmpzLllConstRef(context)) ==
                       CertifiedLllRoute::ineligible &&
               fmpz_mat_equal(matrix.raw(), expected.raw()) &&
               fmpz_mat_equal(transform.raw(), expected_transform.raw());
    };
    FmpzMat single(1, 1);
    fmpz_set_si(fmpz_mat_entry(single.raw(), 0, 0), 9);
    FmpzMat gram(3, 3);
    for (slong i = 0; i < 3; ++i) {
        for (slong j = 0; j < 3; ++j) {
            fmpz_set_si(fmpz_mat_entry(gram.raw(), i, j),
                        i == j ? 10 + i : 3);
        }
    }
    silex::flint::FmpzLll exact_config(0.99, 0.51, GRAM, EXACT);
    if (!ineligible(single, config) || !ineligible(gram, exact_config)) {
        return 1;
    }
    FmpzMat wide(2, 3);
    fmpz_set_si(fmpz_mat_entry(wide.raw(), 0, 0), 5);
    fmpz_set_si(fmpz_mat_entry(wide.raw(), 1, 1), 7);
    FmpzMat wide_transform(2, 2);
    return lll_first_phase_certified(wide.raw(), wide_transform.raw(),
                                     config_ref) ==
                           CertifiedLllRoute::ineligible
                   ? 0
                   : 1;
}

int main() {
    return test_init_set_swap_basis() != 0 || test_hnf_and_transform() != 0 ||
                   test_contains() != 0 ||
                   test_contains_matches_hnf_reference() != 0 ||
                   test_sum_intersection_index() != 0 ||
                   test_intersection_matches_smith_route() != 0 ||
                   test_saturate() != 0 || test_lll_reduce() != 0 ||
                   test_lll_basis_only_contract() != 0 ||
                   test_lll_certifier_boundaries() != 0 ||
                   test_lll_certified_routes() != 0 ||
                   test_certified_fmpz_lll_matches_fmpz_lll() != 0 ||
                   test_integer_gram_rejects_non_positive_definite() != 0 ||
                   test_gram_interval_screen_is_sufficient() != 0 ||
                   test_certified_fmpz_lll_gram_matches_fmpz_lll() != 0 ||
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
