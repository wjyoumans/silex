#include <silex/diagnostics.hpp>
#include <silex/fmpz_smat.hpp>

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/nmod_mat.hpp>

#include "fmpz_smat/fmpz_smat_internal.hpp"

#include <flint/fmpz_mat.h>

#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace {

void set_entry_si(fmpz_mat_t matrix, slong row, slong col, slong value) noexcept {
    fmpz_set_si(fmpz_mat_entry(matrix, row, col), value);
}

void set_entry_si(silex::flint::FmpzMat& matrix,
                  slong row,
                  slong col,
                  slong value) noexcept {
    silex::flint::fmpz_set_si(
            silex::flint::fmpz_mat_entry(matrix, row, col), value);
}

bool entry_is_si(const fmpz_mat_t matrix, slong row, slong col, slong value) noexcept {
    return fmpz_equal_si(fmpz_mat_entry(matrix, row, col), value) != 0;
}

bool fmpz_is_si(const fmpz_t value, slong expected) noexcept {
    return fmpz_equal_si(value, expected) != 0;
}

slong trimmed_nonzero_rows(silex::flint::FmpzMatConstRef matrix) noexcept {
    slong rows = silex::flint::fmpz_mat_nrows(matrix);
    while (rows > 0 && ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    return rows;
}

bool hnf_context_matches_prefix(silex::fmpz_smat::HnfContext& context,
                                silex::flint::FmpzMatConstRef rows,
                                slong used_rows) noexcept {
    const slong nrows = silex::flint::fmpz_mat_nrows(rows);
    const slong ncols = silex::flint::fmpz_mat_ncols(rows);
    if (used_rows < 0 || used_rows > nrows) {
        return false;
    }

    silex::flint::FmpzMatConstWindow prefix(rows, 0, 0, used_rows, ncols);
    silex::flint::FmpzMat hnf_full(used_rows, ncols);
    ::fmpz_mat_hnf(hnf_full.raw(), prefix.raw());

    const slong hnf_rows =
            trimmed_nonzero_rows(silex::flint::FmpzMatConstRef(hnf_full));
    bool ok = hnf_rows == context.rank();
    if (ok) {
        silex::flint::FmpzMatConstWindow hnf_window(
                silex::flint::FmpzMatConstRef(hnf_full), 0, 0, hnf_rows,
                ncols);
        silex::flint::FmpzMat expected(hnf_rows, ncols);
        silex::flint::fmpz_mat_set(
                silex::flint::FmpzMatRef(expected), hnf_window.const_ref());

        silex::flint::FmpzMat got(hnf_rows, ncols);
        ok = context.get_hnf(silex::flint::FmpzMatRef(got)) &&
             silex::flint::fmpz_mat_equal(
                     silex::flint::FmpzMatConstRef(got),
                     silex::flint::FmpzMatConstRef(expected));
    }

    return ok;
}

int test_sparse_row_basic() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::SparseRow copy;
    silex::fmpz_smat::SparseRow other;
    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);

    if (row.length() != 0) {
        return 1;
    }
    row.get_entry(x, 5);
    if (!fmpz_is_zero(x)) {
        return 1;
    }

    fmpz_set_si(x, 4);
    if (!row.set_entry(3, x)) {
        return 1;
    }
    fmpz_set_si(x, -2);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    if (row.length() != 2 || row.col(0) != 1 || row.col(1) != 3) {
        return 1;
    }
    std::span<const slong> columns = row.columns();
    silex::flint::FmpzVecConstRef values = row.values();
    if (columns.size() != 2 || columns[0] != 1 || columns[1] != 3 ||
        values.length() != 2) {
        return 1;
    }
    row.get_value(y, 0);
    if (!fmpz_is_si(y, -2)) {
        return 1;
    }

    fmpz_set_si(x, 7);
    if (!row.set_entry(3, x) || row.length() != 2) {
        return 1;
    }
    row.get_entry(y, 3);
    if (!fmpz_is_si(y, 7)) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(1, x) || row.length() != 1 || row.col(0) != 3) {
        return 1;
    }
    columns = row.columns();
    values = row.values();
    if (columns.size() != 1 || columns[0] != 3 || values.length() != 1) {
        return 1;
    }
    row.get_entry(y, 1);
    if (!fmpz_is_zero(y)) {
        return 1;
    }

    fmpz_set_ui(x, UWORD(1) << 20);
    fmpz_mul_ui(x, x, UWORD(1) << 20);
    if (!row.append_entry(5, x) || row.length() != 2) {
        return 1;
    }
    row.get_entry(y, 5);
    if (!fmpz_equal(y, x)) {
        return 1;
    }

    fmpz_set_si(x, 11);
    if (row.append_entry(4, x) || row.append_entry(-1, x) ||
        row.set_entry(-1, x) || row.length() != 2) {
        return 1;
    }

    if (!copy.set(row)) {
        return 1;
    }
    fmpz_zero(x);
    if (!row.set_entry(3, x) || row.length() != 1 || copy.length() != 2) {
        return 1;
    }
    copy.get_entry(y, 3);
    if (!fmpz_is_si(y, 7)) {
        return 1;
    }

    fmpz_set_si(x, -9);
    if (!other.set_entry(0, x)) {
        return 1;
    }
    copy.swap(other);
    if (copy.length() != 1 || copy.col(0) != 0 || other.length() != 2) {
        return 1;
    }

    silex::fmpz_smat::SparseRow moved(std::move(other));
    silex::fmpz_smat::SparseRow assigned;
    assigned = std::move(moved);
    if (assigned.length() != 2 || assigned.col(0) != 3 ||
        assigned.col(1) != 5) {
        return 1;
    }
    assigned.get_value(y, 0);
    if (!fmpz_is_si(y, 7)) {
        return 1;
    }

    fmpz_clear(y);
    fmpz_clear(x);
    return 0;
}

int test_sparse_matrix_basic() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::SparseRow got;
    silex::fmpz_smat::SparseRow bad;
    silex::fmpz_smat::SparseMat matrix(4);
    silex::fmpz_smat::SparseMat copy(0);
    fmpz_mat_t dense;
    fmpz_mat_t dense2;
    fmpz_t x;
    fmpz_init(x);
    fmpz_mat_init(dense, 2, 4);
    fmpz_mat_init(dense2, 3, 5);

    if (matrix.nrows() != 0 || matrix.ncols() != 4 || matrix.nnz() != 0) {
        return 1;
    }

    fmpz_set_si(x, 5);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    fmpz_set_si(x, -3);
    if (!row.set_entry(3, x) || !matrix.append_row(row)) {
        return 1;
    }
    if (matrix.nrows() != 1 || matrix.nnz() != 2) {
        return 1;
    }
    const silex::fmpz_smat::SparseRow& first_row = matrix.row_ref(0);
    if (first_row.length() != 2 || first_row.col(0) != 1 ||
        first_row.col(1) != 3) {
        return 1;
    }

    matrix.get_entry(x, 0, 1);
    if (!fmpz_is_si(x, 5)) {
        return 1;
    }
    matrix.get_entry(x, 0, 2);
    if (!fmpz_is_zero(x)) {
        return 1;
    }

    matrix.get_row(got, 0);
    if (got.length() != 2) {
        return 1;
    }
    fmpz_zero(x);
    if (!got.set_entry(1, x)) {
        return 1;
    }
    matrix.get_entry(x, 0, 1);
    if (!fmpz_is_si(x, 5)) {
        return 1;
    }
    silex::flint::FmpzMat owned_row = got.to_fmpz_mat_row(4);
    if (fmpz_mat_nrows(owned_row.raw()) != 1 ||
        fmpz_mat_ncols(owned_row.raw()) != 4 ||
        !entry_is_si(owned_row.raw(), 0, 0, 0) ||
        !entry_is_si(owned_row.raw(), 0, 1, 0) ||
        !entry_is_si(owned_row.raw(), 0, 3, -3)) {
        return 1;
    }

    if (!row.set(got) || !matrix.append_row(row)) {
        return 1;
    }
    if (matrix.nrows() != 2 || matrix.nnz() != 3) {
        return 1;
    }

    set_entry_si(dense, 1, 1, 99);
    matrix.get_fmpz_mat(dense);
    if (!entry_is_si(dense, 0, 0, 0) || !entry_is_si(dense, 0, 1, 5) ||
        !entry_is_si(dense, 0, 3, -3) || !entry_is_si(dense, 1, 1, 0) ||
        !entry_is_si(dense, 1, 3, -3)) {
        return 1;
    }
    silex::flint::FmpzMat owned_dense = matrix.to_fmpz_mat();
    if (!entry_is_si(owned_dense.raw(), 0, 0, 0) ||
        !entry_is_si(owned_dense.raw(), 0, 1, 5) ||
        !entry_is_si(owned_dense.raw(), 0, 3, -3) ||
        !entry_is_si(owned_dense.raw(), 1, 1, 0) ||
        !entry_is_si(owned_dense.raw(), 1, 3, -3)) {
        return 1;
    }

    fmpz_set_si(x, 8);
    if (!bad.set_entry(4, x) || matrix.append_row(bad) ||
        matrix.nrows() != 2 || matrix.nnz() != 3) {
        return 1;
    }

    fmpz_set_si(x, 2);
    if (!row.set_entry(0, x) || !matrix.set_row(1, row) ||
        matrix.nnz() != 4 || matrix.set_row(1, bad) ||
        matrix.set_row(-1, row) || matrix.nnz() != 4) {
        return 1;
    }

    if (!copy.set(matrix)) {
        return 1;
    }
    matrix.swap(copy);
    if (matrix.nrows() != 2 || matrix.ncols() != 4 || matrix.nnz() != 4) {
        return 1;
    }

    fmpz_mat_zero(dense2);
    set_entry_si(dense2, 0, 2, 6);
    set_entry_si(dense2, 1, 0, -1);
    set_entry_si(dense2, 2, 4, 9);
    matrix.set_fmpz_mat(dense2);
    if (matrix.nrows() != 3 || matrix.ncols() != 5 || matrix.nnz() != 3) {
        return 1;
    }
    if (!row.set_fmpz_mat_row(dense2, 2) || row.length() != 1 || row.col(0) != 4) {
        return 1;
    }
    set_entry_si(dense2, 1, 1, 99);
    row.get_fmpz_mat_row(dense2, 1);
    if (!entry_is_si(dense2, 1, 1, 0) || !entry_is_si(dense2, 1, 4, 9)) {
        return 1;
    }
    matrix.get_fmpz_mat(dense2);
    if (!entry_is_si(dense2, 0, 2, 6) || !entry_is_si(dense2, 1, 0, -1) ||
        !entry_is_si(dense2, 2, 4, 9)) {
        return 1;
    }

    if (row.set_fmpz_mat_row(dense2, 3) || row.length() != 1 ||
        matrix.append_fmpz_mat_row(dense, 0) ||
        matrix.append_fmpz_mat_row(dense2, -1) || matrix.nrows() != 3 ||
        !matrix.append_fmpz_mat_row(dense2, 1) || matrix.nrows() != 4 ||
        matrix.nnz() != 4) {
        return 1;
    }

    fmpz_mat_clear(dense2);
    fmpz_mat_clear(dense);
    fmpz_clear(x);
    return 0;
}

bool sparse_rows_match(const silex::fmpz_smat::SparseRow& left,
                       const silex::fmpz_smat::SparseRow& right) noexcept {
    if (left.length() != right.length()) {
        return false;
    }
    for (slong i = 0; i < left.length(); ++i) {
        if (left.columns()[static_cast<std::size_t>(i)] !=
                    right.columns()[static_cast<std::size_t>(i)] ||
            fmpz_equal(left.values().data() + i, right.values().data() + i) ==
                    0) {
            return false;
        }
    }
    return true;
}

// Rows built from dense rows are canonical (increasing columns, no stored
// zeros), so comparing their columns and values decides dense-row equality.
// The class-group duplicate-relation check relies on this.
int test_sparse_row_dense_equality_is_canonical() {
    constexpr slong count = 7;
    constexpr slong ncols = 5;
    fmpz_mat_t dense;
    fmpz_mat_init(dense, count, ncols);
    fmpz_mat_zero(dense);
    // Row 0 stays zero.
    set_entry_si(dense, 1, 0, 1);
    set_entry_si(dense, 2, 0, 1);  // differs from row 1 only in the last column
    set_entry_si(dense, 2, 4, 1);
    set_entry_si(dense, 3, 4, 1);  // differs from row 4 only in sign
    set_entry_si(dense, 4, 4, -1);
    set_entry_si(dense, 5, 1, 2);
    set_entry_si(dense, 5, 2, -3);
    set_entry_si(dense, 6, 1, 2);  // differs from row 5 in a multi-limb value
    fmpz_one(fmpz_mat_entry(dense, 6, 2));
    fmpz_mul_2exp(fmpz_mat_entry(dense, 6, 2),
                  fmpz_mat_entry(dense, 6, 2), 100);
    fmpz_neg(fmpz_mat_entry(dense, 6, 2), fmpz_mat_entry(dense, 6, 2));

    silex::fmpz_smat::SparseMat stored(ncols);
    for (slong i = 0; i < count; ++i) {
        if (!stored.append_fmpz_mat_row(dense, i)) {
            return 1;
        }
    }

    fmpz_mat_t left;
    fmpz_mat_t right;
    fmpz_mat_window_init(left, dense, 0, 0, 1, ncols);
    fmpz_mat_window_init(right, dense, 0, 0, 1, ncols);
    int failed = 0;
    silex::fmpz_smat::SparseRow candidate;
    for (slong i = 0; i < count && failed == 0; ++i) {
        if (!candidate.set_fmpz_mat_row(dense, i)) {
            failed = 1;
            break;
        }
        for (slong j = 0; j < count; ++j) {
            fmpz_mat_window_clear(left);
            fmpz_mat_window_clear(right);
            fmpz_mat_window_init(left, dense, i, 0, i + 1, ncols);
            fmpz_mat_window_init(right, dense, j, 0, j + 1, ncols);
            const bool dense_equal = fmpz_mat_equal(left, right) != 0;
            if (dense_equal != (i == j) ||
                sparse_rows_match(stored.row_ref(j), candidate) !=
                        dense_equal) {
                failed = 1;
                break;
            }
        }
    }
    fmpz_mat_window_clear(right);
    fmpz_mat_window_clear(left);

    // An equal dense row in another matrix gives the same sparse row.
    fmpz_mat_t copy;
    fmpz_mat_init(copy, 1, ncols);
    fmpz_set(fmpz_mat_entry(copy, 0, 1), fmpz_mat_entry(dense, 6, 1));
    fmpz_set(fmpz_mat_entry(copy, 0, 2), fmpz_mat_entry(dense, 6, 2));
    if (failed == 0 && (!candidate.set_fmpz_mat_row(copy, 0) ||
                        !sparse_rows_match(stored.row_ref(6), candidate) ||
                        sparse_rows_match(stored.row_ref(5), candidate))) {
        failed = 1;
    }
    fmpz_mat_clear(copy);
    fmpz_mat_clear(dense);
    return failed;
}

int test_sparse_matrix_kernels() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::SparseMat matrix(4);
    silex::fmpz_smat::SparseMat transpose(0);
    silex::fmpz_smat::SparseMat transpose_twice(0);
    silex::fmpz_smat::SparseMat zero(3);
    fmpz_mat_t dense;
    fmpz_mat_t dense_t;
    fmpz_mat_t dense_tt;
    fmpz_mat_t dense_square;
    fmpz_mat_t xmat;
    fmpz_mat_t ymat;
    fmpz_mat_t square;
    fmpz_mat_t square_expected;
    fmpz_mat_t expected;
    fmpz_mat_t row_dense;
    fmpz_mat_t empty_out;
    fmpz_t x;
    fmpz_t dot;

    fmpz_init(x);
    fmpz_init(dot);
    fmpz_mat_init(dense, 3, 4);
    fmpz_mat_init(dense_t, 4, 3);
    fmpz_mat_init(dense_tt, 3, 4);
    fmpz_mat_init(dense_square, 4, 4);
    fmpz_mat_init(xmat, 4, 2);
    fmpz_mat_init(ymat, 3, 2);
    fmpz_mat_init(square, 4, 2);
    fmpz_mat_init(square_expected, 4, 2);
    fmpz_mat_init(expected, 3, 2);
    fmpz_mat_init(row_dense, 2, 4);
    fmpz_mat_init(empty_out, 0, 2);

    fmpz_set_ui(x, 1);
    fmpz_mul_2exp(x, x, 80);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, -7);
    if (!row.set_entry(3, x) || !matrix.append_row(row)) {
        return 1;
    }
    fmpz_zero(x);
    if (!row.set_entry(0, x) || !row.set_entry(3, x) || !matrix.append_row(row)) {
        return 1;
    }
    fmpz_set_si(x, 5);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    fmpz_set_si(x, -2);
    if (!row.set_entry(2, x) || !matrix.append_row(row)) {
        return 1;
    }

    set_entry_si(xmat, 0, 0, 3);
    set_entry_si(xmat, 0, 1, -1);
    set_entry_si(xmat, 1, 0, 4);
    set_entry_si(xmat, 1, 1, 2);
    set_entry_si(xmat, 2, 0, -6);
    set_entry_si(xmat, 2, 1, 5);
    set_entry_si(xmat, 3, 0, 8);
    set_entry_si(xmat, 3, 1, -3);

    matrix.get_fmpz_mat(dense);
    fmpz_mat_mul(expected, dense, xmat);
    matrix.mul_fmpz_mat(ymat, xmat);
    if (!fmpz_mat_equal(ymat, expected) || !entry_is_si(ymat, 1, 0, 0) ||
        !entry_is_si(ymat, 1, 1, 0)) {
        return 1;
    }

    set_entry_si(row_dense, 1, 0, 11);
    set_entry_si(row_dense, 1, 1, 4);
    set_entry_si(row_dense, 1, 2, -6);
    set_entry_si(row_dense, 1, 3, 3);
    matrix.get_row(row, 2);
    row.dot_fmpz_mat_row(dot, row_dense, 1);
    if (!fmpz_is_si(dot, 32)) {
        return 1;
    }
    row.dot_fmpz_mat_row(fmpz_mat_entry(row_dense, 1, 1), row_dense, 1);
    if (!entry_is_si(row_dense, 1, 1, 32) || !entry_is_si(row_dense, 1, 2, -6)) {
        return 1;
    }

    matrix.transpose(transpose);
    if (transpose.nrows() != 4 || transpose.ncols() != 3 ||
        transpose.nnz() != matrix.nnz()) {
        return 1;
    }
    transpose.get_fmpz_mat(dense_t);
    fmpz_mat_transpose(dense_tt, dense_t);
    if (!fmpz_mat_equal(dense_tt, dense)) {
        return 1;
    }

    transpose.transpose(transpose_twice);
    transpose_twice.get_fmpz_mat(dense_tt);
    if (!fmpz_mat_equal(dense_tt, dense)) {
        return 1;
    }

    matrix.transpose(matrix);
    if (matrix.nrows() != 4 || matrix.ncols() != 3) {
        return 1;
    }
    matrix.get_fmpz_mat(dense_t);
    if (!entry_is_si(dense_t, 3, 0, -7) || !entry_is_si(dense_t, 1, 2, 5) ||
        !entry_is_si(dense_t, 2, 2, -2)) {
        return 1;
    }
    matrix.transpose(matrix);
    matrix.get_fmpz_mat(dense_tt);
    if (!fmpz_mat_equal(dense_tt, dense)) {
        return 1;
    }

    set_entry_si(dense_square, 0, 0, 1);
    set_entry_si(dense_square, 1, 1, 2);
    set_entry_si(dense_square, 2, 2, -1);
    set_entry_si(dense_square, 3, 3, 3);
    matrix.set_fmpz_mat(dense_square);
    set_entry_si(square, 0, 0, 2);
    set_entry_si(square, 0, 1, -1);
    set_entry_si(square, 1, 0, 3);
    set_entry_si(square, 1, 1, 4);
    set_entry_si(square, 2, 0, -5);
    set_entry_si(square, 2, 1, 6);
    set_entry_si(square, 3, 0, 7);
    set_entry_si(square, 3, 1, -8);
    fmpz_mat_mul(square_expected, dense_square, square);
    matrix.mul_fmpz_mat(square, square);
    if (!fmpz_mat_equal(square, square_expected)) {
        return 1;
    }

    fmpz_mat_zero(expected);
    zero.mul_fmpz_mat(empty_out, expected);
    if (fmpz_mat_nrows(empty_out) != 0 || fmpz_mat_ncols(empty_out) != 2 ||
        !fmpz_mat_is_zero(empty_out)) {
        return 1;
    }

    fmpz_clear(dot);
    fmpz_clear(x);
    fmpz_mat_clear(empty_out);
    fmpz_mat_clear(row_dense);
    fmpz_mat_clear(expected);
    fmpz_mat_clear(square_expected);
    fmpz_mat_clear(square);
    fmpz_mat_clear(ymat);
    fmpz_mat_clear(xmat);
    fmpz_mat_clear(dense_square);
    fmpz_mat_clear(dense_tt);
    fmpz_mat_clear(dense_t);
    fmpz_mat_clear(dense);
    return 0;
}

int test_sparse_matrix_nmod() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::SparseMat matrix(4);
    silex::fmpz_smat::SparseMat zero(3);
    nmod_mat_t mod_matrix;
    nmod_mat_t empty;
    fmpz_t x;
    slong rank = -1;

    fmpz_init(x);
    nmod_mat_init(mod_matrix, 3, 4, 5);
    nmod_mat_init(empty, 0, 3, 7);

    fmpz_set_ui(x, 1);
    fmpz_mul_2exp(x, x, 80);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, -7);
    if (!row.set_entry(3, x) || !matrix.append_row(row)) {
        return 1;
    }
    fmpz_zero(x);
    if (!row.set_entry(0, x) || !row.set_entry(3, x) || !matrix.append_row(row)) {
        return 1;
    }
    fmpz_set_si(x, 12);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    fmpz_set_si(x, -13);
    if (!row.set_entry(2, x) || !matrix.append_row(row)) {
        return 1;
    }

    nmod_mat_entry(mod_matrix, 1, 1) = 4;
    matrix.get_nmod_mat(mod_matrix);
    if (nmod_mat_entry(mod_matrix, 0, 0) != 1 ||
        nmod_mat_entry(mod_matrix, 0, 3) != 3 ||
        nmod_mat_entry(mod_matrix, 1, 1) != 0 ||
        nmod_mat_entry(mod_matrix, 2, 1) != 2 ||
        nmod_mat_entry(mod_matrix, 2, 2) != 2) {
        return 1;
    }

    nmod_mat_entry(mod_matrix, 1, 0) = 4;
    matrix.get_row(row, 0);
    row.get_nmod_mat_row(mod_matrix, 1);
    if (nmod_mat_entry(mod_matrix, 1, 0) != 1 ||
        nmod_mat_entry(mod_matrix, 1, 1) != 0 ||
        nmod_mat_entry(mod_matrix, 1, 3) != 3) {
        return 1;
    }

    zero.get_nmod_mat(empty);
    if (nmod_mat_nrows(empty) != 0 || nmod_mat_ncols(empty) != 3) {
        return 1;
    }

    if (!matrix.rank_mod_prime_ui(&rank, 5) || rank != 2) {
        return 1;
    }
    rank = 99;
    if (matrix.rank_mod_prime_ui(&rank, 4) || rank != 99 ||
        matrix.rank_mod_prime_ui(&rank, 1) || rank != 99) {
        return 1;
    }

    nmod_mat_clear(empty);
    nmod_mat_clear(mod_matrix);
    fmpz_clear(x);
    return 0;
}

int test_mod_rank_context() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::SparseMat matrix(2);
    silex::fmpz_smat::SparseMat matrix5(5);
    silex::fmpz_smat::ModRankContext context;
    fmpz_mat_t dense;
    fmpz_mat_t wrong_width;
    fmpz_t x;
    fmpz_t y;
    slong rank = -1;
    slong rank2 = -1;
    bool independent = true;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_mat_init(dense, 3, 2);
    fmpz_mat_init(wrong_width, 1, 3);

    if (context.rank() != 0 || context.ambient_dim() != 0 || context.prime() != 0) {
        return 1;
    }
    if (context.add_row(&independent, row) || !independent) {
        return 1;
    }
    if (context.set_prime_ui(2, 4) || context.ambient_dim() != 0 ||
        context.set_prime_ui(-1, 3) || context.ambient_dim() != 0) {
        return 1;
    }

    if (!context.set_prime_ui(2, 3) || context.ambient_dim() != 2 ||
        context.prime() != 3 || context.rank() != 0) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(0, x) || !row.set_entry(1, x)) {
        return 1;
    }
    independent = true;
    if (!context.add_row(&independent, row) || independent || context.rank() != 0) {
        return 1;
    }

    fmpz_set_si(x, 1);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    independent = false;
    if (!context.add_row(&independent, row) || !independent || context.rank() != 1) {
        return 1;
    }

    fmpz_set_si(x, 2);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    independent = true;
    if (!context.add_row(&independent, row) || independent || context.rank() != 1) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, 3);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    independent = true;
    if (!context.add_row(&independent, row) || independent || context.rank() != 1) {
        return 1;
    }

    fmpz_set_si(x, 1);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    independent = false;
    if (!context.add_row(&independent, row) || !independent || context.rank() != 2) {
        return 1;
    }

    context.reset();
    if (context.rank() != 0 || context.ambient_dim() != 2) {
        return 1;
    }

    set_entry_si(dense, 0, 0, 2);
    set_entry_si(dense, 0, 1, 0);
    set_entry_si(dense, 1, 0, 0);
    set_entry_si(dense, 1, 1, 2);
    set_entry_si(dense, 2, 0, 5);
    set_entry_si(dense, 2, 1, -7);

    independent = false;
    if (!context.add_fmpz_mat_row(&independent, dense, 0) || !independent ||
        !context.add_fmpz_mat_row(&independent, dense, 1) || !independent ||
        !context.add_fmpz_mat_row(&independent, dense, 2) || independent ||
        context.rank() != 2) {
        return 1;
    }

    independent = true;
    if (context.add_fmpz_mat_row(&independent, dense, 3) || !independent ||
        context.add_fmpz_mat_row(&independent, wrong_width, 0) || !independent ||
        context.rank() != 2) {
        return 1;
    }

    context.reset();
    matrix.set_fmpz_mat(dense);
    if (!matrix.rank_mod_prime_ui(&rank, 3)) {
        return 1;
    }
    matrix.get_row(row, 0);
    if (!context.add_row(&independent, row)) {
        return 1;
    }
    matrix.get_row(row, 1);
    if (!context.add_row(&independent, row)) {
        return 1;
    }
    rank2 = context.rank();
    if (rank2 != rank) {
        return 1;
    }

    if (!context.set_prime_ui(2, 2)) {
        return 1;
    }
    matrix.get_row(row, 0);
    if (!context.add_row(&independent, row) || independent) {
        return 1;
    }
    matrix.get_row(row, 1);
    if (!context.add_row(&independent, row) || independent || context.rank() != 0) {
        return 1;
    }

    fmpz_set_ui(y, 1);
    fmpz_mul_2exp(y, y, 90);
    fmpz_add_ui(y, y, 1);
    fmpz_zero(x);
    matrix = silex::fmpz_smat::SparseMat(2);
    if (!row.set_entry(1, x) || !row.set_entry(0, y) || !matrix.append_row(row)) {
        return 1;
    }
    fmpz_set_si(x, -10);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, 11);
    if (!row.set_entry(1, x) || !matrix.append_row(row)) {
        return 1;
    }
    if (!matrix.rank_mod_prime_ui(&rank, 5) || !context.set_prime_ui(2, 5)) {
        return 1;
    }
    matrix.get_row(row, 0);
    if (!context.add_row(&independent, row)) {
        return 1;
    }
    matrix.get_row(row, 1);
    if (!context.add_row(&independent, row) || context.rank() != rank) {
        return 1;
    }

    if (!context.set_prime_ui(5, 7)) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(0, x) || !row.set_entry(1, x) || !row.set_entry(2, x) ||
        !row.set_entry(3, x) || !row.set_entry(4, x)) {
        return 1;
    }
    fmpz_set_si(x, 1);
    if (!row.set_entry(3, x) || !matrix5.append_row(row)) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(3, x)) {
        return 1;
    }
    fmpz_set_si(x, -2);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    fmpz_set_si(x, 5);
    if (!row.set_entry(4, x) || !matrix5.append_row(row)) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(1, x) || !row.set_entry(4, x)) {
        return 1;
    }
    fmpz_set_si(x, 9);
    if (!row.set_entry(0, x) || !matrix5.append_row(row)) {
        return 1;
    }

    fmpz_set_si(x, 6);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, 3);
    if (!row.set_entry(2, x) || !matrix5.append_row(row)) {
        return 1;
    }

    fmpz_set_si(x, 1);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    fmpz_set_si(x, 4);
    if (!row.set_entry(3, x) || !matrix5.append_row(row)) {
        return 1;
    }

    for (slong i = 0; i < matrix5.nrows(); ++i) {
        matrix5.get_row(row, i);
        if (!context.add_row(&independent, row)) {
            return 1;
        }

        silex::fmpz_smat::SparseMat prefix(5);
        for (slong j = 0; j <= i; ++j) {
            matrix5.get_row(row, j);
            if (!prefix.append_row(row)) {
                return 1;
            }
        }
        if (!prefix.rank_mod_prime_ui(&rank, 7) || context.rank() != rank) {
            return 1;
        }
    }
    if (context.rank() != 5) {
        return 1;
    }

    fmpz_clear(y);
    fmpz_clear(x);
    fmpz_mat_clear(wrong_width);
    fmpz_mat_clear(dense);
    return 0;
}

int test_mod_rank_context_row_is_independent() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::SparseRow other;
    silex::fmpz_smat::ModRankContext context;
    silex::flint::Fmpz x;
    bool independent = false;

    if (context.row_is_independent(&independent, row)) {
        return 1;
    }

    if (!context.set_prime_ui(3, 5)) {
        return 1;
    }

    independent = true;
    if (!context.row_is_independent(&independent, row) || independent ||
        context.rank() != 0) {
        return 1;
    }

    silex::flint::fmpz_set_si(x, 1);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    independent = false;
    if (!context.row_is_independent(&independent, row) || !independent ||
        context.rank() != 0) {
        return 1;
    }

    if (!context.add_row(&independent, row) || !independent ||
        context.rank() != 1) {
        return 1;
    }

    silex::flint::fmpz_set_si(x, 6);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    independent = true;
    if (!context.row_is_independent(&independent, row) || independent ||
        context.rank() != 1) {
        return 1;
    }

    silex::flint::fmpz_set_si(x, 3);
    if (!other.set_entry(0, x)) {
        return 1;
    }
    silex::flint::fmpz_set_si(x, 1);
    if (!other.set_entry(1, x)) {
        return 1;
    }
    independent = false;
    if (!context.row_is_independent(&independent, other) || !independent ||
        context.rank() != 1) {
        return 1;
    }
    if (!context.add_row(&independent, other) || !independent ||
        context.rank() != 2) {
        return 1;
    }

    silex::flint::fmpz_set_si(x, 1);
    if (!other.set_entry(3, x)) {
        return 1;
    }
    independent = true;
    if (context.row_is_independent(&independent, other) || !independent ||
        context.rank() != 2) {
        return 1;
    }

    return 0;
}

int test_hnf_context() {
    silex::fmpz_smat::SparseRow row;
    silex::fmpz_smat::HnfContext context;
    fmpz_mat_t hnf;
    fmpz_mat_t dense;
    fmpz_mat_t wrong_width;
    fmpz_t x;
    fmpz_t index;
    bool independent = false;
    bool index_refined = false;

    fmpz_mat_init(hnf, 2, 2);
    fmpz_mat_init(dense, 3, 2);
    fmpz_mat_init(wrong_width, 1, 3);
    fmpz_init(x);
    fmpz_init(index);

    if (context.rank() != 0 || context.ambient_dim() != 0 ||
        !context.reset(2, 0) || context.rank() != 0 ||
        context.ambient_dim() != 2 || context.has_mod_rank() ||
        context.reset(-1, 0) || context.ambient_dim() != 2 ||
        context.reset(3, 4) || context.ambient_dim() != 2) {
        return 1;
    }

    fmpz_set_si(x, 2);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    independent = false;
    if (!context.add_row(&independent, row) || !independent ||
        context.rank() != 1 || context.full_rank_index(index)) {
        return 1;
    }

    fmpz_zero(x);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    independent = true;
    if (!context.add_row(&independent, row) || independent ||
        context.rank() != 1) {
        return 1;
    }

    fmpz_set_si(x, 3);
    if (!row.set_entry(1, x)) {
        return 1;
    }
    independent = false;
    if (!context.add_row(&independent, row) || !independent ||
        context.rank() != 2 || !context.full_rank_index(index) ||
        !fmpz_is_si(index, 6) || !context.get_hnf(hnf) ||
        !entry_is_si(hnf, 0, 0, 2) || !entry_is_si(hnf, 0, 1, 0) ||
        !entry_is_si(hnf, 1, 0, 0) || !entry_is_si(hnf, 1, 1, 3)) {
        return 1;
    }

    fmpz_set_si(x, 1);
    if (!row.set_entry(2, x)) {
        return 1;
    }
    independent = true;
    if (context.add_row(&independent, row) || !independent ||
        context.rank() != 2) {
        return 1;
    }

    if (!context.reset(1, 0)) {
        return 1;
    }
    fmpz_mat_clear(hnf);
    fmpz_mat_init(hnf, 1, 1);
    fmpz_zero(x);
    if (!row.set_entry(1, x) || !row.set_entry(2, x)) {
        return 1;
    }
    fmpz_set_si(x, 2);
    index_refined = true;
    if (!row.set_entry(0, x) ||
        !context.add_row(&independent, &index_refined, row) ||
        !independent || index_refined || !context.full_rank_index(index) ||
        !fmpz_is_si(index, 2)) {
        return 1;
    }
    fmpz_set_si(x, 1);
    index_refined = false;
    if (!row.set_entry(0, x) ||
        !context.add_row(&independent, &index_refined, row) ||
        independent || !index_refined || context.rank() != 1 ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 1) ||
        !context.get_hnf(hnf) || !entry_is_si(hnf, 0, 0, 1)) {
        return 1;
    }
    index_refined = true;
    if (!context.add_row(&independent, &index_refined, row) ||
        independent || index_refined || context.rank() != 1 ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 1)) {
        return 1;
    }

    if (!context.reset(2, 3) || !context.has_mod_rank()) {
        return 1;
    }
    fmpz_mat_clear(hnf);
    fmpz_mat_init(hnf, 2, 2);
    fmpz_set_si(x, 1);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_zero(x);
    if (!row.set_entry(1, x) || !context.add_row(&independent, row) ||
        !independent) {
        return 1;
    }
    fmpz_zero(x);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, 3);
    if (!row.set_entry(1, x) || !context.add_row(&independent, row) ||
        !independent || context.rank() != 2 ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 3)) {
        return 1;
    }
    fmpz_set_si(x, 6);
    if (!row.set_entry(1, x) || !context.add_row(&independent, row) ||
        independent || context.rank() != 2 || !context.get_hnf(hnf) ||
        !entry_is_si(hnf, 0, 0, 1) || !entry_is_si(hnf, 1, 1, 3)) {
        return 1;
    }
    silex::fmpz_smat::HnfContext copied;
    if (!copied.set(context) || !copied.has_mod_rank() ||
        !copied.full_rank_index(index) || !fmpz_is_si(index, 3)) {
        return 1;
    }
    fmpz_zero(x);
    if (!row.set_entry(0, x)) {
        return 1;
    }
    fmpz_set_si(x, 1);
    if (!row.set_entry(1, x) || !copied.add_row(&independent, row) ||
        independent || copied.rank() != 2 ||
        !copied.full_rank_index(index) || !fmpz_is_si(index, 1) ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 3)) {
        return 1;
    }

    if (!context.reset(2, 0)) {
        return 1;
    }
    set_entry_si(dense, 0, 0, 2);
    set_entry_si(dense, 0, 1, 0);
    set_entry_si(dense, 1, 0, 0);
    set_entry_si(dense, 1, 1, 5);
    set_entry_si(dense, 2, 0, 4);
    set_entry_si(dense, 2, 1, 10);
    if (!context.add_fmpz_mat_row(&independent, dense, 0) || !independent ||
        !context.add_fmpz_mat_row(&independent, dense, 1) || !independent ||
        !context.add_fmpz_mat_row(&independent, dense, 2) || independent ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 10)) {
        return 1;
    }
    bool reduces_to_zero = false;
    if (!context.fmpz_mat_row_reduces_to_zero(
                reduces_to_zero, silex::flint::FmpzMatConstRef(dense), 2) ||
        !reduces_to_zero || context.rank() != 2 ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 10)) {
        return 1;
    }
    set_entry_si(dense, 2, 0, 1);
    set_entry_si(dense, 2, 1, 0);
    if (!context.fmpz_mat_row_reduces_to_zero(
                reduces_to_zero, silex::flint::FmpzMatConstRef(dense), 2) ||
        reduces_to_zero || context.rank() != 2 ||
        !context.full_rank_index(index) || !fmpz_is_si(index, 10)) {
        return 1;
    }
    independent = true;
    if (context.add_fmpz_mat_row(&independent, dense, 3) || !independent ||
        context.add_fmpz_mat_row(&independent, wrong_width, 0) ||
        !independent ||
        context.fmpz_mat_row_reduces_to_zero(
                reduces_to_zero,
                silex::flint::FmpzMatConstRef(wrong_width), 0)) {
        return 1;
    }

    fmpz_clear(index);
    fmpz_clear(x);
    fmpz_mat_clear(wrong_width);
    fmpz_mat_clear(dense);
    fmpz_mat_clear(hnf);
    return 0;
}

int test_hnf_context_transform() {
    silex::fmpz_smat::HnfContext context;
    silex::flint::FmpzMat basis(2, 2);
    silex::flint::FmpzMat hnf(2, 2);
    silex::flint::FmpzMat transform(2, 2);
    silex::flint::FmpzMat product(2, 2);
    bool independent = false;

    silex::flint::fmpz_set_si(
            silex::flint::fmpz_mat_entry(basis, 0, 0), 2);
    silex::flint::fmpz_set_si(
            silex::flint::fmpz_mat_entry(basis, 0, 1), 1);
    silex::flint::fmpz_set_si(
            silex::flint::fmpz_mat_entry(basis, 1, 1), 3);

    if (!context.reset(2, 0) ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(basis), 0) ||
        !independent ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(basis), 1) ||
        !independent ||
        !context.get_hnf_rows(silex::flint::FmpzMatRef(hnf)) ||
        !context.get_hnf_transform(silex::flint::FmpzMatRef(transform))) {
        return 1;
    }

    silex::flint::fmpz_mat_mul(silex::flint::FmpzMatRef(product),
                               silex::flint::FmpzMatConstRef(transform),
                               silex::flint::FmpzMatConstRef(basis));
    if (!silex::flint::fmpz_mat_equal(
                silex::flint::FmpzMatConstRef(product),
                silex::flint::FmpzMatConstRef(hnf))) {
        return 1;
    }

    silex::flint::FmpzMat dependent(1, 2);
    silex::flint::fmpz_set_si(
            silex::flint::fmpz_mat_entry(dependent, 0, 0), 1);
    if (!context.add_fmpz_mat_row(
                &independent, silex::flint::FmpzMatConstRef(dependent), 0) ||
        independent ||
        !context.get_hnf_rows(silex::flint::FmpzMatRef(hnf)) ||
        !context.get_hnf_transform(silex::flint::FmpzMatRef(transform))) {
        return 1;
    }
    silex::flint::fmpz_mat_mul(silex::flint::FmpzMatRef(product),
                               silex::flint::FmpzMatConstRef(transform),
                               silex::flint::FmpzMatConstRef(hnf));
    if (!silex::flint::fmpz_mat_equal(
                silex::flint::FmpzMatConstRef(product),
                silex::flint::FmpzMatConstRef(hnf))) {
        return 1;
    }

    silex::fmpz_smat::HnfContext copied;
    return (!copied.set(context) ||
            !copied.get_hnf_rows(silex::flint::FmpzMatRef(hnf)) ||
            !copied.get_hnf_transform(silex::flint::FmpzMatRef(transform)))
            ? 1
            : 0;
}

int test_hnf_context_index_before_transform() {
    silex::fmpz_smat::HnfContext context;
    silex::flint::FmpzMat basis(2, 2);
    silex::flint::FmpzMat hnf(2, 2);
    silex::flint::FmpzMat transform(2, 2);
    silex::flint::FmpzMat product(2, 2);
    silex::flint::Fmpz index;
    bool independent = false;

    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 0, 1, 1);
    set_entry_si(basis, 1, 1, 3);

    if (!context.reset(2, 0) ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(basis), 0) ||
        !independent ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(basis), 1) ||
        !independent ||
        !context.full_rank_index(silex::flint::FmpzRef(index)) ||
        !fmpz_is_si(index.raw(), 6) ||
        !context.get_hnf(silex::flint::FmpzMatRef(hnf)) ||
        !context.get_hnf_transform(silex::flint::FmpzMatRef(transform))) {
        return 1;
    }

    silex::flint::fmpz_mat_mul(silex::flint::FmpzMatRef(product),
                               silex::flint::FmpzMatConstRef(transform),
                               silex::flint::FmpzMatConstRef(basis));
    return silex::flint::fmpz_mat_equal(
                   silex::flint::FmpzMatConstRef(product),
                   silex::flint::FmpzMatConstRef(hnf)) != 0
            ? 0
            : 1;
}

int test_hnf_context_incremental_refresh() {
    silex::fmpz_smat::HnfContext context;
    silex::flint::FmpzMat rows(5, 3);
    silex::flint::Fmpz index;
    bool independent = false;

    set_entry_si(rows, 0, 0, 2);
    set_entry_si(rows, 0, 1, 1);
    set_entry_si(rows, 1, 1, 3);
    set_entry_si(rows, 1, 2, 1);
    set_entry_si(rows, 2, 0, 4);
    set_entry_si(rows, 2, 1, 2);
    set_entry_si(rows, 3, 0, 1);
    set_entry_si(rows, 3, 1, -1);
    set_entry_si(rows, 3, 2, 5);
    set_entry_si(rows, 4, 0, 3);
    set_entry_si(rows, 4, 1, 2);
    set_entry_si(rows, 4, 2, -1);

    if (!context.reset(3, 0) ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 0) ||
        !independent || !hnf_context_matches_prefix(context, rows, 1) ||
        context.full_rank_index(silex::flint::FmpzRef(index))) {
        return 1;
    }

    if (!context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 1) ||
        !independent || !hnf_context_matches_prefix(context, rows, 2) ||
        context.full_rank_index(silex::flint::FmpzRef(index))) {
        return 1;
    }

    if (!context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 2) ||
        independent || !hnf_context_matches_prefix(context, rows, 3)) {
        return 1;
    }

    if (!context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 3) ||
        !independent || !hnf_context_matches_prefix(context, rows, 4) ||
        !context.full_rank_index(silex::flint::FmpzRef(index))) {
        return 1;
    }

    if (!context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 4) ||
        independent || !hnf_context_matches_prefix(context, rows, 5) ||
        !context.full_rank_index(silex::flint::FmpzRef(index))) {
        return 1;
    }

    return 0;
}

int test_hnf_context_deferred_dependent_refresh() {
    silex::fmpz_smat::HnfContext context;
    silex::flint::FmpzMat rows(4, 2);
    silex::flint::Fmpz index;
    bool independent = false;

    set_entry_si(rows, 0, 0, 2);
    set_entry_si(rows, 1, 1, 5);
    set_entry_si(rows, 2, 0, 1);
    set_entry_si(rows, 3, 1, 1);

    if (!context.reset(2, 0) ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 0) ||
        !independent ||
        !context.add_fmpz_mat_row(&independent,
                                  silex::flint::FmpzMatConstRef(rows), 1) ||
        !independent ||
        !context.full_rank_index(silex::flint::FmpzRef(index)) ||
        !fmpz_is_si(index.raw(), 10)) {
        return 1;
    }

    if (!context.add_fmpz_mat_row_defer_dependent(
                &independent, silex::flint::FmpzMatConstRef(rows), 2) ||
        independent) {
        return 1;
    }

    silex::fmpz_smat::HnfContext copied;
    if (!copied.set(context) ||
        !copied.full_rank_index(silex::flint::FmpzRef(index)) ||
        !fmpz_is_si(index.raw(), 5)) {
        return 1;
    }

    if (!context.full_rank_index(silex::flint::FmpzRef(index)) ||
        !fmpz_is_si(index.raw(), 5) ||
        !context.add_fmpz_mat_row_defer_dependent(
                &independent, silex::flint::FmpzMatConstRef(rows), 3) ||
        independent ||
        !context.full_rank_index(silex::flint::FmpzRef(index)) ||
        !fmpz_is_si(index.raw(), 1)) {
        return 1;
    }

    silex::flint::FmpzMat hnf(2, 2);
    silex::flint::FmpzMat transform(2, 2);
    silex::flint::FmpzMat product(2, 2);
    if (!context.get_hnf(silex::flint::FmpzMatRef(hnf)) ||
        !entry_is_si(hnf.raw(), 0, 0, 1) ||
        !entry_is_si(hnf.raw(), 0, 1, 0) ||
        !entry_is_si(hnf.raw(), 1, 0, 0) ||
        !entry_is_si(hnf.raw(), 1, 1, 1)) {
        return 1;
    }

    if (!context.get_hnf_transform(silex::flint::FmpzMatRef(transform))) {
        return 1;
    }
    silex::flint::fmpz_mat_mul(silex::flint::FmpzMatRef(product),
                               silex::flint::FmpzMatConstRef(transform),
                               silex::flint::FmpzMatConstRef(hnf));
    if (!silex::flint::fmpz_mat_equal(
                silex::flint::FmpzMatConstRef(product),
                silex::flint::FmpzMatConstRef(hnf))) {
        return 1;
    }

    return 0;
}

// Deterministic generator for the randomized HNF row-insertion checks.
ulong next_random(ulong& state) noexcept {
    state = state * 6364136223846793005UL + 1442695040888963407UL;
    return state >> 33;
}

// Sets `row` of `matrix` to zero, then to a sparse row shaped like a
// class-group relation: a few nonzero entries of at most max_abs in absolute
// value at random columns from first_column on.
void fill_relation_like_row(silex::flint::FmpzMat& matrix,
                            slong row,
                            slong first_column,
                            slong max_abs,
                            ulong& state) noexcept {
    const slong cols = silex::flint::fmpz_mat_ncols(matrix);
    for (slong j = 0; j < cols; ++j) {
        ::fmpz_zero(::fmpz_mat_entry(matrix.raw(), row, j));
    }
    const slong span = cols - first_column;
    if (span <= 0) {
        return;
    }
    const slong terms = 1 + static_cast<slong>(next_random(state) % 6);
    for (slong t = 0; t < terms; ++t) {
        const slong col = first_column +
                          static_cast<slong>(next_random(state) %
                                             static_cast<ulong>(span));
        const slong value =
                static_cast<slong>(next_random(state) %
                                   static_cast<ulong>(2 * max_abs + 1)) -
                max_abs;
        ::fmpz_set_si(::fmpz_mat_entry(matrix.raw(), row, col), value);
    }
}

// Sets `row` of `matrix` to left_coeff * row `left` + right_coeff * row
// `right`: a row dependent on earlier rows, in their lattice.
void set_combination_row(silex::flint::FmpzMat& matrix,
                         slong row,
                         slong left,
                         slong left_coeff,
                         slong right,
                         slong right_coeff) noexcept {
    const slong cols = silex::flint::fmpz_mat_ncols(matrix);
    for (slong j = 0; j < cols; ++j) {
        fmpz* entry = ::fmpz_mat_entry(matrix.raw(), row, j);
        ::fmpz_mul_si(entry, ::fmpz_mat_entry(matrix.raw(), left, j),
                      left_coeff);
        ::fmpz_addmul_si(entry, ::fmpz_mat_entry(matrix.raw(), right, j),
                         right_coeff);
    }
}

// Stacks the trimmed fmpz_mat_hnf of rows [0, old_rows) of `rows` over rows
// [old_rows, nrows), folds the new rows in with fmpz_mat_hnf_add_rows, and
// compares the result with fmpz_mat_hnf of all rows: the same nonzero rows,
// then zero rows.
bool hnf_add_rows_matches_full_hnf(silex::flint::FmpzMatConstRef rows,
                                   slong old_rows) noexcept {
    const slong nrows = silex::flint::fmpz_mat_nrows(rows);
    const slong ncols = silex::flint::fmpz_mat_ncols(rows);
    silex::flint::FmpzMatConstWindow old_window(rows, 0, 0, old_rows, ncols);
    silex::flint::FmpzMat old_hnf(old_rows, ncols);
    ::fmpz_mat_hnf(old_hnf.raw(), old_window.raw());
    const slong old_rank =
            trimmed_nonzero_rows(silex::flint::FmpzMatConstRef(old_hnf));

    const slong stacked_rows = old_rank + nrows - old_rows;
    silex::flint::FmpzMat stacked(stacked_rows, ncols);
    for (slong i = 0; i < old_rank; ++i) {
        for (slong j = 0; j < ncols; ++j) {
            ::fmpz_set(::fmpz_mat_entry(stacked.raw(), i, j),
                       ::fmpz_mat_entry(old_hnf.raw(), i, j));
        }
    }
    for (slong i = old_rows; i < nrows; ++i) {
        for (slong j = 0; j < ncols; ++j) {
            ::fmpz_set(::fmpz_mat_entry(stacked.raw(), old_rank + i - old_rows,
                                        j),
                       ::fmpz_mat_entry(rows.raw(), i, j));
        }
    }

    silex::flint::FmpzMat expected(nrows, ncols);
    ::fmpz_mat_hnf(expected.raw(), rows.raw());
    const slong expected_rank =
            trimmed_nonzero_rows(silex::flint::FmpzMatConstRef(expected));

    slong rank = -1;
    if (!silex::detail::fmpz_mat_hnf_add_rows(
                rank, silex::flint::FmpzMatRef(stacked), old_rank) ||
        rank != expected_rank ||
        rank != ::fmpz_mat_rank(rows.raw())) {
        return false;
    }
    for (slong i = 0; i < stacked_rows; ++i) {
        for (slong j = 0; j < ncols; ++j) {
            const fmpz* got = ::fmpz_mat_entry(stacked.raw(), i, j);
            const bool equal =
                    i < rank ? ::fmpz_equal(got, ::fmpz_mat_entry(
                                                         expected.raw(), i,
                                                         j)) != 0
                             : ::fmpz_is_zero(got) != 0;
            if (!equal) {
                return false;
            }
        }
    }
    return true;
}

int test_hnf_add_rows_random_small() {
    ulong state = 20261012UL;
    for (slong trial = 0; trial < 400; ++trial) {
        const slong cols = 1 + static_cast<slong>(next_random(state) % 8);
        const slong old_rows = static_cast<slong>(next_random(state) % 9);
        const slong new_rows = static_cast<slong>(next_random(state) % 5);
        const slong rows = old_rows + new_rows;
        const slong max_abs = 1 + static_cast<slong>(next_random(state) % 40);
        silex::flint::FmpzMat matrix(rows, cols);
        for (slong i = 0; i < rows; ++i) {
            // Shifted supports make later rows insert pivots in front of
            // earlier ones; some rows are dependent or zero.
            const slong kind = static_cast<slong>(next_random(state) % 6);
            if (kind == 0 && i >= 2) {
                set_combination_row(matrix, i, i - 1, 3, i - 2, -2);
            } else if (kind == 1 && i >= 1) {
                set_combination_row(matrix, i, i - 1, 5, i - 1, 0);
            } else if (kind == 2) {
                fill_relation_like_row(matrix, i, cols, max_abs, state);
            } else {
                const slong first = static_cast<slong>(
                        next_random(state) % static_cast<ulong>(cols));
                fill_relation_like_row(matrix, i, first, max_abs, state);
            }
        }
        if (!hnf_add_rows_matches_full_hnf(
                    silex::flint::FmpzMatConstRef(matrix), old_rows)) {
            return 1;
        }
    }

    // Malformed retained rows are rejected and left unchanged.
    silex::flint::FmpzMat bad(2, 2);
    set_entry_si(bad, 0, 0, -1);
    set_entry_si(bad, 1, 1, 1);
    silex::flint::FmpzMat bad_copy(2, 2);
    silex::flint::fmpz_mat_set(silex::flint::FmpzMatRef(bad_copy),
                               silex::flint::FmpzMatConstRef(bad));
    slong rank = -1;
    if (silex::detail::fmpz_mat_hnf_add_rows(
                rank, silex::flint::FmpzMatRef(bad), 1) ||
        silex::detail::fmpz_mat_hnf_add_rows(
                rank, silex::flint::FmpzMatRef(bad), 3) ||
        !silex::flint::fmpz_mat_equal(
                silex::flint::FmpzMatConstRef(bad),
                silex::flint::FmpzMatConstRef(bad_copy))) {
        return 1;
    }

    // Pivots out of order (row 0 pivots in column 1, row 1 in column 0) and
    // a zero row inside the claimed prefix are rejected, matrix unchanged.
    silex::flint::FmpzMat unordered(3, 2);
    set_entry_si(unordered, 0, 1, 1);
    set_entry_si(unordered, 1, 0, 1);
    set_entry_si(unordered, 2, 0, 1);
    silex::flint::FmpzMat gap(3, 3);
    set_entry_si(gap, 0, 0, 1);
    set_entry_si(gap, 2, 2, 1);
    for (silex::flint::FmpzMat* matrix : {&unordered, &gap}) {
        const slong prefix = 2;
        silex::flint::FmpzMat copy(silex::flint::fmpz_mat_nrows(
                                           silex::flint::FmpzMatConstRef(
                                                   *matrix)),
                                   silex::flint::fmpz_mat_ncols(
                                           silex::flint::FmpzMatConstRef(
                                                   *matrix)));
        silex::flint::fmpz_mat_set(silex::flint::FmpzMatRef(copy),
                                   silex::flint::FmpzMatConstRef(*matrix));
        if (silex::detail::fmpz_mat_hnf_add_rows(
                    rank, silex::flint::FmpzMatRef(*matrix), prefix) ||
            !silex::flint::fmpz_mat_equal(
                    silex::flint::FmpzMatConstRef(*matrix),
                    silex::flint::FmpzMatConstRef(copy))) {
            return 1;
        }
    }
    return 0;
}

// Rows dumped from the first refreshes of HnfContext with a retained HNF
// in the cyclotomic_64_phi32 class-unit instance (161 columns): the retained
// HNF rows followed by the newly accepted row.  Entries are (row, col, value)
// triples of the nonzero entries.
struct CorpusHnfEntry { slong row; slong col; slong value; };
constexpr CorpusHnfEntry corpus_hnf_refresh_0[] = {
    {0, 0, 1}, {0, 110, -1}, {0, 149, 1}, {1, 2, 1},
    {1, 156, -1}, {1, 157, 1}, {2, 3, 1}, {2, 82, -1},
    {2, 156, -1}, {2, 157, 1}, {3, 5, 1}, {3, 82, 1},
    {3, 141, 1}, {3, 148, -1}, {4, 6, 1}, {4, 25, 1},
    {5, 7, 1}, {5, 141, 1}, {5, 148, -1}, {5, 156, -1},
    {5, 157, 1}, {6, 9, 1}, {6, 109, -1}, {6, 137, 1},
    {7, 15, 1}, {7, 83, 1}, {8, 17, 1}, {8, 84, 1},
    {9, 19, 1}, {9, 103, 1}, {9, 136, 1}, {10, 20, 1},
    {10, 106, 1}, {10, 109, -1}, {11, 24, 1}, {11, 142, -1},
    {11, 158, -1}, {11, 159, 1}, {12, 27, 1}, {12, 136, -1},
    {13, 28, 1}, {13, 156, 1}, {13, 157, -1}, {14, 31, 1},
    {14, 157, 1}, {14, 158, -1}, {15, 32, 1}, {15, 109, 1},
    {15, 141, -1}, {15, 158, -1}, {15, 159, 1}, {16, 37, 1},
    {16, 55, 1}, {16, 156, 1}, {16, 157, -1}, {17, 38, 1},
    {17, 82, -1}, {18, 39, 1}, {18, 156, -1}, {18, 157, 1},
    {19, 40, 1}, {19, 56, 1}, {20, 44, 1}, {20, 156, 1},
    {20, 157, -1}, {20, 158, 1}, {20, 159, -1}, {21, 45, 1},
    {21, 155, 1}, {22, 46, 1}, {22, 147, 1}, {23, 48, 1},
    {23, 52, 1}, {23, 156, 1}, {23, 157, -1}, {24, 49, 1},
    {24, 142, 1}, {24, 156, 1}, {24, 157, -1}, {25, 50, 1},
    {25, 145, 1}, {26, 57, 1}, {26, 142, 1}, {27, 59, 3},
    {27, 127, 1}, {27, 134, -3}, {27, 157, 3}, {27, 158, -2},
    {27, 159, 5}, {28, 60, 1}, {28, 98, 1}, {29, 63, 1},
    {29, 141, 1}, {30, 65, 1}, {30, 146, -1}, {30, 148, 1},
    {30, 156, 1}, {30, 157, -1}, {31, 66, 1}, {31, 72, 1},
    {31, 130, 1}, {32, 68, 1}, {32, 155, 1}, {33, 69, 1},
    {33, 106, 1}, {33, 110, -1}, {33, 149, 1}, {34, 73, 1},
    {34, 156, 1}, {34, 157, -1}, {35, 74, 1}, {35, 146, 1},
    {35, 148, -1}, {36, 76, 1}, {36, 156, 1}, {36, 157, -1},
    {37, 78, 1}, {37, 119, 1}, {37, 156, 1}, {37, 157, -1},
    {38, 79, 1}, {38, 81, 1}, {38, 156, 1}, {38, 157, -1},
    {39, 80, 1}, {39, 136, 1}, {39, 156, 1}, {39, 157, -1},
    {40, 86, 1}, {40, 99, 1}, {40, 141, -1}, {40, 148, 1},
    {40, 156, 1}, {40, 157, -1}, {41, 87, 1}, {41, 148, -1},
    {42, 88, 1}, {42, 141, -1}, {42, 156, 1}, {42, 157, -1},
    {43, 89, 1}, {43, 103, 1}, {43, 156, -1}, {43, 157, 1},
    {44, 90, 1}, {44, 134, 1}, {45, 91, 1}, {45, 110, 1},
    {46, 92, 1}, {46, 109, 1}, {46, 156, 1}, {46, 157, -1},
    {47, 94, 1}, {47, 103, -1}, {47, 156, 1}, {47, 157, -1},
    {48, 95, 1}, {48, 137, 1}, {48, 156, 1}, {48, 157, -1},
    {49, 96, 1}, {49, 107, 1}, {49, 156, 1}, {49, 157, -1},
    {50, 100, 1}, {50, 156, 1}, {50, 157, -1}, {51, 101, 1},
    {51, 151, 1}, {51, 155, -1}, {52, 102, 1}, {52, 109, 1},
    {52, 148, 1}, {52, 156, 1}, {52, 157, -1}, {53, 105, 1},
    {53, 145, 1}, {53, 156, 1}, {53, 157, -1}, {54, 111, 1},
    {54, 153, 1}, {55, 112, 1}, {55, 156, -1}, {55, 157, 1},
    {56, 114, 1}, {56, 133, 1}, {56, 156, 1}, {56, 157, -1},
    {57, 115, 1}, {57, 127, 1}, {57, 156, 1}, {57, 157, -1},
    {58, 116, 1}, {58, 136, 1}, {59, 118, 1}, {59, 141, 1},
    {59, 156, -1}, {59, 157, 1}, {60, 120, 1}, {60, 156, 1},
    {60, 157, -1}, {61, 121, 1}, {61, 134, 1}, {61, 156, -1},
    {62, 122, 1}, {62, 148, 1}, {62, 156, 1}, {62, 157, -1},
    {63, 123, 1}, {63, 155, 1}, {63, 156, 1}, {63, 157, -1},
    {64, 124, 1}, {64, 156, -1}, {64, 157, 1}, {65, 125, 1},
    {65, 148, 1}, {66, 126, 1}, {66, 146, 1}, {67, 128, 1},
    {67, 142, -1}, {67, 156, 1}, {68, 129, 1}, {68, 156, 1},
    {68, 157, -1}, {69, 131, 1}, {69, 141, -1}, {69, 148, 1},
    {69, 156, 1}, {69, 157, -1}, {70, 135, 1}, {70, 145, 1},
    {70, 156, 1}, {70, 157, -1}, {71, 152, 1}, {71, 156, 1},
    {71, 157, -1}, {72, 160, 1}, {73, 9, 1}, {73, 19, -1},
    {73, 156, 1}, {73, 157, -1}, {73, 160, 2},
};
constexpr CorpusHnfEntry corpus_hnf_refresh_1[] = {
    {0, 0, 1}, {0, 110, -1}, {0, 149, 1}, {1, 2, 1},
    {1, 156, -1}, {1, 157, 1}, {2, 3, 1}, {2, 82, -1},
    {2, 156, -1}, {2, 157, 1}, {3, 5, 1}, {3, 82, 1},
    {3, 148, -1}, {3, 157, -1}, {4, 6, 1}, {4, 25, 1},
    {5, 7, 1}, {5, 148, -1}, {5, 156, -1}, {6, 9, 1},
    {6, 137, 1}, {6, 156, 1}, {6, 157, -2}, {7, 14, 1},
    {7, 155, -1}, {7, 156, 1}, {7, 157, -1}, {8, 15, 1},
    {8, 83, 1}, {9, 17, 1}, {9, 84, 1}, {10, 19, 1},
    {10, 137, 1}, {10, 157, -1}, {11, 20, 1}, {11, 137, 1},
    {11, 148, -1}, {11, 157, -1}, {12, 21, 1}, {12, 75, 1},
    {12, 156, 1}, {12, 157, -1}, {13, 24, 1}, {13, 142, -1},
    {13, 158, -1}, {13, 159, 1}, {14, 27, 1}, {14, 136, -1},
    {15, 28, 1}, {15, 156, 1}, {15, 157, -1}, {16, 31, 1},
    {16, 157, 1}, {16, 158, -1}, {17, 32, 1}, {17, 156, -1},
    {17, 157, 3}, {17, 158, -1}, {17, 159, 1}, {18, 37, 1},
    {18, 55, 1}, {18, 156, 1}, {18, 157, -1}, {19, 38, 1},
    {19, 82, -1}, {20, 39, 1}, {20, 156, -1}, {20, 157, 1},
    {21, 40, 1}, {21, 56, 1}, {22, 44, 1}, {22, 156, 1},
    {22, 157, -1}, {22, 158, 1}, {22, 159, -1}, {23, 45, 1},
    {23, 155, 1}, {24, 46, 1}, {24, 147, 1}, {25, 48, 1},
    {25, 52, 1}, {25, 156, 1}, {25, 157, -1}, {26, 49, 1},
    {26, 142, 1}, {26, 156, 1}, {26, 157, -1}, {27, 50, 1},
    {27, 145, 1}, {28, 57, 1}, {28, 142, 1}, {29, 59, 3},
    {29, 134, -3}, {29, 157, 4}, {29, 158, -2}, {29, 159, 5},
    {30, 60, 1}, {30, 98, 1}, {31, 63, 1}, {31, 157, -1},
    {32, 65, 1}, {32, 146, -1}, {32, 148, 1}, {32, 156, 1},
    {32, 157, -1}, {33, 66, 1}, {33, 72, 1}, {33, 130, 1},
    {34, 68, 1}, {34, 155, 1}, {35, 69, 1}, {35, 110, -1},
    {35, 137, 1}, {35, 148, -1}, {35, 149, 1}, {35, 156, -1},
    {35, 157, 1}, {36, 73, 1}, {36, 156, 1}, {36, 157, -1},
    {37, 74, 1}, {37, 146, 1}, {37, 148, -1}, {38, 76, 1},
    {38, 156, 1}, {38, 157, -1}, {39, 78, 1}, {39, 156, 1},
    {39, 157, -2}, {40, 79, 1}, {40, 81, 1}, {40, 156, 1},
    {40, 157, -1}, {41, 80, 1}, {41, 136, 1}, {41, 156, 1},
    {41, 157, -1}, {42, 86, 1}, {42, 99, 1}, {42, 148, 1},
    {42, 156, 1}, {43, 87, 1}, {43, 148, -1}, {44, 88, 1},
    {44, 156, 1}, {45, 89, 1}, {45, 136, -1}, {45, 137, 1},
    {45, 156, -1}, {46, 90, 1}, {46, 134, 1}, {47, 91, 1},
    {47, 110, 1}, {48, 92, 1}, {48, 157, 1}, {49, 94, 1},
    {49, 136, 1}, {49, 137, -1}, {49, 156, 1}, {50, 95, 1},
    {50, 137, 1}, {50, 156, 1}, {50, 157, -1}, {51, 96, 1},
    {51, 107, 1}, {51, 156, 1}, {51, 157, -1}, {52, 100, 1},
    {52, 156, 1}, {52, 157, -1}, {53, 101, 1}, {53, 151, 1},
    {53, 155, -1}, {54, 102, 1}, {54, 148, 1}, {54, 157, 1},
    {55, 103, 1}, {55, 136, 1}, {55, 137, -1}, {55, 157, 1},
    {56, 105, 1}, {56, 145, 1}, {56, 156, 1}, {56, 157, -1},
    {57, 106, 1}, {57, 137, -1}, {57, 148, 1}, {57, 156, 1},
    {57, 157, -1}, {58, 109, 1}, {58, 156, 1}, {58, 157, -2},
    {59, 111, 1}, {59, 153, 1}, {60, 112, 1}, {60, 156, -1},
    {60, 157, 1}, {61, 114, 1}, {61, 133, 1}, {61, 156, 1},
    {61, 157, -1}, {62, 115, 1}, {62, 156, 1}, {63, 116, 1},
    {63, 136, 1}, {64, 118, 1}, {64, 156, -1}, {65, 119, 1},
    {65, 157, 1}, {66, 120, 1}, {66, 156, 1}, {66, 157, -1},
    {67, 121, 1}, {67, 134, 1}, {67, 156, -1}, {68, 122, 1},
    {68, 148, 1}, {68, 156, 1}, {68, 157, -1}, {69, 123, 1},
    {69, 155, 1}, {69, 156, 1}, {69, 157, -1}, {70, 124, 1},
    {70, 156, -1}, {70, 157, 1}, {71, 125, 1}, {71, 148, 1},
    {72, 126, 1}, {72, 146, 1}, {73, 127, 1}, {73, 157, -1},
    {74, 128, 1}, {74, 142, -1}, {74, 156, 1}, {75, 129, 1},
    {75, 156, 1}, {75, 157, -1}, {76, 131, 1}, {76, 148, 1},
    {76, 156, 1}, {77, 135, 1}, {77, 145, 1}, {77, 156, 1},
    {77, 157, -1}, {78, 141, 1}, {78, 157, 1}, {79, 143, 1},
    {79, 155, 1}, {79, 156, -1}, {79, 157, 1}, {80, 150, 1},
    {80, 156, 2}, {80, 157, -2}, {81, 152, 1}, {81, 156, 1},
    {81, 157, -1}, {82, 160, 1}, {83, 31, 1}, {83, 41, 1},
    {83, 44, 1}, {83, 63, 1}, {83, 89, -1}, {83, 95, 1},
    {83, 123, -1}, {83, 131, -1}, {83, 154, 1}, {83, 159, -1},
    {83, 160, 1},
};
constexpr CorpusHnfEntry corpus_hnf_refresh_2[] = {
    {0, 0, 1}, {0, 159, 12}, {1, 1, 1}, {1, 159, 11},
    {2, 2, 1}, {2, 159, 14}, {3, 3, 1}, {3, 159, 10},
    {4, 4, 1}, {4, 159, 5}, {5, 5, 1}, {5, 159, 8},
    {6, 6, 1}, {6, 159, 2}, {7, 7, 1}, {7, 159, 1},
    {8, 8, 1}, {8, 159, 13}, {9, 9, 1}, {9, 159, 7},
    {10, 10, 1}, {10, 159, 15}, {11, 11, 1}, {11, 159, 1},
    {12, 12, 1}, {12, 159, 9}, {13, 13, 1}, {13, 159, 3},
    {14, 14, 1}, {14, 159, 5}, {15, 15, 1}, {15, 159, 6},
    {16, 16, 1}, {16, 159, 4}, {17, 17, 1}, {17, 159, 11},
    {18, 18, 1}, {18, 159, 7}, {19, 19, 1}, {19, 159, 4},
    {20, 20, 1}, {20, 159, 13}, {21, 21, 1}, {21, 159, 14},
    {22, 22, 1}, {22, 159, 2}, {23, 23, 1}, {23, 159, 12},
    {24, 24, 1}, {24, 159, 9}, {25, 25, 1}, {25, 159, 15},
    {26, 27, 1}, {26, 159, 6}, {27, 28, 1}, {27, 159, 3},
    {28, 29, 1}, {28, 159, 16}, {29, 30, 1}, {29, 159, 10},
    {30, 31, 1}, {30, 159, 16}, {31, 32, 1}, {31, 159, 2},
    {32, 33, 1}, {32, 159, 1}, {33, 34, 1}, {33, 159, 9},
    {34, 37, 1}, {34, 159, 13}, {35, 38, 1}, {35, 159, 13},
    {36, 39, 1}, {36, 159, 14}, {37, 40, 1}, {37, 159, 1},
    {38, 41, 1}, {38, 159, 12}, {39, 44, 1}, {39, 159, 8},
    {40, 45, 1}, {40, 159, 15}, {41, 46, 1}, {41, 159, 2},
    {42, 48, 1}, {42, 159, 9}, {43, 49, 1}, {43, 159, 6},
    {44, 50, 1}, {44, 159, 4}, {45, 51, 1}, {45, 159, 11},
    {46, 52, 1}, {46, 159, 11}, {47, 53, 1}, {47, 159, 6},
    {48, 54, 1}, {48, 159, 10}, {49, 55, 1}, {49, 159, 7},
    {50, 56, 1}, {50, 159, 16}, {51, 57, 1}, {51, 159, 3},
    {52, 58, 1}, {52, 159, 15}, {53, 59, 1}, {53, 159, 14},
    {54, 60, 1}, {54, 159, 7}, {55, 61, 1}, {55, 159, 8},
    {56, 62, 1}, {56, 159, 3}, {57, 63, 1}, {57, 159, 12},
    {58, 64, 1}, {58, 159, 8}, {59, 65, 1}, {59, 159, 9},
    {60, 66, 1}, {60, 159, 14}, {61, 67, 1}, {61, 159, 1},
    {62, 68, 1}, {62, 159, 15}, {63, 69, 1}, {63, 159, 10},
    {64, 70, 1}, {64, 159, 2}, {65, 71, 1}, {65, 159, 16},
    {66, 72, 1}, {66, 159, 2}, {67, 73, 1}, {67, 159, 3},
    {68, 74, 1}, {68, 159, 11}, {69, 75, 1}, {69, 159, 6},
    {70, 76, 1}, {70, 159, 3}, {71, 77, 1}, {71, 159, 12},
    {72, 78, 1}, {72, 159, 15}, {73, 79, 1}, {73, 159, 16},
    {74, 80, 1}, {74, 159, 14}, {75, 81, 1}, {75, 159, 4},
    {76, 82, 1}, {76, 159, 13}, {77, 83, 1}, {77, 159, 11},
    {78, 84, 1}, {78, 159, 6}, {79, 85, 1}, {79, 159, 5},
    {80, 86, 1}, {80, 159, 1}, {81, 87, 1}, {81, 159, 9},
    {82, 88, 1}, {82, 159, 8}, {83, 89, 1}, {83, 159, 7},
    {84, 90, 1}, {84, 159, 13}, {85, 91, 1}, {85, 159, 4},
    {86, 92, 1}, {86, 159, 5}, {87, 93, 1}, {87, 159, 7},
    {88, 94, 1}, {88, 159, 10}, {89, 95, 1}, {89, 159, 12},
    {90, 96, 1}, {90, 159, 4}, {91, 97, 1}, {91, 159, 1},
    {92, 98, 1}, {92, 159, 10}, {93, 99, 1}, {93, 159, 15},
    {94, 100, 1}, {94, 159, 3}, {95, 101, 1}, {95, 159, 9},
    {96, 102, 1}, {96, 159, 13}, {97, 103, 1}, {97, 159, 7},
    {98, 104, 1}, {98, 159, 4}, {99, 105, 1}, {99, 159, 7},
    {100, 106, 1}, {100, 159, 2}, {101, 107, 1}, {101, 159, 16},
    {102, 108, 1}, {102, 159, 6}, {103, 109, 1}, {103, 159, 15},
    {104, 110, 1}, {104, 159, 13}, {105, 111, 1}, {105, 159, 6},
    {106, 112, 1}, {106, 159, 14}, {107, 113, 1}, {107, 159, 12},
    {108, 114, 1}, {108, 159, 10}, {109, 115, 1}, {109, 159, 8},
    {110, 116, 1}, {110, 159, 11}, {111, 117, 1}, {111, 159, 16},
    {112, 118, 1}, {112, 159, 9}, {113, 119, 1}, {113, 159, 5},
    {114, 120, 1}, {114, 159, 3}, {115, 121, 1}, {115, 159, 5},
    {116, 122, 1}, {116, 159, 11}, {117, 123, 1}, {117, 159, 1},
    {118, 124, 1}, {118, 159, 14}, {119, 125, 1}, {119, 159, 8},
    {120, 126, 1}, {120, 159, 2}, {121, 127, 1}, {121, 159, 12},
    {122, 128, 1}, {122, 159, 5}, {123, 129, 1}, {123, 159, 3},
    {124, 130, 1}, {124, 159, 1}, {125, 131, 1}, {125, 159, 16},
    {126, 132, 1}, {126, 159, 2}, {127, 133, 1}, {127, 159, 10},
    {128, 134, 1}, {128, 159, 4}, {129, 135, 1}, {129, 159, 7},
    {130, 136, 1}, {130, 159, 6}, {131, 137, 1}, {131, 159, 8},
    {132, 138, 1}, {132, 159, 8}, {133, 139, 1}, {133, 159, 13},
    {134, 140, 1}, {134, 159, 7}, {135, 141, 1}, {135, 159, 5},
    {136, 142, 1}, {136, 159, 14}, {137, 143, 1}, {137, 159, 12},
    {138, 144, 1}, {138, 159, 4}, {139, 145, 1}, {139, 159, 13},
    {140, 146, 1}, {140, 159, 15}, {141, 147, 1}, {141, 159, 15},
    {142, 148, 1}, {142, 159, 9}, {143, 149, 1}, {143, 159, 1},
    {144, 150, 1}, {144, 159, 6}, {145, 151, 1}, {145, 159, 10},
    {146, 152, 1}, {146, 159, 3}, {147, 153, 1}, {147, 159, 11},
    {148, 154, 1}, {148, 159, 14}, {149, 155, 1}, {149, 159, 2},
    {150, 156, 1}, {150, 159, 9}, {151, 157, 1}, {151, 159, 12},
    {152, 158, 1}, {152, 159, 11}, {153, 159, 17}, {154, 160, 1},
    {155, 26, 1}, {155, 98, -1}, {155, 121, -2}, {155, 127, 1},
    {155, 160, 3},
};
struct CorpusHnfRefresh { const CorpusHnfEntry* entries; std::size_t count; slong old_rows; slong rows; };
constexpr CorpusHnfRefresh corpus_hnf_refreshes[] = {
    {corpus_hnf_refresh_0, 243, 73, 74},
    {corpus_hnf_refresh_1, 269, 83, 84},
    {corpus_hnf_refresh_2, 313, 155, 156},
};

// Matrices taken from real refreshes of the cyclotomic_64_phi32 class-unit
// run: the incremental update must agree with fmpz_mat_hnf of all rows.
int test_hnf_add_rows_corpus_refreshes() {
    for (const CorpusHnfRefresh& refresh : corpus_hnf_refreshes) {
        silex::flint::FmpzMat matrix(refresh.rows, 161);
        for (std::size_t k = 0; k < refresh.count; ++k) {
            const CorpusHnfEntry& entry = refresh.entries[k];
            set_entry_si(matrix, entry.row, entry.col, entry.value);
        }
        if (!hnf_add_rows_matches_full_hnf(
                    silex::flint::FmpzMatConstRef(matrix), refresh.old_rows)) {
            return 1;
        }
    }
    return 0;
}

// Relation-module shape: 161 columns, retained HNFs of rank up to 160 (the
// last column is never used, so the lattice stays rank deficient), and a few
// new rows with entries of at most 5 bits, including dependent rows.
int test_hnf_add_rows_relation_shape() {
    constexpr slong cols = 161;
    constexpr slong rows = 182;
    ulong state = 442UL;
    silex::flint::FmpzMat matrix(rows, cols);
    for (slong i = 0; i < rows; ++i) {
        if (i % 9 == 8) {
            // Dependent rows in the lattice of earlier rows.
            set_combination_row(matrix, i, i - 1, 2, i - 4, -1);
        } else {
            fill_relation_like_row(matrix, i, 0, 15, state);
            // Drop the last column, then give the row an entry of its own
            // column, mostly 1 as for relation rows, so the rank grows with
            // the row count.
            ::fmpz_zero(::fmpz_mat_entry(matrix.raw(), i, cols - 1));
            const slong own = (i * 37) % (cols - 1);
            const ulong pick = next_random(state) % 8;
            ::fmpz_set_si(::fmpz_mat_entry(matrix.raw(), i, own),
                          pick < 5 ? 1 : static_cast<slong>(pick) - 3);
        }
    }
    // Row 100 becomes row 30 + row 31 and rows 30 and 31 are then tripled,
    // so row 100 is rationally dependent but refines the lattice.
    set_combination_row(matrix, 100, 30, 1, 31, 1);
    for (slong j = 0; j < cols; ++j) {
        ::fmpz_mul_si(::fmpz_mat_entry(matrix.raw(), 30, j),
                      ::fmpz_mat_entry(matrix.raw(), 30, j), 3);
        ::fmpz_mul_si(::fmpz_mat_entry(matrix.raw(), 31, j),
                      ::fmpz_mat_entry(matrix.raw(), 31, j), 3);
    }

    for (const slong old_rows : {slong{80}, slong{99}, slong{176},
                                 slong{178}}) {
        for (const slong new_rows : {slong{1}, slong{3}}) {
            silex::flint::FmpzMatConstWindow window(
                    matrix, 0, 0, old_rows + new_rows, cols);
            if (!hnf_add_rows_matches_full_hnf(window.const_ref(),
                                               old_rows)) {
                return 1;
            }
        }
    }
    return 0;
}

#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
struct HnfRouteCounters {
    slong add_rows = 0;
};

void hnf_route_profile_callback(void* user,
                                silex::DiagnosticsModule,
                                silex::ProfileEvent event,
                                const char*,
                                const char* label) noexcept {
    if (event == silex::ProfileEvent::begin_scope && label != nullptr &&
        std::strcmp(label, "fmpz_smat.hnf_context.refresh.hnf.add_rows") ==
                0) {
        ++static_cast<HnfRouteCounters*>(user)->add_rows;
    }
}
#endif

// The context refreshes the retained HNF incrementally after a few accepted
// or deferred rows; every refreshed HNF must be the canonical HNF of all rows
// added so far, through partial rank, refining dependent rows and full rank.
int test_hnf_context_refresh_adds_rows_to_retained_hnf() {
    constexpr slong cols = 40;
    constexpr slong rows = 72;
    ulong state = 447UL;
    silex::flint::FmpzMat matrix(rows, cols);
    for (slong i = 0; i < rows; ++i) {
        // Supports widen leftwards, so new pivots land in front of old ones.
        const slong first = cols - 1 - (i * cols) / rows;
        if (i % 7 == 6) {
            set_combination_row(matrix, i, i - 1, 3, i - 3, -2);
        } else if (i % 11 == 10) {
            // Half of the doubled row i - 2: refines the lattice.
            for (slong j = 0; j < cols; ++j) {
                ::fmpz_divexact_ui(::fmpz_mat_entry(matrix.raw(), i, j),
                                   ::fmpz_mat_entry(matrix.raw(), i - 2, j),
                                   2);
            }
        } else {
            fill_relation_like_row(matrix, i, first, 15, state);
        }
        if (i % 11 == 8) {
            for (slong j = 0; j < cols; ++j) {
                ::fmpz_mul_2exp(::fmpz_mat_entry(matrix.raw(), i, j),
                                ::fmpz_mat_entry(matrix.raw(), i, j), 1);
            }
        }
    }

    bool independent = false;
    bool reduces_to_zero = false;
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    HnfRouteCounters counters;
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(
                    silex::DiagnosticsModule::class_group),
            hnf_route_profile_callback, &counters);
#endif
    for (const ulong prime : {ulong{0}, ulong{3}, ulong{65537}}) {
        for (const bool defer : {false, true}) {
            silex::fmpz_smat::HnfContext context;
            if (!context.reset(cols, prime)) {
                return 1;
            }
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
            context.set_diagnostics(&diagnostics);
#endif
            for (slong i = 0; i < rows; ++i) {
                const bool ok =
                        defer ? context.add_fmpz_mat_row_defer_dependent(
                                        &independent,
                                        silex::flint::FmpzMatConstRef(matrix),
                                        i)
                              : context.add_fmpz_mat_row(
                                        &independent,
                                        silex::flint::FmpzMatConstRef(matrix),
                                        i);
                if (!ok) {
                    return 1;
                }
                // Refresh after most rows, as the relation search does, and
                // let a few rows accumulate before some refreshes.
                if (i % 5 != 4 &&
                    !context.fmpz_mat_row_reduces_to_zero(
                            reduces_to_zero,
                            silex::flint::FmpzMatConstRef(matrix), i)) {
                    return 1;
                }
                if (!reduces_to_zero && i % 5 != 4) {
                    return 1;
                }
                if ((i % 6 == 5 || i == rows - 1) &&
                    !hnf_context_matches_prefix(context, matrix, i + 1)) {
                    return 1;
                }
            }
        }
    }
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The incremental route must have been attempted on refreshes with a
    // retained HNF (the scope is entered per attempt; a failed validation
    // falls back to fmpz_mat_hnf with the same result).
    if (counters.add_rows == 0) {
        return 1;
    }
#endif
    return 0;
}

int test_native_cpp_raii_call_sites() {
    static_assert(!std::is_copy_constructible_v<silex::fmpz_smat::SparseRow>);
    static_assert(!std::is_copy_assignable_v<silex::fmpz_smat::SparseRow>);
    static_assert(!std::is_copy_constructible_v<silex::fmpz_smat::SparseMat>);
    static_assert(!std::is_copy_assignable_v<silex::fmpz_smat::SparseMat>);
    static_assert(!std::is_copy_constructible_v<silex::fmpz_smat::ModRankContext>);
    static_assert(!std::is_copy_assignable_v<silex::fmpz_smat::ModRankContext>);
    static_assert(!std::is_copy_constructible_v<silex::fmpz_smat::HnfContext>);
    static_assert(!std::is_copy_assignable_v<silex::fmpz_smat::HnfContext>);

    silex::fmpz_smat::SparseMat matrix(0);
    silex::fmpz_smat::SparseRow row;
    silex::flint::Fmpz value;
    silex::flint::Fmpz index;
    silex::flint::FmpzMat dense(3, 3);
    silex::flint::FmpzMat exported(3, 3);
    silex::flint::FmpzMat right(3, 2);
    silex::flint::FmpzMat product(3, 2);
    silex::flint::FmpzMat hnf(3, 3);
    silex::flint::NmodMat mod(3, 3, 5);
    slong rank = -1;
    bool independent = false;

    set_entry_si(dense.raw(), 0, 0, 2);
    set_entry_si(dense.raw(), 1, 1, 3);
    set_entry_si(dense.raw(), 2, 0, 1);
    set_entry_si(dense.raw(), 2, 1, 1);
    set_entry_si(dense.raw(), 2, 2, 5);

    matrix.set_fmpz_mat(dense);
    silex::fmpz_smat::SparseMat moved(std::move(matrix));
    silex::fmpz_smat::SparseMat assigned(0);
    assigned = std::move(moved);

    if (assigned.nrows() != 3 || assigned.ncols() != 3 || assigned.nnz() != 5) {
        return 1;
    }
    assigned.get_fmpz_mat(exported);
    if (fmpz_mat_equal(exported.raw(), dense.raw()) == 0) {
        return 1;
    }

    assigned.get_row(row, 2);
    row.get_entry(value, 2);
    if (!fmpz_is_si(value.raw(), 5)) {
        return 1;
    }

    assigned.get_nmod_mat(mod);
    if (nmod_mat_entry(mod.raw(), 0, 0) != 2 ||
        nmod_mat_entry(mod.raw(), 1, 1) != 3 ||
        nmod_mat_entry(mod.raw(), 2, 0) != 1 ||
        nmod_mat_entry(mod.raw(), 2, 2) != 0) {
        return 1;
    }

    rank = 99;
    if (!assigned.rank_mod_prime_ui(&rank, 7) || rank != 3 ||
        assigned.rank_mod_prime_ui(&rank, 4) || rank != 3) {
        return 1;
    }

    set_entry_si(right.raw(), 0, 0, 1);
    set_entry_si(right.raw(), 0, 1, 2);
    set_entry_si(right.raw(), 1, 1, 1);
    set_entry_si(right.raw(), 2, 0, -1);
    set_entry_si(right.raw(), 2, 1, 3);
    assigned.mul_fmpz_mat(product, right);
    if (!entry_is_si(product.raw(), 0, 0, 2) ||
        !entry_is_si(product.raw(), 0, 1, 4) ||
        !entry_is_si(product.raw(), 1, 0, 0) ||
        !entry_is_si(product.raw(), 1, 1, 3) ||
        !entry_is_si(product.raw(), 2, 0, -4) ||
        !entry_is_si(product.raw(), 2, 1, 18)) {
        return 1;
    }

    silex::fmpz_smat::ModRankContext mod_rank;
    if (!mod_rank.set_prime_ui(3, 7)) {
        return 1;
    }
    for (slong i = 0; i < 3; ++i) {
        independent = false;
        if (!mod_rank.add_fmpz_mat_row(&independent, dense, i) ||
            !independent || mod_rank.rank() != i + 1) {
            return 1;
        }
    }
    silex::fmpz_smat::ModRankContext copied_mod_rank;
    if (!copied_mod_rank.set(mod_rank) ||
        copied_mod_rank.ambient_dim() != mod_rank.ambient_dim() ||
        copied_mod_rank.prime() != mod_rank.prime() ||
        copied_mod_rank.rank() != mod_rank.rank()) {
        return 1;
    }
    independent = true;
    if (!copied_mod_rank.add_fmpz_mat_row(&independent, dense, 0) ||
        independent || copied_mod_rank.rank() != mod_rank.rank()) {
        return 1;
    }

    silex::fmpz_smat::HnfContext hnf_context;
    if (!hnf_context.reset(3, 0)) {
        return 1;
    }
    for (slong i = 0; i < 3; ++i) {
        if (!hnf_context.add_fmpz_mat_row(&independent, dense, i)) {
            return 1;
        }
    }
    if (hnf_context.rank() != 3 || !hnf_context.get_hnf(hnf) ||
        !hnf_context.full_rank_index(index) ||
        !fmpz_is_si(index.raw(), 30)) {
        return 1;
    }

    return 0;
}

int check_square_abs_det(std::initializer_list<slong> entries,
                         slong n,
                         silex::detail::SquareAbsDetRoute expected_route,
                         slong expected) {
    silex::flint::FmpzMat matrix(n, n);
    slong k = 0;
    for (slong value : entries) {
        fmpz_set_si(fmpz_mat_entry(matrix.raw(), k / n, k % n), value);
        ++k;
    }
    silex::flint::Fmpz out;
    fmpz_set_si(out.raw(), -99);
    const silex::detail::SquareAbsDetRoute route =
            silex::detail::fmpz_mat_square_abs_det(out, matrix);
    return route == expected_route && fmpz_is_si(out.raw(), expected) ? 0 : 1;
}

int test_square_abs_det_routes() {
    using silex::detail::SquareAbsDetRoute;
    // Upper triangular, with a negative diagonal entry: |(-2) * 3 * 5|.
    if (check_square_abs_det({-2, 5, 7, 0, 3, 1, 0, 0, 5}, 3,
                             SquareAbsDetRoute::diagonal_product, 30) != 0) {
        return 1;
    }
    // Singular upper triangular.
    if (check_square_abs_det({4, 1, 0, 0}, 2,
                             SquareAbsDetRoute::diagonal_product, 0) != 0) {
        return 1;
    }
    // Empty matrix: the empty product.
    if (check_square_abs_det({}, 0, SquareAbsDetRoute::diagonal_product, 1) !=
        0) {
        return 1;
    }
    // Not upper triangular: det [[1, 2], [3, 4]] = -2.
    if (check_square_abs_det({1, 2, 3, 4}, 2,
                             SquareAbsDetRoute::determinant, 2) != 0) {
        return 1;
    }
    // Lower triangular is not taken by the diagonal route; its only
    // nonzero subdiagonal entry is the last one scanned.
    if (check_square_abs_det({2, 0, 0, 0, 3, 0, 0, 7, -5}, 3,
                             SquareAbsDetRoute::determinant, 30) != 0) {
        return 1;
    }
    return 0;
}

}  // namespace

int main() {

    return test_sparse_row_basic() != 0 || test_sparse_matrix_basic() != 0 ||
                   test_sparse_row_dense_equality_is_canonical() != 0 ||
                   test_sparse_matrix_kernels() != 0 ||
                   test_sparse_matrix_nmod() != 0 || test_mod_rank_context() != 0 ||
                   test_mod_rank_context_row_is_independent() != 0 ||
                   test_hnf_context() != 0 ||
                   test_hnf_context_transform() != 0 ||
                   test_hnf_context_index_before_transform() != 0 ||
                   test_hnf_context_incremental_refresh() != 0 ||
                   test_hnf_context_deferred_dependent_refresh() != 0 ||
                   test_hnf_add_rows_random_small() != 0 ||
                   test_hnf_add_rows_relation_shape() != 0 ||
                   test_hnf_add_rows_corpus_refreshes() != 0 ||
                   test_hnf_context_refresh_adds_rows_to_retained_hnf() != 0 ||
                   test_native_cpp_raii_call_sites() != 0 ||
                   test_square_abs_det_routes() != 0
               ? 1
               : 0;
}
