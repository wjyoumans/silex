#include <silex/fmpz_smat.hpp>

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/nmod_mat.hpp>

#include <flint/fmpz_mat.h>

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

}  // namespace

int main() {
    return test_sparse_row_basic() != 0 || test_sparse_matrix_basic() != 0 ||
                   test_sparse_matrix_kernels() != 0 ||
                   test_sparse_matrix_nmod() != 0 || test_mod_rank_context() != 0 ||
                   test_mod_rank_context_row_is_independent() != 0 ||
                   test_hnf_context() != 0 ||
                   test_hnf_context_transform() != 0 ||
                   test_hnf_context_index_before_transform() != 0 ||
                   test_hnf_context_incremental_refresh() != 0 ||
                   test_hnf_context_deferred_dependent_refresh() != 0 ||
                   test_native_cpp_raii_call_sites() != 0
               ? 1
               : 0;
}
