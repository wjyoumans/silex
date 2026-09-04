#include <silex/fmpz_smat.hpp>

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/nmod_mat.hpp>

#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>

#include <vector>

namespace silex::fmpz_smat {
namespace {

slong next_alloc(slong old_alloc, slong min_alloc) noexcept {
    slong new_alloc = (old_alloc > 0) ? 2 * old_alloc : 1;
    if (new_alloc < min_alloc) {
        new_alloc = min_alloc;
    }
    return new_alloc;
}

void nmod_row_clear(ulong* row, slong n) noexcept {
    for (slong i = 0; i < n; ++i) {
        row[i] = 0;
    }
}

slong nmod_row_first_nonzero(const ulong* row, slong n) noexcept {
    for (slong i = 0; i < n; ++i) {
        if (row[i] != 0) {
            return i;
        }
    }
    return -1;
}

void nmod_row_submul(ulong* row,
        const ulong* basis_row,
        slong n,
        ulong c,
        nmod_t mod) noexcept {
    if (c == 0) {
        return;
    }
    for (slong i = 0; i < n; ++i) {
        if (basis_row[i] != 0) {
            const ulong tmp = nmod_mul(c, basis_row[i], mod);
            row[i] = nmod_sub(row[i], tmp, mod);
        }
    }
}

void nmod_row_normalize(ulong* row, slong n, ulong c, nmod_t mod) noexcept {
    for (slong i = 0; i < n; ++i) {
        row[i] = nmod_mul(row[i], c, mod);
    }
}

slong trimmed_nonzero_rows(const fmpz_mat_t matrix) noexcept {
    slong rows = fmpz_mat_nrows(matrix);
    while (rows > 0 && fmpz_mat_is_zero_row(matrix, rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (fmpz_mat_is_zero_row(matrix, i) != 0) {
            return -1;
        }
    }
    return rows;
}

slong fmpz_mat_row_first_nonzero(const fmpz_mat_t matrix,
                                 slong row) noexcept {
    const slong cols = fmpz_mat_ncols(matrix);
    for (slong j = 0; j < cols; ++j) {
        if (::fmpz_is_zero(::fmpz_mat_entry(matrix, row, j)) == 0) {
            return j;
        }
    }
    return -1;
}

bool row_reduces_to_zero_mod_hnf(bool& reduces_to_zero,
                                 const SparseRow& row,
                                 const fmpz_mat_t hnf,
                                 slong rank,
                                 slong ambient_dim) noexcept {
    reduces_to_zero = false;
    if (fmpz_mat_nrows(hnf) != rank ||
        fmpz_mat_ncols(hnf) != ambient_dim ||
        !row.fits_columns(ambient_dim)) {
        return false;
    }

    flint::FmpzMat remainder(1, ambient_dim);
    row.get_fmpz_mat_row(flint::FmpzMatRef(remainder), 0);

    flint::Fmpz quotient;
    for (slong i = 0; i < rank; ++i) {
        const slong pivot_col = fmpz_mat_row_first_nonzero(hnf, i);
        if (pivot_col < 0) {
            return false;
        }

        const fmpz* pivot = ::fmpz_mat_entry(hnf, i, pivot_col);
        fmpz* entry = ::fmpz_mat_entry(remainder.raw(), 0, pivot_col);
        if (::fmpz_is_zero(entry) != 0) {
            continue;
        }
        if (::fmpz_divisible(entry, pivot) == 0) {
            return true;
        }

        ::fmpz_fdiv_q(quotient.raw(), entry, pivot);
        for (slong j = pivot_col; j < ambient_dim; ++j) {
            fmpz* remainder_entry = ::fmpz_mat_entry(remainder.raw(), 0, j);
            const fmpz* hnf_entry = ::fmpz_mat_entry(hnf, i, j);
            if (::fmpz_is_zero(hnf_entry) == 0) {
                ::fmpz_submul(remainder_entry, quotient.raw(), hnf_entry);
            }
        }
    }

    reduces_to_zero = fmpz_mat_is_zero_row(remainder.raw(), 0) != 0;
    return true;
}

}  // namespace

void SparseRow::swap(SparseRow& other) noexcept {
    cols_.swap(other.cols_);
    vals_.swap(other.vals_);
}

bool SparseRow::set(const SparseRow& other) noexcept {
    if (this == &other) {
        return true;
    }
    const slong other_length = other.length();
    if (!assure_alloc(other_length)) {
        return false;
    }
    cols_.resize(static_cast<std::size_t>(other_length));
    for (slong i = 0; i < other_length; ++i) {
        cols_[static_cast<std::size_t>(i)] =
                other.cols_[static_cast<std::size_t>(i)];
        fmpz_set(vals_.data() + i, other.vals_.data() + i);
    }
    return true;
}

std::span<const slong> SparseRow::columns() const noexcept {
    return std::span<const slong>(cols_.data(), cols_.size());
}

flint::FmpzVecConstRef SparseRow::values() const noexcept {
    return flint::FmpzVecConstRef(vals_.data(), length());
}

slong SparseRow::col(slong i) const {
    if (i < 0 || i >= length()) {
        flint_throw(FLINT_ERROR, "SparseRow::col: index out of range\n");
    }
    return cols_[static_cast<std::size_t>(i)];
}

void SparseRow::get_value(flint::FmpzRef out, slong i) const {
    if (i < 0 || i >= length()) {
        flint_throw(FLINT_ERROR, "SparseRow::get_value: index out of range\n");
    }
    flint::fmpz_set(out, vals_.data() + i);
}

void SparseRow::get_entry(flint::FmpzRef out, slong col) const {
    if (col < 0) {
        flint_throw(FLINT_ERROR, "SparseRow::get_entry: column index out of range\n");
    }
    bool found = false;
    const slong pos = find(col, &found);
    if (found) {
        flint::fmpz_set(out, vals_.data() + pos);
    } else {
        flint::fmpz_zero(out);
    }
}

bool SparseRow::set_entry(slong col, flint::FmpzConstRef value) noexcept {
    if (col < 0) {
        return false;
    }

    bool found = false;
    const slong pos = find(col, &found);
    const slong old_length = length();
    if (found) {
        if (flint::fmpz_is_zero(value)) {
            for (slong i = pos; i + 1 < old_length; ++i) {
                fmpz_set(vals_.data() + i, vals_.data() + i + 1);
            }
            cols_.erase(cols_.begin() + pos);
        } else {
            fmpz_set(vals_.data() + pos, value.raw());
        }
        return true;
    }

    if (flint::fmpz_is_zero(value)) {
        return true;
    }

    if (!assure_alloc(old_length + 1)) {
        return false;
    }
    for (slong i = old_length; i > pos; --i) {
        fmpz_set(vals_.data() + i, vals_.data() + i - 1);
    }
    cols_.insert(cols_.begin() + pos, col);
    fmpz_set(vals_.data() + pos, value.raw());
    return true;
}

bool SparseRow::append_entry(slong col, flint::FmpzConstRef value) noexcept {
    const slong old_length = length();
    if (col < 0 ||
        (old_length > 0 &&
         col <= cols_[static_cast<std::size_t>(old_length - 1)])) {
        return false;
    }
    if (flint::fmpz_is_zero(value)) {
        return true;
    }
    if (!assure_alloc(old_length + 1)) {
        return false;
    }
    cols_.push_back(col);
    fmpz_set(vals_.data() + old_length, value.raw());
    return true;
}

bool SparseRow::set_fmpz_mat_row(flint::FmpzMatConstRef matrix, slong row) noexcept {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix)) {
        return false;
    }

    SparseRow tmp;
    for (slong j = 0; j < flint::fmpz_mat_ncols(matrix); ++j) {
        if (!tmp.append_entry(j, flint::fmpz_mat_entry(matrix, row, j))) {
            return false;
        }
    }
    swap(tmp);
    return true;
}

void SparseRow::get_fmpz_mat_row(flint::FmpzMatRef matrix, slong row) const {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix)) {
        flint_throw(FLINT_ERROR, "SparseRow::get_fmpz_mat_row: row index out of range\n");
    }
    const slong ncols = flint::fmpz_mat_ncols(matrix);
    if (!fits(ncols)) {
        flint_throw(FLINT_ERROR, "SparseRow::get_fmpz_mat_row: output has too few columns\n");
    }

    for (slong j = 0; j < ncols; ++j) {
        flint::fmpz_zero(flint::fmpz_mat_entry(matrix, row, j));
    }
    for (slong j = 0; j < length(); ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(
                                matrix, row,
                                cols_[static_cast<std::size_t>(j)]),
                        vals_.data() + j);
    }
}

flint::FmpzMat SparseRow::to_fmpz_mat_row(slong ncols) const {
    if (!fits(ncols)) {
        flint_throw(FLINT_ERROR, "SparseRow::to_fmpz_mat_row: output has too few columns\n");
    }
    flint::FmpzMat matrix(1, ncols);
    get_fmpz_mat_row(flint::FmpzMatRef(matrix), 0);
    return matrix;
}

void SparseRow::dot_fmpz_mat_row(flint::FmpzRef out,
        flint::FmpzMatConstRef matrix,
        slong row) const {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix)) {
        flint_throw(FLINT_ERROR, "SparseRow::dot_fmpz_mat_row: row index out of range\n");
    }
    if (!fits(flint::fmpz_mat_ncols(matrix))) {
        flint_throw(FLINT_ERROR, "SparseRow::dot_fmpz_mat_row: input has too few columns\n");
    }

    flint::Fmpz acc;
    flint::Fmpz tmp;
    // Sparse/dense dot kernel: keep direct FLINT entry access visible.
    for (slong j = 0; j < length(); ++j) {
        fmpz_mul(tmp.raw(), vals_.data() + j,
                 fmpz_mat_entry(matrix.raw(), row,
                                cols_[static_cast<std::size_t>(j)]));
        fmpz_add(acc.raw(), acc.raw(), tmp.raw());
    }
    fmpz_set(out.raw(), acc.raw());
}

void SparseRow::get_nmod_mat_row(flint::NmodMatRef matrix, slong row) const {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::nmod_mat_nrows(matrix)) {
        flint_throw(FLINT_ERROR, "SparseRow::get_nmod_mat_row: row index out of range\n");
    }
    const slong ncols = flint::nmod_mat_ncols(matrix);
    if (!fits(ncols)) {
        flint_throw(FLINT_ERROR, "SparseRow::get_nmod_mat_row: output has too few columns\n");
    }

    for (slong j = 0; j < ncols; ++j) {
        flint::nmod_mat_set_entry(matrix, row, j, 0);
    }
    for (slong j = 0; j < length(); ++j) {
        flint::nmod_mat_set_entry(
                matrix, row, cols_[static_cast<std::size_t>(j)],
                fmpz_get_nmod(vals_.data() + j, matrix.raw()->mod));
    }
}

bool SparseRow::assure_alloc(slong alloc) noexcept {
    if (alloc <= vals_.length()) {
        return true;
    }

    const slong new_alloc = next_alloc(vals_.length(), alloc);
    cols_.reserve(static_cast<std::size_t>(new_alloc));
    flint::FmpzVec new_vals(new_alloc);
    for (slong i = 0; i < length(); ++i) {
        fmpz_set(new_vals.data() + i, vals_.data() + i);
    }
    vals_.swap(new_vals);
    return true;
}

slong SparseRow::find(slong col, bool* found) const noexcept {
    slong lo = 0;
    const slong row_length = length();
    slong hi = row_length;
    while (lo < hi) {
        const slong mid = lo + (hi - lo) / 2;
        if (cols_[static_cast<std::size_t>(mid)] < col) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    *found = (lo < row_length &&
              cols_[static_cast<std::size_t>(lo)] == col);
    return lo;
}

bool SparseRow::fits(slong ncols) const noexcept {
    const slong row_length = length();
    if (row_length == 0) {
        return true;
    }
    return cols_[static_cast<std::size_t>(row_length - 1)] < ncols;
}

SparseMat::SparseMat(slong ncols) : cols_(ncols) {
    if (ncols < 0) {
        flint_throw(FLINT_ERROR, "SparseMat::SparseMat: negative column count\n");
    }
}

void SparseMat::swap(SparseMat& other) noexcept {
    row_data_.swap(other.row_data_);
    std::swap(cols_, other.cols_);
    std::swap(nnz_, other.nnz_);
}

bool SparseMat::set(const SparseMat& other) noexcept {
    if (this == &other) {
        return true;
    }

    SparseMat tmp(other.cols_);
    tmp.row_data_.reserve(other.row_data_.size());
    for (slong i = 0; i < other.nrows(); ++i) {
        tmp.row_data_.emplace_back();
        if (!tmp.row_data_.back().set(
                    other.row_data_[static_cast<std::size_t>(i)])) {
            return false;
        }
    }
    tmp.nnz_ = other.nnz_;
    swap(tmp);
    return true;
}

bool SparseMat::append_row(const SparseRow& row) noexcept {
    if (!row.fits(cols_)) {
        return false;
    }
    row_data_.emplace_back();
    if (!row_data_.back().set(row)) {
        row_data_.pop_back();
        return false;
    }
    nnz_ += row.length();
    return true;
}

bool SparseMat::append_fmpz_mat_row(flint::FmpzMatConstRef matrix, slong row) noexcept {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != cols_) {
        return false;
    }

    SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }
    return append_row(sparse_row);
}

bool SparseMat::set_row(slong row_index, const SparseRow& row) noexcept {
    if (row_index < 0 || row_index >= nrows() || !row.fits(cols_)) {
        return false;
    }

    SparseRow tmp;
    if (!tmp.set(row)) {
        return false;
    }
    SparseRow& target = row_data_[static_cast<std::size_t>(row_index)];
    nnz_ -= target.length();
    target.swap(tmp);
    nnz_ += target.length();
    return true;
}

const SparseRow& SparseMat::row_ref(slong row_index) const {
    if (row_index < 0 || row_index >= nrows()) {
        flint_throw(FLINT_ERROR, "SparseMat::row_ref: row index out of range\n");
    }
    return row_data_[static_cast<std::size_t>(row_index)];
}

void SparseMat::get_row(SparseRow& row, slong row_index) const {
    if (row_index < 0 || row_index >= nrows()) {
        flint_throw(FLINT_ERROR, "SparseMat::get_row: row index out of range\n");
    }
    row.set(row_data_[static_cast<std::size_t>(row_index)]);
}

void SparseMat::get_entry(flint::FmpzRef out, slong row, slong col) const {
    if (row < 0 || row >= nrows() || col < 0 || col >= cols_) {
        flint_throw(FLINT_ERROR, "SparseMat::get_entry: index out of range\n");
    }
    row_data_[static_cast<std::size_t>(row)].get_entry(out, col);
}

void SparseMat::get_fmpz_mat(flint::FmpzMatRef matrix) const {
    if (matrix.raw() == nullptr || flint::fmpz_mat_nrows(matrix) != nrows() ||
        flint::fmpz_mat_ncols(matrix) != cols_) {
        flint_throw(FLINT_ERROR, "SparseMat::get_fmpz_mat: output has wrong dimensions\n");
    }

    flint::fmpz_mat_zero(matrix);
    for (slong i = 0; i < nrows(); ++i) {
        row_data_[static_cast<std::size_t>(i)].get_fmpz_mat_row(matrix, i);
    }
}

flint::FmpzMat SparseMat::to_fmpz_mat() const {
    flint::FmpzMat matrix(nrows(), cols_);
    get_fmpz_mat(matrix);
    return matrix;
}

void SparseMat::set_fmpz_mat(flint::FmpzMatConstRef matrix) noexcept {
    if (matrix.raw() == nullptr) {
        return;
    }
    SparseMat tmp(flint::fmpz_mat_ncols(matrix));
    for (slong i = 0; i < flint::fmpz_mat_nrows(matrix); ++i) {
        tmp.append_fmpz_mat_row(matrix, i);
    }
    swap(tmp);
}

void SparseMat::mul_fmpz_mat(flint::FmpzMatRef out,
        flint::FmpzMatConstRef right) const {
    if (out.raw() == nullptr || right.raw() == nullptr ||
        flint::fmpz_mat_nrows(right) != cols_ ||
        flint::fmpz_mat_nrows(out) != nrows() ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(right)) {
        flint_throw(FLINT_ERROR, "SparseMat::mul_fmpz_mat: incompatible dimensions\n");
    }

    if (out.raw() == right.raw()) {
        flint::FmpzMat tmp(nrows(), flint::fmpz_mat_ncols(right));
        mul_fmpz_mat(tmp, right);
        flint::fmpz_mat_swap(out, tmp);
        return;
    }

    flint::fmpz_mat_zero(out);
    flint::Fmpz tmp;
    // Sparse-times-dense kernel: raw entries avoid hiding the accumulation cost.
    for (slong i = 0; i < nrows(); ++i) {
        const SparseRow& sparse_row = row_data_[static_cast<std::size_t>(i)];
        for (slong j = 0; j < sparse_row.length(); ++j) {
            const slong col = sparse_row.cols_[static_cast<std::size_t>(j)];
            for (slong k = 0; k < flint::fmpz_mat_ncols(right); ++k) {
                fmpz_mul(tmp.raw(), sparse_row.vals_.data() + j,
                        fmpz_mat_entry(right.raw(), col, k));
                fmpz_add(fmpz_mat_entry(out.raw(), i, k),
                        fmpz_mat_entry(out.raw(), i, k), tmp.raw());
            }
        }
    }
}

void SparseMat::transpose(SparseMat& out) const noexcept {
    SparseMat tmp(nrows());
    SparseRow empty;
    for (slong i = 0; i < cols_; ++i) {
        tmp.append_row(empty);
    }

    for (slong i = 0; i < nrows(); ++i) {
        const SparseRow& sparse_row = row_data_[static_cast<std::size_t>(i)];
        for (slong j = 0; j < sparse_row.length(); ++j) {
            const slong col = sparse_row.cols_[static_cast<std::size_t>(j)];
            tmp.row_data_[static_cast<std::size_t>(col)].set_entry(
                    i, sparse_row.vals_.data() + j);
            ++tmp.nnz_;
        }
    }
    out.swap(tmp);
}

void SparseMat::get_nmod_mat(flint::NmodMatRef matrix) const {
    if (matrix.raw() == nullptr || flint::nmod_mat_nrows(matrix) != nrows() ||
        flint::nmod_mat_ncols(matrix) != cols_) {
        flint_throw(FLINT_ERROR, "SparseMat::get_nmod_mat: output has wrong dimensions\n");
    }

    flint::nmod_mat_zero(matrix);
    for (slong i = 0; i < nrows(); ++i) {
        row_data_[static_cast<std::size_t>(i)].get_nmod_mat_row(matrix, i);
    }
}

bool SparseMat::rank_mod_prime_ui(slong* rank, ulong p) const noexcept {
    if (rank == nullptr || !n_is_prime(p)) {
        return false;
    }

    flint::NmodMat matrix(nrows(), cols_, p);
    get_nmod_mat(matrix.raw());
    *rank = nmod_mat_rank(matrix.raw());
    return true;
}

ModRankContext::ModRankContext() noexcept = default;

ModRankContext::~ModRankContext() noexcept {
    ambient_dim_ = 0;
    rank_ = 0;
    prime_ = 0;
    is_set_ = false;
}

void ModRankContext::swap(ModRankContext& other) noexcept {
    basis_.swap(other.basis_);
    pivots_.swap(other.pivots_);
    std::swap(ambient_dim_, other.ambient_dim_);
    std::swap(rank_, other.rank_);
    std::swap(prime_, other.prime_);
    std::swap(is_set_, other.is_set_);
}

bool ModRankContext::set(const ModRankContext& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_set_) {
        ModRankContext empty;
        swap(empty);
        return true;
    }

    ModRankContext copy;
    if (!copy.set_prime_ui(other.ambient_dim_, other.prime_)) {
        return false;
    }
    ::nmod_mat_set(copy.basis_.raw(), other.basis_.raw());
    copy.pivots_ = other.pivots_;
    copy.rank_ = other.rank_;

    swap(copy);
    return true;
}

bool ModRankContext::set_prime_ui(slong ambient_dim, ulong p) noexcept {
    if (ambient_dim < 0 || !n_is_prime(p)) {
        return false;
    }

    flint::NmodMat basis(ambient_dim, ambient_dim, p);
    std::vector<slong> pivots(static_cast<std::size_t>(ambient_dim));

    basis_.swap(basis);

    pivots_.swap(pivots);
    ambient_dim_ = ambient_dim;
    rank_ = 0;
    prime_ = p;
    is_set_ = true;
    return true;
}

void ModRankContext::reset() noexcept {
    flint::nmod_mat_zero(basis_);
    rank_ = 0;
}

bool ModRankContext::add_row(bool* independent, const SparseRow& row) noexcept {
    if (independent == nullptr || !is_set_ || !row.fits_columns(ambient_dim_)) {
        return false;
    }

    const slong n = ambient_dim_;
    nmod_t mod;
    nmod_init(&mod, prime_);
    std::vector<ulong> vector(static_cast<std::size_t>(n));
    nmod_row_clear(vector.data(), n);

    for (slong i = 0; i < row.length(); ++i) {
        vector[static_cast<std::size_t>(
                row.cols_[static_cast<std::size_t>(i)])] =
                fmpz_get_nmod(row.vals_.data() + i, mod);
    }

    // Modular row-reduction kernel: raw row pointers feed FLINT nmod row ops.
    for (slong i = 0; i < rank_; ++i) {
        const slong pivot = pivots_[static_cast<std::size_t>(i)];
        const ulong c = vector[static_cast<std::size_t>(pivot)];
        if (c != 0) {
            ulong* const basis_row = nmod_mat_entry_ptr(basis_.raw(), i, 0);
            nmod_row_submul(vector.data(), basis_row, n, c, mod);
        }
    }

    const slong pivot = nmod_row_first_nonzero(vector.data(), n);
    if (pivot < 0) {
        *independent = false;
        return true;
    }

    const ulong inv = nmod_inv(vector[static_cast<std::size_t>(pivot)], mod);
    nmod_row_normalize(vector.data(), n, inv, mod);

    slong pos = 0;
    while (pos < rank_ &&
           pivots_[static_cast<std::size_t>(pos)] < pivot) {
        ++pos;
    }

    for (slong i = 0; i < rank_; ++i) {
        const ulong c = nmod_mat_entry(basis_.raw(), i, pivot);
        if (c != 0) {
            ulong* const basis_row = nmod_mat_entry_ptr(basis_.raw(), i, 0);
            nmod_row_submul(basis_row, vector.data(), n, c, mod);
        }
    }

    for (slong i = rank_; i > pos; --i) {
        pivots_[static_cast<std::size_t>(i)] =
                pivots_[static_cast<std::size_t>(i - 1)];
        for (slong j = 0; j < n; ++j) {
            nmod_mat_entry(basis_.raw(), i, j) =
                    nmod_mat_entry(basis_.raw(), i - 1, j);
        }
    }

    pivots_[static_cast<std::size_t>(pos)] = pivot;
    for (slong j = 0; j < n; ++j) {
        nmod_mat_entry(basis_.raw(), pos, j) =
                vector[static_cast<std::size_t>(j)];
    }
    ++rank_;
    *independent = true;

    return true;
}

bool ModRankContext::row_is_independent(
        bool* independent,
        const SparseRow& row) const noexcept {
    if (independent == nullptr || !is_set_ ||
        !row.fits_columns(ambient_dim_)) {
        return false;
    }

    const slong n = ambient_dim_;
    nmod_t mod;
    nmod_init(&mod, prime_);
    std::vector<ulong> vector(static_cast<std::size_t>(n));
    nmod_row_clear(vector.data(), n);

    for (slong i = 0; i < row.length(); ++i) {
        vector[static_cast<std::size_t>(
                row.cols_[static_cast<std::size_t>(i)])] =
                fmpz_get_nmod(row.vals_.data() + i, mod);
    }

    for (slong i = 0; i < rank_; ++i) {
        const slong pivot = pivots_[static_cast<std::size_t>(i)];
        const ulong c = vector[static_cast<std::size_t>(pivot)];
        if (c == 0) {
            continue;
        }
        for (slong j = 0; j < n; ++j) {
            const ulong product = nmod_mul(
                    c, nmod_mat_entry(basis_.raw(), i, j), mod);
            vector[static_cast<std::size_t>(j)] = nmod_sub(
                    vector[static_cast<std::size_t>(j)], product, mod);
        }
    }

    *independent = nmod_row_first_nonzero(vector.data(), n) >= 0;
    return true;
}

bool ModRankContext::add_fmpz_mat_row(bool* independent,
        flint::FmpzMatConstRef matrix,
        slong row) noexcept {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }

    SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }
    return add_row(independent, sparse_row);
}

HnfContext::HnfContext() noexcept
    : basis_(0), hnf_rows_(0), pending_dependent_rows_(0) {}

HnfContext::~HnfContext() noexcept {
    ambient_dim_ = 0;
    rank_ = 0;
    hnf_rows_processed_ = 0;
    has_mod_rank_ = false;
    hnf_valid_ = false;
    hnf_transform_valid_ = false;
    index_valid_ = false;
}

void HnfContext::swap(HnfContext& other) noexcept {
    std::swap(diagnostics_, other.diagnostics_);
    basis_.swap(other.basis_);
    hnf_rows_.swap(other.hnf_rows_);
    pending_dependent_rows_.swap(other.pending_dependent_rows_);
    mod_rank_.swap(other.mod_rank_);
    hnf_.swap(other.hnf_);
    hnf_transform_.swap(other.hnf_transform_);
    index_.swap(other.index_);
    std::swap(ambient_dim_, other.ambient_dim_);
    std::swap(rank_, other.rank_);
    std::swap(hnf_rows_processed_, other.hnf_rows_processed_);
    std::swap(has_mod_rank_, other.has_mod_rank_);
    std::swap(hnf_valid_, other.hnf_valid_);
    std::swap(hnf_transform_valid_, other.hnf_transform_valid_);
    std::swap(index_valid_, other.index_valid_);
}

bool HnfContext::set(const HnfContext& other) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.set");
    HnfContext tmp;
    tmp.diagnostics_ = diagnostics_;
    tmp.basis_ = SparseMat(other.ambient_dim_);
    tmp.hnf_rows_ = SparseMat(other.ambient_dim_);
    tmp.pending_dependent_rows_ = SparseMat(other.ambient_dim_);
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.set.sparse_copy");
        if (!tmp.basis_.set(other.basis_) ||
            !tmp.hnf_rows_.set(other.hnf_rows_) ||
            !tmp.pending_dependent_rows_.set(
                    other.pending_dependent_rows_)) {
            return false;
        }
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.set.dense_copy");
        tmp.hnf_ = flint::FmpzMat(flint::fmpz_mat_nrows(other.hnf_),
                                  flint::fmpz_mat_ncols(other.hnf_));
        tmp.hnf_transform_ = flint::FmpzMat(
                flint::fmpz_mat_nrows(other.hnf_transform_),
                flint::fmpz_mat_ncols(other.hnf_transform_));
        flint::fmpz_mat_set(tmp.hnf_, other.hnf_);
        flint::fmpz_mat_set(tmp.hnf_transform_, other.hnf_transform_);
    }
    flint::fmpz_set(tmp.index_, other.index_);
    tmp.ambient_dim_ = other.ambient_dim_;
    tmp.rank_ = other.rank_;
    tmp.hnf_rows_processed_ = other.hnf_rows_processed_;
    tmp.has_mod_rank_ = false;
    tmp.hnf_valid_ = other.hnf_valid_;
    tmp.hnf_transform_valid_ = other.hnf_transform_valid_;
    tmp.index_valid_ = other.index_valid_;

    if (other.has_mod_rank_) {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.set.mod_rank_copy");
        if (!tmp.mod_rank_.set(other.mod_rank_)) {
            return false;
        }
        tmp.has_mod_rank_ = true;
    }

    swap(tmp);
    return true;
}

bool HnfContext::reset(slong ambient_dim, ulong p) noexcept {
    if (ambient_dim < 0) {
        return false;
    }

    HnfContext tmp;
    tmp.diagnostics_ = diagnostics_;
    tmp.basis_ = SparseMat(ambient_dim);
    tmp.hnf_rows_ = SparseMat(ambient_dim);
    tmp.pending_dependent_rows_ = SparseMat(ambient_dim);
    tmp.hnf_ = flint::FmpzMat(0, ambient_dim);
    tmp.hnf_transform_ = flint::FmpzMat(0, 0);
    tmp.ambient_dim_ = ambient_dim;
    tmp.rank_ = 0;
    tmp.hnf_rows_processed_ = 0;
    tmp.has_mod_rank_ = false;
    tmp.hnf_valid_ = false;
    tmp.hnf_transform_valid_ = false;
    tmp.index_valid_ = false;
    flint::fmpz_zero(tmp.index_);

    if (p != 0) {
        if (!tmp.mod_rank_.set_prime_ui(ambient_dim, p)) {
            return false;
        }
        tmp.has_mod_rank_ = true;
    }

    swap(tmp);
    return true;
}

bool HnfContext::reset_precomputed_hnf(flint::FmpzMatConstRef hnf,
                                       ulong p) noexcept {
    if (hnf.raw() == nullptr || flint::fmpz_mat_nrows(hnf) < 0 ||
        flint::fmpz_mat_nrows(hnf) > flint::fmpz_mat_ncols(hnf)) {
        return false;
    }

    HnfContext tmp;
    tmp.diagnostics_ = diagnostics_;
    if (!tmp.reset(flint::fmpz_mat_ncols(hnf), p)) {
        return false;
    }
    tmp.rank_ = flint::fmpz_mat_nrows(hnf);
    tmp.hnf_ = flint::FmpzMat(tmp.rank_, tmp.ambient_dim_);
    if (!tmp.replace_hnf(hnf.raw())) {
        return false;
    }

    swap(tmp);
    return true;
}

bool HnfContext::add_row(bool* independent, const SparseRow& row) noexcept {
    return add_row(independent, nullptr, row);
}

bool HnfContext::add_row(bool* independent,
                         bool* index_refined,
                         const SparseRow& row) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.add_row");
    if (independent == nullptr || !row.fits_columns(ambient_dim_)) {
        return false;
    }
    if (index_refined != nullptr) {
        *index_refined = false;
    }

    if (rank_ >= ambient_dim_) {
        bool refined = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row.refine_dependent");
            if (!refine_dependent(refined, row)) {
                return false;
            }
            if (index_refined != nullptr) {
                *index_refined = refined;
            }
        }
        *independent = false;
        return true;
    }

    bool row_modular_dependent = false;
    if (has_mod_rank_ && mod_rank_.rank() == rank_) {
        bool mod_independent = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row.mod_rank");
            if (!mod_rank_.add_row(&mod_independent, row)) {
                return false;
            }
        }
        if (mod_independent) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.add_row.accept_independent");
                if (!accept_independent(row)) {
                    return false;
                }
            }
            *independent = true;
            return true;
        }
        row_modular_dependent = true;
    }

    if (row_modular_dependent) {
        bool in_current_lattice = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row.hnf_lattice_member");
            if (!row_reduces_to_zero_in_current_lattice(
                        in_current_lattice, row)) {
                return false;
            }
        }
        if (in_current_lattice) {
            *independent = false;
            return true;
        }
    }

    flint::FmpzMat candidate(rank_ + 1, ambient_dim_);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.add_row.rank_candidate_build");
        if (rank_ > 0) {
            flint::FmpzMatWindow basis_window(
                    candidate, 0, 0, rank_, ambient_dim_);
            basis_.get_fmpz_mat(basis_window.ref());
        }
        row.get_fmpz_mat_row(candidate.raw(), rank_);
    }
    slong candidate_rank = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.add_row.rank");
        candidate_rank = fmpz_mat_rank(candidate.raw());
    }

    if (candidate_rank > rank_) {
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row.accept_independent");
            if (!accept_independent(row)) {
                return false;
            }
        }
        *independent = true;
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.add_row.refine_dependent");
        bool refined = false;
        if (!refine_dependent(refined, row)) {
            return false;
        }
        if (index_refined != nullptr) {
            *index_refined = refined;
        }
    }
    *independent = false;
    return true;
}

bool HnfContext::add_fmpz_mat_row(bool* independent,
        flint::FmpzMatConstRef matrix,
        slong row) noexcept {
    return add_fmpz_mat_row(independent, nullptr, matrix, row);
}

bool HnfContext::add_fmpz_mat_row(bool* independent,
        bool* index_refined,
        flint::FmpzMatConstRef matrix,
        slong row) noexcept {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }

    SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }
    return add_row(independent, index_refined, sparse_row);
}

bool HnfContext::add_row_defer_dependent(bool* independent,
                                         const SparseRow& row) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.add_row_defer_dependent");
    if (independent == nullptr || !row.fits_columns(ambient_dim_)) {
        return false;
    }

    if (rank_ >= ambient_dim_) {
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row_defer_dependent.store_pending");
            if (!pending_dependent_rows_.append_row(row)) {
                return false;
            }
        }
        invalidate_index();
        *independent = false;
        return true;
    }

    bool row_modular_dependent = false;
    if (has_mod_rank_ && mod_rank_.rank() == rank_) {
        bool mod_independent = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row_defer_dependent.mod_rank");
            if (!mod_rank_.add_row(&mod_independent, row)) {
                return false;
            }
        }
        if (mod_independent) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.add_row_defer_dependent.accept_independent");
                if (!accept_independent(row)) {
                    return false;
                }
            }
            *independent = true;
            return true;
        }
        row_modular_dependent = true;
    }

    if (row_modular_dependent) {
        bool in_current_lattice = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row_defer_dependent.hnf_lattice_member");
            if (!row_reduces_to_zero_in_current_lattice(
                        in_current_lattice, row)) {
                return false;
            }
        }
        if (in_current_lattice) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.add_row_defer_dependent.store_pending");
                if (!pending_dependent_rows_.append_row(row)) {
                    return false;
                }
            }
            invalidate_index();
            *independent = false;
            return true;
        }
    }

    flint::FmpzMat candidate(rank_ + 1, ambient_dim_);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.add_row_defer_dependent.rank_candidate_build");
        if (rank_ > 0) {
            flint::FmpzMatWindow basis_window(
                    candidate, 0, 0, rank_, ambient_dim_);
            basis_.get_fmpz_mat(basis_window.ref());
        }
        row.get_fmpz_mat_row(flint::FmpzMatRef(candidate), rank_);
    }

    slong candidate_rank = 0;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.add_row_defer_dependent.rank");
        candidate_rank = fmpz_mat_rank(candidate.raw());
    }

    if (candidate_rank > rank_) {
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::class_group,
                    "fmpz_smat.hnf_context.add_row_defer_dependent.accept_independent");
            if (!accept_independent(row)) {
                return false;
            }
        }
        *independent = true;
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.add_row_defer_dependent.store_pending");
        if (!pending_dependent_rows_.append_row(row)) {
            return false;
        }
    }
    invalidate_index();
    *independent = false;
    return true;
}

bool HnfContext::add_fmpz_mat_row_defer_dependent(
        bool* independent,
        flint::FmpzMatConstRef matrix,
        slong row) noexcept {
    if (matrix.raw() == nullptr || row < 0 ||
        row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }

    SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }
    return add_row_defer_dependent(independent, sparse_row);
}

bool HnfContext::get_hnf(flint::FmpzMatRef matrix) noexcept {
    if (matrix.raw() == nullptr || !refresh(false) ||
        flint::fmpz_mat_nrows(matrix) != rank_ ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }
    flint::fmpz_mat_set(matrix, hnf_);
    return true;
}

bool HnfContext::get_hnf_rows(flint::FmpzMatRef matrix) noexcept {
    if (matrix.raw() == nullptr || !refresh(false) ||
        flint::fmpz_mat_nrows(matrix) != rank_ ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }
    hnf_rows_.get_fmpz_mat(matrix);
    return true;
}

bool HnfContext::get_hnf_transform(flint::FmpzMatRef matrix) noexcept {
    if (matrix.raw() == nullptr || !refresh(true) ||
        flint::fmpz_mat_nrows(matrix) != rank_ ||
        flint::fmpz_mat_ncols(matrix) != rank_) {
        return false;
    }
    flint::fmpz_mat_set(matrix, hnf_transform_);
    return true;
}

bool HnfContext::full_rank_index(flint::FmpzRef index) noexcept {
    if (index.raw() == nullptr || rank_ < ambient_dim_ ||
        !refresh(false) || !index_valid_) {
        return false;
    }
    flint::fmpz_set(index, index_);
    return true;
}

bool HnfContext::fmpz_mat_row_reduces_to_zero(
        bool& reduces_to_zero,
        flint::FmpzMatConstRef matrix,
        slong row) noexcept {
    reduces_to_zero = false;
    if (matrix.raw() == nullptr || row < 0 ||
        row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }

    SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }
    return row_reduces_to_zero_in_current_lattice(reduces_to_zero,
                                                  sparse_row);
}

bool HnfContext::refresh(bool require_transform) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.refresh");
    const slong pending_rows = pending_dependent_rows_.nrows();
    if (hnf_valid_ && hnf_rows_processed_ == rank_ && pending_rows == 0 &&
        (!require_transform || hnf_transform_valid_)) {
        return true;
    }

    if (rank_ == 0) {
        hnf_ = flint::FmpzMat(0, ambient_dim_);
        hnf_transform_ = flint::FmpzMat(0, 0);
        hnf_rows_ = SparseMat(ambient_dim_);
        pending_dependent_rows_ = SparseMat(ambient_dim_);
        hnf_valid_ = true;
        hnf_transform_valid_ = true;
        hnf_rows_processed_ = rank_;
        invalidate_index();
        return true;
    }

    const bool rebuild_transform_from_basis =
            require_transform && hnf_valid_ &&
            hnf_rows_processed_ == rank_ && pending_rows == 0 &&
            !hnf_transform_valid_;
    const slong old_rows =
            (hnf_valid_ && !rebuild_transform_from_basis)
                    ? flint::fmpz_mat_nrows(hnf_)
                    : 0;
    const slong first_new_row =
            (hnf_valid_ && !rebuild_transform_from_basis)
                    ? hnf_rows_processed_
                    : 0;
    const slong new_rows = rank_ - first_new_row;
    if (old_rows < 0 || first_new_row < 0 || first_new_row > rank_ ||
        new_rows < 0) {
        return false;
    }

    const slong dense_rows = old_rows + new_rows + pending_rows;
    const bool compute_transform = require_transform && pending_rows == 0;
    flint::FmpzMat dense(dense_rows, ambient_dim_);
    flint::FmpzMat hnf_full(dense_rows, ambient_dim_);
    flint::FmpzMat transform(compute_transform ? dense_rows : 0,
                             compute_transform ? dense_rows : 0);
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.refresh.dense_build");
        slong out_row = 0;
        if (old_rows > 0) {
            flint::FmpzMatWindow old_window(
                    dense, 0, 0, old_rows, ambient_dim_);
            flint::fmpz_mat_set(old_window.ref(),
                                flint::FmpzMatConstRef(hnf_));
            out_row = old_rows;
        }
        for (slong i = first_new_row; i < rank_; ++i) {
            basis_.row_ref(i).get_fmpz_mat_row(
                    flint::FmpzMatRef(dense), out_row);
            ++out_row;
        }
        for (slong i = 0; i < pending_rows; ++i) {
            pending_dependent_rows_.row_ref(i).get_fmpz_mat_row(
                    flint::FmpzMatRef(dense), out_row);
            ++out_row;
        }
        if (out_row != flint::fmpz_mat_nrows(dense)) {
            return false;
        }
    }
    if (compute_transform) {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.refresh.hnf_transform");
        ::fmpz_mat_hnf_transform(hnf_full.raw(), transform.raw(), dense.raw());
    } else {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.refresh.hnf");
        ::fmpz_mat_hnf(hnf_full.raw(), dense.raw());
    }

    const slong rows = trimmed_nonzero_rows(hnf_full.raw());
    if (rows != rank_) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.refresh.compress");
        flint::FmpzMatConstWindow hnf_window(
                hnf_full, 0, 0, rows, ambient_dim_);
        hnf_ = flint::FmpzMat(rows, ambient_dim_);
        flint::fmpz_mat_set(hnf_, hnf_window.const_ref());
        hnf_transform_ =
                (pending_rows > 0 || require_transform)
                        ? flint::FmpzMat(rows, rank_)
                        : flint::FmpzMat(0, 0);
        if (pending_rows > 0) {
            flint::fmpz_mat_one(flint::FmpzMatRef(hnf_transform_));
        } else if (compute_transform && rows > 0) {
            flint::FmpzMatConstWindow transform_window(
                    transform, 0, 0, rows, rank_);
            flint::fmpz_mat_set(hnf_transform_,
                                transform_window.const_ref());
        }

        hnf_rows_ = SparseMat(ambient_dim_);
        hnf_rows_.set_fmpz_mat(hnf_.raw());
        if (pending_rows > 0) {
            basis_ = SparseMat(ambient_dim_);
            basis_.set_fmpz_mat(hnf_.raw());
            hnf_transform_valid_ = true;
        } else {
            hnf_transform_valid_ = require_transform;
        }
        pending_dependent_rows_ = SparseMat(ambient_dim_);
    }
    hnf_valid_ = true;
    hnf_rows_processed_ = rank_;
    update_index_from_hnf();
    if (pending_rows > 0) {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.refresh.mod_rank");
        return rebuild_mod_rank();
    }
    return true;
}

bool HnfContext::accept_independent(const SparseRow& row) noexcept {
    if (!basis_.append_row(row)) {
        return false;
    }
    ++rank_;
    // Match the C fmpz_smat_hnf_ctx and reference's hnfadd direction: keep the
    // existing HNF as a checkpoint and fold in newly accepted rows on refresh.
    invalidate_index();
    hnf_transform_valid_ = false;
    return true;
}

bool HnfContext::row_reduces_to_zero_in_current_lattice(
        bool& reduces_to_zero,
        const SparseRow& row) noexcept {
    reduces_to_zero = false;
    return refresh(false) &&
           row_reduces_to_zero_mod_hnf(reduces_to_zero, row, hnf_.raw(),
                                       rank_, ambient_dim_);
}

bool HnfContext::refine_dependent(bool& index_refined,
                                  const SparseRow& row) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.refine_dependent");
    index_refined = false;

    bool remainder_zero = false;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.refine_dependent.hnf_remainder");
        if (!row_reduces_to_zero_in_current_lattice(remainder_zero, row)) {
            return false;
        }
    }
    if (remainder_zero) {
        return true;
    }

    flint::FmpzMat candidate(rank_ + 1, ambient_dim_);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.refine_dependent.candidate_build");
        if (rank_ > 0) {
            flint::FmpzMatWindow hnf_rows_window(
                    candidate, 0, 0, rank_, ambient_dim_);
            hnf_rows_.get_fmpz_mat(hnf_rows_window.ref());
        }
        row.get_fmpz_mat_row(candidate.raw(), rank_);
    }

    flint::FmpzMat hnf_full(rank_ + 1, ambient_dim_);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.refine_dependent.hnf");
        ::fmpz_mat_hnf(hnf_full.raw(), candidate.raw());
    }

    const slong rows = trimmed_nonzero_rows(hnf_full.raw());
    if (rows != rank_) {
        return false;
    }

    flint::FmpzMatConstWindow hnf_window(hnf_full, 0, 0, rows, ambient_dim_);
    flint::FmpzMat hnf_new(rows, ambient_dim_);
    flint::fmpz_mat_set(hnf_new, hnf_window.const_ref());

    index_refined = fmpz_mat_equal(hnf_new.raw(), hnf_.raw()) == 0;
    bool ok = true;
    if (index_refined) {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::class_group,
                "fmpz_smat.hnf_context.refine_dependent.replace_hnf");
        ok = replace_hnf(hnf_new.raw());
    }
    return ok;
}

bool HnfContext::replace_hnf(const fmpz_mat_t hnf) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                        "fmpz_smat.hnf_context.replace_hnf");
    if (fmpz_mat_nrows(hnf) != rank_ || fmpz_mat_ncols(hnf) != ambient_dim_) {
        return false;
    }
    basis_ = SparseMat(ambient_dim_);
    basis_.set_fmpz_mat(hnf);
    fmpz_mat_set(hnf_.raw(), hnf);
    hnf_transform_ = flint::FmpzMat(rank_, rank_);
    flint::fmpz_mat_one(flint::FmpzMatRef(hnf_transform_));
    hnf_rows_ = SparseMat(ambient_dim_);
    hnf_rows_.set_fmpz_mat(hnf_.raw());
    hnf_valid_ = true;
    hnf_transform_valid_ = true;
    hnf_rows_processed_ = rank_;
    update_index_from_hnf();
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::class_group,
                            "fmpz_smat.hnf_context.replace_hnf.mod_rank");
        return rebuild_mod_rank();
    }
}

void HnfContext::invalidate_index() noexcept {
    flint::fmpz_zero(index_);
    index_valid_ = false;
}

void HnfContext::update_index_from_hnf() noexcept {
    invalidate_index();
    if (rank_ == ambient_dim_ && rank_ > 0) {
        flint::Fmpz det;
        fmpz_mat_det(det.raw(), hnf_.raw());
        flint::fmpz_abs(index_, det);
        index_valid_ = !flint::fmpz_is_zero(index_);
    }
}

bool HnfContext::rebuild_mod_rank() noexcept {
    if (!has_mod_rank_) {
        return true;
    }
    const ulong p = mod_rank_.prime();
    ModRankContext rebuilt;
    if (!rebuilt.set_prime_ui(ambient_dim_, p)) {
        return false;
    }
    bool independent = false;
    for (slong i = 0; i < rank_; ++i) {
        if (!rebuilt.add_fmpz_mat_row(&independent, hnf_.raw(), i)) {
            return false;
        }
    }
    mod_rank_.swap(rebuilt);
    return true;
}

}  // namespace silex::fmpz_smat
