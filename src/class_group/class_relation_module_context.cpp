#include <silex/detail/class_relation_module_context.hpp>

#include <silex/flint/fmpz.hpp>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace silex::detail {
namespace {

bool sparse_rows_equal(const fmpz_smat::SparseRow& left,
                       const fmpz_smat::SparseRow& right) noexcept {
    if (left.length() != right.length()) {
        return false;
    }
    if (!std::equal(left.columns().begin(), left.columns().end(),
                    right.columns().begin())) {
        return false;
    }

    flint::Fmpz left_value;
    flint::Fmpz right_value;
    for (slong i = 0; i < left.length(); ++i) {
        left.get_value(flint::FmpzRef(left_value), i);
        right.get_value(flint::FmpzRef(right_value), i);
        if (!flint::fmpz_equal(flint::FmpzConstRef(left_value),
                               flint::FmpzConstRef(right_value))) {
            return false;
        }
    }
    return true;
}

bool fmpz_mat_is_row(flint::FmpzMatConstRef matrix, slong row) noexcept {
    return matrix.raw() != nullptr && row >= 0 &&
           row < flint::fmpz_mat_nrows(matrix);
}

bool fmpz_mat_has_shape(flint::FmpzMatRef matrix,
                        slong rows,
                        slong cols) noexcept {
    return matrix.raw() != nullptr && flint::fmpz_mat_nrows(matrix) == rows &&
           flint::fmpz_mat_ncols(matrix) == cols;
}

slong trimmed_nonzero_rows(flint::FmpzMatConstRef matrix) noexcept {
    slong rows = flint::fmpz_mat_nrows(matrix);
    while (rows > 0 &&
           ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    for (slong row = 0; row < rows; ++row) {
        if (::fmpz_mat_is_zero_row(matrix.raw(), row) != 0) {
            return -1;
        }
    }
    return rows;
}

slong fmpz_mat_row_first_nonzero_from(const fmpz_mat_t matrix,
                                      slong row,
                                      slong first_column) noexcept {
    const slong cols = ::fmpz_mat_ncols(matrix);
    for (slong j = first_column; j < cols; ++j) {
        if (::fmpz_is_zero(::fmpz_mat_entry(matrix, row, j)) == 0) {
            return j;
        }
    }
    return -1;
}

// Reduces every entry above a pivot into [0, pivot), the normalization of
// FLINT's canonical row HNF (fmpz_mat_hnf).  Reducing row i2 by pivot row i
// changes row i2 only in columns from pivots[i] on, so one pass with i
// increasing leaves every earlier pivot column reduced.
void reduce_above_hnf_pivots(fmpz_mat_t hnf,
                             const std::vector<slong>& pivots) noexcept {
    const slong cols = ::fmpz_mat_ncols(hnf);
    const slong rank = static_cast<slong>(pivots.size());
    flint::Fmpz quotient;
    for (slong i = 1; i < rank; ++i) {
        const slong pivot_col = pivots[static_cast<std::size_t>(i)];
        const fmpz* pivot = ::fmpz_mat_entry(hnf, i, pivot_col);
        for (slong above = 0; above < i; ++above) {
            const fmpz* entry = ::fmpz_mat_entry(hnf, above, pivot_col);
            if (::fmpz_sgn(entry) >= 0 && ::fmpz_cmp(entry, pivot) < 0) {
                continue;
            }
            ::fmpz_fdiv_q(quotient.raw(), entry, pivot);
            for (slong j = pivot_col; j < cols; ++j) {
                const fmpz* pivot_row_entry = ::fmpz_mat_entry(hnf, i, j);
                if (::fmpz_is_zero(pivot_row_entry) == 0) {
                    ::fmpz_submul(::fmpz_mat_entry(hnf, above, j),
                                  quotient.raw(), pivot_row_entry);
                }
            }
        }
    }
}

// Folds rows [hnf_rows, nrows) of `matrix` into the canonical row HNF held in
// its first hnf_rows rows.  On success the first `rank` rows are the canonical
// row HNF of all rows (FLINT's fmpz_mat_hnf normalization: positive pivots in
// strictly increasing columns, entries above a pivot in [0, pivot)) and the
// remaining rows are zero.  This is the row-addition step of Pernet and
// Stein's HNF algorithm as FLINT implements it (FLINT 3.6.0
// src/fmpz_mat/hnf_pernet_stein.c, add_rows): each new row is eliminated
// against the pivot rows with unimodular extended-gcd row operations (a
// plain row subtraction when the pivot divides the entry, a local shortcut
// where FLINT always calls fmpz_xgcd), a surviving row becomes a new pivot row, and entries above the pivots are
// reduced again.  C. Pernet and W. Stein, "Fast computation of Hermite normal
// forms of random integer matrices", J. Number Theory 130 (2010), 1675-1683.
// The canonical row HNF of a lattice is unique, so the result equals
// fmpz_mat_hnf of the stacked rows.
//
// The code below is derived from that FLINT routine:
//
//   Copyright (C) 2014, 2015 Alex J. Best
//   This file is part of FLINT, licensed under the GNU Lesser General Public
//   License, version 3 or (at your option) any later version
//   (LGPL-3.0-or-later); see <https://www.gnu.org/licenses/>.
//
// The adaptation is distributed with Silex under GPL-3.0-or-later.
// Intentional deviations from FLINT: the divisible-entry shortcut above, a
// pivot-index walk with bubble insertion of new pivot rows, and skipped
// reductions when entries are already reduced.
bool hnf_add_rows(slong& rank, fmpz_mat_t matrix, slong hnf_rows) noexcept {
    rank = 0;
    const slong rows = ::fmpz_mat_nrows(matrix);
    const slong cols = ::fmpz_mat_ncols(matrix);
    if (hnf_rows < 0 || hnf_rows > rows || hnf_rows > cols) {
        return false;
    }

    std::vector<slong> pivots;
    pivots.reserve(static_cast<std::size_t>(std::min(rows, cols)));
    for (slong i = 0; i < hnf_rows; ++i) {
        const slong first =
                pivots.empty() ? 0 : pivots.back() + 1;
        const slong pivot_col =
                fmpz_mat_row_first_nonzero_from(matrix, i, 0);
        if (pivot_col < first ||
            ::fmpz_sgn(::fmpz_mat_entry(matrix, i, pivot_col)) <= 0) {
            return false;
        }
        pivots.push_back(pivot_col);
    }
    rank = hnf_rows;

    flint::Fmpz gcd;
    flint::Fmpz pivot_coeff;
    flint::Fmpz row_coeff;
    flint::Fmpz pivot_cofactor;
    flint::Fmpz row_cofactor;
    flint::Fmpz combined;
    for (slong source = hnf_rows; source < rows; ++source) {
        // Rows [rank, source) are zero, so moving the new row up to `rank`
        // keeps the pivot rows contiguous.
        if (source != rank) {
            ::fmpz_mat_swap_rows(matrix, nullptr, rank, source);
        }
        const slong added = rank;
        bool pivot_rows_changed = false;
        slong pivot_index = 0;
        slong lead = fmpz_mat_row_first_nonzero_from(matrix, added, 0);
        while (lead >= 0) {
            while (pivot_index < rank &&
                   pivots[static_cast<std::size_t>(pivot_index)] < lead) {
                ++pivot_index;
            }
            if (pivot_index == rank ||
                pivots[static_cast<std::size_t>(pivot_index)] != lead) {
                break;
            }

            fmpz* pivot = ::fmpz_mat_entry(matrix, pivot_index, lead);
            fmpz* entry = ::fmpz_mat_entry(matrix, added, lead);
            if (::fmpz_divisible(entry, pivot) != 0) {
                ::fmpz_divexact(row_cofactor.raw(), entry, pivot);
                for (slong j = lead; j < cols; ++j) {
                    const fmpz* pivot_row_entry =
                            ::fmpz_mat_entry(matrix, pivot_index, j);
                    if (::fmpz_is_zero(pivot_row_entry) == 0) {
                        ::fmpz_submul(::fmpz_mat_entry(matrix, added, j),
                                      row_cofactor.raw(), pivot_row_entry);
                    }
                }
            } else {
                // [pivot_row; added] <- [s t; -b/g a/g] [pivot_row; added]
                // with s a + t b = g = gcd(a, b): determinant one, new pivot g.
                ::fmpz_xgcd(gcd.raw(), pivot_coeff.raw(), row_coeff.raw(),
                            pivot, entry);
                ::fmpz_divexact(pivot_cofactor.raw(), pivot, gcd.raw());
                ::fmpz_divexact(row_cofactor.raw(), entry, gcd.raw());
                for (slong j = lead; j < cols; ++j) {
                    fmpz* pivot_row_entry =
                            ::fmpz_mat_entry(matrix, pivot_index, j);
                    fmpz* added_entry = ::fmpz_mat_entry(matrix, added, j);
                    ::fmpz_mul(combined.raw(), pivot_coeff.raw(),
                               pivot_row_entry);
                    ::fmpz_addmul(combined.raw(), row_coeff.raw(),
                                  added_entry);
                    ::fmpz_mul(added_entry, added_entry,
                               pivot_cofactor.raw());
                    ::fmpz_submul(added_entry, row_cofactor.raw(),
                                  pivot_row_entry);
                    ::fmpz_swap(pivot_row_entry, combined.raw());
                }
                pivot_rows_changed = true;
            }
            ++pivot_index;
            lead = fmpz_mat_row_first_nonzero_from(matrix, added, lead + 1);
        }

        if (lead < 0) {
            // The new row lies in the span of the pivot rows; it is now zero.
            if (pivot_rows_changed) {
                reduce_above_hnf_pivots(matrix, pivots);
            }
            continue;
        }

        if (::fmpz_sgn(::fmpz_mat_entry(matrix, added, lead)) < 0) {
            for (slong j = lead; j < cols; ++j) {
                fmpz* added_entry = ::fmpz_mat_entry(matrix, added, j);
                ::fmpz_neg(added_entry, added_entry);
            }
        }
        for (slong i = added; i > pivot_index; --i) {
            ::fmpz_mat_swap_rows(matrix, nullptr, i, i - 1);
        }
        pivots.insert(pivots.begin() + pivot_index, lead);
        ++rank;
        reduce_above_hnf_pivots(matrix, pivots);
    }
    return true;
}

}  // namespace

bool ClassRelationModuleContext::reset(slong ambient_dim,
                                       ulong prime) noexcept {
    if (ambient_dim < 0 || prime < 2) {
        return false;
    }

    fmpz_smat::ModRankContext mod_rank;
    if (!mod_rank.set_prime_ui(ambient_dim, prime)) {
        return false;
    }

    fmpz_smat::HnfContext checkpoint;
    checkpoint.set_diagnostics(diagnostics_);
    if (!checkpoint.reset(ambient_dim, prime)) {
        return false;
    }

    mod_rank_ = std::move(mod_rank);
    bas_gens_ = fmpz_smat::SparseMat(ambient_dim);
    rel_gens_ = fmpz_smat::SparseMat(ambient_dim);
    bas_relation_indices_.clear();
    rel_relation_indices_.clear();
    checkpoint_ = std::move(checkpoint);
    ambient_dim_ = ambient_dim;
    checkpoint_basis_count_ = 0;
    checkpoint_extra_count_ = 0;
    prime_ = prime;
    is_set_ = true;
    checkpoint_valid_ = true;
    return true;
}

void ClassRelationModuleContext::clear() noexcept {
    mod_rank_ = fmpz_smat::ModRankContext();
    bas_gens_ = fmpz_smat::SparseMat(0);
    rel_gens_ = fmpz_smat::SparseMat(0);
    bas_relation_indices_.clear();
    rel_relation_indices_.clear();
    checkpoint_ = fmpz_smat::HnfContext();
    ambient_dim_ = 0;
    checkpoint_basis_count_ = 0;
    checkpoint_extra_count_ = 0;
    prime_ = 0;
    is_set_ = false;
    checkpoint_valid_ = false;
    checkpoint_.set_diagnostics(diagnostics_);
}

bool ClassRelationModuleContext::set(
        const ClassRelationModuleContext& other) noexcept {
    ClassRelationModuleContext tmp;
    tmp.set_diagnostics(diagnostics_);
    if (!other.is_set_) {
        swap(tmp);
        return true;
    }
    if (other.bas_relation_indices_.size() !=
                    static_cast<std::size_t>(other.bas_gens_.nrows()) ||
            other.rel_relation_indices_.size() !=
                    static_cast<std::size_t>(other.rel_gens_.nrows()) ||
            other.checkpoint_basis_count_ < 0 ||
            other.checkpoint_basis_count_ > other.basis_count() ||
            other.checkpoint_extra_count_ < 0 ||
            other.checkpoint_extra_count_ > other.extra_count() ||
            other.checkpoint_.ambient_dim() != other.ambient_dim_ ||
            (other.checkpoint_valid_ &&
             (other.checkpoint_basis_count_ != other.basis_count() ||
              other.checkpoint_extra_count_ != other.extra_count())) ||
            !tmp.reset(other.ambient_dim_, other.prime_)) {
        return false;
    }

    if (!tmp.mod_rank_.set(other.mod_rank_) ||
        !tmp.bas_gens_.set(other.bas_gens_) ||
        !tmp.rel_gens_.set(other.rel_gens_)) {
        return false;
    }
    tmp.bas_relation_indices_ = other.bas_relation_indices_;
    tmp.rel_relation_indices_ = other.rel_relation_indices_;
    if (!tmp.checkpoint_.set(other.checkpoint_)) {
        return false;
    }
    tmp.checkpoint_basis_count_ = other.checkpoint_basis_count_;
    tmp.checkpoint_extra_count_ = other.checkpoint_extra_count_;
    tmp.checkpoint_valid_ = other.checkpoint_valid_;

    swap(tmp);
    return true;
}

void ClassRelationModuleContext::swap(
        ClassRelationModuleContext& other) noexcept {
    mod_rank_.swap(other.mod_rank_);
    bas_gens_.swap(other.bas_gens_);
    rel_gens_.swap(other.rel_gens_);
    bas_relation_indices_.swap(other.bas_relation_indices_);
    rel_relation_indices_.swap(other.rel_relation_indices_);
    checkpoint_.swap(other.checkpoint_);
    std::swap(diagnostics_, other.diagnostics_);
    std::swap(ambient_dim_, other.ambient_dim_);
    std::swap(checkpoint_basis_count_, other.checkpoint_basis_count_);
    std::swap(checkpoint_extra_count_, other.checkpoint_extra_count_);
    std::swap(prime_, other.prime_);
    std::swap(is_set_, other.is_set_);
    std::swap(checkpoint_valid_, other.checkpoint_valid_);
}

void ClassRelationModuleContext::set_diagnostics(
        const DiagnosticsContext* diagnostics) noexcept {
    diagnostics_ = diagnostics;
    checkpoint_.set_diagnostics(diagnostics);
}

bool ClassRelationModuleContext::add_fmpz_mat_row(
        ClassRelationModuleAddResult& result,
        flint::FmpzMatConstRef matrix,
        slong row,
        slong relation_index,
        bool always) noexcept {
    ClassRelationModuleAddUndo undo;
    return add_fmpz_mat_row(undo, result, matrix, row, relation_index,
                            always);
}

bool ClassRelationModuleContext::add_fmpz_mat_row(
        ClassRelationModuleAddUndo& undo,
        ClassRelationModuleAddResult& result,
        flint::FmpzMatConstRef matrix,
        slong row,
        slong relation_index,
        bool always) noexcept {
    result = {};
    if (undo.changed || !is_set_ || !fmpz_mat_is_row(matrix, row) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_ ||
        relation_index < 0) {
        return false;
    }

    fmpz_smat::SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }

    return add_sparse_row(undo, result, sparse_row, relation_index, always);
}

bool ClassRelationModuleContext::rollback_add(
        ClassRelationModuleAddUndo& undo) noexcept {
    if (!undo.changed) {
        return true;
    }
    if (!is_set_) {
        return false;
    }
    // The checkpoint must not yet cover the row being removed.
    if (undo.basis_row) {
        if (basis_count() <= 0 || bas_relation_indices_.empty() ||
            checkpoint_basis_count_ >= basis_count() ||
            !mod_rank_.undo_add_row(undo.mod_rank_position,
                                    undo.mod_rank_eliminated,
                                    undo.mod_rank_displaced_pivot)) {
            return false;
        }
        bas_gens_.pop_last_row();
        bas_relation_indices_.pop_back();
    } else {
        if (extra_count() <= 0 || rel_relation_indices_.empty() ||
            checkpoint_extra_count_ >= extra_count()) {
            return false;
        }
        rel_gens_.pop_last_row();
        rel_relation_indices_.pop_back();
    }
    // A checkpoint taken or rebuilt since the add no longer matches the
    // record, so the restored flag must also agree with the counts.
    checkpoint_valid_ = undo.checkpoint_valid &&
                        checkpoint_basis_count_ == basis_count() &&
                        checkpoint_extra_count_ == extra_count();
    undo.changed = false;
    undo.basis_row = false;
    undo.mod_rank_position = -1;
    undo.mod_rank_eliminated.clear();
    return true;
}

bool ClassRelationModuleContext::same_state_as(
        const ClassRelationModuleContext& other) const noexcept {
    if (is_set_ != other.is_set_) {
        return false;
    }
    if (!is_set_) {
        return true;
    }
    if (ambient_dim_ != other.ambient_dim_ || prime_ != other.prime_ ||
        checkpoint_basis_count_ != other.checkpoint_basis_count_ ||
        checkpoint_extra_count_ != other.checkpoint_extra_count_ ||
        checkpoint_valid_ != other.checkpoint_valid_ ||
        bas_relation_indices_ != other.bas_relation_indices_ ||
        rel_relation_indices_ != other.rel_relation_indices_ ||
        bas_gens_.nrows() != other.bas_gens_.nrows() ||
        bas_gens_.nnz() != other.bas_gens_.nnz() ||
        rel_gens_.nrows() != other.rel_gens_.nrows() ||
        rel_gens_.nnz() != other.rel_gens_.nnz() ||
        !mod_rank_.same_state_as(other.mod_rank_)) {
        return false;
    }
    for (slong i = 0; i < bas_gens_.nrows(); ++i) {
        if (!sparse_rows_equal(bas_gens_.row_ref(i),
                               other.bas_gens_.row_ref(i))) {
            return false;
        }
    }
    for (slong i = 0; i < rel_gens_.nrows(); ++i) {
        if (!sparse_rows_equal(rel_gens_.row_ref(i),
                               other.rel_gens_.row_ref(i))) {
            return false;
        }
    }
    return true;
}

bool ClassRelationModuleContext::classify_fmpz_mat_row(
        ClassRelationModuleAddResult& result,
        flint::FmpzMatConstRef matrix,
        slong row,
        bool always) const noexcept {
    result = {};
    if (!is_set_ || !fmpz_mat_is_row(matrix, row) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }

    fmpz_smat::SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }

    return classify_sparse_row(result, sparse_row, always);
}

bool ClassRelationModuleContext::add_sparse_row(
        ClassRelationModuleAddUndo& undo,
        ClassRelationModuleAddResult& result,
        const fmpz_smat::SparseRow& sparse_row,
        slong relation_index,
        bool always) noexcept {
    result = {};
    if (undo.changed || !is_set_ || relation_index < 0 ||
        !sparse_row.fits_columns(ambient_dim_)) {
        return false;
    }

    if (contains_sparse_row(sparse_row)) {
        result.duplicate = true;
        return true;
    }

    // A dependent row leaves mod_rank_ unchanged; an independent one is
    // recorded so a failed append below can be undone.
    bool independent = false;
    if (!mod_rank_.add_row_recording(&independent, sparse_row,
                                     &undo.mod_rank_position,
                                     &undo.mod_rank_eliminated,
                                     &undo.mod_rank_displaced_pivot)) {
        return false;
    }
    result.modular_independent = independent;

    if (independent) {
        if (!bas_gens_.append_row(sparse_row)) {
            // mod_rank_ was just changed by the recorded add; this inverse
            // cannot fail on its own record.
            static_cast<void>(mod_rank_.undo_add_row(
                    undo.mod_rank_position, undo.mod_rank_eliminated,
                    undo.mod_rank_displaced_pivot));
            result = {};
            return false;
        }
        bas_relation_indices_.push_back(relation_index);
        undo.changed = true;
        undo.basis_row = true;
        undo.checkpoint_valid = checkpoint_valid_;
        result.retained = true;
        checkpoint_valid_ = false;
        return true;
    }

    if (always) {
        if (!rel_gens_.append_row(sparse_row)) {
            result = {};
            return false;
        }
        rel_relation_indices_.push_back(relation_index);
        undo.changed = true;
        undo.basis_row = false;
        undo.checkpoint_valid = checkpoint_valid_;
        result.retained = true;
        checkpoint_valid_ = false;
    }

    return true;
}

bool ClassRelationModuleContext::classify_sparse_row(
        ClassRelationModuleAddResult& result,
        const fmpz_smat::SparseRow& sparse_row,
        bool always) const noexcept {
    result = {};
    if (!is_set_ || !sparse_row.fits_columns(ambient_dim_)) {
        return false;
    }

    if (contains_sparse_row(sparse_row)) {
        result.duplicate = true;
        return true;
    }

    bool independent = false;
    if (!mod_rank_.row_is_independent(&independent, sparse_row)) {
        return false;
    }
    result.modular_independent = independent;
    result.retained = independent || always;
    return true;
}

bool ClassRelationModuleContext::contains_fmpz_mat_row(
        flint::FmpzMatConstRef matrix,
        slong row) const noexcept {
    if (!is_set_ || !fmpz_mat_is_row(matrix, row) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }

    fmpz_smat::SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }
    return contains_sparse_row(sparse_row);
}

bool ClassRelationModuleContext::get_basis_rows(
        flint::FmpzMatRef out) const noexcept {
    if (!is_set_ || !fmpz_mat_has_shape(out, basis_count(), ambient_dim_)) {
        return false;
    }
    bas_gens_.get_fmpz_mat(out);
    return true;
}

bool ClassRelationModuleContext::get_extra_rows(
        flint::FmpzMatRef out) const noexcept {
    if (!is_set_ || !fmpz_mat_has_shape(out, extra_count(), ambient_dim_)) {
        return false;
    }
    rel_gens_.get_fmpz_mat(out);
    return true;
}

bool ClassRelationModuleContext::get_all_rows(
        flint::FmpzMatRef out) const noexcept {
    if (!is_set_ || !fmpz_mat_has_shape(out, relation_count(), ambient_dim_)) {
        return false;
    }
    flint::fmpz_mat_zero(out);
    for (slong i = 0; i < basis_count(); ++i) {
        bas_gens_.row_ref(i).get_fmpz_mat_row(out, i);
    }
    for (slong i = 0; i < extra_count(); ++i) {
        rel_gens_.row_ref(i).get_fmpz_mat_row(out, basis_count() + i);
    }
    return true;
}

bool ClassRelationModuleContext::checkpoint_context(
        fmpz_smat::HnfContext& out) noexcept {
    if (!is_set_ || !rebuild_checkpoint()) {
        return false;
    }
    return out.set(checkpoint_);
}

bool ClassRelationModuleContext::take_checkpoint_context(
        fmpz_smat::HnfContext& out) noexcept {
    if (!is_set_ || !rebuild_checkpoint()) {
        return false;
    }

    fmpz_smat::HnfContext replacement;
    replacement.set_diagnostics(diagnostics_);
    if (!replacement.reset(ambient_dim_, prime_)) {
        return false;
    }
    out.swap(checkpoint_);
    out.set_diagnostics(diagnostics_);
    checkpoint_ = std::move(replacement);
    checkpoint_basis_count_ = 0;
    checkpoint_extra_count_ = 0;
    checkpoint_valid_ = false;
    return true;
}

bool ClassRelationModuleContext::checkpoint_hnf_rows(
        flint::FmpzMatRef out) noexcept {
    if (!is_set_ || !rebuild_checkpoint() ||
        !fmpz_mat_has_shape(out, checkpoint_.rank(), ambient_dim_)) {
        return false;
    }
    return checkpoint_.get_hnf_rows(out);
}

bool ClassRelationModuleContext::checkpoint_hnf_transform(
        flint::FmpzMatRef out) noexcept {
    if (!is_set_ || !rebuild_checkpoint() ||
        !fmpz_mat_has_shape(out, checkpoint_.rank(), checkpoint_.rank())) {
        return false;
    }
    return checkpoint_.get_hnf_transform(out);
}

bool ClassRelationModuleContext::checkpoint_index(flint::FmpzRef out) noexcept {
    if (out.raw() == nullptr || !is_set_) {
        return false;
    }
    if (rank() < ambient_dim_) {
        flint::fmpz_zero(out);
        return true;
    }
    if (!rebuild_checkpoint()) {
        return false;
    }
    return checkpoint_.full_rank_index(out);
}

bool ClassRelationModuleContext::contains_sparse_row(
        const fmpz_smat::SparseRow& row) const noexcept {
    for (slong i = 0; i < bas_gens_.nrows(); ++i) {
        if (sparse_rows_equal(bas_gens_.row_ref(i), row)) {
            return true;
        }
    }
    for (slong i = 0; i < rel_gens_.nrows(); ++i) {
        if (sparse_rows_equal(rel_gens_.row_ref(i), row)) {
            return true;
        }
    }
    return false;
}

bool ClassRelationModuleContext::rebuild_checkpoint() noexcept {
    if (!is_set_) {
        return false;
    }
    if (checkpoint_valid_) {
        return true;
    }

    if (checkpoint_basis_count_ < 0 ||
        checkpoint_basis_count_ > basis_count() ||
        checkpoint_extra_count_ < 0 ||
        checkpoint_extra_count_ > extra_count()) {
        return false;
    }

    const slong old_rows = checkpoint_.rank();
    const slong new_basis_rows = basis_count() - checkpoint_basis_count_;
    const slong new_extra_rows = extra_count() - checkpoint_extra_count_;
    const bool full_rank_checkpoint = old_rows == ambient_dim_ && old_rows > 0;
    flint::Fmpz checkpoint_index;
    const bool checkpoint_index_valid =
            full_rank_checkpoint &&
            checkpoint_.full_rank_index(flint::FmpzRef(checkpoint_index)) &&
            flint::fmpz_sgn(flint::FmpzConstRef(checkpoint_index)) > 0;
    if (full_rank_checkpoint && new_basis_rows == 0) {
        // Invariant: checkpoint_ holds the canonical row HNF H of the lattice
        // spanned by the first checkpoint_basis_count_ basis rows and
        // checkpoint_extra_count_ extra rows, as reset_precomputed_hnf left it
        // (identity transform, index |det H|).  This holds only because
        // nothing calls add_row on checkpoint_ in place; an in-place update
        // would leave a state this skip must not reuse.  When H has full rank
        // and every extra row added since lies in the lattice of H, the rows
        // of H and the new rows span that same lattice, and the row HNF of a
        // lattice is unique (the fmpz_mat_hnf contract), so the rebuild below
        // would recompute H and reset checkpoint_ to the state it already
        // holds.  Membership is the exact remainder test of a row modulo the
        // full-rank, upper-triangular H, the test HnfContext::refine_dependent
        // uses to leave its HNF unchanged.  A failed test, or a row outside
        // the lattice, leaves the decision to the full rebuild.
        bool lattice_unchanged = checkpoint_index_valid;
        for (slong i = checkpoint_extra_count_;
             lattice_unchanged && i < extra_count(); ++i) {
            bool reduces_to_zero = false;
            lattice_unchanged =
                    checkpoint_.row_reduces_to_zero_in_current_lattice(
                            reduces_to_zero, rel_gens_.row_ref(i)) &&
                    reduces_to_zero;
        }
        if (lattice_unchanged) {
            checkpoint_extra_count_ = extra_count();
            checkpoint_valid_ = true;
            return true;
        }
    }
    const slong dense_rows = old_rows + new_basis_rows + new_extra_rows;
    flint::FmpzMat dense(dense_rows, ambient_dim_);
    flint::FmpzMat hnf(dense_rows, ambient_dim_);
    slong out_row = 0;
    if (old_rows > 0) {
        flint::FmpzMatWindow old_window(
                dense, 0, 0, old_rows, ambient_dim_);
        if (!checkpoint_.get_hnf_rows(old_window.ref())) {
            return false;
        }
        out_row = old_rows;
    }
    for (slong i = checkpoint_basis_count_; i < basis_count(); ++i) {
        bas_gens_.row_ref(i).get_fmpz_mat_row(
                flint::FmpzMatRef(dense), out_row);
        ++out_row;
    }
    for (slong i = checkpoint_extra_count_; i < extra_count(); ++i) {
        rel_gens_.row_ref(i).get_fmpz_mat_row(
                flint::FmpzMatRef(dense), out_row);
        ++out_row;
    }
    if (out_row != dense_rows) {
        return false;
    }
    slong incremental_rank = -1;
    if (full_rank_checkpoint) {
        // Adding rows enlarges a full-rank lattice, so its new index divides
        // the retained checkpoint index.  This is FLINT's modular-HNF modulus
        // contract and the retained-index direction used by reference check_index
        // and reference hnfadd_i.
        if (!checkpoint_index_valid) {
            return false;
        }
        ::fmpz_mat_hnf_modular(hnf.raw(), dense.raw(),
                               checkpoint_index.raw());
    } else if (old_rows > 0) {
        // The retained checkpoint rows are already a canonical HNF, so only
        // the new rows need eliminating against them.
        if (!hnf_add_rows(incremental_rank, dense.raw(), old_rows)) {
            return false;
        }
        hnf.swap(dense);
    } else {
        ::fmpz_mat_hnf(hnf.raw(), dense.raw());
    }
    const slong hnf_rows = trimmed_nonzero_rows(
            flint::FmpzMatConstRef(hnf));
    if (hnf_rows < 0 ||
        (incremental_rank >= 0 && incremental_rank != hnf_rows)) {
        return false;
    }

    flint::FmpzMat hnf_basis(hnf_rows, ambient_dim_);
    if (hnf_rows > 0) {
        flint::FmpzMatConstWindow hnf_window(
                hnf, 0, 0, hnf_rows, ambient_dim_);
        flint::fmpz_mat_set(flint::FmpzMatRef(hnf_basis),
                            hnf_window.const_ref());
    }
    fmpz_smat::HnfContext checkpoint;
    checkpoint.set_diagnostics(diagnostics_);
    if (!checkpoint.reset_precomputed_hnf(
                flint::FmpzMatConstRef(hnf_basis), prime_)) {
        return false;
    }

    checkpoint_ = std::move(checkpoint);
    checkpoint_basis_count_ = basis_count();
    checkpoint_extra_count_ = extra_count();
    checkpoint_valid_ = true;
    return true;
}

}  // namespace silex::detail
