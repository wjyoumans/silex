#include <silex/detail/class_relation_module_context.hpp>

#include <silex/flint/fmpz.hpp>

#include <algorithm>
#include <cstddef>
#include <utility>

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
    result = {};
    if (!is_set_ || !fmpz_mat_is_row(matrix, row) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_ ||
        relation_index < 0) {
        return false;
    }

    fmpz_smat::SparseRow sparse_row;
    if (!sparse_row.set_fmpz_mat_row(matrix, row)) {
        return false;
    }

    return add_sparse_row(result, sparse_row, relation_index, always);
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
        ClassRelationModuleAddResult& result,
        const fmpz_smat::SparseRow& sparse_row,
        slong relation_index,
        bool always) noexcept {
    result = {};
    if (!is_set_ || relation_index < 0 ||
        !sparse_row.fits_columns(ambient_dim_)) {
        return false;
    }

    if (contains_sparse_row(sparse_row)) {
        result.duplicate = true;
        return true;
    }

    bool independent = false;
    if (!mod_rank_.add_row(&independent, sparse_row)) {
        return false;
    }
    result.modular_independent = independent;

    if (independent) {
        if (!bas_gens_.append_row(sparse_row)) {
            return false;
        }
        bas_relation_indices_.push_back(relation_index);
        result.retained = true;
        checkpoint_valid_ = false;
        return true;
    }

    if (always) {
        if (!rel_gens_.append_row(sparse_row)) {
            return false;
        }
        rel_relation_indices_.push_back(relation_index);
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
    if (old_rows == ambient_dim_ && old_rows > 0) {
        // Adding rows enlarges a full-rank lattice, so its new index divides
        // the retained checkpoint index.  This is FLINT's modular-HNF modulus
        // contract and the retained-index direction used by reference check_index
        // and reference hnfadd_i.
        flint::Fmpz checkpoint_index;
        if (!checkpoint_.full_rank_index(
                    flint::FmpzRef(checkpoint_index)) ||
            flint::fmpz_sgn(flint::FmpzConstRef(checkpoint_index)) <= 0) {
            return false;
        }
        ::fmpz_mat_hnf_modular(hnf.raw(), dense.raw(),
                               checkpoint_index.raw());
    } else {
        ::fmpz_mat_hnf(hnf.raw(), dense.raw());
    }
    const slong hnf_rows = trimmed_nonzero_rows(
            flint::FmpzMatConstRef(hnf));
    if (hnf_rows < 0) {
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
