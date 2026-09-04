#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/fmpz_smat.hpp>

#include <vector>

namespace silex {

struct DiagnosticsContext;

namespace detail {

struct ClassRelationModuleAddResult {
    bool retained = false;
    bool duplicate = false;
    bool modular_independent = false;
};

class ClassRelationModuleContext {
public:
    ClassRelationModuleContext() noexcept = default;
    ~ClassRelationModuleContext() noexcept = default;

    ClassRelationModuleContext(const ClassRelationModuleContext&) = delete;
    ClassRelationModuleContext& operator=(const ClassRelationModuleContext&) =
            delete;
    ClassRelationModuleContext(ClassRelationModuleContext&&) noexcept =
            default;
    ClassRelationModuleContext& operator=(
            ClassRelationModuleContext&&) noexcept = default;

    bool reset(slong ambient_dim, ulong prime) noexcept;
    void clear() noexcept;
    bool set(const ClassRelationModuleContext& other) noexcept;
    void swap(ClassRelationModuleContext& other) noexcept;
    void set_diagnostics(const DiagnosticsContext* diagnostics) noexcept;

    bool is_set() const noexcept { return is_set_; }
    slong ambient_dim() const noexcept { return ambient_dim_; }
    ulong prime() const noexcept { return prime_; }
    slong rank() const noexcept { return mod_rank_.rank(); }
    slong basis_count() const noexcept { return bas_gens_.nrows(); }
    slong extra_count() const noexcept { return rel_gens_.nrows(); }
    slong relation_count() const noexcept {
        return basis_count() + extra_count();
    }

    const std::vector<slong>& basis_relation_indices() const noexcept {
        return bas_relation_indices_;
    }
    const std::vector<slong>& extra_relation_indices() const noexcept {
        return rel_relation_indices_;
    }

    bool add_fmpz_mat_row(ClassRelationModuleAddResult& result,
                          flint::FmpzMatConstRef matrix,
                          slong row,
                          slong relation_index,
                          bool always = true) noexcept;
    bool classify_fmpz_mat_row(ClassRelationModuleAddResult& result,
                               flint::FmpzMatConstRef matrix,
                               slong row,
                               bool always = true) const noexcept;

    bool contains_fmpz_mat_row(flint::FmpzMatConstRef matrix,
                               slong row) const noexcept;

    bool get_basis_rows(flint::FmpzMatRef out) const noexcept;
    bool get_extra_rows(flint::FmpzMatRef out) const noexcept;
    bool get_all_rows(flint::FmpzMatRef out) const noexcept;

    bool checkpoint_context(fmpz_smat::HnfContext& out) noexcept;
    bool take_checkpoint_context(fmpz_smat::HnfContext& out) noexcept;
    bool checkpoint_hnf_rows(flint::FmpzMatRef out) noexcept;
    bool checkpoint_hnf_transform(flint::FmpzMatRef out) noexcept;
    bool checkpoint_index(flint::FmpzRef out) noexcept;

private:
    bool add_sparse_row(ClassRelationModuleAddResult& result,
                        const fmpz_smat::SparseRow& row,
                        slong relation_index,
                        bool always) noexcept;
    bool classify_sparse_row(ClassRelationModuleAddResult& result,
                             const fmpz_smat::SparseRow& row,
                             bool always) const noexcept;
    bool contains_sparse_row(const fmpz_smat::SparseRow& row) const noexcept;
    bool rebuild_checkpoint() noexcept;

    fmpz_smat::ModRankContext mod_rank_;
    fmpz_smat::SparseMat bas_gens_{0};
    fmpz_smat::SparseMat rel_gens_{0};
    std::vector<slong> bas_relation_indices_;
    std::vector<slong> rel_relation_indices_;
    fmpz_smat::HnfContext checkpoint_;
    const DiagnosticsContext* diagnostics_ = nullptr;
    slong ambient_dim_ = 0;
    // Retained rows already represented by the canonical checkpoint.
    slong checkpoint_basis_count_ = 0;
    slong checkpoint_extra_count_ = 0;
    ulong prime_ = 0;
    bool is_set_ = false;
    bool checkpoint_valid_ = false;
};

}  // namespace detail
}  // namespace silex
