#pragma once

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/flint/nmod_mat.hpp>

#include <span>
#include <utility>
#include <vector>

namespace silex {

struct DiagnosticsContext;

namespace detail {
class ClassRelationModuleContext;
}

}  // namespace silex

namespace silex::fmpz_smat {

class SparseRow {
public:
    SparseRow() noexcept = default;
    ~SparseRow() noexcept = default;

    SparseRow(const SparseRow&) = delete;
    SparseRow& operator=(const SparseRow&) = delete;

    SparseRow(SparseRow&& other) noexcept { swap(other); }
    SparseRow& operator=(SparseRow&& other) noexcept {
        if (this != &other) {
            SparseRow tmp;
            tmp.swap(*this);
            swap(other);
        }
        return *this;
    }

    void swap(SparseRow& other) noexcept;
    bool set(const SparseRow& other) noexcept;

    slong length() const noexcept {
        return static_cast<slong>(cols_.size());
    }
    std::span<const slong> columns() const noexcept;
    flint::FmpzVecConstRef values() const noexcept;
    slong col(slong i) const;
    void get_value(flint::FmpzRef out, slong i) const;
    void get_entry(flint::FmpzRef out, slong col) const;
    bool set_entry(slong col, flint::FmpzConstRef value) noexcept;
    bool append_entry(slong col, flint::FmpzConstRef value) noexcept;
    bool set_fmpz_mat_row(flint::FmpzMatConstRef matrix, slong row) noexcept;
    void get_fmpz_mat_row(flint::FmpzMatRef matrix, slong row) const;
    flint::FmpzMat to_fmpz_mat_row(slong ncols) const;
    void dot_fmpz_mat_row(flint::FmpzRef out,
            flint::FmpzMatConstRef matrix,
            slong row) const;
    void get_nmod_mat_row(flint::NmodMatRef matrix, slong row) const;
    bool fits_columns(slong ncols) const noexcept { return fits(ncols); }

private:
    bool assure_alloc(slong alloc) noexcept;
    slong find(slong col, bool* found) const noexcept;
    bool fits(slong ncols) const noexcept;

    std::vector<slong> cols_;
    flint::FmpzVec vals_{0};

    friend class SparseMat;
    friend class ModRankContext;
};

inline void swap(SparseRow& left, SparseRow& right) noexcept { left.swap(right); }

class SparseMat {
public:
    explicit SparseMat(slong ncols);
    ~SparseMat() noexcept = default;

    SparseMat(const SparseMat&) = delete;
    SparseMat& operator=(const SparseMat&) = delete;

    SparseMat(SparseMat&& other) noexcept { swap(other); }
    SparseMat& operator=(SparseMat&& other) noexcept {
        if (this != &other) {
            SparseMat tmp(0);
            tmp.swap(*this);
            swap(other);
        }
        return *this;
    }

    void swap(SparseMat& other) noexcept;
    bool set(const SparseMat& other) noexcept;

    slong nrows() const noexcept {
        return static_cast<slong>(row_data_.size());
    }
    slong ncols() const noexcept { return cols_; }
    slong nnz() const noexcept { return nnz_; }

    bool append_row(const SparseRow& row) noexcept;
    bool append_fmpz_mat_row(flint::FmpzMatConstRef matrix, slong row) noexcept;
    bool set_row(slong row_index, const SparseRow& row) noexcept;
    const SparseRow& row_ref(slong row_index) const;
    void get_row(SparseRow& row, slong row_index) const;
    void get_entry(flint::FmpzRef out, slong row, slong col) const;
    void get_fmpz_mat(flint::FmpzMatRef matrix) const;
    flint::FmpzMat to_fmpz_mat() const;
    void set_fmpz_mat(flint::FmpzMatConstRef matrix) noexcept;
    void mul_fmpz_mat(flint::FmpzMatRef out, flint::FmpzMatConstRef right) const;
    void transpose(SparseMat& out) const noexcept;
    void get_nmod_mat(flint::NmodMatRef matrix) const;
    bool rank_mod_prime_ui(slong* rank, ulong p) const noexcept;

private:
    std::vector<SparseRow> row_data_;
    slong cols_ = 0;
    slong nnz_ = 0;
};

inline void swap(SparseMat& left, SparseMat& right) noexcept { left.swap(right); }

class ModRankContext {
public:
    ModRankContext() noexcept;
    ~ModRankContext() noexcept;

    ModRankContext(const ModRankContext&) = delete;
    ModRankContext& operator=(const ModRankContext&) = delete;

    ModRankContext(ModRankContext&& other) noexcept { swap(other); }
    ModRankContext& operator=(ModRankContext&& other) noexcept {
        if (this != &other) {
            ModRankContext tmp;
            tmp.swap(*this);
            swap(other);
        }
        return *this;
    }

    void swap(ModRankContext& other) noexcept;
    bool set(const ModRankContext& other) noexcept;
    bool set_prime_ui(slong ambient_dim, ulong p) noexcept;
    void reset() noexcept;
    slong rank() const noexcept { return rank_; }
    slong ambient_dim() const noexcept { return ambient_dim_; }
    ulong prime() const noexcept { return prime_; }
    bool is_set() const noexcept { return is_set_; }
    bool add_row(bool* independent, const SparseRow& row) noexcept;
    bool row_is_independent(bool* independent,
                            const SparseRow& row) const noexcept;
    bool add_fmpz_mat_row(bool* independent,
            flint::FmpzMatConstRef matrix,
            slong row) noexcept;

private:
    flint::NmodMat basis_{0, 0, 2};
    std::vector<slong> pivots_;
    slong ambient_dim_ = 0;
    slong rank_ = 0;
    ulong prime_ = 0;
    bool is_set_ = false;
};

inline void swap(ModRankContext& left, ModRankContext& right) noexcept {
    left.swap(right);
}

class HnfContext {
public:
    HnfContext() noexcept;
    ~HnfContext() noexcept;

    HnfContext(const HnfContext&) = delete;
    HnfContext& operator=(const HnfContext&) = delete;

    HnfContext(HnfContext&& other) noexcept : HnfContext() { swap(other); }
    HnfContext& operator=(HnfContext&& other) noexcept {
        if (this != &other) {
            HnfContext tmp;
            tmp.swap(*this);
            swap(other);
        }
        return *this;
    }

    void swap(HnfContext& other) noexcept;
    void set_diagnostics(const DiagnosticsContext* diagnostics) noexcept {
        diagnostics_ = diagnostics;
    }
    const DiagnosticsContext* diagnostics() const noexcept { return diagnostics_; }
    bool set(const HnfContext& other) noexcept;
    bool reset(slong ambient_dim, ulong p) noexcept;
    slong rank() const noexcept { return rank_; }
    slong ambient_dim() const noexcept { return ambient_dim_; }
    bool has_mod_rank() const noexcept { return has_mod_rank_; }
    bool add_row(bool* independent, const SparseRow& row) noexcept;
    bool add_row(bool* independent,
                 bool* index_refined,
                 const SparseRow& row) noexcept;
    bool add_fmpz_mat_row(bool* independent,
            flint::FmpzMatConstRef matrix,
            slong row) noexcept;
    bool add_fmpz_mat_row(bool* independent,
            bool* index_refined,
            flint::FmpzMatConstRef matrix,
            slong row) noexcept;
    bool add_row_defer_dependent(bool* independent,
                                 const SparseRow& row) noexcept;
    bool add_fmpz_mat_row_defer_dependent(bool* independent,
            flint::FmpzMatConstRef matrix,
            slong row) noexcept;
    bool get_hnf(flint::FmpzMatRef matrix) noexcept;
    bool get_hnf_rows(flint::FmpzMatRef matrix) noexcept;
    bool get_hnf_transform(flint::FmpzMatRef matrix) noexcept;
    bool full_rank_index(flint::FmpzRef index) noexcept;
    bool fmpz_mat_row_reduces_to_zero(bool& reduces_to_zero,
                                      flint::FmpzMatConstRef matrix,
                                      slong row) noexcept;

private:
    friend class detail::ClassRelationModuleContext;

    bool refresh(bool require_transform) noexcept;
    // Private checkpoint handoff for canonical, nonzero HNF basis rows.
    bool reset_precomputed_hnf(flint::FmpzMatConstRef hnf,
                               ulong p) noexcept;
    bool accept_independent(const SparseRow& row) noexcept;
    bool row_reduces_to_zero_in_current_lattice(
            bool& reduces_to_zero,
            const SparseRow& row) noexcept;
    bool refine_dependent(bool& index_refined, const SparseRow& row) noexcept;
    bool replace_hnf(const fmpz_mat_t hnf) noexcept;
    void invalidate_index() noexcept;
    void update_index_from_hnf() noexcept;
    bool rebuild_mod_rank() noexcept;

    const DiagnosticsContext* diagnostics_ = nullptr;
    SparseMat basis_;
    SparseMat hnf_rows_;
    SparseMat pending_dependent_rows_;
    ModRankContext mod_rank_;
    flint::FmpzMat hnf_{0, 0};
    flint::FmpzMat hnf_transform_{0, 0};
    flint::Fmpz index_;
    slong ambient_dim_ = 0;
    slong rank_ = 0;
    slong hnf_rows_processed_ = 0;
    bool has_mod_rank_ = false;
    bool hnf_valid_ = false;
    bool hnf_transform_valid_ = false;
    bool index_valid_ = false;
};

inline void swap(HnfContext& left, HnfContext& right) noexcept { left.swap(right); }

}  // namespace silex::fmpz_smat
