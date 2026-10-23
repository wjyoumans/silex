#include "class_group_internal.hpp"

#include <utility>

namespace silex::detail {

namespace {

slong fmpz_mat_trimmed_nonzero_rows(flint::FmpzMatConstRef matrix) noexcept {
    slong rows = flint::fmpz_mat_nrows(matrix);
    while (rows > 0 &&
           ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (::fmpz_mat_is_zero_row(matrix.raw(), i) != 0) {
            return -1;
        }
    }
    return rows;
}

bool fmpz_mat_copy_row(flint::FmpzMatRef out,
                       slong out_row,
                       flint::FmpzMatConstRef in,
                       slong in_row) noexcept {
    if (out_row < 0 || out_row >= flint::fmpz_mat_nrows(out) ||
        in_row < 0 || in_row >= flint::fmpz_mat_nrows(in) ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(in)) {
        return false;
    }
    for (slong j = 0; j < flint::fmpz_mat_ncols(in); ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, out_row, j),
                        flint::fmpz_mat_entry(in, in_row, j));
    }
    return true;
}

bool hnf_finish_workspace_shape_matches(
        const HnfFinishWorkspace& workspace) noexcept {
    return workspace.exact_valid && workspace.generator_count > 0 &&
           workspace.processed_relation_count >= workspace.relation_rank &&
           workspace.relation_rank >= 0 &&
           workspace.relation_rank <= workspace.generator_count &&
           flint::fmpz_mat_nrows(workspace.hnf) ==
                   workspace.processed_relation_count &&
           flint::fmpz_mat_ncols(workspace.hnf) ==
                   workspace.generator_count &&
           fmpz_mat_trimmed_nonzero_rows(
                   flint::FmpzMatConstRef(workspace.hnf)) ==
                   workspace.relation_rank;
}

// The HNF of the old relations followed by new_relations, from the HNF of
// the update matrix: the old HNF's nonzero rows, which span the old
// relations' row lattice, followed by the new relations.  That matrix spans
// the row lattice of all the relations, so its HNF, padded with zero rows,
// is theirs by uniqueness of the HNF of a row lattice (Cohen, Theorem
// 2.4.3).
bool append_hnf_finish_workspace(
        flint::FmpzMat& hnf_out,
        slong& rank_out,
        const HnfFinishWorkspace& workspace,
        flint::FmpzMatConstRef new_relations,
        slong expected_rank) noexcept {
    const slong old_relations = workspace.processed_relation_count;
    const slong old_rank = workspace.relation_rank;
    const slong new_relation_count =
            flint::fmpz_mat_nrows(new_relations);
    const slong generator_count = flint::fmpz_mat_ncols(new_relations);
    if (!hnf_finish_workspace_shape_matches(workspace) ||
        generator_count != workspace.generator_count ||
        new_relation_count <= 0 || expected_rank < old_rank ||
        expected_rank > generator_count) {
        return false;
    }

    if (old_relations > WORD_MAX - new_relation_count ||
        old_rank > WORD_MAX - new_relation_count) {
        return false;
    }
    const slong relation_count = old_relations + new_relation_count;
    const slong update_rows = old_rank + new_relation_count;
    flint::FmpzMat update_input(update_rows, generator_count);
    for (slong row = 0; row < old_rank; ++row) {
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(update_input), row,
                    flint::FmpzMatConstRef(workspace.hnf), row)) {
            return false;
        }
    }
    for (slong row = 0; row < new_relation_count; ++row) {
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(update_input), old_rank + row,
                    new_relations, row)) {
            return false;
        }
    }

    flint::FmpzMat update_hnf(update_rows, generator_count);
    ::fmpz_mat_hnf(update_hnf.raw(), update_input.raw());
    const slong new_rank = fmpz_mat_trimmed_nonzero_rows(
            flint::FmpzMatConstRef(update_hnf));
    if (new_rank < 0 || new_rank != expected_rank) {
        return false;
    }

    flint::FmpzMat candidate_hnf(relation_count, generator_count);
    for (slong row = 0; row < new_rank; ++row) {
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(candidate_hnf), row,
                    flint::FmpzMatConstRef(update_hnf), row)) {
            return false;
        }
    }

    hnf_out = std::move(candidate_hnf);
    rank_out = new_rank;
    return true;
}

// The n x (g + n) input of the finish kernel's reduction.  Without a usable
// cached extension it is [c * relations | I_n]
// (relation_kernel_extension_input), the input of
// finish_relation_kernel_basis.  With one, of the first m < n relations,
// its first m rows are the cached reduced rows, padded with zeros in the
// identity columns of the new relations, and row i >= m is
// [c * relation_i | e_i].  The cached rows are the rows of
// U * [c * relations_0..m-1 | I_m] for the unimodular U of the earlier
// reduction, so the input is blockdiag(U, I_{n-m}) * [c * relations | I_n]:
// the same lattice, with a unimodular identity part.
bool finish_extension_input(flint::FmpzMat& out,
                            bool& extended,
                            const HnfFinishWorkspace& workspace,
                            flint::FmpzMatConstRef relations,
                            bool use_cache) noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    const slong g = flint::fmpz_mat_ncols(relations);
    if (n <= 0 || g <= 0 || n > WORD_MAX - g) {
        return false;
    }
    const slong m = workspace.finish_extension_relation_count;
    extended = use_cache && workspace.finish_extension_valid && m > 0 &&
               m < n &&
               flint::fmpz_mat_nrows(workspace.finish_extension) == m &&
               flint::fmpz_mat_ncols(workspace.finish_extension) == g + m;

    flint::FmpzMat input(0, 0);
    if (!relation_kernel_extension_input(input, relations)) {
        return false;
    }
    // Rows i < m of [c * relations | I_n] are zero in the columns from
    // g + m on, so overwriting their first g + m entries with the cached
    // rows gives the padded cached rows.
    const slong cached_rows = extended ? m : 0;
    for (slong i = 0; i < cached_rows; ++i) {
        for (slong j = 0; j < g + m; ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(flint::FmpzMatRef(input), i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(
                                    workspace.finish_extension),
                            i, j));
        }
    }
    out = std::move(input);
    return true;
}

// The finish kernel through finish_extension_input and
// relation_kernel_phase_rows; on success the reduced extension replaces
// the workspace's cached one.  `use_cache` false starts from
// [c * relations | I_n].  `extended` reports whether the cached extension
// was used.
bool finish_extension_kernel_rows(flint::FmpzMat& out,
                                  bool& extended,
                                  HnfFinishWorkspace& workspace,
                                  flint::FmpzMatConstRef relations,
                                  slong rank,
                                  bool use_cache) noexcept {
    extended = false;
    flint::FmpzMat extension(0, 0);
    flint::FmpzMat rows(0, 0);
    if (!finish_extension_input(extension, extended, workspace, relations,
                                use_cache) ||
        !relation_kernel_phase_rows(rows, extension, relations, rank)) {
        return false;
    }
    workspace.finish_extension = std::move(extension);
    workspace.finish_extension_relation_count =
            flint::fmpz_mat_nrows(relations);
    workspace.finish_extension_valid = true;
    out = std::move(rows);
    return true;
}

// The finish kernel of the workspace's relations.  When the workspace keeps
// the reduced extension of an earlier relation prefix, the reduction starts
// from it (finish_extension_input); when that input fails the phase or its
// checks, the reduction starts from [c * relations | I_n], whose accepted
// rows are finish_relation_kernel_basis's.  When both fail, or the kernel
// is empty, the rows are hnf_unit_finish_coefficients's, and no extension
// is kept.  Every accepted basis is a basis of the same kernel lattice;
// only the basis depends on the relation counts at which earlier finish
// kernels were taken.
bool hnf_finish_workspace_finish_rows(flint::FmpzMat& out,
                                      HnfFinishWorkspace& workspace,
                                      const ClassGroupContext& context)
        noexcept {
    const slong n = workspace.processed_relation_count;
    const slong rank = workspace.relation_rank;
    if (rank >= 0 && rank < n) {
        flint::FmpzMat relations(n, workspace.generator_count);
        // The rank check of hnf_unit_finish_coefficients: the exact rank
        // of the stored relations must be the context's relation rank,
        // which checked_relation_kernel_rows takes as its precondition.
        if (!context.relations(flint::FmpzMatRef(relations)) ||
            ::fmpz_mat_rank(relations.raw()) != rank) {
            return false;
        }
        bool extended = false;
        if (finish_extension_kernel_rows(
                    out, extended, workspace,
                    flint::FmpzMatConstRef(relations), rank, true)) {
            return true;
        }
        if (extended &&
            finish_extension_kernel_rows(
                    out, extended, workspace,
                    flint::FmpzMatConstRef(relations), rank, false)) {
            return true;
        }
    }

    workspace.finish_extension = flint::FmpzMat(0, 0);
    workspace.finish_extension_relation_count = 0;
    workspace.finish_extension_valid = false;
    return hnf_unit_finish_coefficients(out, context);
}

}  // namespace

bool sync_hnf_finish_workspace(
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept {
    const Order* order = context.parent();
    if (!context.has_factor_base() || order == nullptr ||
        order->parent() == nullptr || context.generator_count() <= 0 ||
        context.relation_count() < context.relation_rank() ||
        context.relation_rank() < 0 ||
        context.relation_rank() > context.generator_count()) {
        return false;
    }

    const slong relation_count = context.relation_count();
    const slong generator_count = context.generator_count();
    if (workspace.exact_valid) {
        if (!hnf_finish_workspace_shape_matches(workspace) ||
            workspace.generator_count != generator_count ||
            workspace.processed_relation_count > relation_count ||
            workspace.relation_rank > context.relation_rank()) {
            return false;
        }
        if (workspace.processed_relation_count == relation_count) {
            return workspace.relation_rank == context.relation_rank();
        }
    }

    flint::FmpzMat candidate_hnf(0, 0);
    slong candidate_rank = -1;
    const bool incremental = workspace.exact_valid;
    if (incremental) {
        flint::FmpzMat new_relations(
                relation_count - workspace.processed_relation_count,
                generator_count);
        if (!class_group_relation_rows(
                    flint::FmpzMatRef(new_relations), context,
                    workspace.processed_relation_count)) {
            return false;
        }
        if (!append_hnf_finish_workspace(
                    candidate_hnf, candidate_rank,
                    workspace, flint::FmpzMatConstRef(new_relations),
                    context.relation_rank())) {
            return false;
        }
    } else {
        flint::FmpzMat relations(relation_count, generator_count);
        if (!context.relations(flint::FmpzMatRef(relations))) {
            return false;
        }
        candidate_hnf = flint::FmpzMat(relation_count, generator_count);
        ::fmpz_mat_hnf(candidate_hnf.raw(), relations.raw());
        candidate_rank = fmpz_mat_trimmed_nonzero_rows(
                flint::FmpzMatConstRef(candidate_hnf));
        if (candidate_rank < 0 ||
            candidate_rank != context.relation_rank()) {
            return false;
        }
    }

    workspace.hnf = std::move(candidate_hnf);
    workspace.generator_count = generator_count;
    workspace.processed_relation_count = relation_count;
    workspace.relation_rank = candidate_rank;
    workspace.exact_valid = true;
    workspace.unit_logs = flint::ArbMat(0, 0);
    workspace.unit_log_precision = 0;
    workspace.unit_logs_valid = false;
    workspace.unit_logs_finish_kernel = false;
    workspace.finish_kernel = flint::FmpzMat(0, 0);
    workspace.finish_kernel_valid = false;
    workspace.witness_kernel = flint::FmpzMat(0, 0);
    workspace.witness_kernel_valid = false;
    if (!incremental) {
        workspace.finish_extension = flint::FmpzMat(0, 0);
        workspace.finish_extension_relation_count = 0;
        workspace.finish_extension_valid = false;
        workspace.relation_logs = flint::ArbMat(0, 0);
        workspace.relation_log_precision = 0;
        workspace.logged_relation_count = 0;
        workspace.relation_logs_valid = false;
    }
    return true;
}

const flint::FmpzMat* hnf_finish_workspace_kernel(
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context,
        bool finish) noexcept {
    if (!sync_hnf_finish_workspace(workspace, context)) {
        return nullptr;
    }

    flint::FmpzMat& kernel =
            finish ? workspace.finish_kernel : workspace.witness_kernel;
    bool& valid = finish ? workspace.finish_kernel_valid
                         : workspace.witness_kernel_valid;
    if (!valid) {
        flint::FmpzMat rows(0, 0);
        const bool have_rows =
                finish ? hnf_finish_workspace_finish_rows(rows, workspace,
                                                          context)
                       : hnf_unit_witness_coefficients(rows, context);
        if (!have_rows ||
            flint::fmpz_mat_nrows(rows) !=
                    workspace.processed_relation_count -
                            workspace.relation_rank ||
            flint::fmpz_mat_ncols(rows) !=
                    workspace.processed_relation_count) {
            return nullptr;
        }
        kernel = std::move(rows);
        valid = true;
    }
    return &kernel;
}

bool hnf_finish_workspace_witness_coefficients(
        flint::FmpzMat& out,
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept {
    const flint::FmpzMat* const kernel =
            hnf_finish_workspace_kernel(workspace, context, false);
    if (kernel == nullptr) {
        return false;
    }
    flint::FmpzMat coefficients(flint::fmpz_mat_nrows(*kernel),
                                flint::fmpz_mat_ncols(*kernel));
    flint::fmpz_mat_set(flint::FmpzMatRef(coefficients),
                        flint::FmpzMatConstRef(*kernel));
    out = std::move(coefficients);
    return true;
}

}  // namespace silex::detail
