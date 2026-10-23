#include "class_group_internal.hpp"

#include <utility>

namespace silex::detail {

namespace {

bool hnf_finish_workspace_shape_matches(
        const HnfFinishWorkspace& workspace) noexcept {
    return workspace.exact_valid && workspace.generator_count > 0 &&
           workspace.processed_relation_count >= workspace.relation_rank &&
           workspace.relation_rank >= 0 &&
           workspace.relation_rank <= workspace.generator_count;
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

    // The exact rank check: the rank of the stored relations must be the
    // context's relation rank, which can be a modular rank or stale when
    // the context's exact row module has not been synced.  fmpz_mat_rank
    // is the number of nonzero rows of the HNF of the relations (unimodular
    // row operations keep the rank, and the nonzero rows of an echelon form
    // are independent), so this is the check the workspace HNF made before,
    // without the HNF.  On failure the workspace is left unchanged.  A sync
    // with no new relations keeps the rank checked here.
    flint::FmpzMat relations(relation_count, generator_count);
    if (!context.relations(flint::FmpzMatRef(relations)) ||
        ::fmpz_mat_rank(relations.raw()) != context.relation_rank()) {
        return false;
    }

    const bool incremental = workspace.exact_valid;
    workspace.generator_count = generator_count;
    workspace.processed_relation_count = relation_count;
    workspace.relation_rank = context.relation_rank();
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
