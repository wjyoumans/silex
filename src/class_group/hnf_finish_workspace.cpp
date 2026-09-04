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
           flint::fmpz_mat_nrows(workspace.relation_transform) ==
                   workspace.processed_relation_count &&
           flint::fmpz_mat_ncols(workspace.relation_transform) ==
                   workspace.processed_relation_count &&
           fmpz_mat_trimmed_nonzero_rows(
                   flint::FmpzMatConstRef(workspace.hnf)) ==
                   workspace.relation_rank;
}

bool append_hnf_finish_workspace(
        flint::FmpzMat& hnf_out,
        flint::FmpzMat& transform_out,
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
    flint::FmpzMat update_transform(update_rows, update_rows);
    ::fmpz_mat_hnf_transform(update_hnf.raw(), update_transform.raw(),
                             update_input.raw());
    const slong new_rank = fmpz_mat_trimmed_nonzero_rows(
            flint::FmpzMatConstRef(update_hnf));
    if (new_rank < 0 || new_rank != expected_rank) {
        return false;
    }

    flint::FmpzMat bridge(update_rows, relation_count);
    for (slong row = 0; row < old_rank; ++row) {
        for (slong col = 0; col < old_relations; ++col) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatRef(bridge), row, col),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(
                                    workspace.relation_transform),
                            row, col));
        }
    }
    for (slong row = 0; row < new_relation_count; ++row) {
        flint::fmpz_one(flint::fmpz_mat_entry(
                flint::FmpzMatRef(bridge), old_rank + row,
                old_relations + row));
    }

    flint::FmpzMat composed(update_rows, relation_count);
    flint::fmpz_mat_mul(flint::FmpzMatRef(composed),
                        flint::FmpzMatConstRef(update_transform),
                        flint::FmpzMatConstRef(bridge));

    flint::FmpzMat candidate_hnf(relation_count, generator_count);
    flint::FmpzMat candidate_transform(relation_count, relation_count);
    for (slong row = 0; row < new_rank; ++row) {
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(candidate_hnf), row,
                    flint::FmpzMatConstRef(update_hnf), row) ||
            !fmpz_mat_copy_row(
                    flint::FmpzMatRef(candidate_transform), row,
                    flint::FmpzMatConstRef(composed), row)) {
            return false;
        }
    }

    slong output_row = new_rank;
    for (slong row = old_rank; row < old_relations; ++row, ++output_row) {
        for (slong col = 0; col < old_relations; ++col) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatRef(candidate_transform),
                            output_row, col),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(
                                    workspace.relation_transform),
                            row, col));
        }
    }
    for (slong row = new_rank; row < update_rows; ++row, ++output_row) {
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(candidate_transform), output_row,
                    flint::FmpzMatConstRef(composed), row)) {
            return false;
        }
    }
    if (output_row != relation_count) {
        return false;
    }

    hnf_out = std::move(candidate_hnf);
    transform_out = std::move(candidate_transform);
    rank_out = new_rank;
    return true;
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
    flint::FmpzMat candidate_transform(0, 0);
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
                    candidate_hnf, candidate_transform, candidate_rank,
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
        candidate_transform = flint::FmpzMat(relation_count,
                                              relation_count);
        ::fmpz_mat_hnf_transform(candidate_hnf.raw(),
                                 candidate_transform.raw(),
                                 relations.raw());
        candidate_rank = fmpz_mat_trimmed_nonzero_rows(
                flint::FmpzMatConstRef(candidate_hnf));
        if (candidate_rank < 0 ||
            candidate_rank != context.relation_rank()) {
            return false;
        }
    }

    workspace.hnf = std::move(candidate_hnf);
    workspace.relation_transform = std::move(candidate_transform);
    workspace.generator_count = generator_count;
    workspace.processed_relation_count = relation_count;
    workspace.relation_rank = candidate_rank;
    workspace.exact_valid = true;
    workspace.unit_logs = flint::ArbMat(0, 0);
    workspace.unit_log_precision = 0;
    workspace.unit_logs_valid = false;
    if (!incremental) {
        workspace.relation_logs = flint::ArbMat(0, 0);
        workspace.relation_log_precision = 0;
        workspace.logged_relation_count = 0;
        workspace.relation_logs_valid = false;
    }
    return true;
}

bool hnf_finish_workspace_witness_coefficients(
        flint::FmpzMat& out,
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept {
    if (!sync_hnf_finish_workspace(workspace, context)) {
        return false;
    }

    const slong relation_count = workspace.processed_relation_count;
    const slong rank = workspace.relation_rank;
    flint::FmpzMat coefficients(relation_count - rank, relation_count);
    for (slong row = rank; row < relation_count; ++row) {
        if (!fmpz_mat_copy_row(
                    flint::FmpzMatRef(coefficients), row - rank,
                    flint::FmpzMatConstRef(workspace.relation_transform),
                    row)) {
            return false;
        }
    }
    out = std::move(coefficients);
    return true;
}

}  // namespace silex::detail
