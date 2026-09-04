#include "order_unit_internal.hpp"

#include "../class_group/class_group_internal.hpp"

#include <utility>

#include <flint/fmpz_mat.h>

#include <silex/archimedean.hpp>

namespace silex::detail {

namespace {

slong trimmed_nonzero_rows(flint::FmpzMatConstRef matrix) noexcept {
    slong rows = flint::fmpz_mat_nrows(matrix);
    while (rows > 0 && ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (::fmpz_mat_is_zero_row(matrix.raw(), i) != 0) {
            return -1;
        }
    }
    return rows;
}

bool hnf_finish_workspace_unit_log_matrix(
        flint::ArbMat& out,
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong precision) noexcept {
    if (!sync_hnf_finish_workspace(workspace, class_group)) {
        return false;
    }

    const slong relation_count = workspace.processed_relation_count;
    const slong unit_columns = relation_count - workspace.relation_rank;
    slong places = 0;
    if (!compact_places(places, embeddings) ||
        flint::arb_mat_nrows_value(out) != unit_columns ||
        flint::arb_mat_ncols_value(out) != places) {
        return false;
    }

    slong first_new_relation = 0;
    if (workspace.relation_logs_valid) {
        if (workspace.relation_log_precision == precision) {
            if (workspace.logged_relation_count < 0 ||
                workspace.logged_relation_count > relation_count ||
                flint::arb_mat_nrows_value(workspace.relation_logs) !=
                        workspace.logged_relation_count ||
                flint::arb_mat_ncols_value(workspace.relation_logs) !=
                        places) {
                return false;
            }
            first_new_relation = workspace.logged_relation_count;
        }
    }

    if (!workspace.relation_logs_valid ||
        workspace.relation_log_precision != precision ||
        first_new_relation != relation_count) {
        flint::ArbMat relation_logs(relation_count, places);
        for (slong row = 0; row < first_new_relation; ++row) {
            for (slong place = 0; place < places; ++place) {
                flint::arb_set(
                        flint::arb_mat_entry_ref(relation_logs, row, place),
                        flint::arb_mat_entry_ref(
                                flint::ArbMatConstRef(
                                        workspace.relation_logs),
                                row, place));
            }
        }

        const NumberField* field = class_group.parent() == nullptr
                ? nullptr
                : class_group.parent()->parent();
        if (field == nullptr) {
            return false;
        }
        flint::ArbVec row_logs(places);
        Element generator(*field);
        for (slong row = first_new_relation; row < relation_count; ++row) {
            if (!class_group.relation_generator(generator, row) ||
                !silex::logarithmic_embedding(
                        flint::ArbVecRef(row_logs), embeddings, generator,
                        LogEmbeddingMode::product, precision)) {
                return false;
            }
            for (slong place = 0; place < places; ++place) {
                flint::arb_set(
                        flint::arb_mat_entry_ref(relation_logs, row, place),
                        flint::ArbConstRef(row_logs.data() + place));
            }
        }

        workspace.relation_logs = std::move(relation_logs);
        workspace.relation_log_precision = precision;
        workspace.logged_relation_count = relation_count;
        workspace.relation_logs_valid = true;
        workspace.unit_logs = flint::ArbMat(0, 0);
        workspace.unit_log_precision = 0;
        workspace.unit_logs_valid = false;
    }

    if (workspace.unit_logs_valid) {
        if (workspace.unit_log_precision != precision ||
            flint::arb_mat_nrows_value(workspace.unit_logs) != unit_columns ||
            flint::arb_mat_ncols_value(workspace.unit_logs) != places) {
            return false;
        }
        flint::arb_mat_set(flint::ArbMatRef(out),
                           flint::ArbMatConstRef(workspace.unit_logs));
        return true;
    }

    flint::ArbMat unit_logs(unit_columns, places);
    for (slong dep = 0; dep < unit_columns; ++dep) {
        const slong transform_row = workspace.relation_rank + dep;
        for (slong relation = 0; relation < relation_count; ++relation) {
            const auto exponent = flint::fmpz_mat_entry(
                    flint::FmpzMatConstRef(workspace.relation_transform),
                    transform_row, relation);
            if (flint::fmpz_is_zero(exponent)) {
                continue;
            }
            for (slong place = 0; place < places; ++place) {
                flint::arb_addmul_fmpz(
                        flint::arb_mat_entry_ref(unit_logs, dep, place),
                        flint::arb_mat_entry_ref(
                                flint::ArbMatConstRef(
                                        workspace.relation_logs),
                                relation, place),
                        exponent, precision);
            }
        }
    }

    workspace.unit_logs = std::move(unit_logs);
    workspace.unit_log_precision = precision;
    workspace.unit_logs_valid = true;
    flint::arb_mat_set(flint::ArbMatRef(out),
                       flint::ArbMatConstRef(workspace.unit_logs));
    return true;
}

}  // namespace

bool hnf_unit_log_matrix(flint::ArbMat& out,
                              const ClassGroupContext& class_group,
                              EmbeddingContext& embeddings,
                              HnfFinishWorkspace* workspace,
                              slong precision) noexcept {
    const Order* order = class_group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !class_group.has_factor_base() ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*field) ||
        precision <= 0 ||
        class_group.relation_count() < class_group.relation_rank() ||
        class_group.generator_count() <= 0 ||
        class_group.relation_rank() > class_group.generator_count()) {
        return false;
    }

    if (workspace != nullptr) {
        return hnf_finish_workspace_unit_log_matrix(
                out, *workspace, class_group, embeddings, precision);
    }

    const slong relation_count = class_group.relation_count();
    const slong generator_count = class_group.generator_count();
    flint::FmpzMat stored(relation_count, generator_count);
    if (!class_group.relations(flint::FmpzMatRef(stored))) {
        return false;
    }

    flint::FmpzMat hnf(relation_count, generator_count);
    flint::FmpzMat transform(relation_count, relation_count);
    ::fmpz_mat_hnf_transform(hnf.raw(), transform.raw(), stored.raw());

    const slong rank = trimmed_nonzero_rows(flint::FmpzMatConstRef(hnf));
    if (rank < 0 || rank != class_group.relation_rank()) {
        return false;
    }

    slong places = 0;
    const slong unit_columns = relation_count - rank;
    if (!compact_places(places, embeddings) ||
        flint::arb_mat_nrows_value(out) != unit_columns ||
        flint::arb_mat_ncols_value(out) != places) {
        return false;
    }

    flint::ArbMat relation_logs(relation_count, places);
    flint::ArbVec row(places);
    Element generator(*field);
    for (slong i = 0; i < relation_count; ++i) {
        if (!class_group.relation_generator(generator, i) ||
            !silex::logarithmic_embedding(flint::ArbVecRef(row), embeddings,
                                          generator,
                                          LogEmbeddingMode::product,
                                          precision)) {
            return false;
        }

        for (slong j = 0; j < places; ++j) {
            arb_set(arb_mat_entry(relation_logs.raw(), i, j),
                    row.data() + j);
        }
    }

    for (slong i = 0; i < unit_columns; ++i) {
        for (slong j = 0; j < places; ++j) {
            flint::arb_zero(flint::arb_mat_entry_ref(out, i, j));
        }
    }

    // reference relation_completion_parameters carries HNF column operations in C and checks
    // Ar = real_i(A) for the zero/HNF-unit columns.  With FLINT's
    // row-oriented hnf = transform * relation_rows, the zero HNF rows are
    // the corresponding retained-relation dependencies.
    for (slong dep = 0; dep < unit_columns; ++dep) {
        const slong transform_row = rank + dep;
        for (slong relation = 0; relation < relation_count; ++relation) {
            const auto exponent = flint::fmpz_mat_entry(
                    flint::FmpzMatConstRef(transform), transform_row,
                    relation);
            if (flint::fmpz_is_zero(exponent)) {
                continue;
            }

            for (slong place = 0; place < places; ++place) {
                flint::arb_addmul_fmpz(
                        flint::arb_mat_entry_ref(out, dep, place),
                        flint::arb_mat_entry_ref(
                                flint::ArbMatConstRef(relation_logs),
                                relation, place),
                        exponent, precision);
            }
        }
    }

    return true;
}

}  // namespace silex::detail
