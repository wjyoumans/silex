#include <silex/detail/class_relation_module_context.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/fmpz_smat.hpp>

#include <initializer_list>
#include <type_traits>
#include <utility>

namespace {
namespace sflint = silex::flint;

void set_entry(sflint::FmpzMat& matrix,
               slong row,
               slong column,
               slong value) noexcept {
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(matrix, row, column), value);
}

bool entry_is(const sflint::FmpzMat& matrix,
              slong row,
              slong column,
              slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(matrix, row, column), value);
}

int test_classification_and_export() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleAddResult result;
    sflint::FmpzMat rows(5, 3);
    sflint::FmpzMat basis(3, 3);
    sflint::FmpzMat extra(1, 3);

    set_entry(rows, 0, 0, 1);
    set_entry(rows, 1, 1, 1);
    set_entry(rows, 2, 2, 1);
    set_entry(rows, 3, 0, 2);
    set_entry(rows, 3, 1, 3);

    if (!context.reset(3, 7) ||
        !context.add_fmpz_mat_row(result, rows, 0, 10) ||
        !result.modular_independent ||
        !context.add_fmpz_mat_row(result, rows, 1, 11) ||
        !result.modular_independent ||
        !context.add_fmpz_mat_row(result, rows, 2, 12) ||
        !result.modular_independent ||
        !context.add_fmpz_mat_row(result, rows, 3, 13) ||
        result.modular_independent || context.rank() != 3 ||
        context.basis_count() != 3 || context.extra_count() != 1 ||
        !context.get_basis_rows(basis) || !context.get_extra_rows(extra)) {
        return 1;
    }

    if (!entry_is(basis, 0, 0, 1) || !entry_is(basis, 1, 1, 1) ||
        !entry_is(basis, 2, 2, 1) || !entry_is(extra, 0, 0, 2) ||
        !entry_is(extra, 0, 1, 3)) {
        return 1;
    }

    if (!context.add_fmpz_mat_row(result, rows, 3, 14) ||
        result.retained || !result.duplicate ||
        context.relation_count() != 4) {
        return 1;
    }
    return 0;
}

int test_checkpoint_matches_hnf_context() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleAddResult result;
    silex::fmpz_smat::HnfContext expected;
    sflint::FmpzMat rows(3, 2);
    sflint::FmpzMat expected_hnf(2, 2);
    sflint::FmpzMat actual_hnf(2, 2);
    sflint::FmpzMat actual_transform(2, 2);
    sflint::Fmpz expected_index;
    sflint::Fmpz actual_index;
    bool independent = false;

    set_entry(rows, 0, 0, 2);
    set_entry(rows, 1, 1, 3);
    set_entry(rows, 2, 0, 1);

    if (!context.reset(2, 5) || !expected.reset(2, 5) ||
        !context.add_fmpz_mat_row(result, rows, 0, 0) ||
        !expected.add_fmpz_mat_row_defer_dependent(&independent, rows, 0) ||
        !context.add_fmpz_mat_row(result, rows, 1, 1) ||
        !expected.add_fmpz_mat_row_defer_dependent(&independent, rows, 1) ||
        !context.add_fmpz_mat_row(result, rows, 2, 2) ||
        !expected.add_fmpz_mat_row_defer_dependent(&independent, rows, 2) ||
        !context.checkpoint_index(actual_index) ||
        !expected.full_rank_index(expected_index) ||
        !context.checkpoint_hnf_rows(actual_hnf) ||
        !context.checkpoint_hnf_transform(actual_transform) ||
        !expected.get_hnf_rows(expected_hnf) ||
        !sflint::fmpz_mat_equal(actual_hnf, expected_hnf) ||
        !entry_is(actual_transform, 0, 0, 1) ||
        !entry_is(actual_transform, 0, 1, 0) ||
        !entry_is(actual_transform, 1, 0, 0) ||
        !entry_is(actual_transform, 1, 1, 1) ||
        !sflint::fmpz_equal(actual_index, expected_index) ||
        !sflint::fmpz_equal_si(actual_index, 3)) {
        return 1;
    }
    return 0;
}

int test_checkpoint_recovers_exact_rank_from_modular_dependence() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleAddResult result;
    silex::fmpz_smat::HnfContext checkpoint;
    sflint::FmpzMat rows(2, 2);
    sflint::Fmpz index;

    set_entry(rows, 0, 0, 5);
    set_entry(rows, 1, 1, 5);

    if (!context.reset(2, 5) ||
        !context.add_fmpz_mat_row(result, rows, 0, 0) ||
        result.modular_independent || !result.retained ||
        !context.add_fmpz_mat_row(result, rows, 1, 1) ||
        result.modular_independent || !result.retained ||
        context.rank() != 0 || !context.checkpoint_context(checkpoint) ||
        checkpoint.rank() != 2 || !checkpoint.full_rank_index(index) ||
        !sflint::fmpz_equal_si(index, 25)) {
        return 1;
    }
    return 0;
}

int test_checkpoint_incrementally_refines_published_lattice() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleAddResult result;
    silex::fmpz_smat::HnfContext first_checkpoint;
    silex::fmpz_smat::HnfContext second_checkpoint;
    sflint::FmpzMat rows(4, 2);
    sflint::FmpzMat hnf(2, 2);
    sflint::Fmpz first_index;
    sflint::Fmpz second_index;

    set_entry(rows, 0, 0, 2);
    set_entry(rows, 1, 1, 3);
    set_entry(rows, 2, 0, 1);
    set_entry(rows, 3, 1, 1);

    if (!context.reset(2, 5) ||
        !context.add_fmpz_mat_row(result, rows, 0, 0) ||
        !context.add_fmpz_mat_row(result, rows, 1, 1) ||
        !context.checkpoint_context(first_checkpoint) ||
        !first_checkpoint.full_rank_index(first_index) ||
        !sflint::fmpz_equal_si(first_index, 6) ||
        !context.add_fmpz_mat_row(result, rows, 2, 2) ||
        !context.add_fmpz_mat_row(result, rows, 3, 3) ||
        !context.checkpoint_context(second_checkpoint) ||
        !second_checkpoint.full_rank_index(second_index) ||
        !second_checkpoint.get_hnf_rows(hnf) ||
        !sflint::fmpz_equal_si(second_index, 1) ||
        !entry_is(hnf, 0, 0, 1) || !entry_is(hnf, 0, 1, 0) ||
        !entry_is(hnf, 1, 0, 0) || !entry_is(hnf, 1, 1, 1)) {
        return 1;
    }
    return 0;
}

// Fills `row` of `rows` with small pseudo-random entries supported on columns
// [first_column, ambient_dim), so earlier rows leave the leading columns
// empty and later rows insert pivots in front of existing ones.
void fill_phase_row(sflint::FmpzMat& rows,
                    slong row,
                    slong first_column,
                    unsigned long& state) noexcept {
    const slong cols = sflint::fmpz_mat_ncols(rows);
    for (slong j = first_column; j < cols; ++j) {
        state = state * 6364136223846793005UL + 1442695040888963407UL;
        const slong value = static_cast<slong>((state >> 33) % 13) - 6;
        set_entry(rows, row, j, value);
    }
}

// The checkpoint HNF must be the canonical row HNF of every retained row,
// however the rows were split across checkpoints.
bool checkpoint_is_canonical_hnf_of_all_rows(
        silex::detail::ClassRelationModuleContext& context) noexcept {
    silex::fmpz_smat::HnfContext checkpoint;
    if (!context.checkpoint_context(checkpoint)) {
        return false;
    }
    const slong cols = context.ambient_dim();
    sflint::FmpzMat all(context.relation_count(), cols);
    sflint::FmpzMat expected(context.relation_count(), cols);
    if (!context.get_all_rows(all)) {
        return false;
    }
    ::fmpz_mat_hnf(expected.raw(), all.raw());
    const slong rank = ::fmpz_mat_rank(all.raw());
    if (checkpoint.rank() != rank) {
        return false;
    }
    sflint::FmpzMat actual(rank, cols);
    if (!checkpoint.get_hnf_rows(actual)) {
        return false;
    }
    for (slong i = 0; i < context.relation_count(); ++i) {
        for (slong j = 0; j < cols; ++j) {
            const bool equal =
                    i < rank
                            ? ::fmpz_equal(
                                      ::fmpz_mat_entry(actual.raw(), i, j),
                                      ::fmpz_mat_entry(expected.raw(), i, j)) != 0
                            : ::fmpz_is_zero(::fmpz_mat_entry(
                                      expected.raw(), i, j)) != 0;
            if (!equal) {
                return false;
            }
        }
    }
    return true;
}

int test_checkpoints_follow_canonical_hnf_across_partial_rank() {
    constexpr slong ambient_dim = 8;
    constexpr slong row_count = 64;
    // Phases widen the support leftwards; the last phase reaches full rank.
    constexpr slong phase_first_column[] = {5, 3, 1, 0};
    sflint::FmpzMat rows(row_count, ambient_dim);
    unsigned long state = 20261010UL;
    for (slong i = 0; i < row_count; ++i) {
        fill_phase_row(rows, i, phase_first_column[(i * 4) / row_count],
                       state);
    }
    // Repeat a few rows, scaled, so dependent rows refine the lattice.
    for (slong i = 6; i < row_count; i += 7) {
        for (slong j = 0; j < ambient_dim; ++j) {
            ::fmpz_mul_si(::fmpz_mat_entry(rows.raw(), i, j),
                          ::fmpz_mat_entry(rows.raw(), i - 5, j), -3);
        }
    }

    for (const ulong prime : {ulong{3}, ulong{65537}}) {
        silex::detail::ClassRelationModuleContext context;
        silex::detail::ClassRelationModuleAddResult result;
        if (!context.reset(ambient_dim, prime)) {
            return 1;
        }
        slong since_checkpoint = 0;
        slong stride = 1;
        for (slong i = 0; i < row_count; ++i) {
            if (!context.add_fmpz_mat_row(result, rows, i, i, true)) {
                return 1;
            }
            if (++since_checkpoint < stride) {
                continue;
            }
            since_checkpoint = 0;
            stride = stride % 4 + 1;
            if (!checkpoint_is_canonical_hnf_of_all_rows(context)) {
                return 1;
            }
        }
        if (!checkpoint_is_canonical_hnf_of_all_rows(context)) {
            return 1;
        }
    }
    return 0;
}

int test_copy_and_move_contracts() {
    static_assert(!std::is_copy_constructible_v<
                  silex::detail::ClassRelationModuleContext>);
    static_assert(!std::is_copy_assignable_v<
                  silex::detail::ClassRelationModuleContext>);
    static_assert(std::is_move_constructible_v<
                  silex::detail::ClassRelationModuleContext>);
    static_assert(std::is_move_assignable_v<
                  silex::detail::ClassRelationModuleContext>);

    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleContext copy;
    silex::detail::ClassRelationModuleAddResult result;
    sflint::FmpzMat rows(2, 1);
    set_entry(rows, 0, 0, 2);
    set_entry(rows, 1, 0, 1);
    if (!context.reset(1, 3) ||
        !context.add_fmpz_mat_row(result, rows, 0, 5) ||
        !context.add_fmpz_mat_row(result, rows, 1, 6) ||
        !copy.set(context)) {
        return 1;
    }

    silex::detail::ClassRelationModuleContext moved(std::move(copy));
    return !moved.is_set() || moved.rank() != 1 ||
                   moved.relation_count() != 2
            ? 1
            : 0;
}

}  // namespace

int main() {
    return test_classification_and_export() != 0 ||
                   test_checkpoint_matches_hnf_context() != 0 ||
                   test_checkpoint_recovers_exact_rank_from_modular_dependence() !=
                           0 ||
                   test_checkpoint_incrementally_refines_published_lattice() !=
                           0 ||
                   test_checkpoints_follow_canonical_hnf_across_partial_rank() !=
                           0 ||
                   test_copy_and_move_contracts() != 0
            ? 1
            : 0;
}
