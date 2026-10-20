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
    // Repeat a few rows, scaled: these lie in the lattice already, so they
    // reduce to zero rows with negative entries.
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

// Checks the published full-rank checkpoint after each step: its index, that
// its HNF is the canonical HNF of every retained row, and that its transform
// is the identity on those HNF rows.
bool full_rank_checkpoint_is(
        silex::detail::ClassRelationModuleContext& context,
        slong expected_index) noexcept {
    sflint::Fmpz index;
    sflint::FmpzMat transform(context.ambient_dim(), context.ambient_dim());
    sflint::FmpzMat one(context.ambient_dim(), context.ambient_dim());
    sflint::fmpz_mat_one(one);
    return context.checkpoint_index(index) &&
           sflint::fmpz_equal_si(index, expected_index) &&
           context.checkpoint_hnf_transform(transform) &&
           sflint::fmpz_mat_equal(transform, one) &&
           checkpoint_is_canonical_hnf_of_all_rows(context);
}

// Once the checkpoint has full rank, extra rows inside its lattice leave it
// unchanged, while a batch with one row outside the lattice (after or before
// rows inside it) refines it to the HNF of all rows.
int test_full_rank_checkpoint_over_lattice_and_new_rows() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleAddResult result;
    sflint::FmpzMat rows(10, 2);
    const slong entries[10][2] = {{2, 0}, {0, 3}, {4, 3}, {2, 6}, {-6, 9},
                                  {2, 3}, {4, 1}, {0, 7}, {1, 0}, {3, 5}};
    for (slong i = 0; i < 10; ++i) {
        set_entry(rows, i, 0, entries[i][0]);
        set_entry(rows, i, 1, entries[i][1]);
    }

    if (!context.reset(2, 5) ||
        !context.add_fmpz_mat_row(result, rows, 0, 0) ||
        !context.add_fmpz_mat_row(result, rows, 1, 1) ||
        !full_rank_checkpoint_is(context, 6) ||
        // One row, then two rows, inside the lattice of index 6.
        !context.add_fmpz_mat_row(result, rows, 2, 2) ||
        !result.retained || !full_rank_checkpoint_is(context, 6) ||
        !context.add_fmpz_mat_row(result, rows, 3, 3) ||
        !context.add_fmpz_mat_row(result, rows, 4, 4) ||
        !full_rank_checkpoint_is(context, 6) ||
        // A row inside, then one outside: index 6 becomes 2.
        !context.add_fmpz_mat_row(result, rows, 5, 5) ||
        !context.add_fmpz_mat_row(result, rows, 6, 6) ||
        !full_rank_checkpoint_is(context, 2) ||
        // A row inside the refined lattice, then a batch whose first row
        // lies outside it.
        !context.add_fmpz_mat_row(result, rows, 7, 7) ||
        !full_rank_checkpoint_is(context, 2) ||
        !context.add_fmpz_mat_row(result, rows, 8, 8) ||
        !context.add_fmpz_mat_row(result, rows, 9, 9) ||
        !full_rank_checkpoint_is(context, 1) ||
        context.relation_count() != 10) {
        return 1;
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

bool add_rows(silex::detail::ClassRelationModuleContext& context,
              const sflint::FmpzMat& rows,
              slong count) noexcept {
    silex::detail::ClassRelationModuleAddResult result;
    for (slong i = 0; i < count; ++i) {
        if (!context.add_fmpz_mat_row(result, rows, i, i)) {
            return false;
        }
    }
    return true;
}

// An in-place add followed by rollback_add restores the module exactly: the
// retained rows and indices, the modular basis (including a middle insertion
// that eliminated a column of an earlier basis row), and the checkpoint
// bookkeeping.  Adds that retain nothing leave the module unchanged.
int test_rollback_restores_state_before_add() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleContext snapshot;
    silex::detail::ClassRelationModuleContext reference;
    silex::detail::ClassRelationModuleAddResult result;
    silex::detail::ClassRelationModuleAddUndo undo;
    silex::fmpz_smat::HnfContext checkpoint;
    sflint::FmpzMat rows(6, 4);

    // r0, r1 have modular pivots 0 and 2; r2 = 2 r0 is an extra row.
    set_entry(rows, 0, 0, 1);
    set_entry(rows, 0, 1, 3);
    set_entry(rows, 0, 3, 4);
    set_entry(rows, 1, 2, 1);
    set_entry(rows, 1, 3, 2);
    set_entry(rows, 2, 0, 2);
    set_entry(rows, 2, 1, 6);
    set_entry(rows, 2, 3, 8);
    // r3 has pivot 1 and is inserted between the two basis rows; it
    // eliminates the entry 3 in column 1 of r0.
    set_entry(rows, 3, 1, 1);
    set_entry(rows, 3, 3, 5);
    // r4 = 3 r0 is modular dependent but not a duplicate.
    set_entry(rows, 4, 0, 3);
    set_entry(rows, 4, 1, 9);
    set_entry(rows, 4, 3, 12);

    if (!context.reset(4, 7) || !add_rows(context, rows, 3) ||
        !context.checkpoint_context(checkpoint) || !snapshot.set(context) ||
        !context.same_state_as(snapshot)) {
        return 1;
    }

    // Independent row, rejected after the add.
    if (!context.add_fmpz_mat_row(undo, result, rows, 3, 3, false) ||
        !result.retained || !result.modular_independent || !undo.changed ||
        !undo.basis_row || undo.mod_rank_position != 1 ||
        context.rank() != 3 || context.same_state_as(snapshot) ||
        !context.rollback_add(undo) || undo.changed ||
        !context.same_state_as(snapshot) || context.rank() != 2 ||
        context.basis_count() != 2 || context.extra_count() != 1) {
        return 1;
    }

    // Modular-dependent row retained as an extra row, then rolled back.
    if (!context.add_fmpz_mat_row(undo, result, rows, 4, 4, true) ||
        !result.retained || result.modular_independent || !undo.changed ||
        undo.basis_row || context.extra_count() != 2 ||
        !context.rollback_add(undo) || !context.same_state_as(snapshot)) {
        return 1;
    }

    // Adds that retain nothing change nothing; rollback is a no-op.
    if (!context.add_fmpz_mat_row(undo, result, rows, 4, 4, false) ||
        result.retained || undo.changed ||
        !context.add_fmpz_mat_row(undo, result, rows, 1, 4, true) ||
        !result.duplicate || undo.changed ||
        !context.rollback_add(undo) || !context.same_state_as(snapshot)) {
        return 1;
    }

    // After the rollbacks the module evolves as if they never happened.
    if (!reference.reset(4, 7) || !add_rows(reference, rows, 3) ||
        !reference.checkpoint_context(checkpoint) ||
        !reference.add_fmpz_mat_row(result, rows, 3, 3) ||
        !context.add_fmpz_mat_row(result, rows, 3, 3) ||
        !context.same_state_as(reference)) {
        return 1;
    }

    // A pending record refuses a second add; a checkpoint rebuilt over the
    // added row refuses the rollback.  Neither changes the module.
    silex::detail::ClassRelationModuleContext after;
    if (!context.add_fmpz_mat_row(undo, result, rows, 4, 4, true) ||
        !undo.changed || !after.set(context) ||
        context.add_fmpz_mat_row(undo, result, rows, 5, 5, true) ||
        !context.same_state_as(after) ||
        !context.checkpoint_context(checkpoint) || !after.set(context) ||
        context.rollback_add(undo) || !context.same_state_as(after)) {
        return 1;
    }
    return 0;
}

// A checkpoint taken between an add and its rollback leaves the restored
// module with an invalid checkpoint whose counts are zero, as if the
// checkpoint had been taken before the add; the rollback never marks a
// checkpoint valid that does not cover the remaining rows.
int test_rollback_after_take_checkpoint_keeps_checkpoint_invalid() {
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleContext reference;
    silex::detail::ClassRelationModuleAddResult result;
    silex::detail::ClassRelationModuleAddUndo undo;
    silex::fmpz_smat::HnfContext checkpoint;
    sflint::FmpzMat rows(4, 4);

    set_entry(rows, 0, 0, 1);
    set_entry(rows, 0, 1, 3);
    set_entry(rows, 0, 3, 4);
    set_entry(rows, 1, 2, 1);
    set_entry(rows, 1, 3, 2);
    set_entry(rows, 2, 0, 2);
    set_entry(rows, 2, 1, 6);
    set_entry(rows, 2, 3, 8);
    set_entry(rows, 3, 1, 1);
    set_entry(rows, 3, 3, 5);

    if (!context.reset(4, 7) || !add_rows(context, rows, 3) ||
        !context.checkpoint_context(checkpoint) ||
        !context.add_fmpz_mat_row(undo, result, rows, 3, 3, false) ||
        !undo.changed || !undo.checkpoint_valid ||
        !context.take_checkpoint_context(checkpoint) ||
        !context.rollback_add(undo)) {
        return 1;
    }
    if (!reference.reset(4, 7) || !add_rows(reference, rows, 3) ||
        !reference.checkpoint_context(checkpoint) ||
        !reference.take_checkpoint_context(checkpoint) ||
        !context.same_state_as(reference)) {
        return 1;
    }
    return 0;
}

// Every in-place add, at whatever modular pivot position it lands, rolls back
// to the exact prior state; the module then takes the row for real.
int test_rollback_after_each_add_of_a_row_stream() {
    constexpr slong kColumns = 6;
    constexpr slong kRows = 24;
    silex::detail::ClassRelationModuleContext context;
    silex::detail::ClassRelationModuleContext snapshot;
    silex::detail::ClassRelationModuleAddResult result;
    sflint::FmpzMat rows(kRows, kColumns);
    ulong state = 12345;
    for (slong i = 0; i < kRows; ++i) {
        for (slong j = 0; j < kColumns; ++j) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            const slong value = static_cast<slong>((state >> 33) % 7) - 3;
            set_entry(rows, i, j, (j * 5 + i) % 4 == 0 ? value : 0);
        }
    }
    if (!context.reset(kColumns, 5)) {
        return 1;
    }
    slong independent_rollbacks = 0;
    for (slong i = 0; i < kRows; ++i) {
        silex::detail::ClassRelationModuleAddUndo undo;
        if (!snapshot.set(context) ||
            !context.add_fmpz_mat_row(undo, result, rows, i, i, true) ||
            !context.rollback_add(undo) || !context.same_state_as(snapshot) ||
            !context.add_fmpz_mat_row(result, rows, i, i)) {
            return 1;
        }
        independent_rollbacks += result.modular_independent ? 1 : 0;
    }
    return independent_rollbacks >= 3 ? 0 : 1;
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
                   test_full_rank_checkpoint_over_lattice_and_new_rows() !=
                           0 ||
                   test_copy_and_move_contracts() != 0 ||
                   test_rollback_restores_state_before_add() != 0 ||
                   test_rollback_after_each_add_of_a_row_stream() != 0 ||
                   test_rollback_after_take_checkpoint_keeps_checkpoint_invalid() !=
                           0
            ? 1
            : 0;
}
