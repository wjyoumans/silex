#include "class_group_internal.hpp"

#include <utility>

#include <flint/fmpz_lll.h>

#include <silex/flint/fmpz_lll.hpp>

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

struct ClassRelationHnfTransform {
    flint::FmpzMat stored{0, 0};
    flint::FmpzMat hnf{0, 0};
    flint::FmpzMat transform{0, 0};
    slong rank = 0;
};

// The context's stored relation matrix (relation_count x generator_count),
// under the preconditions shared by the HNF basis and the unit-row kernel.
bool context_relation_matrix(flint::FmpzMat& out,
                             const ClassGroupContext& context) noexcept {
    const Order* order = context.parent();
    if (!context.has_factor_base() || order == nullptr ||
        order->parent() == nullptr ||
        context.relation_count() < context.relation_rank() ||
        context.generator_count() <= 0 ||
        context.relation_rank() > context.generator_count()) {
        return false;
    }

    flint::FmpzMat stored(context.relation_count(),
                          context.generator_count());
    if (!context.relations(flint::FmpzMatRef(stored))) {
        return false;
    }
    out = std::move(stored);
    return true;
}

bool build_class_relation_hnf_transform(
        ClassRelationHnfTransform& out,
        const ClassGroupContext& context) noexcept {
    flint::FmpzMat stored(0, 0);
    if (!context_relation_matrix(stored, context)) {
        return false;
    }

    const slong relation_count = context.relation_count();
    const slong generator_count = context.generator_count();
    ClassRelationHnfTransform candidate;
    candidate.stored = std::move(stored);
    candidate.hnf = flint::FmpzMat(relation_count, generator_count);
    candidate.transform = flint::FmpzMat(relation_count, relation_count);
    ::fmpz_mat_hnf_transform(candidate.hnf.raw(), candidate.transform.raw(),
                             candidate.stored.raw());

    const slong rank =
            fmpz_mat_trimmed_nonzero_rows(
                    flint::FmpzMatConstRef(candidate.hnf));
    if (rank < 0 || rank != context.relation_rank()) {
        return false;
    }
    candidate.rank = rank;
    out = std::move(candidate);
    return true;
}

bool push_context_relation_witnesses(
        FactoredElement& out,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef coefficients,
        slong row) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*field) || row < 0 ||
        row >= flint::fmpz_mat_nrows(coefficients) ||
        flint::fmpz_mat_ncols(coefficients) != context.relation_count()) {
        return false;
    }

    FactoredElement candidate(*field);
    Element generator(*field);
    if (!candidate.is_defined() || !generator.is_defined() ||
        !candidate.one()) {
        return false;
    }
    for (slong col = 0; col < context.relation_count(); ++col) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(coefficients, row, col);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!context.relation_generator(generator, col) ||
            !multiply_element_power_fmpz(candidate, generator, exponent)) {
            return false;
        }
    }
    candidate.normalize();
    out.swap(candidate);
    return true;
}

// Proves (g_j) = prod_P P^{R_j[P]} for one stored relation, where g_j and R_j
// are the generator and row at index `relation` of the same context matrix
// `stored` that the witnesses are composed from.  Exact fractional-ideal
// equality settles support outside the factor base as well, so a generator
// whose factor-base valuations match but which carries another prime fails.
bool verify_stored_relation_generator(
        const ClassGroupContext& context,
        flint::FmpzMatConstRef stored,
        slong relation) noexcept {
    const Order* order = context.parent();
    const FactorBase* base = context.factor_base();
    if (order == nullptr || order->parent() == nullptr || base == nullptr ||
        relation < 0 || relation >= flint::fmpz_mat_nrows(stored) ||
        flint::fmpz_mat_ncols(stored) != base->length()) {
        return false;
    }

    Element generator(*order->parent());
    flint::FmpzMat row(1, base->length());
    FractionalIdeal principal(*order);
    FractionalIdeal expected(*order);
    return generator.is_defined() && principal.is_defined() &&
           expected.is_defined() &&
           context.relation_generator(generator, relation) &&
           fmpz_mat_copy_row(flint::FmpzMatRef(row), 0, stored, relation) &&
           principal.set_principal(generator, context.diagnostics()) &&
           factor_base_row_ideal(expected, *base,
                                 flint::FmpzMatConstRef(row)) &&
           principal.equal(expected);
}

// The unreduced n x (g + n) matrix [c * relations | I_n].
bool relation_kernel_extension_input(flint::FmpzMat& out,
                                     flint::FmpzMatConstRef relations)
        noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    const slong g = flint::fmpz_mat_ncols(relations);
    if (n <= 0 || g < 0) {
        return false;
    }

    flint::FmpzMat extension(n, g + n);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < g; ++j) {
            flint::fmpz_mul_2exp(
                    flint::fmpz_mat_entry(flint::FmpzMatRef(extension), i, j),
                    flint::fmpz_mat_entry(relations, i, j),
                    relation_kernel_lll_scale_bits);
        }
        flint::fmpz_one(flint::fmpz_mat_entry(flint::FmpzMatRef(extension), i,
                                              g + i));
    }
    out = std::move(extension);
    return true;
}

// FLINT's ULLL truncation threshold (FLINT 3.6.0 lll.c: fmpz_lll is
// fmpz_lll_with_removal_ulll(B, U, 250, NULL, fl)).  For a Z_BASIS input
// with |fmpz_mat_max_bits| at most this, lll_with_removal_ulll.c makes one
// fmpz_lll_wrapper_with_removal_knapsack(B, U, NULL, fl) call, which keeps
// the fmpz_lll_d_with_removal_knapsack output when that phase does not
// return -1 and fmpz_lll_is_reduced(B, fl, 120) accepts it
// (wrapper_with_removal_knapsack.c).
constexpr slong flint_ulll_truncation_bits = 250;

}  // namespace

bool lll_relation_kernel_extension(flint::FmpzMat& out,
                                   flint::FmpzMatConstRef relations) noexcept {
    flint::FmpzMat extension(0, 0);
    if (!relation_kernel_extension_input(extension, relations)) {
        return false;
    }
    // The input rows are independent (identity part), and fmpz_lll applies
    // only unimodular row operations.
    flint::FmpzLll lll;
    ::fmpz_lll(extension.raw(), nullptr, lll.raw());
    out = std::move(extension);
    return true;
}

bool checked_relation_kernel_rows(flint::FmpzMat& out,
                                  flint::FmpzMatConstRef relations,
                                  flint::FmpzMatConstRef reduced_extension,
                                  slong rank) noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    const slong g = flint::fmpz_mat_ncols(relations);
    if (rank < 0 || rank > n ||
        flint::fmpz_mat_nrows(reduced_extension) != n ||
        flint::fmpz_mat_ncols(reduced_extension) != g + n) {
        return false;
    }

    // Each reduced row is (c * v * relations | v) with v a row of the
    // unimodular transform, so any n - rank of its v with v * relations == 0
    // are a basis of the left kernel, whose rank is n - rank.
    const slong kernel_rows = n - rank;
    flint::FmpzMat rows(kernel_rows, n);
    slong found = 0;
    for (slong i = 0; i < n; ++i) {
        bool zero_relation_part = true;
        for (slong j = 0; j < g && zero_relation_part; ++j) {
            zero_relation_part = flint::fmpz_is_zero(
                    flint::fmpz_mat_entry(reduced_extension, i, j));
        }
        if (!zero_relation_part) {
            continue;
        }
        // Check 1: the row count is exactly n - rank.
        if (found == kernel_rows) {
            return false;
        }
        for (slong j = 0; j < n; ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(flint::FmpzMatRef(rows), found, j),
                    flint::fmpz_mat_entry(reduced_extension, i, g + j));
        }
        ++found;
    }
    if (found != kernel_rows) {
        return false;
    }

    // Check 2: every row satisfies v * relations == 0 exactly.
    if (kernel_rows > 0 && g > 0) {
        flint::FmpzMat product(kernel_rows, g);
        flint::fmpz_mat_mul(flint::FmpzMatRef(product),
                            flint::FmpzMatConstRef(rows), relations);
        if (::fmpz_mat_is_zero(product.raw()) == 0) {
            return false;
        }
    }
    out = std::move(rows);
    return true;
}

bool hnf_transform_kernel_rows(flint::FmpzMat& out,
                               flint::FmpzMatConstRef relations,
                               slong rank) noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    const slong g = flint::fmpz_mat_ncols(relations);
    if (rank < 0 || rank > n || rank > g) {
        return false;
    }

    flint::FmpzMat coefficients(n - rank, n);
    if (rank < n) {
        flint::FmpzMat hnf(n, g);
        flint::FmpzMat transform(n, n);
        ::fmpz_mat_hnf_transform(hnf.raw(), transform.raw(), relations.raw());
        if (fmpz_mat_trimmed_nonzero_rows(flint::FmpzMatConstRef(hnf)) !=
            rank) {
            return false;
        }
        for (slong row = rank; row < n; ++row) {
            if (!fmpz_mat_copy_row(flint::FmpzMatRef(coefficients),
                                   row - rank,
                                   flint::FmpzMatConstRef(transform), row)) {
                return false;
            }
        }
    }
    out = std::move(coefficients);
    return true;
}

bool relation_kernel_basis_from_extension(
        flint::FmpzMat& out,
        RelationKernelBasisRoute* route,
        flint::FmpzMatConstRef relations,
        slong rank,
        const flint::FmpzMat* reduced_extension) noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    if (rank < 0 || rank > n) {
        return false;
    }

    flint::FmpzMat checked(0, 0);
    if (reduced_extension != nullptr &&
        checked_relation_kernel_rows(checked, relations,
                                     flint::FmpzMatConstRef(*reduced_extension),
                                     rank)) {
        out = std::move(checked);
        if (route != nullptr) {
            *route = RelationKernelBasisRoute::lll_extended;
        }
        return true;
    }

    // Fallback: the transform is built only here, when the LLL rows are
    // absent or fail a check.
    flint::FmpzMat coefficients(0, 0);
    if (!hnf_transform_kernel_rows(coefficients, relations, rank)) {
        return false;
    }
    out = std::move(coefficients);
    if (route != nullptr) {
        *route = RelationKernelBasisRoute::hnf_transform;
    }
    return true;
}

bool relation_kernel_basis(flint::FmpzMat& out,
                           RelationKernelBasisRoute* route,
                           flint::FmpzMatConstRef relations,
                           slong rank) noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    flint::FmpzMat extension(0, 0);
    const bool have_extension =
            rank >= 0 && rank < n &&
            lll_relation_kernel_extension(extension, relations);
    return relation_kernel_basis_from_extension(
            out, route, relations, rank,
            have_extension ? &extension : nullptr);
}

bool finish_relation_kernel_basis(flint::FmpzMat& out,
                                  RelationKernelBasisRoute* route,
                                  flint::FmpzMatConstRef relations,
                                  slong rank) noexcept {
    const slong n = flint::fmpz_mat_nrows(relations);
    flint::FmpzMat extension(0, 0);
    if (rank < 0 || rank >= n ||
        !relation_kernel_extension_input(extension, relations)) {
        return relation_kernel_basis(out, route, relations, rank);
    }

    flint::FmpzLll lll;
    const slong bits = ::fmpz_mat_max_bits(extension.raw());
    if (bits >= -flint_ulll_truncation_bits &&
        bits <= flint_ulll_truncation_bits) {
        flint::FmpzMat reduced(n, flint::fmpz_mat_ncols(
                                          flint::FmpzMatConstRef(extension)));
        flint::fmpz_mat_set(flint::FmpzMatRef(reduced),
                            flint::FmpzMatConstRef(extension));
        // The first phase of fmpz_lll (see flint_ulll_truncation_bits).  It
        // changes the rows only by row moves and by adding integer
        // multiples of one row to another (lll_d.c, lll_check_babai.c), so
        // the identity part records a unimodular transform, which is all
        // that checked_relation_kernel_rows needs.  FLINT's final
        // reducedness check is not run: the kernel claim rests on the two
        // exact checks, not on the rows being LLL-reduced.
        int status = ::fmpz_lll_d_with_removal_knapsack(reduced.raw(), nullptr,
                                                        nullptr, lll.raw());
        const FinishKernelTestFailure failure =
                FinishKernelTestSeams::failure();
        if (failure != FinishKernelTestFailure::none) {
            // Test-only; see FinishKernelTestSeams.
            flint::fmpz_mat_zero(flint::FmpzMatRef(reduced));
            if (failure == FinishKernelTestFailure::first_phase) {
                status = -1;
            }
            FinishKernelTestSeams::note_injected();
        }
        flint::FmpzMat checked(0, 0);
        if (status != -1 &&
            checked_relation_kernel_rows(checked, relations,
                                         flint::FmpzMatConstRef(reduced),
                                         rank)) {
            out = std::move(checked);
            if (route != nullptr) {
                *route = RelationKernelBasisRoute::lll_double_extended;
            }
            return true;
        }
    }

    // Fallback: relation_kernel_basis on the unchanged input.  fmpz_lll
    // repeats the deterministic phase above before its own later stages.
    ::fmpz_lll(extension.raw(), nullptr, lll.raw());
    return relation_kernel_basis_from_extension(out, route, relations, rank,
                                                &extension);
}

namespace {

// The stored relations and their kernel rows through `basis`
// (relation_kernel_basis or finish_relation_kernel_basis).
template <class KernelBasis>
bool hnf_unit_coefficients(flint::FmpzMat& out,
                           const ClassGroupContext& context,
                           KernelBasis basis) noexcept {
    out = flint::FmpzMat(0, 0);
    flint::FmpzMat stored(0, 0);
    if (!context_relation_matrix(stored, context)) {
        return false;
    }

    // The rank check: the exact rank of the stored relations must be the
    // context's relation rank.  fmpz_mat_rank is the number of nonzero rows
    // of the HNF that the fallback builds, so this is the check that the
    // HNF transform made before, without the transform.
    const slong rank = context.relation_rank();
    if (::fmpz_mat_rank(stored.raw()) != rank) {
        return false;
    }

    flint::FmpzMat coefficients(0, 0);
    if (!basis(coefficients, nullptr, flint::FmpzMatConstRef(stored), rank)) {
        return false;
    }
    out = std::move(coefficients);
    return true;
}

}  // namespace

bool hnf_unit_witness_coefficients(
        flint::FmpzMat& out,
        const ClassGroupContext& context) noexcept {
    return hnf_unit_coefficients(out, context, relation_kernel_basis);
}

bool hnf_unit_finish_coefficients(flint::FmpzMat& out,
                                  const ClassGroupContext& context) noexcept {
    return hnf_unit_coefficients(out, context, finish_relation_kernel_basis);
}

bool hnf_unit_witnesses(std::vector<FactoredElement>& out,
                             const ClassGroupContext& context) noexcept {
    out.clear();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr) {
        return false;
    }

    flint::FmpzMat coefficients(0, 0);
    if (!hnf_unit_witness_coefficients(coefficients, context)) {
        return false;
    }

    const slong rows = flint::fmpz_mat_nrows(coefficients);
    out.reserve(static_cast<std::size_t>(rows));
    for (slong row = 0; row < rows; ++row) {
        FactoredElement candidate(*field);
        if (!candidate.is_defined() ||
            !push_context_relation_witnesses(
                    candidate, context,
                    flint::FmpzMatConstRef(coefficients), row)) {
            out.clear();
            return false;
        }
        out.push_back(std::move(candidate));
    }
    return true;
}

bool class_relation_witnessed_hnf_basis(
        WitnessedClassRelationHnfBasis& out,
        const ClassGroupContext& context) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    ClassRelationHnfTransform hnf_data;
    if (field == nullptr ||
        !build_class_relation_hnf_transform(hnf_data, context) ||
        hnf_data.rank != context.generator_count()) {
        return false;
    }

    const slong rank = hnf_data.rank;
    const slong relation_count = context.relation_count();
    const slong generator_count = context.generator_count();
    WitnessedClassRelationHnfBasis candidate;
    candidate.rows = flint::FmpzMat(rank, generator_count);
    candidate.relation_coefficients =
            flint::FmpzMat(rank, relation_count);
    for (slong row = 0; row < rank; ++row) {
        if (!fmpz_mat_copy_row(flint::FmpzMatRef(candidate.rows), row,
                               flint::FmpzMatConstRef(hnf_data.hnf), row) ||
            !fmpz_mat_copy_row(
                    flint::FmpzMatRef(candidate.relation_coefficients), row,
                    flint::FmpzMatConstRef(hnf_data.transform), row)) {
            return false;
        }
    }

    flint::FmpzMat product(rank, generator_count);
    flint::fmpz_mat_mul(
            flint::FmpzMatRef(product),
            flint::FmpzMatConstRef(candidate.relation_coefficients),
            flint::FmpzMatConstRef(hnf_data.stored));
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(product),
                               flint::FmpzMatConstRef(candidate.rows))) {
        return false;
    }

    // Each witness w_i = prod_j g_j^{C_ij} generates I(rows_i): the ideal map
    // I(v) = prod_P P^{v_P} is a homomorphism, so (g_j) = I(stored_j) for
    // every j with some C_ij != 0, together with the exact identity
    // C * stored == rows checked above, gives
    // (w_i) = prod_j I(stored_j)^{C_ij} = I(rows_i) without expanding w_i.
    for (slong j = 0; j < relation_count; ++j) {
        bool used = false;
        for (slong i = 0; i < rank && !used; ++i) {
            used = !flint::fmpz_is_zero(flint::fmpz_mat_entry(
                    flint::FmpzMatConstRef(candidate.relation_coefficients),
                    i, j));
        }
        if (used && !verify_stored_relation_generator(
                            context, flint::FmpzMatConstRef(hnf_data.stored),
                            j)) {
            return false;
        }
    }

    candidate.witnesses.reserve(static_cast<std::size_t>(rank));
    for (slong i = 0; i < rank; ++i) {
        FactoredElement witness(*field);
        if (!witness.is_defined() ||
            !push_context_relation_witnesses(
                    witness, context,
                    flint::FmpzMatConstRef(candidate.relation_coefficients),
                    i)) {
            return false;
        }
        candidate.witnesses.push_back(std::move(witness));
    }

    out = std::move(candidate);
    return true;
}

}  // namespace silex::detail
