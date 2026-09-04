#include <silex/abelian_group.hpp>

#include <cassert>
#include <utility>

namespace {
namespace sflint = silex::flint;

void set_entry_si(sflint::FmpzMat& matrix,
                  slong row,
                  slong col,
                  slong value) noexcept {
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(matrix, row, col), value);
}

bool entry_is_si(const sflint::FmpzMat& matrix,
                 slong row,
                 slong col,
                 slong value) noexcept {
    return sflint::fmpz_equal_si(sflint::fmpz_mat_entry(matrix, row, col),
                                 value);
}

bool entry_is_zero(const sflint::FmpzMat& matrix,
                   slong row,
                   slong col) noexcept {
    return sflint::fmpz_is_zero(sflint::fmpz_mat_entry(matrix, row, col));
}

sflint::FmpzConstRef entry_const(const sflint::FmpzMat& matrix,
                                 slong row,
                                 slong col) noexcept {
    return sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(matrix), row, col);
}

bool matrix_is_zero(const sflint::FmpzMat& matrix) noexcept {
    for (slong i = 0; i < sflint::fmpz_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < sflint::fmpz_mat_ncols(matrix); ++j) {
            if (!entry_is_zero(matrix, i, j)) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_is_identity(const sflint::FmpzMat& matrix) noexcept {
    if (sflint::fmpz_mat_nrows(matrix) != sflint::fmpz_mat_ncols(matrix)) {
        return false;
    }
    for (slong i = 0; i < sflint::fmpz_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < sflint::fmpz_mat_ncols(matrix); ++j) {
            const slong expected = i == j ? 1 : 0;
            if (!entry_is_si(matrix, i, j, expected)) {
                return false;
            }
        }
    }
    return true;
}

bool invariants_are_si(const silex::FiniteAbelianGroup& group,
                       const slong* values,
                       slong length) noexcept {
    if (group.invariant_count() != length) {
        return false;
    }
    sflint::Fmpz value;
    for (slong i = 0; i < length; ++i) {
        if (!group.invariant(sflint::FmpzRef(value), i) ||
            !sflint::fmpz_equal_si(value, values[i])) {
            return false;
        }
    }
    return true;
}

bool generator_relation_identity(const silex::FiniteAbelianGroup& group) noexcept {
    const slong rank = group.invariant_count();
    const slong m = group.relation_count();
    const slong n = group.generator_count();
    sflint::FmpzMat relations(m, n);
    sflint::FmpzMat generators(rank, n);
    sflint::FmpzMat combinations(rank, m);
    sflint::FmpzMat left(rank, n);
    sflint::Fmpz inv;

    if (!group.relations(sflint::FmpzMatRef(relations)) ||
        !group.invariant_generator_matrix(sflint::FmpzMatRef(generators)) ||
        !group.invariant_generator_relation_matrix(
                sflint::FmpzMatRef(combinations))) {
        return false;
    }
    auto owned_relations = group.relations();
    auto owned_generators = group.invariant_generator_matrix();
    if (!owned_relations.has_value() || !owned_generators.has_value() ||
        !sflint::fmpz_mat_equal(*owned_relations, relations) ||
        !sflint::fmpz_mat_equal(*owned_generators, generators)) {
        return false;
    }

    sflint::fmpz_mat_mul(sflint::FmpzMatRef(left),
                         sflint::FmpzMatConstRef(combinations),
                         sflint::FmpzMatConstRef(relations));
    for (slong i = 0; i < rank; ++i) {
        if (!group.invariant(sflint::FmpzRef(inv), i)) {
            return false;
        }
        auto owned_invariant = group.invariant(i);
        if (!owned_invariant.has_value() ||
            !sflint::fmpz_equal(*owned_invariant, inv)) {
            return false;
        }
        for (slong j = 0; j < n; ++j) {
            sflint::Fmpz expected;
            sflint::fmpz_mul(sflint::FmpzRef(expected),
                             entry_const(generators, i, j),
                             sflint::FmpzConstRef(inv));
            if (!sflint::fmpz_equal(entry_const(left, i, j),
                                    sflint::FmpzConstRef(expected))) {
                return false;
            }
        }
    }
    return true;
}

bool relation_kernel_identity(const silex::FiniteAbelianGroup& group) noexcept {
    const slong kernels = group.relation_kernel_count();
    const slong m = group.relation_count();
    const slong n = group.generator_count();
    sflint::FmpzMat relations(m, n);
    sflint::FmpzMat kernel(kernels, m);
    sflint::FmpzMat kernel_row(1, m);
    sflint::FmpzMat product(kernels, n);

    if (!group.relations(sflint::FmpzMatRef(relations)) ||
        !group.relation_kernel_matrix(sflint::FmpzMatRef(kernel))) {
        return false;
    }
    auto owned_kernel = group.relation_kernel_matrix();
    if (!owned_kernel.has_value() ||
        !sflint::fmpz_mat_equal(*owned_kernel, kernel)) {
        return false;
    }
    for (slong i = 0; i < kernels; ++i) {
        if (!group.relation_kernel_row(sflint::FmpzMatRef(kernel_row), i)) {
            return false;
        }
        for (slong j = 0; j < m; ++j) {
            if (!sflint::fmpz_equal(
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(kernel_row), 0, j),
                        sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(kernel), i, j))) {
                return false;
            }
        }
    }

    sflint::fmpz_mat_mul(sflint::FmpzMatRef(product),
                         sflint::FmpzMatConstRef(kernel),
                         sflint::FmpzMatConstRef(relations));
    return matrix_is_zero(product);
}

int test_basic() {
    silex::FiniteAbelianGroup group;
    sflint::FmpzMat zero_relations(0, 0);
    sflint::FmpzMat relations(2, 2);
    sflint::FmpzMat generators(0, 0);
    sflint::FmpzMat row(1, 2);
    sflint::Fmpz order;
    const slong invs[2] = {2, 4};

    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(zero_relations)));
    assert(group.is_defined());
    assert(group.relation_count() == 0);
    assert(group.generator_count() == 0);
    assert(group.invariant_count() == 0);
    assert(group.relation_kernel_count() == 0);
    assert(invariants_are_si(group, nullptr, 0));
    assert(group.invariant_generator_matrix(sflint::FmpzMatRef(generators)));
    assert(sflint::fmpz_mat_nrows(generators) == 0);
    assert(sflint::fmpz_mat_ncols(generators) == 0);
    assert(group.relation_kernel_matrix(sflint::FmpzMatRef(generators)));
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_is_one(order));

    set_entry_si(relations, 0, 0, 2);
    set_entry_si(relations, 1, 1, 4);
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(group.invariant_count() == 2);
    assert(invariants_are_si(group, invs, 2));
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 8));

    sflint::FmpzMat invariant_generators(2, 2);
    assert(group.invariant_generator_matrix(
            sflint::FmpzMatRef(invariant_generators)));
    assert(matrix_is_identity(invariant_generators));

    set_entry_si(row, 0, 0, 2);
    assert(group.reduce(sflint::FmpzMatRef(row)));
    assert(matrix_is_zero(row));
    set_entry_si(row, 0, 1, 4);
    assert(group.reduce(sflint::FmpzMatRef(row)));
    assert(matrix_is_zero(row));

    return 0;
}

int test_copy_swap_failure_preserves_output() {
    silex::FiniteAbelianGroup group;
    silex::FiniteAbelianGroup copy;
    silex::FiniteAbelianGroup other;
    sflint::FmpzMat relations(1, 1);
    sflint::FmpzMat other_relations(1, 1);
    sflint::FmpzMat bad(1, 2);
    sflint::FmpzMat row(1, 2);
    sflint::Fmpz order;

    set_entry_si(relations, 0, 0, 5);
    set_entry_si(other_relations, 0, 0, 7);
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(other.set_relation_matrix(sflint::FmpzMatConstRef(other_relations)));

    assert(group.set(group));
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 5));

    assert(copy.set(group));
    assert(copy.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 5));

    swap(group, other);
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 7));
    assert(other.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 5));

    set_entry_si(bad, 0, 0, 2);
    assert(!other.set_relation_matrix(sflint::FmpzMatConstRef(bad)));
    assert(other.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 5));

    set_entry_si(row, 0, 0, 11);
    set_entry_si(row, 0, 1, 13);
    assert(!other.reduce(sflint::FmpzMatRef(row)));
    assert(entry_is_si(row, 0, 0, 11));
    assert(entry_is_si(row, 0, 1, 13));

    return 0;
}

int test_element_reduce_and_coordinates() {
    silex::FiniteAbelianGroup group;
    silex::FiniteAbelianGroup unset;
    sflint::FmpzMat relations(3, 2);
    sflint::FmpzMat row(1, 2);
    sflint::FmpzMat coords(1, 2);
    sflint::FmpzMat bad(1, 1);

    set_entry_si(relations, 0, 0, 2);
    set_entry_si(relations, 1, 1, 4);
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));

    set_entry_si(row, 0, 0, 5);
    set_entry_si(row, 0, 1, 7);
    assert(group.reduce(sflint::FmpzMatRef(row)));
    assert(entry_is_si(row, 0, 0, 1));
    assert(entry_is_si(row, 0, 1, 3));
    assert(group.reduce(sflint::FmpzMatRef(row)));
    assert(entry_is_si(row, 0, 0, 1));
    assert(entry_is_si(row, 0, 1, 3));

    set_entry_si(row, 0, 0, 5);
    set_entry_si(row, 0, 1, 7);
    assert(group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                       sflint::FmpzMatConstRef(row)));
    assert(entry_is_si(coords, 0, 0, 1));
    assert(entry_is_si(coords, 0, 1, 3));

    set_entry_si(row, 0, 0, -1);
    set_entry_si(row, 0, 1, -1);
    assert(group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                       sflint::FmpzMatConstRef(row)));
    assert(entry_is_si(coords, 0, 0, 1));
    assert(entry_is_si(coords, 0, 1, 3));

    set_entry_si(coords, 0, 0, 88);
    set_entry_si(coords, 0, 1, 89);
    assert(!unset.invariant_coordinates(sflint::FmpzMatRef(coords),
                                        sflint::FmpzMatConstRef(row)));
    assert(!group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                        sflint::FmpzMatConstRef(bad)));
    assert(entry_is_si(coords, 0, 0, 88));
    assert(entry_is_si(coords, 0, 1, 89));

    return 0;
}

int test_nondiagonal_and_generator_relations() {
    silex::FiniteAbelianGroup group;
    sflint::FmpzMat relations(3, 2);
    sflint::FmpzMat generators(1, 2);
    sflint::FmpzMat row(1, 2);
    sflint::FmpzMat coords(1, 1);
    sflint::Fmpz order;
    const slong inv_6[1] = {6};
    const slong inv_2_6[2] = {2, 6};

    set_entry_si(relations, 0, 0, 2);
    set_entry_si(relations, 1, 1, 3);
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(group.invariant_count() == 1);
    assert(invariants_are_si(group, inv_6, 1));
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 6));
    assert(generator_relation_identity(group));
    assert(group.relation_kernel_count() == 1);
    assert(relation_kernel_identity(group));

    assert(group.invariant_generator_matrix(sflint::FmpzMatRef(generators)));
    sflint::fmpz_set(sflint::fmpz_mat_entry(row, 0, 0),
                     entry_const(generators, 0, 0));
    sflint::fmpz_set(sflint::fmpz_mat_entry(row, 0, 1),
                     entry_const(generators, 0, 1));
    assert(group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                       sflint::FmpzMatConstRef(row)));
    assert(entry_is_si(coords, 0, 0, 1));

    for (slong j = 0; j < 2; ++j) {
        sflint::Fmpz tmp;
        sflint::fmpz_mul(sflint::FmpzRef(tmp),
                         entry_const(row, 0, j),
                         sflint::FmpzConstRef(order));
        sflint::fmpz_set(sflint::fmpz_mat_entry(row, 0, j),
                         sflint::FmpzConstRef(tmp));
    }
    assert(group.reduce(sflint::FmpzMatRef(row)));
    assert(matrix_is_zero(row));

    sflint::fmpz_mat_zero(sflint::FmpzMatRef(relations));
    set_entry_si(relations, 0, 0, 4);
    set_entry_si(relations, 1, 1, 6);
    set_entry_si(relations, 2, 0, 2);
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(group.invariant_count() == 2);
    assert(invariants_are_si(group, inv_2_6, 2));
    assert(group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 12));
    assert(generator_relation_identity(group));
    assert(group.relation_kernel_count() == 1);
    assert(relation_kernel_identity(group));

    sflint::FmpzMat hnf_basis(2, 2);
    silex::FiniteAbelianGroup hnf_group;
    set_entry_si(hnf_basis, 0, 0, 2);
    set_entry_si(hnf_basis, 1, 1, 6);
    assert(hnf_group.set_relation_matrix_with_hnf_basis(
            sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(hnf_basis)));
    assert(hnf_group.relation_count() == 3);
    assert(hnf_group.generator_count() == 2);
    assert(hnf_group.invariant_count() == 2);
    assert(invariants_are_si(hnf_group, inv_2_6, 2));
    assert(hnf_group.order(sflint::FmpzRef(order)));
    assert(sflint::fmpz_equal_si(order, 12));
    assert(generator_relation_identity(hnf_group));
    assert(hnf_group.relation_kernel_count() == 1);
    assert(relation_kernel_identity(hnf_group));

    return 0;
}

int test_trivial_snf_coordinates() {
    silex::FiniteAbelianGroup group;
    sflint::FmpzMat relations(0, 0);
    sflint::FmpzMat row(1, 0);
    sflint::FmpzMat coords(1, 0);

    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                       sflint::FmpzMatConstRef(row)));
    assert(group.reduce(sflint::FmpzMatRef(row)));

    return 0;
}

}  // namespace

int main() {
    assert(test_basic() == 0);
    assert(test_copy_swap_failure_preserves_output() == 0);
    assert(test_element_reduce_and_coordinates() == 0);
    assert(test_nondiagonal_and_generator_relations() == 0);
    assert(test_trivial_snf_coordinates() == 0);
    return 0;
}
