#include <silex/abelian_group.hpp>

#include "abelian_group/abelian_group_internal.hpp"

#include <flint/fmpz_mat.h>

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

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

// True when the rows of `kernel` form a Z-basis of the integer left kernel of
// `relations`: they annihilate it, there are exactly m - rank of them and they
// are independent, and the row lattice is saturated (all Smith invariants of
// `kernel` are 1, equivalently the k x k minors have gcd 1).  A row set that
// spans only a proper sublattice of the kernel, for example a basis with one
// row doubled (|det U| = 2 in the unimodular completion), fails the last test.
bool is_left_kernel_basis(const sflint::FmpzMat& relations,
                          const sflint::FmpzMat& kernel) noexcept {
    const slong m = sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(relations));
    const slong n = sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(relations));
    const slong k = sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(kernel));
    if (sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(kernel)) != m) {
        return false;
    }
    sflint::FmpzMat product(k, n);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(product),
                         sflint::FmpzMatConstRef(kernel),
                         sflint::FmpzMatConstRef(relations));
    if (!matrix_is_zero(product)) {
        return false;
    }
    if (k != m - sflint::fmpz_mat_rank(sflint::FmpzMatConstRef(relations)) ||
        sflint::fmpz_mat_rank(sflint::FmpzMatConstRef(kernel)) != k) {
        return false;
    }
    if (k == 0) {
        return true;
    }
    sflint::FmpzMat snf(k, m);
    ::fmpz_mat_snf(snf.raw(), kernel.raw());
    for (slong i = 0; i < k; ++i) {
        if (!entry_is_si(snf, i, i, 1)) {
            return false;
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
    return matrix_is_zero(product) && is_left_kernel_basis(relations, kernel);
}

// The basis check must reject a kernel row set that is only a sublattice of
// the left kernel, and one that is not the full kernel.
int test_kernel_basis_check_rejects_bad_kernels() {
    // relations = [[2, 0], [0, 3], [2, 3]]: left kernel is spanned by (1,1,-1).
    sflint::FmpzMat relations(3, 2);
    set_entry_si(relations, 0, 0, 2);
    set_entry_si(relations, 1, 1, 3);
    set_entry_si(relations, 2, 0, 2);
    set_entry_si(relations, 2, 1, 3);
    silex::FiniteAbelianGroup group;
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    auto kernel = group.relation_kernel_matrix();
    assert(kernel.has_value());
    assert(is_left_kernel_basis(relations, *kernel));

    // Rebuild the one kernel row doubled: (2, 2, -2).
    sflint::FmpzMat doubled(1, 3);
    set_entry_si(doubled, 0, 0, 2);
    set_entry_si(doubled, 0, 1, 2);
    set_entry_si(doubled, 0, 2, -2);
    assert(sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(*kernel)) == 3);
    assert(!is_left_kernel_basis(relations, doubled));

    sflint::FmpzMat empty(0, 3);
    assert(!is_left_kernel_basis(relations, empty));
    return 0;
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
    assert(silex::detail::FiniteAbelianGroupAccess::
                   set_relation_matrix_with_hnf_basis(
            hnf_group,
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

bool generator_state_equal(const silex::FiniteAbelianGroup& group,
                           const sflint::FmpzMat& generators,
                           const sflint::FmpzMat& coordinates) noexcept {
    const slong n = group.generator_count();
    auto current = group.invariant_generator_matrix();
    if (!current.has_value() ||
        !sflint::fmpz_mat_equal(*current, generators)) {
        return false;
    }
    sflint::FmpzMat unit_row(1, n);
    sflint::FmpzMat coords(1, group.invariant_count());
    for (slong j = 0; j < n; ++j) {
        sflint::fmpz_mat_zero(sflint::FmpzMatRef(unit_row));
        set_entry_si(unit_row, 0, j, 1);
        if (!group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                         sflint::FmpzMatConstRef(unit_row))) {
            return false;
        }
        for (slong k = 0; k < group.invariant_count(); ++k) {
            if (!sflint::fmpz_equal(entry_const(coords, 0, k),
                                    entry_const(coordinates, j, k))) {
                return false;
            }
        }
    }
    return true;
}

bool snapshot_generator_state(const silex::FiniteAbelianGroup& group,
                              sflint::FmpzMat& generators,
                              sflint::FmpzMat& coordinates) noexcept {
    const slong n = group.generator_count();
    sflint::FmpzMat unit_row(1, n);
    sflint::FmpzMat coords(1, group.invariant_count());
    if (!group.invariant_generator_matrix(sflint::FmpzMatRef(generators))) {
        return false;
    }
    for (slong j = 0; j < n; ++j) {
        sflint::fmpz_mat_zero(sflint::FmpzMatRef(unit_row));
        set_entry_si(unit_row, 0, j, 1);
        if (!group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                         sflint::FmpzMatConstRef(unit_row))) {
            return false;
        }
        for (slong k = 0; k < group.invariant_count(); ++k) {
            sflint::fmpz_set(sflint::fmpz_mat_entry(coordinates, j, k),
                             entry_const(coords, 0, k));
        }
    }
    return true;
}

// Invariant generators and coordinates are published from the Smith form of
// the stored HNF basis.  The first witness request takes the HNF of the full
// relation matrix, checks it against that basis, and reruns the Smith form on
// the same basis, so the published state does not move.  Pin that the
// witness request leaves generators and coordinates unchanged on a larger
// group.
int test_witness_request_keeps_generators() {
    constexpr slong n = 6;
    constexpr slong m = 14;
    const slong diagonal[n] = {1, 2, 6, 30, 210, 2310};
    const slong inv[5] = {2, 6, 30, 210, 2310};

    // relations = M * D * W with W unimodular and M containing a unimodular
    // top block, so the group is Z/2 + Z/6 + Z/30 + Z/210 + Z/2310.
    sflint::FmpzMat mix(m, n);
    sflint::FmpzMat scaled(n, n);
    sflint::FmpzMat right(n, n);
    sflint::FmpzMat partial(m, n);
    sflint::FmpzMat relations(m, n);
    for (slong i = 0; i < m; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (i < n) {
                set_entry_si(mix, i, j,
                             j == i ? 1 : (j < i ? (3 * i + j) % 5 - 2 : 0));
            } else {
                set_entry_si(mix, i, j, (7 * i + 3 * j) % 9 - 4);
            }
        }
    }
    for (slong i = 0; i < n; ++i) {
        set_entry_si(scaled, i, i, diagonal[i]);
        for (slong j = 0; j < n; ++j) {
            set_entry_si(right, i, j,
                         j == i ? 1 : (j > i ? (i + 2 * j) % 7 - 3 : 0));
        }
    }
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(partial),
                         sflint::FmpzMatConstRef(mix),
                         sflint::FmpzMatConstRef(scaled));
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(relations),
                         sflint::FmpzMatConstRef(partial),
                         sflint::FmpzMatConstRef(right));

    sflint::FmpzMat hnf_full(m, n);
    sflint::FmpzMat hnf_basis(n, n);
    ::fmpz_mat_hnf(hnf_full.raw(), relations.raw());
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            sflint::fmpz_set(sflint::fmpz_mat_entry(hnf_basis, i, j),
                             entry_const(hnf_full, i, j));
        }
    }
    for (slong i = n; i < m; ++i) {
        for (slong j = 0; j < n; ++j) {
            assert(entry_is_zero(hnf_full, i, j));
        }
    }

    // The public setter forwards the same HNF top block to the internal
    // setter; the second iteration pins the internal path directly.
    for (int with_basis = 0; with_basis < 2; ++with_basis) {
        silex::FiniteAbelianGroup group;
        if (with_basis != 0) {
            assert(silex::detail::FiniteAbelianGroupAccess::
                           set_relation_matrix_with_hnf_basis(
                    group,
                    sflint::FmpzMatConstRef(relations),
                    sflint::FmpzMatConstRef(hnf_basis)));
        } else {
            assert(group.set_relation_matrix(
                    sflint::FmpzMatConstRef(relations)));
        }
        assert(invariants_are_si(group, inv, 5));

        sflint::FmpzMat generators(5, n);
        sflint::FmpzMat coordinates(n, 5);
        assert(snapshot_generator_state(group, generators, coordinates));

        assert(generator_relation_identity(group));
        assert(relation_kernel_identity(group));
        assert(invariants_are_si(group, inv, 5));
        assert(generator_state_equal(group, generators, coordinates));

        assert(generator_relation_identity(group));
        assert(generator_state_equal(group, generators, coordinates));
    }

    return 0;
}

// FLINT's fmpz_mat_snf_transform skips its row-HNF phase when the input is
// already diagonal, so a signed diagonal relation matrix and its HNF basis
// reach the gcd/lcm phase with different signs.  Pin that witness requests on
// the public setter leave generators and coordinates unchanged there too.
bool diagonal_witness_keeps_generators(const slong* diagonal,
                                       slong n,
                                       slong zero_rows,
                                       const slong* inv,
                                       slong inv_count) noexcept {
    sflint::FmpzMat relations(n + zero_rows, n);
    for (slong i = 0; i < n; ++i) {
        set_entry_si(relations, i, i, diagonal[i]);
    }
    silex::FiniteAbelianGroup group;
    if (!group.set_relation_matrix(sflint::FmpzMatConstRef(relations)) ||
        !invariants_are_si(group, inv, inv_count)) {
        return false;
    }

    sflint::FmpzMat generators(inv_count, n);
    sflint::FmpzMat coordinates(n, inv_count);
    return snapshot_generator_state(group, generators, coordinates) &&
           generator_relation_identity(group) &&
           relation_kernel_identity(group) &&
           invariants_are_si(group, inv, inv_count) &&
           generator_state_equal(group, generators, coordinates) &&
           generator_relation_identity(group) &&
           generator_state_equal(group, generators, coordinates);
}

int test_diagonal_witness_keeps_generators() {
    const slong neg_2_3[2] = {-2, 3};
    const slong pos_2_neg_3[2] = {2, -3};
    const slong mixed[3] = {4, -6, 10};
    const slong inv_6[1] = {6};
    const slong inv_2_2_60[3] = {2, 2, 60};

    assert(diagonal_witness_keeps_generators(neg_2_3, 2, 0, inv_6, 1));
    assert(diagonal_witness_keeps_generators(neg_2_3, 2, 1, inv_6, 1));
    assert(diagonal_witness_keeps_generators(pos_2_neg_3, 2, 0, inv_6, 1));
    assert(diagonal_witness_keeps_generators(mixed, 3, 0, inv_2_2_60, 3));
    assert(diagonal_witness_keeps_generators(mixed, 3, 2, inv_2_2_60, 3));
    return 0;
}

// The internal basis setter brings a spanning basis that is not in HNF to
// HNF before its Smith form, so witness requests see the same basis.  A basis
// that does not span the relation lattice makes witness requests fail closed
// without changing the published generators.
int test_internal_basis_setter_witnesses() {
    using silex::detail::FiniteAbelianGroupAccess;
    const slong inv_2_6[2] = {2, 6};
    sflint::FmpzMat relations(3, 2);
    set_entry_si(relations, 0, 0, 4);
    set_entry_si(relations, 1, 1, 6);
    set_entry_si(relations, 2, 0, 2);

    silex::FiniteAbelianGroup reference;
    assert(reference.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    sflint::FmpzMat reference_generators(2, 2);
    sflint::FmpzMat reference_coordinates(2, 2);
    assert(snapshot_generator_state(reference, reference_generators,
                                    reference_coordinates));

    // A unimodular row mix of the HNF basis diag(2, 6).
    sflint::FmpzMat mixed_basis(2, 2);
    set_entry_si(mixed_basis, 0, 0, 2);
    set_entry_si(mixed_basis, 0, 1, 6);
    set_entry_si(mixed_basis, 1, 1, -6);
    silex::FiniteAbelianGroup mixed;
    assert(FiniteAbelianGroupAccess::set_relation_matrix_with_hnf_basis(
            mixed, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(mixed_basis)));
    assert(invariants_are_si(mixed, inv_2_6, 2));
    assert(generator_state_equal(mixed, reference_generators,
                                 reference_coordinates));
    assert(generator_relation_identity(mixed));
    assert(relation_kernel_identity(mixed));
    assert(generator_state_equal(mixed, reference_generators,
                                 reference_coordinates));

    // diag(4, 6) spans an index-2 sublattice of the relation lattice.
    sflint::FmpzMat short_basis(2, 2);
    set_entry_si(short_basis, 0, 0, 4);
    set_entry_si(short_basis, 1, 1, 6);
    silex::FiniteAbelianGroup misuse;
    assert(FiniteAbelianGroupAccess::set_relation_matrix_with_hnf_basis(
            misuse, sflint::FmpzMatConstRef(relations),
            sflint::FmpzMatConstRef(short_basis)));
    const slong count = misuse.invariant_count();
    sflint::FmpzMat generators(count, 2);
    sflint::FmpzMat coordinates(2, count);
    sflint::FmpzMat combinations(count, 3);
    sflint::FmpzMat kernel(misuse.relation_kernel_count(), 3);
    assert(snapshot_generator_state(misuse, generators, coordinates));
    assert(!misuse.invariant_generator_relation_matrix(
            sflint::FmpzMatRef(combinations)));
    assert(!misuse.relation_kernel_matrix(sflint::FmpzMatRef(kernel)));
    assert(misuse.is_defined());
    assert(generator_state_equal(misuse, generators, coordinates));
    return 0;
}

}  // namespace

// invariant_generator_matrix and reduce read the inverse Smith right
// transform, which is filled on first use.  Copies and moves taken before or
// after that first use, and a republished group, must all give the values a
// fresh group gives.
bool generator_matrix_and_reduce_match(const silex::FiniteAbelianGroup& left,
                                       const silex::FiniteAbelianGroup& right) {
    const slong k = left.invariant_count();
    const slong n = left.generator_count();
    if (right.invariant_count() != k || right.generator_count() != n) {
        return false;
    }
    sflint::FmpzMat left_generators(k, n);
    sflint::FmpzMat right_generators(k, n);
    if (!left.invariant_generator_matrix(sflint::FmpzMatRef(left_generators)) ||
        !right.invariant_generator_matrix(
                sflint::FmpzMatRef(right_generators)) ||
        !sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left_generators),
                                sflint::FmpzMatConstRef(right_generators))) {
        return false;
    }
    sflint::FmpzMat left_row(1, n);
    sflint::FmpzMat right_row(1, n);
    for (slong j = 0; j < n; ++j) {
        set_entry_si(left_row, 0, j, 7 * j + 5);
        set_entry_si(right_row, 0, j, 7 * j + 5);
    }
    return left.reduce(sflint::FmpzMatRef(left_row)) &&
           right.reduce(sflint::FmpzMatRef(right_row)) &&
           sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left_row),
                                  sflint::FmpzMatConstRef(right_row));
}

int test_lazy_right_transform_inverse() {
    // Z^2 modulo rows (2, 1), (0, 3), (4, 5): cyclic of order 6, with a
    // nontrivial Smith right transform.
    sflint::FmpzMat relations(3, 2);
    set_entry_si(relations, 0, 0, 2);
    set_entry_si(relations, 0, 1, 1);
    set_entry_si(relations, 1, 1, 3);
    set_entry_si(relations, 2, 0, 4);
    set_entry_si(relations, 2, 1, 5);
    const slong inv_6[1] = {6};

    silex::FiniteAbelianGroup reference;
    silex::FiniteAbelianGroup group;
    silex::FiniteAbelianGroup moved_source;
    assert(reference.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(relations)));
    assert(moved_source.set_relation_matrix(
            sflint::FmpzMatConstRef(relations)));
    assert(invariants_are_si(group, inv_6, 1));

    // Copy and move before any first use.
    silex::FiniteAbelianGroup early_copy;
    assert(early_copy.set(group));
    silex::FiniteAbelianGroup moved(std::move(moved_source));

    // The generator has invariant coordinate 1 and order 6.
    sflint::FmpzMat generators(1, 2);
    sflint::FmpzMat coords(1, 1);
    assert(group.invariant_generator_matrix(sflint::FmpzMatRef(generators)));
    assert(group.invariant_coordinates(sflint::FmpzMatRef(coords),
                                       sflint::FmpzMatConstRef(generators)));
    assert(entry_is_si(coords, 0, 0, 1));
    sflint::FmpzMat row(1, 2);
    sflint::Fmpz six;
    sflint::fmpz_set_si(sflint::FmpzRef(six), 6);
    for (slong j = 0; j < 2; ++j) {
        sflint::fmpz_mul(sflint::fmpz_mat_entry(row, 0, j),
                         entry_const(generators, 0, j),
                         sflint::FmpzConstRef(six));
    }
    assert(group.reduce(sflint::FmpzMatRef(row)));
    assert(matrix_is_zero(row));

    // Copy after first use.
    silex::FiniteAbelianGroup late_copy;
    assert(late_copy.set(group));

    assert(generator_matrix_and_reduce_match(reference, group));
    assert(generator_matrix_and_reduce_match(reference, early_copy));
    assert(generator_matrix_and_reduce_match(reference, moved));
    assert(generator_matrix_and_reduce_match(reference, late_copy));

    // Republishing replaces the cached inverse.
    sflint::FmpzMat other_relations(2, 2);
    set_entry_si(other_relations, 0, 0, 4);
    set_entry_si(other_relations, 0, 1, 2);
    set_entry_si(other_relations, 1, 1, 6);
    silex::FiniteAbelianGroup other_reference;
    assert(other_reference.set_relation_matrix(
            sflint::FmpzMatConstRef(other_relations)));
    assert(group.set_relation_matrix(sflint::FmpzMatConstRef(other_relations)));
    assert(generator_matrix_and_reduce_match(other_reference, group));
    assert(late_copy.set(other_reference));
    assert(generator_matrix_and_reduce_match(other_reference, late_copy));

    // Move assignment carries the cache state with the object.
    moved = std::move(late_copy);
    assert(generator_matrix_and_reduce_match(other_reference, moved));
    group.clear();
    sflint::FmpzMat unused(0, 0);
    assert(!group.invariant_generator_matrix(sflint::FmpzMatRef(unused)));

    return 0;
}

// A cache-taking publish gives the group the plain publish gives, on the
// same relations and basis: same invariants, generators, reduction and
// relation-kernel witnesses (the witness path reruns the Smith form on the
// stored basis and compares it with the published one).
bool cached_publish_matches_plain(
        silex::detail::FiniteAbelianGroupSmithCache& cache,
        const sflint::FmpzMat& relations,
        const sflint::FmpzMat& basis) {
    using silex::detail::FiniteAbelianGroupAccess;
    silex::FiniteAbelianGroup cached;
    silex::FiniteAbelianGroup plain;
    if (!FiniteAbelianGroupAccess::set_relation_matrix_with_hnf_basis(
                cached, sflint::FmpzMatConstRef(relations),
                sflint::FmpzMatConstRef(basis), cache) ||
        !FiniteAbelianGroupAccess::set_relation_matrix_with_hnf_basis(
                plain, sflint::FmpzMatConstRef(relations),
                sflint::FmpzMatConstRef(basis))) {
        return false;
    }
    const slong k = plain.invariant_count();
    if (cached.invariant_count() != k ||
        cached.relation_kernel_count() != plain.relation_kernel_count()) {
        return false;
    }
    for (slong i = 0; i < k; ++i) {
        const auto cached_invariant = cached.invariant(i);
        const auto plain_invariant = plain.invariant(i);
        if (!cached_invariant.has_value() || !plain_invariant.has_value() ||
            !sflint::fmpz_equal(*cached_invariant, *plain_invariant)) {
            return false;
        }
    }
    const slong m = plain.relation_count();
    sflint::FmpzMat cached_kernel(plain.relation_kernel_count(), m);
    sflint::FmpzMat plain_kernel(plain.relation_kernel_count(), m);
    return generator_matrix_and_reduce_match(cached, plain) &&
           cached.relation_kernel_matrix(sflint::FmpzMatRef(cached_kernel)) &&
           plain.relation_kernel_matrix(sflint::FmpzMatRef(plain_kernel)) &&
           sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(cached_kernel),
                                  sflint::FmpzMatConstRef(plain_kernel)) &&
           relation_kernel_identity(cached);
}

// The Smith form cache is keyed on the stored HNF basis alone: a publish with
// more relations but the same basis (given in HNF or not) reuses it, and a
// publish whose basis changed recomputes and replaces it.
int test_smith_cache_reuses_unchanged_basis() {
    silex::detail::FiniteAbelianGroupSmithCache cache;

    // Rows of a lattice with HNF basis [[2, 0, 2], [0, 3, 3], [0, 0, 6]].
    sflint::FmpzMat relations(4, 3);
    const slong rows[4][3] = {{2, 3, 5}, {0, 3, 3}, {0, 0, 6}, {2, 0, 2}};
    for (slong i = 0; i < 4; ++i) {
        for (slong j = 0; j < 3; ++j) {
            set_entry_si(relations, i, j, rows[i][j]);
        }
    }
    sflint::FmpzMat basis(3, 3);
    set_entry_si(basis, 0, 0, 2);
    set_entry_si(basis, 0, 2, 2);
    set_entry_si(basis, 1, 1, 3);
    set_entry_si(basis, 1, 2, 3);
    set_entry_si(basis, 2, 2, 6);
    assert(cached_publish_matches_plain(cache, relations, basis));
    assert(cache.valid && cache.compute_count == 1 && cache.reuse_count == 0);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(cache.hnf_basis),
                                  sflint::FmpzMatConstRef(basis)));

    // One more relation in the same lattice, same basis: reused.
    sflint::FmpzMat more(5, 3);
    for (slong i = 0; i < 4; ++i) {
        for (slong j = 0; j < 3; ++j) {
            set_entry_si(more, i, j, rows[i][j]);
        }
    }
    set_entry_si(more, 4, 0, 4);
    set_entry_si(more, 4, 1, 6);
    set_entry_si(more, 4, 2, 4);
    assert(cached_publish_matches_plain(cache, more, basis));
    assert(cache.compute_count == 1 && cache.reuse_count == 1);

    // A unimodular row mix of the same basis is brought to HNF first, so the
    // key still matches.
    sflint::FmpzMat mixed(3, 3);
    set_entry_si(mixed, 0, 0, 2);
    set_entry_si(mixed, 0, 1, 3);
    set_entry_si(mixed, 0, 2, 5);
    set_entry_si(mixed, 1, 1, 3);
    set_entry_si(mixed, 1, 2, 3);
    set_entry_si(mixed, 2, 2, -6);
    assert(cached_publish_matches_plain(cache, more, mixed));
    assert(cache.compute_count == 1 && cache.reuse_count == 2);

    // A relation that changes the lattice changes the basis: recomputed, and
    // the cache now holds the new basis.
    sflint::FmpzMat grown(6, 3);
    for (slong i = 0; i < 5; ++i) {
        for (slong j = 0; j < 3; ++j) {
            sflint::fmpz_set(sflint::fmpz_mat_entry(grown, i, j),
                             entry_const(more, i, j));
        }
    }
    set_entry_si(grown, 5, 2, 3);
    sflint::FmpzMat grown_basis(3, 3);
    set_entry_si(grown_basis, 0, 0, 2);
    set_entry_si(grown_basis, 0, 2, 2);
    set_entry_si(grown_basis, 1, 1, 3);
    set_entry_si(grown_basis, 2, 2, 3);
    assert(cached_publish_matches_plain(cache, grown, grown_basis));
    assert(cache.compute_count == 2 && cache.reuse_count == 2);
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(cache.hnf_basis),
                                  sflint::FmpzMatConstRef(grown_basis)));

    // The old basis is no longer cached.
    assert(cached_publish_matches_plain(cache, more, basis));
    assert(cache.compute_count == 3 && cache.reuse_count == 2);

    // A basis of another size never matches.
    sflint::FmpzMat small_relations(2, 2);
    set_entry_si(small_relations, 0, 0, 2);
    set_entry_si(small_relations, 1, 1, 6);
    assert(cached_publish_matches_plain(cache, small_relations,
                                        small_relations));
    assert(cache.compute_count == 4 && cache.reuse_count == 2);
    return 0;
}

// The witness rows are the rows of diag(U_H, I) * T that the accessors read,
// where T * relations = [H; 0] is fmpz_mat_hnf_transform and U_H * H * V = S
// is hnf_basis_smith_transform of H (the route the publish and the witness
// path both use; test_unit_pivot_smith_matches_full checks it against
// fmpz_mat_snf_transform).  ensure_left_transform forms only the top
// rows with a non-unit Smith entry; this reference forms the full product and
// checks that every row the accessors return is the same exact integer row.
// The pseudo-random relation matrices give groups whose Smith form mixes
// unit and non-unit entries, or has only unit entries.  Every-entry-non-unit
// (rank == n) is not generated and not asserted.
int test_witness_rows_match_full_left_product() {
    const slong m = 7;
    const slong n = 4;
    ulong state = 12345;
    slong mixed = 0;
    slong trivial = 0;
    slong compared = 0;
    for (slong trial = 0; trial < 60; ++trial) {
        sflint::FmpzMat relations(m, n);
        for (slong i = 0; i < m; ++i) {
            for (slong j = 0; j < n; ++j) {
                state = state * UWORD(6364136223846793005) +
                        UWORD(1442695040888963407);
                const slong value = static_cast<slong>((state >> 33) % 9) - 4;
                // Scale some columns so non-unit invariants are common.
                set_entry_si(relations, i, j,
                             (j >= 2 && trial % 3 != 2) ? 2 * value : value);
            }
        }
        silex::FiniteAbelianGroup group;
        if (!group.set_relation_matrix(
                    sflint::FmpzMatConstRef(relations))) {
            continue;
        }

        sflint::FmpzMat hnf_full(m, n);
        sflint::FmpzMat transform(m, m);
        ::fmpz_mat_hnf_transform(hnf_full.raw(), transform.raw(),
                                 relations.raw());
        sflint::FmpzMat hnf_top(n, n);
        for (slong i = 0; i < n; ++i) {
            for (slong j = 0; j < n; ++j) {
                sflint::fmpz_set(
                        sflint::fmpz_mat_entry(sflint::FmpzMatRef(hnf_top),
                                               i, j),
                        entry_const(hnf_full, i, j));
            }
        }
        sflint::FmpzMat snf(n, n);
        sflint::FmpzMat snf_left(n, n);
        sflint::FmpzMat snf_right(n, n);
        silex::detail::hnf_basis_smith_transform(
                sflint::FmpzMatRef(snf), sflint::FmpzMatRef(snf_left),
                sflint::FmpzMatRef(snf_right),
                sflint::FmpzMatConstRef(hnf_top));
        sflint::FmpzMat transform_top(n, m);
        for (slong i = 0; i < n; ++i) {
            for (slong j = 0; j < m; ++j) {
                sflint::fmpz_set(sflint::fmpz_mat_entry(
                                         sflint::FmpzMatRef(transform_top),
                                         i, j),
                                 entry_const(transform, i, j));
            }
        }
        sflint::FmpzMat full_top(n, m);
        sflint::fmpz_mat_mul(sflint::FmpzMatRef(full_top),
                             sflint::FmpzMatConstRef(snf_left),
                             sflint::FmpzMatConstRef(transform_top));

        const slong rank = group.invariant_count();
        sflint::FmpzMat combinations(rank, m);
        assert(group.invariant_generator_relation_matrix(
                sflint::FmpzMatRef(combinations)));
        slong k = 0;
        for (slong i = 0; i < n; ++i) {
            if (sflint::fmpz_is_pm1(entry_const(snf, i, i))) {
                continue;
            }
            const bool negate = sflint::fmpz_sgn(entry_const(snf, i, i)) < 0;
            for (slong j = 0; j < m; ++j) {
                sflint::Fmpz expected;
                sflint::fmpz_set(sflint::FmpzRef(expected),
                                 entry_const(full_top, i, j));
                if (negate) {
                    sflint::fmpz_neg(sflint::FmpzRef(expected),
                                     sflint::FmpzConstRef(expected));
                }
                assert(sflint::fmpz_equal(entry_const(combinations, k, j),
                                          sflint::FmpzConstRef(expected)));
            }
            ++k;
        }
        assert(k == rank);

        sflint::FmpzMat kernel(m - n, m);
        assert(group.relation_kernel_matrix(sflint::FmpzMatRef(kernel)));
        for (slong i = 0; i < m - n; ++i) {
            for (slong j = 0; j < m; ++j) {
                assert(sflint::fmpz_equal(entry_const(kernel, i, j),
                                          entry_const(transform, n + i, j)));
            }
        }
        assert(generator_relation_identity(group));
        assert(relation_kernel_identity(group));
        ++compared;
        if (rank > 0 && rank < n) {
            ++mixed;
        }
        if (rank == 0) {
            ++trivial;
        }
    }
    assert(compared > 0);
    assert(mixed > 0);
    assert(trivial > 0);
    return 0;
}

namespace {

bool is_unimodular(const sflint::FmpzMat& matrix) noexcept {
    sflint::Fmpz det;
    ::fmpz_mat_det(det.raw(), matrix.raw());
    return sflint::fmpz_is_pm1(sflint::FmpzConstRef(det));
}

// Pseudo-random n x n basis in row Hermite normal form whose diagonal is 1
// except at the positions with nonunit[i] set, where it is 2..7; entries
// above a pivot d_j are drawn from [0, d_j), so rows of unit pivots carry
// nonzero entries in the nonunit columns.
sflint::FmpzMat random_unit_pivot_hnf(slong n,
                                      const std::vector<bool>& nonunit,
                                      ulong& state) noexcept {
    const auto next = [&state]() noexcept {
        state = state * UWORD(6364136223846793005) +
                UWORD(1442695040888963407);
        return state >> 33;
    };
    sflint::FmpzMat basis(n, n);
    for (slong j = 0; j < n; ++j) {
        const slong pivot =
                nonunit[static_cast<std::size_t>(j)]
                        ? static_cast<slong>(2 + next() % 6)
                        : 1;
        set_entry_si(basis, j, j, pivot);
        for (slong i = 0; i < j; ++i) {
            set_entry_si(basis, i, j, static_cast<slong>(next()) % pivot);
        }
    }
    return basis;
}

// U * H * V = S, U and V unimodular, S equal to the Smith form of the full
// fmpz_mat_snf_transform, and every row of H * V a combination of the rows
// of S (column j of H * V is a multiple of S_jj).
bool smith_transform_valid(const sflint::FmpzMat& basis,
                           const sflint::FmpzMat& snf,
                           const sflint::FmpzMat& left,
                           const sflint::FmpzMat& right,
                           const sflint::FmpzMat& full_snf) noexcept {
    const slong n = sflint::fmpz_mat_nrows(basis);
    sflint::FmpzMat basis_right(n, n);
    sflint::FmpzMat product(n, n);
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(basis_right),
                         sflint::FmpzMatConstRef(basis),
                         sflint::FmpzMatConstRef(right));
    sflint::fmpz_mat_mul(sflint::FmpzMatRef(product),
                         sflint::FmpzMatConstRef(left),
                         sflint::FmpzMatConstRef(basis_right));
    if (!sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(product),
                                sflint::FmpzMatConstRef(snf)) ||
        !sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(snf),
                                sflint::FmpzMatConstRef(full_snf)) ||
        !is_unimodular(left) || !is_unimodular(right)) {
        return false;
    }
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            if (!sflint::fmpz_divisible(
                        sflint::FmpzConstRef(sflint::fmpz_mat_entry(
                                sflint::FmpzMatConstRef(basis_right), i, j)),
                        entry_const(snf, j, j))) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

// hnf_basis_smith_transform splits off the unit pivots of an HNF basis and
// runs fmpz_mat_snf_transform on the rest.  On random HNF bases with
// identity blocks it must give the Smith form of the full
// fmpz_mat_snf_transform with valid unimodular transforms, and a group
// published from the basis must have those invariants.  Without unit pivots
// it is fmpz_mat_snf_transform of the whole basis, transforms included.
int test_unit_pivot_smith_matches_full() {
    using silex::detail::FiniteAbelianGroupAccess;
    ulong state = 98765;
    slong reduced = 0;
    slong whole = 0;
    slong trivial = 0;
    for (slong trial = 0; trial < 120; ++trial) {
        const slong n = 1 + trial % 12;
        std::vector<bool> nonunit(static_cast<std::size_t>(n), false);
        slong nonunit_count = 0;
        for (slong i = 0; i < n; ++i) {
            state = state * UWORD(6364136223846793005) +
                    UWORD(1442695040888963407);
            // Every eighth trial has no unit pivot; every eighth (shifted)
            // has only unit pivots; the rest have about one in four.
            const bool pick = trial % 8 == 0
                                      ? true
                                      : (trial % 8 != 4 &&
                                         (state >> 33) % 4 == 0);
            nonunit[static_cast<std::size_t>(i)] = pick;
            nonunit_count += pick ? 1 : 0;
        }
        const sflint::FmpzMat basis =
                random_unit_pivot_hnf(n, nonunit, state);
        assert(::fmpz_mat_is_in_hnf(basis.raw()) != 0);

        sflint::FmpzMat full_snf(n, n);
        sflint::FmpzMat full_left(n, n);
        sflint::FmpzMat full_right(n, n);
        sflint::fmpz_mat_snf_transform(sflint::FmpzMatRef(full_snf),
                                       sflint::FmpzMatRef(full_left),
                                       sflint::FmpzMatRef(full_right),
                                       sflint::FmpzMatConstRef(basis));
        sflint::FmpzMat snf(n, n);
        sflint::FmpzMat left(n, n);
        sflint::FmpzMat right(n, n);
        silex::detail::hnf_basis_smith_transform(
                sflint::FmpzMatRef(snf), sflint::FmpzMatRef(left),
                sflint::FmpzMatRef(right), sflint::FmpzMatConstRef(basis));
        assert(smith_transform_valid(basis, snf, left, right, full_snf));
        assert(smith_transform_valid(basis, full_snf, full_left, full_right,
                                     full_snf));
        if (nonunit_count == n) {
            assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(left),
                                          sflint::FmpzMatConstRef(full_left)));
            assert(sflint::fmpz_mat_equal(
                    sflint::FmpzMatConstRef(right),
                    sflint::FmpzMatConstRef(full_right)));
            ++whole;
        } else if (nonunit_count == 0) {
            assert(matrix_is_identity(snf));
            ++trivial;
        } else {
            ++reduced;
        }

        // Publish from the basis with two extra relation rows (sums of basis
        // rows) and compare the invariants with the full Smith form.
        sflint::FmpzMat relations(n + 2, n);
        for (slong i = 0; i < n; ++i) {
            for (slong j = 0; j < n; ++j) {
                sflint::fmpz_set(sflint::fmpz_mat_entry(relations, i, j),
                                 entry_const(basis, i, j));
                sflint::fmpz_add(
                        sflint::fmpz_mat_entry(relations, n + i % 2, j),
                        entry_const(relations, n + i % 2, j),
                        entry_const(basis, i, j));
            }
        }
        silex::FiniteAbelianGroup group;
        assert(FiniteAbelianGroupAccess::set_relation_matrix_with_hnf_basis(
                group, sflint::FmpzMatConstRef(relations),
                sflint::FmpzMatConstRef(basis)));
        slong k = 0;
        for (slong i = 0; i < n; ++i) {
            if (entry_is_si(full_snf, i, i, 1)) {
                continue;
            }
            sflint::Fmpz value;
            assert(group.invariant(sflint::FmpzRef(value), k));
            assert(sflint::fmpz_equal(sflint::FmpzConstRef(value),
                                      entry_const(full_snf, i, i)));
            ++k;
        }
        assert(k == group.invariant_count());
        assert(generator_relation_identity(group));
        assert(relation_kernel_identity(group));
    }
    assert(reduced > 0 && whole > 0 && trivial > 0);
    return 0;
}

int main() {
    assert(test_basic() == 0);
    assert(test_kernel_basis_check_rejects_bad_kernels() == 0);
    assert(test_copy_swap_failure_preserves_output() == 0);
    assert(test_element_reduce_and_coordinates() == 0);
    assert(test_nondiagonal_and_generator_relations() == 0);
    assert(test_trivial_snf_coordinates() == 0);
    assert(test_witness_request_keeps_generators() == 0);
    assert(test_diagonal_witness_keeps_generators() == 0);
    assert(test_internal_basis_setter_witnesses() == 0);
    assert(test_lazy_right_transform_inverse() == 0);
    assert(test_smith_cache_reuses_unchanged_basis() == 0);
    assert(test_witness_rows_match_full_left_product() == 0);
    assert(test_unit_pivot_smith_matches_full() == 0);
    return 0;
}
