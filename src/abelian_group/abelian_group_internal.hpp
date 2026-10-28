#pragma once

#include <silex/abelian_group.hpp>

namespace silex::detail {

// Smith form of the last HNF basis published through a cache-taking
// set_relation_matrix_with_hnf_basis call: the stored basis (after it was
// brought to Hermite normal form), the Smith form S and the right transform
// V that hnf_basis_smith_transform returned for it.  S and V are a function
// of the basis alone, so a later publish whose stored basis is equal takes
// them from here instead of running hnf_basis_smith_transform again.  Only
// results that published a group are recorded.
struct FiniteAbelianGroupSmithCache {
    flint::FmpzMat hnf_basis{0, 0};
    flint::FmpzMat snf{0, 0};
    flint::FmpzMat right_transform{0, 0};
    bool valid = false;
    // Number of successful publishes that took S and V from the cache or
    // kept those of a target group with an equal stored basis, and number
    // that ran hnf_basis_smith_transform and recorded them; read by tests.
    slong reuse_count = 0;
    slong compute_count = 0;
};

// Smith form S, left transform U and right transform V of a nonsingular square integer
// matrix with U * basis * V = S, U and V unimodular; all four are n x n and
// the outputs must not alias the input.  Unit pivots (columns of the basis
// equal to a unit vector, which for a basis in Hermite normal form are the
// columns with diagonal entry 1) are split off first and only the remaining
// block reaches fmpz_mat_snf_transform; see the definition for the source and
// the argument.  S is the Smith form fmpz_mat_snf_transform returns for the
// basis; U and V are valid transforms but can differ from the ones it
// returns.  Every Smith form of a FiniteAbelianGroup basis comes from here.
void hnf_basis_smith_transform(flint::FmpzMatRef snf,
                               flint::FmpzMatRef left,
                               flint::FmpzMatRef right,
                               flint::FmpzMatConstRef basis) noexcept;

// Noninstalled hook for publishing a presentation from a precomputed basis of
// the relation lattice, so callers that already maintain an HNF checkpoint
// (ClassGroupContext) avoid a second HNF of the full relation matrix.
//
// Precondition: hnf_basis is an n x n basis of the row lattice of relations.
// Only its shape is checked here; a basis in another form is first brought to
// Hermite normal form.  A basis that does not span the relation lattice
// publishes invariants of that basis, and witness requests then fail closed
// because the HNF of the full relation matrix differs from the stored basis.
// For ClassGroupContext::publish_presentation the precondition is invariant
// (I) of the relation-lattice audit in the comment above that function
// (src/class_group/class_group.cpp), which walks every mutation path of the
// relation rows and the row-module checkpoint; an expensive-level debug
// check there asserts that every relation row lies in the basis lattice.
//
// The overload taking a FiniteAbelianGroupSmithCache publishes exactly the
// group the plain overload publishes; it reuses the cached Smith form and
// right transform when the stored basis equals the cached basis, and records
// them otherwise (see also the in-place case below).
class FiniteAbelianGroupAccess {
public:
    static bool set_relation_matrix_with_hnf_basis(
            FiniteAbelianGroup& group,
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis) noexcept {
        return group.set_relation_matrix_with_hnf_basis(relations, hnf_basis);
    }

    static bool set_relation_matrix_with_hnf_basis(
            FiniteAbelianGroup& group,
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis,
            FiniteAbelianGroupSmithCache& cache) noexcept {
        return set_relation_matrix_with_hnf_basis_impl(
                group, relations, hnf_basis, &cache);
    }

    // Publishes into `group` exactly the group the cache overload above
    // publishes from the same relations and basis, but takes the storage of
    // `relations` and `hnf_basis` instead of copying them; both are left in
    // an unspecified valid state.  On failure `group` is unchanged.
    //
    // Both cache-taking entries publish in place when `group` already holds
    // a defined group whose stored basis equals the new stored basis: only
    // its relation rows are replaced.  Its Smith form, right transform,
    // invariants and (if already filled) inverse right transform are the
    // ones this publish would produce, since they depend on the stored basis
    // alone, and its lazily filled left transform, which depends on the
    // relation rows, is dropped.  Such a publish counts as a reuse and
    // leaves the cached basis as it was.
    static bool take_relation_matrix_with_hnf_basis(
            FiniteAbelianGroup& group,
            flint::FmpzMat& relations,
            flint::FmpzMat& hnf_basis,
            FiniteAbelianGroupSmithCache& cache) noexcept;

private:
    friend class silex::FiniteAbelianGroup;

    static bool publish_owned_relation_matrix(
            FiniteAbelianGroup& group,
            flint::FmpzMat& relations,
            flint::FmpzMat& hnf_basis,
            FiniteAbelianGroupSmithCache* cache) noexcept;

    static bool set_relation_matrix_with_hnf_basis_impl(
            FiniteAbelianGroup& group,
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis,
            FiniteAbelianGroupSmithCache* cache) noexcept;
};

}  // namespace silex::detail
