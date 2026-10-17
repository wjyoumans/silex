#pragma once

#include <silex/abelian_group.hpp>

namespace silex::detail {

// Smith form of the last HNF basis published through a cache-taking
// set_relation_matrix_with_hnf_basis call: the stored basis (after it was
// brought to Hermite normal form), the Smith form S and the right transform
// V that fmpz_mat_snf_transform returned for it.  S and V are a function of
// the basis alone, so a later publish whose stored basis is equal takes them
// from here instead of running fmpz_mat_snf_transform again.  Only results
// that published a group are recorded.
struct FiniteAbelianGroupSmithCache {
    flint::FmpzMat hnf_basis{0, 0};
    flint::FmpzMat snf{0, 0};
    flint::FmpzMat right_transform{0, 0};
    bool valid = false;
    // Number of successful publishes that took S and V from the cache, and
    // number that ran fmpz_mat_snf_transform and recorded them; read by
    // tests.
    slong reuse_count = 0;
    slong compute_count = 0;
};

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
// them otherwise.
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

private:
    friend class silex::FiniteAbelianGroup;

    static bool set_relation_matrix_with_hnf_basis_impl(
            FiniteAbelianGroup& group,
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis,
            FiniteAbelianGroupSmithCache* cache) noexcept;
};

}  // namespace silex::detail
