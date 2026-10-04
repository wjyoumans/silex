#pragma once

#include <silex/abelian_group.hpp>

namespace silex::detail {

// Noninstalled hook for publishing a presentation from a precomputed basis of
// the relation lattice, so callers that already maintain an HNF checkpoint
// (ClassGroupContext) avoid a second HNF of the full relation matrix.
//
// Precondition: hnf_basis is an n x n basis of the row lattice of relations.
// Only its shape is checked here; a basis in another form is first brought to
// Hermite normal form.  A basis that does not span the relation lattice
// publishes invariants of that basis, and witness requests then fail closed
// because the HNF of the full relation matrix differs from the stored basis.
class FiniteAbelianGroupAccess {
public:
    static bool set_relation_matrix_with_hnf_basis(
            FiniteAbelianGroup& group,
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis) noexcept {
        return group.set_relation_matrix_with_hnf_basis(relations, hnf_basis);
    }
};

}  // namespace silex::detail
