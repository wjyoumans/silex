#pragma once

#include <silex/element.hpp>

namespace silex::detail {

bool ensure_parent(Element& out, const NumberField& field) noexcept;

// A root in `field` of the cyclotomic polynomial Phi_m (m >= 3), that is, a
// primitive m-th root of unity, found by Hensel lifting from a good prime as
// in reference v0.38.6 `_roots_hensel(Phi_m, max_roots = 1, is_normal = true,
// root_bound = ones)` (see src/element/element.cpp).  Requires a monic
// integral defining polynomial and `root` bound to `field`.  Writes `root`
// only after checking its exact order m.  False means only that no root was
// found: it covers absence and unsupported or unresolved cases alike, and is
// never a proof that `field` has no primitive m-th root of unity.  Callers
// must fail closed on false, never conclude absence from it.
bool cyclotomic_root_hensel(Element& root,
                            const NumberField& field,
                            ulong m) noexcept;

}  // namespace silex::detail
