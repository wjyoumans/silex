#pragma once

#include <silex/element.hpp>

#include <flint/nf.h>

#include <cstdint>

namespace silex::detail {

// Screening hash of an element's stored nf_elem coefficients; `field` is
// the element's parent FLINT field.  It only routes a search: equal stored
// values hash equally, distinct ones may collide, so every match must be
// confirmed with Element::equal.  Two equal elements that hash apart (a
// noncanonical stored form) are merely looked up separately.  Used by
// FactoredElement's distinct-factor scan and the relation-unit factored log
// memo.
std::uint64_t factor_screen_hash(const Element& factor,
                                 const nf_struct* field) noexcept;

}  // namespace silex::detail
