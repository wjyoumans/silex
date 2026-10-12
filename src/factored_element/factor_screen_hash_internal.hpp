#pragma once

#include <silex/element.hpp>
#include <silex/factored_element.hpp>

#include <flint/flint.h>
#include <flint/nf.h>

#include <cstddef>
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

// Internal access to the per-factor tags a FactoredElement keeps beside its
// factor list.  Each stored factor carries an origin and its screen hash.
// A factor stored by push or set_element gets a fresh origin, never reused
// in the process; a factor copied from another factored element (a product,
// quotient, inverse, power, root, a copy-on-write copy, or push_factor)
// keeps its origin.  Stored factors are never modified, so two factors with
// the same origin are copies of one element and Element::equal holds between
// them: an equal origin confirms equality without comparing coefficients.
// Different origins say nothing; such factors may still be equal.
class FactoredElementAccess {
public:
    // factor_screen_hash of factor `index` of `value` (with the factor's
    // parent field), computed once per stored factor and carried with its
    // copies.  `index` must be below value.length().
    static std::uint64_t factor_hash(const FactoredElement& value,
                                     std::size_t index) noexcept;

    // The nonzero origin of factor `index` of `value`.  `index` must be
    // below value.length().
    static std::uint64_t factor_origin(const FactoredElement& value,
                                       std::size_t index) noexcept;

    // out.push(source.factors()[index].factor, exponent), with the same
    // result and failures, except that the copied factor keeps its origin
    // and screen hash.  False when `index` is not below source.length().
    static bool push_factor(FactoredElement& out,
                            const FactoredElement& source,
                            std::size_t index,
                            slong exponent) noexcept;

private:
    friend class silex::FactoredElement;

    // The checks and copy of FactoredElement::push, storing the copy with
    // `origin` and `hash`; origin 0 gives the factor a fresh origin, and
    // hash 0 leaves the screen hash to be computed when first asked for.
    static bool push_tagged(FactoredElement& out,
                            const Element& factor,
                            slong exponent,
                            std::uint64_t origin,
                            std::uint64_t hash) noexcept;
};

}  // namespace silex::detail
