#pragma once

#include <silex/element.hpp>
#include <silex/factored_element.hpp>

#include <flint/flint.h>
#include <flint/nf.h>

#include <cstddef>
#include <cstdint>
#include <vector>

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

// Element::equal between stored factors of different origins (see
// FactoredElementAccess), remembered by the unordered pair of origins so a
// repeated comparison of the same pair is not made again.
//
// The cached result never goes stale.  An origin is minted once per process
// by an atomic counter and never reused (fresh_factor_origin in
// factored_element.cpp), only a copy of the same stored factor carries it
// (push_factor and the factor-list copies), and stored factors are never
// modified; so every factor with a given origin has the same coefficients
// and the same parent field, and Element::equal, which depends only on
// those, gives one answer for a pair of origins.  A true and a false answer
// are both kept.  A factor with origin 0 (stored from a bare Element, as in
// the relation-unit log memo) is compared every time and nothing is kept.
//
// The owner bounds the cache's lifetime; it holds one entry per pair it has
// compared, so it grows with the comparisons made, never past them.  Not
// thread-safe: one owner uses it at a time.
class FactorOriginPairCache {
public:
    // left.equal(right), for factors whose origins are left_origin and
    // right_origin.  Equal nonzero origins give true without a comparison.
    bool equal(const Element& left,
               std::uint64_t left_origin,
               const Element& right,
               std::uint64_t right_origin) noexcept;

    // Number of origin pairs held, and number of calls answered from them.
    std::size_t size() const noexcept {
        return count_;
    }
    std::size_t hits() const noexcept {
        return hits_;
    }

private:
    // An empty slot has low == 0; a held pair has 0 < low <= high.
    struct Slot {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
        bool equal = false;
    };

    // Open addressing with linear probing over a power-of-two number of
    // slots, kept at most half full.
    std::vector<Slot> slots_;
    std::size_t count_ = 0;
    std::size_t hits_ = 0;
};

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

    // value.logarithmic_embedding(out, embeddings, mode, precision,
    // diagnostics), with the same result and failures; the comparisons
    // that merge equal factors of different origins go through `pairs`
    // when it is not null.
    static bool logarithmic_embedding(
            const FactoredElement& value,
            flint::ArbVecRef out,
            EmbeddingContext& embeddings,
            LogEmbeddingMode mode,
            slong precision,
            const DiagnosticsContext* diagnostics,
            FactorOriginPairCache* pairs) noexcept;

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
