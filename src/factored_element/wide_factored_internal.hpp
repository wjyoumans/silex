#pragma once

#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpz_mat.hpp>

#include <span>
#include <vector>

namespace silex::detail {

// A batch of formal products prod_k bases[k]^exponents[i, k] whose exponents
// are arbitrary-size integers.  Row i of `exponents` is element i and column k
// belongs to bases[k]; the bases are pairwise distinct field elements.
//
// Source trace: The reference system's src/Misc/FactoredElem.jl stores a FacElem as a
// dictionary from base to ZZRingElem exponent, so `_transform`
// (FactoredElem.jl:339-365) never overflows.  The second reference, version 2.17.3,
// src/basemath/buch2.c:1214-1240 `vec_chinese_units` keeps its compact units
// in the same batch shape: a vector of bases and an integer exponent matrix
// (bnf_get_sunits), with unused bases removed by `ZM_remove_unused`
// (buch2.c:1104-1118).  Silex's public FactoredElement keeps word-size
// exponents, so this type is private storage for unit generators whose
// exponents do not fit in slong.
struct WideFactoredUnits {
    std::vector<Element> bases;
    flint::FmpzMat exponents{0, 0};
};

bool wide_factored_units_copy(WideFactoredUnits& out,
                              const WideFactoredUnits& input) noexcept;

// True when some exponent does not fit in slong.
bool wide_factored_units_exceed_slong(const WideFactoredUnits& units) noexcept;

// Sets out to the units prod_j units[j]^transform[i, j], i = 0..rank-1, the
// unit-basis transform of reference `_transform` (FactoredElem.jl:339-365)
// with integer exponents of any size.  Equal bases of the inputs are merged
// as in a FacElem dictionary, the new exponent matrix is transform * S where
// row j of S is the merged exponent row of units[j], and bases whose column
// is zero are removed (reference ZM_remove_unused).  Fails, leaving out
// unchanged, on mismatched sizes or parents.
bool wide_factored_units_transform(WideFactoredUnits& out,
                                   std::span<const FactoredElement> units,
                                   flint::FmpzMatConstRef transform) noexcept;

}  // namespace silex::detail
