#pragma once

#include <silex/element.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpq.hpp>

#include <flint/nf.h>
#include <flint/nf_elem.h>

namespace silex::detail {

bool ensure_parent(Element& out, const NumberField& field) noexcept;

// True when z has exact multiplicative order m, that is, when z is a root of
// Phi_m: z^m = 1 and z^(m/l) != 1 for every prime l | m.  This is the
// certificate that w is at least m, used both to publish a cyclotomic
// Hensel root (src/element/hensel.cpp) and a candidate generator in
// src/unit/roots_of_unity.cpp.
bool has_exact_order(const Element& z, ulong m) noexcept;

// Entry points into the private Hensel-lifting and root-reconstruction
// machinery of `src/element/hensel.cpp` (see notes/reviews/T-043-reviewer.md
// N2): exact square/power roots over Q and quadratic fields, pure-power
// Hensel lifting for `Element::is_square` / `Element::is_power`, the generic
// trace and norm via the multiplication matrix, and the cyclotomic
// root-of-unity finder.  Moved out of `element.cpp` as one block because the
// cyclotomic search reuses the pure-power helpers.

// Trace and norm of `element` in `field` (degree `degree`) via its
// multiplication matrix, for backends without a direct computation.
void generic_trace(flint::FmpqRef out,
                   const nf_struct* field,
                   slong degree,
                   const nf_elem_t element) noexcept;
void generic_norm(flint::FmpqRef out,
                  const nf_struct* field,
                  slong degree,
                  const nf_elem_t element) noexcept;

// Exact square root of `input` when it is a rational constant. False means
// unsupported (input is not a rational constant); `is_square` is set only
// on success.
bool is_square_rational_constant(bool& is_square,
                                 Element& root,
                                 const Element& input) noexcept;

// Exact square root of `input` in a quadratic field Q(sqrt(radicand)).
bool is_square_quadratic(bool& is_square,
                         Element& root,
                         const Element& input,
                         const fmpz_t radicand) noexcept;

// Exact `exponent`-th root of `input` when it is a rational constant.
bool is_power_rational_constant(bool& is_power,
                                Element& root,
                                const Element& input,
                                flint::FmpzConstRef exponent) noexcept;

// Hensel-lifted `exponent`-th (respectively square) root of `input` in a
// monic integral field of degree < 10, verified exactly before publication.
bool pure_power_hensel_root(bool& is_power,
                            Element& root,
                            const Element& input,
                            slong exponent,
                            const DiagnosticsContext* diagnostics) noexcept;
bool pure_square_hensel_root(bool& is_square,
                             Element& root,
                             const Element& input,
                             const DiagnosticsContext* diagnostics) noexcept;

// True when a residue field already disproves `input` having an
// `exponent`-th (respectively square) root, before Hensel lifting is tried.
bool pure_power_residue_disproves(bool& is_power,
                                  const Element& input,
                                  slong exponent,
                                  const DiagnosticsContext* diagnostics) noexcept;
bool pure_square_residue_disproves(bool& is_square,
                                   const Element& input,
                                   const DiagnosticsContext* diagnostics) noexcept;

// A root in `field` of the cyclotomic polynomial Phi_m (m >= 3), that is, a
// primitive m-th root of unity, found by Hensel lifting from a good prime as
// in reference v0.38.6 `_roots_hensel(Phi_m, max_roots = 1, is_normal = true,
// root_bound = ones)` (see src/element/hensel.cpp).  Requires a monic
// integral defining polynomial and `root` bound to `field`.  Writes `root`
// only after checking its exact order m.  False means only that no root was
// found: it covers absence and unsupported or unresolved cases alike, and is
// never a proof that `field` has no primitive m-th root of unity.  Callers
// must fail closed on false, never conclude absence from it.
bool cyclotomic_root_hensel(Element& root,
                            const NumberField& field,
                            ulong m) noexcept;

}  // namespace silex::detail
