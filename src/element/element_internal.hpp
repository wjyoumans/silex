#pragma once

#include <silex/element.hpp>
#include <silex/factored.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpq.hpp>

#include <flint/nf.h>
#include <flint/nf_elem.h>

#include <atomic>
#include <span>

namespace silex::detail {

// Degree from which pure_power_residue_disproves runs its scan of
// kHighDegreeDisproofPrimes admissible primes (src/element/hensel.cpp), which
// pure_power_hensel_root also runs before lifting, and from which
// factored_power_residue_disproves screens factored leftovers.  It is the
// degree below which pure_power_hensel_root used to be ported; the
// pure-power lift now runs in every degree, and these residue routes keep
// their limit.
inline constexpr slong kHighDegreeResidueScanDegree = 10;

// Degree from which the power tests of two saturation consumers do not lift
// inside a SaturationLiftLimit scope: the square tests of grh relation
// saturation with units (saturation_relation_with_units in
// src/class_group/class_group.cpp, reached from
// ClassGroupContext::saturate_local_once_with_units_) and the kernel-row
// root tests of the local saturation proofs
// (OrderUnitGroup::prove_local_saturated_ and
// prove_local_saturated_stable_in_place_ in src/order_unit/saturation.cpp,
// through detail::kernel_row_root).  There pure_power_hensel_root returns
// "unsupported" in fields of this degree and higher, and the test is decided
// only by the residue disproof, as it was before the pure-power lift ran
// from degree 10 on.  It keeps those consumers' degree-10 limit, which the
// public Element and FactoredElement square and power tests no longer have;
// the limit avoids lifting candidates the residue scan could not disprove,
// which on high-degree saturation costs time without changing the published
// results.
inline constexpr slong kSaturationLiftDegree = 10;

class TestHookAccess;

// The exponents a SaturationLiftLimit scope covers.  Ordered: a scope never
// narrows the coverage of an enclosing one.
enum class SaturationLiftExponents : unsigned char {
    none,
    // Exponent 2 only (the grh square tests).
    squares,
    // Every exponent (the proof's kernel-row root tests).
    all,
};

// Scoped caller-side option of the saturation power tests (see
// kSaturationLiftDegree).  While an object of this class lives on a thread,
// pure_power_hensel_root on that thread returns "unsupported" for a covered
// exponent in fields of degree kSaturationLiftDegree and higher; other
// threads, and exponents the scope does not cover, are unaffected.
// "Unsupported" is a result every caller already handles (a failed lift
// returns it too), so answers that are given stay exact: "true" only after
// the exact root check, "false" only from a residue disproof.  Scopes nest;
// the destructor restores the previous state.
class SaturationLiftLimit {
public:
    explicit SaturationLiftLimit(SaturationLiftExponents exponents) noexcept
        : previous_(active_) {
        if (exponents > active_) {
            active_ = exponents;
        }
    }
    ~SaturationLiftLimit() { active_ = previous_; }
    SaturationLiftLimit(const SaturationLiftLimit&) = delete;
    SaturationLiftLimit& operator=(const SaturationLiftLimit&) = delete;

    // Whether a scope on the calling thread covers `exponent`.
    static bool covers(slong exponent) noexcept {
        return active_ == SaturationLiftExponents::all ||
               (active_ == SaturationLiftExponents::squares && exponent == 2);
    }

    // Counts one power test whose lift a scope skipped
    // (pure_power_hensel_root).  Only tests read the count, through
    // TestHookAccess.
    static void count_skipped_lift() noexcept {
        skipped_lifts_.fetch_add(1, std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;

    SaturationLiftExponents previous_;
    inline static thread_local SaturationLiftExponents active_ =
            SaturationLiftExponents::none;
    // Power tests whose lift a scope skipped, on any thread, since the
    // process started.
    inline static std::atomic<slong> skipped_lifts_{0};
};

// Switch that only detail::TestHookAccess (test/test_hook_access.hpp) can
// set.  It is off by default and no production code sets it.  While it is
// on, pure_power_hensel_root returns "unsupported" in fields of degree
// kHighDegreeResidueScanDegree and higher, as it did before the lift ran
// there, so tests can still reach the paths that follow an undecided root of
// a true power (for example the fail-closed saturation proofs).  A relaxed
// atomic, so a reader never races with a test setting it.
class ElementPowerTestSeams {
public:
    static bool high_degree_lift_off() noexcept {
        return high_degree_lift_off_.load(std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<bool> high_degree_lift_off_{false};
};

bool ensure_parent(Element& out, const NumberField& field) noexcept;

// True when z has exact multiplicative order m, that is, when z is a root of
// Phi_m: z^m = 1 and z^(m/l) != 1 for every prime l | m.  This is the
// certificate that w is at least m, used both to publish a cyclotomic
// Hensel root (src/element/hensel.cpp) and a candidate generator in
// src/unit/roots_of_unity.cpp.
bool has_exact_order(const Element& z, ulong m) noexcept;

// Per-thread counts of the Hensel data memo in src/element/hensel.cpp
// (norm-change constants of the Friedrich-Fieker lifting exponent, the
// power-root reconstruction lattice, the field embedding data of the
// conjugate bounds, per-radicand lifting exponents, f mod p setups of the
// prime search, the distinct-degree factorizations of f mod p of the residue
// disproof, 1 / f'(theta), whole pure_power_hensel_root answers, and whole
// answers of the degree >= 10 residue disproof scan).
// Diagnostic only: the memo never changes a computed value.
struct HenselDataCacheCounters {
    slong norm_change_constants_hits = 0;
    slong norm_change_constants_misses = 0;
    slong reconstruction_data_hits = 0;
    slong reconstruction_data_misses = 0;
    slong conjugate_bounds_field_hits = 0;
    slong conjugate_bounds_field_misses = 0;
    slong lifting_exponent_hits = 0;
    slong lifting_exponent_misses = 0;
    slong prime_setup_hits = 0;
    slong prime_setup_misses = 0;
    slong residue_setup_hits = 0;
    slong residue_setup_misses = 0;
    slong inverse_derivative_hits = 0;
    slong inverse_derivative_misses = 0;
    slong power_root_answer_hits = 0;
    slong power_root_answer_misses = 0;
    slong residue_disproof_answer_hits = 0;
    slong residue_disproof_answer_misses = 0;
};
HenselDataCacheCounters hensel_data_cache_counters() noexcept;

// Empties every memo counted above on the calling thread, so a test can
// compare a cold computation with a memoized one.  Diagnostic only.
void clear_hensel_data_caches() noexcept;

// Empties only the whole-answer memo of pure_power_hensel_root on the calling
// thread, so a test can reach the memos below it with repeated radicands.
// Diagnostic only.
void clear_power_root_answer_cache() noexcept;

// Empties only the whole-answer memo of the degree >= 10 residue disproof
// scan of pure_power_residue_disproves on the calling thread, so a test can
// reach the residue field setup memo below it with repeated inputs.
// Diagnostic only.
void clear_residue_disproof_answer_cache() noexcept;

// Entry points into the private Hensel-lifting and root-reconstruction
// machinery of `src/element/hensel.cpp`: exact square/power roots over Q
// and quadratic fields, pure-power Hensel lifting for `Element::is_square` /
// `Element::is_power`, the generic trace via the multiplication matrix and
// norm via FLINT, and the cyclotomic root-of-unity finder.  Moved out of
// `element.cpp` as one block because the cyclotomic search reuses the
// pure-power helpers.

// Trace of `element` in `field` (degree `degree`) via its multiplication
// matrix, and norm via FLINT's `nf_elem_norm`, for backends without a
// direct computation.
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

// Whether pure_power_hensel_root first looks for a non-power residue at
// further primes before lifting (see src/element/hensel.cpp,
// pre_lift_residue_disproves).  `disabled` exists only so tests can compare
// against the lift-only answers; production callers use the default.
enum class PreLiftResidueDisproof {
    enabled,
    disabled
};

// A caller-supplied way to compute N(input) for the exact verification of
// pure_power_hensel_root, for callers that can do so more cheaply than
// input.norm(), for example from a factored form of input whose factors are
// small.  `compute(norm, context)` writes exactly N(input) to `norm` and
// returns true, or returns false to fall back to input.norm().  It is called
// at most once per pure_power_hensel_root call, and only when a lifted
// candidate reaches the norm check.  A wrong norm can only turn a found root
// into an undecided result, never give a wrong definite answer, because a
// root is published only after exact powering; debug checks at
// DebugLevel::expensive compare it with input.norm().
struct PowerInputNormSource {
    bool (*compute)(flint::Fmpq& norm, const void* context) noexcept = nullptr;
    const void* context = nullptr;
};

// Hensel-lifted `exponent`-th (respectively square) root of `input` in a
// monic integral field of any degree, verified exactly before publication.
// A definite "not a power" comes only from an exact residue disproof.
// `input_norm_source`, when not null, supplies N(input) to the verification;
// it is not part of the whole-answer memo key, since it does not change the
// answer.
bool pure_power_hensel_root(
        bool& is_power,
        Element& root,
        const Element& input,
        slong exponent,
        const DiagnosticsContext* diagnostics,
        PreLiftResidueDisproof pre_lift = PreLiftResidueDisproof::enabled,
        const PowerInputNormSource* input_norm_source = nullptr) noexcept;

// Element::is_power(is_power, root, exponent, diagnostics) for `input`, with
// an optional source of N(input) passed to pure_power_hensel_root.
// Element::is_power is this with a null `input_norm_source`.
bool element_is_power(bool& is_power,
                      Element& root,
                      const Element& input,
                      flint::FmpzConstRef exponent,
                      const PowerInputNormSource* input_norm_source,
                      const DiagnosticsContext* diagnostics) noexcept;

// Whether input.norm() succeeds and equals `input_norm`: the debug
// cross-check of a norm supplied by a PowerInputNormSource.
bool power_input_norm_matches(const Element& input,
                              const flint::Fmpq& input_norm) noexcept;

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

// pure_power_residue_disproves in degree 10 and higher for the product
// prod factor^exponent of `factors` (nonzero elements of `parent`, positive
// exponents), computed from the residues of the factors without forming the
// product.  True (with `is_power` false) only where
// pure_power_residue_disproves on the evaluated product would return true; a
// false return decides nothing.  It examines only the first admissible prime
// of that scan (a resource bound; see kFactoredScreenAdmissiblePrimes in
// src/element/hensel.cpp), and abandons at the first prime where some factor
// is not p-integral, since the product's admissibility there is unknown.
bool factored_power_residue_disproves(
        bool& is_power,
        const NumberField& parent,
        std::span<const FactorPower<Element>> factors,
        slong exponent,
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
