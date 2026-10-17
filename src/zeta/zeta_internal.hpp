#pragma once

#include <silex/diagnostics.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/zeta.hpp>

#include <cstddef>
#include <vector>

namespace silex {
class Element;
class FactorBase;
}

namespace silex::detail {

// What the Belabas-Friedman evaluation knows about reading residue degrees
// of a maximal order from its field's defining polynomial: the reduced
// coefficients, the order discriminant, and, for an order other than the
// equation order, its index over the equation order.
struct BfResidueDegreeContext {
    const fmpq_poly_struct* defining_polynomial = nullptr;
    std::vector<flint::Fmpz> defining_coefficients;
    slong defining_degree = -1;
    flint::Fmpz order_discriminant;
    flint::Fmpz equation_index;
    const FactorBase* residue_degree_base = nullptr;
    bool has_defining_polynomial_route = false;
    bool has_order_discriminant = false;
    bool check_equation_index = false;
    bool has_complete_residue_degree_base = false;
};

// The residue degrees of the prime ideals of a maximal order above a
// rational prime p, read from the distinct-degree factorisation of the
// field's defining polynomial modulo p (Dedekind-Kummer), without
// decomposing p.  It is the route the Belabas-Friedman evaluation takes
// (bf_residue_degrees_from_defining_polynomial_ui in zeta.cpp), under the
// same gate: the defining polynomial is monic integral, and p divides
// neither the order discriminant nor the index of the equation order in
// the order.  For every other prime the caller must decompose p.
class DefiningPolynomialResidueDegrees {
public:
    // Prepares the route for the maximal order `order`.  Returns false,
    // leaving every later query false, when the route is unavailable for
    // it (not maximal, defining polynomial not monic integral, or the
    // discriminant or index cannot be computed).
    bool init(const Order& order) noexcept;
    // Sets out to the residue degrees above p, one entry per prime ideal,
    // and returns true when the route applies to p.  Returns false when
    // init has not succeeded, the gate rejects p, or the computation
    // fails; out is then unspecified.  `order` must be the order init was
    // called with.
    bool residue_degrees(std::vector<slong>& out,
                         const Order& order,
                         ulong p) noexcept;

private:
    BfResidueDegreeContext context_;
    std::vector<slong> degree_scratch_;
    bool ready_ = false;
};

struct ZetaBfResidueDegreeCacheEntry {
    ulong p = 0;
    std::size_t offset = 0;
    std::size_t length = 0;
};

struct ZetaBfLinearFactorCountEntry {
    ulong p = 0;
    slong count = 0;
};

// Caches, per prime `p`, the residue degrees of the primes of an order above
// `p` (used by the Belabas-Friedman per-prime term in `bf_term`).  Entries
// are keyed only by `p`: there is no order (or computation) identity stored
// alongside them.
//
// Decision (2026-09-26): the cache intentionally carries no order key.  Each
// order, and each computation that needs one, must own and pass its own
// `ZetaBfResidueDegreeCache` instance; callers must never reuse one cache
// across two different orders.  The two owning callers already follow this
// by pairing the cache with an explicit order and checking it before reuse:
//
// - `RelationFactorBasePlan`
//   (class_group/relation_factor_base_plan_internal.hpp) stores its `Order`
//   next to the cache, and `relation_factor_base_plan_residue_degrees`
//   returns the cache only when that order matches.
// - `AnalyticClassRegulatorCache` (order_unit/compute_internal.hpp) tracks
//   the order it last cached for and resets the cache whenever a different
//   order is seen (`reset_for_order_if_needed_`).
//
// The analytic-finish and validation paths (class_group/analytic_finish.cpp,
// order_unit/validation.cpp) only forward a cache pointer already owned this
// way; they do not own one themselves.  A caller that needs to cache residue
// degrees for more than one order must construct one
// `ZetaBfResidueDegreeCache` per order.
//
// The same per-order cache also holds, sorted by `p`, the number of linear
// factors of the defining polynomial modulo `p` (the number of degree-one
// primes above `p`) for the primes where the Belabas-Friedman sums need only
// that count: the Theorem 7 cutoff scan and `bf_term` at every cutoff.  A
// count is valid only for primes passing
// `bf_residue_degree_context_uses_defining_polynomial`, which is checked
// again before every read.  It follows the same one-order rule as the
// residue degrees.
struct ZetaBfResidueDegreeCache {
    std::vector<ZetaBfResidueDegreeCacheEntry> entries;
    std::vector<slong> residue_degrees;
    std::size_t lookup_hint = 0;
    std::vector<ZetaBfLinearFactorCountEntry> linear_factor_counts;
    std::size_t linear_factor_count_hint = 0;
};

// Audit data of the Belabas-Friedman evaluation behind a default-route
// zeta product: the truncation error bound added to the log residue (valid
// under GRH, Belabas-Friedman 2015, Theorem 1), the cutoff used, the
// maximum cutoff the default route allows, and the requested and working
// precisions.  The default route does not require the error bound to meet
// the precision target, so `cutoff` may equal `max_cutoff`.  The grh
// index-one route (zeta_class_regulator_product_index_one_with_diagnostics)
// fills the same fields, with the error bound of Belabas-Friedman 2015,
// Theorem 7, below 1/4.
struct ZetaBfRouteAudit {
    flint::Arb error_bound;
    ulong cutoff = 0;
    ulong max_cutoff = 0;
    slong requested_precision = 0;
    slong work_precision = 0;
};

// `out` is the GRH generation bound min(BDF, Bach), or Bach's bound alone
// when the BDF criterion cannot be computed or decided (for example for a
// non-monic defining polynomial).  `bach_selected`, when given, is set on
// success to whether Bach's bound gave `out`: strictly below the BDF bound,
// or the BDF bound unavailable; on a tie the BDF bound is reported.
bool grh_factor_base_bound_with_diagnostics(
        flint::FmpzRef out,
        const Order& order,
        const DiagnosticsContext* diagnostics,
        bool* bach_selected = nullptr) noexcept;

// Outcome of the Belabas-Diaz y Diaz-Friedman search behind the GRH
// generation bound.
enum class BdfSearchOutcome {
    // The criterion holds at the bound written to `out`.
    found,
    // The criterion failed at a cutoff at least the stop bound, so the
    // bound the full search would find exceeds the stop bound.
    stopped,
    // The criterion could not be computed or decided; `out` is unchanged.
    failed,
};

// The BDF search that `grh_factor_base_bound_with_diagnostics` runs with
// Bach's bound as `stop_bound`; `stop_bound` may be null (no stop).
// `evaluations`, when given, receives the number of evaluations of the
// criterion's right side.  On `found`, `out` is the BDF bound; otherwise
// `out` is unchanged.
BdfSearchOutcome grh_bdf_factor_base_bound_search(
        flint::FmpzRef out,
        const Order& order,
        const flint::Fmpz* stop_bound,
        const DiagnosticsContext* diagnostics,
        slong* evaluations = nullptr) noexcept;

class TestHookAccess;

// Test entry for the batched linear factor counts of the Belabas-Friedman
// prime loop; only detail::TestHookAccess (test/test_hook_access.hpp) can
// call it.  For each prime p in [first, end), in increasing order,
// `primes` receives p, `batch_counts` the number of linear factors of the
// defining polynomial of `order` modulo p from the batched word pass, and
// `per_prime_counts` the same number from the per-prime nmod path; each is
// -1 where that path does not apply to p.  Returns false when `order` is not
// a maximal order with a basis.
class ZetaBfTestEntry {
    friend class TestHookAccess;
    static bool linear_factor_counts(std::vector<ulong>& primes,
                                     std::vector<slong>& batch_counts,
                                     std::vector<slong>& per_prime_counts,
                                     const Order& order,
                                     ulong first,
                                     ulong end) noexcept;
};

// `unconditional`, when given, is set to true exactly when the product came
// from an unconditional route: degree one, or the quadratic L(1, chi) route
// (explicit quadratic backend, maximal order, |D| fitting a ulong).  It is
// false when the product came from the Belabas-Friedman fallback, whose error
// bound assumes GRH for zeta_K and zeta_Q (Belabas-Friedman 2015, Theorem 1).
// `bf_audit`, when given, receives the audit data of that fallback; it is
// written only on success on the fallback route and left unchanged
// otherwise.  `accept_grh_residue`, for a caller serving a `grh` request,
// skips the quadratic L(1, chi) route and takes the Belabas-Friedman residue
// in every degree above one (`unconditional` is then false); the default
// keeps the unconditional route first.
bool zeta_class_regulator_product_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        bool* unconditional = nullptr,
        ZetaBfRouteAudit* bf_audit = nullptr,
        bool accept_grh_residue = false) noexcept;

// True when `zeta_class_regulator_product` has an unconditional route for
// `order`: degree one, or the quadratic L(1, chi) route preconditions.  The
// quadratic route can still fall back to Belabas-Friedman if the L-value
// evaluation fails, so a caller that holds the value should use the route
// reported by `zeta_class_regulator_product_with_diagnostics` instead.
bool zeta_unconditional_route_available(const Order& order) noexcept;

bool zeta_class_regulator_product_bf_audit_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

bool zeta_class_regulator_product_validation_with_diagnostics(
        flint::ArbRef out,
        flint::ArbRef error_bound,
        ulong& cutoff,
        slong& work_precision,
        const Order& order,
        ulong max_cutoff,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

// The Belabas-Friedman hR for the grh analytic index-one test: the log
// residue f_K(X) at the least cutoff X (a multiple of 9 up to the default
// route's maximum) whose Belabas-Friedman 2015 Theorem 7 bound, with
// sigma = 1 + 1/sqrt(log D), is below 1/4 - 2^-20, so that the log residue
// has radius below 1/4, enough for the index-one test (Belabas-Friedman
// 2015, p. 1, footnote 1: below (1/2) log 2).  Valid under GRH for zeta_K
// and zeta_Q, in every degree above one, including a quadratic order that
// also has the unconditional L(1, chi) route: its only caller serves grh
// requests, which accept a GRH-conditional hR.  `audit` receives the
// Theorem 7 error bound, the cutoff, the default route's maximum cutoff,
// and the precisions, and is written only on success.  Fails, leaving the
// outputs unchanged, in degree one or when no cutoff up to the maximum
// reaches the target; callers then use the default route.
bool zeta_class_regulator_product_index_one_with_diagnostics(
        flint::ArbRef out,
        ZetaBfRouteAudit& audit,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        const FactorBase* residue_degree_base = nullptr,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr) noexcept;

// Test entry: the Belabas-Friedman 2015 Theorem 7 bound used above, at
// `cutoff` (a multiple of 9, at least 18), with the prime-ideal sum over
// the norms below `cutoff`.
bool zeta_bf_theorem7_error_bound(flint::ArbRef out,
                                  const Order& order,
                                  ulong cutoff,
                                  slong precision) noexcept;

bool zeta_bf_audit_cutoff_available(const Order& order,
                                    ulong max_cutoff,
                                    slong precision) noexcept;

bool class_regulator_product_estimate_with_diagnostics(
        flint::ArbRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics,
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        flint::Fmpz* torsion_order = nullptr,
        Element* torsion_generator = nullptr) noexcept;

// The character of the quadratic L(1, chi) route: the Kronecker symbol
// (D / n), which for a fundamental discriminant D is the real primitive
// character modulo |D|.
int quadratic_character(flint::FmpzConstRef discriminant, ulong n) noexcept;

// Test controls and counters of the Taylor stepping of the incomplete gamma
// values Gamma(a, z) in `quadratic_dirichlet_l1`.  The counters are indexed
// by shift = 2a: the sums use a = 1/2 and a = 0 for real D, a = 1 and
// a = 1/2 for imaginary D.
struct QuadraticAfeStepAudit {
    // False turns stepping off, so that every value is evaluated directly,
    // as the evaluation did before stepping was added.
    bool stepping = true;
    // Blocks started, and values taken from a block (an anchor carried
    // over from the previous block included).
    ulong blocks[3] = {0, 0, 0};
    ulong stepped_values[3] = {0, 0, 0};
};

// L(1, chi_D) for a fundamental discriminant D with |D| >= 3 fitting a
// ulong, by the approximate functional equation with chi_D evaluated as a
// Kronecker symbol.  Returns false when D is outside that range or the
// enclosure is not finite; D is not checked to be fundamental.  `audit`,
// when given, is used and updated as described above.
bool quadratic_dirichlet_l1(flint::ArbRef out,
                            flint::FmpzConstRef discriminant,
                            slong precision,
                            QuadraticAfeStepAudit* audit = nullptr) noexcept;

// L(1, chi_D) as `quadratic_dirichlet_l1`, with the same two sums, tail
// bound and per-term working precisions, but sized by an absolute tolerance
// instead of FLINT's mapping from a target precision: the first sum is
// truncated at tolerance 2^log2_tolerance and the second at FLINT's multiple
// of it (about 2^log2_tolerance sqrt(|D| / pi)), and both sums and their
// combination are accumulated at `precision`.  The result is a rigorous
// enclosure for any such tolerance and precision; only its radius depends on
// them.  Returns false when D is outside the range of
// `quadratic_dirichlet_l1`, `precision` is not positive, `log2_tolerance` is
// positive, or the enclosure is not finite.
bool quadratic_dirichlet_l1_with_tolerance(
        flint::ArbRef out,
        flint::FmpzConstRef discriminant,
        slong log2_tolerance,
        slong precision,
        QuadraticAfeStepAudit* audit = nullptr) noexcept;

// A fault forced on every block the stepper of `quadratic_dirichlet_l1`
// plans, after planning and before the block starts, so that the start must
// reject the block (for tests of its fail-closed paths): the Cauchy radius
// R = 2^r raised to at least 2 z_n, R lowered to at most half the largest
// step v_max of the block, or the truncation bound replaced by infinity once
// it is computed.  Planning never produces these states.
enum class QuadraticAfeStepFault : unsigned char {
    none,
    radius_above_z,
    radius_below_step,
    infinite_truncation,
};

struct QuadraticAfeGammaStepCount {
    ulong blocks = 0;
    ulong stepped_values = 0;
    // Planned blocks whose start was rejected, by the check that failed:
    // R = 2^r below z_n, v_max below R, a finite truncation bound.
    ulong rejected_radius = 0;
    ulong rejected_step = 0;
    ulong rejected_truncation = 0;
};

// Test entry for the stepper of `quadratic_dirichlet_l1` alone: writes
// Gamma(shift / 2, pi n^2 / q) for n = first, ..., first + out.length() - 1
// to `out`, taking the branches of that sum at every index.  A block is
// planned with anchor tolerance 2^(target - guard) and anchored at a value
// evaluated directly at precision prec + guard; values inside a block have
// tolerance about 2^target at precision `prec`, and the last index of a
// block becomes the next anchor.  Blocks are at most `max_block` long
// (0: the default cap).  An index no block covers is evaluated directly.
// `fault` is forced on every planned block before it starts.
// Returns false on invalid arguments (shift > 2, q or first zero, an empty
// `out`, or a negative precision or guard).
bool quadratic_afe_gamma_stepped(
        flint::ArbVec& out,
        QuadraticAfeGammaStepCount& count,
        ulong shift,
        ulong q,
        ulong first,
        ulong max_block,
        double target,
        slong prec,
        slong guard,
        QuadraticAfeStepFault fault = QuadraticAfeStepFault::none) noexcept;

struct QuadraticAfeGammaStepBlock {
    ulong block = 0;
    slong terms = 0;
    slong radius_exponent = 0;
};

// Test entry for one block of that stepper: plans and starts a block
// anchored at n with anchor tolerance 2^target and coefficients at `prec`,
// the anchor evaluated directly at `prec`, and writes its length B, its
// number of Taylor terms M, the exponent r of its Cauchy radius R = 2^r and
// its truncation bound (Johansson 2021, Lemma 5, eq. 27) to `truncation`.
// Returns false on invalid arguments (shift > 2, q or n zero, `prec` not
// positive) or when no block starts at n.
bool quadratic_afe_gamma_step_block(QuadraticAfeGammaStepBlock& block,
                                    mag_struct* truncation,
                                    ulong shift,
                                    ulong q,
                                    ulong n,
                                    double target,
                                    slong prec) noexcept;

}  // namespace silex::detail
