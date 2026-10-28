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
// equation order, its index over the equation order.  `cyclotomic_index`
// is n when the defining polynomial is exactly the cyclotomic polynomial
// Phi_n (detail::cyclotomic_defining_polynomial_index, detected once when
// the context is set up), and 0 otherwise.  `real_cyclotomic_index` is n
// when it is instead exactly psi_n, the minimal polynomial of
// 2 cos(2 pi / n) that defines the maximal real subfield of Q(zeta_n)
// (detail::real_cyclotomic_defining_polynomial_index), and 0 otherwise.
struct BfResidueDegreeContext {
    const fmpq_poly_struct* defining_polynomial = nullptr;
    std::vector<flint::Fmpz> defining_coefficients;
    slong defining_degree = -1;
    ulong cyclotomic_index = 0;
    ulong real_cyclotomic_index = 0;
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
// field's defining polynomial modulo p (Dedekind-Kummer), or from its
// closed form when that polynomial is exactly Phi_n or psi_n, without
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
    // called with.  Degrees from a distinct-degree factorization are kept in
    // the order data (OrderAccess::defining_polynomial_residue_degrees,
    // order/order_internal.hpp), where every route object and every
    // Belabas-Friedman evaluation over the same order data reads them.
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
// Theorem 7, below 1/4 (or below 1/q for a log-radius denominator q > 4).
// That route also records its log-radius denominator q in
// `theorem7_log_radius_denominator` and, in `theorem7_prime_sum`, the
// Theorem 7 prime-ideal sum at `cutoff` (sum over prime ideals of norm
// below `cutoff` of log N / (N^sigma - 1)); the denominator is zero on the
// default route.
struct ZetaBfRouteAudit {
    flint::Arb error_bound;
    ulong cutoff = 0;
    ulong max_cutoff = 0;
    slong requested_precision = 0;
    slong work_precision = 0;
    flint::Arb theorem7_prime_sum;
    ulong theorem7_log_radius_denominator = 0;
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
// defining polynomial of `order` modulo p from the batched word pass,
// `per_prime_counts` the same number from the per-prime path (the word
// kernel above degree 6 where it applies, then the nmod path), and
// `nmod_counts` from the per-prime nmod path alone; each is -1 where that
// path does not apply to p.  All three run with the cyclotomic closed form
// turned off, so a cyclotomic defining polynomial still exercises the word
// kernel and the nmod path.  Returns false when `order` is not a maximal
// order with a basis.
class ZetaBfTestEntry {
    friend class TestHookAccess;
    static bool linear_factor_counts(std::vector<ulong>& primes,
                                     std::vector<slong>& batch_counts,
                                     std::vector<slong>& per_prime_counts,
                                     std::vector<slong>& nmod_counts,
                                     const Order& order,
                                     ulong first,
                                     ulong end) noexcept;
    // The linear factor counts the Belabas-Friedman prime loop uses, with
    // the cyclotomic closed forms, against the word kernel and the nmod path
    // without them.  `cyclotomic_index` receives the context's Phi_n index
    // (0 when the defining polynomial is not exactly a cyclotomic
    // polynomial) and `real_cyclotomic_index` its psi_n index (0 when it is
    // not exactly the real-subfield polynomial psi_n).
    // `batch_size` receives the number of primes the word-kernel batch
    // covers for the defining polynomial (0 for Phi_n and psi_n).
    // For each prime p in [first, end), in increasing order, `primes`
    // receives p, `counts` the count of the loop's path (batch, then
    // per-prime), `word_counts` the word kernel's batched count and
    // `nmod_counts` the per-prime nmod count without the word kernel, each
    // -1 where that path does
    // not apply to p.  Returns false when `order` is not a maximal order
    // with a basis.
    static bool cyclotomic_linear_factor_counts(
            ulong& cyclotomic_index,
            ulong& real_cyclotomic_index,
            std::size_t& batch_size,
            std::vector<ulong>& primes,
            std::vector<slong>& counts,
            std::vector<slong>& word_counts,
            std::vector<slong>& nmod_counts,
            const Order& order,
            ulong first,
            ulong end) noexcept;
    // The residue degrees the defining-polynomial route gives, with the
    // cyclotomic closed forms, against the nmod distinct-degree path
    // without them.  `cyclotomic_index` and `real_cyclotomic_index` are as
    // for cyclotomic_linear_factor_counts.  For each prime p in
    // [first, end), in increasing order, `primes` receives p, `degrees` the
    // route's degrees and `nmod_degrees` the nmod path's, each empty where
    // the route does not apply to p.  Returns false when `order` is not a
    // maximal order with a basis.
    static bool cyclotomic_residue_degrees(
            ulong& cyclotomic_index,
            ulong& real_cyclotomic_index,
            std::vector<ulong>& primes,
            std::vector<std::vector<slong>>& degrees,
            std::vector<std::vector<slong>>& nmod_degrees,
            const Order& order,
            ulong first,
            ulong end) noexcept;
    // The word kernel's unreduced square of x^k mod f, with symmetric
    // squaring in `symmetric` and the full i, j loop in `full`, for
    // `power` of 2 to 20 entries, each at most the kernel's prime limit for
    // that degree minus two.  Returns false otherwise.
    static bool word_square(std::vector<ulong>& symmetric,
                            std::vector<ulong>& full,
                            const std::vector<ulong>& power) noexcept;
    // x^p mod f and the number of roots of f modulo p from the word kernel,
    // with symmetric and with full squaring, for a monic `f` (low to high)
    // of degree 2 to 20 with coefficients reduced modulo the prime p below
    // the kernel's prime limit for the degree.  Returns false otherwise.
    static bool word_frobenius(std::vector<ulong>& symmetric,
                               std::vector<ulong>& full,
                               slong& symmetric_count,
                               slong& full_count,
                               const std::vector<ulong>& f,
                               ulong p) noexcept;
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

// The log-radius denominator of the coarse grh index-one value: its log
// residue has radius below 1/4.
inline constexpr ulong kZetaBfIndexOneLogRadiusDenominator = 4;

// The exclusive limit on the log-radius denominator q of a grh index-one
// value: its Theorem 7 tail target is 1/q - 2^-20, which is not positive
// for q >= 2^20, so no cutoff can reach it.
inline constexpr ulong kZetaBfIndexOneLogRadiusDenominatorLimit = UWORD(1)
                                                                  << 20;

// The Belabas-Friedman hR for the grh analytic index-one test: the log
// residue f_K(X) at the least cutoff X (a multiple of 9 up to the default
// route's maximum) whose Belabas-Friedman 2015 Theorem 7 bound, with
// sigma = 1 + 1/sqrt(log D), is below 1/q - 2^-20, q =
// `log_radius_denominator`, so that the log residue has radius below 1/q.
// With q = 4 that is enough for the index-one test (Belabas-Friedman
// 2015, p. 1, footnote 1: below (1/2) log 2); a larger q narrows the
// value for the grh index bound (analytic_index_bound_for_validation), and
// a q below 4 fails.  Valid under GRH for zeta_K
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
        ZetaBfResidueDegreeCache* residue_degree_cache = nullptr,
        ulong log_radius_denominator =
                kZetaBfIndexOneLogRadiusDenominator) noexcept;

// True only when it is proven that the grh index-one route with
// log-radius denominator `log_radius_denominator`
// (zeta_class_regulator_product_index_one_with_diagnostics) fails because
// no cutoff up to the default maximum has a Theorem 7 bound below its
// tail target 1/q - 2^-20.  `index_one_audit` is the audit of a successful
// call of that route on the same order with a smaller denominator (the
// coarse value, q = 4).  Below its cutoff X_c no bound was below the
// larger target of that call, and the scan recomputes those bounds
// identically; from X_c to the maximum X_m the bound is at least
//   2.324 (2 sigma - 1) / (sqrt(X_m) log(3 X_m)) (1 + 3.88 / log(X_m/9))
//       * (delta_0 - 2 (S(X_c) + T)),
// where T bounds the sum over prime ideals of norm in [X_c, X_m) with at
// most n ideals of norm p and floor(n/f) of norm p^f above each prime p,
// the primes taken in short blocks at the value of the block start.  The
// test evaluates that bound in ball arithmetic and compares its lower
// endpoint with the target's upper endpoint.  False means only "not
// proven": the route may still fail.
//
// The proof below X_c is not a stand-alone certificate: it assumes the
// Theorem 7 scan reproduces its balls bit-exactly.  That holds only because
// the scan adds identical balls whatever the residue-degree cache holds, so
// both residue-degree paths (cached and uncached) must issue the same
// additions.  A change to the scan or to bf_theorem7_add_prime that breaks
// this also breaks the skip.
bool zeta_bf_index_one_target_unreachable(
        const Order& order,
        const ZetaBfRouteAudit& index_one_audit,
        ulong log_radius_denominator,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

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

// The values chi(n) = (D / n) of `quadratic_character` for n = 1, 2, ...,
// filled on demand by a segmented sieve and kept, so that both sums of the
// approximate functional equation read them from one table.  The Kronecker
// symbol (D / n) is completely multiplicative in n > 0 (Cohen, "A Course in
// Computational Algebraic Number Theory", GTM 138, Section 1.4.2:
// Definition 1.4.8 of the symbol and Theorem 1.4.9(2)), so a
// composite n with smallest prime factor p has chi(n) = chi(p) chi(n / p);
// `quadratic_character` is called only for n = 1 and the primes n, and the
// table returns exactly its values.  A value is two bits, so the table
// allocates at most max_entries / 4 + 1 bytes (about 16 MiB by default;
// while it grows, a reallocation briefly holds the previous buffer too) plus
// one segment of scratch; past max_entries, `value` calls
// `quadratic_character` directly.
// `discriminant` must outlive the table.
class QuadraticCharacterTable {
public:
    static constexpr ulong default_max_entries = ulong(1) << 26;
    static constexpr ulong default_segment_length = ulong(1) << 16;

    explicit QuadraticCharacterTable(
            flint::FmpzConstRef discriminant,
            ulong max_entries = default_max_entries,
            ulong segment_length = default_segment_length) noexcept;

    QuadraticCharacterTable(const QuadraticCharacterTable&) = delete;
    QuadraticCharacterTable& operator=(const QuadraticCharacterTable&) =
            delete;

    // (D / n); n = 0 is not stored and takes the symbol directly.
    int value(ulong n) noexcept {
        if (n - 1 < filled_) {
            return stored(n);
        }
        return value_slow(n);
    }

    // Values held: chi(n) is stored for 1 <= n <= filled().
    ulong filled() const noexcept {
        return filled_;
    }

    // Bytes allocated for the stored values.
    ulong allocated_bytes() const noexcept {
        return packed_.capacity();
    }

    // Calls of `quadratic_character` made by the table.
    ulong symbol_calls() const noexcept {
        return symbol_calls_;
    }

private:
    int stored(ulong n) const noexcept {
        const unsigned code =
                (static_cast<unsigned>(packed_[n >> 2]) >> ((n & 3) * 2)) & 3u;
        return code == 2 ? -1 : static_cast<int>(code);
    }

    void store(ulong n, int chi) noexcept;
    int value_slow(ulong n) noexcept;
    void extend() noexcept;

    flint::FmpzConstRef discriminant_;
    ulong max_entries_;
    ulong segment_length_;
    ulong filled_ = 0;
    ulong symbol_calls_ = 0;
    std::vector<unsigned char> packed_;
    std::vector<ulong> small_primes_;
    ulong small_primes_through_ = 1;
    std::vector<ulong> factor_;
};

// A fault forced on the blocks the stepper of `quadratic_dirichlet_l1`
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

// Test controls and counters of the Taylor stepping of the incomplete gamma
// values Gamma(a, z) in `quadratic_dirichlet_l1`, and of its tail-bound
// skip.  The counters are indexed
// by shift = 2a: the sums use a = 1/2 and a = 0 for real D, a = 1 and
// a = 1/2 for imaginary D.
struct QuadraticAfeStepAudit {
    // False turns stepping off, so that every value is evaluated directly,
    // as the evaluation did before stepping was added.
    bool stepping = true;
    // A fault forced on every block planned in a sum once `fault_after`
    // blocks of that sum have started (0: on every block).
    // `blocks` is cumulative over the audit object's lifetime, so use a
    // fresh audit for each `quadratic_dirichlet_l1` call.
    QuadraticAfeStepFault fault = QuadraticAfeStepFault::none;
    ulong fault_after = 0;
    // Blocks started, and values taken from a block (an anchor carried
    // over from the previous block included).
    ulong blocks[3] = {0, 0, 0};
    ulong stepped_values[3] = {0, 0, 0};
    // Planned blocks whose start was rejected, by the check that failed:
    // R = 2^r below z_n, v_max below R, a finite truncation bound.
    ulong rejected_radius[3] = {0, 0, 0};
    ulong rejected_step[3] = {0, 0, 0};
    ulong rejected_truncation[3] = {0, 0, 0};
    // False evaluates the tail bound at every index of a sum, as the
    // evaluation did before the tail-bound skip; true skips the indices that
    // a lower bound of the tail bound certifies to come before the stop.
    bool tail_skip = true;
    // The index at which each sum stopped (its first index with a tail bound
    // below the tolerance), and the number of tail bounds it evaluated.
    ulong stop_index[3] = {0, 0, 0};
    ulong tail_bounds[3] = {0, 0, 0};
    // The largest index the tail-bound skip passed over, and the number of
    // lower-bound evaluations its search made (both 0 with tail_skip false).
    ulong skip_index[3] = {0, 0, 0};
    ulong skip_evals[3] = {0, 0, 0};
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

struct QuadraticAfeGammaStepCount {
    ulong blocks = 0;
    ulong stepped_values = 0;
    // Planned blocks whose start was rejected, by the check that failed:
    // R = 2^r below z_n, v_max below R, a finite truncation bound.
    ulong rejected_radius = 0;
    ulong rejected_step = 0;
    ulong rejected_truncation = 0;
};

struct QuadraticAfeGammaStepBlock {
    ulong block = 0;
    slong terms = 0;
    slong radius_exponent = 0;
};

// Test entries of the zeta sources; only detail::TestHookAccess
// (test/test_hook_access.hpp) can call them.
class ZetaTestEntry {
    friend class TestHookAccess;

    // Test entry: the Belabas-Friedman 2015 Theorem 7 bound used by
    // `zeta_class_regulator_product_with_diagnostics`, at
    // `cutoff` (a multiple of 9, at least 18), with the prime-ideal sum over
    // the norms below `cutoff`.
    static bool bf_theorem7_error_bound(flint::ArbRef out,
                                        const Order& order,
                                        ulong cutoff,
                                        slong precision) noexcept;

    // Test entry: `grh_factor_base_bound_with_diagnostics` with `stop_bound`
    // in place of Bach's bound, both as the stop bound of the BDF search and
    // as the bound the BDF bound is compared with.  A stop bound below the
    // BDF bound reaches the case where the search stops (`stopped`) and the
    // stop bound is used without the note of a BDF failure; with Bach's
    // actual bound, real fields reach it only at enormous discriminants.
    static bool grh_factor_base_bound_with_stop(
            flint::FmpzRef out,
            const Order& order,
            const flint::Fmpz& stop_bound,
            const DiagnosticsContext* diagnostics,
            bool* bach_selected) noexcept;

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
    static bool quadratic_afe_gamma_stepped(
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

    // Test entry for one block of that stepper: plans and starts a block
    // anchored at n with anchor tolerance 2^target and coefficients at `prec`,
    // the anchor evaluated directly at `prec`, and writes its length B, its
    // number of Taylor terms M, the exponent r of its Cauchy radius R = 2^r and
    // its truncation bound (Johansson 2021, Lemma 5, eq. 27) to `truncation`.
    // Returns false on invalid arguments (shift > 2, q or n zero, `prec` not
    // positive) or when no block starts at n.
    static bool quadratic_afe_gamma_step_block(QuadraticAfeGammaStepBlock& block,
                                               mag_struct* truncation,
                                               ulong shift,
                                               ulong q,
                                               ulong n,
                                               double target,
                                               slong prec) noexcept;
};

}  // namespace silex::detail
