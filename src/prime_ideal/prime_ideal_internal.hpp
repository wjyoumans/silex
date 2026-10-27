#pragma once

#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/prime_ideal.hpp>

#include <atomic>
#include <vector>

namespace silex::detail {

class TestHookAccess;

// Switch that only detail::TestHookAccess (test/test_hook_access.hpp) can
// set.  Off by default and never set by production code; one relaxed atomic
// load per ideal build of unramified_prime_data_from_kernel.
class PrimeIdealTestSeams {
public:
    // True while armed: the next ideal builds of an unramified index-divisor
    // split (ideal_from_mod_rows for a split prime, or pO for an inert p)
    // fail as an allocation failure would.
    static bool ideal_build_fails() noexcept {
        if (!ideal_build_failure_.load(std::memory_order_relaxed)) {
            return false;
        }
        ideal_build_failures_taken_.fetch_add(1, std::memory_order_relaxed);
        return true;
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<bool> ideal_build_failure_{false};
    inline static std::atomic<slong> ideal_build_failures_taken_{0};
};

enum class RetainedQuadraticPrimeKind {
    inert,
    ramified,
    split,
};

// True when prime ideals built now over `order` store their residue
// polynomial in the integral generator omega of a maximal quadratic-backend
// order (basis [1, omega]).  Otherwise the residue polynomial is in the
// NumberField generator alpha.  Every construction route selects its
// convention with this predicate and records the result in the prime; code
// reading an existing prime must use
// PrimeIdealAccess::residue_uses_integral_generator instead, because the
// order's maximality record can change after the prime was built.  See
// PrimeIdeal::residue_polynomial.
//
// Source trace: reference `base2.c:modprinit` builds the residue map from
// the prime's own stored data (for p not dividing the index, the polynomial
// generator and its factor mod p); `quadgen`/`quadpoly` supply the integral
// generator omega used by the direct maximal-quadratic route.
bool residue_polynomial_uses_quadratic_integral_generator(
        const Order& order) noexcept;

// Writes `element` as numerator(v) / denominator, where v is the variable of
// `prime`'s stored residue polynomial (alpha, or omega as recorded in the
// prime when it was built) and numerator has integer coefficients.  On the
// alpha convention the element's alpha-polynomial is used directly: this is
// the image of its order coordinates under the order basis matrix, not the
// coordinates themselves.  Fails, leaving outputs unspecified, when p divides
// the denominator or, on the alpha convention, when p divides the order basis
// denominator.  The alpha convention requires alpha integral (monic integral
// defining polynomial).  The equation and maximal orders on which residue
// data is produced then contain Z[alpha], and for them the failure is exactly
// p | [O : Z[alpha]].  Reducing numerator modulo (p, residue polynomial)
// and multiplying by denominator^{-1} mod p gives the residue class of
// `element` modulo P.
bool residue_variable_numerator(flint::FmpzPoly& numerator,
                                flint::Fmpz& denominator,
                                const PrimeIdeal& prime,
                                const Element& element) noexcept;

// Closed-form factorization of a monic quadratic x^2 + b x + c modulo a
// prime p that fits a word (Cohen, GTM 138, Proposition 5.1.4 for the
// splitting type; Algorithm 1.5.1 for the square root, via FLINT's
// n_sqrtmod).  `reduced` is the polynomial modulo p = the modulus of
// `context`.  For odd p, with D = b^2 - 4c mod p: D = 0 gives (x - r)^2 with
// r = -b/2; a nonresidue D gives the irreducible polynomial itself; a nonzero
// square D = s^2 gives x - r1 and x - r2 with r = (-b +- s)/2.  For p = 2 the
// roots are read off f(0) and f(1).  Canonical order: a split prime's two
// linear factors are listed by increasing root r1 < r2 in [0, p).  `out` is
// overwritten with monic factors and exponents in the layout of
// fmpz_mod_poly_factor (whose order for split primes is not specified).
// Returns false, leaving `out` unspecified, when p is not a word, the
// polynomial is not monic of degree 2, or the computed square root does not
// check.  p must be prime; that is not checked here.
bool monic_quadratic_factor_word(flint::FmpzModPolyFactor& out,
                                 flint::FmpzModPolyConstRef reduced,
                                 const flint::FmpzModCtx& context) noexcept;

// fmpz_mod_poly_factor(out, reduced, context), computed with word-size
// Frobenius powers, gcds and divisions where FLINT would take its
// Cantor-Zassenhaus route (5 * bits(p) + deg <= 75, so p < 2^14); the
// equal-degree splitting stays FLINT's fmpz_mod_poly_factor_equal_deg.  The
// factors, exponents and their order are those of fmpz_mod_poly_factor (see
// the definition).  p must be prime; that is not checked here.
void factor_mod_prime(flint::FmpzModPolyFactor& out,
                      flint::FmpzModPolyConstRef reduced,
                      const flint::FmpzModCtx& context) noexcept;

// The factorization of `reduced` modulo the prime modulus of `context`:
// monic_quadratic_factor_word when it applies, factor_mod_prime
// otherwise.  Every route of decompose_prime that factors a quadratic
// modulo p (omega's minimal polynomial, or a monic quadratic defining
// polynomial on the Dedekind-Kummer routes) factors through it, so a split
// word-size prime's two primes come in the same canonical order on each
// route; larger primes keep FLINT's order.
void factor_with_quadratic_closed_form(
        flint::FmpzModPolyFactor& out,
        flint::FmpzModPolyConstRef reduced,
        const flint::FmpzModCtx& context) noexcept;

class MaximalQuadraticPrimeAccess {
public:
    static bool set_from_integral_generator_factor(
            PrimeIdeal& out,
            const Order& order,
            flint::FmpzConstRef p,
            flint::FmpzModPolyConstRef factor,
            slong ramification_index,
            const Element& integral_generator,
            const flint::FmpzModCtx& context) noexcept;

    static bool set_first_degree_one_prime(
            PrimeIdeal& out,
            RetainedQuadraticPrimeKind& kind,
            const Order& order,
            flint::FmpzConstRef p,
            const DiagnosticsContext* diagnostics) noexcept;
};

// The route decompose_prime(primes, order, p, max_residue_degree) takes.
// One classifier in prime_ideal.cpp chooses it, and decompose_prime,
// PrimeIdealAccess::decompose_prime_containing and
// decompose_prime_residue_degrees all route by it.
enum class DecomposePrimeRoute {
    // Returns false without decomposing: no basis or field, a negative
    // bound, p not prime, a defining polynomial that is not monic integral,
    // or a non-equation order not known to be maximal.
    refused,
    // A degree-one order: (p) itself.
    degree_one,
    // A maximal quadratic-backend order with basis [1, omega]: the factors
    // of the minimal polynomial of omega modulo p (closed form for word-size
    // p, see monic_quadratic_factor_word).
    quadratic_integral_generator,
    // max_residue_degree == 1, an equation or maximal order, and the
    // defining polynomial f squarefree modulo p: the roots of f modulo p.
    squarefree_roots,
    // A maximal non-equation order with p prime to [O : Z[theta]]: the
    // factors of f modulo p (Dedekind-Kummer; closed form for a quadratic
    // f and word-size p, see factor_with_quadratic_closed_form).
    nonindex_factor,
    // A maximal non-equation order with p dividing [O : Z[theta]].
    index_divisor,
    // The equation order: the factors of f modulo p (Dedekind-Kummer;
    // closed form for a quadratic f and word-size p).
    equation_factor,
};

// What decompose_prime_residue_degrees does on an index-divisor cache miss
// for a p not dividing the discriminant, when the bounded or recursive split
// of O/pO gives the decomposition.
enum class IndexDivisorMiss : unsigned char {
    // Build each prime's ideal and store the decomposition in the order's
    // index-divisor cache, as decompose_prime does, so a later
    // decompose_prime or decompose_prime_residue_degrees of p reads it back.
    store_decomposition,
    // Read only the residue degrees, in the same order: no ideal is built
    // and nothing is stored.  When the characteristic polynomial of one of
    // the first two order basis elements is squarefree modulo a word-size p
    // and its factors there all have one degree f, the degrees are read off
    // that factorization (Dedekind-Kummer) without forming O/pO; otherwise
    // they are read off the split.  For a caller that keeps the
    // degrees of p itself and does not expect decompose_prime of p to
    // follow; a later decompose_prime of p computes the same decomposition
    // from scratch.
    degrees_only,
};

// Residue degrees of the primes above p in `order`, in the order
// decompose_prime(primes, order, p) publishes them, with the same success.
// For an index divisor of a maximal non-equation order it reads them from the
// order's index-divisor decomposition cache (OrderAccess::
// index_divisor_decomposition) without rebuilding the stored ideals, and on a
// miss for an unramified p computes the quotient O/pO once for the bounded
// and recursive splits; `miss` says whether it then builds and stores the
// decomposition decompose_prime would store.  Every other case calls
// decompose_prime.  Leaves `out` unchanged on failure.  For callers that need
// only the decomposition type, such as the Bach-Friedman Euler factors in
// zeta.cpp.  `diagnostics` receives the profile events of the monogenic
// generator search this function runs on that miss (not those of a call it
// delegates to decompose_prime).
//
// With IndexDivisorMiss::degrees_only the result (its success, or the order
// of the degrees) can differ from decompose_prime's only where building an
// ideal of the split fails, since decompose_prime then moves on to its next
// route, or where the split itself would fail on an allocation and the
// basis-element route above already gave the degrees.  That build (ideal_from_mod_rows, or pO for an inert p) fails only
// when an allocation fails: the rows pO plus the kept components always span
// a full-rank ideal of the order (see ideal_from_mod_rows).
bool decompose_prime_residue_degrees(
        std::vector<slong>& out,
        const Order& order,
        flint::FmpzConstRef p,
        IndexDivisorMiss miss = IndexDivisorMiss::store_decomposition,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

// Noninstalled prime-ideal hooks: constructors that take caller-supplied
// data, and read-only access to recorded residue-map state.
class PrimeIdealAccess {
public:
    // The integral ideal of `prime`, without the copy get_ideal makes.  Its
    // HNF determines the prime ideal: two primes of one order are the same
    // prime exactly when these compare equal under Ideal::equal, whatever
    // residue data each construction recorded.
    static const Ideal& ideal(const PrimeIdeal& prime) noexcept {
        return prime.ideal_;
    }

    // Whether `left` and `right`, two primes of one order, are the same
    // prime ideal: their integral ideals compare equal under Ideal::equal,
    // which compares canonical HNFs.  e and f then agree too, provided the
    // order is maximal and the recorded e and f are correct (a prime ideal
    // of a maximal order has well-defined e and f).  Unlike PrimeIdeal::
    // equal, this ignores the residue data each construction recorded (the
    // residue polynomial, its variable convention and whether one was
    // recorded at all).  The one same-prime test shared by the factor-base
    // lookup and the class-group partial-relation pairing.
    static bool same_prime_ideal(const PrimeIdeal& left,
                                 const PrimeIdeal& right) noexcept {
        return ideal(left).equal(ideal(right));
    }

    // The residue-polynomial convention recorded when `prime` was built.
    static bool residue_uses_integral_generator(
            const PrimeIdeal& prime) noexcept;

    // Whether p does not divide the denominator of `prime`'s order basis in
    // the alpha power basis (for an order containing Z[alpha]: p does not
    // divide [O : Z[alpha]]).  Cached in the prime after the first call.
    static bool order_basis_is_p_integral(const PrimeIdeal& prime) noexcept;

    // Builds the degree-one prime (p, theta - root) of an equation or
    // maximal order.  It requires f(root) = 0 mod p and the defining
    // polynomial f squarefree mod p, a sufficient condition for
    // (p, theta - root) to be a degree-one prime with e = 1 (p then does not
    // divide the index of Z[theta] and Dedekind-Kummer applies).  The
    // condition is not necessary; when it fails, the call fails and leaves
    // `out` unchanged even if (p, theta - root) is such a prime.
    static bool set_degree_one_prime_ideal_from_root(
            PrimeIdeal& out,
            const Order& order,
            flint::FmpzConstRef p,
            flint::FmpzConstRef root) noexcept;

    // Builds into `out` (defined over `order`) the prime ideal
    // (p, g(theta)) of the irreducible factor `factor` of the defining
    // polynomial modulo p with multiplicity `ramification_index`, as
    // decompose_prime builds it on its Dedekind-Kummer routes: ideal, residue
    // polynomial g and centered Kummer generator.  Shared by those routes and
    // decompose_prime_containing so both build the same prime.
    static bool set_defining_polynomial_factor_prime(
            PrimeIdeal& out,
            const Order& order,
            flint::FmpzConstRef p,
            const fmpz_mod_poly_struct* factor,
            slong ramification_index,
            const flint::FmpzModCtx& ctx,
            const DiagnosticsContext* diagnostics) noexcept;

    // The prime ideals above the rational prime p that can contain the
    // nonzero element `element` of a maximal order, in the order
    // decompose_prime(.., p, 0) lists them and built as it builds them.
    // When decompose_prime would take a Dedekind-Kummer route (p does not
    // divide [O : Z[theta]], and the order is not on the quadratic
    // integral-generator route), the primes above p are P_i = (p, g_i(theta))
    // for the monic irreducible factors g_i of f modulo p.  Writing
    // element = A(theta) / d with d coprime to p, element lies in P_i exactly
    // when g_i divides A modulo p, because P_i meets Z[theta] in
    // p Z[theta] + g_i(theta) Z[theta] (Cohen, GTM 138, Theorem 4.8.13 and
    // its proof).  Only those P_i are built; every omitted P_i has valuation
    // 0 at `element`.  On every other route the result is the whole
    // decomposition.  Returns false on failure; on the Dedekind-Kummer route
    // only the primes it builds can fail.
    static bool decompose_prime_containing(
            std::vector<PrimeIdeal>& out,
            const OrderElement& element,
            flint::FmpzConstRef p,
            const DiagnosticsContext* diagnostics) noexcept;

    // A key for the prime ideal above p containing the nonzero element
    // `element` of a maximal order, when there is only one and it has
    // residue degree 1.  When decompose_prime_containing(.., element, p)
    // takes its Dedekind-Kummer route and G = gcd(f, A) modulo p has degree
    // 1 (element = A(theta) / d as there), sets `found` and writes to `root`
    // the r in [0, p) with G = lead (x - r).  The factors g_i of f modulo p
    // dividing A are the irreducible factors of G, so x - r is the only one:
    // decompose_prime_containing then returns the single prime
    // P_r = (p, theta - r), built from that factor, and it is the only prime
    // above p containing `element`.  For two elements with roots r and s,
    // P_r = P_s exactly when r = s: if r != s and P_r = P_s = P, then
    // r - s = (theta - s) - (theta - r) lies in P meet Z = pZ, which two
    // distinct residues in [0, p) do not allow.  So comparing roots decides
    // what PrimeIdealAccess::same_prime_ideal decides on the two built
    // primes.  Otherwise `found` is false (any other route, or G of another
    // degree).  Returns false on failure.
    static bool single_containing_linear_prime_root(
            bool& found,
            flint::FmpzRef root,
            const OrderElement& element,
            flint::FmpzConstRef p) noexcept;

    // `norm_vp` must be the exact p-adic valuation of the absolute norm of
    // `element`; the result is wrong otherwise.
    static bool prime_ideal_valuation_with_norm_vp(
            slong& out,
            const PrimeIdeal& prime,
            const OrderElement& element,
            slong norm_vp,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;

private:
    friend class TestHookAccess;

    // Test hook: values `element` at `prime` by coordinate divisibility
    // with the anti-uniformizer taken from the HNF of P, the route used for
    // primes without a Kummer quotient map, for every prime of a maximal
    // order (Kummer primes included), so tests can compare it with the
    // power-containment fallback.  For tests only.
    static bool anti_uniformizer_valuation(
            slong& out,
            const PrimeIdeal& prime,
            const OrderElement& element,
            slong bound) noexcept;

    // Test hook: splits O/(p, I_p) for a prime p dividing the index of the
    // maximal order `order` and writes its prime ideals in the order
    // decompose_prime publishes them.  It skips the separating-element search
    // and uses the bounded idempotent enumeration (`recursive_split` false)
    // or the recursive split (true), so tests can compare the two on one
    // algebra.  Returns false if the chosen method fails.  For tests only.
    static bool index_divisor_prime_ideals(
            std::vector<Ideal>& out,
            const Order& order,
            flint::FmpzConstRef p,
            bool recursive_split) noexcept;

    // Test hook: the route decompose_prime(.., order, p, max_residue_degree)
    // takes, from the classifier decompose_prime uses.  For tests only.
    static DecomposePrimeRoute decompose_prime_route_for_max_degree(
            const Order& order,
            flint::FmpzConstRef p,
            slong max_residue_degree) noexcept;

};

}  // namespace silex::detail
