#pragma once

#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/prime_ideal.hpp>

#include <vector>

namespace silex::detail {

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

// Residue degrees of the primes above p in `order`, in the order
// decompose_prime(primes, order, p) publishes them, with the same success.
// For an index divisor of a maximal non-equation order it reads them from the
// order's index-divisor decomposition cache (OrderAccess::
// index_divisor_decomposition) without rebuilding the stored ideals, and on a
// miss for an unramified p stores the decomposition decompose_prime would
// store, computing the quotient O/pO once for the bounded and recursive
// splits.  Every other case calls decompose_prime.  Leaves `out` unchanged on
// failure.  For callers that need only the decomposition type, such as the
// Bach-Friedman Euler factors in zeta.cpp.
bool decompose_prime_residue_degrees(std::vector<slong>& out,
                                     const Order& order,
                                     flint::FmpzConstRef p) noexcept;

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

    // Test hook: runs the exact power-containment valuation fallback (the
    // largest i <= `bound` with `element` in P^i) directly, bypassing the
    // norm, Hensel and coordinate-divisibility routes.  For tests only.
    static bool valuation_by_power_containment_for_testing(
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
    static bool index_divisor_prime_ideals_for_testing(
            std::vector<Ideal>& out,
            const Order& order,
            flint::FmpzConstRef p,
            bool recursive_split) noexcept;

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

    // `norm_vp` must be the exact p-adic valuation of the absolute norm of
    // `element`; the result is wrong otherwise.
    static bool prime_ideal_valuation_with_norm_vp(
            slong& out,
            const PrimeIdeal& prime,
            const OrderElement& element,
            slong norm_vp,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;
};

}  // namespace silex::detail
