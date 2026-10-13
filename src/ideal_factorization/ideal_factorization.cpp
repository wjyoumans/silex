#include <silex/ideal_factorization.hpp>

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_element.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mat.h>
#include <flint/longlong.h>
#include <flint/ulong_extras.h>

#include <limits>
#include <utility>
#include <vector>

#include "ideal_factorization_internal.hpp"
#include "../prime_ideal/prime_ideal_internal.hpp"

namespace silex {
namespace {

// Profiling wrapper so the valuation callers appear as separate scopes.
bool scoped_prime_ideal_valuation(
        [[maybe_unused]] const char* label,
        slong& valuation,
        const PrimeIdeal& prime,
        const OrderElement& element,
        slong norm_valuation,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal, label);
    return detail::PrimeIdealAccess::prime_ideal_valuation_with_norm_vp(
            valuation, prime, element, norm_valuation, diagnostics);
}

using IdealFactor = FactorPower<PrimeIdeal>;

bool copy_append(std::vector<IdealFactor>& factors,
                 const Order& parent,
                 const PrimeIdeal& prime,
                 slong exponent) noexcept {
    if (exponent <= 0 || !same_order_parent(prime.parent(), &parent) ||
        !prime.has_prime_data()) {
        return false;
    }

    IdealFactor entry{PrimeIdeal(parent), exponent};
    if (!entry.factor.is_defined() || !entry.factor.set(prime)) {
        return false;
    }
    factors.emplace_back(std::move(entry));
    return true;
}

bool reconstruct_factors(Ideal& out,
                         const Order& parent,
                         const std::vector<IdealFactor>& factors) noexcept {
    if (!same_order_parent(out.parent(), &parent)) {
        return false;
    }

    Ideal accumulator(parent);
    if (!accumulator.is_defined() || !accumulator.one()) {
        return false;
    }

    Ideal base(parent);
    Ideal product(parent);
    Ideal square(parent);
    if (!base.is_defined() || !product.is_defined() || !square.is_defined()) {
        return false;
    }

    for (const IdealFactor& factor : factors) {
        if (!same_order_parent(factor.factor.parent(), &parent) ||
            factor.exponent <= 0 || !factor.factor.get_ideal(base)) {
            return false;
        }

        slong exponent = factor.exponent;
        while (exponent > 0) {
            if ((exponent & 1) != 0) {
                if (!product.multiply(accumulator, base)) {
                    return false;
                }
                accumulator.swap(product);
            }

            exponent >>= 1;
            if (exponent > 0) {
                if (!square.multiply(base, base)) {
                    return false;
                }
                base.swap(square);
            }
        }
    }

    out.swap(accumulator);
    return true;
}

bool factor_over_base_reconstruction_requested(
        const DiagnosticsContext* diagnostics) noexcept {
#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    return debug_check_enabled(diagnostics, DiagnosticsModule::ideal,
                               DebugLevel::expensive);
#else
    (void)diagnostics;
    return false;
#endif
}

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
bool factor_over_base_reconstruction_matches(
        const Ideal& ideal,
        const Order& parent,
        const std::vector<IdealFactor>& factors) noexcept {
    Ideal reconstructed(parent);
    return reconstruct_factors(reconstructed, parent, factors) &&
           reconstructed.equal(ideal);
}

bool required_prime_factorization_matches_reference(
        bool direct_ok,
        bool direct_matches,
        const Ideal& ideal,
        const FactorBase& base,
        const PrimeIdeal& required_prime) noexcept {
    const Order* parent = ideal.parent();
    if (parent == nullptr || !same_order_parent(base.parent(), parent) ||
        !same_order_parent(required_prime.parent(), parent)) {
        return !direct_ok;
    }

    IdealFactorization factorization(*parent);
    PrimeIdeal factor(*parent);
    if (!factorization.is_defined() || !factor.is_defined() ||
        !factorization.factor(ideal)) {
        return !direct_ok;
    }

    flint::Fmpz required_p;
    flint::Fmpz factor_p;
    if (!required_prime.rational_prime(flint::FmpzRef(required_p))) {
        return !direct_ok;
    }
    bool reference_matches = true;
    slong required_exponent = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = 0;
        if (!factorization.prime(factor, i) ||
            !factorization.exponent(exponent, i)) {
            return false;
        }
        // Matched as ideals, as base.contains matches: the factorization
        // may record residue data other than required_prime's for the same
        // prime.
        if (detail::PrimeIdealAccess::ideal(factor).equal(
                    detail::PrimeIdealAccess::ideal(required_prime))) {
            required_exponent += exponent;
        } else if (!base.contains(factor)) {
            bool admissible = false;
            if (exponent <= 0 ||
                !factor.rational_prime(flint::FmpzRef(factor_p)) ||
                !detail::required_prime_witness_admits_lacking_prime(
                        admissible, base, flint::FmpzConstRef(factor_p),
                        flint::FmpzConstRef(required_p)) ||
                !admissible) {
                reference_matches = false;
            }
        }
    }
    reference_matches = reference_matches && required_exponent == 1;
    return direct_ok && direct_matches == reference_matches;
}
#endif

bool required_prime_norm_support_is_admissible(
        bool& admissible,
        flint::FmpzConstRef norm,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.required_prime_norm_factor");
    admissible = false;
    if (flint::fmpz_sgn(norm) <= 0 ||
        required_prime.residue_degree() <= 0) {
        return false;
    }

    flint::Fmpz required_rational_prime;
    flint::Fmpz rational_prime;
    flint::FmpzFactor factorization;
    if (!required_prime.rational_prime(
                flint::FmpzRef(required_rational_prime))) {
        return false;
    }
    flint::fmpz_factor(flint::FmpzFactorRef(factorization), norm);

    bool required_support_present = false;
    for (slong i = 0;
         i < flint::fmpz_factor_num(
                     flint::FmpzFactorConstRef(factorization));
         ++i) {
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(rational_prime),
                flint::FmpzFactorConstRef(factorization), i);
        if (flint::fmpz_equal(
                    flint::FmpzConstRef(rational_prime),
                    flint::FmpzConstRef(required_rational_prime))) {
            required_support_present =
                    flint::fmpz_factor_exp(
                            flint::FmpzFactorConstRef(factorization), i) >=
                    static_cast<ulong>(required_prime.residue_degree());
            continue;
        }

        slong block_index = -1;
        if (!base.rational_prime_block_index_for_prime(
                    block_index,
                    flint::FmpzConstRef(rational_prime))) {
            return true;
        }
    }

    admissible = required_support_present;
    return true;
}

bool factor_base_order_element_rational_valuation(
        slong& remaining,
        fmpz_mat_struct* exponents,
        const FactorBase& base,
        const OrderElement& element,
        flint::FmpzConstRef rational_prime,
        slong norm_valuation,
        const DiagnosticsContext* diagnostics) noexcept {
    if (norm_valuation < 0 || remaining < 0 || remaining > norm_valuation ||
        (exponents != nullptr &&
         (::fmpz_mat_nrows(exponents) != 1 ||
          ::fmpz_mat_ncols(exponents) != base.length()))) {
        return false;
    }
    if (remaining == 0) {
        return true;
    }

    slong block_index = -1;
    flint::Fmpz block_prime;
    slong block_length = 0;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.rational_valuation_block_lookup");
        if (!base.rational_prime_block_index_for_prime(block_index,
                                                       rational_prime)) {
            return true;
        }
        if (!base.rational_prime_block_data(flint::FmpzRef(block_prime),
                                            block_length, block_index) ||
            !flint::fmpz_equal(flint::FmpzConstRef(block_prime),
                               rational_prime)) {
            return false;
        }
    }

    for (slong offset = 0; offset < block_length && remaining > 0; ++offset) {
        slong index = -1;
        if (!base.rational_prime_block_index(index, block_index, offset)) {
            return false;
        }
        const PrimeIdeal* prime = base.prime_at(index);
        if (prime == nullptr || prime->residue_degree() <= 0) {
            return false;
        }

        slong valuation = -1;
        if (!scoped_prime_ideal_valuation(
                    "ideal_factorization.valuation_caller_prime_block",
                    valuation, *prime, element, norm_valuation, diagnostics) ||
            valuation < 0) {
            return false;
        }
        if (valuation == 0) {
            continue;
        }

        const slong residue_degree = prime->residue_degree();
        if (valuation > remaining / residue_degree) {
            remaining = -1;
            return true;
        }
        if (exponents != nullptr) {
            ::fmpz_set_si(::fmpz_mat_entry(exponents, 0, index), valuation);
        }
        remaining -= valuation * residue_degree;
    }

    return true;
}

bool scalar_factor_over_base_with_required_prime_direct(
        bool& matches,
        flint::FmpzConstRef scalar,
        const FactorBase& base,
        const PrimeIdeal& required_prime) noexcept {
    matches = false;
    const Order* order = base.parent();
    const slong degree = order == nullptr ? 0 : order->degree();
    if (order == nullptr || order->parent() == nullptr || degree <= 0 ||
        !order->is_maximal() || flint::fmpz_is_zero(scalar) ||
        !same_order_parent(required_prime.parent(), order) ||
        !required_prime.has_prime_data() || base.contains(required_prime) ||
        required_prime.ramification_index() <= 0 ||
        required_prime.residue_degree() <= 0 ||
        required_prime.ramification_index() > degree / required_prime.residue_degree()) {
        return false;
    }

    // Source: base4.c:Q_nffactor/prV_e_muls gives v_P(m)=e_P*v_p(m).
    // buch2.c:divide_p_elt/can_factor requires complete norm accounting.
    // Cancel the positive v_q(m) at each q: sum(e_P*f_P) must be degree,
    // except that a q other than the required prime's may leave its one
    // lacking e = 1 ideal (required_prime_witness_admits_lacking_prime).
    if (required_prime.ramification_index() != 1) return true;
    flint::Fmpz required_p, q, absolute;
    if (!required_prime.rational_prime(flint::FmpzRef(required_p))) return false;
    flint::fmpz_abs(absolute, scalar);
    flint::FmpzFactor factors;
    flint::fmpz_factor(flint::FmpzFactorRef(factors), flint::FmpzConstRef(absolute));
    bool required_seen = false;
    for (slong i = 0; i < flint::fmpz_factor_num(flint::FmpzFactorConstRef(factors)); ++i) {
        flint::fmpz_factor_get_fmpz(q, flint::FmpzFactorConstRef(factors), i);
        slong remaining = degree;
        if (flint::fmpz_equal(flint::FmpzConstRef(q), flint::FmpzConstRef(required_p))) {
            required_seen = true;
            if (flint::fmpz_factor_exp(flint::FmpzFactorConstRef(factors), i) != 1)
                return true;
            remaining -= required_prime.residue_degree();
        }
        slong block = -1;
        if (base.rational_prime_block_index_for_prime(block, flint::FmpzConstRef(q))) {
            flint::Fmpz block_p, factor_p;
            slong length = 0;
            if (!base.rational_prime_block_data(block_p, length, block) ||
                !flint::fmpz_equal(flint::FmpzConstRef(block_p), flint::FmpzConstRef(q)))
                return false;
            for (slong offset = 0; offset < length; ++offset) {
                slong index = -1;
                if (!base.rational_prime_block_index(index, block, offset)) return false;
                const PrimeIdeal* prime = base.prime_at(index);
                if (prime == nullptr || !same_order_parent(prime->parent(), order) ||
                    !prime->has_prime_data() || !prime->rational_prime(factor_p) ||
                    !flint::fmpz_equal(flint::FmpzConstRef(factor_p), flint::FmpzConstRef(q)))
                    return false;
                const slong e = prime->ramification_index();
                const slong f = prime->residue_degree();
                if (e <= 0 || f <= 0 || e > remaining / f) return false;
                remaining -= e * f;
            }
        }
        if (remaining != 0) {
            // v_Q(m) = e_Q v_q(m) > 0 at every Q | q, so the prime ideals
            // above q outside the base carry valuation.  Only the single
            // e = 1 ideal above a q other than the required prime's own may.
            bool admissible = false;
            if (!detail::required_prime_witness_admits_lacking_prime(
                        admissible, base, flint::FmpzConstRef(q),
                        flint::FmpzConstRef(required_p)))
                return false;
            if (!admissible) return true;
        }
    }
    matches = required_seen;
    return true;
}

// Retain the element while factoring its norm and evaluating prime-ideal
// valuations, using the source factor-generation and prime-division logic.
bool order_element_factor_over_base_with_required_prime_direct(
        bool& matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics,
        detail::RequiredPrimeDirectAudit* audit) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.required_prime_direct");
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::ideal,
                        "ideal.required_prime_direct.stage_entered");
    matches = false;
    const Order* order = element.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !order->is_maximal() || element.equal_si(0) ||
        !same_order_parent(base.parent(), order) ||
        !same_order_parent(required_prime.parent(), order) ||
        !required_prime.has_prime_data() ||
        required_prime.residue_degree() <= 0 || base.contains(required_prime)) {
        return false;
    }

    Element ambient(*order->parent());
    flint::Fmpq norm;
    flint::Fmpz norm_abs;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.required_prime_direct_norm");
        if (!ambient.is_defined() || !element.get_element(ambient)) {
            return false;
        }
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.required_prime_direct_norm_compute");
        if (!ambient.norm(flint::FmpqRef(norm)) ||
            fmpz_is_one(fmpq_denref(norm.raw())) == 0) {
            return false;
        }
    }
    fmpz_abs(norm_abs.raw(), fmpq_numref(norm.raw()));
    if (fmpz_is_zero(norm_abs.raw()) != 0) {
        return false;
    }

    flint::Fmpz required_rational_prime;
    flint::Fmpz rational_prime;
    flint::FmpzFactor factorization;
    if (!required_prime.rational_prime(
                flint::FmpzRef(required_rational_prime))) {
        return false;
    }
    {
        // Source: buch2.c:can_factor returns 0 when Z_ppo(N, F->prodZ), the
        // part of N prime to every factor-base rational prime, is not 1,
        // before absZ_factor. The required prime's rational prime is removed
        // too. A nonunit cofactor has a rational prime with no base or
        // required prime ideal above it, which the accounting loop below
        // would reject after the full factorization. If a block cannot be
        // read, the prefilter is skipped and the full path decides. Such an
        // element is now a completed non-match even where the full path
        // would have returned false on an internal error at a smaller prime
        // (block mismatch, failed valuation, exponent overflow), so the
        // honesty witness searches take it as final instead of falling back
        // to the full factorization.  It is no witness either way: its norm
        // has a rational prime with no base or required prime ideal above it.
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.required_prime_direct_cofactor");
        flint::Fmpz cofactor;
        flint::fmpz_set(flint::FmpzRef(cofactor),
                        flint::FmpzConstRef(norm_abs));
        bool cofactor_complete =
                fmpz_cmp_ui(required_rational_prime.raw(), 2) >= 0;
        if (cofactor_complete) {
            ::fmpz_remove(cofactor.raw(), cofactor.raw(),
                          required_rational_prime.raw());
        }
        slong block_length = 0;
        const slong block_count = base.rational_prime_block_count();
        for (slong block = 0; cofactor_complete && block < block_count &&
                              fmpz_is_one(cofactor.raw()) == 0;
             ++block) {
            if (!base.rational_prime_block_data(
                        flint::FmpzRef(rational_prime), block_length, block) ||
                fmpz_cmp_ui(rational_prime.raw(), 2) < 0) {
                cofactor_complete = false;
                break;
            }
            ::fmpz_remove(cofactor.raw(), cofactor.raw(),
                          rational_prime.raw());
        }
        if (cofactor_complete && fmpz_is_one(cofactor.raw()) == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::ideal,
                    "ideal.required_prime_direct.reject_cofactor");
            if (audit != nullptr) {
                audit->rejected_by_cofactor = true;
            }
            return true;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.required_prime_direct_fmpz_factor");
        flint::fmpz_factor(flint::FmpzFactorRef(factorization),
                           flint::FmpzConstRef(norm_abs));
    }

    bool required_prime_seen = false;
    for (slong i = 0;
         i < flint::fmpz_factor_num(
                     flint::FmpzFactorConstRef(factorization));
         ++i) {
        const ulong exponent = flint::fmpz_factor_exp(
                flint::FmpzFactorConstRef(factorization), i);
        if (exponent >
            static_cast<ulong>(std::numeric_limits<slong>::max())) {
            return false;
        }
        const slong norm_valuation = static_cast<slong>(exponent);
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(rational_prime),
                flint::FmpzFactorConstRef(factorization), i);

        slong remaining_valuation = norm_valuation;
        if (flint::fmpz_equal(
                    flint::FmpzConstRef(rational_prime),
                    flint::FmpzConstRef(required_rational_prime))) {
            required_prime_seen = true;
            slong required_valuation = -1;
            if (!scoped_prime_ideal_valuation(
                    "ideal_factorization.valuation_caller_required_prime",
                        required_valuation, required_prime, element,
                        norm_valuation, diagnostics)) {
                return false;
            }
            if (required_valuation != 1) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::ideal,
                        "ideal.required_prime_direct.reject_required_valuation");
                return true;
            }
            remaining_valuation -= required_prime.residue_degree();
            if (remaining_valuation < 0) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::ideal,
                        "ideal.required_prime_direct.reject_required_overshoot");
                return true;
            }
        }

        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.required_prime_direct_accounting");
        slong unaccounted = remaining_valuation;
        if (!factor_base_order_element_rational_valuation(
                    unaccounted, nullptr, base, element,
                    flint::FmpzConstRef(rational_prime), norm_valuation,
                    diagnostics)) {
            return false;
        }
        if (unaccounted > 0) {
            // v_q(N(alpha)) = sum over Q | q of f_Q v_Q(alpha), so a positive
            // remainder is carried by the primes above q outside the base.
            // When that is the single e = 1 ideal above a rational prime q
            // other than the required prime's, the principal ideal qO_K puts
            // its class in the subgroup of the base ideals above q, and the
            // witness still places the required class in the base subgroup.
            bool admissible = false;
            if (!detail::required_prime_witness_admits_lacking_prime(
                        admissible, base, flint::FmpzConstRef(rational_prime),
                        flint::FmpzConstRef(required_rational_prime))) {
                return false;
            }
            if (admissible) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::ideal,
                        "ideal.required_prime_direct.lacking_prime");
                unaccounted = 0;
            }
        }
        if (unaccounted != 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::ideal,
                    "ideal.required_prime_direct.reject_accounted");
            return true;
        }
    }

    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::ideal,
            required_prime_seen
                    ? "ideal.required_prime_direct.accepted"
                    : "ideal.required_prime_direct.reject_required_absent");
    matches = required_prime_seen;
    return true;
}

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
bool order_element_required_prime_matches_reference(
        bool direct_ok,
        bool direct_matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* order = element.parent();
    if (order == nullptr) {
        return !direct_ok;
    }

    Ideal principal(*order);
    bool reference_matches = false;
    const bool reference_ok =
            principal.is_defined() && principal.set_principal(element) &&
            detail::ideal_factor_over_base_with_required_prime(
                    reference_matches, principal, base, required_prime,
                    diagnostics);
    return direct_ok == reference_ok &&
           (!direct_ok || direct_matches == reference_matches);
}
#endif

bool ideal_factor_over_base_direct(flint::FmpzMatRef exponents,
                                   bool& smooth,
                                   const Ideal& ideal,
                                   const FactorBase& base,
                                   const DiagnosticsContext* diagnostics,
                                   const PrimeIdeal* required_prime = nullptr) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.factor_over_base_direct");
    const Order* parent = ideal.parent();
    if (parent == nullptr || !same_order_parent(base.parent(), parent) ||
        !ideal.has_hnf() ||
        !parent->is_maximal() || flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length()) {
        return false;
    }
    if (required_prime != nullptr &&
        (!same_order_parent(required_prime->parent(), parent) ||
         !required_prime->has_prime_data() || base.contains(*required_prime))) {
        return false;
    }

    flint::Fmpz norm;
    flint::Fmpz remaining_norm;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.factor_over_base_norm");
        if (!ideal.norm(flint::FmpzRef(norm))) {
            return false;
        }
    }
    fmpz_set(remaining_norm.raw(), norm.raw());
    flint::fmpz_mat_zero(exponents);

    if (required_prime != nullptr) {
        bool norm_support_is_admissible = false;
        if (!required_prime_norm_support_is_admissible(
                    norm_support_is_admissible,
                    flint::FmpzConstRef(norm), base, *required_prime,
                    diagnostics)) {
            return false;
        }
        if (!norm_support_is_admissible) {
            smooth = false;
            return true;
        }
    }

    const bool audit_reconstruction =
            factor_over_base_reconstruction_requested(diagnostics);
    std::vector<IdealFactor> factors;
    if (audit_reconstruction) {
        factors.reserve(static_cast<std::size_t>(base.length()));
    }

    flint::Fmpz required_rational_prime;
    if (required_prime != nullptr) {
        slong valuation = -1;
        if (!required_prime->valuation(valuation, ideal, diagnostics)) {
            return false;
        }
        if (valuation != 1) {
            smooth = false;
            return true;
        }

        flint::Fmpz required_norm;
        const slong residue_degree = required_prime->residue_degree();
        if (residue_degree <= 0 ||
            !required_prime->rational_prime(
                    flint::FmpzRef(required_rational_prime))) {
            return false;
        }
        flint::fmpz_pow_ui(flint::FmpzRef(required_norm),
                           flint::FmpzConstRef(required_rational_prime),
                           static_cast<ulong>(residue_degree));
        if (!flint::fmpz_divisible(flint::FmpzConstRef(remaining_norm),
                                   flint::FmpzConstRef(required_norm))) {
            return false;
        }
        flint::fmpz_divexact(flint::FmpzRef(remaining_norm),
                             flint::FmpzConstRef(remaining_norm),
                             flint::FmpzConstRef(required_norm));
        if (audit_reconstruction &&
            !copy_append(factors, *parent, *required_prime, 1)) {
            return false;
        }
    }

    flint::Fmpz p;
    bool is_smooth = true;
    for (slong block_index = 0;
         block_index < base.rational_prime_block_count() &&
         fmpz_is_one(remaining_norm.raw()) == 0;
         ++block_index) {
        slong length = 0;
        if (!base.rational_prime_block_data(flint::FmpzRef(p), length,
                                            block_index)) {
            return false;
        }

        slong remaining_vp = 0;
        while (fmpz_divisible(remaining_norm.raw(), p.raw()) != 0) {
            fmpz_divexact(remaining_norm.raw(), remaining_norm.raw(), p.raw());
            ++remaining_vp;
        }

        if (remaining_vp == 0) {
            continue;
        }

        for (slong offset = 0; offset < length && remaining_vp > 0; ++offset) {
            slong index = -1;
            if (!base.rational_prime_block_index(index, block_index,
                                                 offset)) {
                return false;
            }
            const PrimeIdeal* prime = base.prime_at(index);
            if (prime == nullptr) {
                return false;
            }

            slong valuation = -1;
            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                                    "ideal.factor_over_base_valuation");
                if (!prime->valuation(valuation, ideal, diagnostics)) {
                    return false;
                }
            }
            if (valuation > 0) {
                fmpz_set_si(flint::fmpz_mat_entry(exponents, 0, index).raw(),
                            valuation);
                if (audit_reconstruction &&
                    !copy_append(factors, *parent, *prime, valuation)) {
                    return false;
                }
                remaining_vp -= valuation * prime->residue_degree();
            }
        }

        if (remaining_vp > 0 && required_prime != nullptr) {
            // A required-prime witness may leave the remainder at the one
            // e = 1 ideal above p that the base lacks, for p other than the
            // required prime's rational prime; see
            // required_prime_witness_admits_lacking_prime.
            PrimeIdeal lacking(*parent);
            bool admissible = false;
            if (!lacking.is_defined() ||
                !detail::required_prime_witness_admits_lacking_prime(
                        admissible, base, flint::FmpzConstRef(p),
                        flint::FmpzConstRef(required_rational_prime),
                        audit_reconstruction ? &lacking : nullptr)) {
                return false;
            }
            if (admissible) {
                if (audit_reconstruction) {
                    slong valuation = -1;
                    if (!lacking.valuation(valuation, ideal, diagnostics) ||
                        valuation <= 0 ||
                        !copy_append(factors, *parent, lacking, valuation)) {
                        return false;
                    }
                }
                remaining_vp = 0;
            }
        }

        if (remaining_vp != 0) {
            is_smooth = false;
            break;
        }
    }

    if (is_smooth && fmpz_is_one(remaining_norm.raw()) == 0) {
        is_smooth = false;
    }

    if (is_smooth && audit_reconstruction) {
        SILEX_DEBUG_CHECK(diagnostics, DiagnosticsModule::ideal,
                          DebugLevel::expensive,
                          "ideal factor-over-base reconstruction",
                          factor_over_base_reconstruction_matches(
                                  ideal, *parent, factors));
    }

    smooth = is_smooth;
    return true;
}

bool checked_add(slong& out, slong left, slong right) noexcept {
    if ((right > 0 && left > std::numeric_limits<slong>::max() - right) ||
        (right < 0 && left < std::numeric_limits<slong>::min() - right)) {
        return false;
    }
    out = left + right;
    return true;
}

bool checked_mul_nonneg(slong& out, slong left, slong right) noexcept {
    if (left < 0 || right < 0 ||
        (left != 0 && right > std::numeric_limits<slong>::max() / left)) {
        return false;
    }
    out = left * right;
    return true;
}

slong strip_rational_prime(flint::Fmpz& value, const flint::Fmpz& p) noexcept {
    slong count = 0;
    while (fmpz_divisible(value.raw(), p.raw()) != 0) {
        fmpz_divexact(value.raw(), value.raw(), p.raw());
        ++count;
    }
    return count;
}

// v_P(N / d) = v_P(N) - e(P/p) v_p(d) for P above p.
bool fractional_prime_valuation(slong& out,
                                const PrimeIdeal& prime,
                                const Ideal& numerator,
                                slong denominator_vp,
                                const DiagnosticsContext* diagnostics) noexcept {
    slong numerator_valuation = -1;
    slong denominator_valuation = 0;
    if (!prime.valuation(numerator_valuation, numerator, diagnostics) ||
        numerator_valuation < 0 || prime.ramification_index() <= 0 ||
        !checked_mul_nonneg(denominator_valuation, prime.ramification_index(),
                            denominator_vp)) {
        return false;
    }
    out = numerator_valuation - denominator_valuation;
    return true;
}

// Factors I = N / d over the base from valuations of I itself, so a base
// with an incomplete rational-prime block still accepts, for example,
// I = P^{-1} over {P}.  Per rational prime p the exact norm identity
// sum_{P | p} f(P) v_P(I) = v_p(N(I)) is checked over the base primes.
// It proves v_P(I) = 0 at the non-base P | p when p does not divide d
// (then every v_P(I) = v_P(N) >= 0) or when the base holds every prime
// above p (sum e(P) f(P) = n).  Otherwise the non-base valuations can
// cancel in the norm, so they are computed from the decomposition of p.
// The stored denominator is minimal (set_integral_den removes
// gcd(content(N), d)), so a prime p | d has some P | p with v_P(I) < 0 and
// a p | d without a base block leaves I nonsmooth.
bool fractional_ideal_factor_over_base_direct(
        flint::FmpzMatRef exponents,
        bool& smooth,
        const Ideal& numerator,
        flint::FmpzConstRef denominator,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.frac_factor_over_base_direct");
    const Order* parent = numerator.parent();
    if (parent == nullptr || !same_order_parent(base.parent(), parent) ||
        !numerator.has_hnf() || !parent->is_maximal() ||
        fmpz_sgn(denominator.raw()) <= 0 ||
        flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length()) {
        return false;
    }
    const slong degree = parent->degree();

    flint::Fmpz remaining_norm;
    flint::Fmpz remaining_denominator;
    if (!numerator.norm(flint::FmpzRef(remaining_norm))) {
        return false;
    }
    fmpz_set(remaining_denominator.raw(), denominator.raw());
    flint::fmpz_mat_zero(exponents);

    flint::Fmpz p;
    bool is_smooth = true;
    for (slong block_index = 0;
         block_index < base.rational_prime_block_count() &&
         (fmpz_is_one(remaining_norm.raw()) == 0 ||
          fmpz_is_one(remaining_denominator.raw()) == 0);
         ++block_index) {
        slong length = 0;
        if (!base.rational_prime_block_data(flint::FmpzRef(p), length,
                                            block_index)) {
            return false;
        }

        const slong numerator_vp = strip_rational_prime(remaining_norm, p);
        const slong denominator_vp =
                strip_rational_prime(remaining_denominator, p);
        if (numerator_vp == 0 && denominator_vp == 0) {
            continue;
        }

        // v_p(N(I)) = v_p(N(N)) - n v_p(d).
        slong denominator_norm_vp = 0;
        if (!checked_mul_nonneg(denominator_norm_vp, degree,
                                denominator_vp)) {
            return false;
        }
        const slong target = numerator_vp - denominator_norm_vp;

        slong accounted = 0;
        slong covered_degree = 0;
        for (slong offset = 0; offset < length; ++offset) {
            slong index = -1;
            if (!base.rational_prime_block_index(index, block_index,
                                                 offset)) {
                return false;
            }
            const PrimeIdeal* prime = base.prime_at(index);
            if (prime == nullptr || prime->residue_degree() <= 0) {
                return false;
            }

            slong valuation = 0;
            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                                    "ideal.frac_factor_over_base_valuation");
                if (!fractional_prime_valuation(valuation, *prime, numerator,
                                                denominator_vp,
                                                diagnostics)) {
                    return false;
                }
            }
            slong local_degree = 0;
            if (!checked_mul_nonneg(local_degree, prime->ramification_index(),
                                    prime->residue_degree()) ||
                !checked_add(covered_degree, covered_degree, local_degree)) {
                return false;
            }
            if (valuation == 0) {
                continue;
            }
            fmpz_set_si(flint::fmpz_mat_entry(exponents, 0, index).raw(),
                        valuation);
            const slong magnitude = valuation < 0 ? -valuation : valuation;
            slong contribution = 0;
            if (!checked_mul_nonneg(contribution, magnitude,
                                    prime->residue_degree()) ||
                !checked_add(accounted, accounted,
                             valuation < 0 ? -contribution : contribution)) {
                return false;
            }
        }

        if (denominator_vp == 0 || covered_degree == degree) {
            if (accounted != target) {
                is_smooth = false;
                break;
            }
            continue;
        }

        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                            "ideal.frac_factor_over_base_incomplete_block");
        PrimeIdealList decomposed;
        if (!decompose_prime(decomposed, *parent, flint::FmpzConstRef(p), 0,
                             diagnostics)) {
            return false;
        }
        for (slong i = 0; i < decomposed.size() && is_smooth; ++i) {
            const PrimeIdeal* prime = decomposed.at(i);
            if (prime == nullptr) {
                return false;
            }
            if (base.contains(*prime)) {
                continue;
            }
            slong valuation = 0;
            if (!fractional_prime_valuation(valuation, *prime, numerator,
                                            denominator_vp, diagnostics)) {
                return false;
            }
            if (valuation != 0) {
                is_smooth = false;
            }
        }
        if (!is_smooth) {
            break;
        }
        if (accounted != target) {
            // Every prime above p was valued; the norm identity must hold.
            return false;
        }
    }

    if (is_smooth && (fmpz_is_one(remaining_norm.raw()) == 0 ||
                      fmpz_is_one(remaining_denominator.raw()) == 0)) {
        is_smooth = false;
    }

    smooth = is_smooth;
    return true;
}

}  // namespace

IdealFactorization::IdealFactorization(const Order& parent) noexcept {
    define(parent);
}

IdealFactorization::~IdealFactorization() noexcept = default;

IdealFactorization::IdealFactorization(IdealFactorization&& other) noexcept {
    swap(other);
}

IdealFactorization& IdealFactorization::operator=(
        IdealFactorization&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void IdealFactorization::swap(IdealFactorization& other) noexcept {
    parent_.swap(other.parent_);
    factors_.swap(other.factors_);
}

void IdealFactorization::clear() noexcept {
    factors_.clear();
    parent_.clear();
}

bool IdealFactorization::define(const Order& parent) noexcept {
    if (!parent.has_basis()) {
        return false;
    }

    clear();
    parent_ = parent;
    return true;
}

bool IdealFactorization::set(const IdealFactorization& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    IdealFactorization copy(other.parent_);
    if (!copy.is_defined()) {
        return false;
    }

    for (const IdealFactor& factor : other.factors_) {
        if (!copy_append(copy.factors_,
                         copy.parent_,
                         factor.factor,
                         factor.exponent)) {
            return false;
        }
    }

    swap(copy);
    return true;
}

bool IdealFactorization::is_defined() const noexcept {
    return parent_.has_basis();
}

const Order* IdealFactorization::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong IdealFactorization::length() const noexcept {
    return is_defined() ? static_cast<slong>(factors_.size()) : 0;
}

bool IdealFactorization::factor(const Ideal& ideal) noexcept {
    const Order* order = ideal.parent();
    if (order == nullptr || !ideal.has_hnf() || !order->is_maximal()) {
        return false;
    }

    IdealFactorization candidate(*order);
    flint::Fmpz norm;
    flint::Fmpz p;
    flint::FmpzFactor norm_factorization;
    if (!candidate.is_defined() || !ideal.norm(flint::FmpzRef(norm))) {
        return false;
    }

    fmpz_factor(norm_factorization.raw(), norm.raw());
    for (slong i = 0; i < norm_factorization.raw()->num; ++i) {
        fmpz_factor_get_fmpz(p.raw(), norm_factorization.raw(), i);

        PrimeIdealList decomposed;
        if (!decompose_prime(decomposed, *order, flint::FmpzConstRef(p))) {
            return false;
        }

        for (slong j = 0; j < decomposed.size(); ++j) {
            const PrimeIdeal* prime = decomposed.at(j);
            slong valuation = -1;
            if (prime == nullptr || !prime->valuation(valuation, ideal)) {
                return false;
            }
            if (valuation > 0 &&
                !copy_append(candidate.factors_, *order, *prime, valuation)) {
                return false;
            }
        }
    }

    Ideal reconstructed(*order);
    if (!candidate.reconstruct(reconstructed) || !reconstructed.equal(ideal)) {
        return false;
    }

    swap(candidate);
    return true;
}

bool IdealFactorization::prime(PrimeIdeal& out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length() ||
        !same_order_parent(out.parent(), &parent_)) {
        return false;
    }
    return out.set(factors_[static_cast<std::size_t>(index)].factor);
}

bool IdealFactorization::exponent(slong& out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length()) {
        return false;
    }
    out = factors_[static_cast<std::size_t>(index)].exponent;
    return true;
}

bool IdealFactorization::reconstruct(Ideal& out) const noexcept {
    if (!is_defined()) {
        return false;
    }
    return reconstruct_factors(out, parent_, factors_);
}

bool ideal_is_smooth(bool& smooth,
                     const Ideal& ideal,
                     const FactorBase& base) noexcept {
    if (ideal.parent() == nullptr ||
        !same_order_parent(ideal.parent(), base.parent())) {
        return false;
    }

    flint::FmpzMat exponents(1, base.length());
    return ideal_factor_over_base_direct(flint::FmpzMatRef(exponents), smooth,
                                         ideal, base, nullptr);
}

namespace detail {

bool required_prime_witness_admits_lacking_prime(
        bool& admissible,
        const FactorBase& base,
        flint::FmpzConstRef q,
        flint::FmpzConstRef required_rational_prime,
        PrimeIdeal* lacking) noexcept {
    admissible = false;
    const Order* order = base.parent();
    const slong degree = order == nullptr ? 0 : order->degree();
    if (order == nullptr || degree <= 0 || !order->is_maximal() ||
        flint::fmpz_cmp_ui(q, 2) < 0 ||
        (lacking != nullptr && !same_order_parent(lacking->parent(), order))) {
        return false;
    }
    // The required prime's own rational prime is excluded: its valuations
    // are accounted for by the required ideal and the base alone.
    if (flint::fmpz_equal(q, required_rational_prime)) {
        return true;
    }
    slong block = -1;
    if (!base.rational_prime_block_index_for_prime(block, q)) {
        return true;
    }

    flint::Fmpz block_p;
    flint::Fmpz factor_p;
    slong length = 0;
    if (!base.rational_prime_block_data(flint::FmpzRef(block_p), length,
                                        block) ||
        !flint::fmpz_equal(flint::FmpzConstRef(block_p), q) || length <= 0) {
        return false;
    }
    slong missing_degree = degree;
    for (slong offset = 0; offset < length; ++offset) {
        slong index = -1;
        if (!base.rational_prime_block_index(index, block, offset)) {
            return false;
        }
        const PrimeIdeal* prime = base.prime_at(index);
        if (prime == nullptr || !same_order_parent(prime->parent(), order) ||
            !prime->has_prime_data() ||
            !prime->rational_prime(flint::FmpzRef(factor_p)) ||
            !flint::fmpz_equal(flint::FmpzConstRef(factor_p), q)) {
            return false;
        }
        const slong e = prime->ramification_index();
        const slong f = prime->residue_degree();
        if (e <= 0 || f <= 0 || e > missing_degree / f) {
            return false;
        }
        missing_degree -= e * f;
    }
    if (missing_degree == 0) {
        return true;
    }
    // sum e_i f_i = n over the primes above q (Cohen, GTM 138, Section
    // 4.8): a remainder of 1 is exactly one prime ideal, with e = f = 1.
    if (missing_degree == 1 && lacking == nullptr) {
        admissible = true;
        return true;
    }

    PrimeIdealList primes;
    if (!decompose_prime(primes, *order, q)) {
        return false;
    }
    const PrimeIdeal* outside = nullptr;
    slong outside_count = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const PrimeIdeal* prime = primes.at(i);
        if (prime == nullptr) {
            return false;
        }
        if (!base.contains(*prime)) {
            outside = prime;
            ++outside_count;
        }
    }
    if (outside_count != 1 || outside->ramification_index() != 1) {
        return true;
    }
    if (lacking != nullptr && !lacking->set(*outside)) {
        return false;
    }
    admissible = true;
    return true;
}

bool order_element_factor_over_base_with_one_large_prime(
        OneLargePrimeFactorStatus& status,
        flint::FmpzMatRef exponents,
        PrimeIdeal& large_prime,
        const OrderElement& element,
        flint::FmpqConstRef norm,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics) noexcept {
    // reference `_factor!` and reference `can_factor` retain an integral element while
    // distributing its known norm over factor-base prime ideals.  Apply that
    // same exact valuation contract while permitting the one non-base prime
    // needed by the existing Silex partial-relation table.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.order_element_factor_one_large_prime");
    status = OneLargePrimeFactorStatus::no_candidate;
    const Order* order = element.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !order->is_maximal() || element.equal_si(0) ||
        !same_order_parent(base.parent(), order) ||
        !same_order_parent(large_prime.parent(), order) ||
        ::fmpz_mat_nrows(exponents.raw()) != 1 ||
        ::fmpz_mat_ncols(exponents.raw()) != base.length() ||
        ::fmpz_is_one(fmpq_denref(norm.raw())) == 0 ||
        ::fmpz_is_zero(fmpq_numref(norm.raw())) != 0) {
        return false;
    }

    flint::Fmpz norm_abs;
    flint::Fmpz rational_prime;
    flint::FmpzFactor norm_factorization;
    flint::FmpzMat candidate_row(1, base.length());
    PrimeIdeal candidate_large_prime(*order);
    if (!candidate_large_prime.is_defined()) {
        return false;
    }
    ::fmpz_abs(norm_abs.raw(), fmpq_numref(norm.raw()));
    flint::fmpz_mat_zero(flint::FmpzMatRef(candidate_row));
    flint::fmpz_factor(flint::FmpzFactorRef(norm_factorization),
                       flint::FmpzConstRef(norm_abs));

    bool have_large_prime = false;
    slong large_exponent = 0;
    const slong factor_count = flint::fmpz_factor_num(
            flint::FmpzFactorConstRef(norm_factorization));
    for (slong i = 0; i < factor_count; ++i) {
        const ulong exponent = flint::fmpz_factor_exp(
                flint::FmpzFactorConstRef(norm_factorization), i);
        if (exponent > static_cast<ulong>(std::numeric_limits<slong>::max())) {
            return false;
        }
        const slong norm_valuation = static_cast<slong>(exponent);
        slong remaining = norm_valuation;
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(rational_prime),
                flint::FmpzFactorConstRef(norm_factorization), i);

        if (!factor_base_order_element_rational_valuation(
                    remaining, candidate_row.raw(), base, element,
                    flint::FmpzConstRef(rational_prime), norm_valuation,
                    diagnostics)) {
            return false;
        }
        if (remaining < 0) {
            return true;
        }
        if (remaining == 0) {
            continue;
        }

        // Only the primes above the rational prime that can contain the
        // element are built: every other one has valuation 0, which the
        // loop below skips.
        std::vector<PrimeIdeal> decomposed;
        if (!detail::PrimeIdealAccess::decompose_prime_containing(
                    decomposed, element, flint::FmpzConstRef(rational_prime),
                    diagnostics)) {
            return false;
        }
        for (const PrimeIdeal& candidate : decomposed) {
            const PrimeIdeal* prime = &candidate;
            if (base.contains(*prime)) {
                continue;
            }

            slong valuation = -1;
            if (!scoped_prime_ideal_valuation(
                    "ideal_factorization.valuation_caller_decomposed_cofactor",
                        valuation, *prime, element, norm_valuation,
                        diagnostics) ||
                valuation < 0 || prime->residue_degree() <= 0) {
                return false;
            }
            if (valuation == 0) {
                continue;
            }

            const slong residue_degree = prime->residue_degree();
            if (valuation > remaining / residue_degree) {
                return true;
            }
            remaining -= valuation * residue_degree;

            // A second sighting of the large prime, compared as an ideal,
            // adds to its exponent; any other prime outside the base ends
            // the search.
            if (!have_large_prime) {
                if (!candidate_large_prime.set(*prime)) {
                    return false;
                }
                have_large_prime = true;
                large_exponent = valuation;
            } else if (detail::PrimeIdealAccess::ideal(candidate_large_prime)
                               .equal(detail::PrimeIdealAccess::ideal(
                                       *prime))) {
                if (valuation >
                    std::numeric_limits<slong>::max() - large_exponent) {
                    return false;
                }
                large_exponent += valuation;
            } else {
                return true;
            }
        }
        if (remaining != 0) {
            return true;
        }
    }

    if (!have_large_prime || large_exponent != 1) {
        return true;
    }
    if (!large_prime.set(candidate_large_prime)) {
        return false;
    }
    flint::fmpz_mat_set(exponents, flint::FmpzMatConstRef(candidate_row));
    status = OneLargePrimeFactorStatus::found;
    return true;
}

bool ideal_factor_over_base_with_required_prime(
        bool& matches,
        const Ideal& ideal,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics) noexcept {
    matches = false;
    if (ideal.parent() == nullptr ||
        !same_order_parent(ideal.parent(), base.parent()) ||
        !same_order_parent(ideal.parent(), required_prime.parent())) {
        return false;
    }

    flint::FmpzMat exponents(1, base.length());
    bool direct_matches = false;
    const bool direct_ok = ideal_factor_over_base_direct(
            flint::FmpzMatRef(exponents), direct_matches, ideal, base,
            diagnostics, &required_prime);

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    SILEX_DEBUG_CHECK(
            diagnostics, DiagnosticsModule::ideal, DebugLevel::expensive,
            "ideal required-prime factor-over-base differential",
            required_prime_factorization_matches_reference(
                    direct_ok, direct_matches, ideal, base, required_prime));
#endif

    if (!direct_ok) {
        return false;
    }
    matches = direct_matches;
    return true;
}

bool scalar_factor_over_base_with_required_prime(
        bool& matches,
        flint::FmpzConstRef scalar,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.scalar_required_prime");
    const bool ok = scalar_factor_over_base_with_required_prime_direct(
            matches, scalar, base, required_prime);
#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    const auto agrees_with_reference = [&]() noexcept {
        if (!ok) return true; // Unsupported evaluations belong to the caller's fallback.
        const Order& order = *base.parent();
        Element ambient(*order.parent());
        OrderElement element(order);
        Ideal principal(order);
        return ambient.set_fmpz(scalar) && element.set_element(ambient) &&
               principal.set_principal(element) &&
               required_prime_factorization_matches_reference(
                       ok, matches, principal, base, required_prime);
    };
    SILEX_DEBUG_CHECK(diagnostics, DiagnosticsModule::ideal, DebugLevel::expensive,
                     "scalar required-prime differential", agrees_with_reference());
#endif
    return ok;
}

bool order_element_factor_over_base_with_required_prime(
        bool& matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics,
        RequiredPrimeDirectAudit* audit) noexcept {
    bool direct_matches = false;
    const bool direct_ok =
            order_element_factor_over_base_with_required_prime_direct(
                    direct_matches, element, base, required_prime, diagnostics,
                    audit);

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    SILEX_DEBUG_CHECK(
            diagnostics, DiagnosticsModule::ideal, DebugLevel::expensive,
            "order-element required-prime factor-over-base differential",
            order_element_required_prime_matches_reference(
                    direct_ok, direct_matches, element, base, required_prime,
                    diagnostics));
#endif

    if (!direct_ok) {
        matches = false;
        return false;
    }
    matches = direct_matches;
    return true;
}

bool required_prime_norm_prefilter_init(
        RequiredPrimeNormPrefilter& prefilter,
        const FactorBase& base,
        const PrimeIdeal& required_prime) noexcept {
    prefilter.ready = false;
    prefilter.degree = 0;
    prefilter.word_table_ready = false;
    // The preconditions of order_element_factor_over_base_with_required_prime_direct
    // for an element of base.parent(); when one fails the direct path reports
    // an evaluation failure, which the prefilter leaves to it.
    const Order* order = base.parent();
    if (order == nullptr || order->parent() == nullptr ||
        !order->is_maximal() ||
        !same_order_parent(required_prime.parent(), order) ||
        !required_prime.has_prime_data() ||
        required_prime.residue_degree() <= 0 ||
        base.contains(required_prime)) {
        return false;
    }
    const slong n = order->degree();
    if (n < 2 || n > 4) {
        return false;
    }

    // The direct cofactor test removes these primes from |N|; it is skipped
    // when one cannot be read, and so is the prefilter.
    flint::Fmpz p;
    slong block_length = 0;
    if (!required_prime.rational_prime(flint::FmpzRef(p)) ||
        fmpz_cmp_ui(p.raw(), 2) < 0) {
        return false;
    }
    ::fmpz_set(prefilter.smooth_product.raw(), p.raw());
    const slong block_count = base.rational_prime_block_count();
    for (slong block = 0; block < block_count; ++block) {
        if (!base.rational_prime_block_data(flint::FmpzRef(p), block_length,
                                            block) ||
            fmpz_cmp_ui(p.raw(), 2) < 0) {
            return false;
        }
        ::fmpz_mul(prefilter.smooth_product.raw(),
                   prefilter.smooth_product.raw(), p.raw());
    }

    // Multiplication table of the order basis w_0, ..., w_(n-1).
    OrderElement left(*order);
    OrderElement right(*order);
    OrderElement product(*order);
    flint::FmpzMat unit(1, n);
    flint::FmpzMat coordinates(1, n);
    flint::FmpzMat table(n * n, n);
    if (!left.is_defined() || !right.is_defined() || !product.is_defined()) {
        return false;
    }
    for (slong i = 0; i < n; ++i) {
        ::fmpz_mat_zero(unit.raw());
        ::fmpz_one(fmpz_mat_entry(unit.raw(), 0, i));
        if (!left.set_coordinates(flint::FmpzMatConstRef(unit))) {
            return false;
        }
        for (slong j = i; j < n; ++j) {
            ::fmpz_mat_zero(unit.raw());
            ::fmpz_one(fmpz_mat_entry(unit.raw(), 0, j));
            if (!right.set_coordinates(flint::FmpzMatConstRef(unit)) ||
                !product.multiply(left, right) ||
                !product.get_coordinates(flint::FmpzMatRef(coordinates))) {
                return false;
            }
            for (slong r = 0; r < n; ++r) {
                const fmpz* entry = fmpz_mat_entry(coordinates.raw(), 0, r);
                ::fmpz_set(fmpz_mat_entry(table.raw(), i * n + j, r), entry);
                ::fmpz_set(fmpz_mat_entry(table.raw(), j * n + i, r), entry);
            }
        }
    }
    if (n == 2) {
        bool fits = true;
        for (slong k = 0; k < 8 && fits; ++k) {
            const fmpz* entry = fmpz_mat_entry(table.raw(), k / 2, k % 2);
            fits = ::fmpz_bits(entry) <=
                   static_cast<flint_bitcnt_t>(FLINT_BITS / 2 - 1);
            prefilter.word_table[k] = fits ? ::fmpz_get_si(entry) : 0;
        }
        prefilter.word_table_ready = fits;
    }
    prefilter.multiplication_table = std::move(table);
    prefilter.multiplication_matrix = flint::FmpzMat(n, n);
    prefilter.degree = n;
    prefilter.ready = true;
    return true;
}

bool required_prime_norm_prefilter_norm(
        flint::FmpzRef out,
        RequiredPrimeNormPrefilter& prefilter,
        flint::FmpzMatConstRef coordinates) noexcept {
    const slong n = prefilter.degree;
    if (!prefilter.ready || ::fmpz_mat_nrows(coordinates.raw()) != 1 ||
        ::fmpz_mat_ncols(coordinates.raw()) != n) {
        return false;
    }
    // Source: the norm of a is the determinant of the Q-linear map x -> a x
    // (Cohen, A Course in Computational Algebraic Number Theory, GTM 138,
    // Section 4.3), here in the order basis, where its matrix is integral.
    // Column j holds the coordinates of a w_j = sum_i a_i w_i w_j, as in
    // reference 2.17.3 src/basemath/base3.c:zk_multable. The determinant is
    // exact integer arithmetic, so the value equals Element::norm.
    fmpz_mat_struct* matrix = prefilter.multiplication_matrix.raw();
    const fmpz_mat_struct* table = prefilter.multiplication_table.raw();
    const fmpz_mat_struct* a = coordinates.raw();
    for (slong r = 0; r < n; ++r) {
        for (slong j = 0; j < n; ++j) {
            fmpz* entry = fmpz_mat_entry(matrix, r, j);
            ::fmpz_zero(entry);
            for (slong i = 0; i < n; ++i) {
                ::fmpz_addmul(entry, fmpz_mat_entry(a, 0, i),
                              fmpz_mat_entry(table, i * n + j, r));
            }
        }
    }
    ::fmpz_mat_det(out.raw(), matrix);
    return true;
}

bool required_prime_norm_prefilter_rejects(
        bool& rejected,
        RequiredPrimeNormPrefilter& prefilter,
        flint::FmpzMatConstRef coordinates) noexcept {
    return required_prime_norm_prefilter_rejects_word(rejected, prefilter,
                                                      coordinates) ||
           required_prime_norm_prefilter_rejects_fmpz(rejected, prefilter,
                                                      coordinates);
}

bool required_prime_norm_prefilter_rejects_word(
        bool& rejected,
        const RequiredPrimeNormPrefilter& prefilter,
        flint::FmpzMatConstRef coordinates) noexcept {
    rejected = false;
    const fmpz_mat_struct* a = coordinates.raw();
    if (!prefilter.ready || prefilter.degree != 2 ||
        !prefilter.word_table_ready || ::fmpz_mat_nrows(a) != 1 ||
        ::fmpz_mat_ncols(a) != 2) {
        return false;
    }
    // At most required_prime_norm_prefilter_word_limit = 2^(FLINT_BITS / 2
    // - 1) - 1 in absolute value exactly when at most FLINT_BITS / 2 - 1 bits.
    constexpr flint_bitcnt_t limit_bits = FLINT_BITS / 2 - 1;
    const fmpz* a0_entry = fmpz_mat_entry(a, 0, 0);
    const fmpz* a1_entry = fmpz_mat_entry(a, 0, 1);
    if (::fmpz_bits(a0_entry) > limit_bits ||
        ::fmpz_bits(a1_entry) > limit_bits) {
        return false;
    }
    // The fmpz path below in words: the same multiplication matrix and
    // 2 x 2 determinant, exact because every coordinate and table entry is at
    // most 2^(FLINT_BITS / 2 - 1) - 1 in absolute value, so each entry
    // a_0 t_0 + a_1 t_1 is below 2^(FLINT_BITS - 1) and each product of two
    // entries below 2^(2 FLINT_BITS - 2) in absolute value.
    const slong a0 = ::fmpz_get_si(a0_entry);
    const slong a1 = ::fmpz_get_si(a1_entry);
    const slong* t = prefilter.word_table;
    // Entry (r, j) is sum_i a_i t[(i * 2 + j) * 2 + r].
    const slong m00 = a0 * t[0] + a1 * t[4];
    const slong m01 = a0 * t[2] + a1 * t[6];
    const slong m10 = a0 * t[1] + a1 * t[5];
    const slong m11 = a0 * t[3] + a1 * t[7];
    ulong left_hi = 0, left_lo = 0, right_hi = 0, right_lo = 0;
    smul_ppmm(left_hi, left_lo, m00, m11);
    smul_ppmm(right_hi, right_lo, m01, m10);
    ulong det_hi = 0, det_lo = 0;
    sub_ddmmss(det_hi, det_lo, left_hi, left_lo, right_hi, right_lo);
    if (static_cast<slong>(det_hi) < 0) {
        // |det| < 2^(2 FLINT_BITS - 1), so negation does not overflow.
        det_lo = -det_lo;
        det_hi = ~det_hi + (det_lo == 0 ? 1 : 0);
    }
    if (det_hi != 0) {
        return false;
    }
    if (det_lo == 0) {
        return true;
    }
    // The Z_ppo loop of the fmpz path below on words; gcd(N, P) is
    // gcd(N, P mod N).
    ulong cofactor = det_lo;
    ulong divisor =
            n_gcd(cofactor, ::fmpz_fdiv_ui(prefilter.smooth_product.raw(),
                                           cofactor));
    while (divisor != 1) {
        cofactor /= divisor;
        divisor = n_gcd(cofactor, divisor);
    }
    rejected = cofactor != 1;
    return true;
}

bool required_prime_norm_prefilter_rejects_fmpz(
        bool& rejected,
        RequiredPrimeNormPrefilter& prefilter,
        flint::FmpzMatConstRef coordinates) noexcept {
    rejected = false;
    if (!required_prime_norm_prefilter_norm(
                flint::FmpzRef(prefilter.norm), prefilter, coordinates)) {
        return false;
    }
    if (::fmpz_is_zero(prefilter.norm.raw()) != 0) {
        return true;
    }
    // Source: buch2.c:can_factor rejects when Z_ppo(N, F->prodZ) is not 1;
    // reference 2.17.3 src/basemath/base4.c:Z_ppo divides out gcd(x, f)
    // until it is 1. The cofactor is 1 exactly when the direct path's
    // cofactor, |N| with every factor-base rational prime and the required
    // rational prime removed, is 1.
    ::fmpz_abs(prefilter.cofactor.raw(), prefilter.norm.raw());
    ::fmpz_gcd(prefilter.divisor.raw(), prefilter.cofactor.raw(),
               prefilter.smooth_product.raw());
    while (::fmpz_is_one(prefilter.divisor.raw()) == 0) {
        ::fmpz_divexact(prefilter.cofactor.raw(), prefilter.cofactor.raw(),
                        prefilter.divisor.raw());
        ::fmpz_gcd(prefilter.divisor.raw(), prefilter.cofactor.raw(),
                   prefilter.divisor.raw());
    }
    rejected = ::fmpz_is_one(prefilter.cofactor.raw()) == 0;
    return true;
}

}  // namespace detail

bool ideal_factor_over_base(flint::FmpzMatRef exponents,
                            const Ideal& ideal,
                            const FactorBase& base,
                            const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.factor_over_base");
    if (ideal.parent() == nullptr ||
        !same_order_parent(ideal.parent(), base.parent()) ||
        flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length()) {
        return false;
    }

    bool smooth = false;
    flint::FmpzMat candidate(1, base.length());
    if (!ideal_factor_over_base_direct(flint::FmpzMatRef(candidate), smooth,
                                       ideal, base, diagnostics) ||
        !smooth) {
        return false;
    }

    flint::fmpz_mat_set(exponents, flint::FmpzMatConstRef(candidate));
    return true;
}

bool ideal_factor_over_base(flint::FmpzMatRef exponents,
                            const FractionalIdeal& ideal,
                            const FactorBase& base,
                            const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::ideal,
                        "ideal.frac_factor_over_base");
    if (ideal.parent() == nullptr ||
        !same_order_parent(ideal.parent(), base.parent()) ||
        flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length()) {
        return false;
    }

    Ideal numerator(*ideal.parent());
    flint::Fmpz denominator;
    if (!ideal.get_integral_den(numerator, flint::FmpzRef(denominator))) {
        return false;
    }
    if (fmpz_is_one(denominator.raw()) != 0) {
        return ideal_factor_over_base(exponents, numerator, base, diagnostics);
    }

    bool smooth = false;
    flint::FmpzMat candidate(1, base.length());
    if (!fractional_ideal_factor_over_base_direct(
                flint::FmpzMatRef(candidate), smooth, numerator,
                flint::FmpzConstRef(denominator), base, diagnostics) ||
        !smooth) {
        return false;
    }

    flint::fmpz_mat_set(exponents, flint::FmpzMatConstRef(candidate));
    return true;
}

}  // namespace silex
