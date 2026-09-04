#include <silex/ideal_factorization.hpp>

#include <silex/diagnostics.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpq_poly.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mat.h>

#include <limits>
#include <utility>
#include <vector>

#include "ideal_factorization_internal.hpp"

namespace silex {
namespace {

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

    bool reference_matches = true;
    slong required_exponent = 0;
    for (slong i = 0; i < factorization.length(); ++i) {
        slong exponent = 0;
        if (!factorization.prime(factor, i) ||
            !factorization.exponent(exponent, i)) {
            return false;
        }
        if (factor.equal(required_prime)) {
            required_exponent += exponent;
        } else if (!base.contains(factor)) {
            reference_matches = false;
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
    if (!base.rational_prime_block_index_for_prime(block_index,
                                                   rational_prime)) {
        return true;
    }

    flint::Fmpz block_prime;
    slong block_length = 0;
    if (!base.rational_prime_block_data(flint::FmpzRef(block_prime),
                                        block_length, block_index) ||
        !flint::fmpz_equal(flint::FmpzConstRef(block_prime), rational_prime)) {
        return false;
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
        if (!detail::prime_ideal_valuation_with_norm_vp(
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

bool factor_base_accounts_for_order_element_rational_valuation(
        bool& accounted,
        const FactorBase& base,
        const OrderElement& element,
        flint::FmpzConstRef rational_prime,
        slong norm_valuation,
        slong remaining_valuation,
        const DiagnosticsContext* diagnostics) noexcept {
    slong remaining = remaining_valuation;
    if (!factor_base_order_element_rational_valuation(
                remaining, nullptr, base, element, rational_prime,
                norm_valuation, diagnostics)) {
        return false;
    }
    accounted = remaining == 0;
    return true;
}

// Source factor-generation, factor-admission, and prime-division logic retains
// the element while factoring its norm and evaluating prime-ideal
// valuations.  Keep that representation boundary while using Silex/FLINT
// primitives for the exact operations.
bool order_element_factor_over_base_with_required_prime_direct(
        bool& matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics) noexcept {
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
    if (!ambient.is_defined() || !element.get_element(ambient) ||
        !ambient.norm(flint::FmpqRef(norm)) ||
        fmpz_is_one(fmpq_denref(norm.raw())) == 0) {
        return false;
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
    flint::fmpz_factor(flint::FmpzFactorRef(factorization),
                       flint::FmpzConstRef(norm_abs));

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
            if (!detail::prime_ideal_valuation_with_norm_vp(
                        required_valuation, required_prime, element,
                        norm_valuation, diagnostics)) {
                return false;
            }
            if (required_valuation != 1) {
                return true;
            }
            remaining_valuation -= required_prime.residue_degree();
            if (remaining_valuation < 0) {
                return true;
            }
        }

        bool accounted = false;
        if (!factor_base_accounts_for_order_element_rational_valuation(
                    accounted, base, element,
                    flint::FmpzConstRef(rational_prime), norm_valuation,
                    remaining_valuation, diagnostics)) {
            return false;
        }
        if (!accounted) {
            return true;
        }
    }

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

    if (required_prime != nullptr) {
        slong valuation = -1;
        if (!required_prime->valuation(valuation, ideal, diagnostics)) {
            return false;
        }
        if (valuation != 1) {
            smooth = false;
            return true;
        }

        flint::Fmpz p;
        flint::Fmpz required_norm;
        const slong residue_degree = required_prime->residue_degree();
        if (residue_degree <= 0 ||
            !required_prime->rational_prime(flint::FmpzRef(p))) {
            return false;
        }
        flint::fmpz_pow_ui(flint::FmpzRef(required_norm),
                           flint::FmpzConstRef(p),
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

bool set_element_fmpz(Element& out, flint::FmpzConstRef value) noexcept {
    if (!out.is_defined()) {
        return false;
    }

    flint::FmpqPoly polynomial;
    flint::Fmpq coefficient;
    fmpq_set_fmpz(coefficient.raw(), value.raw());
    flint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool set_principal_fmpz(Ideal& out, flint::FmpzConstRef value) noexcept {
    const Order* parent = out.parent();
    if (parent == nullptr || parent->parent() == nullptr) {
        return false;
    }

    Element element(*parent->parent());
    OrderElement order_element(*parent);
    return set_element_fmpz(element, value) &&
           order_element.set_element(element) &&
           out.set_principal(order_element);
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

        PrimeIdealList decomposed;
        if (!decompose_prime(decomposed, *order,
                             flint::FmpzConstRef(rational_prime), 0,
                             diagnostics)) {
            return false;
        }
        for (slong j = 0; j < decomposed.size(); ++j) {
            const PrimeIdeal* prime = decomposed.at(j);
            if (prime == nullptr) {
                return false;
            }
            if (base.contains(*prime)) {
                continue;
            }

            slong valuation = -1;
            if (!prime_ideal_valuation_with_norm_vp(
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

            if (!have_large_prime) {
                if (!candidate_large_prime.set(*prime)) {
                    return false;
                }
                have_large_prime = true;
                large_exponent = valuation;
            } else if (candidate_large_prime.equal(*prime)) {
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

bool order_element_factor_over_base_with_required_prime(
        bool& matches,
        const OrderElement& element,
        const FactorBase& base,
        const PrimeIdeal& required_prime,
        const DiagnosticsContext* diagnostics) noexcept {
    bool direct_matches = false;
    const bool direct_ok =
            order_element_factor_over_base_with_required_prime_direct(
                    direct_matches, element, base, required_prime, diagnostics);

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
    Ideal denominator_ideal(*ideal.parent());
    flint::Fmpz denominator;
    flint::FmpzMat numerator_row(1, base.length());
    flint::FmpzMat denominator_row(1, base.length());
    flint::FmpzMat candidate(1, base.length());

    if (!ideal.get_integral_den(numerator, flint::FmpzRef(denominator)) ||
        !ideal_factor_over_base(flint::FmpzMatRef(numerator_row),
                                numerator, base, diagnostics)) {
        return false;
    }

    if (fmpz_is_one(denominator.raw()) != 0) {
        flint::fmpz_mat_zero(flint::FmpzMatRef(denominator_row));
    } else if (!set_principal_fmpz(denominator_ideal,
                                  flint::FmpzConstRef(denominator)) ||
               !ideal_factor_over_base(flint::FmpzMatRef(denominator_row),
                                       denominator_ideal, base, diagnostics)) {
        return false;
    }

    for (slong i = 0; i < base.length(); ++i) {
        fmpz_sub(flint::fmpz_mat_entry(candidate, 0, i).raw(),
                 flint::fmpz_mat_entry(numerator_row, 0, i).raw(),
                 flint::fmpz_mat_entry(denominator_row, 0, i).raw());
    }

    flint::fmpz_mat_set(exponents, flint::FmpzMatConstRef(candidate));
    return true;
}

}  // namespace silex
