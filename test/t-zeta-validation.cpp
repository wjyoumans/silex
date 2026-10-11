#include "order_unit/class_unit_transaction_internal.hpp"
#include "order_unit/compute_internal.hpp"
#include "zeta/zeta_internal.hpp"

#include "test_support.hpp"

#include <silex/class_group.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/zeta.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <optional>
#include <vector>

#include <flint/ulong_extras.h>

namespace {
namespace sflint = silex::flint;

struct FieldSetup {
    silex::NumberField field;
    silex::Order maximal_order;
};

FieldSetup setup_from_coefficients(const slong* coefficients,
                                   slong degree) noexcept {
    assert(coefficients != nullptr);
    assert(degree > 0);

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, degree, 1);
    for (slong i = 0; i < degree; ++i) {
        if (coefficients[i] != 0) {
            sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
        }
    }

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(equation_order));
    assert(setup.maximal_order.is_maximal());
    return setup;
}

FieldSetup cubic_blocker() noexcept {
    const slong coefficients[] = {-5, -2, 0};
    return setup_from_coefficients(coefficients, 3);
}

FieldSetup quadratic_specialized_control() noexcept {
    const slong coefficients[] = {-5, 0};
    return setup_from_coefficients(coefficients, 2);
}

FieldSetup quartic_blocker() noexcept {
    const slong coefficients[] = {4, 1, 3, -8};
    return setup_from_coefficients(coefficients, 4);
}

void zeta_validation_target(sflint::Arb& out, slong precision) noexcept {
    sflint::Fmpq target;
    sflint::fmpq_set_si(target, 6931, 20000);
    sflint::arb_set_fmpq(out, target, precision);
}

void zeta_validation_tail_target(sflint::Arb& out, slong precision) noexcept {
    sflint::Arb validation_target;
    sflint::Arb numerical_reserve;
    zeta_validation_target(validation_target, precision);
    sflint::arb_one(numerical_reserve);
    sflint::arb_mul_2exp_si(numerical_reserve, numerical_reserve, -20);
    sflint::arb_sub(out, validation_target, numerical_reserve, precision);
}

bool radius_lt(const sflint::Arb& value,
               const sflint::Arb& target) noexcept {
    sflint::Arb radius;
    sflint::arb_get_rad_arb(radius, value);
    return sflint::arb_lt(radius, target);
}

bool call_validation(sflint::Arb& product,
                           sflint::Arb& error_bound,
                           ulong& cutoff,
                           slong& work_precision,
                           const silex::Order& order,
                           ulong max_cutoff,
                           slong precision) noexcept {
    return silex::detail::
            zeta_class_regulator_product_validation_with_diagnostics(
                    sflint::ArbRef(product), sflint::ArbRef(error_bound),
                    cutoff, work_precision, order, max_cutoff, precision,
                    nullptr);
}

bool check_failed_call_preserves_outputs(const silex::Order& order,
                                         ulong max_cutoff) noexcept {
    sflint::Arb product;
    sflint::Arb error_bound;
    sflint::arb_set_si(product, 17);
    sflint::arb_set_si(error_bound, 19);
    ulong cutoff = 23;
    slong work_precision = 29;

    return !call_validation(product, error_bound, cutoff,
                                  work_precision, order, max_cutoff, 128) &&
           ::arb_equal_si(product.raw(), 17) != 0 &&
           ::arb_equal_si(error_bound.raw(), 19) != 0 && cutoff == 23 &&
           work_precision == 29;
}

bool check_producer(const FieldSetup& setup,
                    ulong expected_cutoff,
                    ulong strict_max_cutoff) noexcept {
    const ulong preceding_cutoff = expected_cutoff - 9;
    if (!check_failed_call_preserves_outputs(setup.maximal_order,
                                             preceding_cutoff)) {
        return false;
    }

    sflint::Arb product;
    sflint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;
    if (!call_validation(product, error_bound, cutoff, work_precision,
                               setup.maximal_order, expected_cutoff, 128) ||
        cutoff != expected_cutoff || work_precision != 64 ||
        !sflint::arb_is_finite(product) ||
        !sflint::arb_is_positive(product) ||
        !sflint::arb_is_finite(error_bound) ||
        !sflint::arb_is_positive(error_bound)) {
        return false;
    }

    sflint::Arb tail_target;
    zeta_validation_tail_target(tail_target, 128);
    if (!sflint::arb_lt(error_bound, tail_target)) {
        return false;
    }

    sflint::Arb log_product;
    sflint::Arb validation_target;
    sflint::arb_log(log_product, product, 128);
    zeta_validation_target(validation_target, 128);
    if (!sflint::arb_is_finite(log_product) ||
        !radius_lt(log_product, validation_target)) {
        return false;
    }

    auto strict = silex::zeta_class_regulator_product_bf_audit(
            setup.maximal_order, strict_max_cutoff, 128);
    return strict.has_value() &&
           sflint::arb_overlaps(product, strict->value);
}

bool check_cache_separation(const FieldSetup& setup) noexcept {
    constexpr ulong max_cutoff = 20000;
    constexpr slong precision = 128;
    silex::detail::AnalyticClassRegulatorCache cache;
    cache.configure_validation(max_cutoff);
    if (!cache.ensure(setup.maximal_order, precision + 32, nullptr, nullptr,
                      true) ||
        !cache.has_validation(max_cutoff, precision + 32) ||
        cache.has_validation(max_cutoff, precision) ||
        !cache.ensure(setup.maximal_order, precision, nullptr, nullptr, true) ||
        !cache.has_validation(max_cutoff, precision) ||
        cache.has_bf_audit(max_cutoff, precision)) {
        return false;
    }

    sflint::Arb validation_snapshot;
    sflint::arb_set(sflint::ArbRef(validation_snapshot), cache.value());
    if (!cache.ensure_bf_audit(setup.maximal_order, max_cutoff, precision,
                               nullptr) ||
        !cache.has_bf_audit(max_cutoff, precision) ||
        ::arb_equal(validation_snapshot.raw(), cache.value().raw()) == 0) {
        return false;
    }

    auto strict = silex::zeta_class_regulator_product_bf_audit(
            setup.maximal_order, max_cutoff, precision);
    return strict.has_value() && cache.bf_cutoff() == strict->cutoff &&
           cache.bf_work_precision() == strict->work_precision &&
           ::arb_equal(cache.bf_error_bound().raw(),
                       strict->error_bound.raw()) != 0 &&
           sflint::arb_overlaps(cache.bf_value(), strict->value) &&
           sflint::arb_overlaps(cache.value(), validation_snapshot);
}

bool configure_native_options(silex::ClassGroupComputeOptions& options,
                              sflint::Fmpz& factor_base_bound,
                              const silex::Order& order) noexcept {
    options = silex::ClassGroupComputeOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;
    return true;
}

bool check_transaction_failure_leaves_outputs_unset(
        const FieldSetup& setup) noexcept {
    silex::ClassGroupComputeOptions options;
    sflint::Fmpz factor_base_bound;
    if (!configure_native_options(options, factor_base_bound,
                                  setup.maximal_order)) {
        return false;
    }
    options.max_candidates = 0;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport audit;
    const bool computed = silex::detail::compute_class_unit_transaction(
            units, class_group, setup.maximal_order,
            sflint::FmpzConstRef(factor_base_bound), options, 128, audit);
    // A computation that fails without reaching the honesty checkpoint keeps
    // the generic reason.
    return !computed &&
           audit.failure_stage == silex::detail::ClassUnitStage::total &&
           audit.failure_reason != nullptr &&
           std::strcmp(audit.failure_reason,
                       "class_unit_computation_failed") == 0 &&
           !audit.factor_base_honesty_failure.recorded &&
           !audit.final_result_published &&
           audit.class_group_certification ==
                   silex::CertificationMode::unknown &&
           audit.unit_group_certification ==
                   silex::CertificationMode::unknown &&
           !class_group.has_factor_base() && !class_group.has_presentation() &&
           !units.is_set();
}



// Residue degrees of the primes of `order` above `p`, sorted, from the
// general prime-decomposition routine.
bool decomposition_residue_degrees(std::vector<slong>& out,
                                   const silex::Order& order,
                                   ulong p) noexcept {
    sflint::Fmpz prime;
    sflint::fmpz_set_ui(sflint::FmpzRef(prime), p);
    silex::PrimeIdealList primes;
    if (!silex::decompose_prime(primes, order,
                                sflint::FmpzConstRef(prime))) {
        return false;
    }
    out.clear();
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* ideal = primes.at(i);
        if (ideal == nullptr || ideal->residue_degree() <= 0) {
            return false;
        }
        out.push_back(ideal->residue_degree());
    }
    std::sort(out.begin(), out.end());
    return !out.empty();
}

void cache_store(silex::detail::ZetaBfResidueDegreeCache& cache,
                 ulong p,
                 const std::vector<slong>& degrees) noexcept {
    const std::size_t offset = cache.residue_degrees.size();
    cache.residue_degrees.insert(cache.residue_degrees.end(),
                                 degrees.begin(), degrees.end());
    cache.entries.push_back(silex::detail::ZetaBfResidueDegreeCacheEntry{
            p, offset, degrees.size()});
}

bool call_bf_audit(sflint::Arb& product,
                   sflint::Arb& error_bound,
                   ulong& cutoff,
                   slong& work_precision,
                   const silex::Order& order,
                   silex::detail::ZetaBfResidueDegreeCache* cache) noexcept {
    return silex::detail::
            zeta_class_regulator_product_bf_audit_with_diagnostics(
                    sflint::ArbRef(product), sflint::ArbRef(error_bound),
                    cutoff, work_precision, order, 20000, 128, nullptr,
                    nullptr, cache);
}

struct ReferenceDecomposition {
    ulong p;
    slong degrees[4];  // sorted residue degrees, zero-terminated
};

// Residue-degree differential test.  The BF term normally gets residue
// degrees from the defining polynomial mod p (with fallbacks for primes
// dividing the index).  Prefilling the residue-degree cache with degrees
// from decompose_prime for every prime below the cutoff replaces that path,
// so the two BF results must agree exactly.  decompose_prime itself is
// anchored to the external reference GP 2.17.4,
// `idealprimedec(nfinit(P), p)[j].f`.
bool check_residue_degree_differential(
        const FieldSetup& setup,
        const ReferenceDecomposition* anchors,
        std::size_t anchor_count,
        ulong perturbed_prime) noexcept {
    const silex::Order& order = setup.maximal_order;
    std::vector<slong> degrees;
    for (std::size_t i = 0; i < anchor_count; ++i) {
        std::vector<slong> expected;
        for (slong degree : anchors[i].degrees) {
            if (degree == 0) {
                break;
            }
            expected.push_back(degree);
        }
        if (!decomposition_residue_degrees(degrees, order, anchors[i].p) ||
            degrees != expected) {
            return false;
        }
    }

    sflint::Arb reference;
    sflint::Arb reference_error;
    ulong reference_cutoff = 0;
    slong reference_precision = 0;
    if (!call_bf_audit(reference, reference_error, reference_cutoff,
                       reference_precision, order, nullptr)) {
        return false;
    }

    silex::detail::ZetaBfResidueDegreeCache cache;
    silex::detail::ZetaBfResidueDegreeCache perturbed;
    bool perturbed_stored = false;
    for (ulong p = 2; p < reference_cutoff; p = n_nextprime(p, 1)) {
        if (!decomposition_residue_degrees(degrees, order, p)) {
            return false;
        }
        cache_store(cache, p, degrees);
        if (p == perturbed_prime) {
            // Move all of the residue degree into one prime ideal.
            slong total = 0;
            for (slong degree : degrees) {
                total += degree;
            }
            if (degrees.size() == 1) {
                return false;
            }
            cache_store(perturbed, p, std::vector<slong>{total});
            perturbed_stored = true;
        } else {
            cache_store(perturbed, p, degrees);
        }
    }
    if (!perturbed_stored) {
        return false;
    }

    sflint::Arb candidate;
    sflint::Arb candidate_error;
    ulong candidate_cutoff = 0;
    slong candidate_precision = 0;
    if (!call_bf_audit(candidate, candidate_error, candidate_cutoff,
                       candidate_precision, order, &cache) ||
        candidate_cutoff != reference_cutoff ||
        candidate_precision != reference_precision ||
        ::arb_equal(candidate.raw(), reference.raw()) == 0 ||
        ::arb_equal(candidate_error.raw(), reference_error.raw()) == 0) {
        return false;
    }

    // Sensitivity check: a cache with one wrong decomposition is consumed
    // and changes the BF value, so the comparison above is not vacuous.
    if (!call_bf_audit(candidate, candidate_error, candidate_cutoff,
                       candidate_precision, order, &perturbed)) {
        return false;
    }
    return ::arb_equal(candidate.raw(), reference.raw()) == 0;
}

bool check_residue_degree_differentials() noexcept {
    // x^3 - x - 1, disc -23.
    const slong cubic_coefficients[] = {-1, -1, 0};
    const ReferenceDecomposition cubic_anchors[] = {
            {2, {3}},       {3, {3}},    {5, {1, 2}},  {7, {1, 2}},
            {11, {1, 2}},   {13, {3}},   {23, {1, 1}}, {503, {1, 2}},
    };
    // x^3 - x^2 - 2x - 8 (Dedekind): 2 divides the index of Z[x] and
    // splits completely, so the defining polynomial mod 2 is unusable.
    const slong dedekind_coefficients[] = {-8, -2, -1};
    const ReferenceDecomposition dedekind_anchors[] = {
            {2, {1, 1, 1}}, {3, {3}},    {5, {1, 2}},  {7, {3}},
            {11, {3}},      {13, {3}},   {23, {3}},    {503, {1, 1}},
    };
    // x^4 - 4x^2 + 2: 2 is totally ramified.
    const slong quartic_coefficients[] = {2, 0, -4, 0};
    const ReferenceDecomposition quartic_anchors[] = {
            {2, {1}},       {3, {4}},    {5, {4}},     {7, {2, 2}},
            {11, {4}},      {13, {4}},   {23, {2, 2}}, {503, {2, 2}},
    };

    return check_residue_degree_differential(
                   setup_from_coefficients(cubic_coefficients, 3),
                   cubic_anchors, 8, 5) &&
           check_residue_degree_differential(
                   setup_from_coefficients(dedekind_coefficients, 3),
                   dedekind_anchors, 8, 2) &&
           check_residue_degree_differential(
                   setup_from_coefficients(quartic_coefficients, 4),
                   quartic_anchors, 8, 7);
}


struct BdfStopCase {
    slong degree;
    slong coefficients[4];
    slong discriminant;
    // Bach's bound ceil(c log^2 |D|), c = 6 in degree two and 12 otherwise
    // (Bach 1990, Theorem 4 and the remark following it), computed
    // independently of Silex.
    slong bach;
    // The full BDF search: its bound, the smallest integer cutoff where
    // the criterion holds (recomputed independently of Silex at 60 digits
    // from the same sum over prime powers below the cutoff), and
    // its number of evaluations: doubling from 2 reaches the first power of
    // two 2^j at least the bound in j evaluations, and the bisection of
    // (2^(j - 1), 2^j] takes j - 1 more.
    slong bdf;
    slong evaluations;
    // Largest doubling cutoff where the criterion fails, 2^(j - 1).
    slong last_failing_doubling;
};

// The BDF search behind the GRH generation bound stops once the criterion
// fails at a cutoff at least Bach's bound.  The criterion's right side does
// not decrease as the cutoff grows, so the full search would return a bound
// above that cutoff, and the published bound min(BDF, Bach) is unchanged.
bool check_bdf_search_stop() noexcept {
    // The criterion fails at the first cutoff, 2, for every field, so the
    // doubling fails at least once before it holds.
    const BdfStopCase cases[] = {
            {3, {200, 1, 0, 0}, -1080004, 2317, 19, 9, 16},
            {2, {4, 1, 0, 0}, -15, 45, 5, 5, 4},
            {3, {1000003, 1, 0, 0}, -27000162000247, 11478, 472, 17, 256},
            {2, {1000003, 1, 0, 0}, -4000011, 1387, 123, 13, 64},
            {4, {10007, 1, 3, 0}, 256422632363321, 13210, 550, 19, 512},
    };
    for (const BdfStopCase& c : cases) {
        const FieldSetup setup = setup_from_coefficients(c.coefficients,
                                                         c.degree);
        sflint::Fmpz discriminant;
        if (!setup.maximal_order.discriminant(
                    sflint::FmpzRef(discriminant)) ||
            !sflint::fmpz_equal_si(discriminant, c.discriminant)) {
            return false;
        }

        // The full search, without a stop bound.
        sflint::Fmpz bdf;
        slong evaluations = 0;
        if (silex::detail::grh_bdf_factor_base_bound_search(
                    sflint::FmpzRef(bdf), setup.maximal_order, nullptr,
                    nullptr, &evaluations) !=
                    silex::detail::BdfSearchOutcome::found ||
            !sflint::fmpz_equal_si(bdf, c.bdf) ||
            evaluations != c.evaluations) {
            return false;
        }

        // The published bound is min(BDF, Bach), as before; Bach is
        // reported only when strictly smaller.
        sflint::Fmpz bound;
        bool bach_selected = c.bach >= c.bdf;
        if (!silex::detail::grh_factor_base_bound_with_diagnostics(
                    sflint::FmpzRef(bound), setup.maximal_order, nullptr,
                    &bach_selected) ||
            !sflint::fmpz_equal_si(bound, FLINT_MIN(c.bdf, c.bach)) ||
            bach_selected != (c.bach < c.bdf)) {
            return false;
        }

        // With Bach's bound as the stop, these fields never fail the
        // criterion at a cutoff that large, so the search runs in full.
        sflint::Fmpz stop;
        sflint::Fmpz stopped_bdf;
        sflint::fmpz_set_si(sflint::FmpzRef(stop), c.bach);
        slong stop_evaluations = 0;
        if (silex::detail::grh_bdf_factor_base_bound_search(
                    sflint::FmpzRef(stopped_bdf), setup.maximal_order, &stop,
                    nullptr, &stop_evaluations) !=
                    silex::detail::BdfSearchOutcome::found ||
            !sflint::fmpz_equal_si(stopped_bdf, c.bdf) ||
            stop_evaluations != c.evaluations) {
            return false;
        }

        // A stop at the BDF bound itself never stops: every cutoff where the
        // criterion fails lies below the bound the search returns.  (The
        // literal rule of BDF 2008, Section 3, which stops once the doubled
        // cutoff passes the stop, would stop at the first power of two past
        // the bound.)
        sflint::fmpz_set_si(sflint::FmpzRef(stop), c.bdf);
        sflint::fmpz_zero(sflint::FmpzRef(stopped_bdf));
        if (silex::detail::grh_bdf_factor_base_bound_search(
                    sflint::FmpzRef(stopped_bdf), setup.maximal_order, &stop,
                    nullptr, &stop_evaluations) !=
                    silex::detail::BdfSearchOutcome::found ||
            !sflint::fmpz_equal_si(stopped_bdf, c.bdf) ||
            stop_evaluations != c.evaluations) {
            return false;
        }

        // A stop at or below a cutoff where the doubling fails ends the
        // search there, after one evaluation per doubling step, and leaves
        // `out` unchanged; the bound the full search returns is above the
        // stop, so min(BDF, stop) is the stop.  This is how the GRH bound now
        // reaches Bach's bound without the full search.  With Bach's actual
        // bound as the stop it happens only for t(K) > 12, or t(K) > 6 in
        // degree 2 where Silex uses 6 log^2 D (BDF 2008, Section 1 after
        // Definition 1.2, Section 3, Theorem 4.3), that is for enormous
        // discriminants, so the stop here is set below Bach's bound.
        slong doublings = 0;
        for (slong t = c.last_failing_doubling; t > 1; t /= 2) {
            ++doublings;
        }
        const slong stops[] = {2, c.last_failing_doubling};
        for (const slong stop_value : stops) {
            sflint::fmpz_set_si(sflint::FmpzRef(stop), stop_value);
            sflint::fmpz_set_si(sflint::FmpzRef(stopped_bdf), 7);
            stop_evaluations = 0;
            if (silex::detail::grh_bdf_factor_base_bound_search(
                        sflint::FmpzRef(stopped_bdf), setup.maximal_order,
                        &stop, nullptr, &stop_evaluations) !=
                        silex::detail::BdfSearchOutcome::stopped ||
                !sflint::fmpz_equal_si(stopped_bdf, 7) ||
                stop_evaluations != (stop_value == 2 ? 1 : doublings) ||
                stop_evaluations >= c.evaluations ||
                stop_value >= c.bdf) {
                return false;
            }
        }

        // A stop between the last failing doubling and the bound is reached
        // in the bisection: it ends at the failing midpoint bdf - 1, after
        // the doublings and the bisection steps up to that midpoint, which is the
        // last evaluation of the full search.  It is
        // absent when no cutoff lies strictly between them.
        const slong bisection_stop = c.bdf - 1;
        if (bisection_stop > c.last_failing_doubling) {
            slong expected = doublings + 1;
            slong low = c.last_failing_doubling;
            slong high = 2 * c.last_failing_doubling;
            for (;;) {
                const slong middle = low + (high - low) / 2;
                ++expected;
                if (middle >= c.bdf) {
                    high = middle;
                } else if (middle >= bisection_stop) {
                    break;
                } else {
                    low = middle;
                }
            }
            sflint::fmpz_set_si(sflint::FmpzRef(stop), bisection_stop);
            sflint::fmpz_set_si(sflint::FmpzRef(stopped_bdf), 7);
            stop_evaluations = 0;
            if (silex::detail::grh_bdf_factor_base_bound_search(
                        sflint::FmpzRef(stopped_bdf), setup.maximal_order,
                        &stop, nullptr, &stop_evaluations) !=
                        silex::detail::BdfSearchOutcome::stopped ||
                !sflint::fmpz_equal_si(stopped_bdf, 7) ||
                stop_evaluations != expected ||
                stop_evaluations > c.evaluations) {
                return false;
            }
        }
    }
    return true;
}

// The grh index-one route sets its cutoff by Belabas-Friedman 2015,
// Theorem 7 (p. 12), with sigma = 1 + 1/sqrt(log D).  Reference values are
// a floating-point gp evaluation of the same formula over the prime ideals
// of norm below X (idealprimedec), to four decimals.
bool theorem7_bound_near(const silex::Order& order,
                         ulong cutoff,
                         double expected) noexcept {
    sflint::Arb bound;
    if (!silex::detail::zeta_bf_theorem7_error_bound(
                sflint::ArbRef(bound), order, cutoff, 128) ||
        !sflint::arb_is_finite(bound)) {
        return false;
    }
    sflint::Arb window;
    sflint::arb_set_d(window, expected);
    ::mag_set_d(arb_radref(window.raw()), 1e-4);
    return sflint::arb_contains(window, bound);
}

bool check_theorem7_bound_reference() noexcept {
    const slong simplest_cubic[] = {-1, -9, -6};
    const slong octic[] = {-1, -1, 0, 0, 0, 0, 0, 0};
    const FieldSetup cubic = setup_from_coefficients(simplest_cubic, 3);
    const FieldSetup x8 = setup_from_coefficients(octic, 8);
    sflint::Arb unused;
    return theorem7_bound_near(cubic.maximal_order, 279, 0.3293) &&
           theorem7_bound_near(cubic.maximal_order, 603, 0.1737) &&
           theorem7_bound_near(x8.maximal_order, 279, 0.4289) &&
           theorem7_bound_near(x8.maximal_order, 603, 0.2206) &&
           // Theorem 7 needs X > 9; the BF term needs a multiple of 9.
           !silex::detail::zeta_bf_theorem7_error_bound(
                   sflint::ArbRef(unused), cubic.maximal_order, 9, 128) &&
           !silex::detail::zeta_bf_theorem7_error_bound(
                   sflint::ArbRef(unused), cubic.maximal_order, 280, 128);
}

// The index-one value: the least multiple of 9 whose Theorem 7 bound is
// below 1/4 - 2^-20 (gp: 387 for the cubic, 522 for x^8 - x - 1, 531 for
// the multiquadratic octic; 81 for x^2 - 5), a log radius below 1/4, and an
// enclosure
// that meets the default route's hR.  A quadratic order with the
// unconditional L(1, chi) route has one too, as grh requests take it in
// every degree above one; its default hR is the L(1, chi) value, an
// independent check.  Degree one has none.
//
// Independent oracle for 81 (gp 2.17.4, realprecision 60): with
// nf = nfinit(x^2 - 5), n = 2, r = 2, D = 5, sigma = 1 + 1/sqrt(log D),
// s2 = 2 sigma - 1, delta0 the part of delta (Belabas-Friedman 2015,
// Theorem 7, p. 12; also at src/zeta/zeta.cpp:2926-2937) that does not
// depend on X,
//   delta0 = log D + (C/2 + 1 - log(4 pi)/2)/s2 + 2/(sigma - 1) + 2/sigma
//            - n (log 2 pi - psi(sigma))
//            - r (psi((sigma + 1)/2) - psi(sigma/2))/2,
// with C Euler's constant, psi the digamma function, and r the number of
// real places, S(X) = sum over idealprimedec prime ideals of norm N < X
// of log N / (N^sigma - 1), d = delta0 - 2 S(X), t = log(X/9),
// beta = (1/2 + 1/t)(X/9) log((X/9 + 1)/(X/9 - 1))/2, and
//   b(X) = 2.324 s2 / (sqrt(X) log(3X))
//          * (d (1 + 3.88/t) + 4.26 (n - 1) beta / (s2 sqrt(X))),
// the least X = 18, 27, ... with d > 0 and b(X) < 1/4 - 2^-20 is
// X = 81, b(81) = 0.22147147143...; the same script gives 387
// (b = 0.249044782...) and 522 (b = 0.249374529...) for the cubic and
// x^8 - x - 1, matching the values above.  Because the search is the least
// multiple of 9, X = 72 fails with b(72) = 0.25265, the minimality that the
// preceding-bound check below tests in Silex.  Silex scans from 72 (its
// minimum cutoff 70 rounded up to a multiple of 9); the oracle's wider range
// also rules out X = 18, ..., 63 (b = 2.13, 0.98, 0.61, 0.45, 0.36, 0.29).
bool check_index_one_product(const FieldSetup& setup,
                             ulong expected_cutoff) noexcept {
    sflint::Arb product;
    silex::detail::ZetaBfRouteAudit audit;
    if (!silex::detail::zeta_class_regulator_product_index_one_with_diagnostics(
                sflint::ArbRef(product), audit, setup.maximal_order, 128,
                nullptr) ||
        audit.cutoff != expected_cutoff || audit.max_cutoff != 20000 ||
        audit.requested_precision != 128 || audit.work_precision < 64 ||
        !sflint::arb_is_finite(audit.error_bound) ||
        !sflint::arb_is_positive(audit.error_bound)) {
        return false;
    }

    sflint::Arb quarter;
    sflint::Arb log_product;
    sflint::arb_set_d(quarter, 0.25);
    sflint::arb_log(log_product, product, 128);
    if (!sflint::arb_lt(audit.error_bound, quarter) ||
        !sflint::arb_is_finite(log_product) ||
        !radius_lt(log_product, quarter)) {
        return false;
    }
    sflint::Arb preceding_bound;
    sflint::Arb tail_target;
    sflint::arb_set_d(tail_target, 0.25 - 1.0 / 1048576.0);
    if (!silex::detail::zeta_bf_theorem7_error_bound(
                sflint::ArbRef(preceding_bound), setup.maximal_order,
                expected_cutoff - 9, 64) ||
        sflint::arb_lt(preceding_bound, tail_target)) {
        return false;
    }

    sflint::Arb default_product;
    return silex::zeta_class_regulator_product(
                   sflint::ArbRef(default_product), setup.maximal_order,
                   128) &&
           sflint::arb_overlaps(product, default_product);
}

bool check_index_one_products() noexcept {
    const slong simplest_cubic[] = {-1, -9, -6};
    const slong octic[] = {-1, -1, 0, 0, 0, 0, 0, 0};
    const slong multiquadratic[] = {1, 0, -3, 0, 8, 0, -3, 0};
    const slong linear[] = {0};
    const FieldSetup quadratic = quadratic_specialized_control();
    if (!silex::detail::zeta_unconditional_route_available(
                quadratic.maximal_order)) {
        return false;
    }
    const FieldSetup rationals = setup_from_coefficients(linear, 1);
    sflint::Arb product;
    sflint::arb_set_ui(product, 7);
    silex::detail::ZetaBfRouteAudit audit;
    audit.cutoff = 5;
    return check_index_one_product(
                   setup_from_coefficients(simplest_cubic, 3), 387) &&
           check_index_one_product(setup_from_coefficients(octic, 8), 522) &&
           check_index_one_product(
                   setup_from_coefficients(multiquadratic, 8), 531) &&
           check_index_one_product(quadratic, 81) &&
           !silex::detail::
                   zeta_class_regulator_product_index_one_with_diagnostics(
                           sflint::ArbRef(product), audit,
                           rationals.maximal_order, 128, nullptr) &&
           audit.cutoff == 5 && ::arb_equal_si(product.raw(), 7) != 0;
}

// The analytic cache serves the index-one value until the next ensure(),
// which returns the default value; the index-one value is not offered once
// the default value is cached, nor when the validation enclosure is
// enabled.
bool check_index_one_cache(const FieldSetup& setup,
                           ulong index_one_cutoff) noexcept {
    silex::detail::AnalyticClassRegulatorCache cache;
    if (!cache.ensure_index_one(setup.maximal_order, 128) ||
        cache.value_unconditional()) {
        return false;
    }
    const silex::detail::ZetaBfRouteAudit* audit = cache.value_bf_audit();
    if (audit == nullptr || audit->cutoff != index_one_cutoff) {
        return false;
    }
    sflint::Arb index_one_value;
    sflint::arb_set(sflint::ArbRef(index_one_value), cache.value());

    if (!cache.ensure(setup.maximal_order, 128)) {
        return false;
    }
    audit = cache.value_bf_audit();
    if (audit == nullptr || audit->cutoff == index_one_cutoff ||
        ::arb_equal(cache.value().raw(), index_one_value.raw()) != 0 ||
        cache.ensure_index_one(setup.maximal_order, 128) ||
        cache.value_bf_audit() != audit) {
        return false;
    }

    silex::detail::AnalyticClassRegulatorCache validation_cache;
    validation_cache.configure_validation(20000);
    return !validation_cache.ensure_index_one(setup.maximal_order, 128);
}

// A grh transaction accepted by the index-one test records the audit of
// the index-one value, not of the default route.  This holds for x^2 - d
// as well, although that presentation also has the unconditional L(1, chi)
// route.
bool check_grh_transaction_records_index_one_audit(
        const FieldSetup& setup) noexcept {
    silex::ClassGroupComputeOptions options;
    sflint::Fmpz factor_base_bound;
    if (!configure_native_options(options, factor_base_bound,
                                  setup.maximal_order)) {
        return false;
    }
    options.requested_certification = silex::CertificationMode::grh;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport report;
    if (!silex::detail::compute_class_unit_transaction(
                units, class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128,
                report)) {
        return false;
    }
    const std::optional<silex::ClassGroupZetaBfProofRecord> record =
            class_group.zeta_bf_proof_record();

    sflint::Arb product;
    silex::detail::ZetaBfRouteAudit audit;
    return record.has_value() &&
           silex::detail::
                   zeta_class_regulator_product_index_one_with_diagnostics(
                           sflint::ArbRef(product), audit,
                           setup.maximal_order, 128, nullptr) &&
           record->cutoff == audit.cutoff &&
           record->max_cutoff == audit.max_cutoff &&
           record->requested_precision == audit.requested_precision &&
           record->work_precision == audit.work_precision &&
           ::arb_equal(record->error_bound.raw(), audit.error_bound.raw()) !=
                   0;
}

// When the coarse value does not show index one, the transaction discards
// it and records the audit of the default route.  The degree-12 field
// x^12 - x^11 - x^9 + ... + 1 has a coarse cutoff of 2547 but candidate
// index bound above one at the first test.
bool check_grh_transaction_refines_to_default_audit(
        const FieldSetup& setup) noexcept {
    constexpr ulong coarse_cutoff = 2547;
    silex::detail::AnalyticClassRegulatorCache cache;
    if (!cache.ensure_index_one(setup.maximal_order, 128)) {
        return false;
    }
    const silex::detail::ZetaBfRouteAudit* audit = cache.value_bf_audit();
    if (audit == nullptr || audit->cutoff != coarse_cutoff ||
        !cache.ensure(setup.maximal_order, 128)) {
        return false;
    }
    audit = cache.value_bf_audit();
    if (audit == nullptr || audit->cutoff == coarse_cutoff) {
        return false;
    }
    const ulong default_cutoff = audit->cutoff;

    silex::ClassGroupComputeOptions options;
    sflint::Fmpz factor_base_bound;
    if (!configure_native_options(options, factor_base_bound,
                                  setup.maximal_order)) {
        return false;
    }
    options.requested_certification = silex::CertificationMode::grh;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport report;
    if (!silex::detail::compute_class_unit_transaction(
                units, class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128,
                report)) {
        return false;
    }
    const std::optional<silex::ClassGroupZetaBfProofRecord> record =
            class_group.zeta_bf_proof_record();
    return record.has_value() && record->cutoff == default_cutoff &&
           record->cutoff != coarse_cutoff;
}

// The grh default-route product (Belabas-Friedman residue at the default
// maximum cutoff), with its audit.
bool call_grh_default_product(
        sflint::Arb& product,
        silex::detail::ZetaBfRouteAudit& audit,
        const silex::Order& order,
        silex::detail::ZetaBfResidueDegreeCache* cache) noexcept {
    bool unconditional = true;
    return silex::detail::zeta_class_regulator_product_with_diagnostics(
                   sflint::ArbRef(product), order, 128, nullptr, nullptr,
                   cache, &unconditional, &audit, true) &&
           !unconditional;
}

bool same_grh_default_product(
        const sflint::Arb& product,
        const silex::detail::ZetaBfRouteAudit& audit,
        const sflint::Arb& reference,
        const silex::detail::ZetaBfRouteAudit& reference_audit) noexcept {
    return audit.cutoff == reference_audit.cutoff &&
           audit.work_precision == reference_audit.work_precision &&
           ::arb_equal(product.raw(), reference.raw()) != 0 &&
           ::arb_equal(audit.error_bound.raw(),
                       reference_audit.error_bound.raw()) != 0;
}

// The linear-factor counts stored in the residue-degree cache are the
// number of degree-one primes above p (checked against decompose_prime),
// a later Belabas-Friedman evaluation consumes them (a wrong count changes
// the value) and gets the uncached value from them, and a grh ensure() after
// ensure_index_one() on the shared cache gives the value and audit of a
// fresh cache.
bool check_linear_factor_count_cache(const FieldSetup& setup) noexcept {
    const silex::Order& order = setup.maximal_order;
    sflint::Arb reference;
    silex::detail::ZetaBfRouteAudit reference_audit;
    if (!call_grh_default_product(reference, reference_audit, order,
                                  nullptr)) {
        return false;
    }

    silex::detail::ZetaBfResidueDegreeCache cache;
    for (int pass = 0; pass < 2; ++pass) {
        sflint::Arb candidate;
        silex::detail::ZetaBfRouteAudit candidate_audit;
        if (!call_grh_default_product(candidate, candidate_audit, order,
                                      &cache) ||
            !same_grh_default_product(candidate, candidate_audit, reference,
                                      reference_audit) ||
            cache.linear_factor_counts.size() < 2) {
            return false;
        }
    }

    std::vector<slong> degrees;
    const std::size_t count_size = cache.linear_factor_counts.size();
    const std::size_t checked[] = {0, 1, count_size / 2, count_size - 1};
    for (std::size_t index : checked) {
        const silex::detail::ZetaBfLinearFactorCountEntry& entry =
                cache.linear_factor_counts[index];
        if (!decomposition_residue_degrees(degrees, order, entry.p) ||
            entry.count != static_cast<slong>(std::count(
                                   degrees.begin(), degrees.end(), 1))) {
            return false;
        }
    }

    silex::detail::ZetaBfResidueDegreeCache perturbed = cache;
    silex::detail::ZetaBfLinearFactorCountEntry& last =
            perturbed.linear_factor_counts.back();
    last.count = last.count == 0 ? 1 : last.count - 1;
    sflint::Arb perturbed_value;
    silex::detail::ZetaBfRouteAudit perturbed_audit;
    if (!call_grh_default_product(perturbed_value, perturbed_audit, order,
                                  &perturbed) ||
        ::arb_equal(perturbed_value.raw(), reference.raw()) != 0) {
        return false;
    }

    silex::detail::AnalyticClassRegulatorCache shared;
    if (!shared.ensure_index_one(order, 128) ||
        !shared.ensure(order, 128, nullptr, nullptr, false, true)) {
        return false;
    }
    const silex::detail::ZetaBfRouteAudit* shared_audit =
            shared.value_bf_audit();
    sflint::Arb shared_value;
    sflint::arb_set(sflint::ArbRef(shared_value), shared.value());
    if (shared_audit == nullptr ||
        !same_grh_default_product(shared_value, *shared_audit, reference,
                                  reference_audit)) {
        return false;
    }

    // Every count of the cache from the grh run above, changed.  Seeding it
    // into an AnalyticClassRegulatorCache is observable only through the
    // value that a later grh ensure() computes from it.
    silex::detail::ZetaBfResidueDegreeCache wrong = cache;
    for (silex::detail::ZetaBfLinearFactorCountEntry& entry :
         wrong.linear_factor_counts) {
        entry.count = entry.count == 0 ? 1 : entry.count - 1;
    }

    // A grh ensure() reads the shared cache (a cache hit): seeded wrong
    // counts change its value, so ensure() does not pass nullptr.
    silex::detail::AnalyticClassRegulatorCache seeded;
    if (wrong.entries.empty() ||
        !seeded.seed_bf_residue_degree_cache(order, wrong, nullptr) ||
        !seeded.ensure(order, 128, nullptr, nullptr, false, true) ||
        ::arb_equal(seeded.value().raw(), reference.raw()) != 0) {
        return false;
    }

    // The seed check rejects counts out of order, below 0 or above the
    // degree, and a rejected seed leaves the cache as it was.
    silex::detail::ZetaBfResidueDegreeCache unsorted = cache;
    std::swap(unsorted.linear_factor_counts[0],
              unsorted.linear_factor_counts[1]);
    silex::detail::ZetaBfResidueDegreeCache negative = cache;
    negative.linear_factor_counts[0].count = -1;
    silex::detail::ZetaBfResidueDegreeCache too_large = cache;
    too_large.linear_factor_counts[0].count = order.degree() + 1;
    silex::detail::AnalyticClassRegulatorCache rejecting;
    if (rejecting.seed_bf_residue_degree_cache(order, unsorted, nullptr) ||
        rejecting.seed_bf_residue_degree_cache(order, negative, nullptr) ||
        rejecting.seed_bf_residue_degree_cache(order, too_large, nullptr) ||
        !rejecting.ensure(order, 128, nullptr, nullptr, false, true) ||
        ::arb_equal(rejecting.value().raw(), reference.raw()) == 0) {
        return false;
    }

    // Switching to a different order clears the counts: the wrong counts of
    // `order` must not reach the value of another order.
    const FieldSetup other = cubic_blocker();
    sflint::Arb other_reference;
    silex::detail::ZetaBfRouteAudit other_audit;
    silex::detail::AnalyticClassRegulatorCache switched;
    return call_grh_default_product(other_reference, other_audit,
                                    other.maximal_order, nullptr) &&
           switched.seed_bf_residue_degree_cache(order, wrong, nullptr) &&
           switched.ensure(other.maximal_order, 128, nullptr, nullptr, false,
                           true) &&
           ::arb_equal(switched.value().raw(), other_reference.raw()) != 0;
}
}  // namespace

int main() {
    const FieldSetup quadratic = quadratic_specialized_control();
    const FieldSetup cubic = cubic_blocker();
    const FieldSetup quartic = quartic_blocker();
    const slong simplest_cubic_a6[] = {-1, -9, -6};
    const slong d12_coefficients[] = {1, 0, -2, 0, 2, 1, 1, 1, 2, -1, 0, -1};
    const slong x2_minus_210[] = {-210, 0};

    return !check_failed_call_preserves_outputs(quadratic.maximal_order,
                                                20000) ||
                   !check_producer(cubic, 1008, 20000) ||
                   !check_producer(quartic, 1971, 30000) ||
                   !check_cache_separation(cubic) ||
                   !check_transaction_failure_leaves_outputs_unset(cubic) ||
                   !check_transaction_failure_leaves_outputs_unset(quartic) ||
                   !check_residue_degree_differentials() ||
                   !check_bdf_search_stop() ||
                   !check_theorem7_bound_reference() ||
                   !check_index_one_products() ||
                   !check_index_one_cache(
                           setup_from_coefficients(simplest_cubic_a6, 3),
                           387) ||
                   !check_grh_transaction_records_index_one_audit(
                           setup_from_coefficients(simplest_cubic_a6, 3)) ||
                   !check_grh_transaction_records_index_one_audit(
                           setup_from_coefficients(x2_minus_210, 2)) ||
                   !check_grh_transaction_refines_to_default_audit(
                           setup_from_coefficients(d12_coefficients, 12)) ||
                   !check_linear_factor_count_cache(
                           setup_from_coefficients(d12_coefficients, 12))
            ? 1
            : 0;
}
