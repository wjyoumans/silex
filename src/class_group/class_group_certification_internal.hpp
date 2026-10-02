#pragma once

#include <silex/class_group.hpp>

namespace silex::detail {

class ClassGroupCertificationAccess {
public:
    static bool exact_imaginary_quadratic_class_order_for_run(
            flint::FmpzRef out,
            ClassGroupContext& context,
            flint::FmpzConstRef discriminant) noexcept;

    static bool try_certify_imaginary_quadratic_from_exact_order(
            ClassGroupContext& context,
            CertificationMode requested,
            flint::FmpzConstRef discriminant,
            flint::FmpzConstRef exact_order) noexcept;

    static bool record_factor_base_honesty_proof(
            ClassGroupContext& context,
            flint::FmpzConstRef required_bound) noexcept;

    // Exact imaginary-quadratic index h_cand / h for a rank-zero maximal
    // order.  Under a proven request, true with index one means the class
    // group and `units` have been published proven (this requires the
    // torsion of `units` to be the computed torsion).  False leaves the
    // certification state as it was, including the factor-base generation
    // check a proven request runs first.  Index > 1, or any index under a
    // grh request, publishes no label and is only a saturation target or
    // acceptance test for the caller; a proven request then keeps its
    // generation check.
    static bool rank_zero_quadratic_class_index_bound(
            flint::FmpzRef out,
            ClassGroupContext& context,
            OrderUnitGroup& units,
            CertificationMode requested) noexcept;

    static bool saturate_relations_for_index_bound_with_units(
            bool& changed,
            bool& saturated,
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept;

    // The gates below trust caller-supplied proof data: an analytic
    // class-regulator product hR or a saturation index bound.  They are
    // internal so that only Silex routes that computed that data themselves
    // (for example from the zeta function) can publish `proven` through them.
    //
    // An analytic hR publishes `proven` only when it is unconditional:
    // degree one, or the quadratic L(1, chi) route.  A Belabas-Friedman hR
    // assumes GRH (Belabas-Friedman 2015, Theorem 1) and never publishes
    // `proven` on its own.  `hr_unconditional` states the route that
    // produced the supplied value (for example
    // `AnalyticClassRegulatorCache::value_unconditional()`); it is honoured
    // only for an order that has an unconditional route at all.  There is
    // deliberately no overload without it: whether a value is unconditional
    // depends on the route that computed it, not on the order, because the
    // quadratic route falls back to Belabas-Friedman when L(1, chi) fails.
    static bool try_certify_class_unit_with_units(
            ClassGroupContext& context,
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            slong precision,
            bool hr_unconditional) noexcept;

    // Records the analytic index-one check that accepted a `grh`-requested
    // class/unit pair.  The record is informational: it sets
    // `analytic_class_regulator_status` to `verified` with the
    // conditionality of the hR that was used (`hr_unconditional` for the
    // quadratic L(1, chi) route, GRH for a Belabas-Friedman hR; degree one
    // takes the exact route and never reaches this record) and never
    // changes the class-group or unit certification labels, the
    // unit/regulator proof states, or any other proof record.  A Belabas-Friedman check recorded here never serves as
    // a proof component.  The caller must have accepted the pair by the
    // analytic index-one test against that hR; the exact imaginary-quadratic
    // route, which uses no analytic value, records nothing.  An
    // unconditional record stored here can later let an explicit
    // try_certify_with_units(proven) promote the class group (see
    // analytic_class_regulator_certification() in class_group.hpp); this
    // call itself never promotes.
    static void record_grh_acceptance_analytic_check(
            ClassGroupContext& context,
            bool hr_unconditional) noexcept;

    // Publishes the grh labels on a completed class/unit pair.  GRH is
    // assumed only by the analytic hR check; the torsion does not depend on
    // it, and a wrong torsion order w would make the grh unit label wrong.
    // So the torsion of `units` must be the torsion Silex computes for the
    // order (order_unit_torsion_is_computed, one rank_zero_torsion);
    // otherwise this fails closed: false, both labels unchanged.
    static bool publish_grh_labels(
            ClassGroupContext& context,
            OrderUnitGroup& units,
            const DiagnosticsContext* diagnostics) noexcept;

    // Adds the Belabas-Friedman audit data (error bound, cutoff, and
    // precisions) of the hR used by record_grh_acceptance_analytic_check()
    // when that hR was the default zeta route's Belabas-Friedman fallback.
    // Like that record, it never changes certification labels.
    static bool record_grh_acceptance_bf_audit(
            ClassGroupContext& context,
            flint::ArbConstRef error_bound,
            ulong cutoff,
            ulong max_cutoff,
            slong requested_precision,
            slong work_precision) noexcept;

    static bool try_analytic_index_bound_with_units(
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            flint::FmpzConstRef aux_prime_bound,
            slong precision,
            bool hr_unconditional) noexcept;

    // Unconditional class-group proof by saturation: with factor-base
    // generation checked and `units` proven, prove relation saturation at
    // every prime dividing the candidate class number (local discrete-log
    // proofs over degree-one auxiliary primes up to `aux_prime_bound`) and
    // publish `proven`.  Fails and leaves the context unchanged otherwise.
    static bool try_prove_class_order_saturation_with_units(
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            flint::FmpzConstRef aux_prime_bound) noexcept;

    static bool try_prove_relation_saturation_index_bound_with_units(
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            flint::FmpzConstRef aux_prime_bound) noexcept;

    // The stored ell-local saturation proof record for `ell`, for tests and
    // audits; false when none is stored.
    static bool relation_saturation_proof_record(
            const ClassGroupContext& context,
            flint::FmpzConstRef ell,
            ProofState& status,
            slong& rank,
            slong& target_rank,
            slong& local_primes) noexcept;

    // A Belabas-Friedman audit records a GRH-conditional analytic check; it
    // publishes `proven` only when saturation has already been proven (or in
    // degree one, where hR = 1 exactly).
    static bool try_certify_class_unit_with_bf_audit(
            ClassGroupContext& context,
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            flint::ArbConstRef error_bound,
            ulong cutoff,
            ulong max_cutoff,
            slong requested_precision,
            slong work_precision) noexcept;

private:
    // Stores a verified relation-saturation record, backed by an
    // exact-class-order proof record, at every prime dividing the verified
    // exact class order `exact_order`.  A prime that already carries a
    // verified ell-local proof keeps it.
    static bool record_exact_class_order_saturation(
            ClassGroupContext& context,
            flint::FmpzConstRef exact_order) noexcept;
};

}  // namespace silex::detail
