#pragma once

#include <silex/class_group.hpp>

#include "../factor_base/generation_different_route_internal.hpp"

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

    // Records GRH-conditional factor-base generation for a base that was
    // just built to contain every prime ideal of norm at most its build
    // bound (the grh policy's norm-bounded base), when that build bound is
    // at least `grh_bound`, the GRH generation bound min(BDF, Bach).  Under
    // GRH such a base generates the class group (Bach 1990, Thm 4 and the
    // remark following it, p. 376; Belabas-Diaz y Diaz-Friedman 2008,
    // Cor 2.2).  The record is separate
    // from the unconditional generation status, which stays as built and is
    // never set to `verified` here; it is written only at build time and
    // never by the honesty receipt.  Returns false, recording nothing, when
    // the build bound is below `grh_bound` or the context has no base.
    // `basis` is the theorem that gave `grh_bound`, bdf or bach; any other
    // value records nothing and returns false.  It is reported by
    // factor_base_generation_basis() while generation is not also
    // unconditional.
    static bool record_grh_factor_base_generation(
            ClassGroupContext& context,
            flint::FmpzConstRef grh_bound,
            ClassGroupGenerationBasis basis) noexcept;

    // True when the factor base carries the GRH generation record above.
    // Only a grh class/unit request writes it; readers that need
    // unconditional generation use factor_base_generation_status().
    static bool grh_factor_base_generation_covered(
            const ClassGroupContext& context) noexcept;

    // The theorem stored with that record, bdf or bach, and none when there
    // is no record.  Unlike factor_base_generation_basis(), it reports the
    // record even when generation is also unconditional.
    static ClassGroupGenerationBasis grh_factor_base_generation_basis(
            const ClassGroupContext& context) noexcept;

    // How the recorded unconditional generation requirement
    // (factor_base_generation_bound()) treats the class of the different:
    // index_one or conductor_primes when it is Zimmert's twin-class bound,
    // none when it is factor_base_class_group_bound or nothing is recorded.
    static GenerationDifferentRoute factor_base_generation_different_route(
            const ClassGroupContext& context) noexcept;

    // Exact imaginary-quadratic index h_cand / h for a rank-zero maximal
    // order, under a proven request only; any other request returns false
    // without listing reduced forms (a grh request takes h from the GRH
    // relation lattice and the analytic index test).  True with index one
    // means the class group and `units` have been published proven (this
    // requires the torsion of `units` to be the computed torsion).  False
    // leaves the certification state as it was, including the factor-base
    // generation check it runs first.  Index > 1 publishes no label and is
    // only a saturation target for the caller, which keeps the generation
    // check.
    static bool rank_zero_quadratic_class_index_bound(
            flint::FmpzRef out,
            ClassGroupContext& context,
            OrderUnitGroup& units,
            CertificationMode requested) noexcept;

    // Bounded saturation at the primes of `index_bound`, a bound on
    // h_cand / h.  `index_bound_unconditional` says whether that bound holds
    // without GRH; when it does not, relation_saturation_status() becomes
    // `verified` only if verified ell-records cover every p | h_cand.
    static bool saturate_relations_for_index_bound_with_units(
            bool& changed,
            bool& saturated,
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            bool index_bound_unconditional,
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
    // route (proven requests only) uses no analytic value.  An
    // unconditional record stored here can later let an explicit
    // try_certify_with_units(proven) promote the class group (see
    // analytic_class_regulator_certification() in class_group.hpp); this
    // call itself never promotes.
    //
    // It also records relation and unit completeness of the accepted pair
    // (class_unit_completeness_certification() in class_group.hpp):
    // `proven` with basis unconditional_analytic for an unconditional hR,
    // `grh` with basis belabas_friedman otherwise (Belabas-Friedman 2015,
    // Theorem 1, which assumes RH for zeta_K and zeta_Q).
    static void record_grh_acceptance_analytic_check(
            ClassGroupContext& context,
            bool hr_unconditional) noexcept;

    // Publishes the grh labels on a completed class/unit pair.  GRH is
    // assumed only by the analytic hR check; the torsion does not depend on
    // it, and a wrong torsion order w would make the grh unit label wrong.
    // So the torsion of `units` must be the torsion Silex computes for the
    // order (order_unit_torsion_is_computed, one rank_zero_torsion);
    // otherwise this fails closed: false, both labels unchanged.  When the
    // class group is already `proven` (the exact degree-one route of a grh
    // request), completeness is recorded `proven` with basis
    // unconditional_certification before the relabel.
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
    // degree one, where hR = 1 exactly).  Class/unit validation, the only
    // caller, uses it only in degree greater than two; the degree-one
    // branch is defensive and matches try_certify_class_unit_with_zeta_bf.
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
    friend class silex::ClassGroupContext;

    // The per-prime part shared by the two saturation proofs above: proves
    // relation saturation at each prime of `required_ells`, in order, by
    // the ell-local discrete-log proof and marks it verified.  With
    // `skip_verified_ells` a prime that already carries a verified
    // ell-record is kept instead.  Then requires the records to cover every
    // prime dividing h_cand and sets the saturation status `verified`.
    //
    // It checks none of its preconditions; the caller must ensure that
    //   - the context has a presentation;
    //   - `units` is set and proven on the same parent order;
    //   - factor-base generation is checked up to the required bound;
    //   - the unit and regulator proof states are verified, or the caller
    //     sets them to verified in the same transaction;
    //   - the caller holds a CertificationTransaction_ and returns on false,
    //     so that the transaction rolls back every change made here.
    // `required_ells` need not contain every prime dividing h_cand for
    // soundness, since the coverage check above rejects a missing prime.
    // A prime p | h_cand left out of `required_ells` makes the call fail
    // unless an earlier proof already left both a verified ell-record and
    // a verified proof record at p, which the coverage check accepts (see
    // relation_saturation_ell_verified_).  Both callers pass every
    // p | h_cand.
    static bool prove_relation_saturation_ells(
            ClassGroupContext& context,
            const OrderUnitGroup& units,
            const std::vector<flint::Fmpz>& required_ells,
            flint::FmpzConstRef aux_prime_bound,
            bool skip_verified_ells) noexcept;

    // Stores a verified relation-saturation record, backed by an
    // exact-class-order proof record, at every prime dividing the verified
    // exact class order `exact_order`.  A prime that already carries a
    // verified ell-local proof keeps it; any other proof record at that
    // prime (failed or not checked) is overwritten as verified.
    static bool record_exact_class_order_saturation(
            ClassGroupContext& context,
            flint::FmpzConstRef exact_order) noexcept;
    // Publishes the class group `proven` from h_cand equal to the exact
    // class order `exact_order`, with generation verified: stores the
    // per-prime saturation records, requires them to cover h_cand, and sets
    // the certification and the saturation, unit, and regulator proof
    // states.  The caller holds the certification transaction and commits
    // it only on success.
    static bool publish_exact_class_order_proof(
            ClassGroupContext& context,
            flint::FmpzConstRef exact_order) noexcept;
};

}  // namespace silex::detail
