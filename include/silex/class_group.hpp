#pragma once

#include <memory>
#include <optional>
#include <vector>

#include <flint/flint.h>

#include <silex/abelian_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/detail/class_relation_module_context.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/fmpz_smat.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>
#include <silex/status.hpp>

namespace silex {

class OrderUnitGroup;
class ClassGroupContext;
class FactoredElement;

namespace detail {
struct ClassUnitTransactionContext;
struct ClassGroupRelationOptions;
struct ClassGroupContextStorage;
struct RelationSaturationCandidateFilter;
class ClassGroupCertificationAccess;
class ClassGroupRelationAccess;
class ClassGroupFinishAccess;
class ClassGroupRelationSearchAccess;
class ClassUnitTransactionAccess;
}  // namespace detail

struct ClassGroupFactorBaseGenerationRecord {
    flint::Fmpz p;
    ProofState status = ProofState::not_checked;
};

struct ClassGroupRelationSaturationRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
};

struct ClassGroupZetaBfProofRecord {
    ulong cutoff = 0;
    ulong max_cutoff = 0;
    slong requested_precision = 0;
    slong work_precision = 0;
    flint::Arb error_bound;
};

// The result that establishes that the factor base generates the class
// group (ClassGroupContext::factor_base_generation_basis()).
// - minkowski_type: unconditional; every prime ideal up to the
//   Minkowski-type bound factor_base_class_group_bound() is in the base or
//   was checked to lie in its span.
// - zimmert_twin_class: unconditional; every prime ideal up to Zimmert's
//   twin-class bound (Zimmert 1981, Korollar 1), which every ideal class K
//   meets in K or in D K^(-1) for the different D, is in the base or was
//   checked to lie in its span.  These primes generate the class group with
//   no hypothesis on the class of D, for any defining polynomial
//   (GB-kor1-generation; see "Class groups and order units" in the
//   algorithms and sources reference).  Used only when it is smaller than
//   the Minkowski-type bound.
// - bdf: under the Riemann hypothesis for the L-functions of the nontrivial
//   characters of the class group (Belabas-Diaz y Diaz-Friedman 2008,
//   Corollary 2.2).
// - bach: under the Riemann hypothesis for zeta_K and those L-functions
//   (Bach 1990, Theorem 4 and the remark following it).
// - structure_certificate: unconditional; imaginary quadratic `proven`
//   class groups only, with no generation bound.  The published invariant
//   generators g_j, products of factor-base primes with invariants d_j,
//   satisfy, in exact reduced binary quadratic form arithmetic: the product
//   of the d_j is the exact class number h, g_j^(d_j) is principal, and for
//   every prime ell dividing h no nonzero F_ell-combination of the
//   g_j^(d_j / ell) with ell | d_j is principal.  The map from the direct sum
//   of the Z/d_j to the class group is then injective between groups of
//   order h, so the g_j, hence the factor base, generate the class group.
enum class ClassGroupGenerationBasis {
    none = 0,
    minkowski_type = 1,
    bdf = 2,
    bach = 3,
    zimmert_twin_class = 4,
    structure_certificate = 5
};

// The result that establishes relation and unit completeness of a published
// class/unit pair (ClassGroupContext::class_unit_completeness_basis()).
// Completeness is relative to factor-base generation: it says the relations
// and units are complete for a generating factor base.
// - belabas_friedman: analytic index one against a Belabas-Friedman hR, under
//   the Riemann hypothesis for zeta_K and zeta_Q (Belabas-Friedman 2015,
//   Theorem 1).
// - unconditional_analytic: analytic index one against an unconditional hR
//   (the quadratic L(1, chi) route).
// - exact_class_number: the candidate class number equals the exact
//   imaginary-quadratic class number (`grh` imaginary-quadratic pairs with
//   |D| < 2^34; above that the pair is accepted by the analytic index-one
//   test).
// - unconditional_certification: the class group and units were certified
//   `proven` (saturation or an unconditional analytic check, see
//   relation_saturation_status() and analytic_class_regulator_certification()),
//   including the exact degree-one route of a `grh` request.  A `proven`
//   class group certified without a unit pair also reports this basis.
enum class ClassUnitCompletenessBasis {
    none = 0,
    belabas_friedman = 1,
    unconditional_analytic = 2,
    exact_class_number = 3,
    unconditional_certification = 4
};

// Resource ceilings for the deterministic standalone relation candidate.
// Successful candidates always retain CertificationMode::unknown.
struct ClassGroupCandidateOptions {
    slong max_candidates = WORD_MAX;
    slong max_relations = WORD_MAX;
    const DiagnosticsContext* diagnostics = nullptr;
};

// Proof resources for the atomic paired class/unit transaction.
struct ClassGroupComputeOptions {
    slong max_candidates = WORD_MAX;
    slong max_relations = WORD_MAX;
    slong relation_saturation_aux_prime_bound = 0;
    slong relation_saturation_max_appends_per_ell = 0;
    ulong zeta_bf_max_cutoff = 0;
    CertificationMode requested_certification = CertificationMode::unknown;
    const DiagnosticsContext* diagnostics = nullptr;
};

enum class ClassGroupRelationSource {
    Unknown = 0,
    Search = 1,
    RandomProduct = 2,
    Supplied = 3,
    Saturation = 4,
    LargePrimeMatch = 5,
    ClassGenerator = 6
};

class ClassGroupContext {
public:
    ClassGroupContext() noexcept;
    explicit ClassGroupContext(const Order& order) noexcept;
    ~ClassGroupContext() noexcept;

    ClassGroupContext(const ClassGroupContext&) = delete;
    ClassGroupContext& operator=(const ClassGroupContext&) = delete;

    ClassGroupContext(ClassGroupContext&& other) noexcept;
    ClassGroupContext& operator=(ClassGroupContext&& other) noexcept;

    void swap(ClassGroupContext& other) noexcept;
    void clear() noexcept;
    bool define(const Order& order) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    void set_diagnostics(const DiagnosticsContext* diagnostics) noexcept;
    const DiagnosticsContext* diagnostics() const noexcept;
    slong analytic_finish_precision() const noexcept;
    bool analytic_finish_product(flint::ArbRef out,
                                      slong& precision) const noexcept;

    bool has_factor_base() const noexcept;
    const FactorBase* factor_base() const noexcept;
    slong generator_count() const noexcept;
    bool factor_base_prime(PrimeIdeal& out, slong index) const noexcept;
    bool factor_base_prime_is_principal(bool& out,
                                        slong index) const noexcept;
    bool factor_base_prime_is_hnf_covered(bool& out,
                                          slong index) noexcept;
    bool build_factor_base(flint::FmpzConstRef bound) noexcept;
    bool set_factor_base(const FactorBase& base) noexcept;

    bool append_relation(const Relation& relation) noexcept;
    bool append_relation(const Relation& relation,
                         ClassGroupRelationSource source) noexcept;
    bool try_append_generator_relation(
            const Element& generator,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept;
    bool try_append_generator_relation(
            bool& partial_throttle_exit,
            const Element& generator,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept;
    slong relation_count() const noexcept;
    slong relation_rank() const noexcept;
    slong skipped_dependent_relation_count() const noexcept;
    bool relation_source(ClassGroupRelationSource& out,
                         slong index) const noexcept;
    slong relation_source_count(ClassGroupRelationSource source) const noexcept;
    bool relations(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::FmpzMat> relations() const noexcept;
    bool relation_generator(Element& out, slong index) const noexcept;
    const Element* relation_generator_at(slong index) const noexcept;

    bool set_relation_matrix(const RelationMatrix& matrix) noexcept;
    bool publish_presentation() noexcept;
    bool has_presentation() const noexcept;
    CertificationMode certification_status() const noexcept;
    ProofState factor_base_generation_status() const noexcept;
    ProofState factor_base_generation_checked_status() const noexcept;
    bool factor_base_build_bound(flint::FmpzRef out) const noexcept;
    bool factor_base_generation_bound(flint::FmpzRef out) const noexcept;
    bool factor_base_generation_checked_bound(
            flint::FmpzRef out) const noexcept;
    slong factor_base_generation_record_count() const noexcept;
    bool factor_base_generation_record(flint::FmpzRef p,
                                       ProofState& status,
                                       slong index) const noexcept;
    std::optional<ClassGroupFactorBaseGenerationRecord>
    factor_base_generation_record(slong index) const noexcept;
    bool check_factor_base_generation_bound(
            flint::FmpzConstRef required_bound) noexcept;
    // Where factor-base generation rests, for the current presentation (`unknown`
    // when there is none, even if a factor base exists):
    // `proven` when factor_base_generation_status() is `verified` (basis
    // minkowski_type, or zimmert_twin_class when the recorded requirement is
    // the twin-class bound), or when a `proven` imaginary quadratic class
    // group was certified by its structure certificate (basis
    // structure_certificate; factor_base_generation_status() then stays as
    // built, since no generation bound was checked), `grh` when a `grh` class/unit request built the base
    // to contain every prime ideal up to the GRH generation bound
    // min(BDF, Bach), or Bach's bound alone when the BDF criterion cannot be
    // computed or decided, and generation is not also unconditional (basis
    // bdf or bach, whichever theorem gave that bound), and `unknown`
    // otherwise (basis none).  GRH generation never changes
    // factor_base_generation_status(), which keeps meaning unconditional
    // generation, and never counts toward a `proven` label.
    CertificationMode factor_base_generation_certification() const noexcept;
    ClassGroupGenerationBasis factor_base_generation_basis() const noexcept;
    // The bound of that theorem: the recorded requirement
    // (factor_base_generation_bound(), the Minkowski-type or twin-class
    // bound) for `proven`, the GRH generation
    // bound before the factor-base size doubling for `grh`.  False when the
    // certification is `unknown` or the basis is structure_certificate,
    // which has no bound.
    bool factor_base_generation_certification_bound(
            flint::FmpzRef out) const noexcept;
    // Where relation and unit completeness rests for the published pair
    // (or for a `proven` class group alone, which also reports `proven` /
    // unconditional_certification),
    // relative to factor-base generation: `proven` for a `proven` class
    // group (basis unconditional_certification) and for a `grh` pair
    // accepted through the exact imaginary-quadratic class number (|D| <
    // 2^34), an unconditional quadratic L(1, chi) hR (a `proven` request),
    // or the exact degree-one route;
    // `grh` for a `grh` pair accepted through a Belabas-Friedman hR (every
    // degree >= 3 field, a degree-two field outside the explicit quadratic
    // backend, every quadratic field of a `grh` request (including
    // x^2 - d) whose acceptance is not the exact route, and a quadratic field
    // with |D| >= 2^34 or whose L(1, chi) evaluation fell back to
    // Belabas-Friedman); `unknown` otherwise.  A
    // `grh` result whose generation and completeness are both `proven` keeps
    // the `grh` label.
    CertificationMode class_unit_completeness_certification() const noexcept;
    ClassUnitCompletenessBasis class_unit_completeness_basis() const noexcept;
    // Whether the relations of the current presentation are known to be
    // saturated.  The routes that set `verified` do not all mean the same
    // coverage:
    // - the saturation proofs (try_prove_relation_saturation_with_units and
    //   the saturation route of a `proven` try_certify_with_units), the
    //   exact imaginary-quadratic route, the internal index-bound proof
    //   try_prove_relation_saturation_index_bound_with_units_, and the
    //   rank-zero exact index-one route
    //   (rank_zero_quadratic_class_index_bound) set it only once every
    //   prime p | h_cand has a verified relation_saturation_record(); the
    //   last two require relation_saturation_covers_class_order_();
    // - try_certify_trivial_quotient, and try_certify_quadratic with
    //   h_cand = 1, set it with no prime required and no record kept;
    // - the unconditional analytic index-one routes (degree one or the
    //   quadratic L(1, chi) route) of a `proven` try_certify_with_units,
    //   including its promotion from an unconditional check stored
    //   earlier (a `grh` run on x^2 - d no longer stores one) and the internal
    //   fallback try_analytic_index_bound_with_units_ it runs when the
    //   direct index-one check fails, of try_certify_class_unit_with_zeta,
    //   and, in degree one only (where hR = 1 exactly), of
    //   try_certify_class_unit_with_zeta_bf and the internal
    //   ClassGroupCertificationAccess::try_certify_class_unit_with_bf_audit
    //   set it once they certify `proven`, because index one with
    //   generation checked means saturation at every prime.  They clear
    //   the per-prime records unless verified records already cover every
    //   p | h_cand, which are kept, so the status can be `verified` with no
    //   records;
    // - saturate_relations_bounded_with_units sets `verified` when its
    //   bounded search removed every prime from the candidate index and
    //   `unavailable` otherwise; its records describe that search, not an
    //   ell-local proof.  Its bound is h_cand itself, which bounds the index
    //   of the computed relation lattice in the full one unconditionally
    //   (that index equals h_cand / h only once generation is known).
    //   Internal class/unit validation uses the
    //   index-bounded variant
    //   saturate_relations_bounded_for_index_with_units_ (the public
    //   method forwards to it), which removes only the primes of a supplied
    //   bound on h_cand / h.  It sets `verified` once those primes are
    //   removed and either the bound holds unconditionally or verified
    //   records cover every p | h_cand, and `unavailable` otherwise.  The
    //   pass keeps records only for primes of the bound, so coverage holds
    //   only once the resulting h_cand is 1 and never rescues a conditional
    //   bound.  The
    //   validation pass supplies the exact imaginary-quadratic index (under
    //   a `grh` request only for |D| < 2^34), which is unconditional when
    //   generation is verified up to the Minkowski-type bound (always so
    //   under a `proven` request) and GRH-conditional otherwise;
    //   a Belabas-Friedman analytic index bound only decides whether the
    //   pass runs and is never the bound it saturates against.
    // The bounded search leaves the status unchanged when the candidate
    // index is 1 or has no prime factors (it returns early).
    // No route publishes `proven` from this status: each `proven` route
    // rests on verified records for every p | h_cand, on h_cand = 1, on the
    // exact class number, or on an unconditional analytic check.
    ProofState relation_saturation_status() const noexcept;
    // `verified` once an analytic class-regulator check (hR from the zeta
    // function) has passed for the current presentation.  The check is
    // unconditional only for degree one and for the quadratic L(1, chi)
    // route of a `proven` request; a Belabas-Friedman hR (every degree >= 3
    // field, a degree-two field outside the explicit quadratic backend, every
    // quadratic field of a `grh` request including x^2 - d, and a quadratic
    // field with |D| >= 2^34 or whose L(1, chi) evaluation fails and falls
    // back to Belabas-Friedman) assumes GRH.  Use
    // analytic_class_regulator_certification() to read which one was used.
    // A `grh` class/unit transaction accepted by the analytic index-one test
    // records that check here as well; the record is
    // informational and leaves the `grh` labels unchanged.  With a
    // Belabas-Friedman hR it also records that evaluation's audit data
    // (zeta_bf_proof_record()).  The exact imaginary-quadratic (|D| < 2^34)
    // and degree-one `grh` routes use no analytic check and record none.
    ProofState analytic_class_regulator_status() const noexcept;
    // `proven` for an unconditional analytic check, `grh` for a
    // Belabas-Friedman (GRH-conditional) check, and `unknown` when no
    // analytic check is recorded.  Once an unconditional check is recorded,
    // a later Belabas-Friedman audit of the same presentation keeps the
    // `proven` label.  This is the conditionality of the analytic check
    // only, not of the result: a `grh` transaction accepted through an
    // unconditional hR (quadratic L(1, chi), no longer taken by a `grh`
    // request) reports `proven` here while certification_status() stays
    // `grh`.  The `grh` transaction itself
    // never promotes a label.  The one exception to "the record never
    // promotes" is a later explicit try_certify_with_units(units, proven)
    // with proven `units`: a stored `proven` record then counts as the
    // unconditional analytic check and can promote the class group to
    // `proven`.  A `grh` run on x^2 - d stores a Belabas-Friedman record,
    // not a `proven` one, so for it the call computes L(1, chi) itself at a
    // usable precision.  A GRH-conditional check never
    // serves as a proof component and never makes the class group or the
    // units `proven`: a `proven` label then requires factor-base
    // generation, proven units, and verified relation saturation at every
    // prime dividing the candidate class number.
    CertificationMode analytic_class_regulator_certification() const noexcept;
    // Audit of a Belabas-Friedman hR evaluation.  Its error bound is valid
    // under GRH for zeta_K and zeta_Q (Belabas-Friedman 2015, Theorem 1 or
    // Theorem 7), so it is a GRH-conditional record, not an unconditional
    // proof.  A record stored by an accepted `grh` class/unit transaction
    // comes from the coarse index-one evaluation (Theorem 7 error bound, cutoff
    // chosen so that the absolute error of log(hR) is below 1/4, enough to
    // decide h_cand R_cand/hR < 2) only when that coarse value itself shows
    // index one; any other outcome (no cutoff up to 20007, or a coarse value
    // that does not show index one) is recomputed and recorded from the default
    // zeta route.  Such a coarse value is then discarded, so a field can pay
    // up to part of one extra analytic evaluation.  A field whose default
    // route fails can be accepted when the coarse value shows index one.  Its error bound need not meet the target for the requested
    // precision (unlike `try_certify_class_unit_with_zeta_bf`).
    ProofState zeta_bf_proof_status() const noexcept;
    bool zeta_bf_proof_record(ulong& cutoff,
                              ulong& max_cutoff,
                              slong& requested_precision,
                              slong& work_precision,
                              flint::ArbRef error_bound) const noexcept;
    std::optional<ClassGroupZetaBfProofRecord>
    zeta_bf_proof_record() const noexcept;
    slong relation_saturation_record_count() const noexcept;
    bool relation_saturation_record(flint::FmpzRef ell,
                                    ProofState& status,
                                    slong index) const noexcept;
    std::optional<ClassGroupRelationSaturationRecord>
    relation_saturation_record(slong index) const noexcept;
    ProofState unit_proof_status() const noexcept;
    ProofState regulator_proof_status() const noexcept;
    bool try_certify_quadratic(CertificationMode requested) noexcept;
    bool try_certify_trivial_quotient(CertificationMode requested) noexcept;
    // A `proven` request succeeds through factor-base generation and proven
    // units plus either an unconditional analytic hR (degree one or the
    // quadratic L(1, chi) route) or verified relation saturation at every
    // prime dividing the candidate class number.  A Belabas-Friedman hR,
    // which assumes GRH, is never used to publish `proven`, so
    // `zeta_bf_max_cutoff` is ignored for a `proven` request (and no other
    // request evaluates hR); the overload is kept for source compatibility.
    // On degree >= 3 a `proven` request runs the saturation proof, which
    // computes discrete logarithms at auxiliary primes for every prime
    // dividing the candidate class number.
    bool try_certify_with_units(const OrderUnitGroup& units,
                                CertificationMode requested,
                                slong precision) noexcept;
    bool try_certify_with_units(const OrderUnitGroup& units,
                                CertificationMode requested,
                                slong precision,
                                ulong zeta_bf_max_cutoff) noexcept;
    // Publish `proven` for the class group and units from an analytic hR
    // only when that hR is unconditional (degree one or the quadratic
    // L(1, chi) route).  In degree one the Belabas-Friedman entry point uses
    // the exact value hR = 1 and evaluates no series.  With a
    // Belabas-Friedman hR in any other degree, which assumes GRH, they
    // succeed only if `units` is already proven and the class group already
    // has an unconditional proof component for this presentation: relation
    // saturation proven at every prime dividing the candidate class number, or
    // an unconditional analytic check recorded earlier (degree one or the
    // quadratic L(1, chi) route, which a `grh` run on x^2 - d no longer
    // stores); then they only record the GRH-conditional analytic check.
    // Otherwise they fail and leave both objects unchanged.
    bool try_certify_class_unit_with_zeta(OrderUnitGroup& units,
                                          slong precision) noexcept;
    bool try_certify_class_unit_with_zeta_bf(OrderUnitGroup& units,
                                             ulong max_cutoff,
                                             slong precision) noexcept;
    bool try_prove_relation_saturation_with_units(
            const OrderUnitGroup& units,
            flint::FmpzConstRef ell,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool saturate_relations_bounded_with_units(
            bool& changed,
            bool& saturated,
            const OrderUnitGroup& units,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept;
    bool presentation(FiniteAbelianGroup& out) const noexcept;
    std::optional<FiniteAbelianGroup> presentation() const noexcept;
    slong invariant_count() const noexcept;
    bool invariant(flint::FmpzRef out, slong index) const noexcept;
    bool invariants(flint::FmpzVecRef out) const noexcept;
    bool order(flint::FmpzRef out) const noexcept;
    bool invariant_generator_matrix(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::Fmpz> invariant(slong index) const noexcept;
    std::optional<flint::Fmpz> order() const noexcept;
    std::optional<flint::FmpzMat> invariant_generator_matrix() const noexcept;
    bool invariant_generator(FractionalIdeal& out, slong index) const noexcept;
    bool ideal_class_coordinates(flint::FmpzMatRef out,
                                 const FractionalIdeal& ideal) const noexcept;
    bool invariant_generator_power_witness(FactoredElement& out,
                                           slong index) const noexcept;
    slong relation_kernel_unit_count() const noexcept;
    bool relation_kernel_unit(FactoredElement& out, slong index) const noexcept;
    bool compute_candidate(const Order& order,
                           flint::FmpzConstRef factor_base_bound,
                           const ClassGroupCandidateOptions& options = {})
            noexcept;

private:
    friend class OrderUnitGroup;
    friend class detail::ClassGroupCertificationAccess;
    friend class detail::ClassGroupRelationAccess;
    friend class detail::ClassGroupFinishAccess;
    friend class detail::ClassGroupRelationSearchAccess;
    friend class detail::ClassUnitTransactionAccess;

    enum class RelationAppendOutcome {
        none,
        rank,
        index,
        kernel,
        skipped_dependent
    };

    enum class DependentRelationPolicy {
        native,
        skip,
        keep,
        keep_nonduplicate
    };

    enum class SpecializedRelationBackendStatus {
        unavailable,
        succeeded,
        failed
    };

    // Gates that trust caller-supplied proof data (an analytic hR, an index
    // bound, or an exact norm).  They are reachable only through the
    // noninstalled detail access classes; the public API certifies only from
    // data it computes itself.
    bool try_append_generator_relation_with_norm_(
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpqConstRef norm,
            ClassGroupRelationSource source) noexcept;
    bool try_certify_class_unit_with_units_(
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            slong precision,
            bool hr_unconditional) noexcept;
    bool try_prove_relation_saturation_index_bound_with_units_(
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool try_analytic_index_bound_with_units_(
            const OrderUnitGroup& units,
            flint::ArbConstRef analytic_class_regulator_product,
            flint::FmpzConstRef aux_prime_bound,
            slong precision,
            bool hr_unconditional) noexcept;
    bool try_prove_class_order_saturation_with_units_(
            const OrderUnitGroup& units,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool prove_relation_saturation_dlog_ell_(
            const OrderUnitGroup& units,
            flint::FmpzConstRef ell,
            flint::FmpzConstRef aux_prime_bound) noexcept;
    bool try_promote_proven_certification_() noexcept;
    bool try_certify_analytic_class_regulator_(
            const OrderUnitGroup& units,
            flint::ArbConstRef analytic_hR,
            slong precision,
            bool hr_unconditional) noexcept;
    bool record_analytic_class_unit_regulator_(
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_hR,
            slong precision,
            bool hr_unconditional) noexcept;
    bool try_certify_analytic_class_unit_regulator_(
            OrderUnitGroup& units,
            flint::ArbConstRef analytic_hR,
            slong precision,
            bool hr_unconditional) noexcept;
    void record_analytic_class_regulator_check_(bool unconditional) noexcept;
    bool record_zeta_bf_audit_(flint::ArbConstRef error_bound,
                               ulong cutoff,
                               ulong max_cutoff,
                               slong requested_precision,
                               slong work_precision) noexcept;
    bool relation_row_refines_(const Relation& relation) noexcept;
    bool sync_row_module_checkpoint_() noexcept;
    bool append_relation_with_outcome_(RelationAppendOutcome& outcome,
                                       const Relation& relation,
                                       ClassGroupRelationSource source)
            noexcept;
    bool append_relation_with_outcome_(
            RelationAppendOutcome& outcome,
            const Relation& relation,
            ClassGroupRelationSource source,
            DependentRelationPolicy dependent_policy) noexcept;
    SpecializedRelationBackendStatus
    run_maximal_imaginary_quadratic_relation_backend_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool build_relation_factor_base_(flint::FmpzConstRef bound) noexcept;
    bool build_search_factor_base_(
            flint::FmpzConstRef bound,
            bool incomplete = false,
            FactorBase* norm_bounded_base = nullptr) noexcept;
    void reset_partial_relations_() noexcept;
    void configure_partial_relations_(
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool try_append_integral_generator_relation_(
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            ClassGroupRelationSource source) noexcept;
    bool try_partial_relation_(bool& partial_throttle_exit,
                               const Element& generator,
                               const fmpq* known_norm = nullptr,
                               const fmpz_mat_struct*
                                       known_integral_coordinates = nullptr,
                               const fmpz_poly_struct*
                                       known_integral_polynomial = nullptr)
            noexcept;
    bool saturate_local_once_(bool& changed,
                              const PrimeIdeal& prime,
                              flint::FmpzConstRef ell) noexcept;
    bool saturate_local_(bool& changed,
                         bool& index_cleared,
                         flint::FmpzConstRef ell,
                         flint::FmpzConstRef aux_prime_bound,
                         slong max_appends) noexcept;
    bool saturate_relations_bounded_(bool& changed,
                                     bool& saturated,
                                     flint::FmpzConstRef aux_prime_bound,
                                     slong max_appends_per_ell,
                                     slong max_appends_total) noexcept;
    bool try_auto_relation_saturation_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool* changed_out = nullptr) noexcept;
    bool run_native_experimental_relation_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event) noexcept;
    bool run_relation_search_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event,
            bool clamp_relation_kernel_units_to_rank) noexcept;
    bool run_relation_production_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event) noexcept;
    bool run_relation_production_prepass_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event) noexcept;
    bool run_lll_relation_route_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options,
            bool emit_norm_prefilter_profile_event,
            bool accept_tentative_presentation = false,
            bool extension_only = false) noexcept;
    bool pivot_info_(flint::FmpzRef h,
                           std::vector<slong>& pivots) noexcept;
    bool compute_tentative_candidate_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool extend_tentative_relations_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool extend_lll_relation_slice_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool saturate_local_once_with_units_(
            bool& changed,
            const OrderUnitGroup& units,
            bool include_torsion,
            const PrimeIdeal& prime,
            flint::FmpzConstRef ell,
            detail::RelationSaturationCandidateFilter* candidate_filter)
            noexcept;
    bool saturate_local_with_units_(bool& changed,
                                    bool& index_cleared,
                                    const OrderUnitGroup& units,
                                    flint::FmpzConstRef ell,
                                    flint::FmpzConstRef aux_prime_bound,
                                    slong max_appends) noexcept;
    bool saturate_relations_bounded_for_index_with_units_(
            bool& changed,
            bool& saturated,
            const OrderUnitGroup& units,
            flint::FmpzConstRef index_bound,
            bool index_bound_unconditional,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept;
    bool extend_relation_kernel_units_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;
    bool compute_relation_candidate_(
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const detail::ClassGroupRelationOptions& options) noexcept;

    void reset_certification_metadata_() noexcept;
    void reset_factor_base_generation_check_() noexcept;
    bool mark_factor_base_generation_check_(
            flint::FmpzConstRef checked_bound,
            ProofState status) noexcept;
    bool record_factor_base_generation_(
            flint::FmpzConstRef build_bound) noexcept;
    bool record_factor_base_honesty_proof_(
            flint::FmpzConstRef required_bound) noexcept;
    bool relation_saturation_proof_prereqs_verified_() const noexcept;
    bool record_relation_saturation_proof_(flint::FmpzConstRef ell,
                                           ProofState status,
                                           slong rank,
                                           slong target_rank,
                                           slong local_primes) noexcept;
    bool relation_saturation_proof_verified_(
            flint::FmpzConstRef ell) const noexcept;
    bool mark_relation_saturation_verified_(
            flint::FmpzConstRef ell) noexcept;
    bool mark_relation_saturation_(flint::FmpzConstRef ell,
                                   ProofState status) noexcept;
    bool relation_saturation_ell_verified_(
            flint::FmpzConstRef ell) const noexcept;
    bool relation_saturation_covers_class_order_() const noexcept;
    bool relation_saturation_proof_complete_() const noexcept;
    struct CertificationSnapshot_;
    class CertificationTransaction_;
    void save_certification_state_(
            CertificationSnapshot_& out) const noexcept;
    void restore_certification_state_(
            CertificationSnapshot_& saved) noexcept;
    bool ensure_private_storage_() noexcept;
    Order parent_;
    FactorBase base_;
    RelationMatrix relations_;
    std::vector<ClassGroupRelationSource> relation_sources_;
    std::vector<char> relation_basis_flags_;
    Relation generator_relation_scratch_;
    // Allocated transactionally by define(); null only while undefined.
    std::unique_ptr<detail::ClassGroupContextStorage> private_storage_;
    FiniteAbelianGroup quotient_;
    fmpz_smat::HnfContext row_module_;
    detail::ClassRelationModuleContext relation_module_;
    slong relation_rank_ = 0;
    bool row_module_synced_ = true;
    slong skipped_dependent_relations_ = 0;
    slong relation_kernel_units_target_ = 0;
    DependentRelationPolicy generator_relation_policy_ =
            DependentRelationPolicy::native;
    CertificationMode certification_ = CertificationMode::unknown;
    flint::Fmpz factor_base_build_bound_;
    flint::Fmpz factor_base_generation_bound_;
    flint::Fmpz factor_base_generation_checked_bound_;
    ProofState factor_base_generation_status_ = ProofState::not_checked;
    ProofState factor_base_generation_checked_status_ =
            ProofState::not_checked;
    ProofState relation_saturation_status_ = ProofState::not_checked;
    ProofState analytic_class_regulator_status_ = ProofState::not_checked;
    // True when the recorded analytic check used a Belabas-Friedman hR,
    // whose error bound assumes GRH.
    bool analytic_class_regulator_assumes_grh_ = false;
    ProofState zeta_bf_status_ = ProofState::not_checked;
    ulong zeta_bf_cutoff_ = 0;
    ulong zeta_bf_max_cutoff_ = 0;
    slong zeta_bf_requested_precision_ = 0;
    slong zeta_bf_work_precision_ = 0;
    flint::Arb zeta_bf_error_bound_;
    const DiagnosticsContext* diagnostics_ = nullptr;
    // Driver-scoped and non-owning; intentionally not swapped with results.
    detail::ClassUnitTransactionContext* class_unit_transaction_context_ =
            nullptr;
    ProofState unit_proof_status_ = ProofState::not_checked;
    ProofState regulator_proof_status_ = ProofState::not_checked;
    bool has_base_ = false;
};

inline void swap(ClassGroupContext& left,
                 ClassGroupContext& right) noexcept {
    left.swap(right);
}

}  // namespace silex
