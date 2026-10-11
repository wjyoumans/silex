#pragma once

#include <silex/class_group.hpp>
#include <silex/element.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_field.hpp>

#include "relation_completion_scheduler_internal.hpp"
#include "../factor_base/generation_different_route_internal.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace silex::detail {

struct FactorBaseHonestyFailure;

struct PartialRelationEntry {
    explicit PartialRelationEntry(const Order& order) noexcept
        : prime(order),
          generator(*order.parent()) {
    }

    PrimeIdeal prime;
    Element generator;
    flint::FmpzMat row{0, 0};
    flint::Fmpz rational_prime;
    flint::FmpzPoly residue_key;
    bool has_residue_key = false;
    bool has_row = false;
};

struct DirectPartialResidueBlockCacheEntry {
    DirectPartialResidueBlockCacheEntry(
            slong block_index,
            flint::FmpzConstRef rational_prime) noexcept
        : context(rational_prime.raw()),
          input(context),
          remainder(context),
          block_index(block_index) {
        flint::fmpz_set(flint::FmpzRef(this->rational_prime),
                        rational_prime);
    }

    flint::Fmpz rational_prime;
    flint::FmpzModCtx context;
    flint::FmpzModPoly input;
    flint::FmpzModPoly remainder;
    std::vector<flint::FmpzModPoly> residue_polynomials;
    std::vector<slong> prime_indices;
    std::vector<slong> residue_degrees;
    slong block_index = -1;
};

// Test hook: runs the direct large-prime partial-relation factorization
// (the residue-polynomial branch of reference Rel_add.jl `_factor!`) on
// `generator` over `base` with a fresh residue-block cache.  Returns true
// and writes the factor-base row and leftover rational prime when that
// branch finds a single large prime; returns false otherwise.  When given,
// `known_norm` must be the exact field norm of `generator` and
// `known_integral_coordinates` its maximal-order coordinates.  For tests only.
bool direct_partial_relation_row_for_testing(
        flint::FmpzMat& out_row,
        flint::Fmpz& out_rational_large_prime,
        const Element& generator,
        const FactorBase& base,
        flint::FmpzConstRef factor_base_bound,
        const fmpq* known_norm = nullptr,
        const fmpz_mat_struct* known_integral_coordinates = nullptr) noexcept;

struct FactorBaseGenerationRecord {
    flint::Fmpz p;
    ProofState status = ProofState::not_checked;
};

struct RelationSaturationRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
};

// What a relation-saturation proof record at `ell` rests on.
//
// - `ell_local_dlog`: the ell-local test of
//   ClassGroupContext::prove_relation_saturation_dlog_ell_; `rank`,
//   `target_rank`, and `local_primes` describe that computation.  When the
//   local target rank is 0 the proof is vacuous and no computation ran;
//   the record then stores zero for all three.
// - `exact_class_order`: the exact imaginary-quadratic route.  With
//   generation verified and h_cand equal to the exact class number h, the
//   computed relation lattice has index h_cand / h = 1 in the full one, so it
//   is saturated at every prime; `rank`, `target_rank`, and `local_primes`
//   are zero because no ell-local computation ran.
enum class RelationSaturationProofBasis : unsigned char {
    ell_local_dlog,
    exact_class_order,
};

struct RelationSaturationProofRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
    slong rank = 0;
    slong target_rank = 0;
    slong local_primes = 0;
    RelationSaturationProofBasis basis =
            RelationSaturationProofBasis::ell_local_dlog;
};

struct RelationAdmissionCache {
    slong target_relation_count = 0;
    slong relation_surplus = 0;
    slong relation_count_before_init = 0;
    slong retained_relation_count = 0;
    slong trivial_relation_count = 0;
    slong duplicate_relation_count = 0;
    slong skipped_dependent_relation_count = 0;
    slong missing = 0;
    std::vector<ulong> modular_basis;
    std::vector<ulong> relation_row_hashes;
};

struct RelationSearchControlState {
    RelationCompletionState completion;
    slong checkpoint_relation_count = -1;
    slong done_small = 0;
    slong small_fail = 0;
    slong fail_limit = 0;
    slong candidates_tried = 0;
    slong accepted_relations = 0;
    slong checkpoint_extra_relation_requests = 0;
    slong analytic_extra_relation_requests = 0;
    slong analytic_finish_precision = 160;
    slong analytic_finish_product_precision = 0;
    slong analytic_precision_doublings = 0;
    slong analytic_precision_inconclusive_count = 0;
    slong bounded_full_factor_base_restart_requests = 0;
    slong honesty_restart_requests = 0;
    slong squash_index = 0;
    slong finish_unit_log_rotation = 0;
    bool bounded_full_factor_base_restart_pending = false;
    bool bounded_full_factor_base_restart_allow_past_half = false;
    bool finish_unit_log_rotation_active = false;
    bool finish_full_rank_relation_active = false;
    bool honesty_checked = false;
    bool analytic_finish_product_valid = false;
    // reference relation_completion_parameters retains nfrootsof1 from compute_inverse_residue through
    // fundamental-unit publication. Keep the same private finish lifetime.
    bool analytic_finish_torsion_valid = false;
    flint::Arb analytic_finish_product;
    flint::Fmpz analytic_finish_torsion_order;
    Element analytic_finish_torsion_generator;
};

struct RelationFinishState {
    RelationAdmissionCache cache;
    SubfactorBaseSchedule subfactor_base_schedule;
    RelationSearchControlState route_state;
    std::vector<slong> small_multiplier;
    flint::Fmpz factor_base_bound;
    slong target_relation_kernel_units = 0;
};

struct HnfFinishWorkspace {
    flint::FmpzMat hnf{0, 0};
    flint::FmpzMat relation_transform{0, 0};
    flint::ArbMat relation_logs{0, 0};
    flint::ArbMat unit_logs{0, 0};
    slong generator_count = 0;
    slong processed_relation_count = 0;
    slong relation_rank = 0;
    slong relation_log_precision = 0;
    slong logged_relation_count = 0;
    slong unit_log_precision = 0;
    bool exact_valid = false;
    bool relation_logs_valid = false;
    bool unit_logs_valid = false;

    void clear() noexcept {
        hnf = flint::FmpzMat(0, 0);
        relation_transform = flint::FmpzMat(0, 0);
        relation_logs = flint::ArbMat(0, 0);
        unit_logs = flint::ArbMat(0, 0);
        generator_count = 0;
        processed_relation_count = 0;
        relation_rank = 0;
        relation_log_precision = 0;
        logged_relation_count = 0;
        unit_log_precision = 0;
        exact_valid = false;
        relation_logs_valid = false;
        unit_logs_valid = false;
    }
};

struct RandomIdealSearchState {
    std::vector<slong> base_indices;
    std::vector<slong> exponents;
    std::vector<Ideal> base_ideals;
    std::vector<Ideal> inverse_base_numerators;
    std::vector<flint::Fmpz> inverse_base_denominators;
    Ideal random_ideal;
    flint::Fmpz lower_bound;
    flint::Fmpz upper_bound;
    ulong random_state = UWORD(0);
    bool defined = false;
};

struct ReducedIdealLatticeCacheEntry {
    ReducedIdealLatticeCacheEntry(flint::FmpzMatConstRef ideal_hnf,
                                  flint::FmpzMatConstRef reduced_basis)
        : ideal_hnf(flint::fmpz_mat_nrows(ideal_hnf),
                    flint::fmpz_mat_ncols(ideal_hnf)),
          reduced_basis(flint::fmpz_mat_nrows(reduced_basis),
                        flint::fmpz_mat_ncols(reduced_basis)) {
        flint::fmpz_mat_set(this->ideal_hnf, ideal_hnf);
        flint::fmpz_mat_set(this->reduced_basis, reduced_basis);
    }

    flint::FmpzMat ideal_hnf;
    flint::FmpzMat reduced_basis;
};

struct ReducedIdealLatticeCache {
    // Match the existing source-backed post-finite total candidate budget.
    static constexpr std::size_t capacity = 64;

    std::vector<ReducedIdealLatticeCacheEntry> entries;
    std::size_t next_eviction = 0;

    void clear() noexcept {
        entries.clear();
        next_eviction = 0;
    }
};

// Residue-field quotient logarithms reused by local class-relation
// saturation with units (`saturate_local_once_with_units_`).  Every value is
// the exact `ResidueFieldQuotientLog::apply` image of one number-field
// element at one (prime, ell), so it stays valid while that element is
// unchanged.  The saturation pass checks the relation generators and unit
// factorizations against `elements` by value on every use and clears the
// cache on any mismatch; the cache never decides which primes or kernel rows
// are tried.
struct LocalSaturationDlogEntry {
    enum class State : unsigned char {
        unknown,
        unit,          // nonzero residue; `logs[i]` holds its quotient log
        zero,          // the residue is zero
        unreducible,   // the element does not reduce modulo the prime
        failed,        // the quotient log of a nonzero residue failed
    };

    flint::Fmpz ell;
    ResidueFieldQuotientLog quotient_log;
    ResidueFieldElement image;
    // Indexed like LocalSaturationDlogCache::elements.
    std::vector<State> states;
    std::vector<flint::Fmpz> logs;
};

struct LocalSaturationDlogCache {
    struct UnitTerm {
        slong element = 0;
        slong exponent = 0;
    };

    // Distinct relation generators and unit-factor bases.
    std::vector<Element> elements;
    std::vector<std::uint64_t> element_hashes;
    // Relation generator i is elements[relation_elements[i]].
    std::vector<slong> relation_elements;
    // Free unit generator u is the product of elements[t.element]^t.exponent
    // over unit_terms[u], factor by factor in its stored order.
    std::vector<std::vector<UnitTerm>> unit_terms;
    bool units_recorded = false;
    std::vector<LocalSaturationDlogEntry> entries;

    void clear() noexcept {
        elements.clear();
        element_hashes.clear();
        relation_elements.clear();
        unit_terms.clear();
        units_recorded = false;
        entries.clear();
    }
};

struct ClassGroupContextStorage {
    std::vector<PartialRelationEntry> partial_relations;
    std::vector<DirectPartialResidueBlockCacheEntry>
            direct_partial_residue_blocks;
    slong partial_relations_max = 0;
    slong partial_relations_nonproductive_streak = 0;
    bool use_partial_relations = false;
    bool partial_relations_configured = false;

    RandomIdealSearchState random_ideal_search_state;
    ReducedIdealLatticeCache native_post_finite_reduced_lattices;
    std::optional<RelationFinishState> finish_state;
    HnfFinishWorkspace hnf_finish_workspace;
    slong analytic_finish_precision = 0;
    flint::Arb analytic_finish_product;
    slong analytic_finish_product_precision = 0;
    bool analytic_finish_product_valid = false;

    std::vector<FactorBaseGenerationRecord> factor_base_generation_records;
    // GRH-conditional factor-base generation (grh class/unit requests only):
    // set when the factor base was built to contain every prime ideal of
    // norm at most `grh_generation_bound`, the GRH generation bound
    // min(BDF, Bach).  Kept apart from the unconditional generation status,
    // which it never sets.  Written only when the base is built, like
    // factor_base_generation_status_, so certification transactions do not
    // restore it.
    flint::Fmpz grh_generation_bound;
    bool grh_generation_covered = false;
    // The theorem that gave grh_generation_bound (bdf or bach).
    ClassGroupGenerationBasis grh_generation_basis =
            ClassGroupGenerationBasis::none;
    // How the unconditional generation requirement (the recorded
    // factor_base_generation_bound_, from proven_generation_requirement)
    // treats the class of the different: none for the
    // factor_base_class_group_bound requirement, kor1_generation for
    // Zimmert's twin-class requirement.
    // Written with that bound when the base is built, so certification
    // transactions do not restore it.
    GenerationDifferentRoute generation_different_route =
            GenerationDifferentRoute::none;
    // Relation and unit completeness of a `grh`-requested pair, written
    // when the pair is accepted (record_grh_acceptance_analytic_check,
    // record_grh_acceptance_exact_class_number) or relabelled
    // (publish_grh_labels).  Read only while the published label is `grh`;
    // cleared with the other certification metadata and restored by
    // certification transactions.
    CertificationMode class_unit_completeness = CertificationMode::unknown;
    ClassUnitCompletenessBasis class_unit_completeness_basis =
            ClassUnitCompletenessBasis::none;
    // Set when the published imaginary quadratic class group was certified
    // `proven` by its structure certificate, which proves that the factor
    // base generates the class group without a generation bound; the
    // generation status stays as built.  Read by
    // factor_base_generation_certification() (basis structure_certificate).
    // It describes the published presentation, so it is cleared with the
    // other certification metadata and restored by certification
    // transactions.
    bool generation_structure_certified = false;
    std::vector<RelationSaturationRecord> relation_saturation_records;
    std::vector<RelationSaturationProofRecord>
            relation_saturation_proof_records;
    LocalSaturationDlogCache local_saturation_dlog_cache;
    // LLL-reduced basis of the relation kernel of the published
    // presentation, read by ClassGroupContext::relation_kernel_unit.
    // Rebuilt on first use after the presentation changes; every write of
    // the presentation resets the flag.
    flint::FmpzMat reduced_relation_kernel{0, 0};
    bool reduced_relation_kernel_valid = false;
    // Integral Gram-Schmidt data of reduced_relation_kernel (Cohen, GTM 138,
    // Algorithm 2.6.7): a 1 x (m + 1) row d with d_0 = 1 and
    // d_{j+1} = d_j |b*_j|^2, and an m x m matrix lambda with
    // lambda(k, j) = d_{j+1} mu(k, j) for j < k.  Read by the
    // relation-saturation beta size reduction; built on first use and reset
    // together with the reduced kernel.
    flint::FmpzMat reduced_relation_kernel_gram_d{0, 0};
    flint::FmpzMat reduced_relation_kernel_gram_lambda{0, 0};
    bool reduced_relation_kernel_gram_valid = false;

    void invalidate_reduced_relation_kernel() noexcept {
        reduced_relation_kernel = flint::FmpzMat(0, 0);
        reduced_relation_kernel_valid = false;
        reduced_relation_kernel_gram_d = flint::FmpzMat(0, 0);
        reduced_relation_kernel_gram_lambda = flint::FmpzMat(0, 0);
        reduced_relation_kernel_gram_valid = false;
    }

    void clear() noexcept {
        partial_relations.clear();
        direct_partial_residue_blocks.clear();
        partial_relations_max = 0;
        partial_relations_nonproductive_streak = 0;
        use_partial_relations = false;
        partial_relations_configured = false;

        random_ideal_search_state.base_indices.clear();
        random_ideal_search_state.exponents.clear();
        random_ideal_search_state.base_ideals.clear();
        random_ideal_search_state.inverse_base_numerators.clear();
        random_ideal_search_state.inverse_base_denominators.clear();
        random_ideal_search_state.random_ideal.clear();
        flint::fmpz_zero(
                flint::FmpzRef(random_ideal_search_state.lower_bound));
        flint::fmpz_zero(
                flint::FmpzRef(random_ideal_search_state.upper_bound));
        random_ideal_search_state.random_state = UWORD(0);
        random_ideal_search_state.defined = false;
        native_post_finite_reduced_lattices.clear();
        finish_state.reset();
        hnf_finish_workspace.clear();
        analytic_finish_precision = 0;
        flint::arb_zero(flint::ArbRef(analytic_finish_product));
        analytic_finish_product_precision = 0;
        analytic_finish_product_valid = false;

        factor_base_generation_records.clear();
        flint::fmpz_zero(flint::FmpzRef(grh_generation_bound));
        grh_generation_covered = false;
        grh_generation_basis = ClassGroupGenerationBasis::none;
        generation_different_route = GenerationDifferentRoute::none;
        class_unit_completeness = CertificationMode::unknown;
        class_unit_completeness_basis = ClassUnitCompletenessBasis::none;
        generation_structure_certified = false;
        relation_saturation_records.clear();
        relation_saturation_proof_records.clear();
        local_saturation_dlog_cache.clear();
        invalidate_reduced_relation_kernel();
    }
};

class ClassGroupRelationSearchAccess {
public:
    static slong relation_kernel_units_target(
            const ClassGroupContext& context) noexcept {
        return context.relation_kernel_units_target_;
    }

    static void set_relation_kernel_units_target(
            ClassGroupContext& context,
            slong target) noexcept {
        context.relation_kernel_units_target_ = target;
    }

    static bool build_relation_factor_base(
            ClassGroupContext& context,
            flint::FmpzConstRef bound) noexcept {
        return context.build_relation_factor_base_(bound);
    }

    static const RelationMatrix& relations(
            const ClassGroupContext& context) noexcept {
        return context.relations_;
    }

    static void record_skipped_dependent(
            ClassGroupContext& context) noexcept {
        ++context.skipped_dependent_relations_;
    }

    static bool append_search_relation(
            bool& skipped_dependent,
            ClassGroupContext& context,
            const Relation& relation) noexcept {
        ClassGroupContext::RelationAppendOutcome outcome =
                ClassGroupContext::RelationAppendOutcome::none;
        if (!context.append_relation_with_outcome_(
                    outcome, relation, ClassGroupRelationSource::Search,
                    ClassGroupContext::DependentRelationPolicy::keep)) {
            return false;
        }
        skipped_dependent = outcome ==
                ClassGroupContext::RelationAppendOutcome::skipped_dependent;
        return true;
    }

    static bool defer_native_goal_publication(
            const ClassGroupContext& context) noexcept;

    // The unconditional generation requirement recorded with the installed
    // factor base: proven_generation_requirement for the context's order,
    // evaluated once when the base was built
    // (record_factor_base_generation_).  False, leaving `out` unchanged,
    // when nothing is recorded or `order` is not the context's order; the
    // caller then computes the requirement itself.
    static bool recorded_generation_requirement(
            flint::FmpzRef out,
            const ClassGroupContext& context,
            const Order& order) noexcept {
        if (!context.has_factor_base() ||
            !same_order_parent(&context.parent_, &order) ||
            flint::fmpz_is_zero(flint::FmpzConstRef(
                    context.factor_base_generation_bound_))) {
            return false;
        }
        flint::fmpz_set(out, flint::FmpzConstRef(
                                     context.factor_base_generation_bound_));
        return true;
    }

    // Record, or clear, in the attached class/unit transaction report the
    // unwitnessed required prime that ended a candidate attempt; no-ops
    // without a transaction.  Clearing also clears the report's
    // factor-base size-limit flag, so each attempt starts clean.
    static void record_factor_base_honesty_failure(
            ClassGroupContext& context,
            FactorBaseHonestyFailure&& failure) noexcept;
    static void clear_factor_base_honesty_failure(
            ClassGroupContext& context) noexcept;

    // Whether the attached class/unit transaction keeps a passed honesty
    // scan (TransactionFactorBaseHonesty) of a factor base equal to `base`
    // at the same active and required bounds.  False without a transaction.
    static bool transaction_factor_base_honest(
            const ClassGroupContext& context,
            const FactorBase& base,
            flint::FmpzConstRef active_bound,
            flint::FmpzConstRef required_bound) noexcept;
    // Keep, in the attached class/unit transaction, a passed honesty scan of
    // `base` at these bounds, replacing any earlier one; no-op without a
    // transaction.  Nothing is kept if the base cannot be copied.
    static void record_transaction_factor_base_honest(
            ClassGroupContext& context,
            const FactorBase& base,
            flint::FmpzConstRef active_bound,
            flint::FmpzConstRef required_bound) noexcept;

    // The class-group candidate computation and extension entry points the
    // class/unit transaction calls for each attempt.
    static bool compute_relation_candidate(
            ClassGroupContext& context,
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const ClassGroupRelationOptions& options) noexcept {
        return context.compute_relation_candidate_(
                order, factor_base_bound, options);
    }
    static bool extend_relation_kernel_units(
            ClassGroupContext& context,
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const ClassGroupRelationOptions& options) noexcept {
        return context.extend_relation_kernel_units_(
                order, factor_base_bound, options);
    }

    static bool sync_row_module_checkpoint(
            ClassGroupContext& context) noexcept {
        return context.sync_row_module_checkpoint_();
    }

    static slong relation_kernel_row_count(
            const ClassGroupContext& context) noexcept {
        return relation_kernel_row_count_from_dimensions(
                context.generator_count(), context.relation_count(),
                context.relation_rank());
    }

    static slong relation_kernel_row_count_from_dimensions(
            slong generators,
            slong relations,
            slong rank) noexcept {
        return rank == generators && relations > generators
                ? relations - generators
                : 0;
    }

    static std::size_t direct_partial_residue_block_cache_size(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? 0
                : context.private_storage_
                          ->direct_partial_residue_blocks.size();
    }

    // `norm` must be the exact field norm of `generator` and
    // `integral_coordinates` its coordinates.
    static bool try_append_integral_generator_relation(
            ClassGroupContext& context,
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            ClassGroupRelationSource source) noexcept {
        return context.try_append_integral_generator_relation_(
                partial_throttle_exit, generator, integral_coordinates, norm,
                integral_coordinate_polynomial, source);
    }

    // `norm` must be the exact field norm of `generator`.
    static bool try_append_generator_relation_with_norm(
            ClassGroupContext& context,
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpqConstRef norm,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept {
        return context.try_append_generator_relation_with_norm_(
                partial_throttle_exit, generator, norm, source);
    }

    static ReducedIdealLatticeCache* native_post_finite_reduced_lattices(
            ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_
                           ->native_post_finite_reduced_lattices;
    }
};

}  // namespace silex::detail
