#pragma once

#include <silex/class_group.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>

#include "relation_completion_scheduler_internal.hpp"

#include <optional>
#include <vector>

namespace silex::detail {

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

struct FactorBaseGenerationRecord {
    flint::Fmpz p;
    ProofState status = ProofState::not_checked;
};

struct RelationSaturationRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
};

struct RelationSaturationProofRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
    slong rank = 0;
    slong target_rank = 0;
    slong local_primes = 0;
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
    std::vector<FractionalIdeal> base_fractional_ideals;
    std::vector<FractionalIdeal> inverse_base_ideals;
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
    std::vector<RelationSaturationRecord> relation_saturation_records;
    std::vector<RelationSaturationProofRecord>
            relation_saturation_proof_records;

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
        random_ideal_search_state.base_fractional_ideals.clear();
        random_ideal_search_state.inverse_base_ideals.clear();
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
        relation_saturation_records.clear();
        relation_saturation_proof_records.clear();
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

    static ReducedIdealLatticeCache* native_post_finite_reduced_lattices(
            ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_
                           ->native_post_finite_reduced_lattices;
    }
};

}  // namespace silex::detail
