#pragma once

#include "class_group_storage_internal.hpp"

#include <vector>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;

}  // namespace silex::detail

namespace silex::detail::relation_search {

struct NormPrefilter;

bool publish_and_check_compute_goal(
        ClassGroupContext& context,
        slong target_relation_kernel_units,
        CertificationMode requested_certification) noexcept;

bool collect_integral_ideal_lattice_relations_den(
        ClassGroupContext& context,
        const Ideal& ideal,
        flint::FmpzConstRef den,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        ulong restart_seed,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::RelationAdmissionCache* admission_cache = nullptr,
        bool random_relation = false,
        slong* factor_attempts = nullptr,
        bool* factor_attempt_limit_reached = nullptr,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept;

bool collect_integral_ideal_lattice_relations(
        ClassGroupContext& context,
        const Ideal& ideal,
        slong ideal_search_radius,
        NormPrefilter* norm_prefilter,
        ulong restart_seed,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::ReducedIdealLatticeCache* reduced_lattice_cache = nullptr)
        noexcept;

bool collect_order_small_lll_relations(
        ClassGroupContext& context,
        const Order& order,
        const Ideal& order_ideal,
        NormPrefilter* norm_prefilter,
        slong target_relation_kernel_units,
        const ClassGroupRelationOptions& options,
        slong& candidates_tried,
        slong& accepted_relations,
        bool& goal_reached,
        bool& partial_throttle_exit,
        ClassGroupRelationSource source = ClassGroupRelationSource::Search,
        slong max_good_relations = WORD(-1)) noexcept;

bool collect_private_direct_t2_relations_den(
        ClassGroupContext& context,
        const Ideal& ideal,
        flint::FmpzConstRef den,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong max_relations_per_ideal,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        detail::RelationAdmissionCache& admission_cache,
        bool random_relation,
        detail::OrderMinkowskiEmbeddingCache* t2_embedding_cache) noexcept;

bool build_nonprincipal_indices(std::vector<slong>& nonprincipal,
                                ClassGroupContext& context,
                                const FactorBase& base) noexcept;

bool build_hnf_covered_flags(std::vector<char>& covered,
                             ClassGroupContext& context,
                             const FactorBase& base) noexcept;

bool build_uncovered_indices_from_flags(
        std::vector<slong>& uncovered,
        const std::vector<char>& hnf_covered,
        const std::vector<slong>& nonprincipal) noexcept;

slong max_slong_value(slong left, slong right) noexcept;
slong min_slong_value(slong left, slong right) noexcept;

ulong relation_search_phase_seed(const ClassGroupContext& context,
                                 slong factor_base_length,
                                 ulong phase,
                                 slong phase_restarts) noexcept;

ulong next_relation_random_exponent(ulong& state) noexcept;

}  // namespace silex::detail::relation_search
