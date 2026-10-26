#pragma once

#include "class_group_storage_internal.hpp"

#include <vector>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;

}  // namespace silex::detail

namespace silex::detail::relation_search {

struct NormPrefilter;
struct IdealLatticeRoundedEmbeddingCache;

// Torsion orbits of small-LLL candidates.
//
// In a totally complex field with w > 2 roots of unity, let Z be the
// order-coordinate matrix of multiplication by a primitive w-th root of
// unity zeta (row j holds the coordinates of zeta * b_j).  Candidate
// coordinates u and u * Z^k describe alpha and zeta^k * alpha.  These
// generate the same ideal, so they have the same |N(alpha)| (the norm
// prefilter decides on that alone) and, when alpha is smooth, the same
// factor-base exponent row.  The candidate loop uses this only to skip a
// candidate whose full test is known to end in the keep_nonduplicate
// duplicate-row skip of append_relation_with_outcome_; see
// try_small_lll_orbit_skip (lll_relation_search.cpp).  Fields with a real
// place have w = 2 and are left inactive.
struct SmallLllTorsionOrbit {
    bool active = false;
    slong degree = 0;
    slong torsion_order = 0;
    // Z^k = numerators[k - 1] / denominators[k - 1], k = 1..w-1, exactly.
    std::vector<flint::FmpzMat> numerators;
    std::vector<flint::Fmpz> denominators;
    // key_rows[k * degree + j] = (Z^k c)_j mod a fixed prime for a fixed
    // vector c, k = 0..w-1, so that the key c . (u Z^k) of an image is
    // u . (Z^k c).  Keys only index exact checks; they decide nothing.
    std::vector<ulong> key_rows;
    // Candidates the loop skipped through this data.
    slong skipped_candidates = 0;
};

// Leaves `out` inactive (the loop then tests every candidate in full) when
// the field is not totally complex, has only +-1, or any step fails; a
// failure here never fails the relation route.
void build_small_lll_torsion_orbit(SmallLllTorsionOrbit& out,
                                   const Order& order,
                                   const NormPrefilter* norm_prefilter)
        noexcept;

// Relation-kernel dimension counted toward a kernel-unit target: the
// kernel dimension (relation_kernel_unit_count, or the deferred
// relation-kernel row count) less the Galois-image relations stored without
// raising the modular rank.  Kernel images are left out of the target by
// design (see relation_candidate.cpp, the Galois-image step); off the
// cyclotomic image route both equal the plain counts.
slong target_relation_kernel_unit_count(
        const ClassGroupContext& context) noexcept;
slong target_relation_kernel_row_count(
        const ClassGroupContext& context) noexcept;

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
        slong max_good_relations = WORD(-1),
        IdealLatticeRoundedEmbeddingCache* lll_cache = nullptr,
        SmallLllTorsionOrbit* torsion_orbit = nullptr) noexcept;

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
