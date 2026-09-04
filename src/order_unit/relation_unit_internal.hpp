#pragma once

#include <optional>
#include <vector>

#include <silex/class_group.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order_unit.hpp>

namespace silex::detail {

inline constexpr slong kUnitCandidateSortPrecision = 32;
inline constexpr slong kTorsionStartPrecision = 16;
inline constexpr slong kRelAddStartPrecision = 32;
inline constexpr slong kReduceModUnitsMaxPrecision = 10000;

enum class RelationTorsionStatus {
    torsion,
    non_torsion,
    inconclusive
};

struct RelationUnitCandidateBatchState {
    std::vector<slong> relations_used;
    slong extra_relation_count = 0;
    ulong random_state = UWORD(0x9e3779b97f4a7c15);
};

struct ReduceUnitsInverseCacheEntry {
    slong precision = 0;
    slong rank = 0;
    flint::ArbMat inverse;

    ReduceUnitsInverseCacheEntry(slong p, slong r) noexcept
        : precision(p),
          rank(r),
          inverse(r, r) {
    }

    ReduceUnitsInverseCacheEntry(
            ReduceUnitsInverseCacheEntry&&) noexcept = default;
    ReduceUnitsInverseCacheEntry& operator=(
            ReduceUnitsInverseCacheEntry&&) noexcept = default;

    ReduceUnitsInverseCacheEntry(
            const ReduceUnitsInverseCacheEntry&) = delete;
    ReduceUnitsInverseCacheEntry& operator=(
            const ReduceUnitsInverseCacheEntry&) = delete;
};

struct RelationUnitExtractionState {
    slong torsion_precision = kTorsionStartPrecision;
    slong rel_add_precision = kRelAddStartPrecision;
    slong cutoff_inverse_precision = 0;
    slong cutoff_inverse_rank = 0;
    std::optional<flint::ArbMat> cutoff_inverse;
    std::vector<ReduceUnitsInverseCacheEntry>
            reduce_mod_units_inverses;
    std::vector<FactoredElement> reduce_mod_units_log_generators;
    RelationUnitCandidateBatchState unit_candidate_batch;
    flint::Fmpz validation_index_bound;
    flint::Arb expected_regulator;
    bool has_expected_regulator = false;
    bool saturation_at_two_done = false;
    bool expand_reduce_mod_units_log_generators = false;
    bool reduce_mod_units_log_generator_expansion_attempted = false;
    bool reduce_mod_units_log_generators_expanded = false;

    const Element* expanded_reduce_mod_units_log_generator(
            slong index,
            slong rank) const noexcept {
        if (!reduce_mod_units_log_generators_expanded || index < 0 ||
            rank <= 0 ||
            reduce_mod_units_log_generators.size() !=
                    static_cast<std::size_t>(rank) ||
            index >= rank) {
            return nullptr;
        }

        const FactoredElement& generator =
                reduce_mod_units_log_generators[
                        static_cast<std::size_t>(index)];
        slong exponent = 0;
        const Element* factor = generator.factor(0);
        if (generator.length() != 1 || factor == nullptr ||
            !generator.exponent(exponent, 0) || exponent != 1) {
            return nullptr;
        }
        return factor;
    }

    void clear_dependent_unit_cache() noexcept {
        cutoff_inverse.reset();
        cutoff_inverse_precision = 0;
        cutoff_inverse_rank = 0;
        reduce_mod_units_inverses.clear();
        reduce_mod_units_log_generators.clear();
        reduce_mod_units_log_generator_expansion_attempted = false;
        reduce_mod_units_log_generators_expanded = false;
    }
};

bool unit_candidate_torsion_status(
        RelationTorsionStatus& status,
        RelationUnitExtractionState& extraction_state,
        const FactoredElement& candidate,
        EmbeddingContext& embeddings) noexcept;
bool add_dependent_unit(bool& changed,
                              OrderUnitGroup& group,
                              const FactoredElement& candidate,
                              EmbeddingContext& embeddings,
                              RelationUnitExtractionState& extraction_state,
                              slong precision) noexcept;
bool reduce_stored_relation_units(bool& changed,
                               OrderUnitGroup& group,
                               EmbeddingContext& embeddings,
                               slong precision) noexcept;
bool reduce_relation_units_modulo(std::vector<FactoredElement>& values,
                            const OrderUnitGroup& group,
                            EmbeddingContext& embeddings,
                            slong start_precision,
                            RelationUnitExtractionState* extraction_state =
                                    nullptr) noexcept;

}  // namespace silex::detail
