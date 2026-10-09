#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <silex/class_group.hpp>
#include <silex/embedding.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order_unit.hpp>

namespace silex::detail {

inline constexpr slong kUnitCandidateSortPrecision = 32;
inline constexpr slong kTorsionStartPrecision = 16;
inline constexpr slong kRelAddStartPrecision = 32;
// Resource cap on the precision doubling of dependent-unit recovery; it is
// not a mathematical bound.  Reconstructing a denominator below 2^63 needs
// coordinate radii below 2^-127, so the cap leaves ample room for the
// log-embedding loss of large compact exponents.  Reaching it fails the
// search, and add_dependent_unit then leaves the group unchanged.
inline constexpr slong kRelAddMaxPrecision = WORD(1) << 16;
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

struct FactoredLogMemoValue {
    LogEmbeddingMode mode = LogEmbeddingMode::plain;
    slong precision = 0;
    std::uint64_t cache_version = 0;
    flint::ArbVec logs;

    explicit FactoredLogMemoValue(slong places) noexcept
        : logs(places) {
    }

    FactoredLogMemoValue(FactoredLogMemoValue&&) noexcept = default;
    FactoredLogMemoValue& operator=(FactoredLogMemoValue&&) noexcept =
            default;

    FactoredLogMemoValue(const FactoredLogMemoValue&) = delete;
    FactoredLogMemoValue& operator=(const FactoredLogMemoValue&) = delete;
};

// A memoized candidate, stored compactly: factor i is
// FactoredLogMemo::factors[factor_indices[i]] with exponent exponents[i].
struct FactoredLogMemoEntry {
    std::uint64_t fingerprint = 0;
    std::vector<std::size_t> factor_indices;
    std::vector<slong> exponents;
    std::vector<FactoredLogMemoValue> results;
};

// Logarithmic embeddings of relation-unit candidates, kept across the
// candidate batches and validation rounds of one extraction.  The same
// candidate (the same factors and exponents, in the same order) is embedded
// again by the reduction, the sort, and the torsion test, and again when a
// later round rebuilds it from the same relations.  Candidates repeat a few
// distinct factors many times, so each distinct factor is stored once.
//
// A hit returns exactly the value FactoredElement::logarithmic_embedding
// would compute again.  That call reads every factor log from the embedding
// context's cache and sums them in a fixed order; a result is stored only
// when the cache version did not change during the call, and returned only
// while the version is still the same.  Every factor log then reads back
// the same cached value, so the sum, its precision loop, and the side
// effects on the context (none: no factor is recomputed) are unchanged.
struct FactoredLogMemo {
    std::vector<Element> factors;
    std::vector<std::uint64_t> factor_hashes;
    std::vector<FactoredLogMemoEntry> entries;
    slong hits = 0;
    slong stores = 0;

    void clear() noexcept {
        factors.clear();
        factor_hashes.clear();
        entries.clear();
        hits = 0;
        stores = 0;
    }
};

// FactoredElement::logarithmic_embedding through `memo` (a direct call when
// `memo` is null).  Same arguments, result, and failure as the direct call.
bool memoized_factored_log_embedding(
        flint::ArbVecRef out,
        FactoredLogMemo* memo,
        const FactoredElement& value,
        EmbeddingContext& embeddings,
        LogEmbeddingMode mode,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

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
    FactoredLogMemo factored_log_memo;
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

bool conj_log_cutoff_inverse(const flint::ArbMat*& out,
                             RelationUnitExtractionState& extraction_state,
                             const OrderUnitGroup& group,
                             EmbeddingContext& embeddings,
                             slong& precision) noexcept;

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
// True when the reduce-mod-units precision ladder should stop on the
// coordinate matrix `rows`: some entry is a ball of radius at most 0.1 that
// contains a half-integer, and every entry that defeats the approximate
// rounding at `precision` is a ball of radius at most 0.1.
bool coordinate_rows_stop_at_half_integer_tie(const flint::ArbMat& rows,
                                              slong precision) noexcept;
bool reduce_relation_units_modulo(std::vector<FactoredElement>& values,
                            const OrderUnitGroup& group,
                            EmbeddingContext& embeddings,
                            slong start_precision,
                            RelationUnitExtractionState* extraction_state =
                                    nullptr) noexcept;

}  // namespace silex::detail
