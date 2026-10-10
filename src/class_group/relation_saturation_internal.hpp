#pragma once

#include <cstddef>
#include <vector>

#include <silex/class_group.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex {

class EmbeddingContext;
class OrderUnitGroup;

namespace detail {

struct RelationUnitExtractionState;

struct SaturationCandidateSearchResult {
    slong auxiliary_prime_count = 0;
    slong degree_one_prime_count = 0;
    slong quotient_log_row_count = 0;
    slong skipped_quotient_log_row_count = 0;
    slong collapse_after_quotient_log_rows = 0;
    slong final_stable_count = 0;
    slong final_candidate_rows = 0;
    slong final_candidate_cols = 0;
    bool collapsed_to_zero = false;
};

struct SaturationCandidateProcessingResult {
    bool empty_candidate_space = false;
    bool wasted = false;
    bool saw_power = false;
    bool saw_non_power = false;
    bool saw_class_relation = false;
    bool saw_dependent_unit = false;
    slong inspected_columns = 0;
    slong tested_power_candidates = 0;
    slong structural_power_candidates = 0;
    slong compact_power_candidates = 0;
    slong appended_class_relations = 0;
    slong appended_dependent_units = 0;
    slong simplified_relation_count = 0;
    slong simplified_generator_count = 0;
    slong candidate_rows = 0;
    slong candidate_cols = 0;
    slong auxiliary_prime_count = 0;
    SaturationCandidateSearchResult candidate_search;
};

bool relation_basis_extra_indices(std::vector<slong>& basis_indices,
                                        std::vector<slong>& extra_indices,
                                        const ClassGroupContext& context)
        noexcept;
bool saturation_relation_data(std::vector<FactoredElement>& relations,
                                    flint::FmpzMatRef relation_rows,
                                    const ClassGroupContext& context) noexcept;
bool simplified_saturation_relation_data(
        std::vector<FactoredElement>& relations,
        flint::FmpzMat& relation_rows,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        slong n,
        slong precision) noexcept;
bool simplified_saturation_context(
        ClassGroupContext& simplified,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        const RelationUnitExtractionState& unit_state,
        slong n,
        slong precision) noexcept;
bool saturation_candidates(
        flint::FmpzMat& candidates,
        slong& auxiliary_prime_count,
        const ClassGroupContext& context,
        slong n,
        double stable,
        SaturationCandidateSearchResult* search_result = nullptr)
        noexcept;
bool saturation_process_candidates_once(
        SaturationCandidateProcessingResult& result,
        ClassGroupContext& context,
        flint::FmpzMatConstRef candidates,
        slong n,
        bool append_class_relations,
        OrderUnitGroup* units = nullptr,
        EmbeddingContext* embeddings = nullptr,
        RelationUnitExtractionState* unit_state = nullptr,
        slong precision = 0,
        slong max_class_relation_appends = WORD_MAX) noexcept;
bool saturation_process_candidates_once(
        SaturationCandidateProcessingResult& result,
        const ClassGroupContext& input_context,
        ClassGroupContext& append_context,
        flint::FmpzMatConstRef candidates,
        slong n,
        bool append_class_relations,
        OrderUnitGroup* units = nullptr,
        EmbeddingContext* embeddings = nullptr,
        RelationUnitExtractionState* unit_state = nullptr,
        slong precision = 0,
        slong max_class_relation_appends = WORD_MAX) noexcept;
// The degree-one quotient-log row of saturation_candidates_from_relations at
// (q, theta - root), from the exact fmpz evaluation and from the word-size
// path used for q < 2^64; exposed so tests can compare them.
bool degree_one_quotient_log_row_general(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& relations,
        flint::FmpzConstRef q,
        flint::FmpzConstRef root,
        slong n) noexcept;
bool degree_one_quotient_log_row_word_size(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& relations,
        flint::FmpzConstRef q,
        flint::FmpzConstRef root,
        slong n) noexcept;
// The factors of the exponent-reduced inputs of the degree-one residue
// characters, as rational polynomials, copied once so that each auxiliary
// prime reads them in place.
struct DegreeOneResidueFactor {
    flint::FmpqPoly polynomial;
    slong exponent = 0;
};

using DegreeOneResidueInput = std::vector<DegreeOneResidueFactor>;

// Class-relation saturation with units: the candidate filter from the
// stacked residue characters of reference
// `Saturate.jl:compute_candidates_for_saturate`, over the relation
// generators, free units and (when include_torsion) the torsion generator.
// The stack grows one degree-one prime at a time, only when a candidate that
// the current stack does not rule out is about to take its power test, and
// stops by the reference rule with `stable` starting at 1.5 and doubling
// after a failed power test.  See class_group.cpp.
struct RelationSaturationCandidateFilter {
    bool started = false;
    // False when the filter does not apply or failed; every candidate is
    // then power tested.
    bool available = false;
    // The reference stopping rule holds: the kernel dimension has stayed the
    // same for more than `stable * width` degree-one primes.
    bool stopped = false;
    slong ell = 0;
    slong width = 0;
    double stable = 1.5;
    slong stable_count = 1;
    slong stacked_prime_count = 0;
    const Order* order = nullptr;
    const DiagnosticsContext* diagnostics = nullptr;
    std::vector<FactoredElement> inputs;
    std::vector<DegreeOneResidueInput> prepared_inputs;
    bool prepared = false;
    flint::Fmpz discriminant;
    flint::Fmpz equation_index;
    bool equation_index_known = false;
    flint::Fmpz q;
    bool q_loaded = false;
    std::vector<flint::Fmpz> roots;
    std::size_t next_root = 0;
    // The stacked characters modulo ell as independent rows, each reduced
    // against the earlier ones (pivot entry 1); the kernel dimension is
    // width - rank.
    std::vector<std::vector<ulong>> rows;
    std::vector<slong> pivots;
};

bool relation_saturation_candidate_filter_start(
        RelationSaturationCandidateFilter& filter,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        bool include_torsion,
        flint::FmpzConstRef ell) noexcept;
bool relation_saturation_candidate_filter_grow(
        RelationSaturationCandidateFilter& filter) noexcept;
bool relation_saturation_candidate_filter_rules_out(
        bool& ruled_out,
        const RelationSaturationCandidateFilter& filter,
        flint::FmpzMatConstRef kernel,
        slong row) noexcept;
bool relation_saturation_candidate_filter_kernel_empty(
        const RelationSaturationCandidateFilter& filter) noexcept;
void relation_saturation_candidate_filter_failed_power_test(
        RelationSaturationCandidateFilter& filter) noexcept;
bool saturate_class_unit_context(
        bool& success,
        SaturationCandidateProcessingResult& last_result,
        ClassGroupContext& context,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        RelationUnitExtractionState& unit_state,
        slong n,
        double stable,
        slong precision) noexcept;

}  // namespace detail
}  // namespace silex
