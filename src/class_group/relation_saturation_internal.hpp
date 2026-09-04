#pragma once

#include <vector>

#include <silex/class_group.hpp>
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
