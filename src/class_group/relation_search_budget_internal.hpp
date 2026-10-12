#pragma once

namespace silex::detail {

// The search budget, if any, that ended a slice of the LLL relation route
// (ClassGroupContext::run_lll_relation_route_) without reaching its goal:
// the candidate cap (ClassGroupRelationOptions::max_candidates, counted per
// route call) or the relation cap (max_relations, counted over all
// relations of the context).  Both caps are unlimited by default, so only a
// caller that sets them can see a value other than none.  A class/unit
// transaction that fails because a cap stopped the search reports
// relation_search_budget_exhausted with the cap's name.
enum class RelationSearchBudgetCap : unsigned char {
    none,
    max_candidates,
    max_relations,
};

// The option name of `cap` ("max_candidates" or "max_relations"), or
// nullptr for none.
inline const char* relation_search_budget_cap_name(
        RelationSearchBudgetCap cap) noexcept {
    switch (cap) {
    case RelationSearchBudgetCap::max_candidates:
        return "max_candidates";
    case RelationSearchBudgetCap::max_relations:
        return "max_relations";
    case RelationSearchBudgetCap::none:
        break;
    }
    return nullptr;
}

}  // namespace silex::detail
