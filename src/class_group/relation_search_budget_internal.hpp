#pragma once

#include <flint/flint.h>

namespace silex::detail {

// The search budget, if any, that ended a relation search without reaching
// its goal: the candidate cap (ClassGroupRelationOptions::max_candidates,
// counted per route call) or the relation cap (max_relations, counted over
// all relations of the context).  Both caps are unlimited by default, so
// only a caller that sets them can see a value other than none.  A
// class/unit transaction that fails because a cap stopped the search
// reports relation_search_budget_exhausted with the cap's name.
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

// The cap reached by a relation search that has tried `candidates_tried`
// candidates and holds `accepted_relations` relations: the candidate cap
// when it is reached, else the relation cap when it is reached, else none.
// Callers ask only when the search ended short of its goal, so a reached
// cap is what stopped it.
inline RelationSearchBudgetCap relation_search_budget_cap_reached(
        slong candidates_tried,
        slong max_candidates,
        slong accepted_relations,
        slong max_relations) noexcept {
    if (candidates_tried >= max_candidates) {
        return RelationSearchBudgetCap::max_candidates;
    }
    if (accepted_relations >= max_relations) {
        return RelationSearchBudgetCap::max_relations;
    }
    return RelationSearchBudgetCap::none;
}

}  // namespace silex::detail
