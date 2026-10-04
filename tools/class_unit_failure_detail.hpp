#pragma once

#include "order_unit/class_unit_transaction_internal.hpp"

#include <cstring>
#include <ostream>

namespace silex_tools {

// Writes the value of the class/unit instance JSON field `failure_detail`,
// the one place for structured detail about a failure: an object whose keys
// depend on `failure_reason`, or null when that reason carries no detail.
// Only factor_base_honesty_unwitnessed carries detail: the unwitnessed
// required prime (`p`, null if it does not fit a signed machine word, and
// `residue_degree`), its witness search route (`search`, "lattice" or
// "t2"), the stage of its last search against the last stage allowed, and
// that stage's caps under the class-group detail log's names: `radius`,
// `twists` and `random_tries` on the lattice route, `random_tries`,
// `factor_attempts` and `element_steps` on the T2 route.  Nested lines are
// indented for a top-level field.
inline void write_class_unit_failure_detail_json(
        std::ostream& out,
        const silex::detail::ClassUnitTransactionReport& report) {
    const silex::detail::FactorBaseHonestyFailure& failure =
            report.factor_base_honesty_failure;
    if (report.failure_reason == nullptr ||
        std::strcmp(report.failure_reason,
                    "factor_base_honesty_unwitnessed") != 0 ||
        !failure.recorded) {
        out << "null";
        return;
    }
    const auto field = [&out](const char* name, slong value) {
        out << ",\n    \"" << name << "\": " << static_cast<long long>(value);
    };
    out << "{\n    \"p\": ";
    if (failure.p < 0) {
        out << "null";
    } else {
        out << static_cast<long long>(failure.p);
    }
    field("residue_degree", failure.residue_degree);
    out << ",\n    \"search\": \""
        << (failure.direct_witness_search ? "t2" : "lattice") << "\"";
    field("stage", failure.stage);
    field("max_stage", failure.max_stage);
    if (failure.direct_witness_search) {
        field("random_tries", failure.random_tries);
        field("factor_attempts", failure.factor_attempts);
        field("element_steps", failure.element_steps);
    } else {
        field("radius", failure.radius);
        field("twists", failure.max_twists);
        field("random_tries", failure.random_tries);
    }
    out << "\n  }";
}

}  // namespace silex_tools
