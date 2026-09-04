#pragma once

#include "class_group_storage_internal.hpp"

namespace silex::detail {

bool rebuild_relation_factor_base_and_replay(
        ClassGroupContext& context,
        RelationAdmissionCache& cache,
        const Order& order,
        flint::FmpzConstRef next_bound,
        slong add_need) noexcept;

}  // namespace silex::detail
