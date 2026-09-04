#pragma once

#include "class_group_storage_internal.hpp"

#include <vector>

namespace silex::detail {

bool initialize_relation_admission_cache(
        ClassGroupContext& context,
        RelationAdmissionCache& cache,
        const Order& order,
        slong add_need) noexcept;

bool try_admit_relation(
        ClassGroupContext& context,
        bool& retained,
        RelationAdmissionCache& cache,
        const Relation& relation,
        bool in_random_relation,
        bool factor_base_verified = false) noexcept;

bool try_admit_deferred_integral_relation(
        ClassGroupContext& context,
        bool& retained,
        RelationAdmissionCache& cache,
        Relation& relation,
        flint::FmpzMatConstRef integral_coordinates,
        bool in_random_relation,
        bool factor_base_verified = false) noexcept;

}  // namespace silex::detail

namespace silex::detail::relation_search {

bool begin_selected_pivot_recollection(
        detail::RelationAdmissionCache& cache,
        const std::vector<slong>& indices,
        slong n) noexcept;

void end_selected_pivot_recollection(
        detail::RelationAdmissionCache& cache,
        const std::vector<slong>& indices,
        slong n) noexcept;

}  // namespace silex::detail::relation_search
