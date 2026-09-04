#include "class_group/relation_completion_scheduler_internal.hpp"
#include "class_group/class_group_internal.hpp"
#include "test_support.hpp"

#include <silex/class_group.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order.hpp>
#include <silex/relation.hpp>

#include <initializer_list>
#include <vector>

namespace {
namespace sflint = silex::flint;
namespace scheduler = silex::detail::relation_search;

using Schedule = silex::detail::SubfactorBaseSchedule;
using CompletionState = silex::detail::RelationCompletionState;

struct FieldSetup {
    silex::NumberField field;
    silex::Order order;
};

FieldSetup cubic_x3_minus_2() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.order = silex::Order(setup.field);
    if (!setup.order.maximal_order(equation_order) ||
        !setup.order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

FieldSetup degree_one() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    setup.order = silex::test::equation_order(setup.field);
    return setup;
}

FieldSetup real_quadratic() noexcept {
    FieldSetup setup;
    setup.field = silex::test::quadratic_field(2);
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.order = silex::Order(setup.field);
    if (!setup.order.maximal_order(equation_order) ||
        !setup.order.is_maximal()) {
        return FieldSetup{};
    }
    return setup;
}

bool set_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

bool build_schedule(Schedule& schedule,
                    silex::FactorBase& base,
                    slong bound_value,
                    slong product_bound_value = 1) noexcept {
    sflint::Fmpz bound;
    sflint::Fmpz product_bound;
    return set_si(bound, bound_value) &&
           set_si(product_bound, product_bound_value) &&
           base.build_prime_ideal_norm_bounded(
                   sflint::FmpzConstRef(bound)) &&
           scheduler::build_initial_subfactor_base(
                   schedule, base,
                   sflint::FmpzConstRef(product_bound));
}

bool schedules_equal(const Schedule& left,
                     const Schedule& right) noexcept {
    return left.permutation == right.permutation &&
           left.subfactor_base == right.subfactor_base &&
           left.excluded_indices == right.excluded_indices &&
           left.dependency_growth_threshold ==
                   right.dependency_growth_threshold &&
           left.dependency_rotation_threshold ==
                   right.dependency_rotation_threshold;
}

bool build_relation_context(
        silex::ClassGroupContext& context,
        const silex::NumberField& field,
        std::initializer_list<slong> generators) noexcept {
    sflint::Fmpz bound;
    if (!set_si(bound, 2) ||
        !context.build_factor_base(sflint::FmpzConstRef(bound))) {
        return false;
    }

    const silex::FactorBase* base = context.factor_base();
    if (base == nullptr) {
        return false;
    }
    silex::Relation relation(*base);
    silex::Element generator(field);
    for (const slong value : generators) {
        if (!generator.set_si(value) ||
            !relation.set_generator(generator) ||
            !context.append_relation(relation)) {
            return false;
        }
    }
    return true;
}

bool decision_matches(
        const scheduler::RelationNeedDecision& decision,
        slong need,
        slong class_defect,
        slong kernel_defect,
        bool goal_reached,
        scheduler::RelationNeedReason reason) noexcept {
    return decision.need == need &&
           decision.class_relation_defect == class_defect &&
           decision.relation_kernel_unit_defect == kernel_defect &&
           decision.goal_reached == goal_reached &&
           decision.reason == reason;
}

int test_initial_subfactor_base_schedule() {
    FieldSetup setup = cubic_x3_minus_2();
    if (!setup.order.is_defined()) {
        return 1;
    }

    silex::FactorBase base(setup.order);
    Schedule schedule;
    if (!build_schedule(schedule, base, 30) || base.length() != 8 ||
        schedule.permutation !=
                std::vector<slong>({2, 4, 5, 0, 1, 6, 3, 7}) ||
        schedule.subfactor_base != std::vector<slong>({2, 4, 5}) ||
        schedule.excluded_indices !=
                std::vector<char>({1, 1, 0, 1, 0, 0, 0, 0}) ||
        schedule.dependency_growth_threshold != 48 ||
        schedule.dependency_rotation_threshold != 4) {
        return 1;
    }

    silex::FactorBase strict_cutoff_base(setup.order);
    Schedule strict_cutoff_schedule;
    if (!build_schedule(
                strict_cutoff_schedule, strict_cutoff_base, 30, 935) ||
        strict_cutoff_schedule.permutation !=
                std::vector<slong>({2, 4, 5, 6, 0, 1, 3, 7}) ||
        strict_cutoff_schedule.subfactor_base !=
                std::vector<slong>({2, 4, 5, 6}) ||
        strict_cutoff_schedule.excluded_indices !=
                std::vector<char>({1, 1, 0, 1, 0, 0, 0, 0}) ||
        strict_cutoff_schedule.dependency_growth_threshold != 64 ||
        strict_cutoff_schedule.dependency_rotation_threshold != 6) {
        return 1;
    }

    silex::FactorBase short_base(setup.order);
    Schedule short_schedule;
    return !build_schedule(short_schedule, short_base, 12) ||
                   short_base.length() != 4 ||
                   short_schedule.permutation !=
                           std::vector<slong>({2, 3, 0, 1}) ||
                   short_schedule.subfactor_base !=
                           std::vector<slong>({2, 3}) ||
                   short_schedule.excluded_indices !=
                           std::vector<char>({1, 1, 0, 0}) ||
                   short_schedule.dependency_growth_threshold != 32 ||
                   short_schedule.dependency_rotation_threshold != 3
            ? 1
            : 0;
}

int test_subfactor_base_changes() {
    using Change = scheduler::SubfactorBaseChange;
    using AdvanceDecision = scheduler::SubfactorBaseAdvanceDecision;

    FieldSetup setup = cubic_x3_minus_2();
    silex::FactorBase base(setup.order);
    Schedule initial;
    if (!setup.order.is_defined() || !build_schedule(initial, base, 30)) {
        return 1;
    }

    const std::vector<slong> useful = {4, 0, 2};
    Schedule changed = initial;
    if (!scheduler::change_subfactor_base(
                changed, base, useful, Change::rotate) ||
        changed.subfactor_base != std::vector<slong>({4, 2, 6}) ||
        changed.dependency_growth_threshold != 48 ||
        changed.dependency_rotation_threshold != 4 ||
        !scheduler::change_subfactor_base(
                changed, base, useful, Change::increase) ||
        changed.subfactor_base != std::vector<slong>({4, 2, 6, 7}) ||
        changed.dependency_growth_threshold != 64 ||
        changed.dependency_rotation_threshold != 6) {
        return 1;
    }

    silex::FactorBase short_base(setup.order);
    Schedule short_schedule;
    if (!build_schedule(short_schedule, short_base, 12)) {
        return 1;
    }
    const Schedule saved_short_schedule = short_schedule;
    if (scheduler::change_subfactor_base(
                short_schedule, short_base, {}, Change::increase) ||
        !schedules_equal(short_schedule, saved_short_schedule)) {
        return 1;
    }

    Schedule rotation_schedule = initial;
    CompletionState rotation_state;
    rotation_state.dependent_trials =
            rotation_schedule.dependency_rotation_threshold - 2;
    if (scheduler::advance_subfactor_base_schedule(
                rotation_schedule, rotation_state, base, useful,
                false, false) != AdvanceDecision::ready ||
        rotation_state.dependent_trials != 3 ||
        !schedules_equal(rotation_schedule, initial) ||
        scheduler::advance_subfactor_base_schedule(
                rotation_schedule, rotation_state, base, useful,
                false, false) != AdvanceDecision::ready ||
        rotation_state.dependent_trials != 4 ||
        rotation_schedule.subfactor_base !=
                std::vector<slong>({4, 2, 6})) {
        return 1;
    }

    Schedule growth_schedule = initial;
    CompletionState growth_state;
    growth_state.dependent_trials =
            growth_schedule.dependency_growth_threshold - 1;
    const std::vector<slong> stable_useful = {2, 4, 5};
    if (scheduler::advance_subfactor_base_schedule(
                growth_schedule, growth_state, base, stable_useful,
                false, false) != AdvanceDecision::ready ||
        growth_state.dependent_trials != 48 ||
        growth_state.subfactor_base_trials != 0 ||
        !schedules_equal(growth_schedule, initial) ||
        scheduler::advance_subfactor_base_schedule(
                growth_schedule, growth_state, base, stable_useful,
                false, false) != AdvanceDecision::ready ||
        growth_state.dependent_trials != 0 ||
        growth_state.subfactor_base_trials != 1 ||
        growth_schedule.subfactor_base !=
                std::vector<slong>({2, 4, 5, 6}) ||
        growth_schedule.dependency_growth_threshold != 64 ||
        growth_schedule.dependency_rotation_threshold != 6) {
        return 1;
    }

    Schedule guarded_schedule = initial;
    CompletionState guarded_state;
    guarded_state.dependent_trials =
            guarded_schedule.dependency_growth_threshold;
    guarded_state.subfactor_base_trials = 3;
    if (scheduler::advance_subfactor_base_schedule(
                guarded_schedule, guarded_state, base, useful,
                true, true) !=
                    AdvanceDecision::request_guarded_factor_base_rebuild ||
        guarded_state.dependent_trials != 49 ||
        guarded_state.subfactor_base_trials != 4 ||
        !schedules_equal(guarded_schedule, initial)) {
        return 1;
    }

    CompletionState limited_growth_state;
    limited_growth_state.dependent_trials =
            short_schedule.dependency_growth_threshold;
    if (scheduler::advance_subfactor_base_schedule(
                short_schedule, limited_growth_state, short_base, {},
                false, true) !=
                    AdvanceDecision::request_factor_base_rebuild_to_limit ||
        limited_growth_state.dependent_trials != 33 ||
        limited_growth_state.subfactor_base_trials != 1 ||
        !schedules_equal(short_schedule, saved_short_schedule)) {
        return 1;
    }

    short_schedule = saved_short_schedule;
    CompletionState failed_growth_state;
    failed_growth_state.dependent_trials =
            short_schedule.dependency_growth_threshold;
    if (scheduler::advance_subfactor_base_schedule(
                short_schedule, failed_growth_state, short_base, {},
                false, false) != AdvanceDecision::increase_failed ||
        failed_growth_state.dependent_trials != 33 ||
        failed_growth_state.subfactor_base_trials != 1 ||
        !schedules_equal(short_schedule, saved_short_schedule)) {
        return 1;
    }

    const std::vector<slong> insufficient_rotation = {2, 0, 1};
    short_schedule = saved_short_schedule;
    CompletionState limited_rotation_state;
    limited_rotation_state.dependent_trials =
            short_schedule.dependency_rotation_threshold - 1;
    if (scheduler::advance_subfactor_base_schedule(
                short_schedule, limited_rotation_state, short_base,
                insufficient_rotation, false, true) !=
                    AdvanceDecision::request_factor_base_rebuild_to_limit ||
        limited_rotation_state.dependent_trials != 3 ||
        limited_rotation_state.subfactor_base_trials != 0 ||
        !schedules_equal(short_schedule, saved_short_schedule)) {
        return 1;
    }

    short_schedule = saved_short_schedule;
    CompletionState failed_rotation_state;
    failed_rotation_state.dependent_trials =
            short_schedule.dependency_rotation_threshold - 1;
    return scheduler::advance_subfactor_base_schedule(
                   short_schedule, failed_rotation_state, short_base,
                   insufficient_rotation, false, false) !=
                           AdvanceDecision::rotation_failed ||
                   failed_rotation_state.dependent_trials != 3 ||
                   failed_rotation_state.subfactor_base_trials != 0 ||
                   !schedules_equal(short_schedule, saved_short_schedule)
            ? 1
            : 0;
}

int test_relation_need_decisions() {
    using Reason = scheduler::RelationNeedReason;

    FieldSetup setup = degree_one();
    silex::ClassGroupContext empty_context(setup.order);
    if (!build_relation_context(empty_context, setup.field, {})) {
        return 1;
    }
    const silex::FactorBase* empty_base = empty_context.factor_base();
    scheduler::RelationNeedDecision decision;
    if (empty_base == nullptr ||
        !scheduler::compute_relation_need(
                decision, empty_context, *empty_base, 3) ||
        !decision_matches(decision, 1, 1, 0, false,
                          Reason::class_relation_defect) ||
        empty_context.has_presentation()) {
        return 1;
    }

    silex::ClassGroupContext one_relation_context(setup.order);
    if (!build_relation_context(one_relation_context, setup.field, {4})) {
        return 1;
    }
    const silex::FactorBase* one_relation_base =
            one_relation_context.factor_base();
    if (one_relation_base == nullptr ||
        !scheduler::compute_relation_need(
                decision, one_relation_context, *one_relation_base, 3) ||
        !decision_matches(decision, 1, 0, 3, false,
                          Reason::relation_kernel_unit_defect) ||
        !one_relation_context.has_presentation()) {
        return 1;
    }

    silex::ClassGroupContext two_relation_context(setup.order);
    if (!build_relation_context(
                two_relation_context, setup.field, {4, 2})) {
        return 1;
    }
    const silex::FactorBase* two_relation_base =
            two_relation_context.factor_base();
    if (two_relation_base == nullptr ||
        !scheduler::compute_relation_need(
                decision, two_relation_context, *two_relation_base, 1) ||
        !decision_matches(decision, 0, 0, 0, true, Reason::none) ||
        !scheduler::compute_relation_need(
                decision, two_relation_context, *two_relation_base, 2) ||
        !decision_matches(decision, 1, 0, 1, false,
                          Reason::relation_kernel_unit_defect)) {
        return 1;
    }
    return 0;
}

int test_need_updates_and_targets() {
    CompletionState completion;
    completion.need = 3;
    completion.old_need = 3;
    completion.dependent_trials = 7;
    completion.subfactor_base_trials = 5;
    scheduler::update_relation_need(completion, 3);
    if (completion.need != 3 || completion.old_need != 3 ||
        completion.dependent_trials != 7 ||
        completion.subfactor_base_trials != 5) {
        return 1;
    }
    scheduler::update_relation_need(completion, 2);
    if (completion.need != 2 || completion.old_need != 2 ||
        completion.dependent_trials != 0 ||
        completion.subfactor_base_trials != 5) {
        return 1;
    }
    completion.dependent_trials = 9;
    scheduler::update_relation_need(completion, 2);
    if (completion.dependent_trials != 9) {
        return 1;
    }
    scheduler::update_relation_need(completion, 0);
    if (completion.need != 0 || completion.old_need != 0 ||
        completion.dependent_trials != 0 ||
        completion.subfactor_base_trials != 5) {
        return 1;
    }

    FieldSetup quadratic = real_quadratic();
    silex::detail::ClassGroupRelationOptions options;
    slong target = 99;
    if (!quadratic.order.is_defined() ||
        !scheduler::relation_kernel_unit_target(
                target, quadratic.order, options) ||
        target != 1) {
        return 1;
    }
    options.target_relation_kernel_units = 4;
    if (!scheduler::relation_kernel_unit_target(
                target, quadratic.order, options) ||
        target != 4) {
        return 1;
    }
    options.target_relation_kernel_units = -1;
    target = 99;
    if (scheduler::relation_kernel_unit_target(
                target, quadratic.order, options) ||
        target != 99) {
        return 1;
    }
    options.target_relation_kernel_units = 0;
    const silex::Order undefined_order;
    if (scheduler::relation_kernel_unit_target(
                target, undefined_order, options) ||
        target != 99) {
        return 1;
    }

    FieldSetup linear = degree_one();
    silex::ClassGroupContext context(linear.order);
    if (!build_relation_context(context, linear.field, {4, 2}) ||
        context.relation_count() != 2) {
        return 1;
    }

    completion = CompletionState{};
    completion.need = 0;
    options.max_relations = 5;
    if (scheduler::relation_collection_target_count(
                context, options, completion) != 2) {
        return 1;
    }
    completion.need = 2;
    if (scheduler::relation_collection_target_count(
                context, options, completion) != 4) {
        return 1;
    }
    completion.need = 10;
    if (scheduler::relation_collection_target_count(
                context, options, completion) != 5) {
        return 1;
    }
    options.max_relations = 1;
    if (scheduler::relation_collection_target_count(
                context, options, completion) != 2) {
        return 1;
    }
    completion.need = WORD_MAX;
    options.max_relations = WORD_MAX;
    return scheduler::relation_collection_target_count(
                   context, options, completion) != WORD_MAX
            ? 1
            : 0;
}

}  // namespace

int main() {
    return test_initial_subfactor_base_schedule() != 0 ||
                   test_subfactor_base_changes() != 0 ||
                   test_relation_need_decisions() != 0 ||
                   test_need_updates_and_targets() != 0
            ? 1
            : 0;
}
