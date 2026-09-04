#include "class_group/relation_factor_base_plan_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "test_support.hpp"
#include "zeta/zeta_internal.hpp"

#include <silex/class_group.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order.hpp>

#include <array>
#include <cassert>
#include <cstddef>

namespace {
namespace sflint = silex::flint;

struct FieldSetup {
    silex::NumberField field;
    silex::Order maximal_order;
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
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(equation_order));
    assert(setup.maximal_order.is_maximal());
    return setup;
}

FieldSetup degree_one() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    setup.maximal_order = silex::test::equation_order(setup.field);
    assert(setup.maximal_order.is_defined());
    assert(setup.maximal_order.degree() == 1);
    return setup;
}

FieldSetup imaginary_quadratic_minus_three() noexcept {
    FieldSetup setup;
    setup.field = silex::test::quadratic_field(-3);
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(equation_order));
    assert(setup.maximal_order.is_maximal());
    return setup;
}

FieldSetup quintic_discriminant_negative_401370255() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, -7);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -6);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 3);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(equation_order));
    assert(setup.maximal_order.is_maximal());
    return setup;
}

bool cache_shape_is_valid(
        const silex::detail::ZetaBfResidueDegreeCache& cache,
        slong field_degree) noexcept {
    if (cache.entries.empty() || cache.residue_degrees.empty()) {
        return false;
    }

    ulong previous_prime = 0;
    for (const auto& entry : cache.entries) {
        if (entry.p <= previous_prime || entry.length == 0 ||
            entry.offset > cache.residue_degrees.size() ||
            entry.length > cache.residue_degrees.size() - entry.offset) {
            return false;
        }
        previous_prime = entry.p;
        for (std::size_t i = 0; i < entry.length; ++i) {
            const slong degree = cache.residue_degrees[entry.offset + i];
            if (degree <= 0 || degree > field_degree) {
                return false;
            }
        }
    }
    return true;
}

bool caches_are_equal(
        const silex::detail::ZetaBfResidueDegreeCache& left,
        const silex::detail::ZetaBfResidueDegreeCache& right) noexcept {
    if (left.entries.size() != right.entries.size() ||
        left.residue_degrees != right.residue_degrees ||
        left.lookup_hint != right.lookup_hint) {
        return false;
    }
    for (std::size_t i = 0; i < left.entries.size(); ++i) {
        if (left.entries[i].p != right.entries[i].p ||
            left.entries[i].offset != right.entries[i].offset ||
            left.entries[i].length != right.entries[i].length) {
            return false;
        }
    }
    return true;
}

int test_exact_working_bounds_and_order_scoped_cache() {
    FieldSetup cubic = cubic_x3_minus_2();
    silex::detail::RelationFactorBasePlan cubic_plan;
    if (!silex::detail::relation_search::build_relation_factor_base_plan(
                cubic_plan, cubic.maximal_order, nullptr) ||
        !cubic_plan.valid ||
        !sflint::fmpz_equal_si(
                sflint::FmpzConstRef(cubic_plan.working_bound), 17)) {
        return 1;
    }

    auto* cubic_cache =
            silex::detail::relation_factor_base_plan_residue_degrees(
                    cubic_plan, cubic.maximal_order);
    if (cubic_cache != &cubic_plan.residue_degrees ||
        !cache_shape_is_valid(*cubic_cache, cubic.maximal_order.degree())) {
        return 1;
    }

    FieldSetup quintic = quintic_discriminant_negative_401370255();
    silex::detail::RelationFactorBasePlan quintic_plan;
    sflint::Fmpz public_bound;
    if (!silex::detail::relation_search::build_relation_factor_base_plan(
                quintic_plan, quintic.maximal_order, nullptr) ||
        !quintic_plan.valid ||
        !sflint::fmpz_equal_si(
                sflint::FmpzConstRef(quintic_plan.working_bound), 33) ||
        !silex::factor_base_class_group_bound(
                sflint::FmpzRef(public_bound), quintic.maximal_order) ||
        !sflint::fmpz_equal_si(
                sflint::FmpzConstRef(public_bound), 1539)) {
        return 1;
    }

    if (silex::detail::relation_factor_base_plan_residue_degrees(
                cubic_plan, quintic.maximal_order) != nullptr) {
        return 1;
    }

    silex::ClassGroupContext context(cubic.maximal_order);
    silex::detail::ClassUnitTransactionReport audit;
    audit.reset();
    silex::detail::ClassUnitTransactionContext run_context{
            .audit = audit,
            .relation_factor_base_plan = &cubic_plan,
            .imaginary_quadratic_exact_order_discriminant = {},
            .imaginary_quadratic_exact_order = {},
            .imaginary_quadratic_exact_order_valid = false};
    silex::detail::ClassUnitTransactionAccess::set_run_context(
            context, &run_context);
    const bool same_order =
            silex::detail::ClassUnitTransactionAccess::
                    relation_factor_base_plan_residue_degrees(
                            context, cubic.maximal_order) == cubic_cache;
    const bool crossed_order =
            silex::detail::ClassUnitTransactionAccess::
                    relation_factor_base_plan_residue_degrees(
                            context, quintic.maximal_order) != nullptr;
    silex::detail::ClassUnitTransactionAccess::set_run_context(context, nullptr);

    silex::ClassGroupContext foreign_context(quintic.maximal_order);
    silex::detail::ClassUnitTransactionAccess::set_run_context(
            foreign_context, &run_context);
    const bool crossed_context_parent =
            silex::detail::ClassUnitTransactionAccess::
                    relation_factor_base_plan_residue_degrees(
                            foreign_context, cubic.maximal_order) != nullptr;
    silex::detail::ClassUnitTransactionAccess::set_run_context(
            foreign_context, nullptr);
    return same_order && !crossed_order && !crossed_context_parent ? 0 : 1;
}

int test_degree_one_and_failed_build_preservation() {
    FieldSetup linear = degree_one();
    silex::detail::RelationFactorBasePlan linear_plan;
    const auto* linear_cache =
            silex::detail::relation_factor_base_plan_residue_degrees(
                    linear_plan, linear.maximal_order);
    if (linear_cache != nullptr ||
        !silex::detail::relation_search::build_relation_factor_base_plan(
                linear_plan, linear.maximal_order, nullptr) ||
        !linear_plan.valid ||
        !sflint::fmpz_is_one(
                sflint::FmpzConstRef(linear_plan.working_bound)) ||
        !linear_plan.residue_degrees.entries.empty() ||
        !linear_plan.residue_degrees.residue_degrees.empty() ||
        silex::detail::relation_factor_base_plan_residue_degrees(
                linear_plan, linear.maximal_order) !=
                &linear_plan.residue_degrees) {
        return 1;
    }

    FieldSetup cubic = cubic_x3_minus_2();
    silex::detail::RelationFactorBasePlan plan;
    if (!silex::detail::relation_search::build_relation_factor_base_plan(
                plan, cubic.maximal_order, nullptr)) {
        return 1;
    }
    sflint::Fmpz saved_bound;
    sflint::fmpz_set(sflint::FmpzRef(saved_bound),
                     sflint::FmpzConstRef(plan.working_bound));
    const silex::Order saved_order = plan.order;
    const auto saved_cache = plan.residue_degrees;
    const bool saved_valid = plan.valid;

    const silex::Order undefined_order;
    return silex::detail::relation_search::build_relation_factor_base_plan(
                   plan, undefined_order, nullptr) ||
                   plan.valid != saved_valid ||
                   !plan.order.has_same_data(saved_order) ||
                   !sflint::fmpz_equal(
                           sflint::FmpzConstRef(plan.working_bound),
                           sflint::FmpzConstRef(saved_bound)) ||
                   !caches_are_equal(plan.residue_degrees, saved_cache)
            ? 1
            : 0;
}

int test_producer_cache_matches_fresh_bf_evaluation() {
    FieldSetup setup = cubic_x3_minus_2();
    silex::detail::RelationFactorBasePlan plan;
    if (!silex::detail::relation_search::build_relation_factor_base_plan(
                plan, setup.maximal_order, nullptr)) {
        return 1;
    }

    sflint::Arb cached_product;
    sflint::Arb cached_error;
    ulong cached_cutoff = 0;
    slong cached_precision = 0;
    auto consumed_cache = plan.residue_degrees;
    if (!silex::detail::zeta_class_regulator_product_bf_audit_with_diagnostics(
                sflint::ArbRef(cached_product),
                sflint::ArbRef(cached_error), cached_cutoff,
                cached_precision, setup.maximal_order, 20000, 128, nullptr,
                nullptr, &consumed_cache)) {
        return 1;
    }

    sflint::Arb fresh_product;
    sflint::Arb fresh_error;
    ulong fresh_cutoff = 0;
    slong fresh_precision = 0;
    silex::detail::ZetaBfResidueDegreeCache fresh_cache;
    return !silex::detail::zeta_class_regulator_product_bf_audit_with_diagnostics(
                   sflint::ArbRef(fresh_product),
                   sflint::ArbRef(fresh_error), fresh_cutoff,
                   fresh_precision, setup.maximal_order, 20000, 128, nullptr,
                   nullptr, &fresh_cache) ||
                   cached_cutoff != fresh_cutoff ||
                   cached_precision != fresh_precision ||
                   !sflint::arb_overlaps(cached_product, fresh_product) ||
                   !sflint::arb_overlaps(cached_error, fresh_error)
            ? 1
            : 0;
}

int test_restart_bound_progression() {
    FieldSetup setup = cubic_x3_minus_2();
    sflint::Fmpz limit;
    if (!silex::detail::relation_search::relation_factor_base_restart_limit(
                limit, setup.maximal_order) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(limit), 87)) {
        return 1;
    }

    sflint::Fmpz current;
    sflint::Fmpz next;
    sflint::fmpz_set_si(sflint::FmpzRef(current), 17);
    constexpr std::array<slong, 7> expected = {
            21, 25, 29, 33, 37, 41, 45};
    for (const slong value : expected) {
        if (!silex::detail::relation_search::next_relation_factor_base_bound(
                    next, sflint::FmpzConstRef(current),
                    sflint::FmpzConstRef(limit)) ||
            !sflint::fmpz_equal_si(sflint::FmpzConstRef(next), value)) {
            return 1;
        }
        sflint::fmpz_set(sflint::FmpzRef(current),
                         sflint::FmpzConstRef(next));
    }
    if (silex::detail::relation_search::next_relation_factor_base_bound(
                next, sflint::FmpzConstRef(current),
                sflint::FmpzConstRef(limit))) {
        return 1;
    }

    sflint::fmpz_set_si(sflint::FmpzRef(current), 85);
    if (!silex::detail::relation_search::
                 next_relation_factor_base_bound_to_limit(
                         next, sflint::FmpzConstRef(current),
                         sflint::FmpzConstRef(limit)) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(next), 87)) {
        return 1;
    }

    sflint::fmpz_set_si(sflint::FmpzRef(limit), 100);
    sflint::fmpz_set_si(sflint::FmpzRef(current), 1);
    if (!silex::detail::relation_search::
                 next_relation_factor_base_bound_to_limit(
                         next, sflint::FmpzConstRef(current),
                         sflint::FmpzConstRef(limit)) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(next), 2)) {
        return 1;
    }

    sflint::fmpz_set_si(sflint::FmpzRef(current), 8);
    if (!silex::detail::relation_search::
                 next_relation_factor_base_bound_to_limit(
                         next, sflint::FmpzConstRef(current),
                         sflint::FmpzConstRef(limit)) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(next), 13)) {
        return 1;
    }

    sflint::fmpz_set_si(sflint::FmpzRef(current), 50);
    if (!silex::detail::relation_search::next_relation_factor_base_bound(
                next, sflint::FmpzConstRef(current),
                sflint::FmpzConstRef(limit)) ||
        !sflint::fmpz_equal_si(sflint::FmpzConstRef(next), 55)) {
        return 1;
    }

    sflint::fmpz_set_si(sflint::FmpzRef(current), 51);
    sflint::fmpz_set_si(sflint::FmpzRef(next), 999);
    return silex::detail::relation_search::next_relation_factor_base_bound(
                   next, sflint::FmpzConstRef(current),
                   sflint::FmpzConstRef(limit)) ||
                   !sflint::fmpz_equal_si(sflint::FmpzConstRef(next), 999)
            ? 1
            : 0;
}

int test_discriminant_minus_three_exact_plan() {
    FieldSetup setup = imaginary_quadratic_minus_three();
    silex::detail::RelationFactorBasePlan plan;
    sflint::Fmpz public_generation_bound;
    sflint::Fmpz generic_restart_limit;
    if (!silex::detail::relation_search::
                 build_maximal_imaginary_quadratic_factor_base_plan(
                         plan, setup.maximal_order, nullptr) ||
        !plan.valid || !plan.order.has_same_data(setup.maximal_order) ||
        !silex::factor_base_class_group_bound(
                sflint::FmpzRef(public_generation_bound),
                setup.maximal_order) ||
        sflint::fmpz_cmp(
                sflint::FmpzConstRef(plan.working_bound),
                sflint::FmpzConstRef(public_generation_bound)) < 0 ||
        !silex::detail::relation_search::relation_factor_base_restart_limit(
                generic_restart_limit, setup.maximal_order) ||
        sflint::fmpz_cmp(
                sflint::FmpzConstRef(plan.working_bound),
                sflint::FmpzConstRef(generic_restart_limit)) <= 0) {
        return 1;
    }
    return 0;
}

}  // namespace

int main() {
    return test_exact_working_bounds_and_order_scoped_cache() != 0 ||
                   test_degree_one_and_failed_build_preservation() != 0 ||
                   test_producer_cache_matches_fresh_bf_evaluation() != 0 ||
                   test_restart_bound_progression() != 0 ||
                   test_discriminant_minus_three_exact_plan() != 0
            ? 1
            : 0;
}
