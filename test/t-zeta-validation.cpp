#include "order_unit/class_unit_transaction_internal.hpp"
#include "order_unit/compute_internal.hpp"
#include "zeta/zeta_internal.hpp"

#include "test_support.hpp"

#include <silex/class_group.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>
#include <silex/zeta.hpp>

#include <cassert>

namespace {
namespace sflint = silex::flint;

struct FieldSetup {
    silex::NumberField field;
    silex::Order maximal_order;
};

FieldSetup setup_from_coefficients(const slong* coefficients,
                                   slong degree) noexcept {
    assert(coefficients != nullptr);
    assert(degree > 0);

    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, degree, 1);
    for (slong i = 0; i < degree; ++i) {
        if (coefficients[i] != 0) {
            sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
        }
    }

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

FieldSetup cubic_blocker() noexcept {
    const slong coefficients[] = {-5, -2, 0};
    return setup_from_coefficients(coefficients, 3);
}

FieldSetup quadratic_specialized_control() noexcept {
    const slong coefficients[] = {-5, 0};
    return setup_from_coefficients(coefficients, 2);
}

FieldSetup quartic_blocker() noexcept {
    const slong coefficients[] = {4, 1, 3, -8};
    return setup_from_coefficients(coefficients, 4);
}

void zeta_validation_target(sflint::Arb& out, slong precision) noexcept {
    sflint::Fmpq target;
    sflint::fmpq_set_si(target, 6931, 20000);
    sflint::arb_set_fmpq(out, target, precision);
}

void zeta_validation_tail_target(sflint::Arb& out, slong precision) noexcept {
    sflint::Arb validation_target;
    sflint::Arb numerical_reserve;
    zeta_validation_target(validation_target, precision);
    sflint::arb_one(numerical_reserve);
    sflint::arb_mul_2exp_si(numerical_reserve, numerical_reserve, -20);
    sflint::arb_sub(out, validation_target, numerical_reserve, precision);
}

bool radius_lt(const sflint::Arb& value,
               const sflint::Arb& target) noexcept {
    sflint::Arb radius;
    sflint::arb_get_rad_arb(radius, value);
    return sflint::arb_lt(radius, target);
}

bool call_validation(sflint::Arb& product,
                           sflint::Arb& error_bound,
                           ulong& cutoff,
                           slong& work_precision,
                           const silex::Order& order,
                           ulong max_cutoff,
                           slong precision) noexcept {
    return silex::detail::
            zeta_class_regulator_product_validation_with_diagnostics(
                    sflint::ArbRef(product), sflint::ArbRef(error_bound),
                    cutoff, work_precision, order, max_cutoff, precision,
                    nullptr);
}

bool check_failed_call_preserves_outputs(const silex::Order& order,
                                         ulong max_cutoff) noexcept {
    sflint::Arb product;
    sflint::Arb error_bound;
    sflint::arb_set_si(product, 17);
    sflint::arb_set_si(error_bound, 19);
    ulong cutoff = 23;
    slong work_precision = 29;

    return !call_validation(product, error_bound, cutoff,
                                  work_precision, order, max_cutoff, 128) &&
           ::arb_equal_si(product.raw(), 17) != 0 &&
           ::arb_equal_si(error_bound.raw(), 19) != 0 && cutoff == 23 &&
           work_precision == 29;
}

bool check_producer(const FieldSetup& setup,
                    ulong expected_cutoff,
                    ulong strict_max_cutoff) noexcept {
    const ulong preceding_cutoff = expected_cutoff - 9;
    if (!check_failed_call_preserves_outputs(setup.maximal_order,
                                             preceding_cutoff)) {
        return false;
    }

    sflint::Arb product;
    sflint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;
    if (!call_validation(product, error_bound, cutoff, work_precision,
                               setup.maximal_order, expected_cutoff, 128) ||
        cutoff != expected_cutoff || work_precision != 64 ||
        !sflint::arb_is_finite(product) ||
        !sflint::arb_is_positive(product) ||
        !sflint::arb_is_finite(error_bound) ||
        !sflint::arb_is_positive(error_bound)) {
        return false;
    }

    sflint::Arb tail_target;
    zeta_validation_tail_target(tail_target, 128);
    if (!sflint::arb_lt(error_bound, tail_target)) {
        return false;
    }

    sflint::Arb log_product;
    sflint::Arb validation_target;
    sflint::arb_log(log_product, product, 128);
    zeta_validation_target(validation_target, 128);
    if (!sflint::arb_is_finite(log_product) ||
        !radius_lt(log_product, validation_target)) {
        return false;
    }

    auto strict = silex::zeta_class_regulator_product_bf_audit(
            setup.maximal_order, strict_max_cutoff, 128);
    return strict.has_value() &&
           sflint::arb_overlaps(product, strict->value);
}

bool check_cache_separation(const FieldSetup& setup) noexcept {
    constexpr ulong max_cutoff = 20000;
    constexpr slong precision = 128;
    silex::detail::AnalyticClassRegulatorCache cache;
    cache.configure_validation(max_cutoff);
    if (!cache.ensure(setup.maximal_order, precision + 32, nullptr, nullptr,
                      true) ||
        !cache.has_validation(max_cutoff, precision + 32) ||
        cache.has_validation(max_cutoff, precision) ||
        !cache.ensure(setup.maximal_order, precision, nullptr, nullptr, true) ||
        !cache.has_validation(max_cutoff, precision) ||
        cache.has_bf_audit(max_cutoff, precision)) {
        return false;
    }

    sflint::Arb validation_snapshot;
    sflint::arb_set(sflint::ArbRef(validation_snapshot), cache.value());
    if (!cache.ensure_bf_audit(setup.maximal_order, max_cutoff, precision,
                               nullptr) ||
        !cache.has_bf_audit(max_cutoff, precision) ||
        ::arb_equal(validation_snapshot.raw(), cache.value().raw()) == 0) {
        return false;
    }

    auto strict = silex::zeta_class_regulator_product_bf_audit(
            setup.maximal_order, max_cutoff, precision);
    return strict.has_value() && cache.bf_cutoff() == strict->cutoff &&
           cache.bf_work_precision() == strict->work_precision &&
           ::arb_equal(cache.bf_error_bound().raw(),
                       strict->error_bound.raw()) != 0 &&
           sflint::arb_overlaps(cache.bf_value(), strict->value) &&
           sflint::arb_overlaps(cache.value(), validation_snapshot);
}

bool configure_native_options(silex::ClassGroupComputeOptions& options,
                              sflint::Fmpz& factor_base_bound,
                              const silex::Order& order) noexcept {
    options = silex::ClassGroupComputeOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;
    return true;
}

bool check_transaction_failure_leaves_outputs_unset(
        const FieldSetup& setup) noexcept {
    silex::ClassGroupComputeOptions options;
    sflint::Fmpz factor_base_bound;
    if (!configure_native_options(options, factor_base_bound,
                                  setup.maximal_order)) {
        return false;
    }
    options.max_candidates = 0;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport audit;
    const bool computed = silex::detail::compute_class_unit_transaction(
            units, class_group, setup.maximal_order,
            sflint::FmpzConstRef(factor_base_bound), options, 128, audit);
    return !computed &&
           audit.failure_stage != silex::detail::ClassUnitStage::none &&
           audit.failure_reason != nullptr &&
           !audit.final_result_published &&
           audit.class_group_certification ==
                   silex::CertificationMode::unknown &&
           audit.unit_group_certification ==
                   silex::CertificationMode::unknown &&
           !class_group.has_factor_base() && !class_group.has_presentation() &&
           !units.is_set();
}


}  // namespace

int main() {
    const FieldSetup quadratic = quadratic_specialized_control();
    const FieldSetup cubic = cubic_blocker();
    const FieldSetup quartic = quartic_blocker();

    return !check_failed_call_preserves_outputs(quadratic.maximal_order,
                                                20000) ||
                   !check_producer(cubic, 1008, 20000) ||
                   !check_producer(quartic, 1971, 30000) ||
                   !check_cache_separation(cubic) ||
                   !check_transaction_failure_leaves_outputs_unset(cubic) ||
                   !check_transaction_failure_leaves_outputs_unset(quartic)
            ? 1
            : 0;
}
