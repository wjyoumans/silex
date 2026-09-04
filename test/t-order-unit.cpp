#include <silex/order_unit.hpp>

#include <silex/class_group.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_field.hpp>
#include <silex/signature.hpp>
#include <silex/zeta.hpp>

#include "order_unit/order_unit_internal.hpp"
#include "order_unit/relation_unit_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <limits>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

int test_unit_extraction_cache_clear() {
    silex::NumberField field = silex::test::quadratic_field(2);
    silex::detail::RelationUnitExtractionState state;
    state.cutoff_inverse.emplace(1, 1);
    state.cutoff_inverse_precision = 64;
    state.cutoff_inverse_rank = 1;
    state.reduce_mod_units_inverses.emplace_back(64, 1);
    state.reduce_mod_units_log_generators.emplace_back(field);
    assert(state.reduce_mod_units_log_generators.back().one());
    state.expand_reduce_mod_units_log_generators = true;
    state.reduce_mod_units_log_generator_expansion_attempted = true;
    state.reduce_mod_units_log_generators_expanded = true;

    state.clear_dependent_unit_cache();
    assert(!state.cutoff_inverse.has_value());
    assert(state.cutoff_inverse_precision == 0);
    assert(state.cutoff_inverse_rank == 0);
    assert(state.reduce_mod_units_inverses.empty());
    assert(state.reduce_mod_units_log_generators.empty());
    assert(state.expand_reduce_mod_units_log_generators);
    assert(!state.reduce_mod_units_log_generator_expansion_attempted);
    assert(!state.reduce_mod_units_log_generators_expanded);
    return 0;
}

void set_zmat_entry_si(sflint::FmpzMat& matrix,
                       slong row,
                       slong column,
                       slong value) noexcept {
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(matrix), row, column),
            value);
}

void set_zmat_entry_str(sflint::FmpzMat& matrix,
                        slong row,
                        slong column,
                        const char* value) noexcept {
    assert(sflint::fmpz_set_str(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(matrix), row, column),
            value));
}

void set_arb_mat_entry_si_2exp(sflint::ArbMat& matrix,
                               slong row,
                               slong column,
                               slong mantissa,
                               slong exponent) noexcept {
    arb_struct* entry = arb_mat_entry(matrix.raw(), row, column);
    ::arf_set_si(arb_midref(entry), mantissa);
    ::arf_mul_2exp_si(arb_midref(entry), arb_midref(entry), exponent);
    ::mag_zero(arb_radref(entry));
}

void add_arb_mat_entry_2exp(sflint::ArbMat& matrix,
                            slong row,
                            slong column,
                            slong exponent) noexcept {
    sflint::Arf term;
    ::arf_set_ui(term.raw(), 1);
    ::arf_mul_2exp_si(term.raw(), term.raw(), exponent);
    arb_struct* entry = arb_mat_entry(matrix.raw(), row, column);
    ::arf_add(arb_midref(entry), arb_midref(entry), term.raw(),
              ARF_PREC_EXACT, ARF_RND_NEAR);
    ::mag_zero(arb_radref(entry));
}

void set_arb_mat_entry_interval_si_2exp(
        sflint::ArbMat& matrix,
        slong row,
        slong column,
        slong lower_mantissa,
        slong lower_exponent,
        slong upper_mantissa,
        slong upper_exponent) noexcept {
    sflint::Arf lower;
    sflint::Arf upper;
    ::arf_set_si(lower.raw(), lower_mantissa);
    ::arf_mul_2exp_si(lower.raw(), lower.raw(), lower_exponent);
    ::arf_set_si(upper.raw(), upper_mantissa);
    ::arf_mul_2exp_si(upper.raw(), upper.raw(), upper_exponent);
    ::arb_set_interval_arf(arb_mat_entry(matrix.raw(), row, column),
                           lower.raw(), upper.raw(), 256);
}

void set_arb_mat_entry_huge_2exp(sflint::ArbMat& matrix,
                                 slong row,
                                 slong column,
                                 const char* exponent) noexcept {
    sflint::Fmpz power;
    assert(sflint::fmpz_set_str(sflint::FmpzRef(power), exponent));
    arb_struct* entry = arb_mat_entry(matrix.raw(), row, column);
    ::arf_one(arb_midref(entry));
    ::arf_mul_2exp_fmpz(arb_midref(entry), arb_midref(entry), power.raw());
    ::mag_zero(arb_radref(entry));
}

void assert_regulator_pivot_rows(
        const sflint::ArbMat& matrix,
        const std::vector<slong>& expected) noexcept {
    std::vector<slong> actual;
    assert(silex::detail::regulator_pivot_rows_for_testing(
            actual, sflint::ArbMatConstRef(matrix), 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(actual == expected);
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(matrix), row, col),
            value);
}

bool check_unit_proof_record(const silex::OrderUnitGroup& group,
                             slong index,
                             slong expected_ell,
                             silex::ProofState expected_status,
                             slong expected_aux_bound,
                             slong expected_local_primes,
                             bool expected_changed) noexcept {
    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    silex::ProofState status = silex::ProofState::not_checked;
    slong local_primes = -1;
    bool changed = !expected_changed;
    return group.unit_proof_record(sflint::FmpzRef(ell), status,
                                   sflint::FmpzRef(aux_bound), local_primes,
                                   changed, index) &&
           sflint::fmpz_equal_si(ell, expected_ell) &&
           status == expected_status &&
           sflint::fmpz_equal_si(aux_bound, expected_aux_bound) &&
           (expected_local_primes >= 0
                    ? local_primes == expected_local_primes
                    : local_primes > 0) &&
           changed == expected_changed;
}

bool has_unit_proof_record(const silex::OrderUnitGroup& group,
                           slong expected_ell,
                           silex::ProofState expected_status,
                           bool expected_changed) noexcept {
    for (slong i = 0; i < group.unit_proof_record_count(); ++i) {
        sflint::Fmpz ell;
        sflint::Fmpz aux_bound;
        silex::ProofState status = silex::ProofState::not_checked;
        slong local_primes = 0;
        bool changed = !expected_changed;
        if (group.unit_proof_record(sflint::FmpzRef(ell), status,
                                    sflint::FmpzRef(aux_bound), local_primes,
                                    changed, i) &&
            sflint::fmpz_equal_si(ell, expected_ell) &&
            status == expected_status && local_primes > 0 &&
            changed == expected_changed) {
            return true;
        }
    }
    return false;
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField cubic_field(slong c1, slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quartic_field(slong c3,
                                 slong c2,
                                 slong c1,
                                 slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, c3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, c2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quintic_field(slong c3,
                                 slong c2,
                                 slong c1,
                                 slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, c3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, c2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::OrderUnitGroup real_quadratic_unit_group() noexcept {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);

    silex::OrderUnitGroup group;
    assert(group.compute(order));
    return group;
}

bool order_value_has_power_one(const silex::OrderElement& generator,
                               slong exponent) noexcept {
    const silex::Order* order = generator.parent();
    if (order == nullptr || order->parent() == nullptr) {
        return false;
    }

    silex::Element value(*order->parent());
    silex::Element power(*order->parent());
    sflint::Fmpz exp;
    sflint::fmpz_set_si(sflint::FmpzRef(exp), exponent);
    return generator.get_element(value) &&
           power.pow_fmpz(value, sflint::FmpzConstRef(exp)) &&
           power.equal_si(1);
}

bool check_rank_zero_group(const silex::OrderUnitGroup& group,
                           const silex::Order& order,
                           slong expected_order) noexcept {
    if (!group.is_defined() || !group.is_set() ||
        !silex::same_order_parent(group.parent(), &order) ||
        group.free_rank() != 0 ||
        group.certification_status() != silex::CertificationMode::proven) {
        return false;
    }

    sflint::Fmpz torsion_order;
    silex::OrderElement generator(order);
    sflint::Arb regulator;
    return group.torsion_order(sflint::FmpzRef(torsion_order)) &&
           sflint::fmpz_equal_si(torsion_order, expected_order) &&
           group.torsion_generator(generator) &&
           order_value_has_power_one(generator, expected_order) &&
           group.regulator(sflint::ArbRef(regulator)) &&
           sflint::arb_is_one(regulator);
}

bool set_real_quadratic_unit(silex::Element& out) noexcept {
    const silex::NumberField* field = out.parent();
    if (field == nullptr) {
        return false;
    }

    silex::Element theta(*field);
    return theta.gen() && out.add_si(theta, 1);
}

bool check_first_free_generator(const silex::OrderUnitGroup& group,
                                const silex::Element& expected) noexcept {
    const silex::NumberField* field = expected.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           expanded.equal(expected);
}

bool check_first_free_generator_square(const silex::OrderUnitGroup& group,
                                       const silex::Element& expected_square)
        noexcept {
    const silex::NumberField* field = expected_square.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    silex::Element square(*field);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           square.multiply(expanded, expanded) &&
           square.equal(expected_square);
}

bool check_first_free_generator_is_order_unit(
        const silex::OrderUnitGroup& group,
        const silex::Order& order) noexcept {
    const silex::NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    silex::OrderElement order_element(order);
    silex::Ideal principal(order);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           order_element.set_element(expanded) &&
           principal.set_principal(order_element) &&
           principal.is_one();
}

bool set_quadratic_coeffs(silex::Element& out,
                          slong constant_num,
                          ulong constant_den,
                          slong linear_num,
                          ulong linear_den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;

    sflint::fmpq_set_si(coefficient, constant_num, constant_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    sflint::fmpq_set_si(coefficient, linear_num, linear_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool set_integral_quadratic_coeffs_str(silex::Element& out,
                                       const char* constant,
                                       const char* linear) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::Fmpz integer;

    if (!sflint::fmpz_set_str(sflint::FmpzRef(integer), constant)) {
        return false;
    }
    sflint::fmpq_set_fmpz(coefficient, sflint::FmpzConstRef(integer));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    if (!sflint::fmpz_set_str(sflint::FmpzRef(integer), linear)) {
        return false;
    }
    sflint::fmpq_set_fmpz(coefficient, sflint::FmpzConstRef(integer));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool degree_one_residue_root(sflint::Fmpz& root,
                             const silex::PrimeIdeal& prime) noexcept {
    sflint::Fmpz p;
    sflint::FmpzPoly residue;
    sflint::Fmpz constant;
    sflint::Fmpz leading;
    sflint::Fmpz inverse_leading;
    if (prime.residue_degree() != 1 ||
        !prime.rational_prime(sflint::FmpzRef(p)) ||
        !prime.residue_polynomial(sflint::FmpzPolyRef(residue)) ||
        ::fmpz_poly_degree(residue.raw()) != 1) {
        return false;
    }

    ::fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
    ::fmpz_poly_get_coeff_fmpz(leading.raw(), residue.raw(), 1);
    ::fmpz_mod(leading.raw(), leading.raw(), p.raw());
    if (::fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) == 0) {
        return false;
    }
    ::fmpz_neg(root.raw(), constant.raw());
    ::fmpz_mod(root.raw(), root.raw(), p.raw());
    ::fmpz_mul(root.raw(), root.raw(), inverse_leading.raw());
    ::fmpz_mod(root.raw(), root.raw(), p.raw());
    return true;
}

bool check_real_quadratic_group(const silex::OrderUnitGroup& group,
                                const silex::Order& order,
                                const silex::Element& expected_unit) noexcept {
    if (!group.is_defined() || !group.is_set() ||
        !silex::same_order_parent(group.parent(), &order) ||
        group.free_rank() != 1 ||
        group.certification_status() != silex::CertificationMode::proven ||
        !check_first_free_generator(group, expected_unit)) {
        return false;
    }

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    return group.torsion_order(sflint::FmpzRef(torsion_order)) &&
           sflint::fmpz_equal_si(torsion_order, 2) &&
           group.regulator(sflint::ArbRef(regulator)) &&
           sflint::arb_is_positive(regulator);
}

int test_degree_one_rank_zero_group() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderUnitGroup group;
    assert(!group.is_defined());
    assert(!group.is_set());
    assert(group.free_rank() == -1);
    assert(group.certification_status() == silex::CertificationMode::unknown);

    assert(group.compute(order));
    assert(check_rank_zero_group(group, order, 2));

    silex::OrderElement generator(order);
    assert(group.torsion_generator(generator));
    assert(generator.equal_si(-1));
    return 0;
}

int test_compute_real_quadratic_proven() {
    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Order sqrt2_order;
    sqrt2_order = silex::test::equation_order(sqrt2);
    assert(sqrt2_order.is_maximal());
    silex::Element epsilon2(sqrt2);
    assert(set_quadratic_coeffs(epsilon2, 1, 1, 1, 1));

    silex::OrderUnitGroup group;
    assert(group.compute(sqrt2_order));
    assert(check_real_quadratic_group(group, sqrt2_order, epsilon2));

    silex::NumberField sqrt5 = quadratic_field(5);
    silex::Order equation;
    silex::Order maximal(sqrt5);
    equation = silex::test::equation_order(sqrt5);
    assert(!equation.is_maximal());
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    silex::Element epsilon5(sqrt5);
    assert(set_quadratic_coeffs(epsilon5, 1, 2, 1, 2));

    silex::OrderUnitGroup half_integral;
    assert(half_integral.compute(maximal));
    assert(check_real_quadratic_group(half_integral, maximal, epsilon5));

    silex::NumberField sqrt17345 = quadratic_field(17345);
    silex::Order equation17345;
    silex::Order maximal17345(sqrt17345);
    equation17345 = silex::test::equation_order(sqrt17345);
    assert(!equation17345.is_maximal());
    assert(maximal17345.maximal_order(equation17345));
    assert(maximal17345.is_maximal());
    silex::Element epsilon17345(sqrt17345);
    assert(set_integral_quadratic_coeffs_str(
            epsilon17345, "301977958012", "2292915721"));

    silex::OrderUnitGroup large_half_integral;
    assert(large_half_integral.compute(maximal17345));
    assert(check_real_quadratic_group(
            large_half_integral, maximal17345, epsilon17345));

    return 0;
}

int test_compute_with_class_group_rank_zero() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 64;
    options.max_relations = 16;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 80));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(class_group.relation_saturation_record_count() == 0);
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(!class_group.try_certify_with_units(
            units, silex::CertificationMode::grh, 80));
    assert(check_rank_zero_group(units, order, 2));

    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 1));

    silex::ClassGroupContext proven_class_group;
    silex::OrderUnitGroup proven_units;
    options.requested_certification = silex::CertificationMode::proven;
    assert(proven_units.compute_with_class_group(
            proven_class_group, order, sflint::FmpzConstRef(bound),
            options, 80));
    assert(proven_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.relation_saturation_record_count() == 0);
    assert(proven_class_group.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.unit_proof_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(check_rank_zero_group(proven_units, order, 2));

    return 0;
}

int test_compute_with_class_group_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 160));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(units.is_set());
    assert(silex::same_order_parent(units.parent(), &order));
    assert(units.free_rank() == 1);
    silex::Element epsilon(field);
    assert(set_quadratic_coeffs(epsilon, 1, 1, 1, 1));
    assert(check_real_quadratic_group(units, order, epsilon));

    silex::ClassGroupContext proven_class_group;
    silex::OrderUnitGroup proven_units;
    options.requested_certification = silex::CertificationMode::proven;
    assert(proven_units.compute_with_class_group(
            proven_class_group, order, sflint::FmpzConstRef(bound),
            options, 160));
    assert(proven_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.unit_proof_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(check_real_quadratic_group(proven_units, order, epsilon));

    silex::NumberField sqrt5 = quadratic_field(5);
    silex::Order equation;
    silex::Order maximal(sqrt5);
    equation = silex::test::equation_order(sqrt5);
    assert(maximal.maximal_order(equation));
    silex::Element epsilon5(sqrt5);
    assert(set_quadratic_coeffs(epsilon5, 1, 2, 1, 2));

    sflint::Fmpz sqrt5_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(sqrt5_bound), maximal));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(sqrt5_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(sqrt5_bound), 2);
    }
    silex::ClassGroupComputeOptions sqrt5_options;
    sqrt5_options.max_candidates = 256;
    sqrt5_options.max_relations = 48;
    sqrt5_options.requested_certification =
            silex::CertificationMode::proven;

    silex::ClassGroupContext sqrt5_class_group;
    silex::OrderUnitGroup sqrt5_units;
    assert(sqrt5_units.compute_with_class_group(
            sqrt5_class_group, maximal, sflint::FmpzConstRef(sqrt5_bound),
            sqrt5_options, 160));
    assert(sqrt5_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(sqrt5_class_group.relation_kernel_unit_count() >=
           sqrt5_units.free_rank());
    assert(sqrt5_class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(check_real_quadratic_group(sqrt5_units, maximal, epsilon5));

    return 0;
}

int test_compute_with_class_group_real_quadratic_210() {
    silex::NumberField field = quadratic_field(210);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::Element expected_unit(field);
    assert(set_quadratic_coeffs(expected_unit, 29, 1, 2, 1));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 10000;
    options.max_relations = 1000;
    options.requested_certification = silex::CertificationMode::proven;

    auto check_pair = [&](const silex::ClassGroupContext& class_group,
                          const silex::OrderUnitGroup& units) {
        assert(class_group.has_presentation());
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.regulator_proof_status() ==
               silex::ProofState::verified);

        sflint::Fmpz class_order;
        sflint::Fmpz invariant;
        assert(class_group.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 4));
        assert(class_group.invariant_count() == 2);
        assert(class_group.invariant(sflint::FmpzRef(invariant), 0));
        assert(sflint::fmpz_equal_si(invariant, 2));
        assert(class_group.invariant(sflint::FmpzRef(invariant), 1));
        assert(sflint::fmpz_equal_si(invariant, 2));
        assert(check_real_quadratic_group(units, order, expected_unit));
    };

    auto compute_pair = [&](silex::ClassGroupContext& class_group,
                            silex::OrderUnitGroup& units,
                            slong requested_bound, slong precision) {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, requested_bound));
        assert(units.compute_with_class_group(
                class_group, order, sflint::FmpzConstRef(bound), options,
                precision));
        check_pair(class_group, units);
    };

    const slong requested_bounds[] = {10, 19, 20, 50};
    for (const slong requested_bound : requested_bounds) {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        compute_pair(class_group, units, requested_bound, 128);
    }

    silex::ClassGroupContext preserved_class_group;
    silex::OrderUnitGroup preserved_units;
    compute_pair(preserved_class_group, preserved_units, 20, 192);

    const slong relation_count_before =
            preserved_class_group.relation_count();
    const auto relations_before = preserved_class_group.relations();
    const auto invariant_generators_before =
            preserved_class_group.invariant_generator_matrix();
    const auto regulator_before = preserved_units.regulator();
    const auto torsion_order_before = preserved_units.torsion_order();
    assert(relations_before.has_value());
    assert(invariant_generators_before.has_value());
    assert(regulator_before.has_value());
    assert(torsion_order_before.has_value());

    silex::ClassGroupComputeOptions constrained_options = options;
    constrained_options.max_candidates = 0;
    constrained_options.max_relations = 0;
    sflint::Fmpz constrained_bound;
    assert(set_fmpz_si(constrained_bound, 10));
    assert(!preserved_units.compute_with_class_group(
            preserved_class_group, order,
            sflint::FmpzConstRef(constrained_bound), constrained_options,
            128));

    check_pair(preserved_class_group, preserved_units);
    assert(preserved_class_group.relation_count() == relation_count_before);
    const auto relations_after = preserved_class_group.relations();
    const auto invariant_generators_after =
            preserved_class_group.invariant_generator_matrix();
    const auto regulator_after = preserved_units.regulator();
    const auto torsion_order_after = preserved_units.torsion_order();
    assert(relations_after.has_value());
    assert(invariant_generators_after.has_value());
    assert(regulator_after.has_value());
    assert(torsion_order_after.has_value());
    assert(sflint::fmpz_mat_equal(*relations_before, *relations_after));
    assert(sflint::fmpz_mat_equal(*invariant_generators_before,
                                  *invariant_generators_after));
    assert(::arb_equal(regulator_before->raw(), regulator_after->raw()) != 0);
    assert(sflint::fmpz_equal(*torsion_order_before,
                              *torsion_order_after));

    return 0;
}

int test_nonmaximal_quadratic_pair_rejection_preserves_output() {
    silex::NumberField field = quadratic_field(-47);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(!equation.is_maximal());
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), maximal));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext fresh_class_group;
    silex::OrderUnitGroup fresh_units;
    assert(!fresh_units.compute_with_class_group(
            fresh_class_group, equation, sflint::FmpzConstRef(bound),
            options, 128));
    assert(!fresh_class_group.has_factor_base());
    assert(!fresh_class_group.has_presentation());
    assert(!fresh_units.is_set());

    silex::ClassGroupContext preserved_class_group;
    silex::OrderUnitGroup preserved_units;
    assert(preserved_units.compute_with_class_group(
            preserved_class_group, maximal, sflint::FmpzConstRef(bound),
            options, 128));
    assert(preserved_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(preserved_units.certification_status() ==
           silex::CertificationMode::proven);
    assert(check_rank_zero_group(preserved_units, maximal, 2));

    const slong relation_count_before =
            preserved_class_group.relation_count();
    const auto relations_before = preserved_class_group.relations();
    const auto invariant_generators_before =
            preserved_class_group.invariant_generator_matrix();
    const auto class_order_before = preserved_class_group.order();
    const auto regulator_before = preserved_units.regulator();
    const auto torsion_order_before = preserved_units.torsion_order();
    silex::OrderElement torsion_before(maximal);
    assert(relations_before.has_value());
    assert(invariant_generators_before.has_value());
    assert(class_order_before.has_value());
    assert(regulator_before.has_value());
    assert(torsion_order_before.has_value());
    assert(preserved_units.torsion_generator(torsion_before));
    assert(sflint::fmpz_equal_si(*class_order_before, 5));

    assert(!preserved_units.compute_with_class_group(
            preserved_class_group, equation, sflint::FmpzConstRef(bound),
            options, 128));
    assert(silex::same_order_parent(preserved_class_group.parent(),
                                    &maximal));
    assert(silex::same_order_parent(preserved_units.parent(), &maximal));
    assert(preserved_class_group.relation_count() == relation_count_before);
    assert(preserved_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(preserved_units.certification_status() ==
           silex::CertificationMode::proven);

    const auto relations_after = preserved_class_group.relations();
    const auto invariant_generators_after =
            preserved_class_group.invariant_generator_matrix();
    const auto class_order_after = preserved_class_group.order();
    const auto regulator_after = preserved_units.regulator();
    const auto torsion_order_after = preserved_units.torsion_order();
    silex::OrderElement torsion_after(maximal);
    assert(relations_after.has_value());
    assert(invariant_generators_after.has_value());
    assert(class_order_after.has_value());
    assert(regulator_after.has_value());
    assert(torsion_order_after.has_value());
    assert(preserved_units.torsion_generator(torsion_after));
    assert(sflint::fmpz_mat_equal(*relations_before, *relations_after));
    assert(sflint::fmpz_mat_equal(*invariant_generators_before,
                                  *invariant_generators_after));
    assert(sflint::fmpz_equal(*class_order_before, *class_order_after));
    assert(::arb_equal(regulator_before->raw(), regulator_after->raw()) != 0);
    assert(sflint::fmpz_equal(*torsion_order_before,
                              *torsion_order_after));
    assert(torsion_before.equal(torsion_after));

    return 0;
}

int test_compute_with_class_group_quintic_proven() {
    silex::NumberField field = quintic_field(0, 0, -1, -1);
    silex::Order equation;
    equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 1000;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 128));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    // The native reference-style unit extraction can certify this pair without
    // preserving the older incidental surplus of five kernel rows.
    assert(class_group.relation_kernel_unit_count() >= units.free_rank());
    assert(class_group.analytic_class_regulator_status() ==
           silex::ProofState::verified);
    assert(class_group.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(class_group.relation_saturation_status() !=
           silex::ProofState::unavailable);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(class_group.zeta_bf_proof_status() ==
           silex::ProofState::not_checked);

    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_is_one(sflint::FmpzConstRef(class_order)));
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(units.free_rank() == 2);

    return 0;
}

int test_compute_with_class_group_quintic_small_bf_budget_trivial() {
    silex::NumberField field = quintic_field(0, 0, -1, -1);
    silex::Order equation;
    equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 1000;
    options.zeta_bf_max_cutoff = 1;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 128));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(class_group.zeta_bf_proof_status() ==
           silex::ProofState::not_checked);
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(units.free_rank() == 2);

    return 0;
}

int test_compute_with_class_group_failures_preserve_output() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);

    silex::OrderUnitGroup units;
    assert(units.compute(degree_one_order));
    assert(check_rank_zero_group(units, degree_one_order, 2));

    silex::ClassGroupContext class_group(degree_one_order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(class_group.build_factor_base(sflint::FmpzConstRef(bound)));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 64;
    options.max_relations = 16;

    assert(!units.compute_with_class_group(class_group, order,
                                           sflint::FmpzConstRef(bound),
                                           options, 160));
    assert(check_rank_zero_group(units, degree_one_order, 2));
    assert(silex::same_order_parent(class_group.parent(), &degree_one_order));

    assert(!units.compute_with_class_group(class_group, degree_one_order,
                                           sflint::FmpzConstRef(bound),
                                           options, 0));
    assert(check_rank_zero_group(units, degree_one_order, 2));
    assert(silex::same_order_parent(class_group.parent(), &degree_one_order));

    return 0;
}

int test_relation_kernel_rescale_log_matrix_word_rounding() {
    // Source: reference implementation polarit2.c:RgM_rescale_to_int rescales t_REAL
    // matrices using the represented, word-rounded real precision before
    // integer LLL.  Relation-kernel unit selection uses this path for compact
    // log matrices, so keep it aligned with the class-group reference boundary.
    sflint::ArbMat matrix(2, 2);
    sflint::FmpzMat rows32(2, 2);
    sflint::FmpzMat rows64(2, 2);
    sflint::FmpzMat native_rows(2, 2);
    sflint::FmpzMat expected32(2, 2);
    sflint::FmpzMat expected64(2, 2);

    set_arb_mat_entry_si_2exp(matrix, 0, 0, 1, 0);
    add_arb_mat_entry_2exp(matrix, 0, 0, -100);
    set_arb_mat_entry_si_2exp(matrix, 0, 1, 1, -1);
    set_arb_mat_entry_si_2exp(matrix, 1, 0, 1, -66);
    set_arb_mat_entry_si_2exp(matrix, 1, 1, -3, -2);

    set_zmat_entry_str(expected32, 0, 0, "73786976294838206464");
    set_zmat_entry_str(expected32, 0, 1, "36893488147419103232");
    set_zmat_entry_si(expected32, 1, 0, 1);
    set_zmat_entry_str(expected32, 1, 1, "-55340232221128654848");

    set_zmat_entry_str(expected64, 0, 0,
                       "1267650600228229401496703205377");
    set_zmat_entry_str(expected64, 0, 1,
                       "633825300114114700748351602688");
    set_zmat_entry_str(expected64, 1, 0, "17179869184");
    set_zmat_entry_str(expected64, 1, 1,
                       "-950737950171172051122527404032");

    assert(silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows32), sflint::ArbMatConstRef(matrix), 65,
            32));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(rows32),
                                  sflint::FmpzMatConstRef(expected32)));
    assert(silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), 65,
            64));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(rows64),
                                  sflint::FmpzMatConstRef(expected64)));

    assert(!silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), 65,
            16));
    assert(!silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), 0,
            64));
    assert(!silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), -1,
            64));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(rows64),
                                  sflint::FmpzMatConstRef(expected64)));

    sflint::ArbMat empty_matrix(0, 0);
    sflint::FmpzMat empty_rows(0, 0);
    for (const slong word_bits : {slong{32}, slong{64}}) {
        const slong maximum_safe_precision =
                std::numeric_limits<slong>::max() - (word_bits - 1);
        assert(silex::detail::
                       relation_kernel_rescale_log_matrix_for_testing(
                               sflint::FmpzMatRef(empty_rows),
                               sflint::ArbMatConstRef(empty_matrix),
                               maximum_safe_precision, word_bits));
        sflint::FmpzMat overflow_rows(2, 2);
        sflint::fmpz_mat_set(sflint::FmpzMatRef(overflow_rows),
                             sflint::FmpzMatConstRef(expected64));
        assert(!silex::detail::
                        relation_kernel_rescale_log_matrix_for_testing(
                                sflint::FmpzMatRef(overflow_rows),
                                sflint::ArbMatConstRef(matrix),
                                maximum_safe_precision + 1, word_bits));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(overflow_rows),
                sflint::FmpzMatConstRef(expected64)));
    }

    const slong native_word_bits =
            silex::detail::relation_kernel_native_word_bits();
    assert(native_word_bits == 32 || native_word_bits == 64);
    assert(silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(native_rows), sflint::ArbMatConstRef(matrix),
            65, native_word_bits));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(native_rows),
            native_word_bits == 32 ? sflint::FmpzMatConstRef(expected32)
                                   : sflint::FmpzMatConstRef(expected64)));
    return 0;
}

int test_reduced_regulator_reconstruction_correction() {
    constexpr slong precision = 256;

    auto reset_sentinels = [](sflint::Arb& regulator,
                              sflint::FmpzMat& coordinates) noexcept {
        sflint::arb_set_si(regulator, 42);
        coordinates = sflint::FmpzMat(1, 2);
        set_zmat_entry_si(coordinates, 0, 0, 47);
        set_zmat_entry_si(coordinates, 0, 1, -53);
    };
    auto assert_sentinels = [](const sflint::Arb& regulator,
                               const sflint::FmpzMat& coordinates) noexcept {
        assert(::arb_equal_si(regulator.raw(), 42) != 0);
        assert(sflint::fmpz_mat_nrows(coordinates) == 1);
        assert(sflint::fmpz_mat_ncols(coordinates) == 2);
        assert(mat_entry_is_si(coordinates, 0, 0, 47));
        assert(mat_entry_is_si(coordinates, 0, 1, -53));
    };
    auto call = [](sflint::Arb& out,
                   const sflint::ArbMat& coordinates,
                   const sflint::Arb& regulator_multiple,
                   const sflint::Arb& z,
                   sflint::FmpzMat* integer_coordinates) noexcept {
        return silex::detail::reduced_regulator_from_coordinates(
                sflint::ArbRef(out), coordinates,
                sflint::ArbConstRef(regulator_multiple),
                sflint::ArbConstRef(z), integer_coordinates, precision,
                nullptr);
    };

    sflint::ArbMat one_coordinate(1, 1);
    set_arb_mat_entry_si_2exp(one_coordinate, 0, 0, 1, 0);
    sflint::Arb z;
    sflint::arb_set_ui(z, 8);
    sflint::Arb regulator_multiple;
    sflint::Arb regulator_out;
    sflint::FmpzMat integer_out(1, 2);

    sflint::arb_set_ui(regulator_multiple, 1);
    ::arb_mul_2exp_si(regulator_multiple.raw(), regulator_multiple.raw(), -4);
    reset_sentinels(regulator_out, integer_out);
    assert(!call(regulator_out, one_coordinate, regulator_multiple, z,
                 &integer_out));
    assert_sentinels(regulator_out, integer_out);

    sflint::arb_set_ui(regulator_multiple, 1);
    ::arb_mul_2exp_si(regulator_multiple.raw(), regulator_multiple.raw(), -3);
    assert(call(regulator_out, one_coordinate, regulator_multiple, z,
                &integer_out));
    assert(::arb_equal(regulator_out.raw(), regulator_multiple.raw()) != 0);
    assert(sflint::fmpz_mat_nrows(integer_out) == 1);
    assert(sflint::fmpz_mat_ncols(integer_out) == 1);
    assert(mat_entry_is_si(integer_out, 0, 0, 1));

    sflint::arb_set_ui(regulator_multiple, 1);
    ::arb_mul_2exp_si(regulator_multiple.raw(), regulator_multiple.raw(), -2);
    sflint::arb_set_si(regulator_out, 42);
    assert(call(regulator_out, one_coordinate, regulator_multiple, z,
                nullptr));
    assert(::arb_equal(regulator_out.raw(), regulator_multiple.raw()) != 0);

    sflint::Arf lower;
    sflint::Arf upper;
    ::arf_set_ui(lower.raw(), 3);
    ::arf_mul_2exp_si(lower.raw(), lower.raw(), -5);
    ::arf_set_ui(upper.raw(), 5);
    ::arf_mul_2exp_si(upper.raw(), upper.raw(), -5);
    ::arb_set_interval_arf(regulator_multiple.raw(), lower.raw(), upper.raw(),
                           precision);
    reset_sentinels(regulator_out, integer_out);
    assert(!call(regulator_out, one_coordinate, regulator_multiple, z,
                 &integer_out));
    assert_sentinels(regulator_out, integer_out);

    sflint::ArbMat huge_error(1, 1);
    set_arb_mat_entry_si_2exp(huge_error, 0, 0, 2, 0);
    sflint::Fmpz huge_error_exponent;
    sflint::fmpz_set_si(sflint::FmpzRef(huge_error_exponent),
                        std::numeric_limits<slong>::max());
    sflint::fmpz_add_ui(sflint::FmpzRef(huge_error_exponent),
                        sflint::FmpzConstRef(huge_error_exponent), 1);
    assert(!sflint::fmpz_fits_si(
            sflint::FmpzConstRef(huge_error_exponent)));
    ::arb_add_error_2exp_fmpz(
            arb_mat_entry(huge_error.raw(), 0, 0),
            huge_error_exponent.raw());
    assert(::arb_is_finite(arb_mat_entry(huge_error.raw(), 0, 0)) != 0);
    assert(::arb_contains_si(arb_mat_entry(huge_error.raw(), 0, 0), 2) != 0);
    sflint::arb_one(regulator_multiple);
    sflint::arb_one(z);
    reset_sentinels(regulator_out, integer_out);
    assert(!call(regulator_out, huge_error, regulator_multiple, z,
                 &integer_out));
    assert_sentinels(regulator_out, integer_out);

    sflint::Fmpz q;
    sflint::Fmpz k;
    sflint::Fmpz x_denominator;
    sflint::fmpz_one(sflint::FmpzRef(q));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(q), sflint::FmpzConstRef(q), 40);
    sflint::fmpz_add_ui(sflint::FmpzRef(k), sflint::FmpzConstRef(q), 1);
    sflint::fmpz_mul_2exp(sflint::FmpzRef(x_denominator),
                          sflint::FmpzConstRef(q), 1);
    sflint::fmpz_add_ui(sflint::FmpzRef(x_denominator),
                        sflint::FmpzConstRef(x_denominator), 1);

    sflint::Fmpz one;
    sflint::Fmpz two;
    sflint::fmpz_one(sflint::FmpzRef(one));
    sflint::fmpz_set_si(sflint::FmpzRef(two), 2);
    sflint::Fmpq x;
    sflint::Fmpq one_over_q;
    sflint::Fmpq one_over_k;
    sflint::fmpq_set_fmpz_frac(x, sflint::FmpzConstRef(two),
                               sflint::FmpzConstRef(x_denominator));
    sflint::fmpq_set_fmpz_frac(one_over_q, sflint::FmpzConstRef(one),
                               sflint::FmpzConstRef(q));
    sflint::fmpq_set_fmpz_frac(one_over_k, sflint::FmpzConstRef(one),
                               sflint::FmpzConstRef(k));
    sflint::Fmpq crossover;
    ::fmpq_add(crossover.raw(), one_over_q.raw(), one_over_k.raw());
    ::fmpq_div_2exp(crossover.raw(), crossover.raw(), 1);
    sflint::Fmpq midpoint;
    ::fmpq_add(midpoint.raw(), x.raw(), crossover.raw());
    ::fmpq_div_2exp(midpoint.raw(), midpoint.raw(), 1);
    assert(::fmpq_cmp(x.raw(), midpoint.raw()) < 0);
    assert(::fmpq_cmp(midpoint.raw(), crossover.raw()) < 0);

    sflint::ArbMat high_denominator(1, 1);
    arb_struct* high_denominator_entry =
            arb_mat_entry(high_denominator.raw(), 0, 0);
    ::arf_set_fmpq(arb_midref(high_denominator_entry), midpoint.raw(),
                   precision, ARF_RND_NEAR);
    ::mag_set_ui_2exp_si(arb_radref(high_denominator_entry), 1, -79);
    assert(sflint::arb_contains_fmpq(high_denominator_entry, one_over_q));
    assert(sflint::arb_contains_fmpq(high_denominator_entry, one_over_k));
    sflint::arb_set_fmpz(regulator_multiple, sflint::FmpzConstRef(k));
    sflint::arb_set_ui(z, 1);
    ::arb_mul_2exp_si(z.raw(), z.raw(), -1);

    sflint::Fmpq expected_reduced;
    sflint::fmpq_set_fmpz_frac(expected_reduced, sflint::FmpzConstRef(k),
                               sflint::FmpzConstRef(q));
    assert(call(regulator_out, high_denominator, regulator_multiple, z,
                &integer_out));
    assert(sflint::arb_contains_fmpq(sflint::ArbConstRef(regulator_out),
                                     expected_reduced));
    assert(!sflint::arb_contains_si(regulator_out, 1));
    assert(mat_entry_is_si(integer_out, 0, 0, 1));

    ::arb_neg(high_denominator_entry, high_denominator_entry);
    sflint::Fmpq negative_one_over_q;
    sflint::Fmpq negative_one_over_k;
    sflint::fmpq_neg(negative_one_over_q, one_over_q);
    sflint::fmpq_neg(negative_one_over_k, one_over_k);
    assert(sflint::arb_contains_fmpq(high_denominator_entry,
                                     negative_one_over_q));
    assert(sflint::arb_contains_fmpq(high_denominator_entry,
                                     negative_one_over_k));
    assert(call(regulator_out, high_denominator, regulator_multiple, z,
                &integer_out));
    assert(sflint::arb_contains_fmpq(sflint::ArbConstRef(regulator_out),
                                     expected_reduced));
    assert(!sflint::arb_contains_si(regulator_out, 1));
    assert(mat_entry_is_si(integer_out, 0, 0, -1));
    return 0;
}

int test_hnf_regulator_multiple_certification_correction() {
    constexpr slong precision = 256;

    auto call = [](sflint::Arb& out,
                   const sflint::ArbMat& basis) noexcept {
        return silex::detail::hnf_regulator_multiple_from_basis(
                sflint::ArbRef(out), basis, 2, precision);
    };
    auto reset_sentinel = [](sflint::Arb& out) noexcept {
        sflint::arb_set_si(out, 42);
    };
    auto assert_sentinel = [](const sflint::Arb& out) noexcept {
        assert(::arb_equal_si(out.raw(), 42) != 0);
    };
    auto assert_exact_inverse_product = [](const sflint::ArbMat& basis) {
        sflint::ArbMat inverse(2, 2);
        sflint::ArbMat product(2, 2);
        assert(::arb_mat_inv(inverse.raw(), basis.raw(), precision) != 0);
        ::arb_mat_mul(product.raw(), inverse.raw(), basis.raw(), precision);
        for (slong i = 0; i < 2; ++i) {
            for (slong j = 0; j < 2; ++j) {
                assert(::arb_equal_si(
                               arb_mat_entry(product.raw(), i, j),
                               i == j ? 1 : 0) != 0);
            }
        }
    };
    auto set_symmetric_basis = [](sflint::ArbMat& basis,
                                  slong mantissa,
                                  slong exponent) noexcept {
        set_arb_mat_entry_si_2exp(basis, 0, 0, 1, 0);
        set_arb_mat_entry_si_2exp(basis, 1, 0, 1, 0);
        set_arb_mat_entry_si_2exp(basis, 0, 1, -mantissa, exponent);
        set_arb_mat_entry_si_2exp(basis, 1, 1, mantissa, exponent);
    };

    sflint::Arb out;

    sflint::ArbMat passing(2, 2);
    set_symmetric_basis(passing, 1, 0);
    assert_exact_inverse_product(passing);
    reset_sentinel(out);
    assert(call(out, passing));
    assert(::arb_equal_si(out.raw(), 1) != 0);

    // This basis has an exact inverse but kR = 1/2 while its lower-right
    // minor is 1.  Inverse-residual validation alone cannot reject it.
    sflint::ArbMat mismatch(2, 2);
    set_arb_mat_entry_si_2exp(mismatch, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(mismatch, 1, 0, 1, 0);
    set_arb_mat_entry_si_2exp(mismatch, 0, 1, 0, 0);
    set_arb_mat_entry_si_2exp(mismatch, 1, 1, 1, 0);
    assert_exact_inverse_product(mismatch);
    reset_sentinel(out);
    assert(!call(out, mismatch));
    assert_sentinel(out);

    sflint::ArbMat consistent_small(2, 2);
    set_symmetric_basis(consistent_small, 1, -4);
    reset_sentinel(out);
    assert(!call(out, consistent_small));
    assert_sentinel(out);

    sflint::ArbMat threshold(2, 2);
    set_symmetric_basis(threshold, 1, -3);
    reset_sentinel(out);
    assert(call(out, threshold));
    sflint::Arb one_eighth;
    sflint::arb_one(one_eighth);
    ::arb_mul_2exp_si(one_eighth.raw(), one_eighth.raw(), -3);
    assert(::arb_equal(out.raw(), one_eighth.raw()) != 0);

    sflint::ArbMat threshold_overlap(2, 2);
    set_arb_mat_entry_si_2exp(threshold_overlap, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(threshold_overlap, 1, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(
            threshold_overlap, 1, 1, 3, -5, 5, -5);
    ::arb_neg(arb_mat_entry(threshold_overlap.raw(), 0, 1),
              arb_mat_entry(threshold_overlap.raw(), 1, 1));
    reset_sentinel(out);
    assert(!call(out, threshold_overlap));
    assert_sentinel(out);

    // A zero-containing difference ball remains certifiable through its
    // absolute upper bound; it must not be mistaken for threshold overlap.
    sflint::ArbMat rounded(2, 2);
    set_arb_mat_entry_si_2exp(rounded, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(rounded, 1, 0, 1, 0);
    set_arb_mat_entry_si_2exp(rounded, 0, 1, -1, 0);
    ::mag_set_ui_2exp_si(
            arb_radref(arb_mat_entry(rounded.raw(), 0, 1)), 1, -80);
    set_arb_mat_entry_si_2exp(rounded, 1, 1, 1, 0);
    reset_sentinel(out);
    assert(call(out, rounded));
    assert(::arb_contains_si(out.raw(), 1) != 0);
    assert(::arb_is_exact(out.raw()) == 0);
    sflint::Arb one;
    sflint::Arb difference;
    sflint::arb_one(one);
    sflint::arb_sub(difference, one, out, precision);
    assert(::arb_contains_zero(difference.raw()) != 0);
    assert(::arb_is_zero(difference.raw()) == 0);

    return 0;
}

int test_regulator_pivot_exponent_policy() {
    // The source regulator-multiple pivot compares integer exponent bins,
    // preserves the first row on a tie, and accepts only gexpo > -32.
    sflint::ArbMat threshold(1, 1);
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 0, 0);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{-1});
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 1, -32);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{-1});
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 3, -33);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{-1});
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 1, -31);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{0});

    sflint::ArbMat tie(2, 1);
    set_arb_mat_entry_si_2exp(tie, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(tie, 1, 0, 7, -2);
    assert_regulator_pivot_rows(tie, std::vector<slong>{0});
    set_arb_mat_entry_si_2exp(tie, 0, 0, 7, -2);
    set_arb_mat_entry_si_2exp(tie, 1, 0, 1, 0);
    assert_regulator_pivot_rows(tie, std::vector<slong>{0});
    set_arb_mat_entry_si_2exp(tie, 1, 0, 1, 1);
    assert_regulator_pivot_rows(tie, std::vector<slong>{1});

    sflint::ArbMat same_bin_ball(2, 1);
    set_arb_mat_entry_interval_si_2exp(
            same_bin_ball, 0, 0, 5, -2, 3, -1);
    set_arb_mat_entry_si_2exp(same_bin_ball, 1, 0, 7, -2);
    assert_regulator_pivot_rows(same_bin_ball, std::vector<slong>{0});

    sflint::ArbMat huge(2, 1);
    set_arb_mat_entry_huge_2exp(
            huge, 0, 0, "92233720368547758081234567890");
    set_arb_mat_entry_huge_2exp(
            huge, 1, 0, "92233720368547758081234567891");
    assert_regulator_pivot_rows(huge, std::vector<slong>{1});
    return 0;
}

int test_regulator_pivot_arb_ambiguity_and_invalidity() {
    std::vector<slong> pivots{47};
    sflint::ArbMat exact_zero(2, 1);
    set_arb_mat_entry_si_2exp(exact_zero, 0, 0, 0, 0);
    set_arb_mat_entry_si_2exp(exact_zero, 1, 0, 1, 0);
    assert_regulator_pivot_rows(exact_zero, std::vector<slong>{1});

    sflint::ArbMat later_equal(2, 1);
    set_arb_mat_entry_si_2exp(later_equal, 0, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(later_equal, 1, 0, 3, -2, 5, -2);
    assert_regulator_pivot_rows(later_equal, std::vector<slong>{0});

    sflint::ArbMat earlier_equal(2, 1);
    set_arb_mat_entry_interval_si_2exp(earlier_equal, 0, 0, 3, -2, 5, -2);
    set_arb_mat_entry_si_2exp(earlier_equal, 1, 0, 1, 0);
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(earlier_equal), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat later_outranks(2, 1);
    set_arb_mat_entry_si_2exp(later_outranks, 0, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(
            later_outranks, 1, 0, 3, -1, 5, -1);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(later_outranks), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat contains_zero(2, 1);
    set_arb_mat_entry_si_2exp(contains_zero, 0, 0, 1, 0);
    ::arb_zero(arb_mat_entry(contains_zero.raw(), 1, 0));
    ::arb_add_error_2exp_si(arb_mat_entry(contains_zero.raw(), 1, 0), -2);
    assert_regulator_pivot_rows(contains_zero, std::vector<slong>{0});
    ::arb_add_error_2exp_si(arb_mat_entry(contains_zero.raw(), 1, 0), 2);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(contains_zero), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat earlier_contains_zero(2, 1);
    ::arb_zero(arb_mat_entry(earlier_contains_zero.raw(), 0, 0));
    ::arb_add_error_2exp_si(
            arb_mat_entry(earlier_contains_zero.raw(), 0, 0), 0);
    set_arb_mat_entry_si_2exp(earlier_contains_zero, 1, 0, 1, 0);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(earlier_contains_zero), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat threshold_crossing(1, 1);
    set_arb_mat_entry_interval_si_2exp(
            threshold_crossing, 0, 0, 3, -33, 1, -31);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(threshold_crossing), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat below_threshold(1, 1);
    set_arb_mat_entry_interval_si_2exp(
            below_threshold, 0, 0, 3, -35, 5, -34);
    assert_regulator_pivot_rows(below_threshold, std::vector<slong>{-1});

    sflint::ArbMat nonfinite_tail(1, 2);
    set_arb_mat_entry_si_2exp(nonfinite_tail, 0, 0, 1, 0);
    ::arb_pos_inf(arb_mat_entry(nonfinite_tail.raw(), 0, 1));
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(nonfinite_tail), 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(pivots.empty());

    // The first column is certified and pivoted before the second column's
    // threshold-crossing ball makes the result inconclusive.  No partial
    // pivot trace may escape.
    sflint::ArbMat late_ambiguity(2, 2);
    set_arb_mat_entry_si_2exp(late_ambiguity, 0, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(
            late_ambiguity, 1, 1, 3, -33, 1, -31);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(late_ambiguity), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(exact_zero), 0) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(pivots.empty());
    return 0;
}

int test_regulator_pivot_clean_cols_and_mapping() {
    const silex::Signature signature(2, 0);
    std::vector<slong> selected{47};
    std::vector<slong> pivots{47};

    sflint::ArbMat discarded(1, 2);
    set_arb_mat_entry_si_2exp(discarded, 0, 0, 1, -3);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, discarded, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(selected.empty());
    assert((pivots == std::vector<slong>{0}));

    sflint::ArbMat retained(1, 2);
    set_arb_mat_entry_si_2exp(retained, 0, 0, 1, -2);
    set_arb_mat_entry_si_2exp(retained, 0, 1, 1, 0);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, retained, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert((selected == std::vector<slong>{0}));

    sflint::ArbMat ambiguous(1, 2);
    set_arb_mat_entry_interval_si_2exp(ambiguous, 0, 0, 3, -4, 5, -4);
    set_arb_mat_entry_si_2exp(ambiguous, 0, 1, 0, 0);
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, ambiguous, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(selected.empty());
    assert(pivots.empty());

    ::arb_zero(arb_mat_entry(ambiguous.raw(), 0, 0));
    ::arb_add_error_2exp_si(arb_mat_entry(ambiguous.raw(), 0, 0), -2);
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, ambiguous, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(selected.empty());
    assert(pivots.empty());

    // An earlier possible retention is resolved by a later certified entry.
    sflint::ArbMat resolved(1, 2);
    set_arb_mat_entry_interval_si_2exp(resolved, 0, 0, 3, -4, 5, -4);
    set_arb_mat_entry_si_2exp(resolved, 0, 1, 1, 0);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, resolved, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert((selected == std::vector<slong>{0}));

    // Certified retain cannot hide a later nonfinite entry.
    set_arb_mat_entry_si_2exp(resolved, 0, 0, 1, -2);
    ::arb_pos_inf(arb_mat_entry(resolved.raw(), 0, 1));
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, resolved, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(selected.empty());
    assert(pivots.empty());

    // Invalid in a later unit row dominates earlier cleaning uncertainty.
    sflint::ArbMat late_invalid(2, 2);
    set_arb_mat_entry_interval_si_2exp(
            late_invalid, 0, 0, 3, -4, 5, -4);
    ::arb_pos_inf(arb_mat_entry(late_invalid.raw(), 1, 1));
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, late_invalid, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(selected.empty());
    assert(pivots.empty());

    // Two retained originals separated by a discarded row keep source order.
    const silex::Signature signature3(3, 0);
    sflint::ArbMat mapped(3, 3);
    set_arb_mat_entry_si_2exp(mapped, 0, 0, 1, -2);
    set_arb_mat_entry_si_2exp(mapped, 0, 1, 1, 0);
    set_arb_mat_entry_si_2exp(mapped, 1, 0, 1, -3);
    set_arb_mat_entry_si_2exp(mapped, 2, 1, 1, -2);
    set_arb_mat_entry_si_2exp(mapped, 2, 2, 1, 0);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, mapped, signature3, 2, 3, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert((selected == std::vector<slong>{0, 2}));

    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, mapped, silex::Signature(2, 0), 2, 3, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(selected.empty());
    assert(pivots.empty());
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, mapped,
            silex::Signature(std::numeric_limits<slong>::max(),
                             std::numeric_limits<slong>::max()),
            2, 3, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);

    // reference's rank-zero branch bypasses filter_unit_log_columns, so a finite ambiguous log
    // is irrelevant rather than precision-inconclusive.
    sflint::ArbMat rank_zero_ambiguous(1, 1);
    ::arb_zero(arb_mat_entry(rank_zero_ambiguous.raw(), 0, 0));
    ::arb_add_error_2exp_si(
            arb_mat_entry(rank_zero_ambiguous.raw(), 0, 0), -2);
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, rank_zero_ambiguous, silex::Signature(1, 0),
            0, 1, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(selected.empty());
    assert(pivots.empty());
    ::arb_pos_inf(arb_mat_entry(rank_zero_ambiguous.raw(), 0, 0));
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, rank_zero_ambiguous, silex::Signature(1, 0),
            0, 1, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(selected.empty());
    assert(pivots.empty());
    return 0;
}

int test_regulator_pivot_multi_column_elimination() {
    sflint::ArbMat matrix(2, 2);
    set_arb_mat_entry_si_2exp(matrix, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(matrix, 0, 1, 1, 0);
    set_arb_mat_entry_si_2exp(matrix, 1, 0, 7, -2);
    set_arb_mat_entry_si_2exp(matrix, 1, 1, 3, 0);
    sflint::ArbMat input_snapshot(2, 2);
    ::arb_mat_set(input_snapshot.raw(), matrix.raw());
    assert_regulator_pivot_rows(matrix, std::vector<slong>{0, 1});
    assert(::arb_mat_equal(input_snapshot.raw(), matrix.raw()) != 0);
    assert_regulator_pivot_rows(matrix, std::vector<slong>{0, 1});
    assert(::arb_mat_equal(input_snapshot.raw(), matrix.raw()) != 0);
    return 0;
}

void assert_regulator_pivot_route_state_preserved(
        const silex::detail::RegulatorPivotFinishTestResult& result)
        noexcept {
    assert(result.goal_reached_after == result.goal_reached_before);
    assert(result.relation_need_after == result.relation_need_before);
    assert(result.completion_old_need_after ==
           result.completion_old_need_before);
    assert(result.completion_dependent_trials_after ==
           result.completion_dependent_trials_before);
    assert(result.completion_subfactor_base_trials_after ==
           result.completion_subfactor_base_trials_before);
    assert(result.analytic_extra_relation_requests_after ==
           result.analytic_extra_relation_requests_before);
    assert(result.finish_unit_log_rotation_after ==
           result.finish_unit_log_rotation_before);
    assert(result.squash_index_after == result.squash_index_before);
    assert(result.candidates_tried_after == result.candidates_tried_before);
    assert(result.accepted_relations_after ==
           result.accepted_relations_before);
    assert(result.finish_unit_log_rotation_active_after ==
           result.finish_unit_log_rotation_active_before);
    assert(result.finish_full_rank_relation_active_after ==
           result.finish_full_rank_relation_active_before);
    assert(result.relation_control_state_unchanged);
}

int test_regulator_pivot_finish_route_matrix() {
    using Action = silex::detail::RegulatorPivotFinishAction;
    using Outcome = silex::detail::RegulatorPivotOutcome;
    using Path = silex::detail::RegulatorPivotFinishPath;
    for (const Path path : {Path::independent_unit_preprobe,
                            Path::regulator_product}) {
        const auto full =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::success, 3, 3, true, 160, 0);
        assert(full.action == Action::proceed);
        assert(full.all_control_state_unchanged);

        const auto partial =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::success, 1, 3, true, 160, 0);
        assert(partial.action == Action::request_relations);
        assert(!partial.goal_reached_after);
        assert(partial.relation_need_after == 2);
        assert(partial.completion_old_need_after == 2);
        assert(partial.completion_dependent_trials_after == 0);
        assert(partial.completion_subfactor_base_trials_after ==
               partial.completion_subfactor_base_trials_before);
        assert(partial.candidates_tried_after ==
               partial.candidates_tried_before);
        assert(partial.accepted_relations_after ==
               partial.accepted_relations_before);
        assert(partial.analytic_extra_relation_requests_after ==
               partial.analytic_extra_relation_requests_before + 1);
        assert(partial.finish_unit_log_rotation_active_after);
        assert(partial.finish_unit_log_rotation_after ==
               partial.squash_index_before);
        assert(partial.squash_index_after == partial.squash_index_before + 1);
        assert(!partial.finish_full_rank_relation_active_after);
        assert(partial.analytic_finish_precision_after ==
               partial.analytic_finish_precision_before);
        assert(partial.analytic_precision_doublings_after ==
               partial.analytic_precision_doublings_before);
        assert(partial.analytic_precision_inconclusive_after ==
               partial.analytic_precision_inconclusive_before);
        assert(partial.factor_base_restart_requests_after ==
               partial.factor_base_restart_requests_before);
        assert(partial.factor_base_restart_pending_after ==
               partial.factor_base_restart_pending_before);
        assert(partial.factor_base_restart_allow_past_half_after ==
               partial.factor_base_restart_allow_past_half_before);

        const auto precision =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::precision_inconclusive, 0, 3,
                                false, 160, 0);
        assert(precision.action == Action::retry_precision);
        assert_regulator_pivot_route_state_preserved(precision);
        assert(precision.analytic_finish_precision_after == 320);
        assert(precision.analytic_precision_doublings_after ==
               precision.analytic_precision_doublings_before + 1);
        assert(precision.analytic_precision_inconclusive_after ==
               precision.analytic_precision_inconclusive_before + 1);
        assert(precision.factor_base_restart_requests_after ==
               precision.factor_base_restart_requests_before);
        assert(!precision.factor_base_restart_pending_after);
        assert(precision.factor_base_restart_allow_past_half_after ==
               precision.factor_base_restart_allow_past_half_before);

        const auto precision_cap =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::precision_inconclusive, 0, 3,
                                false, 40960, 0);
        assert(precision_cap.action == Action::failed);
        assert(precision_cap.all_control_state_unchanged);

        const auto restart =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::precision_inconclusive, 0, 3,
                                true, 160, 7);
        assert(restart.action == Action::restart_factor_base);
        assert_regulator_pivot_route_state_preserved(restart);
        assert(restart.analytic_finish_precision_after ==
               restart.analytic_finish_precision_before);
        assert(restart.analytic_precision_doublings_after ==
               restart.analytic_precision_doublings_before);
        assert(restart.analytic_precision_inconclusive_after ==
               restart.analytic_precision_inconclusive_before);
        assert(restart.factor_base_restart_requests_after ==
               restart.factor_base_restart_requests_before + 1);
        assert(restart.factor_base_restart_pending_after);
        assert(!restart.factor_base_restart_allow_past_half_after);

        const auto invalid =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::invalid, 0, 3, true, 160, 7);
        assert(invalid.action == Action::failed);
        assert(invalid.all_control_state_unchanged);
    }
    return 0;
}

int test_relation_kernel_units_rank_one() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.relation_kernel_unit_count() >= 6);

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    assert(group.set_relation_kernel_units(equation_order, context,
                                           embeddings, 160));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    assert(group.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
    auto owned_torsion_order = group.torsion_order();
    assert(owned_torsion_order.has_value());
    assert(sflint::fmpz_equal_si(*owned_torsion_order, 2));
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_positive(regulator));
    auto owned_regulator = group.regulator();
    assert(owned_regulator.has_value());
    assert(sflint::arb_is_positive(*owned_regulator));

    return 0;
}

int test_relation_kernel_units_bounded_rank_one() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz factor_bound;
    assert(set_fmpz_si(factor_bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() >= 6);

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    sflint::Fmpz denominator_bound;
    assert(set_fmpz_si(denominator_bound, 2));
    assert(group.set_relation_kernel_units_bounded(
            equation_order, context, embeddings,
            sflint::FmpzConstRef(denominator_bound), 80, 160));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    silex::OrderUnitGroup bound_one;
    assert(set_fmpz_si(denominator_bound, 1));
    assert(bound_one.set_relation_kernel_units_bounded(
            equation_order, context, embeddings,
            sflint::FmpzConstRef(denominator_bound), 80, 160));
    assert(bound_one.is_set());
    assert(bound_one.free_rank() == 1);
    assert(check_first_free_generator_is_order_unit(bound_one,
                                                    equation_order));

    return 0;
}

int test_relation_kernel_units_index_bounded_rank_one() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz factor_bound;
    assert(set_fmpz_si(factor_bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() >= 6);

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    assert(group.set_relation_kernel_units_index_bounded(
            equation_order, context, embeddings, 80, 160));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    return 0;
}

int test_relation_kernel_units_index_bounded_quartic_power_root() {
    silex::NumberField field = quartic_field(-1, -1, 0, 2);
    silex::Order equation_order;
    silex::Order order(field);
    equation_order = silex::test::equation_order(field);
    assert(order.maximal_order(equation_order));

    sflint::Fmpz factor_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_bound), order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_bound), 2);
    }

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() == 6);

    silex::EmbeddingContext embeddings(field);
    sflint::Arb analytic_hR;
    assert(silex::zeta_class_regulator_product(
            sflint::ArbRef(analytic_hR), order, 128));

    silex::OrderUnitGroup direct(order);
    sflint::Fmpz direct_bound;
    assert(direct.set_relation_kernel_units(order, context, embeddings, 128));
    assert(direct.class_regulator_index_bound(
            sflint::FmpzRef(direct_bound), context,
            sflint::ArbConstRef(analytic_hR), 128));
    assert(sflint::fmpz_equal_si(direct_bound, 11));

    silex::OrderUnitGroup index_bounded(order);
    sflint::Fmpz bounded_value;
    assert(index_bounded.set_relation_kernel_units_index_bounded(
            order, context, embeddings, 128, 128));
    assert(index_bounded.class_regulator_index_bound(
            sflint::FmpzRef(bounded_value), context,
            sflint::ArbConstRef(analytic_hR), 128));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(bounded_value),
                              sflint::FmpzConstRef(direct_bound)));

    return 0;
}

int test_relation_kernel_units_index_bounded_cubic_rank_target() {
    silex::NumberField field = cubic_field(-2, -5);
    silex::Order equation_order;
    silex::Order order(field);
    equation_order = silex::test::equation_order(field);
    assert(order.maximal_order(equation_order));

    sflint::Fmpz factor_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_bound), order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_bound), 2);
    }

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() == 1);

    sflint::Fmpz aux_start;
    sflint::Fmpz aux_max;
    assert(set_fmpz_si(aux_start, 2));
    assert(set_fmpz_si(aux_max, 31));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup refined(order);
    bool changed = false;
    bool stable = false;
    assert(refined.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, context, embeddings, 128, 128, 1,
            sflint::FmpzConstRef(aux_start),
            sflint::FmpzConstRef(aux_max), 2));

    sflint::Fmpz index_bound;
    assert(refined.regulator_index_bound(sflint::FmpzRef(index_bound), 128));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(index_bound), 74));

    return 0;
}

int test_relation_kernel_units_index_bounded_saturated() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz factor_bound;
    assert(set_fmpz_si(factor_bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() >= 6);

    sflint::Fmpz start;
    sflint::Fmpz max;
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    bool changed = true;
    bool stable = true;
    assert(group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, equation_order, context, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 3));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    return 0;
}

int test_relation_kernel_units_rank_zero_delegates() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    assert(group.set_relation_kernel_units(order, context, embeddings, 80));
    assert(check_rank_zero_group(group, order, 2));

    silex::OrderUnitGroup bounded;
    assert(bounded.set_relation_kernel_units_bounded(
            order, context, embeddings, sflint::FmpzConstRef(bound), 40, 80));
    assert(check_rank_zero_group(bounded, order, 2));

    silex::OrderUnitGroup index_bounded;
    assert(index_bounded.set_relation_kernel_units_index_bounded(
            order, context, embeddings, 40, 80));
    assert(check_rank_zero_group(index_bounded, order, 2));

    silex::OrderUnitGroup saturated;
    bool changed = true;
    bool stable = false;
    assert(saturated.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, context, embeddings, 40, 80, 1,
            sflint::FmpzConstRef(bound), sflint::FmpzConstRef(bound), 2));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(saturated, order, 2));

    return 0;
}

int test_relation_kernel_unit_failures_preserve_output() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);

    silex::OrderUnitGroup group;
    assert(group.compute(degree_one_order));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    silex::ClassGroupContext insufficient(order);
    assert(insufficient.build_factor_base(sflint::FmpzConstRef(bound)));
    silex::EmbeddingContext embeddings(field);
    assert(!group.set_relation_kernel_units(order, insufficient, embeddings,
                                            160));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    sflint::Fmpz denominator_bound;
    assert(set_fmpz_si(denominator_bound, 0));
    assert(!group.set_relation_kernel_units_bounded(
            order, insufficient, embeddings,
            sflint::FmpzConstRef(denominator_bound), 80, 160));
    assert(set_fmpz_si(denominator_bound, 2));
    assert(!group.set_relation_kernel_units_bounded(
            order, insufficient, embeddings,
            sflint::FmpzConstRef(denominator_bound), 160, 80));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    assert(!group.set_relation_kernel_units_index_bounded(
            order, insufficient, embeddings, 80, 160));
    assert(!group.set_relation_kernel_units_index_bounded(
            order, insufficient, embeddings, 160, 80));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    silex::NumberField other_field = quadratic_field(3);
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);
    silex::ClassGroupContext wrong_parent(other_order);
    silex::EmbeddingContext other_embeddings(other_field);
    assert(!group.set_relation_kernel_units(order, wrong_parent, embeddings,
                                            160));
    assert(!group.set_relation_kernel_units(order, insufficient,
                                            other_embeddings, 160));
    assert(!group.set_relation_kernel_units(order, insufficient, embeddings,
                                            0));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    sflint::Fmpz start;
    sflint::Fmpz max;
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));
    bool changed = true;
    bool stable = false;
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(changed);
    assert(!stable);
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, other_embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 160, 80, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 0,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(set_fmpz_si(start, 1));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(set_fmpz_si(start, 5));
    assert(set_fmpz_si(max, 4));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 0));
    assert(changed);
    assert(!stable);
    assert(check_rank_zero_group(group, degree_one_order, 2));

    return 0;
}

int test_saturate_row_real_quadratic_square() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));

    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::FmpzMat row(1, 1);
    sflint::Fmpz ell;
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    assert(set_fmpz_si(ell, 2));

    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    assert(saturated.saturate_row(changed, group,
                                  sflint::FmpzMatConstRef(row), 0,
                                  sflint::FmpzConstRef(ell), embeddings,
                                  128));
    assert(changed);
    assert(saturated.free_rank() == 1);
    assert(saturated.certification_status() ==
           silex::CertificationMode::unknown);
    assert(check_first_free_generator_square(saturated, epsilon2));

    return 0;
}

int test_saturate_row_no_root_and_divisible_copy() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::FmpzMat row(1, 1);
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);

    silex::OrderUnitGroup out(order);
    bool changed = true;
    assert(out.saturate_row(changed, group, sflint::FmpzMatConstRef(row), 0,
                            sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!changed);
    assert(check_first_free_generator(out, epsilon));

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 2);
    changed = true;
    assert(group.saturate_row(changed, group, sflint::FmpzMatConstRef(row), 0,
                              sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!changed);
    assert(check_first_free_generator(group, epsilon));

    return 0;
}

int test_saturate_row_rank_zero_and_failures() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);

    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));
    assert(check_rank_zero_group(rank_zero, degree_one_order, 2));

    sflint::FmpzMat zero_width(1, 0);
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    assert(copied.saturate_row(changed, rank_zero,
                               sflint::FmpzMatConstRef(zero_width), 0,
                               sflint::FmpzConstRef(ell),
                               degree_one_embeddings, 80));
    assert(!changed);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::FmpzMat row(1, 1);
    sflint::FmpzMat bad_width(1, 2);
    sflint::Fmpz not_prime;
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    assert(set_fmpz_si(not_prime, 4));
    assert(!copied.saturate_row(changed, group, sflint::FmpzMatConstRef(row),
                                1, sflint::FmpzConstRef(ell), embeddings,
                                128));
    assert(!copied.saturate_row(changed, group,
                                sflint::FmpzMatConstRef(bad_width), 0,
                                sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!copied.saturate_row(changed, group, sflint::FmpzMatConstRef(row),
                                0, sflint::FmpzConstRef(not_prime),
                                embeddings, 128));
    assert(!copied.saturate_row(changed, group, sflint::FmpzMatConstRef(row),
                                0, sflint::FmpzConstRef(ell), embeddings, 0));
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    return 0;
}

int test_residue_dlog_kernel_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 5));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());
    sflint::FmpzMat kernel(0, 0);
    assert(group.residue_dlog_kernel(kernel, prime_span,
                                     sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(kernel) == 1);
    assert(sflint::fmpz_mat_nrows(kernel) == 1);
    assert(mat_entry_is_si(kernel, 0, 0, 1));

    return 0;
}

int test_residue_dlog_proof_kernel_torsion() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::OrderUnitGroup group;
    assert(group.compute(order));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 3));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    sflint::FmpzMat kernel(1, 1);
    assert(group.residue_dlog_proof_kernel(
            kernel, silex::PrimeIdealSpan(primes.at(0), primes.size()),
            sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 1);

    assert(set_fmpz_si(ell, 3));
    assert(group.residue_dlog_proof_kernel(
            kernel, silex::PrimeIdealSpan(primes.at(0), primes.size()),
            sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 0);

    return 0;
}

int test_rank_one_degree_one_root_image_matches_prime_ideal() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));

    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::FactoredElement stored_generator(field);
    assert(group.free_generator(stored_generator, 0));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(ell, 2));

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p), 1));
    assert(primes.size() > 0);

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);

        bool via_prime = false;
        bool via_root = false;
        bool via_nmod_root = false;
        sflint::Fmpz root;
        assert(degree_one_residue_root(root, *prime));
        assert(silex::detail::saturation_proof_prime_known_rank_one_free_generator_nonzero(
                via_prime, group, stored_generator, *prime,
                sflint::FmpzConstRef(ell)));
        assert(silex::detail::rank_one_factored_image_nonzero_at_degree_one_root(
                via_root, stored_generator, sflint::FmpzConstRef(p),
                sflint::FmpzConstRef(root), sflint::FmpzConstRef(ell)));
        assert(silex::detail::rank_one_factored_image_nonzero_at_degree_one_root_nmod(
                via_nmod_root, stored_generator, sflint::fmpz_get_ui(
                                                        sflint::FmpzConstRef(p)),
                sflint::fmpz_get_ui(sflint::FmpzConstRef(root)),
                sflint::FmpzConstRef(ell)));
        assert(via_root == via_prime);
        assert(via_nmod_root == via_prime);
        ++checked;
    }
    assert(checked == primes.size());
    return 0;
}

int test_direct_degree_one_proof_column_matches_residue_field() {
    silex::OrderUnitGroup group = real_quadratic_unit_group();
    const silex::Order* order = group.parent();
    assert(order != nullptr);
    if (order == nullptr) {
        return 1;
    }

    sflint::Fmpz p;
    sflint::Fmpz ell;
    const bool prime_set = set_fmpz_si(p, 7) && set_fmpz_si(ell, 2);
    assert(prime_set);
    if (!prime_set) {
        return 1;
    }

    silex::PrimeIdealList primes;
    const bool decomposed = silex::decompose_prime(
            primes, *order, sflint::FmpzConstRef(p), 1);
    assert(decomposed);
    assert(primes.size() > 0);
    if (!decomposed || primes.size() == 0) {
        return 1;
    }

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        if (prime == nullptr) {
            return 1;
        }

        sflint::FmpzMat generic_column(0, 0);
        sflint::FmpzMat direct_column(0, 0);
        const bool columns_match =
                silex::detail::saturation_proof_prime_column(
                        generic_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                silex::detail::saturation_proof_prime_column_direct_degree_one(
                        direct_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(generic_column),
                        sflint::FmpzMatConstRef(direct_column));
        assert(columns_match);
        if (!columns_match) {
            return 1;
        }
        ++checked;
    }
    assert(checked == primes.size());
    return checked == primes.size() ? 0 : 1;
}

bool generic_residue_dlog_column_for_test(sflint::FmpzMat& out,
                                          const silex::OrderUnitGroup& group,
                                          const silex::PrimeIdeal& prime,
                                          sflint::FmpzConstRef ell) noexcept {
    const silex::Order* order = group.parent();
    const silex::NumberField* field =
            order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr) {
        return false;
    }

    silex::ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    sflint::FmpzMat candidate(rank, 1);
    silex::ResidueFieldElement image(residue_field);
    silex::ResidueFieldQuotientLog quotient_log(residue_field);
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        silex::FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !quotient_log.apply(sflint::fmpz_mat_entry(candidate, i, 0),
                                image)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

int test_direct_degree_one_residue_dlog_matrix_matches_residue_field() {
    silex::OrderUnitGroup group = real_quadratic_unit_group();
    const silex::Order* order = group.parent();
    assert(order != nullptr);
    if (order == nullptr) {
        return 1;
    }

    sflint::Fmpz p;
    sflint::Fmpz ell;
    const bool prime_set = set_fmpz_si(p, 7) && set_fmpz_si(ell, 2);
    assert(prime_set);
    if (!prime_set) {
        return 1;
    }

    silex::PrimeIdealList primes;
    const bool decomposed = silex::decompose_prime(
            primes, *order, sflint::FmpzConstRef(p), 1);
    assert(decomposed);
    assert(primes.size() > 0);
    if (!decomposed || primes.size() == 0) {
        return 1;
    }

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        if (prime == nullptr) {
            return 1;
        }

        sflint::FmpzMat generic_column(0, 0);
        sflint::FmpzMat direct_column(0, 0);
        sflint::FmpzMat saturation_column(0, 0);
        const bool columns_match =
                generic_residue_dlog_column_for_test(
                        generic_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                silex::detail::residue_dlog_matrix_direct_degree_one(
                        direct_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                silex::detail::saturation_prime_column(
                        saturation_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(generic_column),
                        sflint::FmpzMatConstRef(direct_column)) &&
                sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(generic_column),
                        sflint::FmpzMatConstRef(saturation_column));
        assert(columns_match);
        if (!columns_match) {
            return 1;
        }
        ++checked;
    }
    assert(checked == primes.size());
    return checked == primes.size() ? 0 : 1;
}

int test_half_integral_direct_degree_one_residue_dlog_matches_residue_field() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup group;
    assert(group.compute(maximal));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 11));
    assert(set_fmpz_si(ell, 2));

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(
            primes, maximal, sflint::FmpzConstRef(p), 1));
    assert(primes.size() == 2);

    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);

        sflint::FmpzMat generic_column(0, 0);
        sflint::FmpzMat direct_column(0, 0);
        sflint::FmpzMat saturation_column(0, 0);
        assert(generic_residue_dlog_column_for_test(
                generic_column, group, *prime,
                sflint::FmpzConstRef(ell)));
        assert(silex::detail::residue_dlog_matrix_direct_degree_one(
                direct_column, group, *prime,
                sflint::FmpzConstRef(ell)));
        assert(silex::detail::saturation_prime_column(
                saturation_column, group, *prime,
                sflint::FmpzConstRef(ell)));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(generic_column),
                sflint::FmpzMatConstRef(direct_column)));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(generic_column),
                sflint::FmpzMatConstRef(saturation_column)));
    }

    return 0;
}

int test_select_saturation_proof_primes() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 31));

    silex::PrimeIdealList selected;
    sflint::FmpzMat kernel(1, 1);
    bool certified = false;
    assert(group.select_saturation_proof_primes(
            selected, certified, kernel, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(bound)));
    assert(certified);
    assert(selected.size() > 0);
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 2);
    for (slong i = 0; i < selected.size(); ++i) {
        assert(selected.at(i) != nullptr);
        assert(silex::same_order_parent(selected.at(i)->parent(), &order));
        assert(selected.at(i)->residue_degree() == 1);
    }

    assert(set_fmpz_si(bound, 2));
    certified = false;
    sflint::FmpzMat preserved_kernel(1, 1);
    silex::PrimeIdealList preserved;
    assert(group.select_saturation_proof_primes(
            preserved, certified, preserved_kernel,
            sflint::FmpzConstRef(ell), sflint::FmpzConstRef(bound)) == false);
    assert(!certified);
    assert(!preserved.is_defined());
    assert(sflint::fmpz_mat_nrows(preserved_kernel) == 1);
    assert(sflint::fmpz_mat_ncols(preserved_kernel) == 1);

    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));
    assert(set_fmpz_si(ell, 3));
    assert(set_fmpz_si(bound, 5));
    assert(rank_zero.select_saturation_proof_primes(
            preserved, certified, preserved_kernel,
            sflint::FmpzConstRef(ell), sflint::FmpzConstRef(bound)));
    assert(certified);
    assert(preserved.size() == 0);
    assert(sflint::fmpz_mat_nrows(preserved_kernel) == 0);
    assert(sflint::fmpz_mat_ncols(preserved_kernel) == 0);

    return 0;
}

int test_unit_proof_records() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::OrderUnitGroup group;
    assert(group.compute(order));
    assert(group.unit_proof_record_count() == 0);

    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(aux_bound, 11));
    assert(!group.unit_proof_verified(sflint::FmpzConstRef(ell)));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_local_saturated(
            status, changed, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 80));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.unit_proof_record_count() == 1);
    assert(proved.unit_proof_verified(sflint::FmpzConstRef(ell)));
    assert(check_unit_proof_record(proved, 0, 2,
                                   silex::ProofState::verified, 11, -1,
                                   false));
    auto proof_record = proved.unit_proof_record(0);
    assert(proof_record.has_value());
    assert(sflint::fmpz_equal_si(proof_record->ell, 2));
    assert(proof_record->status == silex::ProofState::verified);
    assert(sflint::fmpz_equal_si(proof_record->aux_prime_bound, 11));
    assert(proof_record->local_primes > 0);
    assert(!proof_record->changed);

    assert(set_fmpz_si(aux_bound, 13));
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_local_saturated(
            status, changed, proved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 80));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.unit_proof_record_count() == 1);
    assert(check_unit_proof_record(proved, 0, 2,
                                   silex::ProofState::verified, 13, -1,
                                   false));

    assert(set_fmpz_si(ell, 3));
    assert(set_fmpz_si(aux_bound, 17));
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_local_saturated(
            status, changed, proved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 80));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.unit_proof_record_count() == 2);
    assert(check_unit_proof_record(proved, 1, 3,
                                   silex::ProofState::verified, 17, 0,
                                   false));

    silex::OrderUnitGroup copied(order);
    assert(copied.set(proved));
    assert(copied.unit_proof_record_count() == 2);
    assert(copied.unit_proof_verified(sflint::FmpzConstRef(ell)));

    silex::OrderUnitGroup empty(order);
    assert(empty.compute(order));
    assert(empty.unit_proof_record_count() == 0);
    copied.swap(empty);
    assert(copied.unit_proof_record_count() == 0);
    assert(empty.unit_proof_record_count() == 2);
    assert(empty.unit_proof_verified(sflint::FmpzConstRef(ell)));

    assert(set_fmpz_si(ell, 5));
    assert(!copied.unit_proof_verified(sflint::FmpzConstRef(ell)));
    slong local_primes = 17;
    status = silex::ProofState::verified;
    changed = true;
    assert(!copied.unit_proof_record(sflint::FmpzRef(ell), status,
                                     sflint::FmpzRef(aux_bound),
                                     local_primes, changed, 99));
    assert(sflint::fmpz_equal_si(ell, 5));
    assert(sflint::fmpz_equal_si(aux_bound, 17));
    assert(status == silex::ProofState::verified);
    assert(local_primes == 17);
    assert(changed);
    assert(!copied.unit_proof_record(99).has_value());

    return 0;
}

int test_prove_local_saturated_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order,
                           silex::FactoredElementSpan(primitive_generators,
                                                      1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 31));

    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_local_saturated(
            status, changed, group, sflint::FmpzConstRef(ell), 2,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(check_first_free_generator(proved, epsilon));
    assert(proved.certification_status() == silex::CertificationMode::unknown);
    assert(check_unit_proof_record(proved, 0, 2,
                                   silex::ProofState::verified, 31, -1,
                                   false));

    silex::Element epsilon2(field);
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement square_generators[] = {std::move(square)};
    silex::OrderUnitGroup square_group(order);
    assert(square_group.set_units(order,
                                  silex::FactoredElementSpan(square_generators,
                                                             1),
                                  embeddings, 128));

    silex::OrderUnitGroup mutated(order);
    status = silex::ProofState::not_checked;
    changed = false;
    assert(mutated.prove_local_saturated(
            status, changed, square_group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::unavailable);
    assert(changed);
    assert(check_first_free_generator_square(mutated, epsilon2));
    assert(check_unit_proof_record(mutated, 0, 2,
                                   silex::ProofState::unavailable, 31, -1,
                                   true));

    assert(set_fmpz_si(bound, 2));
    silex::OrderUnitGroup unavailable(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(unavailable.prove_local_saturated(
            status, changed, square_group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::unavailable);
    assert(!changed);
    assert(check_first_free_generator(unavailable, epsilon2));
    assert(check_unit_proof_record(unavailable, 0, 2,
                                   silex::ProofState::unavailable, 2, 0,
                                   false));

    sflint::Fmpz not_prime;
    assert(set_fmpz_si(not_prime, 4));
    status = silex::ProofState::verified;
    changed = true;
    assert(!unavailable.prove_local_saturated(
            status, changed, square_group, sflint::FmpzConstRef(not_prime), 1,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(changed);

    return 0;
}

int test_regulator_index_bound() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz bound;
    assert(rank_zero.regulator_index_bound(sflint::FmpzRef(bound), 128));
    assert(sflint::fmpz_equal_si(bound, 1));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup full;
    assert(full.compute(order));
    assert(full.regulator_index_bound(sflint::FmpzRef(bound), 160));
    assert(sflint::fmpz_sgn(sflint::FmpzConstRef(bound)) > 0);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(square)};
    silex::OrderUnitGroup subgroup(order);
    assert(subgroup.set_units(order, silex::FactoredElementSpan(generators, 1),
                              embeddings, 160));
    sflint::Fmpz subgroup_bound;
    assert(subgroup.regulator_index_bound(sflint::FmpzRef(subgroup_bound),
                                          160));
    assert(fmpz_cmp(subgroup_bound.raw(), bound.raw()) >= 0);

    sflint::Fmpz sentinel;
    assert(set_fmpz_si(sentinel, 42));
    silex::OrderUnitGroup unset(order);
    assert(!unset.regulator_index_bound(sflint::FmpzRef(sentinel), 160));
    assert(sflint::fmpz_equal_si(sentinel, 42));
    assert(!rank_zero.regulator_index_bound(sflint::FmpzRef(sentinel), 0));
    assert(sflint::fmpz_equal_si(sentinel, 42));

    return 0;
}

int test_class_regulator_index_bound() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderUnitGroup units;
    assert(units.compute(order));

    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        context = silex::ClassGroupContext(order);
        assert(context.is_defined());
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

        const silex::FactorBase* base = context.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element alpha(field);
        assert(alpha.set_si(-2));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        assert(context.has_presentation());
    };

    silex::ClassGroupContext class_group;
    prepare_context(class_group);

    sflint::Arb product;
    assert(units.class_regulator_product(sflint::ArbRef(product),
                                         class_group, 128));
    assert(sflint::arb_contains_si(product, 1));

    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb analytic_hR;
    sflint::arb_set_fmpq(analytic_hR, half, 128);

    sflint::Fmpz index_bound;
    assert(units.class_regulator_index_bound(
            sflint::FmpzRef(index_bound), class_group,
            sflint::ArbConstRef(analytic_hR), 128));
    assert(sflint::fmpz_equal_si(index_bound, 2));

    sflint::Fmpz sentinel;
    assert(set_fmpz_si(sentinel, 42));
    sflint::Arb zero;
    sflint::arb_zero(zero);
    assert(!units.class_regulator_index_bound(
            sflint::FmpzRef(sentinel), class_group,
            sflint::ArbConstRef(zero), 128));
    assert(sflint::fmpz_equal_si(sentinel, 42));
    assert(!units.class_regulator_index_bound(
            sflint::FmpzRef(sentinel), class_group,
            sflint::ArbConstRef(analytic_hR), 0));
    assert(sflint::fmpz_equal_si(sentinel, 42));

    return 0;
}

int test_class_regulator_index_bound_interval_boundary() {
    sflint::Arb candidate_product;
    sflint::Arb analytic_product;
    sflint::Arb analytic_error;
    sflint::arb_set_si(candidate_product, 2);
    sflint::arb_set_si(analytic_product, 2);
    sflint::arb_set_si(analytic_error, 1);
    sflint::arb_add_error(analytic_product, analytic_error);
    assert(sflint::arb_contains_si(analytic_product, 1));
    assert(sflint::arb_contains_si(analytic_product, 3));

    sflint::Fmpz bound;
    assert(silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(analytic_product), 128, nullptr));
    assert(sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) >= 0);

    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb narrow_analytic_product;
    sflint::Arb narrow_error;
    sflint::arb_set_si(narrow_analytic_product, 2);
    sflint::arb_set_fmpq(narrow_error, half, 128);
    sflint::arb_add_error(narrow_analytic_product, narrow_error);
    assert(silex::detail::
                   class_regulator_index_is_one_from_candidate_product(
                           sflint::ArbConstRef(candidate_product),
                           sflint::ArbConstRef(narrow_analytic_product), 128,
                           nullptr));
    assert(silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(narrow_analytic_product), 128, nullptr));
    assert(sflint::fmpz_equal_si(bound, 2));

    return 0;
}

int test_class_unit_regulator_certification() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;

    silex::ClassGroupContext class_group;
    assert(class_group.compute_candidate(order, sflint::FmpzConstRef(bound),
                                         options));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::unknown);

    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup units(order);
    assert(units.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 256));
    assert(units.certification_status() == silex::CertificationMode::unknown);

    sflint::Arb analytic_hR;
    assert(units.class_regulator_product(sflint::ArbRef(analytic_hR),
                                         class_group, 256));
    assert(class_group.try_certify_class_unit_with_units(
            units, sflint::ArbConstRef(analytic_hR), 256));
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.analytic_class_regulator_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(units.certification_status() == silex::CertificationMode::proven);

    return 0;
}

int test_prove_index_bound() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 5));
    silex::OrderUnitGroup proved_rank_zero(degree_one_order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved_rank_zero.prove_index_bound(
            status, changed, rank_zero, 1, sflint::FmpzConstRef(aux_bound), 1,
            degree_one_embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved_rank_zero.certification_status() ==
           silex::CertificationMode::proven);
    assert(proved_rank_zero.unit_proof_record_count() == 0);

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup full_rank(order);
    assert(full_rank.set_units(order,
                               silex::FactoredElementSpan(
                                       primitive_generators, 1),
                               embeddings, 128));

    assert(set_fmpz_si(aux_bound, 31));
    silex::OrderUnitGroup proved(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_index_bound(
            status, changed, full_rank, 1, sflint::FmpzConstRef(aux_bound), 1,
            embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    assert(check_first_free_generator(proved, epsilon));
    assert(has_unit_proof_record(proved, 2, silex::ProofState::verified,
                                 false));

    silex::Element epsilon2(field);
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement square_generators[] = {std::move(square)};
    silex::OrderUnitGroup square_group(order);
    assert(square_group.set_units(order,
                                  silex::FactoredElementSpan(square_generators,
                                                             1),
                                  embeddings, 128));

    silex::OrderUnitGroup restarted(order);
    status = silex::ProofState::not_checked;
    changed = false;
    assert(restarted.prove_index_bound(
            status, changed, square_group, 1, sflint::FmpzConstRef(aux_bound),
            1, embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(changed);
    assert(restarted.certification_status() ==
           silex::CertificationMode::proven);
    assert(check_first_free_generator_square(restarted, epsilon2));
    assert(has_unit_proof_record(restarted, 2, silex::ProofState::verified,
                                 false));

    silex::OrderUnitGroup restart_limited(order);
    status = silex::ProofState::not_checked;
    changed = false;
    assert(restart_limited.prove_index_bound(
            status, changed, square_group, 1, sflint::FmpzConstRef(aux_bound),
            0, embeddings, 128));
    assert(status == silex::ProofState::unavailable);
    assert(changed);
    assert(restart_limited.certification_status() ==
           silex::CertificationMode::unknown);
    assert(check_first_free_generator_square(restart_limited, epsilon2));
    assert(has_unit_proof_record(restart_limited, 2,
                                 silex::ProofState::unavailable, true));

    sflint::Fmpz bad_bound;
    assert(set_fmpz_si(bad_bound, 1));
    status = silex::ProofState::verified;
    changed = true;
    assert(!proved.prove_index_bound(
            status, changed, full_rank, 1, sflint::FmpzConstRef(bad_bound), 1,
            embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(changed);
    assert(check_first_free_generator(proved, epsilon));

    return 0;
}

int test_saturate_index_bounded() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 5));
    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    bool stable = false;
    assert(copied.saturate_index_bounded(
            changed, stable, rank_zero, degree_one_embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 2, 128));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(square)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::OrderUnitGroup saturated(order);
    changed = false;
    stable = true;
    assert(saturated.saturate_index_bounded(
            changed, stable, group, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator_square(saturated, epsilon2));
    assert(saturated.certification_status() ==
           silex::CertificationMode::unknown);

    silex::OrderUnitGroup zero_pass(order);
    changed = true;
    stable = true;
    assert(zero_pass.saturate_index_bounded(
            changed, stable, group, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 0, 128));
    assert(!changed);
    assert(!stable);
    assert(check_first_free_generator(zero_pass, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup no_progress(order);
    assert(no_progress.set_units(order,
                                 silex::FactoredElementSpan(
                                         primitive_generators, 1),
                                 embeddings, 128));
    silex::OrderUnitGroup preserved(order);
    assert(preserved.set(no_progress));
    assert(set_fmpz_si(aux_bound, 2));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator(preserved, epsilon));

    sflint::Fmpz bad_bound;
    assert(set_fmpz_si(bad_bound, 1));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(bad_bound), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator(preserved, epsilon));

    return 0;
}

int test_saturate_index_bounded_adaptive() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz start;
    sflint::Fmpz max;
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 8));
    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    bool stable = false;
    assert(copied.saturate_index_bounded_adaptive(
            changed, stable, rank_zero, degree_one_embeddings, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2, 128));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(square)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));
    silex::OrderUnitGroup adaptive(order);
    changed = false;
    stable = true;
    assert(adaptive.saturate_index_bounded_adaptive(
            changed, stable, group, embeddings, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator_square(adaptive, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup no_progress(order);
    assert(no_progress.set_units(order,
                                 silex::FactoredElementSpan(
                                         primitive_generators, 1),
                                 embeddings, 128));
    silex::OrderUnitGroup preserved(order);
    assert(preserved.set(no_progress));
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 2));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded_adaptive(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator(preserved, epsilon));

    sflint::Fmpz bad_start;
    assert(set_fmpz_si(bad_start, 1));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded_adaptive(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(bad_start), sflint::FmpzConstRef(max), 3,
            128));
    assert(changed);
    assert(!stable);

    return 0;
}

int test_saturate_local_once_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 5));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    assert(saturated.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(changed);
    assert(saturated.free_rank() == 1);
    assert(check_first_free_generator_square(saturated, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup already_saturated(order);
    assert(already_saturated.set_units(
            order, silex::FactoredElementSpan(primitive_generators, 1),
            embeddings, 128));

    silex::OrderUnitGroup copied(order);
    changed = true;
    assert(copied.saturate_local_once(
            changed, already_saturated, prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(!changed);
    assert(check_first_free_generator(copied, epsilon));

    return 0;
}

int test_saturate_local_once_rank_zero_and_failures() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 5));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, degree_one_order,
                                  sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    sflint::FmpzMat kernel(5, 5);
    assert(rank_zero.residue_dlog_kernel(kernel, prime_span,
                                         sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 0);

    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    assert(copied.saturate_local_once(
            changed, rank_zero, prime_span, sflint::FmpzConstRef(ell),
            degree_one_embeddings, 80));
    assert(!changed);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::NumberField other_field = quadratic_field(3);
    silex::Order order;
    silex::Order other_order;
    order = silex::test::equation_order(field);
    other_order = silex::test::equation_order(other_field);
    silex::EmbeddingContext embeddings(field);
    silex::EmbeddingContext other_embeddings(other_field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::PrimeIdealList other_primes;
    assert(silex::decompose_prime(other_primes, other_order,
                                  sflint::FmpzConstRef(p)));
    assert(other_primes.size() >= 1);
    silex::PrimeIdealSpan other_prime_span(other_primes.at(0), 1);

    sflint::Fmpz not_prime;
    assert(set_fmpz_si(not_prime, 4));
    assert(!copied.saturate_local_once(
            changed, group, silex::PrimeIdealSpan(),
            sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!copied.saturate_local_once(
            changed, group, other_prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(!copied.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(not_prime),
            embeddings, 128));
    assert(!copied.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(ell),
            other_embeddings, 128));
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    return 0;
}

int test_select_saturation_primes() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 5));

    silex::PrimeIdealList selected;
    assert(group.select_saturation_primes(
            selected, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound)));
    assert(selected.size() == 1);
    assert(selected.at(0) != nullptr);
    assert(silex::same_order_parent(selected.at(0)->parent(), &order));

    assert(group.select_saturation_primes(
            selected, sflint::FmpzConstRef(ell), 0,
            sflint::FmpzConstRef(bound)));
    assert(selected.size() == 0);

    silex::PrimeIdealList preserved;
    assert(group.select_saturation_primes(
            preserved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound)));
    assert(set_fmpz_si(bound, 2));
    assert(!group.select_saturation_primes(
            preserved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound)));
    assert(preserved.size() == 1);

    return 0;
}

int test_saturate_bounded_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 5));

    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    bool stable = false;
    assert(saturated.saturate_bounded(
            changed, stable, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, embeddings, 128));
    assert(changed);
    assert(stable);
    assert(saturated.free_rank() == 1);
    assert(check_first_free_generator_square(saturated, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup already(order);
    assert(already.set_units(order,
                             silex::FactoredElementSpan(primitive_generators,
                                                        1),
                             embeddings, 128));

    silex::OrderUnitGroup unchanged(order);
    changed = true;
    stable = false;
    assert(unchanged.saturate_bounded(
            changed, stable, already, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, embeddings, 128));
    assert(!changed);
    assert(stable);
    assert(check_first_free_generator(unchanged, epsilon));

    return 0;
}

int test_saturate_bounded_rank_zero_zero_pass_partial_and_failures() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 5));

    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    bool stable = false;
    assert(copied.saturate_bounded(
            changed, stable, rank_zero, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, degree_one_embeddings, 80));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    changed = true;
    stable = true;
    assert(copied.saturate_bounded(
            changed, stable, rank_zero, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 0, degree_one_embeddings, 80));
    assert(!changed);
    assert(!stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    silex::Element epsilon4(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    assert(epsilon4.multiply(epsilon2, epsilon2));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon4));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::OrderUnitGroup partial(order);
    changed = false;
    stable = true;
    assert(partial.saturate_bounded(
            changed, stable, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 1, embeddings, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator_square(partial, epsilon4));

    silex::OrderUnitGroup preserved(degree_one_order);
    assert(preserved.set(rank_zero));
    assert(set_fmpz_si(bound, 2));
    changed = false;
    stable = false;
    assert(!preserved.saturate_bounded(
            changed, stable, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, embeddings, 128));
    assert(check_rank_zero_group(preserved, degree_one_order, 2));

    return 0;
}

int test_supplied_real_quadratic_units() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element unit(field);
    assert(set_real_quadratic_unit(unit));
    silex::FactoredElement generator(field);
    assert(generator.set_element(unit));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));
    assert(group.is_defined());
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    assert(group.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_positive(regulator));
    assert(check_first_free_generator(group, unit));

    silex::OrderUnitGroup copy;
    assert(copy.set(group));
    assert(copy.free_rank() == 1);
    assert(copy.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_positive(regulator));
    assert(check_first_free_generator(copy, unit));

    silex::OrderUnitGroup moved(std::move(copy));
    assert(moved.free_rank() == 1);
    assert(check_first_free_generator(moved, unit));
    assert(!copy.is_defined());

    return 0;
}

int test_supplied_rank_zero_units() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(), embeddings,
                           80));
    assert(group.is_defined());
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &order));
    assert(group.free_rank() == 0);
    assert(group.certification_status() == silex::CertificationMode::unknown);

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    assert(group.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_one(regulator));

    return 0;
}

int test_supplied_unit_failures_preserve_output() {
    silex::NumberField field = quadratic_field(2);
    silex::NumberField other_field = quadratic_field(3);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::EmbeddingContext other_embeddings(other_field);

    silex::Element unit(field);
    assert(set_real_quadratic_unit(unit));
    silex::FactoredElement compact_unit(field);
    assert(compact_unit.set_element(unit));
    silex::FactoredElement units[] = {std::move(compact_unit)};

    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(units, 1),
                           embeddings, 128));
    sflint::Arb regulator;
    assert(group.regulator(sflint::ArbRef(regulator)));

    silex::Element two(field);
    assert(two.set_si(2));
    silex::FactoredElement nonunit(field);
    assert(nonunit.set_element(two));
    silex::FactoredElement nonunits[] = {std::move(nonunit)};
    assert(!group.set_units(order, silex::FactoredElementSpan(nonunits, 1),
                            embeddings, 128));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator(group, unit));

    silex::Element minus_one(field);
    assert(minus_one.set_si(-1));
    silex::FactoredElement torsion(field);
    assert(torsion.set_element(minus_one));
    silex::FactoredElement torsion_generators[] = {std::move(torsion)};
    assert(!group.set_units(
            order, silex::FactoredElementSpan(torsion_generators, 1),
            embeddings, 128));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator(group, unit));

    assert(!group.set_units(order, silex::FactoredElementSpan(), embeddings,
                            128));
    assert(!group.set_units(order, silex::FactoredElementSpan(units, 1),
                            embeddings, 0));
    assert(!group.set_units(order, silex::FactoredElementSpan(units, 1),
                            other_embeddings, 128));

    silex::Element foreign_unit(other_field);
    assert(set_real_quadratic_unit(foreign_unit));
    silex::FactoredElement foreign(other_field);
    assert(foreign.set_element(foreign_unit));
    silex::FactoredElement foreign_generators[] = {std::move(foreign)};
    assert(!group.set_units(
            order, silex::FactoredElementSpan(foreign_generators, 1),
            embeddings, 128));
    assert(group.is_set());
    assert(group.free_rank() == 1);
    assert(check_first_free_generator(group, unit));

    return 0;
}

int test_imaginary_quadratic_torsion() {
    silex::NumberField qi = quadratic_field(-1);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(qi);

    silex::OrderUnitGroup group;
    assert(group.compute(equation_order));
    assert(check_rank_zero_group(group, equation_order, 4));

    return 0;
}

int test_nonmaximal_torsion_filtering() {
    silex::NumberField qsqrt_minus3 = quadratic_field(-3);
    silex::Order equation_order;
    silex::Order maximal_order(qsqrt_minus3);
    equation_order = silex::test::equation_order(qsqrt_minus3);
    assert(maximal_order.maximal_order(equation_order));

    silex::OrderUnitGroup nonmaximal;
    assert(nonmaximal.compute(equation_order));
    assert(check_rank_zero_group(nonmaximal, equation_order, 2));

    silex::OrderUnitGroup maximal;
    assert(maximal.compute(maximal_order));
    assert(check_rank_zero_group(maximal, maximal_order, 6));

    return 0;
}

int test_positive_rank_failure_preserves_output() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);

    silex::OrderUnitGroup group;
    assert(group.compute(degree_one_order));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    silex::NumberField real_quadratic = quadratic_field(5);
    silex::Order real_order;
    real_order = silex::test::equation_order(real_quadratic);
    assert(!real_order.is_maximal());
    assert(!group.compute(real_order));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    return 0;
}

int test_group_keeps_parent_order_alive() {
    silex::OrderUnitGroup group = real_quadratic_unit_group();
    assert(group.is_defined());
    assert(group.is_set());
    assert(group.parent() != nullptr);
    assert(group.parent()->parent() != nullptr);
    assert(group.free_rank() == 1);

    silex::OrderElement torsion(*group.parent());
    assert(group.torsion_generator(torsion));

    silex::FactoredElement generator(*group.parent()->parent());
    silex::Element expanded(*group.parent()->parent());
    silex::OrderElement order_element(*group.parent());
    silex::Ideal principal(*group.parent());
    assert(group.free_generator(generator, 0));
    assert(generator.evaluate(expanded));
    assert(order_element.set_element(expanded));
    assert(principal.set_principal(order_element));
    assert(principal.is_one());

    return 0;
}

int test_set_move_and_access_failures() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderUnitGroup group;
    assert(group.compute(order));

    silex::OrderUnitGroup copy;
    assert(copy.set(group));
    assert(check_rank_zero_group(copy, order, 2));

    silex::OrderUnitGroup moved(std::move(copy));
    assert(check_rank_zero_group(moved, order, 2));
    assert(!copy.is_defined());

    sflint::Fmpz order_out;
    silex::NumberField other_field = quadratic_field(-1);
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);
    silex::OrderElement wrong_generator(other_order);
    assert(!moved.torsion_generator(wrong_generator));
    assert(moved.torsion_order(sflint::FmpzRef(order_out)));

    silex::FactoredElement free_generator(field);
    assert(!moved.free_generator(free_generator, 0));

    moved.clear();
    assert(!moved.is_defined());
    assert(moved.define(order));
    assert(moved.is_defined());
    assert(!moved.is_set());

    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Order sqrt2_order;
    sqrt2_order = silex::test::equation_order(sqrt2);
    silex::Element epsilon(sqrt2);
    assert(set_quadratic_coeffs(epsilon, 1, 1, 1, 1));

    silex::OrderUnitGroup positive_rank;
    assert(positive_rank.compute(sqrt2_order));
    assert(check_real_quadratic_group(positive_rank, sqrt2_order, epsilon));

    silex::OrderUnitGroup positive_rank_copy;
    assert(positive_rank_copy.set(positive_rank));
    assert(check_real_quadratic_group(positive_rank_copy, sqrt2_order,
                                      epsilon));

    silex::OrderUnitGroup positive_rank_move;
    positive_rank_move = std::move(positive_rank_copy);
    assert(check_real_quadratic_group(positive_rank_move, sqrt2_order,
                                      epsilon));
    assert(!positive_rank_copy.is_defined());

    return 0;
}

}  // namespace

int main() {
    test_unit_extraction_cache_clear();
    test_degree_one_rank_zero_group();
    test_compute_real_quadratic_proven();
    test_compute_with_class_group_rank_zero();
    test_compute_with_class_group_real_quadratic();
    test_compute_with_class_group_real_quadratic_210();
    test_nonmaximal_quadratic_pair_rejection_preserves_output();
    test_compute_with_class_group_quintic_proven();
    test_compute_with_class_group_quintic_small_bf_budget_trivial();
    test_compute_with_class_group_failures_preserve_output();
    test_relation_kernel_rescale_log_matrix_word_rounding();
    test_reduced_regulator_reconstruction_correction();
    test_hnf_regulator_multiple_certification_correction();
    test_regulator_pivot_exponent_policy();
    test_regulator_pivot_arb_ambiguity_and_invalidity();
    test_regulator_pivot_clean_cols_and_mapping();
    test_regulator_pivot_multi_column_elimination();
    test_regulator_pivot_finish_route_matrix();
    test_relation_kernel_units_rank_one();
    test_relation_kernel_units_bounded_rank_one();
    test_relation_kernel_units_index_bounded_rank_one();
    test_relation_kernel_units_index_bounded_quartic_power_root();
    test_relation_kernel_units_index_bounded_cubic_rank_target();
    test_relation_kernel_units_index_bounded_saturated();
    test_relation_kernel_units_rank_zero_delegates();
    test_relation_kernel_unit_failures_preserve_output();
    test_saturate_row_real_quadratic_square();
    test_saturate_row_no_root_and_divisible_copy();
    test_saturate_row_rank_zero_and_failures();
    test_residue_dlog_kernel_real_quadratic();
    test_residue_dlog_proof_kernel_torsion();
    test_rank_one_degree_one_root_image_matches_prime_ideal();
    test_direct_degree_one_proof_column_matches_residue_field();
    test_direct_degree_one_residue_dlog_matrix_matches_residue_field();
    test_half_integral_direct_degree_one_residue_dlog_matches_residue_field();
    test_select_saturation_proof_primes();
    test_unit_proof_records();
    test_prove_local_saturated_real_quadratic();
    test_regulator_index_bound();
    test_class_regulator_index_bound();
    test_class_regulator_index_bound_interval_boundary();
    test_class_unit_regulator_certification();
    test_prove_index_bound();
    test_saturate_index_bounded();
    test_saturate_index_bounded_adaptive();
    test_saturate_local_once_real_quadratic();
    test_saturate_local_once_rank_zero_and_failures();
    test_select_saturation_primes();
    test_saturate_bounded_real_quadratic();
    test_saturate_bounded_rank_zero_zero_pass_partial_and_failures();
    test_supplied_real_quadratic_units();
    test_supplied_rank_zero_units();
    test_supplied_unit_failures_preserve_output();
    test_imaginary_quadratic_torsion();
    test_nonmaximal_torsion_filtering();
    test_positive_rank_failure_preserves_output();
    test_group_keeps_parent_order_alive();
    test_set_move_and_access_failures();
    return 0;
}
