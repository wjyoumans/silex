#include <silex/zeta.hpp>

#include <silex/embedding.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/unit.hpp>

#include "test_support.hpp"

#include <cassert>

namespace {
namespace sflint = silex::flint;

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

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quintic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -6);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, -7);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

bool radius_lt_2exp_si(const sflint::Arb& value, slong exponent) noexcept {
    sflint::Arb radius;
    sflint::Arb target;
    sflint::arb_get_rad_arb(radius, value);
    sflint::arb_one(target);
    sflint::arb_mul_2exp_si(target, target, exponent);
    return sflint::arb_lt(radius, target);
}

int test_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    sflint::Arb hR;
    sflint::Arb residue;
    sflint::Arb log_residue;
    assert(silex::zeta_class_regulator_product(sflint::ArbRef(hR),
                                               order, 128));
    assert(sflint::arb_contains_si(hR, 1));
    assert(silex::zeta_residue(sflint::ArbRef(residue), order, 128));
    assert(sflint::arb_contains_si(residue, 1));
    assert(silex::zeta_log_residue(sflint::ArbRef(log_residue),
                                   order, 128));
    assert(sflint::arb_contains_zero(log_residue));

    ulong cutoff = 1;
    slong work_precision = 0;
    sflint::Arb error_bound;
    assert(silex::zeta_class_regulator_product_bf_audit(
            sflint::ArbRef(hR), sflint::ArbRef(error_bound), cutoff,
            work_precision, order, 20000, 128));
    auto hR_audit = silex::zeta_class_regulator_product_bf_audit(
            order, 20000, 128);
    assert(hR_audit.has_value());
    assert(sflint::arb_contains_si(hR, 1));
    assert(sflint::arb_is_zero(error_bound));
    assert(cutoff == 0);
    assert(work_precision == 128);
    assert(sflint::arb_contains_si(hR_audit->value, 1));
    assert(sflint::arb_is_zero(hR_audit->error_bound));
    assert(hR_audit->cutoff == 0);
    assert(hR_audit->work_precision == 128);
    assert(!silex::zeta_log_residue_bf(sflint::ArbRef(log_residue),
                                       order, 20000, 128));

    return 0;
}

int test_imaginary_quadratic() {
    silex::NumberField field = quadratic_field(-47);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Arb hR;
    sflint::Arb residue;
    sflint::Arb log_residue;
    assert(silex::zeta_residue(sflint::ArbRef(residue), maximal, 192));
    assert(sflint::arb_is_positive(residue));
    assert(silex::zeta_log_residue(sflint::ArbRef(log_residue),
                                   maximal, 192));
    assert(sflint::arb_is_finite(log_residue));
    assert(silex::zeta_class_regulator_product(sflint::ArbRef(hR),
                                               maximal, 192));
    assert(sflint::arb_contains_si(hR, 5));

    assert(!silex::zeta_class_regulator_product(sflint::ArbRef(hR),
                                                equation, 192));
    assert(!silex::zeta_residue(sflint::ArbRef(residue), equation, 192));
    assert(!silex::zeta_log_residue(sflint::ArbRef(log_residue),
                                    equation, 192));

    return 0;
}

int test_real_quadratic_regulator_overlap() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::Element epsilon(field);
    assert(silex::quadratic_fundamental_unit(epsilon, field));

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.refine(256));
    silex::Element units[] = {std::move(epsilon)};
    sflint::Arb regulator;
    assert(silex::unit_regulator(sflint::ArbRef(regulator), embeddings,
                                  silex::ElementSpan(units, 1), 256));

    sflint::Arb hR;
    assert(silex::zeta_class_regulator_product(sflint::ArbRef(hR),
                                               maximal, 256));
    assert(sflint::arb_overlaps(hR, regulator));

    return 0;
}

int test_cubic_bf() {
    silex::NumberField field = cubic_field();
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Arb hR;
    sflint::Arb residue;
    sflint::Arb log_residue;
    sflint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;

    assert(silex::zeta_log_residue(sflint::ArbRef(log_residue),
                                   maximal, 128));
    assert(sflint::arb_is_finite(log_residue));
    assert(radius_lt_2exp_si(log_residue, -2));
    assert(silex::zeta_residue(sflint::ArbRef(residue), maximal, 128));
    assert(sflint::arb_is_positive(residue));
    assert(silex::zeta_class_regulator_product(sflint::ArbRef(hR),
                                               maximal, 128));
    assert(sflint::arb_is_positive(hR));

    assert(!silex::zeta_log_residue_bf(sflint::ArbRef(log_residue),
                                       maximal, 72, 128));
    assert(silex::zeta_log_residue_bf(sflint::ArbRef(log_residue),
                                      maximal, 20000, 128));
    assert(sflint::arb_is_finite(log_residue));
    assert(radius_lt_2exp_si(log_residue, -2));
    auto log_audit = silex::zeta_log_residue_bf_audit(maximal, 20000, 128);
    assert(log_audit.has_value());
    assert(sflint::arb_is_finite(log_audit->value));
    assert(sflint::arb_is_positive(log_audit->error_bound));
    assert(log_audit->cutoff >= 70);
    assert(log_audit->cutoff <= 20007);
    assert(log_audit->work_precision >= 192);

    assert(silex::zeta_residue_bf_audit(
            sflint::ArbRef(residue), sflint::ArbRef(error_bound),
            cutoff, work_precision, maximal, 20000, 128));
    auto residue_audit = silex::zeta_residue_bf_audit(maximal, 20000, 128);
    assert(residue_audit.has_value());
    assert(sflint::arb_is_positive(residue));
    assert(sflint::arb_is_positive(error_bound));
    assert(cutoff >= 70);
    assert(cutoff <= 20007);
    assert(work_precision >= 192);
    assert(sflint::arb_is_positive(residue_audit->value));
    assert(sflint::arb_is_positive(residue_audit->error_bound));
    assert(residue_audit->cutoff >= 70);
    assert(residue_audit->cutoff <= 20007);
    assert(residue_audit->work_precision >= 192);

    cutoff = 0;
    work_precision = 0;
    assert(silex::zeta_class_regulator_product_bf_audit(
            sflint::ArbRef(hR), sflint::ArbRef(error_bound), cutoff,
            work_precision, maximal, 20000, 128));
    auto product_audit = silex::zeta_class_regulator_product_bf_audit(
            maximal, 20000, 128);
    assert(product_audit.has_value());
    assert(sflint::arb_is_positive(hR));
    assert(sflint::arb_is_positive(error_bound));
    assert(cutoff >= 70);
    assert(cutoff <= 20007);
    assert(work_precision >= 192);
    assert(sflint::arb_is_positive(product_audit->value));
    assert(sflint::arb_is_positive(product_audit->error_bound));
    assert(product_audit->cutoff >= 70);
    assert(product_audit->cutoff <= 20007);
    assert(product_audit->work_precision >= 192);

    return 0;
}

int test_quintic_bf() {
    silex::NumberField field = quintic_field();
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Arb log_residue;
    assert(silex::zeta_log_residue(sflint::ArbRef(log_residue),
                                   maximal, 128));
    assert(sflint::arb_is_finite(log_residue));

    return 0;
}

}  // namespace

int main() {
    test_degree_one();
    test_imaginary_quadratic();
    test_real_quadratic_regulator_overlap();
    test_cubic_bf();
    test_quintic_bf();
    return 0;
}
