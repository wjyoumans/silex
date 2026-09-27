#include <silex/zeta.hpp>

#include <silex/embedding.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/unit.hpp>

#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <initializer_list>

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

silex::NumberField field_from_coefficients(const slong* coefficients,
                                           slong length) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i < length; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::Order maximal_order_of(const silex::NumberField& field) noexcept {
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    return maximal;
}

bool arb_identical(const sflint::Arb& left, const sflint::Arb& right) noexcept {
    return ::arb_equal(left.raw(), right.raw()) != 0;
}

// Ground-truth hR values come from an external reference, GP 2.17.4, not
// from Silex:
//
//     default(realprecision, 38);
//     b = bnfinit(P, 1); bnfcertify(b)   \\ returns 1 for every P below
//     b.no * b.reg
//
// The decimal strings are the printed 38-digit values; the ball radius
// 1e-30 covers their rounding.  All fields have degree >= 3, so Silex uses
// the Belabas-Friedman route, whose error bound assumes GRH.
struct HrGroundTruth {
    slong coefficients[8];  // constant term first
    slong length;
    const char* hR;
};

bool ground_truth_ball(sflint::Arb& out, const char* decimal) noexcept {
    char buffer[96];
    std::size_t length = 0;
    while (decimal[length] != '\0' && length + 16 < sizeof(buffer)) {
        buffer[length] = decimal[length];
        ++length;
    }
    if (decimal[length] != '\0') {
        return false;
    }
    const char suffix[] = " +/- 1e-30";
    for (std::size_t i = 0; i < sizeof(suffix); ++i) {
        buffer[length + i] = suffix[i];
    }
    return ::arb_set_str(out.raw(), buffer, 256) == 0;
}

void set_audit_sentinels(sflint::Arb& value,
                         sflint::Arb& error_bound,
                         ulong& cutoff,
                         slong& work_precision) noexcept {
    sflint::arb_set_si(value, 12345);
    sflint::arb_set_si(error_bound, -777);
    cutoff = 424242;
    work_precision = 31337;
}

bool audit_sentinels_unchanged(const sflint::Arb& value,
                               const sflint::Arb& error_bound,
                               ulong cutoff,
                               slong work_precision) noexcept {
    return ::arb_equal_si(value.raw(), 12345) != 0 &&
           ::arb_equal_si(error_bound.raw(), -777) != 0 &&
           cutoff == 424242 && work_precision == 31337;
}

// Calls every audit API in its output-buffer form, expects failure, and
// checks that none of the caller's outputs changed; also checks that the
// owned-return forms report failure.
void assert_audits_fail_atomically(const silex::Order& order,
                                   ulong max_cutoff,
                                   slong precision,
                                   bool check_product) noexcept {
    sflint::Arb value;
    sflint::Arb error_bound;
    ulong cutoff = 0;
    slong work_precision = 0;

    set_audit_sentinels(value, error_bound, cutoff, work_precision);
    assert(!silex::zeta_log_residue_bf_audit(
            sflint::ArbRef(value), sflint::ArbRef(error_bound), cutoff,
            work_precision, order, max_cutoff, precision));
    assert(audit_sentinels_unchanged(value, error_bound, cutoff,
                                     work_precision));
    assert(!silex::zeta_log_residue_bf_audit(order, max_cutoff, precision)
                    .has_value());

    set_audit_sentinels(value, error_bound, cutoff, work_precision);
    assert(!silex::zeta_residue_bf_audit(
            sflint::ArbRef(value), sflint::ArbRef(error_bound), cutoff,
            work_precision, order, max_cutoff, precision));
    assert(audit_sentinels_unchanged(value, error_bound, cutoff,
                                     work_precision));
    assert(!silex::zeta_residue_bf_audit(order, max_cutoff, precision)
                    .has_value());

    if (check_product) {
        set_audit_sentinels(value, error_bound, cutoff, work_precision);
        assert(!silex::zeta_class_regulator_product_bf_audit(
                sflint::ArbRef(value), sflint::ArbRef(error_bound), cutoff,
                work_precision, order, max_cutoff, precision));
        assert(audit_sentinels_unchanged(value, error_bound, cutoff,
                                         work_precision));
        assert(!silex::zeta_class_regulator_product_bf_audit(
                        order, max_cutoff, precision)
                        .has_value());
    }

    sflint::arb_set_si(value, 12345);
    assert(!silex::zeta_log_residue_bf(sflint::ArbRef(value), order,
                                       max_cutoff, precision));
    assert(::arb_equal_si(value.raw(), 12345) != 0);
    assert(!silex::zeta_residue_bf(sflint::ArbRef(value), order,
                                   max_cutoff, precision));
    assert(::arb_equal_si(value.raw(), 12345) != 0);
    if (check_product) {
        assert(!silex::zeta_class_regulator_product_bf(
                sflint::ArbRef(value), order, max_cutoff, precision));
        assert(::arb_equal_si(value.raw(), 12345) != 0);
    }
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

int test_audit_failure_leaves_outputs_unchanged() {
    // Degree one: the log-residue and residue BF routes are undefined and
    // fail; the product audit succeeds there (tested in test_degree_one).
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order = silex::test::equation_order(degree_one);
    assert_audits_fail_atomically(degree_one_order, 20000, 128, false);

    // Non-maximal order: x^3 - x^2 - 2x - 8 has index 2 over Z[x].
    const slong dedekind_coefficients[] = {-8, -2, -1, 1};
    silex::NumberField dedekind =
            field_from_coefficients(dedekind_coefficients, 4);
    silex::Order dedekind_equation = silex::test::equation_order(dedekind);
    assert_audits_fail_atomically(dedekind_equation, 20000, 128, true);

    silex::NumberField cubic = cubic_field();
    silex::Order maximal = maximal_order_of(cubic);
    // The cutoff cap is too small to reach the BF target.
    assert_audits_fail_atomically(maximal, 72, 128, true);
    // A zero cutoff cap and a nonpositive precision are invalid.
    assert_audits_fail_atomically(maximal, 0, 128, true);
    assert_audits_fail_atomically(maximal, 20000, 0, true);

    return 0;
}

int test_max_cutoff_near_uword_max() {
    silex::NumberField field = cubic_field();
    silex::Order maximal = maximal_order_of(field);

    auto log_reference = silex::zeta_log_residue_bf_audit(maximal, 20000, 128);
    auto residue_reference = silex::zeta_residue_bf_audit(maximal, 20000, 128);
    auto product_reference =
            silex::zeta_class_regulator_product_bf_audit(maximal, 20000, 128);
    assert(log_reference.has_value());
    assert(residue_reference.has_value());
    assert(product_reference.has_value());

    // The BF target is reached far below either cap, so an unbounded cap
    // selects the same cutoff and value as the 20000 cap.  UWORD_MAX and
    // UWORD_MAX - 5 are not multiples of 9; rounding them up used to wrap.
    for (ulong max_cutoff : {UWORD_MAX, UWORD_MAX - 5}) {
        auto log_audit =
                silex::zeta_log_residue_bf_audit(maximal, max_cutoff, 128);
        assert(log_audit.has_value());
        assert(log_audit->cutoff == log_reference->cutoff);
        assert(log_audit->work_precision == log_reference->work_precision);
        assert(arb_identical(log_audit->value, log_reference->value));
        assert(arb_identical(log_audit->error_bound,
                             log_reference->error_bound));

        auto residue_audit =
                silex::zeta_residue_bf_audit(maximal, max_cutoff, 128);
        assert(residue_audit.has_value());
        assert(residue_audit->cutoff == residue_reference->cutoff);
        assert(residue_audit->work_precision ==
               residue_reference->work_precision);
        assert(arb_identical(residue_audit->value, residue_reference->value));
        assert(arb_identical(residue_audit->error_bound,
                             residue_reference->error_bound));

        auto product_audit = silex::zeta_class_regulator_product_bf_audit(
                maximal, max_cutoff, 128);
        assert(product_audit.has_value());
        assert(product_audit->cutoff == product_reference->cutoff);
        assert(product_audit->work_precision ==
               product_reference->work_precision);
        assert(arb_identical(product_audit->value, product_reference->value));
        assert(arb_identical(product_audit->error_bound,
                             product_reference->error_bound));

        sflint::Arb value;
        assert(silex::zeta_log_residue_bf(sflint::ArbRef(value), maximal,
                                          max_cutoff, 128));
        assert(arb_identical(value, log_reference->value));
        assert(silex::zeta_residue_bf(sflint::ArbRef(value), maximal,
                                      max_cutoff, 128));
        assert(arb_identical(value, residue_reference->value));
        assert(silex::zeta_class_regulator_product_bf(
                sflint::ArbRef(value), maximal, max_cutoff, 128));
        assert(arb_identical(value, product_reference->value));
    }

    return 0;
}

int test_aliased_out_error_bound_publishes_value_last() {
    // A caller that passes the same Arb for `out` and `error_bound` is
    // nonsensical but legal. Every BF audit buffer form must leave that
    // Arb holding the value on success, like `zeta_log_residue_bf_audit`
    // (which already writes the error bound, then `out`, into the caller's
    // buffers): `out` is always the last write.
    silex::NumberField field = cubic_field();
    silex::Order maximal = maximal_order_of(field);
    ulong cutoff = 0;
    slong work_precision = 0;

    sflint::Arb log_residue;
    sflint::Arb log_error_bound;
    assert(silex::zeta_log_residue_bf_audit(
            sflint::ArbRef(log_residue), sflint::ArbRef(log_error_bound),
            cutoff, work_precision, maximal, 20000, 128));
    // The value and the error bound must differ, or the test would pass
    // vacuously.
    assert(!arb_identical(log_residue, log_error_bound));

    sflint::Arb log_aliased;
    ulong aliased_cutoff = 0;
    slong aliased_work_precision = 0;
    assert(silex::zeta_log_residue_bf_audit(
            sflint::ArbRef(log_aliased), sflint::ArbRef(log_aliased),
            aliased_cutoff, aliased_work_precision, maximal, 20000, 128));
    assert(arb_identical(log_aliased, log_residue));

    sflint::Arb residue;
    sflint::Arb residue_error_bound;
    assert(silex::zeta_residue_bf_audit(
            sflint::ArbRef(residue), sflint::ArbRef(residue_error_bound),
            cutoff, work_precision, maximal, 20000, 128));
    assert(!arb_identical(residue, residue_error_bound));

    sflint::Arb residue_aliased;
    assert(silex::zeta_residue_bf_audit(
            sflint::ArbRef(residue_aliased), sflint::ArbRef(residue_aliased),
            aliased_cutoff, aliased_work_precision, maximal, 20000, 128));
    assert(arb_identical(residue_aliased, residue));

    sflint::Arb product;
    sflint::Arb product_error_bound;
    assert(silex::zeta_class_regulator_product_bf_audit(
            sflint::ArbRef(product), sflint::ArbRef(product_error_bound),
            cutoff, work_precision, maximal, 20000, 128));
    assert(!arb_identical(product, product_error_bound));

    sflint::Arb product_aliased;
    assert(silex::zeta_class_regulator_product_bf_audit(
            sflint::ArbRef(product_aliased), sflint::ArbRef(product_aliased),
            aliased_cutoff, aliased_work_precision, maximal, 20000, 128));
    assert(arb_identical(product_aliased, product));

    return 0;
}

int test_hR_ground_truth() {
    const HrGroundTruth cases[] = {
            // x^3 - x - 1, disc -23, h = 1
            {{-1, -1, 0, 1}, 4,
             "0.28119957432296184651205076406787829979"},
            // x^3 - x^2 - 2x - 8 (Dedekind), disc -503, index 2, h = 1
            {{-8, -2, -1, 1}, 4,
             "7.0273467933610955236852570948239393112"},
            // x^4 - 4x^2 + 2, totally real, disc 2048, h = 1
            {{2, 0, -4, 0, 1}, 5,
             "2.4417950066199157657221421142843130099"},
            // x^5 - 7x^4 - 6x^3 - 3x^2 - 3x + 3, disc -401370255, h = 1
            {{3, -3, -3, -6, -7, 1}, 6,
             "734.65403300200295476728211517178502268"},
            // x^3 - x^2 - 2x + 1, totally real cyclic, disc 49, h = 1
            {{1, -2, -1, 1}, 4,
             "0.52545468212257238833882604544832450954"},
            // x^4 + 1 (Phi_8), disc 256, w = 8, h = 1
            {{1, 0, 0, 0, 1}, 5,
             "1.7627471740390860504652186499595846181"},
            // x^4 - x^2 + 1 (Phi_12), disc 144, w = 12, h = 1
            {{1, 0, -1, 0, 1}, 5,
             "1.3169578969248167086250463473079684440"},
            // x^6 - x^3 + 1 (Phi_9), disc -19683, w = 18, h = 1.  Before
            // T-038 Silex used w = 6 here and returned hR / 3.
            {{1, 0, 0, -1, 0, 0, 1}, 7,
             "3.3971498025847701145790480197497869449"},
            // Phi_7, disc -16807, w = 14, h = 1.  Before T-038 the w search
            // failed here.
            {{1, 1, 1, 1, 1, 1, 1}, 7,
             "2.1018187284902895533553041817932980382"},
    };

    for (const HrGroundTruth& entry : cases) {
        silex::NumberField field = field_from_coefficients(
                entry.coefficients, entry.length);
        silex::Order maximal = maximal_order_of(field);

        sflint::Arb truth;
        assert(ground_truth_ball(truth, entry.hR));

        // Degree >= 3 uses the Belabas-Friedman route, whose error bound
        // assumes GRH; containment of the certified reference value checks the
        // BF sum, the error bound, and the residue-to-hR conversion.
        sflint::Arb hR;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(hR),
                                                   maximal, 128));
        assert(sflint::arb_contains(hR, truth));

        // The quintic needs a cutoff above 20000 to reach the BF target,
        // so the audit runs without a practical cap.
        auto audit = silex::zeta_class_regulator_product_bf_audit(
                maximal, UWORD_MAX, 128);
        assert(audit.has_value());
        assert(sflint::arb_contains(audit->value, truth));
    }

    return 0;
}

}  // namespace

int main() {
    test_degree_one();
    test_imaginary_quadratic();
    test_real_quadratic_regulator_overlap();
    test_cubic_bf();
    test_quintic_bf();
    test_audit_failure_leaves_outputs_unchanged();
    test_max_cutoff_near_uword_max();
    test_aliased_out_error_bound_publishes_value_last();
    test_hR_ground_truth();
    return 0;
}
