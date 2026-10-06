#include "zeta/zeta_internal.hpp"

#include <silex/class_group.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_unit.hpp>
#include <silex/flint/acb.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/dirichlet.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>

#include <flint/acb_dirichlet.h>
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

namespace {
namespace sflint = silex::flint;

bool radius_lt_2exp_si(const sflint::Arb& value, slong exponent) noexcept {
    sflint::Arb radius;
    sflint::Arb target;
    sflint::arb_get_rad_arb(radius, value);
    sflint::arb_one(target);
    sflint::arb_mul_2exp_si(target, target, exponent);
    return sflint::arb_lt(radius, target);
}

int kronecker(const sflint::Fmpz& discriminant, ulong n) noexcept {
    sflint::Fmpz zn;
    sflint::fmpz_set_ui(sflint::FmpzRef(zn), n);
    return sflint::fmpz_kronecker(sflint::FmpzConstRef(discriminant),
                                  sflint::FmpzConstRef(zn));
}

// Value of the FLINT Dirichlet character as -1, 0, or 1 when it is real.
int dirichlet_real_value(const sflint::DirichletGroup& group,
                         const sflint::DirichletChar& character,
                         ulong n) noexcept {
    const ulong got = sflint::dirichlet_chi(group, character, n);
    if (got == DIRICHLET_CHI_NULL) {
        return 0;
    }
    if (got == 0) {
        return 1;
    }
    return group.exponent() % 2 == 0 && got == group.exponent() / 2 ? -1
                                                                    : 2;
}

// For a fundamental discriminant D, the Kronecker character used by the
// quadratic L(1, chi) route is the real primitive FLINT Dirichlet character
// modulo |D| on every residue, and its L(1, chi) enclosure overlaps FLINT's
// `acb_dirichlet_l_fmpq`.
void check_against_flint_dirichlet(slong d, slong precision) {
    sflint::Fmpz discriminant;
    sflint::fmpz_set_si(sflint::FmpzRef(discriminant), d);
    const ulong q = static_cast<ulong>(d < 0 ? -d : d);

    sflint::DirichletGroup group(q);
    assert(group.is_initialized());
    sflint::DirichletChar candidate(group);
    sflint::DirichletChar character(group);
    bool found = false;
    for (ulong j = 0; j < group.character_count() && !found; ++j) {
        sflint::dirichlet_char_index(candidate, group, j);
        if (!sflint::dirichlet_char_is_real(group, candidate) ||
            !sflint::dirichlet_char_is_primitive(group, candidate)) {
            continue;
        }
        bool ok = true;
        for (ulong n = 1; n <= q && ok; ++n) {
            ok = dirichlet_real_value(group, candidate, n) ==
                 silex::detail::quadratic_character(
                         sflint::FmpzConstRef(discriminant), n);
        }
        if (ok) {
            sflint::dirichlet_char_set(character, group, candidate);
            found = true;
        }
    }
    assert(found);

    sflint::Arb value;
    assert(silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(value), sflint::FmpzConstRef(discriminant),
            precision));
    assert(sflint::arb_is_positive(value));
    assert(radius_lt_2exp_si(value, -(precision - 16)));

    sflint::Fmpq one;
    sflint::Acb reference;
    sflint::fmpq_one(one);
    sflint::acb_dirichlet_l_fmpq(reference, sflint::FmpqConstRef(one),
                                 group.raw(), character.raw(), precision);
    assert(::arb_overlaps(value.raw(), acb_realref(reference.raw())) != 0);
}

// Wrapper check: `quadratic_character` returns the Kronecker symbol (D/n) on
// sampled n for a large D, and that symbol is periodic modulo |D| and
// completely multiplicative there.  This checks the wrapper only; that (D/.)
// is the primitive character of the L-function is established by the
// small-D comparison with FLINT's Dirichlet characters above and by the
// class-number and regulator checks below.
void check_large_character(const sflint::Fmpz& discriminant, ulong q) {
    const auto chi = [&discriminant](ulong n) noexcept {
        return silex::detail::quadratic_character(
                sflint::FmpzConstRef(discriminant), n);
    };
    ulong n = 1;
    for (ulong i = 0; i < 2000; ++i) {
        assert(chi(n) == kronecker(discriminant, n));
        assert(chi(n + q) == chi(n));
        const ulong m = 2 * i + 3;
        const ulong k = n % 1000003;
        assert(chi(m * k) == chi(m) * chi(k));
        n = (n * 6364136223846793005UL + 1442695040888963407UL) % (4 * q) + 1;
    }
    assert(chi(q) == 0);
}

// Both sums of the evaluation for D used Taylor stepping of their
// incomplete gamma values Gamma(a, z): a = 1 and a = 1/2 for imaginary D,
// a = 1/2 and a = 0 for real D (counters indexed by shift = 2a).
void assert_stepping_used(const silex::detail::QuadraticAfeStepAudit& audit,
                          slong d) {
    const ulong unused = d < 0 ? 0 : 2;
    for (ulong shift = 0; shift < 3; ++shift) {
        if (shift == unused) {
            assert(audit.blocks[shift] == 0);
            assert(audit.stepped_values[shift] == 0);
        } else {
            assert(audit.blocks[shift] > 0);
            assert(audit.stepped_values[shift] > audit.blocks[shift]);
        }
    }
}

// 64-bit L(1, chi_D) enclosure for a large fundamental discriminant D, with
// the elapsed time printed.  Both sums use stepping at these sizes.
void large_l1(sflint::Arb& value, const sflint::Fmpz& discriminant, slong d) {
    const auto start = std::chrono::steady_clock::now();
    silex::detail::QuadraticAfeStepAudit audit;
    assert(silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(value), sflint::FmpzConstRef(discriminant), 64,
            &audit));
    assert_stepping_used(audit, d);
    const double seconds = std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - start)
                                   .count();
    std::printf("L(1, chi_D) for D = %ld: %.2f s\n", d, seconds);
    assert(sflint::arb_is_positive(value));
    assert(radius_lt_2exp_si(value, -48));
}

// Imaginary D < -4: the L(1, chi) enclosure gives
// h(D) = sqrt(|D|) L(1, chi) / pi = class_number, taken from an independent
// class-number computation.  The enclosure is returned in value.
void check_large_imaginary_discriminant(sflint::Arb& value, slong d,
                                        slong class_number) {
    const ulong q = static_cast<ulong>(-d);
    sflint::Fmpz discriminant;
    sflint::fmpz_set_si(sflint::FmpzRef(discriminant), d);
    check_large_character(discriminant, q);

    large_l1(value, discriminant, d);

    sflint::Arb h;
    sflint::Arb t;
    sflint::arb_sqrt_ui(t, q, 128);
    sflint::arb_mul(h, value, t, 128);
    sflint::arb_const_pi(t, 128);
    sflint::arb_div(h, h, t, 128);
    assert(sflint::arb_contains_si(h, class_number));
    assert(radius_lt_2exp_si(h, -2));
}

// The same imaginary D at a higher precision: the enclosure has the
// precision's radius, still gives h(D) = class_number, and overlaps the
// 64-bit enclosure low.  This is an agreement test for the
// precision-dependent stepping parameters, not a check of the truncation
// bound.
void check_large_imaginary_discriminant_at(const sflint::Arb& low, slong d,
                                           slong class_number,
                                           slong precision) {
    const ulong q = static_cast<ulong>(-d);
    sflint::Fmpz discriminant;
    sflint::fmpz_set_si(sflint::FmpzRef(discriminant), d);

    sflint::Arb value;
    assert(silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(value), sflint::FmpzConstRef(discriminant),
            precision));
    assert(sflint::arb_is_positive(value));
    assert(radius_lt_2exp_si(value, -(precision - 16)));
    assert(::arb_overlaps(value.raw(), low.raw()) != 0);

    sflint::Arb h;
    sflint::Arb t;
    sflint::arb_sqrt_ui(t, q, precision);
    sflint::arb_mul(h, value, t, precision);
    sflint::arb_const_pi(t, precision);
    sflint::arb_div(h, h, t, precision);
    assert(sflint::arb_contains_si(h, class_number));
    assert(radius_lt_2exp_si(h, -2));
}

// Real D > 0 with class number one: the L(1, chi) enclosure contains
// 2 R / sqrt(D), with h = 1 and the regulator R (a decimal ball) taken from
// an independent class-group and regulator computation.
void check_large_real_discriminant(slong d, const char* regulator) {
    const ulong q = static_cast<ulong>(d);
    sflint::Fmpz discriminant;
    sflint::fmpz_set_si(sflint::FmpzRef(discriminant), d);
    check_large_character(discriminant, q);

    sflint::Arb value;
    large_l1(value, discriminant, d);

    sflint::Arb expected;
    sflint::Arb t;
    assert(::arb_set_str(expected.raw(), regulator, 128) == 0);
    sflint::arb_mul_2exp_si(expected, expected, 1);
    sflint::arb_sqrt_ui(t, q, 128);
    sflint::arb_div(expected, expected, t, 128);
    assert(::arb_overlaps(value.raw(), expected.raw()) != 0);
}

// Public route: D = -100000020 (radicand -25000005), h(D) = 5056, taken
// from an independent class-number computation.  The previous L(1, chi)
// route searched the Dirichlet group modulo |D| for the Kronecker
// character, O(|D|) work that takes about a minute here and grows linearly;
// the approximate functional equation takes a fraction of a second.  The
// proven paired computation uses the exact imaginary-quadratic class number
// and records no analytic check;
// try_certify_class_unit_with_zeta then evaluates hR through the quadratic
// L(1, chi) route and records an unconditional (`proven`) analytic check.
void check_public_class_unit_with_zeta() {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), -25000005);
    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    silex::Order order = silex::Order::equation_order(field);
    assert(field.is_defined() && order.is_defined() && order.is_maximal());
    sflint::Fmpz discriminant;
    assert(order.discriminant(sflint::FmpzRef(discriminant)));
    assert(sflint::fmpz_equal_si(discriminant, -100000020));

    sflint::Fmpz factor_base_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(factor_base_bound), 2);
    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::proven;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(
            class_group, order, sflint::FmpzConstRef(factor_base_bound),
            options, 128));
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.analytic_class_regulator_certification() ==
           silex::CertificationMode::unknown);
    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 5056));

    const auto start = std::chrono::steady_clock::now();
    assert(class_group.try_certify_class_unit_with_zeta(units, 128));
    const double seconds = std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - start)
                                   .count();
    std::printf("try_certify_class_unit_with_zeta for D = -100000020: "
                "%.2f s\n",
                seconds);
    assert(class_group.analytic_class_regulator_status() ==
           silex::ProofState::verified);
    assert(class_group.analytic_class_regulator_certification() ==
           silex::CertificationMode::proven);
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(units.certification_status() == silex::CertificationMode::proven);
}

// Maximal order of Q(sqrt(radicand)) for a squarefree radicand, checked to
// have discriminant `expected`.
silex::Order quadratic_maximal_order(silex::NumberField& field,
                                     const char* radicand,
                                     const char* expected) {
    sflint::Fmpz r;
    assert(::fmpz_set_str(r.raw(), radicand, 10) == 0);
    field = silex::NumberField::quadratic(sflint::FmpzConstRef(r));
    sflint::Fmpz one;
    sflint::fmpz_one(sflint::FmpzRef(one));
    silex::Order order =
            silex::Order::quadratic_order(field, sflint::FmpzConstRef(one));
    assert(field.is_defined() && order.is_defined() && order.is_maximal());
    sflint::Fmpz discriminant;
    sflint::Fmpz want;
    assert(order.discriminant(sflint::FmpzRef(discriminant)));
    assert(::fmpz_set_str(want.raw(), expected, 10) == 0);
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(discriminant),
                              sflint::FmpzConstRef(want)));
    return order;
}

// The quadratic L(1, chi) route takes only |D| < 2^44 = 17592186044416;
// larger fields fall back to Belabas-Friedman.  The route predicate is
// checked on both sides of the cap and at |D| about 1e12 (the fields of the
// slow tests below), without evaluating an L-value.
void check_route_size_cap() {
    struct Case {
        const char* radicand;
        const char* discriminant;
        bool quadratic_route;
    };
    const Case cases[] = {
            {"-4398046511105", "-17592186044420", false},  // -(2^44 + 4)
            {"17592186044417", "17592186044417", false},   // 2^44 + 1
            {"-4398046511102", "-17592186044408", true},   // -(2^44 - 8)
            {"17592186044413", "17592186044413", true},    // 2^44 - 3
            {"-1000000000039", "-1000000000039", true},
            {"1000000000061", "1000000000061", true},
    };
    for (const Case& c : cases) {
        silex::NumberField field;
        const silex::Order order =
                quadratic_maximal_order(field, c.radicand, c.discriminant);
        assert(silex::detail::zeta_unconditional_route_available(order) ==
               c.quadratic_route);
    }

    // Above the cap the residue comes from the GRH-conditional
    // Belabas-Friedman fallback.
    silex::NumberField field;
    const silex::Order order = quadratic_maximal_order(
            field, "-4398046511105", "-17592186044420");
    sflint::Arb product;
    bool unconditional = true;
    assert(silex::detail::zeta_class_regulator_product_with_diagnostics(
            sflint::ArbRef(product), order, 64, nullptr, nullptr, nullptr,
            &unconditional));
    assert(!unconditional);
    assert(sflint::arb_is_positive(product));
}

bool slow_tests_enabled() noexcept {
    const char* value = std::getenv("SILEX_TEST_SLOW");
    return value != nullptr && value[0] != '\0' &&
           !(value[0] == '0' && value[1] == '\0');
}

// Invalid input fails closed and leaves `out` untouched.
void check_fails_closed(const sflint::Fmpz& discriminant, slong precision) {
    sflint::Arb out;
    ::arb_set_si(out.raw(), 7);
    assert(!silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(out), sflint::FmpzConstRef(discriminant),
            precision));
    sflint::Arb expected;
    ::arb_set_si(expected.raw(), 7);
    assert(::arb_equal(out.raw(), expected.raw()) != 0);
}

void check_failure_paths() {
    sflint::Fmpz d;
    // D = 0, and |D| = 1, 2 (modulus q < 3).
    for (const slong value : {0L, 1L, -1L, 2L, -2L}) {
        sflint::fmpz_set_si(sflint::FmpzRef(d), value);
        check_fails_closed(d, 128);
    }
    // |D| = 2^64 + 3 does not fit a ulong.
    sflint::fmpz_set_ui(sflint::FmpzRef(d), 1);
    sflint::fmpz_mul_2exp(sflint::FmpzRef(d), sflint::FmpzConstRef(d), 64);
    sflint::fmpz_add_ui(sflint::FmpzRef(d), sflint::FmpzConstRef(d), 3);
    sflint::fmpz_neg(sflint::FmpzRef(d), sflint::FmpzConstRef(d));
    check_fails_closed(d, 128);
    // Non-positive precision with a valid discriminant.
    sflint::fmpz_set_si(sflint::FmpzRef(d), -23);
    check_fails_closed(d, 0);
    check_fails_closed(d, -5);
}


// Gamma(shift / 2, pi n^2 / q) from the stepper of the quadratic L(1, chi)
// sums, against FLINT's `arb_hypgeom_gamma_upper` at 256 bits, over several
// chained blocks of the default length and over blocks forced to the
// minimum length 4.
struct GammaStepCase {
    ulong shift;
    ulong q;
    ulong first;
    slong count;
    ulong max_block;
    // A lower bound on the number of blocks started.
    ulong min_blocks;
};

void check_gamma_stepping_against_hypgeom() {
    const double target = -70.0;
    const slong prec = 80;
    const slong guard = 12;
    const GammaStepCase cases[] = {
            // a = 1 with z about 3e-4: default blocks of about 1500 indices.
            {2, UWORD(10000000000), 1000, 4000, 0, 2},
            {2, UWORD(10000000000), 1000, 41, 4, 10},
            // a = 1/2 and a = 0 with z about 3.1, where rho = 1: default
            // blocks of about 1500 indices.
            {1, UWORD(100000000000000), 10000000, 4000, 0, 2},
            {1, UWORD(100000000000000), 10000000, 41, 4, 10},
            {0, UWORD(100000000000000), 10000000, 4000, 0, 2},
            {0, UWORD(100000000000000), 10000000, 41, 4, 10},
            // a = 1/2 and a = 0 with z about 0.13, where rho = 1 / z:
            // default blocks of about 9 indices.
            {1, UWORD(10000000000), 20000, 100, 0, 8},
            {0, UWORD(10000000000), 20000, 100, 0, 8},
    };
    sflint::Arb z;
    sflint::Arb a;
    sflint::Arb reference;
    for (const GammaStepCase& c : cases) {
        sflint::ArbVec values(c.count);
        silex::detail::QuadraticAfeGammaStepCount count;
        assert(silex::detail::quadratic_afe_gamma_stepped(
                values, count, c.shift, c.q, c.first, c.max_block, target,
                prec, guard));
        // Every value after the first, direct anchor comes from a block.
        assert(count.stepped_values == static_cast<ulong>(c.count - 1));
        assert(count.blocks >= c.min_blocks);
        if (c.max_block != 0) {
            assert(count.blocks >=
                   static_cast<ulong>(c.count - 1) / c.max_block);
        }

        ::arb_set_ui(a.raw(), c.shift);
        ::arb_mul_2exp_si(a.raw(), a.raw(), -1);
        for (slong i = 0; i < c.count; ++i) {
            const ulong n = c.first + static_cast<ulong>(i);
            ::arb_const_pi(z.raw(), 256);
            ::arb_mul_ui(z.raw(), z.raw(), n, 256);
            ::arb_mul_ui(z.raw(), z.raw(), n, 256);
            ::arb_div_ui(z.raw(), z.raw(), c.q, 256);
            ::arb_hypgeom_gamma_upper(reference.raw(), a.raw(), z.raw(), 0,
                                      256);
            assert(radius_lt_2exp_si(reference, -200));
            sflint::Arb value;
            ::arb_set(value.raw(), values.data() + i);
            assert(::arb_overlaps(value.raw(), reference.raw()) != 0);
            assert(radius_lt_2exp_si(value, static_cast<slong>(target)));
        }
    }
}

// L(1, chi_D) with stepping against the evaluation with stepping turned off
// (every gamma value evaluated directly, as before stepping was added).  The
// two enclosures are bit-identical in every case below.
struct StepComparisonCase {
    slong discriminant;
    slong precision;
    // Whether both sums use stepping.  Below |D| about 1e10 the blocks are
    // too short for most a, and at |D| about 1e8 at most a = 1 steps.
    bool both_sums_step;
};

void check_stepping_matches_direct(const StepComparisonCase& c) {
    sflint::Fmpz d;
    sflint::fmpz_set_si(sflint::FmpzRef(d), c.discriminant);
    sflint::Arb stepped;
    sflint::Arb direct;
    silex::detail::QuadraticAfeStepAudit on;
    silex::detail::QuadraticAfeStepAudit off;
    off.stepping = false;
    assert(silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(stepped), sflint::FmpzConstRef(d), c.precision,
            &on));
    assert(silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(direct), sflint::FmpzConstRef(d), c.precision,
            &off));
    assert(::arb_equal(stepped.raw(), direct.raw()) != 0);
    for (ulong shift = 0; shift < 3; ++shift) {
        assert(off.blocks[shift] == 0);
        assert(off.stepped_values[shift] == 0);
    }
    if (c.both_sums_step) {
        assert_stepping_used(on, c.discriminant);
    }
}

// Fourteen fundamental discriminants from -1019 to about 1e12 at 64 to 192
// bits.  The cases that take more than about a second (up to about 90 s
// for D = 1000000000061 without stepping) run with SILEX_TEST_SLOW=1.
void check_stepping_before_after(bool slow) {
    const std::initializer_list<StepComparisonCase> fast_cases = {
            {-1019, 64, false},      {-100003, 128, false},
            {100049, 192, false},    {-99999971, 64, false},
            {-111546435, 64, false},
    };
    const std::initializer_list<StepComparisonCase> slow_cases = {
            {99999989, 64, false},         {48612265, 128, false},
            {-10000000019, 64, true},      {-10000000019, 128, true},
            {-10000000019, 192, true},     {100000000057, 64, true},
            {-100000000003, 64, true},     {-1000000000039, 64, true},
            {1000000000061, 64, true},
    };
    for (const StepComparisonCase& c : fast_cases) {
        check_stepping_matches_direct(c);
    }
    if (slow) {
        for (const StepComparisonCase& c : slow_cases) {
            check_stepping_matches_direct(c);
        }
    }
}

}  // namespace

int main() {
    for (const slong d : {5L, 8L, 12L, 13L, -3L, -4L, -7L, -8L, -23L, -47L,
                          1009L, -1019L, 4012L, -2008L}) {
        check_against_flint_dirichlet(d, 128);
    }
    check_against_flint_dirichlet(-100003, 128);
    check_against_flint_dirichlet(100049, 192);

    check_failure_paths();
    check_route_size_cap();
    check_gamma_stepping_against_hypgeom();
    check_stepping_before_after(slow_tests_enabled());

    sflint::Arb imaginary_l1;
    check_large_imaginary_discriminant(imaginary_l1, -10000000019, 39809);
    check_large_imaginary_discriminant_at(imaginary_l1, -10000000019, 39809,
                                          128);
    check_public_class_unit_with_zeta();

    // |D| about 1e11 runs Taylor stepping of the incomplete gamma values over
    // most of both sums, including Gamma(0, z) = E_1(z), which only real D
    // use (D = -10000000019 above covers Gamma(1, z) and Gamma(1/2, z)).
    // h = 1 and R = 468981.557709373150346..., taken from an independent
    // class-group and regulator computation with an unconditional
    // class-group certificate.
    check_large_real_discriminant(100000000057,
                                  "468981.557709373150346155393 +/- 1e-21");

    // |D| about 1e12 is slow; opt in with
    // SILEX_TEST_SLOW=1.
    if (slow_tests_enabled()) {
        check_large_imaginary_discriminant(imaginary_l1, -1000000000039,
                                           1113261);
        // h = 1 and R = 236155.816169219892016..., taken from an independent
        // class-group and regulator computation.
        check_large_real_discriminant(1000000000061,
                                      "236155.816169219892016 +/- 1e-15");
    }
    return 0;
}
