#include "zeta/zeta_internal.hpp"

#include <silex/flint/acb.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/dirichlet.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>

#include <flint/acb_dirichlet.h>
#include <flint/arb.h>

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

// 64-bit L(1, chi_D) enclosure for a large fundamental discriminant D, with
// the elapsed time printed.
void large_l1(sflint::Arb& value, const sflint::Fmpz& discriminant, slong d) {
    const auto start = std::chrono::steady_clock::now();
    assert(silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(value), sflint::FmpzConstRef(discriminant), 64));
    const double seconds = std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - start)
                                   .count();
    std::printf("L(1, chi_D) for D = %ld: %.2f s\n", d, seconds);
    assert(sflint::arb_is_positive(value));
    assert(radius_lt_2exp_si(value, -48));
}

// Imaginary D < -4: the L(1, chi) enclosure gives
// h(D) = sqrt(|D|) L(1, chi) / pi = class_number, taken from an independent
// class-number computation.
void check_large_imaginary_discriminant(slong d, slong class_number) {
    const ulong q = static_cast<ulong>(-d);
    sflint::Fmpz discriminant;
    sflint::fmpz_set_si(sflint::FmpzRef(discriminant), d);
    check_large_character(discriminant, q);

    sflint::Arb value;
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

// Real D = 1000000000061: the L(1, chi) enclosure contains 2 h R / sqrt(D)
// with h = 1 and R = 236155.816169219892016..., taken from an independent
// class-group and regulator computation.
void check_large_real_discriminant() {
    const slong d = 1000000000061;
    const ulong q = 1000000000061UL;
    sflint::Fmpz discriminant;
    sflint::fmpz_set_si(sflint::FmpzRef(discriminant), d);
    check_large_character(discriminant, q);

    sflint::Arb value;
    large_l1(value, discriminant, d);

    sflint::Arb expected;
    sflint::Arb t;
    assert(::arb_set_str(expected.raw(), "236155.816169219892016 +/- 1e-15",
                         128) == 0);
    sflint::arb_mul_2exp_si(expected, expected, 1);
    sflint::arb_sqrt_ui(t, q, 128);
    sflint::arb_div(expected, expected, t, 128);
    assert(::arb_overlaps(value.raw(), expected.raw()) != 0);
}

bool slow_tests_enabled() noexcept {
    const char* value = std::getenv("SILEX_TEST_SLOW");
    return value != nullptr && value[0] != '\0' &&
           !(value[0] == '0' && value[1] == '\0');
}

}  // namespace

int main() {
    for (const slong d : {5L, 8L, 12L, 13L, -3L, -4L, -7L, -8L, -23L, -47L,
                          1009L, -1019L, 4012L, -2008L}) {
        check_against_flint_dirichlet(d, 128);
    }
    check_against_flint_dirichlet(-100003, 128);
    check_against_flint_dirichlet(100049, 192);

    sflint::Arb unused;
    sflint::Fmpz zero;
    assert(!silex::detail::quadratic_dirichlet_l1(
            sflint::ArbRef(unused), sflint::FmpzConstRef(zero), 128));

    check_large_imaginary_discriminant(-10000000019, 39809);

    // |D| about 1e12 takes tens of seconds to minutes; opt in with
    // SILEX_TEST_SLOW=1.
    if (slow_tests_enabled()) {
        check_large_imaginary_discriminant(-1000000000039, 1113261);
        check_large_real_discriminant();
    }
    return 0;
}
