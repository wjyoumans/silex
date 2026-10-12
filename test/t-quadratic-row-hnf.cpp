// Differential tests of the maximal imaginary quadratic witness row ideal:
// the HNF built from Hensel-lifted prime roots and CRT must equal the HNF the
// ideal power/product route builds, with the same norm and denominator.

#include "class_group/maximal_imaginary_quadratic_row_hnf_internal.hpp"
#include "factor_base/factor_base_internal.hpp"

#include <silex/factor_base.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/qfb.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/prime_ideal.hpp>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <cstdio>
#include <limits>

namespace {
namespace sflint = silex::flint;
namespace detail = silex::detail;

int failures = 0;

void check(bool condition, const char* message, slong index) noexcept {
    if (!condition) {
        std::fprintf(stderr, "FAIL [%ld]: %s\n", static_cast<long>(index),
                     message);
        ++failures;
    }
}

ulong next_random(ulong& state) noexcept {
    state = state * UWORD(6364136223846793005) + UWORD(1442695040888963407);
    return state >> 11;
}

struct Fixture {
    silex::NumberField field;
    silex::Order order;
    silex::FactorBase base;
};

// incomplete: the relation backend's own base (one prime per split block);
// otherwise every degree-one prime ideal of norm <= bound (both primes per
// split block; inert primes are left out, as the backend requires).
bool make_fixture(Fixture& fixture, slong radicand, slong bound,
                  bool incomplete) noexcept {
    sflint::Fmpz value;
    sflint::fmpz_set_si(sflint::FmpzRef(value), radicand);
    fixture.field = silex::NumberField::quadratic(sflint::FmpzConstRef(value));
    sflint::fmpz_one(sflint::FmpzRef(value));
    fixture.order = silex::Order::quadratic_order(
            fixture.field, sflint::FmpzConstRef(value));
    if (!fixture.order.is_defined() || !fixture.order.is_maximal() ||
        !fixture.base.define(fixture.order)) {
        return false;
    }
    sflint::fmpz_set_si(sflint::FmpzRef(value), bound);
    return incomplete
            ? detail::FactorBaseBlockAccess::
                      build_maximal_imaginary_quadratic_for_testing(
                              fixture.base, sflint::FmpzConstRef(value))
            : detail::FactorBaseBlockAccess::build_prime_ideal_norm_bounded(
                      fixture.base, sflint::FmpzConstRef(value), false, 0,
                      false, nullptr, nullptr);
}

// Runs both routes on row and checks that they agree.  Returns false when
// either route failed (which is itself recorded as a failure).
bool routes_agree(const Fixture& fixture, const sflint::FmpzMat& row,
                  slong index) noexcept {
    sflint::FmpzMat binary_hnf(2, 2);
    sflint::FmpzMat ideal_hnf(2, 2);
    sflint::Fmpz binary_norm;
    sflint::Fmpz ideal_norm;
    sflint::Fmpz binary_denominator;
    sflint::Fmpz ideal_denominator;
    const bool ideal_ok = detail::maximal_imaginary_quadratic_row_hnf(
            ideal_hnf, ideal_norm, ideal_denominator, fixture.order,
            fixture.base, sflint::FmpzMatConstRef(row),
            detail::QuadraticRowHnfRoute::ideal_product);
    const bool binary_ok = detail::maximal_imaginary_quadratic_row_hnf(
            binary_hnf, binary_norm, binary_denominator, fixture.order,
            fixture.base, sflint::FmpzMatConstRef(row),
            detail::QuadraticRowHnfRoute::binary_form);
    check(ideal_ok, "ideal-product route failed", index);
    check(binary_ok, "binary-form route failed", index);
    if (!ideal_ok || !binary_ok) {
        return false;
    }
    check(fmpz_mat_equal(binary_hnf.raw(), ideal_hnf.raw()) != 0,
          "row HNF differs", index);
    check(fmpz_equal(binary_norm.raw(), ideal_norm.raw()) != 0,
          "row norm differs", index);
    check(fmpz_equal(binary_denominator.raw(), ideal_denominator.raw()) != 0,
          "row denominator differs", index);
    return true;
}

// Counts the split blocks of length one whose factor-base prime is the
// conjugate of the prime form's ideal p Z + ((-b + sqrt(D)) / 2) Z, that is
// whose cached omega root differs from the prime form's shift.
slong reversed_single_prime_blocks(const Fixture& fixture) noexcept {
    sflint::Fmpz discriminant;
    if (!fixture.order.discriminant(sflint::FmpzRef(discriminant))) {
        return -1;
    }
    slong reversed = 0;
    for (slong block = 0; block < fixture.base.rational_prime_block_count();
         ++block) {
        sflint::Fmpz p;
        slong start = 0;
        slong length = 0;
        if (!fixture.base.rational_prime_block(sflint::FmpzRef(p), start,
                                               length, block)) {
            return -1;
        }
        if (length != 1 || fmpz_divisible(discriminant.raw(), p.raw())) {
            continue;
        }
        const silex::PrimeIdeal* prime = fixture.base.prime_at(start);
        const sflint::Fmpz* root =
                prime == nullptr ? nullptr
                                 : detail::linear_residue_root_ptr(*prime);
        if (root == nullptr) {
            return -1;
        }
        sflint::Qfb form;
        sflint::qfb_prime_form(sflint::QfbRef(form),
                               sflint::FmpzConstRef(discriminant),
                               sflint::FmpzConstRef(p));
        sflint::Fmpz shift;
        fmpz_set(shift.raw(), sflint::qfb_b(sflint::QfbConstRef(form)).raw());
        if (fmpz_is_odd(discriminant.raw())) {
            fmpz_add_ui(shift.raw(), shift.raw(), 1);
        }
        fmpz_fdiv_q_2exp(shift.raw(), shift.raw(), 1);
        fmpz_mod(shift.raw(), shift.raw(), p.raw());
        if (!fmpz_equal(shift.raw(), root->raw())) {
            ++reversed;
        }
    }
    return reversed;
}

void random_rows(const Fixture& fixture, ulong seed, slong count,
                 slong max_exponent, ulong density, slong index_base) noexcept {
    ulong state = seed;
    const slong length = fixture.base.length();
    for (slong trial = 0; trial < count; ++trial) {
        sflint::FmpzMat row(1, length);
        for (slong j = 0; j < length; ++j) {
            if (next_random(state) % 8 >= density) {
                continue;
            }
            const slong exponent = static_cast<slong>(
                    next_random(state) %
                    static_cast<ulong>(2 * max_exponent + 1)) -
                    max_exponent;
            fmpz_set_si(fmpz_mat_entry(row.raw(), 0, j), exponent);
        }
        routes_agree(fixture, row, index_base + trial);
    }
}

void check_fixture(slong radicand, slong bound, bool incomplete,
                   ulong seed, slong index_base) noexcept {
    Fixture fixture;
    if (!make_fixture(fixture, radicand, bound, incomplete)) {
        check(false, "fixture construction failed", index_base);
        return;
    }
    const slong length = fixture.base.length();

    // Zero row: the unit ideal.
    {
        sflint::FmpzMat row(1, length);
        routes_agree(fixture, row, index_base);
    }
    // Every single column with exponents of both signs, including exponents
    // deep enough that the Hensel lift takes several Newton steps.
    const slong single_exponents[] = {1, -1, 2, -2, 3, -5, 17, -40, 129};
    for (slong j = 0; j < length; ++j) {
        for (const slong exponent : single_exponents) {
            sflint::FmpzMat row(1, length);
            fmpz_set_si(fmpz_mat_entry(row.raw(), 0, j), exponent);
            routes_agree(fixture, row, index_base + 1 + j);
        }
    }
    // Sparse and dense random rows with small and larger exponents.
    random_rows(fixture, seed, 40, 6, 3, index_base + 100);
    random_rows(fixture, seed + 1, 40, 6, 8, index_base + 200);
    random_rows(fixture, seed + 2, 10, 60, 5, index_base + 300);

    // Out-of-range exponents fail on both routes.
    if (length > 0) {
        sflint::FmpzMat hnf(2, 2);
        sflint::Fmpz norm;
        sflint::Fmpz denominator;
        sflint::FmpzMat row(1, length);
        fmpz_set_si(fmpz_mat_entry(row.raw(), 0, 0),
                    std::numeric_limits<slong>::min());
        for (const auto route : {detail::QuadraticRowHnfRoute::binary_form,
                                 detail::QuadraticRowHnfRoute::ideal_product}) {
            check(!detail::maximal_imaginary_quadratic_row_hnf(
                          hnf, norm, denominator, fixture.order, fixture.base,
                          sflint::FmpzMatConstRef(row), route),
                  "LONG_MIN exponent accepted", index_base + 400);
        }
        fmpz_one(fmpz_mat_entry(row.raw(), 0, 0));
        fmpz_mul_2exp(fmpz_mat_entry(row.raw(), 0, 0),
                      fmpz_mat_entry(row.raw(), 0, 0), 70);
        for (const auto route : {detail::QuadraticRowHnfRoute::binary_form,
                                 detail::QuadraticRowHnfRoute::ideal_product}) {
            check(!detail::maximal_imaginary_quadratic_row_hnf(
                          hnf, norm, denominator, fixture.order, fixture.base,
                          sflint::FmpzMatConstRef(row), route),
                  "exponent beyond slong accepted", index_base + 401);
        }
    }
}

// A complete base (both primes of each split p) with both columns of a
// split prime set, so the row ideal has content m > 1 or cancels to p^k.
void check_both_conjugates(slong radicand, slong bound,
                           slong index_base) noexcept {
    Fixture fixture;
    if (!make_fixture(fixture, radicand, bound, false)) {
        check(false, "complete fixture construction failed", index_base);
        return;
    }
    slong pairs = 0;
    const slong exponent_pairs[][2] = {{1, 1},  {3, 1},  {1, 3},  {-2, -2},
                                       {-3, 1}, {2, -5}, {4, 4},  {-1, 6},
                                       {7, -7}, {0, -3}, {-9, 0}, {12, 5}};
    for (slong block = 0; block < fixture.base.rational_prime_block_count();
         ++block) {
        sflint::Fmpz p;
        slong start = 0;
        slong length = 0;
        slong first = 0;
        slong second = 0;
        if (!fixture.base.rational_prime_block(sflint::FmpzRef(p), start,
                                               length, block) ||
            length != 2 ||
            !fixture.base.rational_prime_block_index(first, block, 0) ||
            !fixture.base.rational_prime_block_index(second, block, 1)) {
            continue;
        }
        ++pairs;
        for (const auto& exponents : exponent_pairs) {
            sflint::FmpzMat row(1, fixture.base.length());
            fmpz_set_si(fmpz_mat_entry(row.raw(), 0, first), exponents[0]);
            fmpz_set_si(fmpz_mat_entry(row.raw(), 0, second), exponents[1]);
            routes_agree(fixture, row, index_base + block);
        }
    }
    check(pairs > 1, "complete base has no split pair", index_base);
}

}  // namespace

int main() {
    // D = -100000007 (odd, 2 split since D = 1 mod 8): the proven benchmark
    // field, with the backend's one-prime-per-split-block base.
    {
        Fixture fixture;
        check(make_fixture(fixture, -100000007, 300, true),
              "fixture construction failed", 0);
        const slong reversed = reversed_single_prime_blocks(fixture);
        // Both orientations of a single split prime must be exercised.
        check(reversed > 0, "no reversed single-prime block", 0);
        check(reversed < fixture.base.rational_prime_block_count() - 1,
              "no form-oriented single-prime block", 0);
    }
    check_fixture(-100000007, 300, true, UWORD(0x7101), 1000);
    check_fixture(-100000007, 120, false, UWORD(0x7102), 2000);
    // D = -23 and -31: odd D with 2 split.
    check_fixture(-23, 200, true, UWORD(0x7103), 3000);
    check_fixture(-31, 200, false, UWORD(0x7104), 4000);
    // D = 4 * -30030: even D with ramified 2, 3, 5, 7, 11, 13.
    check_fixture(-30030, 250, true, UWORD(0x7105), 5000);
    check_fixture(-30030, 150, false, UWORD(0x7106), 6000);
    // D = 4 * -105 = -420 (radicand 3 mod 4): ramified 2, 3, 5, 7.
    check_fixture(-105, 250, true, UWORD(0x7107), 7000);
    check_fixture(-105, 150, false, UWORD(0x7108), 8000);
    // D = -2 * 3 * 5 * 7 * 11 * 13 * 19 * 1009 (radicand 2 mod 4).
    check_fixture(-2 * 3 * 5 * 7 * 11 * 13 * 19 * 1009, 200, true,
                  UWORD(0x7109), 9000);
    check_both_conjugates(-100000007, 120, 10000);
    check_both_conjugates(-23, 200, 11000);
    check_both_conjugates(-30030, 150, 12000);

    if (failures != 0) {
        std::fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    return 0;
}
