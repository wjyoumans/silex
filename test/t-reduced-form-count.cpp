#include "class_group/class_group_internal.hpp"
#include "class_group/reduced_form_count_internal.hpp"

#include <silex/flint/fmpz.hpp>

#include <flint/fmpz.h>
#include <flint/qfb.h>
#include <flint/ulong_extras.h>

#include <cstdio>

namespace {
namespace sflint = silex::flint;
namespace detail = silex::detail;

int failures = 0;

void check(bool condition, const char* message, slong discriminant) noexcept {
    if (!condition) {
        std::fprintf(stderr, "FAIL [D = %ld]: %s\n",
                     static_cast<long>(discriminant), message);
        ++failures;
    }
}

// D = -n is fundamental: n = 3 mod 4 squarefree, or n = 4m with -m = 2 or
// 3 mod 4 and m squarefree, written independently of the library gate.
bool reference_fundamental(ulong n) noexcept {
    if (n % 4 == 3) {
        return n_is_squarefree(n) != 0;
    }
    if (n % 4 != 0) {
        return false;
    }
    const ulong m = n / 4;
    const ulong minus_m_mod_4 = (4 - m % 4) % 4;
    return (minus_m_mod_4 == 2 || minus_m_mod_4 == 3) &&
           n_is_squarefree(m) != 0;
}

// The number of reduced forms FLINT `qfb_reduced_forms` lists for D.
slong listed_reduced_forms(slong discriminant) noexcept {
    qfb* forms = nullptr;
    const slong count = qfb_reduced_forms(&forms, discriminant);
    if (forms != nullptr) {
        qfb_array_clear(&forms, count);
    }
    return count;
}

bool class_number(slong& out, slong discriminant) noexcept {
    sflint::Fmpz d;
    sflint::fmpz_set_si(sflint::FmpzRef(d), discriminant);
    sflint::Fmpz h;
    if (!detail::imaginary_quadratic_class_number(sflint::FmpzRef(h),
                                                  sflint::FmpzConstRef(d))) {
        return false;
    }
    if (::fmpz_fits_si(h.raw()) == 0) {
        return false;
    }
    out = ::fmpz_get_si(h.raw());
    return true;
}

// The count equals FLINT's listing on every fundamental |D| <= 2e5, and the
// gate agrees with the reference definition on every |D| in that range.
void check_against_listing() noexcept {
    slong compared = 0;
    for (ulong n = 1; n <= 200000; ++n) {
        const slong d = -static_cast<slong>(n);
        const bool fundamental = reference_fundamental(n);
        check(detail::imaginary_quadratic_discriminant_is_fundamental(n) ==
                      fundamental,
              "fundamental gate", d);
        ulong count = 0;
        const bool counted =
                detail::imaginary_quadratic_fundamental_reduced_form_count(
                        count, n);
        check(counted == fundamental, "count refuses exactly the "
                                      "non-fundamental D", d);
        if (!fundamental) {
            continue;
        }
        check(counted && count == static_cast<ulong>(listed_reduced_forms(d)),
              "count equals qfb_reduced_forms", d);
        ++compared;
    }
    check(compared == 60800, "number of fundamental D compared", 0);
}

// Fixed rows: the count, and the class number through the shipped switch.
// h(-17592186044415) is the largest fundamental |D| row below the listing
// cutoff 2^44, and D = -(2^44 + 3) is counted above it (the count cutoff is
// 2^54); h(-(2^44 + 3)) = 1148064 is taken from an independent class-number
// computation.
void check_fixed_rows() noexcept {
    struct Row {
        slong discriminant;
        ulong class_number;
    };
    const Row rows[] = {
            {-3, 1},
            {-4, 1},
            {-40000000008L, 45992},
            {-10000000799L, 202707},
            {-1000000000271L, 1499213},
            {-17592186044415L, 3384192},
            {-17592186044419L, 1148064},
    };
    for (const Row& row : rows) {
        const ulong n = static_cast<ulong>(-row.discriminant);
        ulong count = 0;
        check(detail::imaginary_quadratic_fundamental_reduced_form_count(
                      count, n) &&
                      count == row.class_number,
              "fixed count", row.discriminant);
        slong h = 0;
        check(class_number(h, row.discriminant) &&
                      h == static_cast<slong>(row.class_number),
              "fixed class number", row.discriminant);
    }
}

// A non-fundamental D is refused by the count and takes FLINT's listing.
// Above the listing cutoff 2^44 the count also refuses it, so the switch
// passes it to the class number formula, which gives no h(D) for this value
// (the slow test t-zeta-quadratic-l1; not a property of the formula): D = -(2^44 + 23) = -3^2 * 59 * 4457 *
// 7433317.
void check_non_fundamental_fallback() noexcept {
    {
        ulong count = 0;
        check(!detail::imaginary_quadratic_fundamental_reduced_form_count(
                      count, UWORD(17592186044439)),
              "count refuses a non-fundamental D above 2^44",
              -17592186044439L);
    }
    for (const slong d : {-12L, -16L, -40000000012L}) {
        const ulong n = static_cast<ulong>(-d);
        ulong count = 0;
        check(!detail::imaginary_quadratic_fundamental_reduced_form_count(
                      count, n),
              "count refuses a non-fundamental D", d);
        slong h = 0;
        check(class_number(h, d) && h == listed_reduced_forms(d),
              "non-fundamental D gives the listing's count", d);
    }
    slong h = 0;
    check(class_number(h, -40000000012L) && h == 31614,
          "listing's count for D = -40000000012", -40000000012L);
}

}  // namespace

int main() {
    check_against_listing();
    check_fixed_rows();
    check_non_fundamental_fallback();
    if (failures != 0) {
        std::fprintf(stderr, "%d failures\n", failures);
        return 1;
    }
    return 0;
}
