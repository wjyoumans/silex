#include <silex/flint/acb.hpp>
#include <silex/flint/acb_mat.hpp>
#include <silex/flint/acb_poly.hpp>
#include <silex/flint/acb_vec.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/dirichlet.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/fmpz_poly_factor.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/flint/fq.hpp>
#include <silex/flint/nf.hpp>
#include <silex/flint/nf_elem.hpp>
#include <silex/flint/nmod_mat.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/flint/nmod_poly_factor.hpp>
#include <silex/flint/qfb.hpp>

#include <type_traits>
#include <utility>

namespace {

bool equals_si(const fmpz_t value, slong expected) noexcept {
    return fmpz_equal_si(value, expected) != 0;
}

bool equals_frac(const fmpq_t value, slong numerator, ulong denominator) noexcept {
    return fmpz_equal_si(fmpq_numref(value), numerator) != 0 &&
           fmpz_equal_ui(fmpq_denref(value), denominator) != 0;
}

bool equals_poly2(const fmpz_poly_t value, slong constant, slong linear) noexcept {
    return fmpz_poly_length(value) == 2 &&
           fmpz_poly_get_coeff_si(value, 0) == constant &&
           fmpz_poly_get_coeff_si(value, 1) == linear;
}

bool equals_qpoly2(const fmpq_poly_t value, slong constant, slong linear) noexcept {
    fmpq_t coeff;
    fmpq_init(coeff);

    fmpq_poly_get_coeff_fmpq(coeff, value, 0);
    const bool constant_ok = fmpq_equal_si(coeff, constant) != 0;

    fmpq_poly_get_coeff_fmpq(coeff, value, 1);
    const bool linear_ok = fmpq_equal_si(coeff, linear) != 0;

    fmpq_clear(coeff);
    return fmpq_poly_length(value) == 2 && constant_ok && linear_ok;
}

bool equals_mat2(const fmpz_mat_t value,
                 slong a00,
                 slong a01,
                 slong a10,
                 slong a11) noexcept {
    return fmpz_mat_nrows(value) == 2 && fmpz_mat_ncols(value) == 2 &&
           fmpz_equal_si(fmpz_mat_entry(value, 0, 0), a00) != 0 &&
           fmpz_equal_si(fmpz_mat_entry(value, 0, 1), a01) != 0 &&
           fmpz_equal_si(fmpz_mat_entry(value, 1, 0), a10) != 0 &&
           fmpz_equal_si(fmpz_mat_entry(value, 1, 1), a11) != 0;
}

bool equals_qmat2(const fmpq_mat_t value,
                  slong a00,
                  slong a01,
                  slong a10,
                  slong a11) noexcept {
    return fmpq_mat_nrows(value) == 2 && fmpq_mat_ncols(value) == 2 &&
           fmpq_equal_si(fmpq_mat_entry(value, 0, 0), a00) != 0 &&
           fmpq_equal_si(fmpq_mat_entry(value, 0, 1), a01) != 0 &&
           fmpq_equal_si(fmpq_mat_entry(value, 1, 0), a10) != 0 &&
           fmpq_equal_si(fmpq_mat_entry(value, 1, 1), a11) != 0;
}

bool equals_arb_mat2(const arb_mat_t value,
                     slong a00,
                     slong a01,
                     slong a10,
                     slong a11) noexcept {
    return arb_mat_nrows(value) == 2 && arb_mat_ncols(value) == 2 &&
           arb_equal_si(arb_mat_entry(value, 0, 0), a00) != 0 &&
           arb_equal_si(arb_mat_entry(value, 0, 1), a01) != 0 &&
           arb_equal_si(arb_mat_entry(value, 1, 0), a10) != 0 &&
           arb_equal_si(arb_mat_entry(value, 1, 1), a11) != 0;
}

bool equals_acb_mat2(const acb_mat_t value,
                     slong a00,
                     slong a01,
                     slong a10,
                     slong a11) noexcept {
    return acb_mat_nrows(value) == 2 && acb_mat_ncols(value) == 2 &&
           acb_equal_si(acb_mat_entry(value, 0, 0), a00) != 0 &&
           acb_equal_si(acb_mat_entry(value, 0, 1), a01) != 0 &&
           acb_equal_si(acb_mat_entry(value, 1, 0), a10) != 0 &&
           acb_equal_si(acb_mat_entry(value, 1, 1), a11) != 0;
}

bool equals_nmod_mat2(const nmod_mat_t value,
                      ulong a00,
                      ulong a01,
                      ulong a10,
                      ulong a11) noexcept {
    return nmod_mat_nrows(value) == 2 && nmod_mat_ncols(value) == 2 &&
           nmod_mat_get_entry(value, 0, 0) == a00 &&
           nmod_mat_get_entry(value, 0, 1) == a01 &&
           nmod_mat_get_entry(value, 1, 0) == a10 &&
           nmod_mat_get_entry(value, 1, 1) == a11;
}

int test_fmpz() {
    silex::flint::Fmpz a;
    fmpz_set_si(a.raw(), 7);
    if (!equals_si(a.raw(), 7)) {
        return 1;
    }

    silex::flint::Fmpz b;
    fmpz_set_si(b.raw(), -3);
    if (!equals_si(b.raw(), -3)) {
        return 1;
    }

    a.swap(b);
    if (!equals_si(a.raw(), -3) || !equals_si(b.raw(), 7)) {
        return 1;
    }

    silex::flint::Fmpz c(std::move(a));
    if (!equals_si(c.raw(), -3)) {
        return 1;
    }
    fmpz_set_si(a.raw(), 11);
    if (!equals_si(a.raw(), 11)) {
        return 1;
    }

    silex::flint::Fmpz d;
    fmpz_set_si(d.raw(), 19);
    d = std::move(c);
    if (!equals_si(d.raw(), -3)) {
        return 1;
    }
    fmpz_set_si(c.raw(), 23);
    if (!equals_si(c.raw(), 23)) {
        return 1;
    }

    swap(b, d);
    if (!equals_si(b.raw(), -3) || !equals_si(d.raw(), 7)) {
        return 1;
    }

    const silex::flint::Fmpz& const_d = d;
    if (!equals_si(const_d.raw(), 7)) {
        return 1;
    }

    silex::flint::FmpzRef d_ref(d);
    fmpz_set_si(d_ref.raw(), 29);
    if (!equals_si(d.raw(), 29)) {
        return 1;
    }

    const silex::flint::FmpzConstRef d_const_ref(d);
    if (!equals_si(d_const_ref.raw(), 29)) {
        return 1;
    }

    fmpz_t external;
    fmpz_init(external);
    {
        silex::flint::FmpzRef external_ref(external);
        fmpz_set_si(external_ref.raw(), -41);

        silex::flint::FmpzConstRef external_const_ref(external);
        if (!equals_si(external_const_ref.raw(), -41)) {
            fmpz_clear(external);
            return 1;
        }
    }

    fmpz_set_si(external, 43);
    if (!equals_si(external, 43)) {
        fmpz_clear(external);
        return 1;
    }
    fmpz_clear(external);

    return 0;
}

int test_fmpq() {
    silex::flint::Fmpq a;
    fmpq_set_si(a.raw(), 7, 3);
    if (!equals_frac(a.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::Fmpq b;
    fmpq_set_si(b.raw(), -5, 2);
    if (!equals_frac(b.raw(), -5, 2)) {
        return 1;
    }

    a.swap(b);
    if (!equals_frac(a.raw(), -5, 2) || !equals_frac(b.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::Fmpq c(std::move(a));
    if (!equals_frac(c.raw(), -5, 2)) {
        return 1;
    }
    fmpq_set_si(a.raw(), 11, 4);
    if (!equals_frac(a.raw(), 11, 4)) {
        return 1;
    }

    silex::flint::Fmpq d;
    fmpq_set_si(d.raw(), 19, 6);
    d = std::move(c);
    if (!equals_frac(d.raw(), -5, 2)) {
        return 1;
    }
    fmpq_set_si(c.raw(), 23, 8);
    if (!equals_frac(c.raw(), 23, 8)) {
        return 1;
    }

    swap(b, d);
    if (!equals_frac(b.raw(), -5, 2) || !equals_frac(d.raw(), 7, 3)) {
        return 1;
    }

    const silex::flint::Fmpq& const_d = d;
    if (!equals_frac(const_d.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpqRef d_ref(d);
    fmpq_set_si(d_ref.raw(), 29, 10);
    if (!equals_frac(d.raw(), 29, 10)) {
        return 1;
    }

    const silex::flint::FmpqConstRef d_const_ref(d);
    if (!equals_frac(d_const_ref.raw(), 29, 10)) {
        return 1;
    }

    fmpq_t external;
    fmpq_init(external);
    {
        silex::flint::FmpqRef external_ref(external);
        fmpq_set_si(external_ref.raw(), -41, 12);

        silex::flint::FmpqConstRef external_const_ref(external);
        if (!equals_frac(external_const_ref.raw(), -41, 12)) {
            fmpq_clear(external);
            return 1;
        }
    }

    fmpq_set_si(external, 43, 14);
    if (!equals_frac(external, 43, 14)) {
        fmpq_clear(external);
        return 1;
    }
    fmpq_clear(external);

    return 0;
}

int test_fmpz_poly() {
    silex::flint::FmpzPoly a;
    fmpz_poly_set_coeff_si(a.raw(), 0, 7);
    fmpz_poly_set_coeff_si(a.raw(), 1, 3);
    if (!equals_poly2(a.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpzPoly b;
    fmpz_poly_set_coeff_si(b.raw(), 0, -5);
    fmpz_poly_set_coeff_si(b.raw(), 1, 2);
    if (!equals_poly2(b.raw(), -5, 2)) {
        return 1;
    }

    a.swap(b);
    if (!equals_poly2(a.raw(), -5, 2) || !equals_poly2(b.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpzPoly c(std::move(a));
    if (!equals_poly2(c.raw(), -5, 2)) {
        return 1;
    }
    fmpz_poly_set_coeff_si(a.raw(), 0, 11);
    fmpz_poly_set_coeff_si(a.raw(), 1, 4);
    if (!equals_poly2(a.raw(), 11, 4)) {
        return 1;
    }

    silex::flint::FmpzPoly d;
    fmpz_poly_set_coeff_si(d.raw(), 0, 19);
    fmpz_poly_set_coeff_si(d.raw(), 1, 6);
    d = std::move(c);
    if (!equals_poly2(d.raw(), -5, 2)) {
        return 1;
    }
    fmpz_poly_set_coeff_si(c.raw(), 0, 23);
    fmpz_poly_set_coeff_si(c.raw(), 1, 8);
    if (!equals_poly2(c.raw(), 23, 8)) {
        return 1;
    }

    swap(b, d);
    if (!equals_poly2(b.raw(), -5, 2) || !equals_poly2(d.raw(), 7, 3)) {
        return 1;
    }

    const silex::flint::FmpzPoly& const_d = d;
    if (!equals_poly2(const_d.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpzPolyRef d_ref(d);
    fmpz_poly_set_coeff_si(d_ref.raw(), 0, 29);
    fmpz_poly_set_coeff_si(d_ref.raw(), 1, 10);
    if (!equals_poly2(d.raw(), 29, 10)) {
        return 1;
    }

    const silex::flint::FmpzPolyConstRef d_const_ref(d);
    if (!equals_poly2(d_const_ref.raw(), 29, 10)) {
        return 1;
    }

    fmpz_poly_t external;
    fmpz_poly_init(external);
    {
        silex::flint::FmpzPolyRef external_ref(external);
        fmpz_poly_set_coeff_si(external_ref.raw(), 0, -41);
        fmpz_poly_set_coeff_si(external_ref.raw(), 1, 12);

        silex::flint::FmpzPolyConstRef external_const_ref(external);
        if (!equals_poly2(external_const_ref.raw(), -41, 12)) {
            fmpz_poly_clear(external);
            return 1;
        }
    }

    fmpz_poly_set_coeff_si(external, 0, 43);
    fmpz_poly_set_coeff_si(external, 1, 14);
    if (!equals_poly2(external, 43, 14)) {
        fmpz_poly_clear(external);
        return 1;
    }
    fmpz_poly_clear(external);

    return 0;
}

int test_fmpq_poly() {
    silex::flint::FmpqPoly a;
    fmpq_poly_set_coeff_si(a.raw(), 0, 7);
    fmpq_poly_set_coeff_si(a.raw(), 1, 3);
    if (!equals_qpoly2(a.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpqPoly b;
    fmpq_poly_set_coeff_si(b.raw(), 0, -5);
    fmpq_poly_set_coeff_si(b.raw(), 1, 2);
    if (!equals_qpoly2(b.raw(), -5, 2)) {
        return 1;
    }

    a.swap(b);
    if (!equals_qpoly2(a.raw(), -5, 2) || !equals_qpoly2(b.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpqPoly c(std::move(a));
    if (!equals_qpoly2(c.raw(), -5, 2)) {
        return 1;
    }
    fmpq_poly_set_coeff_si(a.raw(), 0, 11);
    fmpq_poly_set_coeff_si(a.raw(), 1, 4);
    if (!equals_qpoly2(a.raw(), 11, 4)) {
        return 1;
    }

    silex::flint::FmpqPoly d;
    fmpq_poly_set_coeff_si(d.raw(), 0, 19);
    fmpq_poly_set_coeff_si(d.raw(), 1, 6);
    d = std::move(c);
    if (!equals_qpoly2(d.raw(), -5, 2)) {
        return 1;
    }
    fmpq_poly_set_coeff_si(c.raw(), 0, 23);
    fmpq_poly_set_coeff_si(c.raw(), 1, 8);
    if (!equals_qpoly2(c.raw(), 23, 8)) {
        return 1;
    }

    swap(b, d);
    if (!equals_qpoly2(b.raw(), -5, 2) || !equals_qpoly2(d.raw(), 7, 3)) {
        return 1;
    }

    const silex::flint::FmpqPoly& const_d = d;
    if (!equals_qpoly2(const_d.raw(), 7, 3)) {
        return 1;
    }

    silex::flint::FmpqPolyRef d_ref(d);
    fmpq_poly_set_coeff_si(d_ref.raw(), 0, 29);
    fmpq_poly_set_coeff_si(d_ref.raw(), 1, 10);
    if (!equals_qpoly2(d.raw(), 29, 10)) {
        return 1;
    }

    const silex::flint::FmpqPolyConstRef d_const_ref(d);
    if (!equals_qpoly2(d_const_ref.raw(), 29, 10)) {
        return 1;
    }

    fmpq_poly_t external;
    fmpq_poly_init(external);
    {
        silex::flint::FmpqPolyRef external_ref(external);
        fmpq_poly_set_coeff_si(external_ref.raw(), 0, -41);
        fmpq_poly_set_coeff_si(external_ref.raw(), 1, 12);

        silex::flint::FmpqPolyConstRef external_const_ref(external);
        if (!equals_qpoly2(external_const_ref.raw(), -41, 12)) {
            fmpq_poly_clear(external);
            return 1;
        }
    }

    fmpq_poly_set_coeff_si(external, 0, 43);
    fmpq_poly_set_coeff_si(external, 1, 14);
    if (!equals_qpoly2(external, 43, 14)) {
        fmpq_poly_clear(external);
        return 1;
    }
    fmpq_poly_clear(external);

    return 0;
}

int test_fmpz_mat() {
    silex::flint::FmpzMat a(2, 2);
    fmpz_set_si(fmpz_mat_entry(a.raw(), 0, 0), 1);
    fmpz_set_si(fmpz_mat_entry(a.raw(), 0, 1), 2);
    fmpz_set_si(fmpz_mat_entry(a.raw(), 1, 0), 3);
    fmpz_set_si(fmpz_mat_entry(a.raw(), 1, 1), 4);
    if (!equals_mat2(a.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::FmpzMat b(2, 2);
    fmpz_set_si(fmpz_mat_entry(b.raw(), 0, 0), -1);
    fmpz_set_si(fmpz_mat_entry(b.raw(), 0, 1), -2);
    fmpz_set_si(fmpz_mat_entry(b.raw(), 1, 0), -3);
    fmpz_set_si(fmpz_mat_entry(b.raw(), 1, 1), -4);

    a.swap(b);
    if (!equals_mat2(a.raw(), -1, -2, -3, -4) || !equals_mat2(b.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::FmpzMat c(std::move(a));
    if (!equals_mat2(c.raw(), -1, -2, -3, -4)) {
        return 1;
    }

    silex::flint::FmpzMat d(1, 1);
    d = std::move(c);
    if (!equals_mat2(d.raw(), -1, -2, -3, -4)) {
        return 1;
    }

    silex::flint::FmpzMatRef d_ref(d);
    fmpz_set_si(fmpz_mat_entry(d_ref.raw(), 0, 0), 9);
    if (fmpz_equal_si(fmpz_mat_entry(d.raw(), 0, 0), 9) == 0) {
        return 1;
    }

    const silex::flint::FmpzMatConstRef d_const_ref(d);
    if (fmpz_equal_si(fmpz_mat_entry(d_const_ref.raw(), 0, 0), 9) == 0) {
        return 1;
    }
    silex::flint::fmpz_one(silex::flint::fmpz_mat_entry(d, 0, 0));
    if (fmpz_equal_si(silex::flint::fmpz_mat_entry(d_const_ref, 0, 0).raw(), 1) == 0) {
        return 1;
    }

    fmpz_mat_t external;
    fmpz_mat_init(external, 2, 2);
    {
        silex::flint::FmpzMatRef external_ref(external);
        fmpz_set_si(silex::flint::fmpz_mat_entry(external_ref, 1, 1).raw(), 41);
        if (fmpz_equal_si(silex::flint::fmpz_mat_entry(external_ref, 1, 1).raw(), 41) ==
                0) {
            fmpz_mat_clear(external);
            return 1;
        }
    }
    fmpz_set_si(fmpz_mat_entry(external, 0, 0), 43);
    if (fmpz_equal_si(fmpz_mat_entry(external, 0, 0), 43) == 0) {
        fmpz_mat_clear(external);
        return 1;
    }
    fmpz_mat_clear(external);

    silex::flint::FmpzMat window_source(3, 3);
    for (slong i = 0; i < 3; ++i) {
        for (slong j = 0; j < 3; ++j) {
            fmpz_set_si(fmpz_mat_entry(window_source.raw(), i, j), 10 * i + j);
        }
    }
    silex::flint::FmpzMat copied_window(2, 2);
    {
        silex::flint::FmpzMatConstWindow center(window_source, 1, 1, 3, 3);
        silex::flint::fmpz_mat_set(copied_window, center.const_ref());
    }
    if (!equals_mat2(copied_window.raw(), 11, 12, 21, 22)) {
        return 1;
    }
    {
        silex::flint::FmpzMatWindow corner(window_source, 0, 0, 2, 2);
        silex::flint::fmpz_mat_zero(corner.ref());
    }
    if (fmpz_equal_si(fmpz_mat_entry(window_source.raw(), 0, 0), 0) == 0 ||
        fmpz_equal_si(fmpz_mat_entry(window_source.raw(), 1, 1), 0) == 0 ||
        fmpz_equal_si(fmpz_mat_entry(window_source.raw(), 2, 2), 22) == 0) {
        return 1;
    }

    return 0;
}

int test_fmpq_mat() {
    silex::flint::FmpqMat a(2, 2);
    fmpq_set_si(fmpq_mat_entry(a.raw(), 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(a.raw(), 0, 1), 2, 1);
    fmpq_set_si(fmpq_mat_entry(a.raw(), 1, 0), 3, 1);
    fmpq_set_si(fmpq_mat_entry(a.raw(), 1, 1), 4, 1);
    if (!equals_qmat2(a.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::FmpqMat b(2, 2);
    fmpq_set_si(fmpq_mat_entry(b.raw(), 0, 0), -1, 1);
    fmpq_set_si(fmpq_mat_entry(b.raw(), 0, 1), -2, 1);
    fmpq_set_si(fmpq_mat_entry(b.raw(), 1, 0), -3, 1);
    fmpq_set_si(fmpq_mat_entry(b.raw(), 1, 1), -4, 1);

    a.swap(b);
    if (!equals_qmat2(a.raw(), -1, -2, -3, -4) ||
        !equals_qmat2(b.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::FmpqMat c(std::move(a));
    if (!equals_qmat2(c.raw(), -1, -2, -3, -4)) {
        return 1;
    }

    silex::flint::FmpqMat d(1, 1);
    d = std::move(c);
    if (!equals_qmat2(d.raw(), -1, -2, -3, -4)) {
        return 1;
    }

    silex::flint::FmpqMatRef d_ref(d);
    fmpq_set_si(fmpq_mat_entry(d_ref.raw(), 0, 0), 9, 1);
    if (fmpq_equal_si(fmpq_mat_entry(d.raw(), 0, 0), 9) == 0) {
        return 1;
    }

    const silex::flint::FmpqMatConstRef d_const_ref(d);
    if (fmpq_equal_si(fmpq_mat_entry(d_const_ref.raw(), 0, 0), 9) == 0) {
        return 1;
    }

    fmpq_mat_t external;
    fmpq_mat_init(external, 2, 2);
    {
        silex::flint::FmpqMatRef external_ref(external);
        fmpq_set_si(fmpq_mat_entry(external_ref.raw(), 1, 1), 41, 1);
        if (fmpq_equal_si(fmpq_mat_entry(external, 1, 1), 41) == 0) {
            fmpq_mat_clear(external);
            return 1;
        }
    }
    fmpq_set_si(fmpq_mat_entry(external, 0, 0), 43, 1);
    if (fmpq_equal_si(fmpq_mat_entry(external, 0, 0), 43) == 0) {
        fmpq_mat_clear(external);
        return 1;
    }
    fmpq_mat_clear(external);

    return 0;
}

int test_arb() {
    silex::flint::Arb a;
    arb_set_si(a.raw(), 7);
    if (arb_equal_si(a.raw(), 7) == 0) {
        return 1;
    }

    silex::flint::Arb b;
    arb_set_si(b.raw(), -3);
    a.swap(b);
    if (arb_equal_si(a.raw(), -3) == 0 || arb_equal_si(b.raw(), 7) == 0) {
        return 1;
    }

    silex::flint::Arb c(std::move(a));
    if (arb_equal_si(c.raw(), -3) == 0) {
        return 1;
    }
    arb_set_si(a.raw(), 11);
    if (arb_equal_si(a.raw(), 11) == 0) {
        return 1;
    }

    silex::flint::Arb d;
    d = std::move(c);
    if (arb_equal_si(d.raw(), -3) == 0) {
        return 1;
    }
    silex::flint::arb_neg(a, d);
    if (arb_equal_si(a.raw(), 3) == 0 || arb_equal_si(d.raw(), -3) == 0) {
        return 1;
    }
    silex::flint::arb_zero(b);
    if (arb_is_zero(b.raw()) == 0) {
        return 1;
    }
    silex::flint::arb_abs(b, d);
    if (arb_equal_si(b.raw(), 3) == 0 || !silex::flint::arb_is_positive(b)) {
        return 1;
    }
    silex::flint::arb_add(a, b, d, 128);
    if (arb_is_zero(a.raw()) == 0) {
        return 1;
    }
    silex::flint::arb_sub(a, b, d, 128);
    if (arb_equal_si(a.raw(), 6) == 0) {
        return 1;
    }
    silex::flint::arb_mul(a, b, b, 128);
    if (arb_equal_si(a.raw(), 9) == 0) {
        return 1;
    }
    silex::flint::arb_sqrt(a, a, 128);
    if (arb_equal_si(a.raw(), 3) == 0) {
        return 1;
    }
    silex::flint::arb_sqr(a, a, 128);
    if (arb_equal_si(a.raw(), 9) == 0) {
        return 1;
    }
    silex::flint::arb_div(a, a, b, 128);
    if (arb_equal_si(a.raw(), 3) == 0) {
        return 1;
    }
    silex::flint::Fmpz factor;
    fmpz_set_si(factor.raw(), 4);
    silex::flint::arb_mul_fmpz(a, b, factor, 128);
    if (arb_equal_si(a.raw(), 12) == 0) {
        return 1;
    }
    silex::flint::arb_set_fmpz(a, factor);
    if (arb_equal_si(a.raw(), 4) == 0) {
        return 1;
    }
    silex::flint::Arf lower_bound;
    silex::flint::Arf upper_bound;
    silex::flint::arb_get_lbound_arf(lower_bound, b, 128);
    silex::flint::arb_get_ubound_arf(upper_bound, b, 128);
    silex::flint::Fmpz bound_int;
    silex::flint::arf_get_fmpz(bound_int, lower_bound, ARF_RND_FLOOR);
    if (!equals_si(bound_int.raw(), 3) || silex::flint::arf_is_nan(upper_bound) ||
        silex::flint::arf_is_inf(upper_bound)) {
        return 1;
    }

    silex::flint::ArbRef d_ref(d);
    arb_set_si(d_ref.raw(), 29);
    if (arb_equal_si(d.raw(), 29) == 0) {
        return 1;
    }

    const silex::flint::ArbConstRef d_const_ref(d);
    if (arb_equal_si(d_const_ref.raw(), 29) == 0) {
        return 1;
    }

    arb_t external;
    arb_init(external);
    {
        silex::flint::ArbRef external_ref(external);
        arb_set_si(external_ref.raw(), -41);

        silex::flint::ArbConstRef external_const_ref(external);
        if (arb_equal_si(external_const_ref.raw(), -41) == 0) {
            arb_clear(external);
            return 1;
        }
    }

    arb_set_si(external, 43);
    if (arb_equal_si(external, 43) == 0) {
        arb_clear(external);
        return 1;
    }
    arb_clear(external);

    return 0;
}

int test_nf() {
    static_assert(!std::is_copy_constructible_v<silex::flint::Nf>);
    static_assert(!std::is_copy_assignable_v<silex::flint::Nf>);
    static_assert(std::is_move_constructible_v<silex::flint::Nf>);
    static_assert(std::is_move_assignable_v<silex::flint::Nf>);

    silex::flint::Nf empty;
    if (empty.is_defined() || empty.raw() != nullptr) {
        return 1;
    }

    silex::flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_si(polynomial.raw(), 0, -2);
    fmpq_poly_set_coeff_si(polynomial.raw(), 2, 1);

    silex::flint::Nf field(polynomial.raw());
    if (!field.is_defined() || fmpq_poly_degree(field.raw()->pol) != 2) {
        return 1;
    }

    silex::flint::Nf moved(std::move(field));
    if (!moved.is_defined() || field.is_defined() || field.raw() != nullptr ||
        fmpq_poly_degree(moved.raw()->pol) != 2) {
        return 1;
    }

    field = std::move(moved);
    if (!field.is_defined() || moved.is_defined() || moved.raw() != nullptr ||
        fmpq_poly_degree(field.raw()->pol) != 2) {
        return 1;
    }

    empty.swap(field);
    if (!empty.is_defined() || field.is_defined() || field.raw() != nullptr ||
        fmpq_poly_degree(empty.raw()->pol) != 2) {
        return 1;
    }

    field.swap(empty);
    if (!field.is_defined() || empty.is_defined() || empty.raw() != nullptr ||
        fmpq_poly_degree(field.raw()->pol) != 2) {
        return 1;
    }

    silex::flint::NfRef field_ref(field);
    if (fmpq_poly_degree(field_ref.raw()->pol) != 2) {
        return 1;
    }

    const silex::flint::NfConstRef field_const_ref(field);
    if (fmpq_poly_degree(field_const_ref.raw()->pol) != 2) {
        return 1;
    }

    nf_t external;
    nf_init(external, polynomial.raw());
    {
        silex::flint::NfConstRef external_const_ref(external);
        if (fmpq_poly_degree(external_const_ref.raw()->pol) != 2) {
            nf_clear(external);
            return 1;
        }
    }
    if (fmpq_poly_degree(external->pol) != 2) {
        nf_clear(external);
        return 1;
    }
    nf_clear(external);

    return 0;
}

int test_arf() {
    silex::flint::Arf a;
    arf_set_si(a.raw(), 7);
    if (arf_equal_si(a.raw(), 7) == 0) {
        return 1;
    }

    silex::flint::Arf b;
    arf_set_si(b.raw(), -3);
    a.swap(b);
    if (arf_equal_si(a.raw(), -3) == 0 || arf_equal_si(b.raw(), 7) == 0) {
        return 1;
    }

    silex::flint::Arf c(std::move(a));
    if (arf_equal_si(c.raw(), -3) == 0) {
        return 1;
    }
    arf_set_si(a.raw(), 11);
    if (arf_equal_si(a.raw(), 11) == 0) {
        return 1;
    }

    silex::flint::ArfRef ref(c);
    arf_set_si(ref.raw(), 29);
    const silex::flint::ArfConstRef const_ref(c);
    if (arf_equal_si(const_ref.raw(), 29) == 0) {
        return 1;
    }
    silex::flint::Arf copied;
    silex::flint::arf_set(copied, const_ref.raw());
    if (silex::flint::arf_is_nan(copied) || silex::flint::arf_is_inf(copied) ||
        silex::flint::arf_get_d(copied, ARF_RND_NEAR) != 29.0) {
        return 1;
    }

    arf_t external;
    arf_init(external);
    {
        silex::flint::ArfRef external_ref(external);
        arf_set_si(external_ref.raw(), -41);
        silex::flint::ArfConstRef external_const_ref(external);
        if (arf_equal_si(external_const_ref.raw(), -41) == 0) {
            arf_clear(external);
            return 1;
        }
    }
    arf_clear(external);

    return 0;
}

int test_acb() {
    silex::flint::Acb a;
    acb_set_si(a.raw(), 7);
    if (acb_equal_si(a.raw(), 7) == 0) {
        return 1;
    }

    silex::flint::Acb b;
    acb_set_si(b.raw(), -3);
    a.swap(b);
    if (acb_equal_si(a.raw(), -3) == 0 || acb_equal_si(b.raw(), 7) == 0) {
        return 1;
    }

    silex::flint::Acb c(std::move(a));
    if (acb_equal_si(c.raw(), -3) == 0) {
        return 1;
    }
    acb_set_si(a.raw(), 11);
    if (acb_equal_si(a.raw(), 11) == 0) {
        return 1;
    }

    silex::flint::AcbRef ref(c);
    acb_set_si(ref.raw(), 29);
    const silex::flint::AcbConstRef const_ref(c);
    if (acb_equal_si(const_ref.raw(), 29) == 0) {
        return 1;
    }

    acb_t external;
    acb_init(external);
    {
        silex::flint::AcbRef external_ref(external);
        acb_set_si(external_ref.raw(), -41);
        silex::flint::AcbConstRef external_const_ref(external);
        if (acb_equal_si(external_const_ref.raw(), -41) == 0) {
            acb_clear(external);
            return 1;
        }
    }
    acb_clear(external);

    return 0;
}

int test_acb_poly() {
    silex::flint::AcbPoly a;
    acb_poly_set_coeff_si(a.raw(), 0, 7);
    acb_poly_set_coeff_si(a.raw(), 1, 3);
    if (acb_poly_length(a.raw()) != 2) {
        return 1;
    }

    silex::flint::AcbPoly b;
    acb_poly_set_coeff_si(b.raw(), 0, -5);
    a.swap(b);
    if (acb_poly_length(a.raw()) != 1 || acb_poly_length(b.raw()) != 2) {
        return 1;
    }

    silex::flint::AcbPoly c(std::move(a));
    if (acb_poly_length(c.raw()) != 1) {
        return 1;
    }
    acb_poly_set_coeff_si(a.raw(), 0, 11);
    if (acb_poly_length(a.raw()) != 1) {
        return 1;
    }

    silex::flint::AcbPolyRef ref(c);
    acb_poly_set_coeff_si(ref.raw(), 2, 29);
    const silex::flint::AcbPolyConstRef const_ref(c);
    if (acb_poly_length(const_ref.raw()) != 3) {
        return 1;
    }

    acb_poly_t external;
    acb_poly_init(external);
    {
        silex::flint::AcbPolyRef external_ref(external);
        acb_poly_set_coeff_si(external_ref.raw(), 0, -41);
        silex::flint::AcbPolyConstRef external_const_ref(external);
        if (acb_poly_length(external_const_ref.raw()) != 1) {
            acb_poly_clear(external);
            return 1;
        }
    }
    acb_poly_clear(external);

    silex::flint::AcbVec points(2);
    silex::flint::AcbVec values(2);
    acb_set_si(points.data() + 0, 1);
    acb_set_si(points.data() + 1, 2);

    silex::flint::AcbPoly polynomial;
    acb_poly_set_coeff_si(polynomial.raw(), 0, 1);
    acb_poly_set_coeff_si(polynomial.raw(), 2, 1);

    silex::flint::AcbPolyEvaluationTree tree;
    if (!tree.build(points.data(), points.length(), 64) ||
        !tree.matches(2, 64) || tree.length() != 2 ||
        tree.precision() != 64 ||
        !tree.evaluate(values.data(),
                       silex::flint::AcbPolyConstRef(polynomial), 64) ||
        acb_equal_si(values.data() + 0, 2) == 0 ||
        acb_equal_si(values.data() + 1, 5) == 0) {
        return 1;
    }
    if (tree.build(nullptr, 2, 64) || !tree.matches(2, 64)) {
        return 1;
    }

    silex::flint::AcbPolyEvaluationTree moved_tree(std::move(tree));
    if (tree.matches(2, 64) || !moved_tree.matches(2, 64) ||
        moved_tree.evaluate(values.data(),
                            silex::flint::AcbPolyConstRef(polynomial), 128)) {
        return 1;
    }

    return 0;
}

int test_acb_vec() {
    silex::flint::AcbVec a(2);
    if (a.length() != 2) {
        return 1;
    }
    acb_set_si(a.data() + 0, 7);
    acb_set_si(a.data() + 1, -3);

    silex::flint::AcbVec b(2);
    b.set_from(silex::flint::AcbVecConstRef(a));
    if (acb_equal_si(b.data() + 0, 7) == 0 ||
        acb_equal_si(b.data() + 1, -3) == 0) {
        return 1;
    }

    silex::flint::AcbVec c(1);
    if (c.set_from(silex::flint::AcbVecConstRef(a))) {
        return 1;
    }

    a.swap(c);
    if (a.length() != 1 || c.length() != 2) {
        return 1;
    }

    silex::flint::AcbVec d(std::move(c));
    if (d.length() != 2 || acb_equal_si(d.data() + 0, 7) == 0) {
        return 1;
    }

    silex::flint::AcbVecRef ref(d);
    acb_set_si(ref.data() + 1, 29);
    const silex::flint::AcbVecConstRef const_ref(d);
    if (const_ref.length() != 2 ||
        acb_equal_si(const_ref.data() + 1, 29) == 0) {
        return 1;
    }

    return 0;
}

int test_acb_mat() {
    silex::flint::AcbMat a(2, 2);
    acb_set_si(acb_mat_entry(a.raw(), 0, 0), 1);
    acb_set_si(acb_mat_entry(a.raw(), 0, 1), 2);
    acb_set_si(acb_mat_entry(a.raw(), 1, 0), 3);
    acb_set_si(acb_mat_entry(a.raw(), 1, 1), 4);
    if (!equals_acb_mat2(a.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::AcbMat b(2, 2);
    acb_set_si(acb_mat_entry(b.raw(), 0, 0), -1);
    acb_set_si(acb_mat_entry(b.raw(), 0, 1), -2);
    acb_set_si(acb_mat_entry(b.raw(), 1, 0), -3);
    acb_set_si(acb_mat_entry(b.raw(), 1, 1), -4);
    a.swap(b);
    if (!equals_acb_mat2(a.raw(), -1, -2, -3, -4) ||
        !equals_acb_mat2(b.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::AcbMat c(std::move(a));
    if (!equals_acb_mat2(c.raw(), -1, -2, -3, -4)) {
        return 1;
    }

    silex::flint::AcbMatRef ref(c);
    acb_set_si(silex::flint::acb_mat_entry_ref(ref, 0, 0).raw(), 29);
    const silex::flint::AcbMatConstRef const_ref(c);
    if (acb_equal_si(
                silex::flint::acb_mat_entry_ref(const_ref, 0, 0).raw(),
                29) == 0) {
        return 1;
    }

    silex::flint::AcbMat identity(2, 2);
    silex::flint::AcbMat rhs(2, 1);
    silex::flint::AcbMat solution(2, 1);
    acb_one(acb_mat_entry(identity.raw(), 0, 0));
    acb_one(acb_mat_entry(identity.raw(), 1, 1));
    acb_set_si(acb_mat_entry(rhs.raw(), 0, 0), 5);
    acb_set_si(acb_mat_entry(rhs.raw(), 1, 0), -7);
    if (silex::flint::acb_mat_solve(solution, identity, rhs, 128) == 0 ||
        acb_equal_si(acb_mat_entry(solution.raw(), 0, 0), 5) == 0 ||
        acb_equal_si(acb_mat_entry(solution.raw(), 1, 0), -7) == 0) {
        return 1;
    }

    return 0;
}

int test_arb_vec() {
    silex::flint::ArbVec a(2);
    if (a.length() != 2) {
        return 1;
    }
    arb_set_si(a.data() + 0, 7);
    arb_set_si(a.data() + 1, -3);

    silex::flint::ArbVec b(2);
    b.set_from(silex::flint::ArbVecConstRef(a));
    if (arb_contains_si(b.data() + 0, 7) == 0 ||
        arb_contains_si(b.data() + 1, -3) == 0) {
        return 1;
    }

    silex::flint::ArbVec c(1);
    if (c.set_from(silex::flint::ArbVecConstRef(a))) {
        return 1;
    }

    a.swap(c);
    if (a.length() != 1 || c.length() != 2) {
        return 1;
    }

    silex::flint::ArbVec d(std::move(c));
    if (d.length() != 2 || arb_contains_si(d.data() + 0, 7) == 0) {
        return 1;
    }

    silex::flint::ArbVecRef ref(d);
    arb_set_si(ref.data() + 1, 29);
    const silex::flint::ArbVecConstRef const_ref(d);
    if (const_ref.length() != 2 ||
        arb_contains_si(const_ref.data() + 1, 29) == 0) {
        return 1;
    }

    return 0;
}

int test_arb_mat() {
    silex::flint::ArbMat a(2, 2);
    arb_set_si(arb_mat_entry(a.raw(), 0, 0), 1);
    arb_set_si(arb_mat_entry(a.raw(), 0, 1), 2);
    arb_set_si(arb_mat_entry(a.raw(), 1, 0), 3);
    arb_set_si(arb_mat_entry(a.raw(), 1, 1), 4);
    if (!equals_arb_mat2(a.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::ArbMat b(2, 2);
    arb_set_si(arb_mat_entry(b.raw(), 0, 0), -1);
    arb_set_si(arb_mat_entry(b.raw(), 0, 1), -2);
    arb_set_si(arb_mat_entry(b.raw(), 1, 0), -3);
    arb_set_si(arb_mat_entry(b.raw(), 1, 1), -4);
    a.swap(b);
    if (!equals_arb_mat2(a.raw(), -1, -2, -3, -4) ||
        !equals_arb_mat2(b.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::ArbMat c(std::move(a));
    if (!equals_arb_mat2(c.raw(), -1, -2, -3, -4)) {
        return 1;
    }

    silex::flint::ArbMatRef ref(c);
    arb_set_si(arb_mat_entry(ref.raw(), 0, 0), 29);
    const silex::flint::ArbMatConstRef const_ref(c);
    if (arb_equal_si(arb_mat_entry(const_ref.raw(), 0, 0), 29) == 0) {
        return 1;
    }
    silex::flint::ArbMat identity(2, 2);
    silex::flint::ArbMat cholesky(2, 2);
    arb_one(arb_mat_entry(identity.raw(), 0, 0));
    arb_one(arb_mat_entry(identity.raw(), 1, 1));
    if (silex::flint::arb_mat_cho(cholesky, identity, 128) == 0) {
        return 1;
    }

    arb_mat_t external;
    arb_mat_init(external, 2, 2);
    {
        silex::flint::ArbMatRef external_ref(external);
        arb_set_si(arb_mat_entry(external_ref.raw(), 1, 1), -41);
        silex::flint::ArbMatConstRef external_const_ref(external);
        if (arb_equal_si(arb_mat_entry(external_const_ref.raw(), 1, 1), -41) == 0) {
            arb_mat_clear(external);
            return 1;
        }
    }
    arb_mat_clear(external);

    return 0;
}

int test_nmod_mat() {
    silex::flint::NmodMat a(2, 2, 17);
    silex::flint::nmod_mat_set_entry(a, 0, 0, 1);
    silex::flint::nmod_mat_set_entry(a, 0, 1, 2);
    silex::flint::nmod_mat_set_entry(a, 1, 0, 3);
    silex::flint::nmod_mat_set_entry(a, 1, 1, 4);
    if (!equals_nmod_mat2(a.raw(), 1, 2, 3, 4)) {
        return 1;
    }

    silex::flint::NmodMat b(2, 2, 17);
    silex::flint::nmod_mat_set_entry(b, 0, 0, 5);
    a.swap(b);
    if (silex::flint::nmod_mat_get_entry(a, 0, 0) != 5 ||
        silex::flint::nmod_mat_get_entry(b, 0, 0) != 1) {
        return 1;
    }

    silex::flint::NmodMat c(std::move(a));
    if (silex::flint::nmod_mat_get_entry(c, 0, 0) != 5) {
        return 1;
    }

    silex::flint::NmodMatRef ref(c);
    silex::flint::nmod_mat_set_entry(ref, 1, 1, 9);
    const silex::flint::NmodMatConstRef const_ref(c);
    if (silex::flint::nmod_mat_get_entry(const_ref, 1, 1) != 9) {
        return 1;
    }

    nmod_mat_t external;
    nmod_mat_init(external, 2, 2, 17);
    {
        silex::flint::NmodMatRef external_ref(external);
        silex::flint::nmod_mat_set_entry(external_ref, 1, 1, 11);
        silex::flint::NmodMatConstRef external_const_ref(external);
        if (silex::flint::nmod_mat_get_entry(external_const_ref, 1, 1) != 11) {
            nmod_mat_clear(external);
            return 1;
        }
    }
    nmod_mat_clear(external);

    return 0;
}

int test_nmod_poly_factor() {
    silex::flint::NmodPoly polynomial(5);
    nmod_poly_set_coeff_ui(polynomial.raw(), 0, 2);
    nmod_poly_set_coeff_ui(polynomial.raw(), 1, 2);
    nmod_poly_set_coeff_ui(polynomial.raw(), 2, 1);
    nmod_poly_set_coeff_ui(polynomial.raw(), 3, 1);
    if (nmod_poly_degree(polynomial.raw()) != 3) {
        return 1;
    }

    silex::flint::NmodPoly other(5);
    nmod_poly_set_coeff_ui(other.raw(), 1, 1);
    polynomial.swap(other);
    if (nmod_poly_degree(polynomial.raw()) != 1 ||
        nmod_poly_degree(other.raw()) != 3) {
        return 1;
    }

    silex::flint::NmodPoly moved(std::move(other));
    if (nmod_poly_degree(moved.raw()) != 3) {
        return 1;
    }

    silex::flint::NmodPolyRef poly_ref(moved);
    nmod_poly_set_coeff_ui(poly_ref.raw(), 4, 0);
    const silex::flint::NmodPolyConstRef poly_const_ref(moved);
    if (nmod_poly_degree(poly_const_ref.raw()) != 3) {
        return 1;
    }

    silex::flint::NmodPolyFactor factorization;
    slong degrees_storage[3] = {};
    slong* degrees = degrees_storage;
    nmod_poly_factor_distinct_deg(factorization.raw(), moved.raw(),
                                  &degrees);
    if (factorization.raw()->num != 2) {
        return 1;
    }

    bool saw_degree_one = false;
    bool saw_degree_two = false;
    const silex::flint::NmodPolyFactorConstRef factorization_const_ref(
            factorization);
    for (slong i = 0; i < factorization_const_ref.raw()->num; ++i) {
        const slong product_degree =
                nmod_poly_degree(factorization_const_ref.raw()->p + i);
        if (degrees_storage[i] == 1 && product_degree == 1) {
            saw_degree_one = true;
        }
        if (degrees_storage[i] == 2 && product_degree == 2) {
            saw_degree_two = true;
        }
    }
    if (!saw_degree_one || !saw_degree_two) {
        return 1;
    }

    silex::flint::NmodPolyFactor moved_factor(std::move(factorization));
    if (moved_factor.raw()->num != 2) {
        return 1;
    }

    silex::flint::NmodPolyFactor other_factor;
    moved_factor.swap(other_factor);
    if (moved_factor.raw()->num != 0 || other_factor.raw()->num != 2) {
        return 1;
    }

    silex::flint::NmodPolyFactorRef factor_ref(other_factor);
    if (factor_ref.raw()->num != 2) {
        return 1;
    }

    return 0;
}

int test_fmpz_factor() {
    silex::flint::Fmpz n;
    fmpz_set_ui(n.raw(), 12);

    silex::flint::FmpzFactor factor;
    fmpz_factor(factor.raw(), n.raw());
    if (factor.raw()->num != 2) {
        return 1;
    }

    silex::flint::FmpzFactorRef ref(factor);
    if (ref.raw()->num != 2) {
        return 1;
    }

    const silex::flint::FmpzFactorConstRef const_ref(factor);
    if (const_ref.raw()->num != 2) {
        return 1;
    }

    fmpz_factor_t external;
    fmpz_factor_init(external);
    {
        silex::flint::FmpzFactorRef external_ref(external);
        fmpz_factor(external_ref.raw(), n.raw());
        silex::flint::FmpzFactorConstRef external_const_ref(external);
        if (external_const_ref.raw()->num != 2) {
            fmpz_factor_clear(external);
            return 1;
        }
    }
    fmpz_factor_clear(external);

    return 0;
}

int test_fmpz_poly_factor() {
    silex::flint::FmpzPoly polynomial;
    fmpz_poly_set_coeff_si(polynomial.raw(), 0, -1);
    fmpz_poly_set_coeff_si(polynomial.raw(), 2, 1);

    silex::flint::FmpzPolyFactor factorization;
    silex::flint::fmpz_poly_factor(
            silex::flint::FmpzPolyFactorRef(factorization),
            silex::flint::FmpzPolyConstRef(polynomial));

    const silex::flint::FmpzPolyFactorConstRef factorization_ref(
            factorization);
    if (silex::flint::fmpz_poly_factor_num(factorization_ref) != 2) {
        return 1;
    }
    if (silex::flint::fmpz_poly_factor_exp(factorization_ref, 0) != 1 ||
        silex::flint::fmpz_poly_factor_exp(factorization_ref, 1) != 1) {
        return 1;
    }
    if (silex::flint::fmpz_poly_factor_poly_degree(
                factorization_ref, 0) != 1 ||
        silex::flint::fmpz_poly_factor_poly_degree(
                factorization_ref, 1) != 1) {
        return 1;
    }

    fmpz_poly_factor_t external;
    fmpz_poly_factor_init(external);
    {
        silex::flint::FmpzPolyFactorRef external_ref(external);
        silex::flint::fmpz_poly_factor(
                external_ref,
                silex::flint::FmpzPolyConstRef(polynomial));
        const silex::flint::FmpzPolyFactorConstRef external_const_ref(
                external);
        if (silex::flint::fmpz_poly_factor_num(external_const_ref) != 2) {
            fmpz_poly_factor_clear(external);
            return 1;
        }
    }
    fmpz_poly_factor_clear(external);

    return 0;
}

int test_fmpz_vec() {
    silex::flint::FmpzVec vector(2);
    if (vector.length() != 2 || vector.data() == nullptr) {
        return 1;
    }
    fmpz_set_si(vector.data() + 0, 3);
    fmpz_set_si(vector.data() + 1, 5);

    silex::flint::FmpzVec other(1);
    fmpz_set_si(other.data(), -7);
    vector.swap(other);
    if (vector.length() != 1 || !equals_si(vector.data(), -7) ||
        other.length() != 2 || !equals_si(other.data() + 1, 5)) {
        return 1;
    }

    silex::flint::FmpzVec moved(std::move(other));
    if (moved.length() != 2 || !equals_si(moved.data(), 3)) {
        return 1;
    }

    silex::flint::FmpzVecRef ref(moved);
    fmpz_set_si(ref.data() + 0, 11);
    const silex::flint::FmpzVecConstRef const_ref(moved);
    if (const_ref.length() != 2 || !equals_si(const_ref.data(), 11)) {
        return 1;
    }

    fmpz* external = _fmpz_vec_init(2);
    {
        silex::flint::FmpzVecRef external_ref(external, 2);
        fmpz_set_si(external_ref.data() + 1, 13);
        silex::flint::FmpzVecConstRef external_const_ref(external, 2);
        if (external_const_ref.length() != 2 ||
            !equals_si(external_const_ref.data() + 1, 13)) {
            _fmpz_vec_clear(external, 2);
            return 1;
        }
    }
    _fmpz_vec_clear(external, 2);

    return 0;
}

int test_fmpz_mod() {
    silex::flint::FmpzModCtx ctx(17UL);
    silex::flint::FmpzModCtxRef ctx_ref(ctx);
    const silex::flint::FmpzModCtxConstRef ctx_const_ref(ctx);

    silex::flint::Fmpz coeff;
    fmpz_set_ui(coeff.raw(), 5);

    silex::flint::FmpzModPoly poly(ctx);
    fmpz_mod_poly_set_coeff_fmpz(poly.raw(), 1, coeff.raw(), ctx.raw());
    if (fmpz_mod_poly_length(poly.raw(), ctx_const_ref.raw()) != 2 ||
        poly.context() != ctx_const_ref.raw()) {
        return 1;
    }

    silex::flint::Fmpz out;
    fmpz_mod_poly_get_coeff_fmpz(out.raw(), poly.raw(), 1, ctx.raw());
    if (!equals_si(out.raw(), 5)) {
        return 1;
    }

    silex::flint::FmpzModPolyRef poly_ref(poly);
    const silex::flint::FmpzModPolyConstRef poly_const_ref(poly);
    if (fmpz_mod_poly_length(poly_ref.raw(), ctx.raw()) != 2 ||
        fmpz_mod_poly_length(poly_const_ref.raw(), ctx.raw()) != 2) {
        return 1;
    }

    {
        silex::flint::FmpzModCtx movable_ctx(19UL);
        const fmpz_mod_ctx_struct* const moved_ctx_address = movable_ctx.raw();
        silex::flint::FmpzModCtx moved_ctx(std::move(movable_ctx));
        if (moved_ctx.raw() != moved_ctx_address) {
            return 1;
        }

        silex::flint::FmpzModPoly moved_poly(moved_ctx);
        fmpz_mod_poly_set_coeff_fmpz(
                moved_poly.raw(), 1, coeff.raw(), moved_ctx.raw());
        silex::flint::FmpzModPoly assigned_poly;
        if (assigned_poly.is_initialized()) {
            return 1;
        }
        assigned_poly = std::move(moved_poly);
        if (!assigned_poly.is_initialized() ||
            assigned_poly.context() != moved_ctx.raw() ||
            fmpz_mod_poly_length(assigned_poly.raw(), moved_ctx.raw()) != 2) {
            return 1;
        }
    }

    silex::flint::FmpzModMat mat(2, 2, ctx);
    fmpz_mod_mat_set_entry(mat.raw(), 0, 1, coeff.raw(), ctx.raw());
    if (fmpz_mod_mat_nrows(mat.raw(), ctx.raw()) != 2 ||
        fmpz_mod_mat_ncols(mat.raw(), ctx.raw()) != 2 ||
        mat.context() != ctx_const_ref.raw()) {
        return 1;
    }

    fmpz_zero(out.raw());
    fmpz_mod_mat_get_entry(out.raw(), mat.raw(), 0, 1, ctx.raw());
    if (!equals_si(out.raw(), 5)) {
        return 1;
    }

    silex::flint::FmpzModMatRef mat_ref(mat);
    const silex::flint::FmpzModMatConstRef mat_const_ref(mat);
    fmpz_mod_mat_get_entry(out.raw(), mat_ref.raw(), 0, 1, ctx.raw());
    if (!equals_si(out.raw(), 5)) {
        return 1;
    }
    fmpz_mod_mat_get_entry(out.raw(), mat_const_ref.raw(), 0, 1, ctx.raw());
    if (!equals_si(out.raw(), 5)) {
        return 1;
    }

    return 0;
}

int test_fq() {
    static_assert(!std::is_copy_constructible_v<silex::flint::FqCtx>);
    static_assert(!std::is_copy_assignable_v<silex::flint::FqCtx>);
    static_assert(std::is_move_constructible_v<silex::flint::FqCtx>);
    static_assert(std::is_move_assignable_v<silex::flint::FqCtx>);
    static_assert(!std::is_copy_constructible_v<silex::flint::Fq>);
    static_assert(!std::is_copy_assignable_v<silex::flint::Fq>);
    static_assert(std::is_move_constructible_v<silex::flint::Fq>);
    static_assert(std::is_move_assignable_v<silex::flint::Fq>);

    silex::flint::FmpzModCtx prime_ctx(5UL);
    silex::flint::FmpzModPoly modulus(prime_ctx);
    fmpz_mod_poly_set_coeff_ui(modulus.raw(), 2, 1, prime_ctx.raw());
    fmpz_mod_poly_set_coeff_ui(modulus.raw(), 0, 2, prime_ctx.raw());

    silex::flint::FqCtx field(modulus, prime_ctx, "z");
    if (!field.is_initialized()) {
        return 1;
    }

    silex::flint::FqCtx moved_field(std::move(field));
    if (!moved_field.is_initialized() || field.is_initialized()) {
        return 1;
    }

    silex::flint::Fq value(moved_field);
    silex::flint::Fq square(moved_field);
    silex::flint::Fq root(moved_field);
    silex::flint::Fq check(moved_field);
    if (!value.is_initialized() || !square.is_initialized() ||
        !root.is_initialized() || !check.is_initialized()) {
        return 1;
    }

    fq_gen(value.raw(), moved_field.raw());
    fq_sqr(square.raw(), value.raw(), moved_field.raw());
    if (fq_sqrt(root.raw(), square.raw(), moved_field.raw()) == 0) {
        return 1;
    }

    silex::flint::Fq moved_root(std::move(root));
    if (!moved_root.is_initialized() || root.is_initialized()) {
        return 1;
    }

    fq_sqr(check.raw(), moved_root.raw(), moved_field.raw());
    return fq_equal(check.raw(), square.raw(), moved_field.raw()) == 0 ? 1 : 0;
}

int test_fmpz_mod_poly_factor() {
    silex::flint::FmpzModCtx ctx(17UL);
    silex::flint::FmpzModPoly poly(ctx);

    silex::flint::Fmpz coeff;
    fmpz_set_ui(coeff.raw(), 1);
    fmpz_mod_poly_set_coeff_fmpz(poly.raw(), 1, coeff.raw(), ctx.raw());
    fmpz_mod_poly_set_coeff_fmpz(poly.raw(), 0, coeff.raw(), ctx.raw());

    silex::flint::FmpzModPolyFactor factor(ctx);
    fmpz_mod_poly_factor_insert(factor.raw(), poly.raw(), 1, ctx.raw());
    if (factor.raw()->num != 1 || factor.context() != ctx.raw()) {
        return 1;
    }

    silex::flint::FmpzModPolyFactor other(ctx);
    factor.swap(other);
    if (factor.raw()->num != 0 || other.raw()->num != 1) {
        return 1;
    }

    silex::flint::FmpzModPolyFactor moved(std::move(other));
    if (moved.raw()->num != 1) {
        return 1;
    }

    silex::flint::FmpzModPolyFactorRef ref(moved);
    const silex::flint::FmpzModPolyFactorConstRef const_ref(moved);
    if (ref.raw()->num != 1 || const_ref.raw()->num != 1) {
        return 1;
    }

    fmpz_mod_poly_factor_t external;
    fmpz_mod_poly_factor_init(external, ctx.raw());
    {
        silex::flint::FmpzModPolyFactorRef external_ref(external);
        fmpz_mod_poly_factor_insert(external_ref.raw(), poly.raw(), 1, ctx.raw());
        silex::flint::FmpzModPolyFactorConstRef external_const_ref(external);
        if (external_const_ref.raw()->num != 1) {
            fmpz_mod_poly_factor_clear(external, ctx.raw());
            return 1;
        }
    }
    fmpz_mod_poly_factor_clear(external, ctx.raw());

    return 0;
}

int test_fmpz_lll() {
    silex::flint::FmpzLll config;
    if (config.raw()->delta <= config.raw()->eta) {
        return 1;
    }

    silex::flint::FmpzLll custom(0.95, 0.51, Z_BASIS, APPROX);
    if (custom.raw()->rt != Z_BASIS || custom.raw()->gt != APPROX) {
        return 1;
    }

    config.swap(custom);
    if (config.raw()->rt != Z_BASIS || custom.raw()->delta <= custom.raw()->eta) {
        return 1;
    }

    silex::flint::FmpzLll moved(std::move(config));
    if (moved.raw()->rt != Z_BASIS) {
        return 1;
    }

    silex::flint::FmpzLllRef ref(moved);
    ref.raw()->delta = 0.75;
    const silex::flint::FmpzLllConstRef const_ref(moved);
    if (const_ref.raw()->delta != 0.75) {
        return 1;
    }

    fmpz_lll_t external;
    fmpz_lll_context_init_default(external);
    {
        silex::flint::FmpzLllRef external_ref(external);
        external_ref.raw()->eta = 0.52;
        silex::flint::FmpzLllConstRef external_const_ref(external);
        if (external_const_ref.raw()->eta != 0.52) {
            return 1;
        }
    }

    return 0;
}

int test_nf_elem() {
    static_assert(!std::is_copy_constructible_v<silex::flint::NfElem>);
    static_assert(!std::is_copy_assignable_v<silex::flint::NfElem>);
    static_assert(std::is_move_constructible_v<silex::flint::NfElem>);
    static_assert(std::is_move_assignable_v<silex::flint::NfElem>);
    static_assert(!std::is_copy_constructible_v<silex::flint::NfElemVec>);
    static_assert(!std::is_copy_assignable_v<silex::flint::NfElemVec>);
    static_assert(!std::is_move_constructible_v<silex::flint::NfElemVec>);
    static_assert(!std::is_move_assignable_v<silex::flint::NfElemVec>);

    silex::flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_si(polynomial.raw(), 0, -2);
    fmpq_poly_set_coeff_si(polynomial.raw(), 2, 1);

    silex::flint::Nf field(polynomial.raw());
    silex::flint::NfElem empty_element;
    if (empty_element.is_defined() || empty_element.raw() != nullptr ||
        empty_element.field() != nullptr) {
        return 1;
    }

    silex::flint::NfElem value(field);
    if (!value.is_defined()) {
        return 1;
    }
    nf_elem_set_si(value.raw(), 3, field.raw());
    if (nf_elem_equal_si(value.raw(), 3, field.raw()) == 0 ||
        value.field() != field.raw()) {
        return 1;
    }

    silex::flint::NfElem moved(std::move(value));
    if (!moved.is_defined() || value.is_defined() || value.raw() != nullptr ||
        nf_elem_equal_si(moved.raw(), 3, field.raw()) == 0) {
        return 1;
    }

    value = std::move(moved);
    if (!value.is_defined() || moved.is_defined() || moved.raw() != nullptr ||
        nf_elem_equal_si(value.raw(), 3, field.raw()) == 0) {
        return 1;
    }

    empty_element.swap(value);
    if (!empty_element.is_defined() || value.is_defined() ||
        value.raw() != nullptr ||
        nf_elem_equal_si(empty_element.raw(), 3, field.raw()) == 0) {
        return 1;
    }

    value.swap(empty_element);
    if (!value.is_defined() || empty_element.is_defined() ||
        empty_element.raw() != nullptr ||
        nf_elem_equal_si(value.raw(), 3, field.raw()) == 0) {
        return 1;
    }

    silex::flint::NfElemRef ref(value);
    nf_elem_set_si(ref.raw(), -5, field.raw());
    const silex::flint::NfElemConstRef const_ref(value);
    if (nf_elem_equal_si(const_ref.raw(), -5, field.raw()) == 0) {
        return 1;
    }

    nf_elem_t external;
    nf_elem_init(external, field.raw());
    {
        silex::flint::NfElemRef external_ref(external);
        nf_elem_set_si(external_ref.raw(), 7, field.raw());
        silex::flint::NfElemConstRef external_const_ref(external);
        if (nf_elem_equal_si(external_const_ref.raw(), 7, field.raw()) == 0) {
            nf_elem_clear(external, field.raw());
            return 1;
        }
    }
    nf_elem_clear(external, field.raw());

    silex::flint::NfElemVec vector(2, field);
    if (vector.length() != 2 || vector.data() == nullptr ||
        vector.field() != field.raw()) {
        return 1;
    }
    nf_elem_set_si(vector.data() + 0, 11, field.raw());
    nf_elem_set_si(vector.data() + 1, -13, field.raw());
    if (nf_elem_equal_si(vector.data() + 0, 11, field.raw()) == 0 ||
        nf_elem_equal_si(vector.data() + 1, -13, field.raw()) == 0) {
        return 1;
    }

    silex::flint::NfElemVecRef vector_ref(vector);
    nf_elem_set_si(vector_ref.data() + 1, 17, field.raw());
    const silex::flint::NfElemVecConstRef vector_const_ref(vector);
    if (vector_const_ref.length() != 2 ||
        nf_elem_equal_si(vector_const_ref.data() + 1, 17,
                         field.raw()) == 0) {
        return 1;
    }

    silex::flint::NfElemVec empty(0, field);
    if (empty.length() != 0 || empty.data() != nullptr) {
        return 1;
    }

    return 0;
}

int test_dirichlet() {
    silex::flint::DirichletGroup group(5);
    if (dirichlet_group_size(group.raw()) != 4) {
        return 1;
    }

    silex::flint::DirichletGroupRef group_ref(group);
    const silex::flint::DirichletGroupConstRef group_const_ref(group);
    if (dirichlet_group_size(group_ref.raw()) != 4 ||
        dirichlet_group_size(group_const_ref.raw()) != 4) {
        return 1;
    }

    silex::flint::DirichletChar principal(group);
    dirichlet_char_one(principal.raw(), group.raw());
    if (dirichlet_char_is_principal(group.raw(), principal.raw()) == 0 ||
        principal.group() != group.raw()) {
        return 1;
    }

    silex::flint::DirichletCharRef char_ref(principal);
    const silex::flint::DirichletCharConstRef char_const_ref(principal);
    if (dirichlet_char_is_principal(group.raw(), char_ref.raw()) == 0 ||
        dirichlet_char_is_principal(group.raw(), char_const_ref.raw()) == 0) {
        return 1;
    }

    return 0;
}

int test_qfb() {
    static_assert(!std::is_copy_constructible_v<silex::flint::Qfb>);
    static_assert(!std::is_copy_assignable_v<silex::flint::Qfb>);
    static_assert(std::is_nothrow_move_constructible_v<silex::flint::Qfb>);
    static_assert(std::is_nothrow_move_assignable_v<silex::flint::Qfb>);

    silex::flint::Fmpz discriminant;
    silex::flint::Fmpz prime;
    silex::flint::fmpz_set_si(silex::flint::FmpzRef(discriminant), -47);
    silex::flint::fmpz_set_ui(silex::flint::FmpzRef(prime), 2);

    silex::flint::Qfb form;
    silex::flint::qfb_prime_form(
            silex::flint::QfbRef(form),
            silex::flint::FmpzConstRef(discriminant),
            silex::flint::FmpzConstRef(prime));
    if (!silex::flint::fmpz_equal_si(
                silex::flint::qfb_a(silex::flint::QfbConstRef(form)), 2)) {
        return 1;
    }

    silex::flint::Qfb square;
    silex::flint::qfb_pow_ui(
            silex::flint::QfbRef(square),
            silex::flint::QfbConstRef(form),
            silex::flint::FmpzConstRef(discriminant), 2);

    silex::flint::Qfb moved(std::move(square));
    silex::flint::Qfb principal;
    silex::flint::qfb_principal_form(
            silex::flint::QfbRef(principal),
            silex::flint::FmpzConstRef(discriminant));
    principal.swap(moved);
    if (!silex::flint::qfb_is_principal_form(
                silex::flint::QfbConstRef(moved),
                silex::flint::FmpzConstRef(discriminant))) {
        return 1;
    }

    qfb_t external;
    qfb_init(external);
    {
        silex::flint::QfbRef external_ref(external);
        silex::flint::qfb_set(external_ref,
                              silex::flint::QfbConstRef(moved));
    }
    const bool external_is_principal =
            qfb_is_principal_form(external, discriminant.raw()) != 0;
    qfb_clear(external);
    return external_is_principal ? 0 : 1;
}

}  // namespace

int main() {
    return test_fmpz() != 0 || test_fmpq() != 0 || test_fmpz_poly() != 0 ||
                   test_fmpq_poly() != 0 || test_fmpz_mat() != 0 ||
                   test_fmpq_mat() != 0 || test_arb() != 0 || test_nf() != 0 ||
                   test_arf() != 0 || test_acb() != 0 || test_acb_poly() != 0 ||
                   test_acb_vec() != 0 || test_acb_mat() != 0 ||
                   test_arb_vec() != 0 ||
                   test_arb_mat() != 0 ||
                   test_nmod_mat() != 0 ||
                   test_nmod_poly_factor() != 0 ||
                   test_fmpz_factor() != 0 ||
                   test_fmpz_poly_factor() != 0 ||
                   test_fmpz_vec() != 0 ||
                   test_fmpz_mod() != 0 ||
                   test_fq() != 0 ||
                   test_fmpz_mod_poly_factor() != 0 || test_fmpz_lll() != 0 ||
                   test_nf_elem() != 0 || test_dirichlet() != 0 ||
                   test_qfb() != 0
               ? 1
               : 0;
}
