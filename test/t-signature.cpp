#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>
#include <silex/signature.hpp>

#include "test_support.hpp"

#include <cassert>

namespace {
namespace sflint = silex::flint;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x2_minus(sflint::FmpqPoly& polynomial, slong a) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -a);
}

void poly_x3_plus_ax_plus_b(sflint::FmpqPoly& polynomial,
                            slong a,
                            slong b) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, a);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, b);
}

silex::NumberField field_by_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// White-box helper: returns a defined field whose FLINT context has been
// replaced by one for the non-squarefree polynomial x^2.  NumberField
// construction rejects x^2, so this is the only way to keep the internal
// non-squarefree failure paths covered.  The field starts as the generic
// (non-quadratic-backend) field defined by x^2 - x - 1 so that no quadratic
// backend data disagrees with the installed polynomial.
silex::NumberField nonsquarefree_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    assert(field.backend_kind() == silex::NumberFieldBackendKind::generic);

    sflint::FmpqPoly square;
    sflint::fmpq_poly_set_coeff_si(square, 2, 1);
    nf_struct* raw = field.raw_flint_field();
    assert(raw != nullptr);
    nf_clear(raw);
    nf_init(raw, square.raw());
    assert(field.degree() == 2);
    return field;
}

void check_square_rejected() noexcept {
    sflint::FmpqPoly square;
    sflint::fmpq_poly_set_coeff_si(square, 2, 1);
    assert(!silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(square)).is_defined());
}

}  // namespace

int main() {
    silex::Signature sig;
    assert(sig.r1() == 0);
    assert(sig.r2() == 0);
    assert(sig.degree() == 0);

    sflint::FmpqPoly polynomial;
    silex::NumberField field;

    poly_x(polynomial);
    field = field_by_polynomial(polynomial);
    assert(signature(sig, field));
    assert(sig.r1() == 1);
    assert(sig.r2() == 0);
    assert(sig.degree() == 1);
    assert(is_totally_real(field));
    assert(!is_totally_complex(field));

    silex::Signature copy;
    copy.set(sig.r1(), sig.r2());
    assert(copy.r1() == 1);
    assert(copy.r2() == 0);

    poly_x2_minus(polynomial, 2);
    field = field_by_polynomial(polynomial);
    assert(sig.compute(field));
    assert(sig.r1() == 2);
    assert(sig.r2() == 0);
    assert(is_totally_real(field));
    assert(!is_totally_complex(field));

    poly_x2_minus(polynomial, -1);
    field = field_by_polynomial(polynomial);
    assert(signature(sig, field));
    assert(sig.r1() == 0);
    assert(sig.r2() == 1);
    assert(!is_totally_real(field));
    assert(is_totally_complex(field));

    poly_x3_plus_ax_plus_b(polynomial, -3, 1);
    field = field_by_polynomial(polynomial);
    assert(signature(sig, field));
    assert(sig.r1() == 3);
    assert(sig.r2() == 0);
    assert(is_totally_real(field));
    assert(!is_totally_complex(field));

    poly_x3_plus_ax_plus_b(polynomial, 0, -2);
    field = field_by_polynomial(polynomial);
    assert(signature(sig, field));
    assert(sig.r1() == 1);
    assert(sig.r2() == 1);
    assert(!is_totally_real(field));
    assert(!is_totally_complex(field));

    field = silex::test::quadratic_field(5);
    assert(signature(sig, field));
    assert(sig.r1() == 2);
    assert(sig.r2() == 0);
    assert(is_totally_real(field));

    field = silex::test::quadratic_field(-5);
    assert(signature(sig, field));
    assert(sig.r1() == 0);
    assert(sig.r2() == 1);
    assert(is_totally_complex(field));

    check_square_rejected();
    silex::NumberField undefined;
    sig.set(2, 0);
    assert(!signature(sig, undefined));
    assert(sig.r1() == 2);
    assert(sig.r2() == 0);
    assert(!sig.compute(undefined));
    assert(sig.r1() == 2);
    assert(sig.r2() == 0);
    assert(!is_totally_real(undefined));
    assert(!is_totally_complex(undefined));

    sig.set(2, 0);
    field = nonsquarefree_field();
    assert(!signature(sig, field));
    assert(sig.r1() == 2);
    assert(sig.r2() == 0);
    assert(!is_totally_real(field));
    assert(!is_totally_complex(field));

    return 0;
}
