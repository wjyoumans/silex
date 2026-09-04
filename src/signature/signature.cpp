#include <silex/signature.hpp>

#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>

#include <silex/flint/fmpz_poly.hpp>

namespace silex {

void Signature::set(slong real_embeddings, slong complex_pairs) noexcept {
    r1_ = real_embeddings;
    r2_ = complex_pairs;
}

bool Signature::compute(const NumberField& field) noexcept {
    return signature(*this, field);
}

bool signature(Signature& out, const NumberField& field) noexcept {
    if (!field.is_defined()) {
        return false;
    }

    flint::Fmpz radicand;
    if (field.backend_kind() == NumberFieldBackendKind::quadratic &&
        field.quadratic_radicand(flint::FmpzRef(radicand))) {
        if (fmpz_sgn(radicand.raw()) > 0) {
            out.set(2, 0);
        } else {
            out.set(0, 1);
        }
        return true;
    }

    const nf_struct* raw_field = field.raw_flint_field();
    if (raw_field == nullptr || fmpq_poly_is_squarefree(raw_field->pol) == 0) {
        return false;
    }

    flint::FmpzPoly numerator;
    slong real_embeddings = 0;
    slong complex_pairs = 0;
    fmpq_poly_get_numerator(numerator.raw(), raw_field->pol);
    fmpz_poly_signature(&real_embeddings, &complex_pairs, numerator.raw());

    out.set(real_embeddings, complex_pairs);
    return true;
}

bool is_totally_real(const NumberField& field) noexcept {
    Signature sig;
    return signature(sig, field) && sig.r2() == 0;
}

bool is_totally_complex(const NumberField& field) noexcept {
    Signature sig;
    return signature(sig, field) && sig.r1() == 0;
}

}  // namespace silex
