#include "cyclotomic_galois_images_internal.hpp"

#include <silex/element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/hom.hpp>

#include "../factor_base/factor_base_internal.hpp"
#include "../unit/unit_internal.hpp"

#include <flint/fmpz_poly.h>
#include <flint/nf.h>
#include <flint/ulong_extras.h>

#include <cstddef>
#include <utility>

namespace silex::detail {

bool cyclotomic_field_index(ulong& n, const NumberField& field) noexcept {
    n = 0;
    const nf_struct* raw = field.raw_flint_field();
    if (raw == nullptr) {
        return false;
    }
    // The defining polynomial must be exactly Phi_n: integral and monic,
    // with no rescaling or change of variable.
    if (!fmpz_is_one(fmpq_poly_denref(raw->pol))) {
        return true;
    }
    flint::FmpzPoly polynomial;
    flint::fmpq_poly_get_numerator(flint::FmpzPolyRef(polynomial),
                                   flint::FmpqPolyConstRef(raw->pol));
    const ulong index = fmpz_poly_is_cyclotomic(polynomial.raw());
    if (index == 0) {
        return true;
    }
    flint::FmpzPoly phi;
    fmpz_poly_cyclotomic(phi.raw(), index);
    if (fmpz_poly_equal(phi.raw(), polynomial.raw())) {
        n = index;
    }
    return true;
}

bool prepare_cyclotomic_galois_images(CyclotomicGaloisImages& out,
                                      const Order& order) noexcept {
    CyclotomicGaloisImages next;
    const NumberField* field = order.parent();
    if (field == nullptr || !order.has_basis()) {
        return false;
    }
    // Share the order's data so Order::has_same_data identifies it later.
    next.order = order;
    next.initialized = true;

    ulong n = 0;
    if (!cyclotomic_field_index(n, *field)) {
        return false;
    }
    bool real_subfield = false;
    if (n == 0) {
        n = real_cyclotomic_defining_polynomial_index(*field);
        real_subfield = n != 0;
    }
    if (n == 0 || field->degree() < 4) {
        out = std::move(next);
        return true;
    }

    Element generator(*field);
    if (!generator.is_defined() || !generator.gen()) {
        return false;
    }
    // On psi_n, the images are the Dickson polynomials D_k(x), from
    // D_0 = 2, D_1 = x and D_(j+1) = x D_j - D_(j-1).
    Element previous(*field);
    Element current(*field);
    if (real_subfield &&
        (!previous.is_defined() || !current.is_defined() ||
         !previous.set_si(2) || !current.set(generator))) {
        return false;
    }
    for (ulong k = 2; real_subfield ? 2 * k < n : k < n; ++k) {
        Element image(*field);
        if (!image.is_defined()) {
            return false;
        }
        if (real_subfield) {
            // current = D_(k-1), previous = D_(k-2); image = D_k.
            Element product(*field);
            if (!product.is_defined() ||
                !product.multiply(generator, current) ||
                !image.subtract(product, previous) ||
                !previous.set(current) || !current.set(image)) {
                return false;
            }
        }
        if (n_gcd(k, n) != 1) {
            continue;
        }
        // Phi_n: sigma_k(x) = x^k.  psi_n: sigma_k(x) = D_k(x).
        // set_generator_image checks f(image) = 0 in K exactly for the
        // defining polynomial f, and image_matrix checks that sigma_k maps
        // the order into itself with integral coordinates.
        if (!real_subfield) {
            flint::Fmpz exponent;
            fmpz_set_ui(exponent.raw(), k);
            if (!image.pow_fmpz(generator,
                                flint::FmpzConstRef(exponent))) {
                return false;
            }
        }
        FieldHom hom(*field, *field);
        OrderHom order_hom(next.order, next.order);
        flint::FmpzMat matrix(0, 0);
        if (!hom.is_defined() || !hom.set_generator_image(image) ||
            !order_hom.is_defined() ||
            !order_hom.set_field_homomorphism(hom) ||
            !order_hom.image_matrix(matrix)) {
            return false;
        }
        next.exponents.push_back(k);
        next.image_matrices.emplace_back(std::move(matrix));
    }

    next.real_subfield = real_subfield;
    next.conductor = n;
    next.available = !next.image_matrices.empty();
    out = std::move(next);
    return true;
}

bool factor_base_rational_prime_blocks_complete(
        bool& complete,
        const FactorBase& base) noexcept {
    complete = false;
    if (!base.is_defined()) {
        return false;
    }
    const slong blocks = base.rational_prime_block_count();
    for (slong block = 0; block < blocks; ++block) {
        bool block_complete = false;
        if (!FactorBaseBlockAccess::rational_prime_block_is_complete(
                    block_complete, base, block)) {
            return false;
        }
        if (!block_complete) {
            return true;
        }
    }
    complete = blocks > 0;
    return true;
}

bool cyclotomic_galois_image_coordinates(
        flint::FmpzMat& out,
        const CyclotomicGaloisImages& images,
        slong index,
        flint::FmpzMatConstRef coordinates) noexcept {
    if (!images.available || index < 0 ||
        index >= static_cast<slong>(images.image_matrices.size())) {
        return false;
    }
    const flint::FmpzMat& matrix =
            images.image_matrices[static_cast<std::size_t>(index)];
    const slong degree = fmpz_mat_nrows(matrix.raw());
    if (fmpz_mat_ncols(coordinates.raw()) != degree) {
        return false;
    }
    flint::FmpzMat product(fmpz_mat_nrows(coordinates.raw()),
                           fmpz_mat_ncols(matrix.raw()));
    fmpz_mat_mul(product.raw(), coordinates.raw(), matrix.raw());
    out = std::move(product);
    return true;
}

}  // namespace silex::detail
