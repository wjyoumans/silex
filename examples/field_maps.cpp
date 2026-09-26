#include <silex/aut.hpp>
#include <silex/element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/hom.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include <iostream>

namespace {
namespace sflint = silex::flint;

int fail(const char* message) {
    std::cerr << "field maps: " << message << "\n";
    return 1;
}

bool set_linear(silex::Element& out, slong constant, slong coefficient) {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

}  // namespace

int main() {
    // field-maps-parents-begin
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    silex::NumberField K =
            silex::NumberField::by_polynomial(sflint::FmpqPolyConstRef(polynomial));
    silex::NumberField L =
            silex::NumberField::by_polynomial(sflint::FmpqPolyConstRef(polynomial));
    if (!K.is_defined() || !L.is_defined() || K.has_same_data(L)) {
        return fail("expected two distinct Q(sqrt(2)) parents");
    }
    silex::NumberField K_copy = K;
    if (!K_copy.has_same_data(K)) {
        return fail("copying a field handle must preserve parent identity");
    }

    silex::Element theta_K(K);
    silex::Element theta_L(L);
    if (!theta_K.gen() || !theta_L.gen()) {
        return fail("generator construction failed");
    }
    // field-maps-parents-end

    // field-maps-hom-begin
    // Parent-only construction does not yet certify a generator image.
    silex::FieldHom cross(K, L);
    if (!cross.is_defined() || cross.has_generator_image()) {
        return fail("unexpected initial field-map state");
    }
    // This checks theta_L^2 - 2 = 0 exactly in L.
    if (!cross.set_generator_image(theta_L) || !cross.has_generator_image() ||
        !cross.is_isomorphism() || cross.is_identity()) {
        return fail("distinct-parent isomorphism certification failed");
    }

    // A nonroot is rejected without replacing the certified image.
    silex::Element one_L(L);
    silex::Element saved_image(L);
    if (!one_L.one()) {
        return fail("constant construction failed");
    }
    if (cross.set_generator_image(one_L)) {
        return fail("a nonroot was accepted as a generator image");
    }
    if (!cross.generator_image(saved_image) || !saved_image.equal(theta_L)) {
        return fail("rejected certification changed the stored image");
    }

    silex::Element x(K);
    silex::Element cross_image(L);
    silex::Element expected_L(L);
    if (!set_linear(x, 3, 4) || !set_linear(expected_L, 3, 4) ||
        !cross.apply(cross_image, x) || !cross_image.equal(expected_L)) {
        return fail("cross-parent image of 3 + 4*theta_K was incorrect");
    }

    // A copied handle shares K's parent, so this map really is the identity.
    silex::FieldHom identity(K, K_copy);
    silex::Element identity_image(K);
    if (!identity.is_defined() || !identity.set_generator_image(theta_K) ||
        !identity.is_identity() || !identity.is_isomorphism() ||
        !identity.apply(identity_image, x) || !identity_image.equal(x)) {
        return fail("same-parent identity failed");
    }
    // field-maps-hom-end

    // field-maps-aut-begin
    silex::FieldAutomorphism identity_aut(K);
    if (!identity_aut.is_defined() || !identity_aut.set_identity() ||
        !identity_aut.is_identity() || !identity_aut.apply(identity_image, x) ||
        !identity_image.equal(x)) {
        return fail("automorphism identity failed");
    }

    silex::FieldAutomorphism conjugation(K);
    silex::Element minus_theta(K);
    silex::Element generator_image(K);
    silex::Element conjugate_x(K);
    silex::Element expected_K(K);
    silex::Element twice(K);
    if (!conjugation.is_defined() || !conjugation.set_quadratic_conjugation() ||
        conjugation.is_identity() || !minus_theta.negate(theta_K) ||
        !conjugation.apply(generator_image, theta_K) ||
        !generator_image.equal(minus_theta)) {
        return fail("quadratic conjugation must send theta_K to -theta_K");
    }
    if (!set_linear(expected_K, 3, -4) ||
        !conjugation.apply(conjugate_x, x) || !conjugate_x.equal(expected_K) ||
        !conjugation.apply(twice, conjugate_x) || !twice.equal(x)) {
        return fail("quadratic conjugation must be an involution");
    }

    silex::FieldHom conjugation_hom;
    silex::Element extracted_image(K);
    if (!conjugation.homomorphism(conjugation_hom) ||
        !conjugation_hom.is_isomorphism() || conjugation_hom.is_identity() ||
        !conjugation_hom.apply(extracted_image, x) ||
        !extracted_image.equal(conjugate_x)) {
        return fail("extracted homomorphism disagrees with conjugation");
    }
    // field-maps-aut-end

    // field-maps-order-begin
    silex::Order order = silex::Order::equation_order(K);
    if (!order.is_defined()) {
        return fail("equation-order construction failed");
    }
    silex::OrderHom order_conjugation(order, order);
    if (!order_conjugation.is_defined() ||
        !order_conjugation.set_field_homomorphism(conjugation_hom) ||
        !order_conjugation.has_field_homomorphism()) {
        return fail("conjugation did not certify integral source-basis images");
    }

    // Rows give images of the source basis (1, theta_K) in the target basis.
    sflint::FmpzMat image_matrix(0, 0);
    if (!order_conjugation.image_matrix(image_matrix) ||
        sflint::fmpz_mat_nrows(image_matrix) != 2 ||
        sflint::fmpz_mat_ncols(image_matrix) != 2) {
        return fail("order image matrix export failed");
    }
    if (!sflint::fmpz_equal_si(sflint::fmpz_mat_entry(image_matrix, 0, 0), 1) ||
        !sflint::fmpz_equal_si(sflint::fmpz_mat_entry(image_matrix, 0, 1), 0) ||
        !sflint::fmpz_equal_si(sflint::fmpz_mat_entry(image_matrix, 1, 0), 0) ||
        !sflint::fmpz_equal_si(sflint::fmpz_mat_entry(image_matrix, 1, 1), -1)) {
        return fail("order image matrix must be diag(1, -1)");
    }
    silex::Element order_image(K);
    if (!order_conjugation.apply(order_image, x) ||
        !order_image.equal(conjugate_x)) {
        return fail("order-map application disagrees with conjugation");
    }
    // field-maps-order-end

    std::cout << "K and L: separately constructed Q(sqrt(2)) parents\n"
              << "theta_K -> theta_L: isomorphism, identity = false\n"
              << "nonroot generator image rejected; certified image preserved\n"
              << "same-parent identity fixes 3 + 4*theta_K\n"
              << "conjugation: 3 + 4*theta_K -> 3 - 4*theta_K\n"
              << "conjugation applied twice fixes 3 + 4*theta_K\n"
              << "OrderHom image matrix for (1, theta_K): [[1, 0], [0, -1]]\n";
    return 0;
}
