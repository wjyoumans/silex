#include <silex/class_group.hpp>
#include <silex/element.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>

#include "class_group/class_group_storage_internal.hpp"
#include "class_group/cyclotomic_galois_images_internal.hpp"
#include "class_group/relation_admission_cache_internal.hpp"
#include "class_group/relation_candidate_internal.hpp"
#include "class_group/relation_search_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <initializer_list>

namespace {
namespace sflint = silex::flint;

silex::NumberField field_from(std::initializer_list<slong> coefficients) {
    sflint::FmpzPoly polynomial;
    slong i = 0;
    for (slong c : coefficients) {
        sflint::fmpz_poly_set_coeff_si(sflint::FmpzPolyRef(polynomial), i,
                                       c);
        ++i;
    }
    return silex::test::field_by_polynomial(
            sflint::FmpzPolyConstRef(polynomial));
}

silex::Order maximal_order_of(const silex::NumberField& field) {
    return silex::test::verified_maximal_order(
            silex::test::equation_order(field));
}

ulong cyclotomic_index(const silex::NumberField& field) {
    ulong n = 99;
    assert(silex::detail::cyclotomic_field_index(n, field));
    return n;
}

// Exactly Phi_n is detected; a field that is cyclotomic only after a change
// of variable, or not cyclotomic at all, is not.
int test_detection() {
    assert(cyclotomic_index(field_from({1, 0, 0, 0, 1})) == 8);
    assert(cyclotomic_index(field_from({1, 1, 1, 1, 1})) == 5);
    assert(cyclotomic_index(field_from({1, -1, 1, -1, 1, -1, 1})) == 14);
    // x^4 + x + 1 is not cyclotomic.
    assert(cyclotomic_index(field_from({1, 1, 0, 0, 1})) == 0);
    // Phi_3(x + 1) = x^2 + 3x + 3 defines Q(zeta_3), but is not Phi_3.
    assert(cyclotomic_index(field_from({3, 3, 1})) == 0);
    // Phi_5(2x) / 16 is monic only over Q; not integral Phi_n.
    assert(cyclotomic_index(field_from({1, 2, 4, 8, 16})) == 0);
    return 0;
}

// sigma_k(x) = x^k on Phi_8, and every image matrix preserves the norm of an
// order element.
int test_image_matrices() {
    silex::NumberField field = field_from({1, 0, 0, 0, 1});
    silex::Order order = maximal_order_of(field);
    silex::detail::CyclotomicGaloisImages images;
    assert(silex::detail::prepare_cyclotomic_galois_images(images, order));
    assert(images.initialized && images.available);
    assert(images.conductor == 8);
    assert(images.exponents.size() == 3);
    assert(images.exponents[0] == 3 && images.exponents[1] == 5 &&
           images.exponents[2] == 7);

    silex::Element generator(field);
    assert(generator.gen());
    sflint::FmpqMat generator_coordinates(1, order.degree());
    assert(order.coordinates(sflint::FmpqMatRef(generator_coordinates),
                             generator));
    sflint::FmpzMat x_row(1, order.degree());
    for (slong j = 0; j < order.degree(); ++j) {
        const fmpq* entry =
                fmpq_mat_entry(generator_coordinates.raw(), 0, j);
        assert(fmpz_is_one(fmpq_denref(entry)));
        fmpz_set(fmpz_mat_entry(x_row.raw(), 0, j), fmpq_numref(entry));
    }

    // alpha = 3 - 2x + x^2 + 5x^3 in coordinates on the order basis.
    sflint::FmpzMat alpha_row(1, order.degree());
    {
        sflint::FmpqPoly polynomial;
        const slong coefficients[] = {3, -2, 1, 5};
        for (slong i = 0; i < 4; ++i) {
            sflint::fmpq_poly_set_coeff_si(sflint::FmpqPolyRef(polynomial),
                                           i, coefficients[i]);
        }
        silex::Element alpha(field);
        assert(alpha.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
        silex::OrderElement alpha_order;
        assert(alpha_order.define(order));
        assert(alpha_order.set_element(alpha));
        assert(alpha_order.get_coordinates(sflint::FmpzMatRef(alpha_row)));
    }

    for (std::size_t i = 0; i < images.exponents.size(); ++i) {
        const slong index = static_cast<slong>(i);
        // The image of x is x^k.
        sflint::FmpzMat image(0, 0);
        assert(silex::detail::cyclotomic_galois_image_coordinates(
                image, images, index, sflint::FmpzMatConstRef(x_row)));
        silex::OrderElement image_element;
        assert(image_element.define(order));
        assert(image_element.set_coordinates(sflint::FmpzMatConstRef(image)));
        silex::Element image_value(field);
        assert(image_element.get_element(image_value));
        silex::Element expected(field);
        sflint::Fmpz exponent;
        fmpz_set_ui(exponent.raw(), images.exponents[i]);
        assert(expected.pow_fmpz(generator, sflint::FmpzConstRef(exponent)));
        assert(image_value.equal(expected));

        // N(sigma_k(alpha)) = N(alpha).
        assert(silex::detail::cyclotomic_galois_image_coordinates(
                image, images, index, sflint::FmpzMatConstRef(alpha_row)));
        silex::OrderElement alpha_element;
        silex::OrderElement alpha_image;
        assert(alpha_element.define(order) && alpha_image.define(order));
        assert(alpha_element.set_coordinates(
                sflint::FmpzMatConstRef(alpha_row)));
        assert(alpha_image.set_coordinates(sflint::FmpzMatConstRef(image)));
        silex::Element alpha_value(field);
        silex::Element alpha_image_value(field);
        assert(alpha_element.get_element(alpha_value));
        assert(alpha_image.get_element(alpha_image_value));
        assert(!alpha_image_value.equal(alpha_value));
        sflint::Fmpq norm;
        sflint::Fmpq image_norm;
        assert(alpha_value.norm(sflint::FmpqRef(norm)));
        assert(alpha_image_value.norm(sflint::FmpqRef(image_norm)));
        assert(fmpq_equal(norm.raw(), image_norm.raw()));
    }

    // Non-cyclotomic and degree-two cyclotomic fields get no images.
    silex::detail::CyclotomicGaloisImages none;
    assert(silex::detail::prepare_cyclotomic_galois_images(
            none, maximal_order_of(field_from({1, 1, 0, 0, 1}))));
    assert(none.initialized && !none.available);
    assert(silex::detail::prepare_cyclotomic_galois_images(
            none, maximal_order_of(field_from({1, 0, 1}))));
    assert(none.initialized && !none.available);
    return 0;
}

// A prime-ideal-norm-bounded base of a Galois field has complete
// rational-prime blocks, since all primes above p share one residue degree.
// In Q(2^(1/3)) the base of norm at most 10 holds the degree-one prime above
// 5 but not the degree-two one (norm 25), so it is not stable.
bool norm_bounded_base_complete(const silex::Order& order, ulong bound) {
    sflint::Fmpz value;
    fmpz_set_ui(value.raw(), bound);
    silex::FactorBase base(order);
    assert(base.build_prime_ideal_norm_bounded(sflint::FmpzConstRef(value)));
    bool complete = false;
    assert(silex::detail::factor_base_rational_prime_blocks_complete(complete,
                                                                     base));
    return complete;
}

int test_factor_base_stability() {
    assert(norm_bounded_base_complete(
            maximal_order_of(field_from({1, 0, 0, 0, 1})), 30));
    assert(norm_bounded_base_complete(
            maximal_order_of(field_from({1, 1, 1, 1, 1})), 30));
    assert(!norm_bounded_base_complete(
            maximal_order_of(field_from({-2, 0, 0, 1})), 10));
    return 0;
}

// Element sum c_i x^i of `field`.
silex::Element element_from(const silex::NumberField& field,
                            std::initializer_list<slong> coefficients) {
    sflint::FmpqPoly polynomial;
    slong i = 0;
    for (slong c : coefficients) {
        sflint::fmpq_poly_set_coeff_si(sflint::FmpqPolyRef(polynomial), i,
                                       c);
        ++i;
    }
    silex::Element element(field);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    return element;
}

// Q(zeta_5) with the factor base of prime-ideal norm at most 30: P5, the
// inert prime 2 (norm 16) and the four primes above 11.  pi = 2 + zeta has
// norm Phi_5(-2) = 11, and its conjugates 2 + zeta^k (k = 2, 3, 4) generate
// the other three primes above 11.  The Galois-image step runs only from
// degree kLargeFieldRelationSearchMinDegree on; unless `crossover` is set,
// the setup lowers that degree for its context to 4 (an internal test seam)
// so that the step runs on Q(zeta_5).
struct Zeta5RelationSetup {
    silex::NumberField field = field_from({1, 1, 1, 1, 1});
    silex::Order order = maximal_order_of(field);
    silex::ClassGroupContext context{order};

    explicit Zeta5RelationSetup(bool crossover = false) {
        sflint::Fmpz bound;
        fmpz_set_ui(bound.raw(), 30);
        silex::FactorBase base(order);
        assert(base.build_prime_ideal_norm_bounded(
                sflint::FmpzConstRef(bound)));
        assert(context.set_factor_base(base));
        assert(context.generator_count() == 6);
        if (!crossover) {
            assert(silex::detail::ClassGroupRelationSearchAccess::
                           set_galois_image_min_degree(context, 4));
        }
    }

    // Tries alpha = pi * sigma_2(pi) as a search candidate of certification
    // `requested` while generator relations are kept unless they duplicate
    // a stored row (the policy of the LLL route), with kernel-unit target
    // `target`.
    void try_pi_sigma2_pi(slong target, bool& goal_reached,
                          silex::CertificationMode requested =
                                  silex::CertificationMode::grh) {
        silex::Element alpha(field);
        assert(alpha.multiply(element_from(field, {2, 1}),
                              element_from(field, {2, 0, 1})));
        silex::OrderElement alpha_order;
        assert(alpha_order.define(order));
        assert(alpha_order.set_element(alpha));
        sflint::FmpzMat coordinates(1, order.degree());
        assert(alpha_order.get_coordinates(sflint::FmpzMatRef(coordinates)));
        using Access = silex::detail::ClassGroupRelationSearchAccess;
        const Access::GeneratorRelationPolicyScope policy(
                context, Access::DependentRelationPolicy::keep_nonduplicate);
        slong candidates_tried = 0;
        slong accepted_relations = context.relation_count();
        bool partial_throttle_exit = false;
        goal_reached = false;
        assert(silex::detail::relation_search::try_coordinate_candidate(
                context, order, coordinates, nullptr,
                silex::ClassGroupRelationSource::Search, target, 100, 100,
                candidates_tried, accepted_relations, requested,
                goal_reached, partial_throttle_exit));
    }

    const silex::detail::CyclotomicGaloisImages& images() {
        const silex::detail::CyclotomicGaloisImages* out =
                silex::detail::ClassGroupRelationSearchAccess::
                        cyclotomic_galois_images(context);
        assert(out != nullptr);
        return *out;
    }
};

// Below degree kLargeFieldRelationSearchMinDegree (30) a grh candidate gets
// no images: alpha alone is stored and no image is tried.  A proven
// candidate gets none at any degree.
int test_no_images_below_crossover_or_proven() {
    static_assert(silex::detail::kLargeFieldRelationSearchMinDegree == 30);
    Zeta5RelationSetup crossover(true);
    assert(silex::detail::ClassGroupRelationSearchAccess::
                   galois_image_min_degree(crossover.context) == 30);
    bool goal_reached = false;
    crossover.try_pi_sigma2_pi(100, goal_reached);
    assert(crossover.images().image_candidates == 0);
    assert(crossover.context.relation_count() == 1);

    Zeta5RelationSetup proven;
    proven.try_pi_sigma2_pi(100, goal_reached,
                            silex::CertificationMode::proven);
    assert(proven.images().image_candidates == 0);
    assert(proven.context.relation_count() == 1);
    return 0;
}

// Before full modular rank an image is kept only when it raises the rank
// (reference add_rel_i with m = NULL, buch2.c:2350-2351; the generator
// path has no surplus allowance).  With P11 = a, b, c, d for pi,
// sigma_2(pi), sigma_4(pi), sigma_3(pi), alpha = pi * sigma_2(pi) gives
// a + b and its images b + c, a + d and c + d; any two of the images span
// the third with a + b, so the third is not kept.
int test_images_before_full_rank() {
    Zeta5RelationSetup setup;
    bool goal_reached = false;
    setup.try_pi_sigma2_pi(100, goal_reached);
    assert(!goal_reached);
    assert(setup.images().image_candidates == 3);
    assert(setup.images().image_relations == 2);
    assert(setup.context.relation_count() == 3);
    assert(setup.context.relation_rank() == 3);
    using Access = silex::detail::ClassGroupRelationSearchAccess;
    assert(Access::galois_image_relations(setup.context) == 2);
    assert(Access::galois_image_kernel_relations(setup.context) == 0);
    return 0;
}

// At full modular rank every stored relation gets its images, and every
// image that is not already known is kept (buch2.c:2294, 2350).  The
// images b + c, a + d and c + d of alpha are new rows, so all three are
// stored; none raises the rank, so all three are kernel images, which the
// kernel-unit target leaves out.
int test_images_at_full_rank() {
    Zeta5RelationSetup setup;
    for (const silex::Element& generator :
         {element_from(setup.field, {2}), element_from(setup.field, {1, -1}),
          element_from(setup.field, {2, 1}),
          element_from(setup.field, {2, 0, 1}),
          element_from(setup.field, {2, 0, 0, 1}),
          element_from(setup.field, {1, -1, -1, -1})}) {
        assert(setup.context.try_append_generator_relation(generator));
    }
    assert(setup.context.relation_count() == 6);
    assert(setup.context.relation_rank() == 6);

    bool goal_reached = false;
    setup.try_pi_sigma2_pi(2, goal_reached);
    assert(setup.images().image_candidates == 3);
    assert(setup.images().image_relations == 3);
    assert(setup.context.relation_count() == 10);
    assert(setup.context.relation_rank() == 6);
    using Access = silex::detail::ClassGroupRelationSearchAccess;
    assert(Access::galois_image_relations(setup.context) == 3);
    assert(Access::galois_image_kernel_relations(setup.context) == 3);
    assert(Access::capped_relation_count(setup.context) == 7);
    return 0;
}

// The admission-cache path (the port of the reference relation cache)
// admits each image as a non-random relation, as reference add_rel_i with
// m = NULL: one that raises the rank, or a dependent one while the
// surplus allowance relsup is open (buch2.c:2351-2361).  Initialization
// stores the rational-prime relations of 5, 2 and 11 (the last is
// a + b + c + d), so after a + b and its image b + c the images a + d and
// c + d are dependent and stored under the allowance.
int test_admission_cache_images() {
    Zeta5RelationSetup setup;
    silex::detail::RelationAdmissionCache cache;
    assert(silex::detail::initialize_relation_admission_cache(
            setup.context, cache, setup.order, 3));
    assert(setup.context.relation_count() == 3);
    assert(cache.missing == 3);
    assert(cache.relation_surplus == 3);

    silex::Element alpha(setup.field);
    assert(alpha.multiply(element_from(setup.field, {2, 1}),
                          element_from(setup.field, {2, 0, 1})));
    silex::OrderElement alpha_order;
    assert(alpha_order.define(setup.order));
    assert(alpha_order.set_element(alpha));
    sflint::FmpzMat coordinates(1, setup.order.degree());
    assert(alpha_order.get_coordinates(sflint::FmpzMatRef(coordinates)));
    sflint::Fmpz den;
    fmpz_one(den.raw());
    slong candidates_tried = 0;
    slong accepted_relations = 0;
    bool goal_reached = false;
    bool partial_throttle_exit = false;
    assert(silex::detail::relation_search::try_coordinate_candidate_den(
            setup.context, setup.order, coordinates,
            sflint::FmpzConstRef(den), nullptr,
            silex::ClassGroupRelationSource::Search, 100, 100, 100,
            candidates_tried, accepted_relations,
            silex::CertificationMode::grh, goal_reached,
            partial_throttle_exit, &cache));
    assert(setup.images().image_candidates == 3);
    assert(setup.images().image_relations == 3);
    assert(setup.context.relation_count() == 7);
    assert(accepted_relations == 4);
    assert(cache.missing == 1);
    assert(cache.relation_surplus == 1);
    return 0;
}

struct ClassUnitRun {
    bool computed = false;
    silex::CertificationMode certification = silex::CertificationMode::unknown;
    slong class_order = 0;
    slong image_candidates = 0;
    slong image_relations = 0;
};

ClassUnitRun class_unit_run(const silex::Order& order,
                            silex::CertificationMode requested) {
    ClassUnitRun run;
    sflint::Fmpz minkowski;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(minkowski),
                                               order));
    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = requested;
    silex::ClassGroupContext context;
    silex::OrderUnitGroup units;
    run.computed = units.compute_with_class_group(
            context, order, sflint::FmpzConstRef(minkowski), options, 128);
    if (!run.computed) {
        return run;
    }
    run.certification = context.certification_status();
    sflint::Fmpz class_order;
    assert(context.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_fits_si(class_order));
    run.class_order = sflint::fmpz_get_si(class_order);
    const silex::detail::CyclotomicGaloisImages* images =
            silex::detail::ClassGroupRelationSearchAccess::
                    cyclotomic_galois_images(context);
    assert(images != nullptr);
    run.image_candidates = images->image_candidates;
    run.image_relations = images->image_relations;
    return run;
}

// A grh class/unit run adds image relations on Q(zeta_64) (degree 32, at
// least the crossover; class group Z/17, Washington, Introduction to
// Cyclotomic Fields, class number tables) and none on Q(zeta_8) (degree 4,
// below it, h = 1).
int test_class_unit_runs() {
    silex::Order small = maximal_order_of(field_from({1, 0, 0, 0, 1}));
    const ClassUnitRun below =
            class_unit_run(small, silex::CertificationMode::grh);
    assert(below.computed);
    assert(below.certification == silex::CertificationMode::grh);
    assert(below.class_order == 1);
    assert(below.image_candidates == 0);

    silex::Order large = maximal_order_of(field_from(
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}));
    const ClassUnitRun above =
            class_unit_run(large, silex::CertificationMode::grh);
    assert(above.computed);
    assert(above.certification == silex::CertificationMode::grh);
    assert(above.class_order == 17);
    assert(above.image_relations > 0);
    assert(above.image_candidates >= above.image_relations);
    return 0;
}

}  // namespace

int main() {
    if (test_detection() != 0 || test_image_matrices() != 0 ||
        test_factor_base_stability() != 0 ||
        test_no_images_below_crossover_or_proven() != 0 ||
        test_images_before_full_rank() != 0 ||
        test_images_at_full_rank() != 0 ||
        test_admission_cache_images() != 0 || test_class_unit_runs() != 0) {
        return 1;
    }
    return 0;
}
