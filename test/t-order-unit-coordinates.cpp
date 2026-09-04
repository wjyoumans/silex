#include <silex/class_group.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>

#include "order_unit/order_unit_internal.hpp"
#include "test_support.hpp"

#include <cassert>

namespace {
namespace sflint = silex::flint;

struct UnitFixture {
    silex::NumberField field;
    silex::Order order;
    silex::OrderUnitGroup units;
    silex::EmbeddingContext embeddings;
};

UnitFixture quadratic_fixture(slong radicand) noexcept {
    UnitFixture out;
    out.field = silex::test::quadratic_field(radicand);
    out.order = silex::test::equation_order(out.field);
    out.embeddings = silex::EmbeddingContext(out.field);
    assert(out.order.is_maximal());
    assert(out.units.compute(out.order));
    assert(out.units.certification_status() ==
           silex::CertificationMode::proven);
    return out;
}

UnitFixture cubic_disc81_fixture() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);

    UnitFixture out;
    out.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation = silex::test::equation_order(out.field);
    out.order = silex::Order(out.field);
    assert(out.order.maximal_order(equation));
    assert(out.order.is_maximal());
    out.embeddings = silex::EmbeddingContext(out.field);

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                out.order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;
    silex::ClassGroupContext class_group;
    assert(out.units.compute_with_class_group(
            class_group, out.order, sflint::FmpzConstRef(bound), options,
            128));
    assert(out.units.free_rank() == 2);
    assert(out.units.certification_status() ==
           silex::CertificationMode::proven);
    return out;
}

silex::FactoredElement unit_from_coordinates(
        const UnitFixture& fixture,
        slong torsion_exponent,
        const slong* free_exponents) noexcept {
    silex::FactoredElement out(fixture.field);
    assert(out.one());

    silex::OrderElement torsion(fixture.order);
    silex::Element torsion_value(fixture.field);
    assert(fixture.units.torsion_generator(torsion));
    assert(torsion.get_element(torsion_value));
    if (torsion_exponent != 0) {
        assert(out.push(torsion_value, torsion_exponent));
    }

    for (slong i = 0; i < fixture.units.free_rank(); ++i) {
        silex::FactoredElement generator(fixture.field);
        silex::Element generator_value(fixture.field);
        assert(fixture.units.free_generator(generator, i));
        assert(generator.evaluate(generator_value));
        if (free_exponents[i] != 0) {
            assert(out.push(generator_value, free_exponents[i]));
        }
    }
    out.normalize();
    return out;
}

void assert_coordinates(
        const silex::detail::OrdinaryUnitCoordinateResult& result,
        const silex::detail::OrdinaryUnitCoordinates& coordinates,
        slong torsion_exponent,
        const slong* free_exponents,
        slong rank) noexcept {
    assert(result.success);
    assert(result.outcome ==
           silex::detail::OrdinaryUnitCoordinateOutcome::verified);
    assert(result.stage ==
           silex::detail::OrdinaryUnitCoordinateStage::none);
    assert(coordinates.defined);
    assert(coordinates.work_precision == result.work_precision);
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(coordinates.torsion_exponent),
            torsion_exponent));
    assert(sflint::fmpz_mat_nrows(coordinates.free_exponents) == 1);
    assert(sflint::fmpz_mat_ncols(coordinates.free_exponents) == rank);
    for (slong i = 0; i < rank; ++i) {
        assert(sflint::fmpz_equal_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatConstRef(
                                coordinates.free_exponents),
                        0, i),
                free_exponents[i]));
    }
}

int test_rank_zero_and_torsion_coordinates() {
    UnitFixture fixture = quadratic_fixture(-1);
    assert(fixture.units.free_rank() == 0);

    silex::OrderElement torsion(fixture.order);
    silex::Element torsion_value(fixture.field);
    assert(fixture.units.torsion_generator(torsion));
    assert(torsion.get_element(torsion_value));
    silex::FactoredElement factored(fixture.field);
    assert(factored.set_element(torsion_value));

    silex::detail::OrdinaryUnitCoordinateResult result;
    silex::detail::OrdinaryUnitCoordinates coordinates;
    assert(silex::detail::ordinary_unit_coordinates(
            result, coordinates, fixture.units, torsion_value,
            fixture.embeddings, 8, 128));
    assert_coordinates(result, coordinates, 1, nullptr, 0);

    silex::detail::OrdinaryUnitCoordinates factored_coordinates;
    assert(silex::detail::ordinary_unit_coordinates(
            result, factored_coordinates, fixture.units, factored,
            fixture.embeddings, 8, 128));
    assert_coordinates(result, factored_coordinates, 1, nullptr, 0);
    assert(fixture.units.certification_status() ==
           silex::CertificationMode::proven);
    return 0;
}

int test_real_quadratic_round_trip_and_nonunit() {
    UnitFixture fixture = quadratic_fixture(2);
    const slong exponents[] = {7};
    silex::FactoredElement target =
            unit_from_coordinates(fixture, 1, exponents);
    silex::Element expanded(fixture.field);
    assert(target.evaluate(expanded));

    silex::detail::OrdinaryUnitCoordinateResult result;
    silex::detail::OrdinaryUnitCoordinates coordinates;
    assert(silex::detail::ordinary_unit_coordinates(
            result, coordinates, fixture.units, target, fixture.embeddings,
            8, 256));
    assert_coordinates(result, coordinates, 1, exponents, 1);

    silex::detail::OrdinaryUnitCoordinates expanded_coordinates;
    assert(silex::detail::ordinary_unit_coordinates(
            result, expanded_coordinates, fixture.units, expanded,
            fixture.embeddings, 8, 256));
    assert_coordinates(result, expanded_coordinates, 1, exponents, 1);

    silex::detail::OrdinaryUnitCoordinates preserved;
    preserved.free_exponents = sflint::FmpzMat(1, 1);
    sflint::fmpz_set_si(sflint::FmpzRef(preserved.torsion_exponent), 99);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatRef(preserved.free_exponents), 0, 0),
            77);
    preserved.work_precision = 123;
    preserved.defined = true;

    silex::Element two(fixture.field);
    assert(two.set_si(2));
    assert(silex::detail::ordinary_unit_coordinates(
            result, preserved, fixture.units, two, fixture.embeddings, 8,
            256));
    assert(result.success);
    assert(result.outcome ==
           silex::detail::OrdinaryUnitCoordinateOutcome::not_unit);
    assert(result.stage ==
           silex::detail::OrdinaryUnitCoordinateStage::none);
    assert(preserved.defined && preserved.work_precision == 123);
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(preserved.torsion_exponent), 99));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(preserved.free_exponents), 0, 0),
            77));

    silex::FactoredElement factored_two(fixture.field);
    assert(factored_two.set_element(two));
    assert(silex::detail::ordinary_unit_coordinates(
            result, preserved, fixture.units, factored_two,
            fixture.embeddings, 8, 256));
    assert(result.outcome ==
           silex::detail::OrdinaryUnitCoordinateOutcome::not_unit);
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(preserved.torsion_exponent), 99));
    return 0;
}

int test_higher_rank_and_precision_exhaustion() {
    UnitFixture fixture = cubic_disc81_fixture();
    const slong expected[] = {2, -3};
    silex::FactoredElement target =
            unit_from_coordinates(fixture, 1, expected);

    silex::detail::OrdinaryUnitCoordinateResult result;
    silex::detail::OrdinaryUnitCoordinates coordinates;
    assert(silex::detail::ordinary_unit_coordinates(
            result, coordinates, fixture.units, target, fixture.embeddings,
            16, 512));
    assert_coordinates(result, coordinates, 1, expected, 2);

    silex::Element theta(fixture.field);
    silex::Element theta_minus_one(fixture.field);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));
    silex::FactoredElement source_fixture(fixture.field);
    assert(source_fixture.one());
    assert(source_fixture.push(theta, 2));
    assert(source_fixture.push(theta_minus_one, -3));
    assert(silex::detail::ordinary_unit_coordinates(
            result, coordinates, fixture.units, source_fixture,
            fixture.embeddings, 16, 512));
    assert(result.outcome ==
           silex::detail::OrdinaryUnitCoordinateOutcome::verified);

    silex::detail::OrdinaryUnitCoordinates preserved;
    preserved.free_exponents = sflint::FmpzMat(1, 1);
    sflint::fmpz_set_si(sflint::FmpzRef(preserved.torsion_exponent), 41);
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatRef(preserved.free_exponents), 0, 0),
            43);
    preserved.defined = true;
    assert(!silex::detail::ordinary_unit_coordinates(
            result, preserved, fixture.units, source_fixture,
            fixture.embeddings, 1, 1));
    assert(!result.success);
    assert(result.outcome ==
           silex::detail::OrdinaryUnitCoordinateOutcome::unknown);
    assert(result.stage ==
           silex::detail::OrdinaryUnitCoordinateStage::precision_exhausted);
    assert(preserved.defined);
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(preserved.torsion_exponent), 41));
    assert(sflint::fmpz_mat_ncols(preserved.free_exponents) == 1);
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(
                    sflint::FmpzMatConstRef(preserved.free_exponents), 0, 0),
            43));
    assert(fixture.units.certification_status() ==
           silex::CertificationMode::proven);
    return 0;
}

}  // namespace

int main() {
    assert(test_rank_zero_and_torsion_coordinates() == 0);
    assert(test_real_quadratic_round_trip_and_nonunit() == 0);
    assert(test_higher_rank_and_precision_exhaustion() == 0);
    return 0;
}
