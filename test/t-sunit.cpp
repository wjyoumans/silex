#include <silex/sunit.hpp>

#include "test_support.hpp"

#include <cassert>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct ProvenQuadraticFixture {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::EmbeddingContext embeddings;
    sflint::Fmpz factor_base_bound;
};

ProvenQuadraticFixture proven_quadratic(slong radicand) noexcept {
    ProvenQuadraticFixture out;
    out.field = silex::test::quadratic_field(radicand);
    const silex::Order equation = silex::test::equation_order(out.field);
    out.order = silex::Order(out.field);
    assert(out.order.maximal_order(equation));
    assert(out.order.is_maximal());
    out.embeddings = silex::EmbeddingContext(out.field);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 4096;
    options.max_relations = 256;
    options.requested_certification = silex::CertificationMode::proven;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(out.factor_base_bound), out.order));
    if (sflint::fmpz_cmp_ui(
                sflint::FmpzConstRef(out.factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(out.factor_base_bound), 2);
    }
    assert(out.units.compute_with_class_group(
            out.class_group, out.order,
            sflint::FmpzConstRef(out.factor_base_bound), options, 192));
    assert(out.class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(out.units.certification_status() ==
           silex::CertificationMode::proven);
    return out;
}

std::vector<silex::PrimeIdeal> first_prime_above(
        const silex::Order& order,
        slong rational_prime) noexcept {
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    silex::PrimeIdealList decomposition;
    assert(silex::decompose_prime(
            decomposition, order, sflint::FmpzConstRef(p)));
    assert(decomposition.size() > 0 && decomposition.at(0) != nullptr);
    std::vector<silex::PrimeIdeal> out;
    out.emplace_back(order);
    assert(out.back().set(*decomposition.at(0)));
    return out;
}

silex::SUnitCoordinates coordinates(const silex::SUnitGroup& group,
                                    slong torsion,
                                    const slong* ordinary,
                                    const slong* nonunit) noexcept {
    silex::SUnitCoordinates out;
    sflint::fmpz_set_si(sflint::FmpzRef(out.torsion_exponent), torsion);
    out.ordinary_free_exponents =
            sflint::FmpzMat(1, group.ordinary_free_rank());
    out.nonunit_exponents = sflint::FmpzMat(1, group.nonunit_rank());
    for (slong i = 0; i < group.ordinary_free_rank(); ++i) {
        sflint::fmpz_set_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(out.ordinary_free_exponents), 0, i),
                ordinary[i]);
    }
    for (slong i = 0; i < group.nonunit_rank(); ++i) {
        sflint::fmpz_set_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(out.nonunit_exponents), 0, i),
                nonunit[i]);
    }
    out.defined = true;
    return out;
}

void assert_coordinates_equal(const silex::SUnitCoordinates& actual,
                              slong torsion,
                              const slong* ordinary,
                              slong ordinary_rank,
                              const slong* nonunit,
                              slong nonunit_rank) noexcept {
    assert(actual.defined);
    assert(sflint::fmpz_equal_si(
            sflint::FmpzConstRef(actual.torsion_exponent), torsion));
    assert(sflint::fmpz_mat_nrows(actual.ordinary_free_exponents) == 1);
    assert(sflint::fmpz_mat_ncols(actual.ordinary_free_exponents) ==
           ordinary_rank);
    assert(sflint::fmpz_mat_nrows(actual.nonunit_exponents) == 1);
    assert(sflint::fmpz_mat_ncols(actual.nonunit_exponents) == nonunit_rank);
    for (slong i = 0; i < ordinary_rank; ++i) {
        assert(sflint::fmpz_equal_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatConstRef(
                                actual.ordinary_free_exponents),
                        0, i),
                ordinary[i]));
    }
    for (slong i = 0; i < nonunit_rank; ++i) {
        assert(sflint::fmpz_equal_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatConstRef(actual.nonunit_exponents), 0,
                        i),
                nonunit[i]));
    }
}

void assert_regulator_formula(const ProvenQuadraticFixture& fixture,
                              const silex::SClassGroup& s_class_group,
                              const silex::SUnitGroup& s_unit_group,
                              const std::vector<silex::PrimeIdeal>& selected,
                              slong precision) noexcept {
    sflint::Arb expected;
    sflint::Arb actual;
    assert(fixture.units.regulator(sflint::ArbRef(expected)));
    if (!selected.empty()) {
        sflint::Fmpz class_order;
        assert(s_class_group.order(sflint::FmpzRef(class_order)));
        sflint::arb_mul_fmpz(expected, expected,
                            sflint::FmpzConstRef(class_order), precision);
        for (const silex::PrimeIdeal& prime : selected) {
            sflint::Fmpz norm;
            sflint::Arb log_norm;
            assert(prime.norm(sflint::FmpzRef(norm)));
            sflint::arb_log_fmpz(log_norm, sflint::FmpzConstRef(norm),
                                 precision);
            sflint::arb_mul(expected, expected, log_norm, precision);
        }
    }
    assert(s_unit_group.regulator(sflint::ArbRef(actual)));
    assert(sflint::arb_overlaps(actual, expected));
}

void assert_public_s_class_witnesses(
        const silex::SClassGroup& group) noexcept {
    const silex::Order* order = group.parent();
    assert(order != nullptr && order->parent() != nullptr);
    sflint::FmpzMat selected_row(1, group.selected_prime_count());
    for (slong i = 0; i < group.invariant_count(); ++i) {
        sflint::Fmpz invariant;
        silex::FractionalIdeal generator(*order);
        silex::FactoredElement witness(*order->parent());
        assert(group.invariant(sflint::FmpzRef(invariant), i));
        assert(group.invariant_generator(generator, i));
        assert(group.invariant_generator_power_witness(witness, i));
        assert(group.invariant_generator_power_selected_exponents(
                sflint::FmpzMatRef(selected_row), i));

        silex::FractionalIdeal expected(*order);
        silex::FractionalIdeal generator_power(*order);
        assert(generator_power.pow_fmpz(generator,
                                        sflint::FmpzConstRef(invariant)));
        assert(expected.set(generator_power));
        for (slong j = 0; j < group.selected_prime_count(); ++j) {
            silex::PrimeIdeal prime(*order);
            silex::Ideal prime_ideal(*order);
            silex::FractionalIdeal prime_fractional(*order);
            silex::FractionalIdeal prime_power(*order);
            silex::FractionalIdeal product(*order);
            assert(group.selected_prime(prime, j));
            assert(prime.get_ideal(prime_ideal));
            assert(prime_fractional.set_integral(prime_ideal));
            assert(prime_power.pow_fmpz(
                    prime_fractional,
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(selected_row), 0, j)));
            assert(product.multiply(expected, prime_power));
            expected.swap(product);
        }

        silex::Element witness_value(*order->parent());
        silex::FractionalIdeal principal(*order);
        assert(witness.evaluate(witness_value));
        assert(principal.set_principal(witness_value));
        assert(principal.equal(expected));
    }
}

int test_empty_s_publication_and_regulator() {
    ProvenQuadraticFixture fixture = proven_quadratic(-5);
    std::vector<silex::PrimeIdeal> selected;
    silex::SUnitComputeOptions options;
    options.regulator_precision = 192;
    silex::SUnitComputeResult result;
    silex::SClassGroup s_class_group;
    silex::SUnitGroup s_unit_group;
    assert(silex::compute_sunit_groups(
            result, s_class_group, s_unit_group, fixture.class_group,
            fixture.units, silex::PrimeIdealSpan(), options));
    assert(result.success && result.stage == silex::SUnitComputeStage::none);
    assert(s_class_group.is_defined() && s_unit_group.is_defined());
    assert(s_class_group.selected_prime_count() == 0);
    assert(s_unit_group.selected_prime_count() == 0);
    assert(s_class_group.invariant_count() == 1);
    sflint::Fmpz class_order;
    assert(s_class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));
    assert(s_unit_group.ordinary_free_rank() == 0);
    assert(s_unit_group.nonunit_rank() == 0);
    assert(s_unit_group.free_rank() == 0);
    assert(s_unit_group.generator_count() == 1);
    assert(s_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(s_unit_group.certification_status() ==
           silex::CertificationMode::proven);
    assert_regulator_formula(fixture, s_class_group, s_unit_group, selected,
                             192);
    assert_public_s_class_witnesses(s_class_group);

    const slong no_exponents[] = {0};
    silex::SUnitCoordinates input =
            coordinates(s_unit_group, 1, no_exponents, no_exponents);
    silex::FactoredElement value(fixture.field);
    assert(s_unit_group.image(value, input));
    silex::SUnitCoordinates recovered;
    silex::SUnitMembershipResult membership;
    assert(s_unit_group.preimage(membership, recovered, value,
                                 fixture.embeddings, 16, 256));
    assert(membership.outcome == silex::SUnitMembershipOutcome::verified);
    assert_coordinates_equal(recovered, 1, no_exponents, 0, no_exponents, 0);
    return 0;
}

int test_default_pair_sunit_round_trips_and_outside_support() {
    ProvenQuadraticFixture fixture = proven_quadratic(2);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    silex::SUnitComputeOptions options;
    options.regulator_precision = 192;
    silex::SUnitComputeResult result;
    silex::SClassGroup s_class_group;
    silex::SUnitGroup s_unit_group;
    assert(silex::compute_sunit_groups(
            result, s_class_group, s_unit_group, fixture.class_group,
            fixture.units,
            silex::PrimeIdealSpan(selected.data(), selected.size()), options));
    assert(s_unit_group.ordinary_free_rank() == 1);
    assert(s_unit_group.nonunit_rank() == 1);
    assert(s_unit_group.free_rank() == 2);
    assert(s_unit_group.generator_count() == 3);
    sflint::Fmpz s_class_order;
    assert(s_class_group.order(sflint::FmpzRef(s_class_order)));
    assert(sflint::fmpz_is_one(s_class_order));
    sflint::FmpzMat valuations(1, 1);
    assert(s_unit_group.nonunit_valuation_matrix(
            sflint::FmpzMatRef(valuations)));
    assert(::fmpz_mat_rank(valuations.raw()) == 1);
    assert_regulator_formula(fixture, s_class_group, s_unit_group, selected,
                             192);

    ulong state = UWORD(0x5eed1234);
    for (slong sample = 0; sample < 12; ++sample) {
        state = state * UWORD(6364136223846793005) + UWORD(1);
        const slong torsion = static_cast<slong>(state & UWORD(1));
        const slong ordinary[] = {
                static_cast<slong>((state >> 8) % UWORD(7)) - 3};
        const slong nonunit[] = {
                static_cast<slong>((state >> 16) % UWORD(7)) - 3};
        silex::SUnitCoordinates input =
                coordinates(s_unit_group, torsion, ordinary, nonunit);
        silex::FactoredElement value(fixture.field);
        assert(s_unit_group.image(value, input));

        silex::SUnitCoordinates recovered;
        silex::SUnitMembershipResult membership;
        assert(s_unit_group.preimage(membership, recovered, value,
                                     fixture.embeddings, 16, 512));
        assert(membership.success);
        assert(membership.outcome ==
               silex::SUnitMembershipOutcome::verified);
        assert_coordinates_equal(recovered, torsion, ordinary, 1, nonunit, 1);
    }

    const slong ordinary_only[] = {2};
    const slong no_nonunit[] = {0};
    silex::SUnitCoordinates ordinary_coordinates =
            coordinates(s_unit_group, 1, ordinary_only, no_nonunit);
    silex::Element ordinary_value(fixture.field);
    assert(s_unit_group.image(ordinary_value, ordinary_coordinates));
    silex::SUnitCoordinates recovered;
    silex::SUnitMembershipResult membership;
    assert(s_unit_group.preimage(membership, recovered, ordinary_value,
                                 fixture.embeddings, 16, 512));
    assert_coordinates_equal(recovered, 1, ordinary_only, 1, no_nonunit, 1);

    const slong sentinel_ordinary[] = {37};
    const slong sentinel_nonunit[] = {41};
    silex::SUnitCoordinates preserved = coordinates(
            s_unit_group, 1, sentinel_ordinary, sentinel_nonunit);
    silex::Element three(fixture.field);
    assert(three.set_si(3));
    assert(s_unit_group.preimage(membership, preserved, three,
                                 fixture.embeddings, 16, 512));
    assert(membership.success);
    assert(membership.outcome ==
           silex::SUnitMembershipOutcome::not_sunit);
    assert(membership.stage == silex::SUnitMembershipStage::residual_unit);
    assert_coordinates_equal(preserved, 1, sentinel_ordinary, 1,
                             sentinel_nonunit, 1);
    return 0;
}

int test_rank_zero_nonunit_and_nontrivial_s_class() {
    ProvenQuadraticFixture killed = proven_quadratic(-5);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(killed.order, 2);
    silex::SUnitComputeOptions options;
    options.regulator_precision = 192;
    silex::SUnitComputeResult result;
    silex::SClassGroup killed_s_class;
    silex::SUnitGroup killed_s_units;
    assert(silex::compute_sunit_groups(
            result, killed_s_class, killed_s_units, killed.class_group,
            killed.units,
            silex::PrimeIdealSpan(selected.data(), selected.size()), options));
    sflint::Fmpz class_order;
    assert(killed_s_class.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_is_one(class_order));
    assert(killed_s_units.ordinary_free_rank() == 0);
    assert(killed_s_units.nonunit_rank() == 1);
    sflint::FmpzMat valuations(1, 1);
    assert(killed_s_units.nonunit_valuation_row(
            sflint::FmpzMatRef(valuations), 0));
    assert(sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(valuations), 0, 0),
            2));
    const slong no_ordinary[] = {0};
    const slong nonunit[] = {3};
    silex::SUnitCoordinates input =
            coordinates(killed_s_units, 1, no_ordinary, nonunit);
    silex::FactoredElement value(killed.field);
    assert(killed_s_units.image(value, input));
    silex::SUnitCoordinates recovered;
    silex::SUnitMembershipResult membership;
    assert(killed_s_units.preimage(membership, recovered, value,
                                   killed.embeddings, 16, 256));
    assert_coordinates_equal(recovered, 1, no_ordinary, 0, nonunit, 1);

    ProvenQuadraticFixture surviving = proven_quadratic(-14);
    selected = first_prime_above(surviving.order, 2);
    silex::SClassGroup surviving_s_class;
    silex::SUnitGroup surviving_s_units;
    assert(silex::compute_sunit_groups(
            result, surviving_s_class, surviving_s_units,
            surviving.class_group, surviving.units,
            silex::PrimeIdealSpan(selected.data(), selected.size()), options));
    assert(surviving_s_class.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));
    assert(surviving_s_class.invariant_count() == 1);
    assert_public_s_class_witnesses(surviving_s_class);
    assert_regulator_formula(surviving, surviving_s_class, surviving_s_units,
                             selected, 192);
    return 0;
}

int test_fail_closed_preserves_publication() {
    ProvenQuadraticFixture fixture = proven_quadratic(2);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    silex::SUnitComputeResult result;
    silex::SClassGroup s_class_group;
    silex::SUnitGroup s_unit_group;
    assert(silex::compute_sunit_groups(
            result, s_class_group, s_unit_group, fixture.class_group,
            fixture.units,
            silex::PrimeIdealSpan(selected.data(), selected.size())));

    silex::ClassGroupCandidateOptions unknown_options;
    unknown_options.max_candidates = 256;
    unknown_options.max_relations = 32;
    sflint::Fmpz unknown_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(unknown_bound), fixture.order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(unknown_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(unknown_bound), 2);
    }
    silex::ClassGroupContext unknown;
    assert(unknown.compute_candidate(fixture.order,
                                     sflint::FmpzConstRef(unknown_bound),
                                     unknown_options));
    assert(unknown.certification_status() ==
           silex::CertificationMode::unknown);
    assert(!silex::compute_sunit_groups(
            result, s_class_group, s_unit_group, unknown, fixture.units,
            silex::PrimeIdealSpan(selected.data(), selected.size())));
    assert(result.stage == silex::SUnitComputeStage::input_validation);
    sflint::Fmpz class_order;
    assert(s_class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_is_one(class_order));
    assert(s_unit_group.free_rank() == 2);

    return 0;
}

}  // namespace

int main() {
    assert(test_empty_s_publication_and_regulator() == 0);
    assert(test_default_pair_sunit_round_trips_and_outside_support() == 0);
    assert(test_rank_zero_nonunit_and_nontrivial_s_class() == 0);
    assert(test_fail_closed_preserves_publication() == 0);
    return 0;
}
