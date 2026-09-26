#include <silex/sunit.hpp>

#include "test_support.hpp"

#include <cassert>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct ProvenFixture {
    silex::NumberField field;
    silex::Order order;
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::EmbeddingContext embeddings;
    sflint::Fmpz factor_base_bound;
};

ProvenFixture proven_field(silex::NumberField field,
                            slong requested_bound = 0,
                            slong max_candidates = 4096,
                            slong max_relations = 256,
                            ulong bf_cutoff = 0) noexcept {
    ProvenFixture out;
    out.field = std::move(field);
    const silex::Order equation = silex::test::equation_order(out.field);
    out.order = silex::Order(out.field);
    assert(out.order.maximal_order(equation));
    assert(out.order.is_maximal());
    out.embeddings = silex::EmbeddingContext(out.field);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = max_candidates;
    options.max_relations = max_relations;
    options.zeta_bf_max_cutoff = bf_cutoff;
    options.requested_certification = silex::CertificationMode::proven;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(out.factor_base_bound), out.order));
    if (sflint::fmpz_cmp_ui(
                sflint::FmpzConstRef(out.factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(out.factor_base_bound), 2);
    }
    if (requested_bound > 0) {
        sflint::fmpz_set_si(sflint::FmpzRef(out.factor_base_bound),
                            requested_bound);
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

ProvenFixture proven_quadratic(slong radicand) noexcept {
    return proven_field(silex::test::quadratic_field(radicand));
}

ProvenFixture proven_polynomial(const slong* coefficients,
                                 slong degree,
                                 slong requested_bound = 0,
                                 slong max_candidates = 5000,
                                 slong max_relations = 500) noexcept {
    sflint::FmpqPoly polynomial;
    for (slong i = 0; i <= degree; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }
    return proven_field(silex::test::field_by_polynomial(
                                sflint::FmpqPolyConstRef(polynomial)),
                         requested_bound, max_candidates, max_relations,
                         20000);
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

void assert_regulator_formula(const ProvenFixture& fixture,
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
    assert(sflint::arb_is_finite(actual));
    assert(sflint::arb_is_positive(actual));
    assert(sflint::arb_overlaps(actual, expected));
    if (selected.empty()) {
        assert(::arb_equal(actual.raw(), expected.raw()));
    }
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
        assert(silex::same_order_parent(generator.parent(), order));
        assert(witness.parent() != nullptr &&
               witness.parent()->has_same_data(*order->parent()));

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

void assert_public_nonunit_witnesses(
        const silex::SClassGroup& s_class_group,
        const silex::SUnitGroup& group,
        const std::vector<silex::PrimeIdeal>& selected,
        slong expected_index) noexcept {
    const silex::Order* order = group.parent();
    assert(order != nullptr && order->parent() != nullptr);
    assert(silex::same_order_parent(s_class_group.parent(), order));
    const slong count = static_cast<slong>(selected.size());
    assert(group.selected_prime_count() == count);
    assert(s_class_group.selected_prime_count() == count);
    assert(group.nonunit_rank() == count);
    sflint::FmpzMat matrix(count, count);
    sflint::Fmpz index;
    assert(group.nonunit_valuation_matrix(sflint::FmpzMatRef(matrix)));
    sflint::fmpz_mat_det(sflint::FmpzRef(index), matrix);
    sflint::fmpz_abs(sflint::FmpzRef(index), sflint::FmpzConstRef(index));
    assert(sflint::fmpz_equal_si(index, expected_index));
    for (slong j = 0; j < count; ++j) {
        silex::PrimeIdeal published(*order);
        assert(group.selected_prime(published, j));
        assert(silex::same_order_parent(published.parent(), order));
        assert(published.equal(selected[static_cast<std::size_t>(j)]));
        assert(s_class_group.selected_prime(published, j));
        assert(published.equal(selected[static_cast<std::size_t>(j)]));
    }
    for (slong i = 0; i < count; ++i) {
        silex::FactoredElement generator(*order->parent());
        silex::Element expanded(*order->parent());
        silex::FractionalIdeal principal(*order);
        silex::FractionalIdeal product(*order);
        sflint::FmpzMat row(1, count);
        assert(group.nonunit_generator(generator, i));
        assert(generator.parent() != nullptr &&
               generator.parent()->has_same_data(*order->parent()));
        assert(generator.evaluate(expanded));
        assert(principal.set_principal(expanded));
        assert(product.one());
        assert(group.nonunit_valuation_row(sflint::FmpzMatRef(row), i));
        for (slong j = 0; j < count; ++j) {
            const auto& prime = selected[static_cast<std::size_t>(j)];
            slong compact_valuation = 0;
            slong expanded_valuation = 0;
            assert(prime.valuation(compact_valuation, generator, nullptr));
            assert(prime.valuation(expanded_valuation, expanded, nullptr));
            assert(compact_valuation == expanded_valuation);
            assert(sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(matrix, i, j),
                    expanded_valuation));
            assert(sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(row, 0, j), expanded_valuation));
            silex::Ideal integral(*order);
            silex::FractionalIdeal fractional(*order);
            silex::FractionalIdeal power(*order);
            silex::FractionalIdeal next(*order);
            assert(prime.get_ideal(integral));
            assert(fractional.set_integral(integral));
            assert(power.pow_fmpz(fractional,
                                   sflint::fmpz_mat_entry(
                                           sflint::FmpzMatConstRef(row), 0, j)));
            assert(next.multiply(product, power));
            product.swap(next);
        }
        // This identity excludes support at every prime outside S.
        assert(principal.equal(product));
    }
}

int test_empty_s_publication_and_regulator() {
    ProvenFixture fixture = proven_quadratic(-5);
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
    ProvenFixture fixture = proven_quadratic(2);
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

struct SUnitSourceProofs {
    silex::ProofState relation_saturation = silex::ProofState::not_checked;
    silex::ProofState units = silex::ProofState::not_checked;
    silex::ProofState regulator = silex::ProofState::not_checked;
};

void assert_split_prime_publication(
        const silex::SClassGroup& s_class_group,
        const silex::SUnitGroup& s_unit_group,
        const SUnitSourceProofs& source_proofs,
        const sflint::FmpzMat* prime_hnfs,
        bool reversed) noexcept {
    const silex::Order* order = s_unit_group.parent();
    assert(order != nullptr && order->parent() != nullptr);
    assert(silex::same_order_parent(s_class_group.parent(), order));
    assert(s_class_group.is_defined() && s_unit_group.is_defined());
    assert(s_class_group.selected_prime_count() == 2);
    assert(s_unit_group.selected_prime_count() == 2);
    assert(s_class_group.invariant_count() == 0);
    sflint::Fmpz class_order;
    assert(s_class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_is_one(class_order));
    assert(s_unit_group.ordinary_free_rank() == 1);
    assert(s_unit_group.nonunit_rank() == 2);
    assert(s_unit_group.free_rank() == 3);
    assert(s_unit_group.generator_count() == 4);

    assert(s_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(s_class_group.source_class_certification() ==
           silex::CertificationMode::proven);
    assert(s_class_group.proof_status() == silex::ProofState::verified);
    assert(s_unit_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(s_unit_group.source_class_certification() ==
           silex::CertificationMode::proven);
    assert(s_unit_group.source_unit_certification() ==
           silex::CertificationMode::proven);
    assert(s_unit_group.source_relation_saturation_status() ==
           source_proofs.relation_saturation);
    assert(s_unit_group.source_unit_proof_status() == source_proofs.units);
    assert(s_unit_group.source_regulator_proof_status() ==
           source_proofs.regulator);
    assert(s_unit_group.proof_status() == silex::ProofState::verified);
    assert(s_unit_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(s_unit_group.regulator_precision() == 192);

    sflint::FmpzMat valuations(2, 2);
    sflint::Fmpz index;
    assert(s_unit_group.nonunit_valuation_matrix(
            sflint::FmpzMatRef(valuations)));
    sflint::fmpz_mat_det(sflint::FmpzRef(index), valuations);
    sflint::fmpz_abs(sflint::FmpzRef(index), sflint::FmpzConstRef(index));
    assert(sflint::fmpz_is_one(index));

    for (slong j = 0; j < 2; ++j) {
        silex::PrimeIdeal class_prime(*order);
        silex::PrimeIdeal unit_prime(*order);
        silex::Ideal ideal(*order);
        sflint::FmpzMat hnf(2, 2);
        assert(s_class_group.selected_prime(class_prime, j));
        assert(s_unit_group.selected_prime(unit_prime, j));
        assert(silex::same_order_parent(unit_prime.parent(), order));
        assert(class_prime.equal(unit_prime));
        assert(unit_prime.get_ideal(ideal));
        assert(ideal.get_hnf(sflint::FmpzMatRef(hnf)));
        assert(sflint::fmpz_mat_equal(hnf, prime_hnfs[reversed ? 1 - j : j]));

        for (slong i = 0; i < 2; ++i) {
            silex::FactoredElement generator(*order->parent());
            silex::Element expanded(*order->parent());
            sflint::FmpzMat row(1, 2);
            assert(s_unit_group.nonunit_generator(generator, i));
            assert(generator.evaluate(expanded));
            assert(s_unit_group.nonunit_valuation_row(
                    sflint::FmpzMatRef(row), i));
            slong compact_valuation = 0;
            slong expanded_valuation = 0;
            assert(unit_prime.valuation(compact_valuation, generator,
                                         nullptr));
            assert(unit_prime.valuation(expanded_valuation, expanded,
                                         nullptr));
            assert(compact_valuation == expanded_valuation);
            assert(sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(valuations, i, j),
                    expanded_valuation));
            assert(sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(row, 0, j), expanded_valuation));
        }
    }
}

void assert_mixed_sunit_round_trips(
        silex::Element& expanded,
        const silex::SUnitGroup& group,
        silex::EmbeddingContext& embeddings) noexcept {
    std::vector<slong> ordinary(
            static_cast<std::size_t>(group.ordinary_free_rank()));
    std::vector<slong> nonunit(
            static_cast<std::size_t>(group.nonunit_rank()));
    for (slong i = 0; i < group.ordinary_free_rank(); ++i) {
        ordinary[static_cast<std::size_t>(i)] = i % 2 == 0 ? i + 2 : -i - 2;
    }
    for (slong i = 0; i < group.nonunit_rank(); ++i) {
        nonunit[static_cast<std::size_t>(i)] = i % 2 == 0 ? -i - 2 : i + 2;
    }
    silex::SUnitCoordinates input =
            coordinates(group, 1, ordinary.data(), nonunit.data());
    silex::FactoredElement compact(*expanded.parent());
    silex::Element evaluated(*expanded.parent());
    assert(group.image(compact, input));
    assert(group.image(expanded, input));
    assert(compact.evaluate(evaluated));
    assert(evaluated.equal(expanded));

    silex::SUnitCoordinates recovered;
    silex::SUnitMembershipResult membership;
    assert(group.preimage(membership, recovered, compact, embeddings, 16,
                           512));
    assert(membership.success);
    assert(membership.outcome == silex::SUnitMembershipOutcome::verified);
    assert(membership.stage == silex::SUnitMembershipStage::none);
    assert_coordinates_equal(recovered, 1, ordinary.data(),
                              group.ordinary_free_rank(), nonunit.data(),
                              group.nonunit_rank());

    silex::SUnitCoordinates expanded_coordinates;
    silex::SUnitMembershipResult expanded_membership;
    assert(group.preimage(expanded_membership, expanded_coordinates,
                           expanded, embeddings, 16, 512));
    assert(expanded_membership.success);
    assert(expanded_membership.outcome ==
           silex::SUnitMembershipOutcome::verified);
    assert(expanded_membership.stage == silex::SUnitMembershipStage::none);
    assert_coordinates_equal(expanded_coordinates, 1, ordinary.data(),
                              group.ordinary_free_rank(), nonunit.data(),
                              group.nonunit_rank());
}

// Match the manifest's (p, beta) ideal, independently of decomposition order.
silex::PrimeIdeal prime_from_witness(const ProvenFixture& fixture,
                                      slong rational_prime,
                                      const slong* beta_numerators,
                                      ulong beta_denominator,
                                      slong expected_e,
                                      slong expected_f) noexcept {
    sflint::FmpqPoly polynomial;
    for (slong i = 0; i < fixture.order.degree(); ++i) {
        sflint::Fmpq coefficient;
        sflint::fmpq_set_si(coefficient, beta_numerators[i], beta_denominator);
        sflint::fmpq_poly_set_coeff_fmpq(polynomial, i,
                                        sflint::FmpqConstRef(coefficient));
    }
    silex::Element beta(fixture.field);
    silex::OrderElement integral_beta(fixture.order);
    silex::OrderElement p_element(fixture.order);
    silex::Ideal beta_ideal(fixture.order);
    silex::Ideal p_ideal(fixture.order);
    silex::Ideal expected(fixture.order);
    assert(beta.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
    assert(integral_beta.set_element(beta));
    assert(p_element.set_si(rational_prime));
    assert(beta_ideal.set_principal(integral_beta));
    assert(p_ideal.set_principal(p_element));
    assert(expected.add(p_ideal, beta_ideal));
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    silex::PrimeIdealList decomposition;
    assert(silex::decompose_prime(decomposition, fixture.order,
                                   sflint::FmpzConstRef(p)));
    silex::PrimeIdeal out(fixture.order);
    slong matches = 0;
    for (slong i = 0; i < decomposition.size(); ++i) {
        const auto* prime = decomposition.at(i);
        assert(prime != nullptr);
        silex::Ideal ideal(fixture.order);
        assert(prime->get_ideal(ideal));
        if (ideal.equal(expected)) {
            assert(out.set(*prime));
            ++matches;
        }
    }
    assert(matches == 1);
    assert(out.ramification_index() == expected_e);
    assert(out.residue_degree() == expected_f);
    return out;
}

void assert_replay_publication(ProvenFixture& fixture,
                                const std::vector<silex::PrimeIdeal>& selected,
                                slong expected_class_order,
                                slong expected_ordinary_rank,
                                slong expected_index,
                                slong expected_torsion = 2) noexcept {
    silex::SUnitComputeOptions options;
    options.regulator_precision = 192;
    silex::SUnitComputeResult result;
    silex::SClassGroup s_class;
    silex::SUnitGroup s_units;
    assert(silex::compute_sunit_groups(
            result, s_class, s_units, fixture.class_group, fixture.units,
            silex::PrimeIdealSpan(selected.data(), selected.size()), options));
    assert(result.success && result.stage == silex::SUnitComputeStage::none);
    assert(result.selected_index == -1);
    assert(silex::same_order_parent(s_class.parent(), &fixture.order));
    assert(silex::same_order_parent(s_units.parent(), &fixture.order));
    assert(s_class.source_class_certification() ==
           silex::CertificationMode::proven);
    assert(s_units.source_class_certification() ==
           silex::CertificationMode::proven);
    assert(s_units.source_unit_certification() ==
           silex::CertificationMode::proven);
    assert(s_units.source_relation_saturation_status() ==
           fixture.class_group.relation_saturation_status());
    assert(s_units.source_unit_proof_status() ==
           fixture.class_group.unit_proof_status());
    assert(s_units.source_regulator_proof_status() ==
           fixture.class_group.regulator_proof_status());
    assert(s_class.proof_status() == silex::ProofState::verified);
    assert(s_units.proof_status() == silex::ProofState::verified);
    assert(s_units.regulator_proof_status() == silex::ProofState::verified);
    assert(s_units.regulator_precision() == 192);
    sflint::Fmpz class_order;
    assert(s_class.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, expected_class_order));
    assert(s_class.invariant_count() == (expected_class_order == 1 ? 0 : 1));
    if (expected_class_order != 1) {
        sflint::Fmpz invariant;
        assert(s_class.invariant(sflint::FmpzRef(invariant), 0));
        assert(sflint::fmpz_equal_si(invariant, expected_class_order));
    }
    assert(s_units.ordinary_free_rank() == expected_ordinary_rank);
    assert(s_units.free_rank() == expected_ordinary_rank +
                                         static_cast<slong>(selected.size()));
    sflint::Fmpz torsion_order;
    assert(s_units.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, expected_torsion));
    assert_public_s_class_witnesses(s_class);
    assert_public_nonunit_witnesses(s_class, s_units, selected, expected_index);
    assert_regulator_formula(fixture, s_class, s_units, selected, 192);
    silex::Element expanded(fixture.field);
    assert_mixed_sunit_round_trips(expanded, s_units, fixture.embeddings);
    if (expected_class_order == 2) {
        // Preserve both a nontrivial quotient witness and the unit publication.
        silex::FractionalIdeal before_generator(fixture.order);
        silex::FactoredElement witness(fixture.field);
        silex::Element before_witness(fixture.field);
        sflint::FmpzMat before_row(1, s_class.selected_prime_count());
        sflint::FmpzMat before_valuations(s_units.nonunit_rank(),
                                          s_units.selected_prime_count());
        sflint::Arb before_regulator;
        assert(s_class.invariant_generator(before_generator, 0));
        assert(s_class.invariant_generator_power_witness(witness, 0));
        assert(witness.evaluate(before_witness));
        assert(s_class.invariant_generator_power_selected_exponents(
                sflint::FmpzMatRef(before_row), 0));
        assert(s_units.nonunit_valuation_matrix(
                sflint::FmpzMatRef(before_valuations)));
        assert(s_units.regulator(sflint::ArbRef(before_regulator)));
        std::vector<silex::PrimeIdeal> duplicate;
        for (slong i = 0; i < 2; ++i) {
            duplicate.emplace_back(fixture.order);
            assert(duplicate.back().set(selected[0]));
        }
        assert(!silex::compute_sunit_groups(
                result, s_class, s_units, fixture.class_group, fixture.units,
                silex::PrimeIdealSpan(duplicate.data(), duplicate.size()),
                options));
        assert(!result.success);
        // S-class context validation rejects duplicates without an indexed
        // diagnostic, before either output publishes.
        assert(result.stage == silex::SUnitComputeStage::s_class_context);
        assert(result.selected_index == -1);
        assert(silex::same_order_parent(s_class.parent(), &fixture.order));
        assert(silex::same_order_parent(s_units.parent(), &fixture.order));
        assert(s_class.invariant_count() == 1);
        assert(s_class.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, expected_class_order));
        silex::FractionalIdeal after_generator(fixture.order);
        silex::Element after_witness(fixture.field);
        sflint::FmpzMat after_row(1, s_class.selected_prime_count());
        sflint::FmpzMat after_valuations(s_units.nonunit_rank(),
                                         s_units.selected_prime_count());
        sflint::Arb after_regulator;
        assert(s_class.invariant_generator(after_generator, 0));
        assert(after_generator.equal(before_generator));
        assert(s_class.invariant_generator_power_witness(witness, 0));
        assert(witness.evaluate(after_witness));
        assert(after_witness.equal(before_witness));
        assert(s_class.invariant_generator_power_selected_exponents(
                sflint::FmpzMatRef(after_row), 0));
        assert(sflint::fmpz_mat_equal(before_row, after_row));
        assert(s_units.nonunit_valuation_matrix(
                sflint::FmpzMatRef(after_valuations)));
        assert(sflint::fmpz_mat_equal(before_valuations, after_valuations));
        assert(s_units.regulator(sflint::ArbRef(after_regulator)));
        assert(::arb_equal(before_regulator.raw(), after_regulator.raw()));
        silex::Element after_image(fixture.field);
        assert_mixed_sunit_round_trips(after_image, s_units, fixture.embeddings);
        assert(after_image.equal(expanded));
        assert_public_s_class_witnesses(s_class);
        assert_public_nonunit_witnesses(s_class, s_units, selected, expected_index);
    }

    std::vector<slong> ordinary(static_cast<std::size_t>(expected_ordinary_rank));
    std::vector<slong> nonunit(selected.size());
    silex::SUnitCoordinates negative_torsion =
            coordinates(s_units, -1, ordinary.data(), nonunit.data());
    assert(s_units.image(expanded, negative_torsion));
    silex::SUnitMembershipResult membership;
    silex::SUnitCoordinates recovered;
    assert(s_units.preimage(membership, recovered, expanded,
                             fixture.embeddings, 16, 512));
    assert(membership.success &&
           membership.outcome == silex::SUnitMembershipOutcome::verified);
    assert_coordinates_equal(recovered, expected_torsion - 1, ordinary.data(),
                              expected_ordinary_rank, nonunit.data(),
                              s_units.nonunit_rank());

    for (auto& exponent : ordinary) {
        exponent = 37;
    }
    for (auto& exponent : nonunit) {
        exponent = 41;
    }
    silex::SUnitCoordinates preserved =
            coordinates(s_units, 1, ordinary.data(), nonunit.data());
    const auto assert_preserved = [&]() {
        assert_coordinates_equal(preserved, 1, ordinary.data(),
                                  expected_ordinary_rank, nonunit.data(),
                                  s_units.nonunit_rank());
    };
    // 13 is outside every selected-prime list in these finite replays.
    assert(expanded.set_si(13));
    assert(s_units.preimage(membership, preserved, expanded,
                             fixture.embeddings, 16, 512));
    assert(membership.success);
    assert(membership.outcome == silex::SUnitMembershipOutcome::not_sunit);
    assert(membership.stage == silex::SUnitMembershipStage::residual_unit);
    assert_preserved();
    assert(expanded.set_si(0));
    // The expanded overload first uses FactoredElement::set_element,
    // which rejects zero before the selected-valuation membership stage.
    assert(!s_units.preimage(membership, preserved, expanded,
                              fixture.embeddings, 16, 512));
    assert(!membership.success);
    assert(membership.outcome == silex::SUnitMembershipOutcome::unknown);
    assert(membership.stage == silex::SUnitMembershipStage::input_validation);
    assert_preserved();
    assert(expanded.set_si(1));
    assert(!s_units.preimage(membership, preserved, expanded,
                              fixture.embeddings, 16, 8));
    assert(!membership.success);
    assert(membership.outcome == silex::SUnitMembershipOutcome::unknown);
    assert(membership.stage == silex::SUnitMembershipStage::input_validation);
    assert_preserved();
    if (expected_index == 3) {
        assert(selected.size() == 1);
        assert(expanded.set_si(2));
        slong valuation = 0;
        assert(selected[0].valuation(valuation, expanded, nullptr));
        assert(valuation == 1);
        assert(s_units.preimage(membership, preserved, expanded,
                                 fixture.embeddings, 16, 512));
        assert(membership.success);
        assert(membership.outcome == silex::SUnitMembershipOutcome::not_sunit);
        assert(membership.stage == silex::SUnitMembershipStage::valuation_solve);
        assert_preserved();
    }
}

int test_manifest_witness_replays() {
    // Coefficients and (p, beta) witnesses follow data/sunit_fields.json.
    const slong cubic_coefficients[] = {-2, 0, 0, 1};
    ProvenFixture cubic = proven_polynomial(cubic_coefficients, 3);
    std::vector<silex::PrimeIdeal> selected;
    assert_replay_publication(cubic, selected, 1, 1, 1);
    const slong cubic_beta[][3] = {
            {0, 1, 0}, {1, 1, 0}, {-1, -2, 1}, {2, 1, 0}, {7, 0, 0}};
    const slong primes[] = {2, 3, 5, 5, 7};
    const slong ramification[] = {3, 3, 1, 1, 1};
    const slong residue_degrees[] = {1, 1, 2, 1, 3};
    for (slong i = 0; i < 5; ++i) {
        selected.push_back(prime_from_witness(cubic, primes[i], cubic_beta[i],
                                              1, ramification[i],
                                              residue_degrees[i]));
    }
    assert_replay_publication(cubic, selected, 1, 1, 1);

    const slong imaginary_coefficients[] = {23, 0, 1};
    ProvenFixture imaginary = proven_polynomial(imaginary_coefficients, 2, 5);
    const slong imaginary_beta[] = {-1, 1};
    selected.clear();
    selected.push_back(prime_from_witness(imaginary, 2, imaginary_beta, 2, 1, 1));
    assert_replay_publication(imaginary, selected, 1, 0, 3);

    const slong real_coefficients[] = {-210, 0, 1};
    ProvenFixture real = proven_polynomial(real_coefficients, 2, 10, 10000, 1000);
    const slong real_beta[] = {0, 1};
    selected.clear();
    selected.push_back(prime_from_witness(real, 2, real_beta, 1, 2, 1));
    assert_replay_publication(real, selected, 2, 1, 2);
    // Fundamental unit 29 + 2*sqrt(210), as in t-order-unit.cpp.
    sflint::Arb expected_regulator;
    sflint::Arb regulator;
    sflint::arb_sqrt_ui(expected_regulator, 210, 192);
    sflint::arb_mul_ui(expected_regulator, expected_regulator, 2, 192);
    sflint::arb_add_ui(expected_regulator, expected_regulator, 29, 192);
    sflint::arb_log(expected_regulator, expected_regulator, 192);
    assert(real.units.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_overlaps(regulator, expected_regulator));

    const slong shifted_coefficients[] = {-1, -1, 1};
    ProvenFixture shifted = proven_polynomial(shifted_coefficients, 2);
    const slong split_beta[][2] = {{-4, 1}, {3, 1}};
    selected.clear();
    for (const auto& beta : split_beta) {
        selected.push_back(prime_from_witness(shifted, 11, beta, 1, 1, 1));
    }
    assert_replay_publication(shifted, selected, 1, 1, 1);
    const slong ramified_beta[] = {2, 1};
    selected.clear();
    selected.push_back(prime_from_witness(shifted, 5, ramified_beta, 1, 2, 1));
    assert_replay_publication(shifted, selected, 1, 1, 1);
    return 0;
}

int test_exact_edge_unit_coordinates() {
    const slong degree_one[] = {0, 1};
    const slong imaginary[][3] = {{3, 0, 1}, {1, -1, 1}, {1, 0, 1}};
    const slong torsion_orders[] = {6, 6, 4};
    const std::vector<silex::PrimeIdeal> empty;
    ProvenFixture rational = proven_polynomial(degree_one, 1);
    assert_replay_publication(rational, empty, 1, 0, 1);
    for (slong i = 0; i < 3; ++i) {
        ProvenFixture fixture = proven_polynomial(imaginary[i], 2);
        assert_replay_publication(fixture, empty, 1, 0, 1, torsion_orders[i]);
    }
    return 0;
}

int test_split_prime_order_lifetime_and_failed_publication() {
    silex::SClassGroup s_class_groups[2];
    silex::SUnitGroup s_unit_groups[2];
    SUnitSourceProofs source_proofs;
    // These snapshots retain no field, order, or selected-prime handles.
    sflint::FmpzMat prime_hnfs[] = {
            sflint::FmpzMat(2, 2), sflint::FmpzMat(2, 2)};
    sflint::FmpzMat valuation_matrices[] = {
            sflint::FmpzMat(2, 2), sflint::FmpzMat(2, 2)};
    sflint::FmpqPoly image_polynomials[2];
    sflint::Arb regulators[2];
    {
        ProvenFixture fixture = proven_quadratic(5);
        source_proofs = {fixture.class_group.relation_saturation_status(),
                         fixture.class_group.unit_proof_status(),
                         fixture.class_group.regulator_proof_status()};
        sflint::Fmpz p;
        sflint::fmpz_set_si(sflint::FmpzRef(p), 11);
        silex::PrimeIdealList decomposition;
        assert(silex::decompose_prime(
                decomposition, fixture.order, sflint::FmpzConstRef(p)));
        assert(decomposition.size() == 2);
        assert(decomposition.at(0) != nullptr &&
               decomposition.at(1) != nullptr);
        assert(!decomposition.at(0)->equal(*decomposition.at(1)));
        for (slong j = 0; j < 2; ++j) {
            silex::Ideal ideal(fixture.order);
            sflint::Fmpz norm;
            assert(decomposition.at(j)->get_ideal(ideal));
            assert(ideal.get_hnf(sflint::FmpzMatRef(prime_hnfs[j])));
            assert(decomposition.at(j)->norm(sflint::FmpzRef(norm)));
            assert(sflint::fmpz_equal_si(norm, 11));
        }

        silex::SUnitComputeOptions options;
        options.regulator_precision = 192;
        options.diagnostics = nullptr;
        assert(fixture.class_group.diagnostics() == nullptr);
        silex::SUnitComputeResult result;
        for (slong order_index = 0; order_index < 2; ++order_index) {
            std::vector<silex::PrimeIdeal> selected;
            for (slong j = 0; j < 2; ++j) {
                selected.emplace_back(fixture.order);
                assert(selected.back().set(*decomposition.at(
                        order_index == 0 ? j : 1 - j)));
            }
            auto& s_class_group = s_class_groups[order_index];
            auto& s_unit_group = s_unit_groups[order_index];
            assert(silex::compute_sunit_groups(
                    result, s_class_group, s_unit_group, fixture.class_group,
                    fixture.units,
                    silex::PrimeIdealSpan(selected.data(), selected.size()),
                    options));
            assert(result.success);
            assert(result.stage == silex::SUnitComputeStage::none);
            assert(result.selected_index == -1);
            for (slong j = 0; j < 2; ++j) {
                silex::PrimeIdeal published(fixture.order);
                assert(s_class_group.selected_prime(published, j));
                assert(published.equal(selected[j]));
                assert(s_unit_group.selected_prime(published, j));
                assert(published.equal(selected[j]));
            }
            assert_split_prime_publication(s_class_group, s_unit_group,
                                            source_proofs, prime_hnfs,
                                            order_index != 0);
            assert_regulator_formula(fixture, s_class_group, s_unit_group,
                                     selected, 192);
            silex::Element value(fixture.field);
            assert_mixed_sunit_round_trips(value, s_unit_group,
                                           fixture.embeddings);
            assert(value.get_fmpq_poly(
                    sflint::FmpqPolyRef(image_polynomials[order_index])));
            assert(s_unit_group.nonunit_valuation_matrix(
                    sflint::FmpzMatRef(valuation_matrices[order_index])));
            assert(s_unit_group.regulator(
                    sflint::ArbRef(regulators[order_index])));
        }

        const silex::NumberField foreign_field =
                silex::test::field_by_polynomial(sflint::FmpqPolyConstRef(
                        fixture.field.raw_flint_field()->pol));
        assert(!foreign_field.has_same_data(fixture.field));
        assert(::fmpq_poly_equal(foreign_field.raw_flint_field()->pol,
                                 fixture.field.raw_flint_field()->pol));
        const silex::Order foreign_equation =
                silex::test::equation_order(foreign_field);
        silex::Order foreign_order(foreign_field);
        assert(foreign_order.maximal_order(foreign_equation));
        assert(!foreign_order.has_same_data(fixture.order));
        std::vector<silex::PrimeIdeal> foreign_selected =
                first_prime_above(foreign_order, 11);
        std::vector<silex::PrimeIdeal> invalid_selected;
        invalid_selected.emplace_back(fixture.order);
        assert(invalid_selected.back().set(*decomposition.at(0)));
        invalid_selected.emplace_back(foreign_order);
        assert(invalid_selected.back().set(foreign_selected[0]));
        for (slong order_index = 0; order_index < 2; ++order_index) {
            assert(!silex::compute_sunit_groups(
                    result, s_class_groups[order_index],
                    s_unit_groups[order_index], fixture.class_group,
                    fixture.units,
                    silex::PrimeIdealSpan(invalid_selected.data(),
                                          invalid_selected.size()),
                    options));
            assert(!result.success);
            assert(result.stage == silex::SUnitComputeStage::input_validation);
            assert(result.selected_index == 1);
        }
    }

    // Only the publications keep the original mathematical parents alive.
    const silex::Order* retained_order = s_unit_groups[0].parent();
    assert(retained_order != nullptr && retained_order->parent() != nullptr);
    const silex::NumberField& retained_field = *retained_order->parent();
    silex::EmbeddingContext embeddings(retained_field);
    for (slong order_index = 0; order_index < 2; ++order_index) {
        auto& s_class_group = s_class_groups[order_index];
        auto& s_unit_group = s_unit_groups[order_index];
        assert(silex::same_order_parent(s_unit_group.parent(), retained_order));
        assert_split_prime_publication(s_class_group, s_unit_group,
                                        source_proofs, prime_hnfs,
                                        order_index != 0);
        sflint::FmpzMat valuations(2, 2);
        sflint::Arb regulator;
        assert(s_unit_group.nonunit_valuation_matrix(
                sflint::FmpzMatRef(valuations)));
        assert(sflint::fmpz_mat_equal(valuations,
                                      valuation_matrices[order_index]));
        assert(s_unit_group.regulator(sflint::ArbRef(regulator)));
        assert(::arb_equal(regulator.raw(), regulators[order_index].raw()));
        silex::Element value(retained_field);
        silex::Element expected(retained_field);
        assert_mixed_sunit_round_trips(value, s_unit_group, embeddings);
        assert(expected.set_fmpq_poly(
                sflint::FmpqPolyConstRef(image_polynomials[order_index])));
        assert(value.equal(expected));

        // Obtain the other basis's coordinates; generators need not permute.
        const auto& other = s_unit_groups[1 - order_index];
        silex::SUnitCoordinates other_coordinates;
        silex::SUnitMembershipResult membership;
        assert(other.preimage(membership, other_coordinates, value,
                               embeddings, 16, 512));
        assert(membership.success);
        assert(membership.outcome == silex::SUnitMembershipOutcome::verified);
        assert(membership.stage == silex::SUnitMembershipStage::none);
        silex::FactoredElement compact(retained_field);
        silex::Element expanded(retained_field);
        silex::Element evaluated(retained_field);
        assert(other.image(compact, other_coordinates));
        assert(other.image(expanded, other_coordinates));
        assert(compact.evaluate(evaluated));
        assert(evaluated.equal(value));
        assert(expanded.equal(value));
    }
    return 0;
}

int test_rank_zero_nonunit_and_nontrivial_s_class() {
    ProvenFixture killed = proven_quadratic(-5);
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

    ProvenFixture surviving = proven_quadratic(-14);
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
    ProvenFixture fixture = proven_quadratic(2);
    std::vector<silex::PrimeIdeal> selected =
            first_prime_above(fixture.order, 2);
    silex::SUnitComputeResult result;
    silex::SClassGroup s_class_group;
    silex::SUnitGroup s_unit_group;
    assert(silex::compute_sunit_groups(
            result, s_class_group, s_unit_group, fixture.class_group,
            fixture.units,
            silex::PrimeIdealSpan(selected.data(), selected.size())));
    sflint::FmpzMat before_valuations(1, 1);
    sflint::Arb before_regulator;
    silex::Element before_image(fixture.field);
    assert(s_unit_group.nonunit_valuation_matrix(
            sflint::FmpzMatRef(before_valuations)));
    assert(s_unit_group.regulator(sflint::ArbRef(before_regulator)));
    assert_mixed_sunit_round_trips(before_image, s_unit_group, fixture.embeddings);

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
    assert(silex::same_order_parent(s_class_group.parent(), &fixture.order));
    assert(silex::same_order_parent(s_unit_group.parent(), &fixture.order));
    assert(s_class_group.proof_status() == silex::ProofState::verified);
    assert(s_unit_group.proof_status() == silex::ProofState::verified);
    assert(s_class_group.invariant_count() == 0);
    assert_public_nonunit_witnesses(s_class_group, s_unit_group, selected, 1);
    sflint::FmpzMat after_valuations(1, 1);
    sflint::Arb after_regulator;
    silex::Element after_image(fixture.field);
    assert(s_unit_group.nonunit_valuation_matrix(
            sflint::FmpzMatRef(after_valuations)));
    assert(sflint::fmpz_mat_equal(before_valuations, after_valuations));
    assert(s_unit_group.regulator(sflint::ArbRef(after_regulator)));
    assert(::arb_equal(before_regulator.raw(), after_regulator.raw()));
    assert_mixed_sunit_round_trips(after_image, s_unit_group, fixture.embeddings);
    assert(after_image.equal(before_image));

    return 0;
}

}  // namespace

int main() {
    assert(test_manifest_witness_replays() == 0);
    assert(test_exact_edge_unit_coordinates() == 0);
    assert(test_empty_s_publication_and_regulator() == 0);
    assert(test_default_pair_sunit_round_trips_and_outside_support() == 0);
    assert(test_split_prime_order_lifetime_and_failed_publication() == 0);
    assert(test_rank_zero_nonunit_and_nontrivial_s_class() == 0);
    assert(test_fail_closed_preserves_publication() == 0);
    return 0;
}
