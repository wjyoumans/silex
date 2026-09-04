#include <benchmark/benchmark.h>

#include <silex/archimedean.hpp>
#include <silex/embedding.hpp>
#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

#include "benchmark_contract.hpp"

namespace {
namespace sflint = silex::flint;

bool set_quadratic_field(silex::NumberField& field) noexcept {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(radicand, 5);
    return field.define_quadratic(radicand);
}

bool set_cubic_field(silex::NumberField& field) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -5);
    return field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial));
}

bool set_quadratic_square(silex::Element& out) noexcept {
    const silex::NumberField* field = out.parent();
    if (field == nullptr) {
        return false;
    }

    silex::Element theta(*field);
    silex::Element unit_like(*field);
    return theta.gen() &&
           unit_like.one() &&
           unit_like.add(unit_like, theta) &&
           out.multiply(unit_like, unit_like);
}

bool set_factored_unit_like(silex::FactoredElement& out,
                            const silex::NumberField& field) noexcept {
    silex::Element theta(field);
    silex::Element unit_like(field);
    silex::Element three(field);
    return theta.gen() &&
           unit_like.one() &&
           unit_like.add(unit_like, theta) &&
           three.set_si(3) &&
           out.one() &&
           out.push(unit_like, 18) &&
           out.push(three, -12);
}

bool set_factored_unit_like_value(silex::Element& out,
                                  const silex::NumberField& field) noexcept {
    silex::Element theta(field);
    silex::Element unit_like(field);
    silex::Element unit_power(field);
    silex::Element three(field);
    silex::Element denominator(field);
    silex::Element inverse_denominator(field);
    sflint::Fmpz exponent;
    if (!theta.gen() || !unit_like.one() ||
        !unit_like.add(unit_like, theta)) {
        return false;
    }
    sflint::fmpz_set_ui(exponent, 18);
    if (!unit_power.pow_fmpz(unit_like, sflint::FmpzConstRef(exponent)) ||
        !three.set_si(3)) {
        return false;
    }
    sflint::fmpz_set_ui(exponent, 12);
    return denominator.pow_fmpz(
                   three, sflint::FmpzConstRef(exponent)) &&
           inverse_denominator.invert(denominator) &&
           out.multiply(unit_power, inverse_denominator);
}

bool square_matches(const silex::Element& root,
                    const silex::Element& expected,
                    const silex::NumberField& field) noexcept {
    silex::Element square(field);
    return square.multiply(root, root) && square.equal(expected);
}

bool factored_square_matches(const silex::FactoredElement& root,
                             const silex::Element& expected,
                             const silex::NumberField& field) noexcept {
    silex::FactoredElement square(field);
    silex::Element evaluated(field);
    return square.pow_si(root, 2) && square.evaluate(evaluated) &&
           evaluated.equal(expected);
}

void BM_nf_elt_is_square_quadratic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_quadratic_field(field)) {
        silex::bench_contract::fail(
                state, "quadratic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::Element value(field);
    silex::Element root(field);
    if (!set_quadratic_square(value)) {
        silex::bench_contract::fail(
                state, "quadratic square setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool is_square = false;
    bool call_ok = true;
    bool all_square = true;
    for (auto _ : state) {
        const bool iteration_ok = value.is_square(is_square, root);
        call_ok = iteration_ok && call_ok;
        all_square = is_square && all_square;
        benchmark::DoNotOptimize(is_square);
    }

    benchmark::ClobberMemory();
    state.counters["base"] = 2.0;
    state.counters["is_power"] = is_square ? 1.0 : 0.0;
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "quadratic square operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!all_square) {
        silex::bench_contract::fail(
                state, "quadratic square fixture was not recognized",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    if (!square_matches(root, value, field)) {
        silex::bench_contract::fail(
                state, "quadratic square root did not re-power to input",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_elt_is_power_quadratic_square(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_quadratic_field(field)) {
        silex::bench_contract::fail(
                state, "quadratic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::Element value(field);
    silex::Element root(field);
    sflint::Fmpz exponent;
    sflint::fmpz_set_ui(exponent, 2);
    if (!set_quadratic_square(value)) {
        silex::bench_contract::fail(
                state, "quadratic square setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool is_power = false;
    bool call_ok = true;
    bool all_power = true;
    for (auto _ : state) {
        const bool iteration_ok = value.is_power(
                is_power, root, sflint::FmpzConstRef(exponent));
        call_ok = iteration_ok && call_ok;
        all_power = is_power && all_power;
        benchmark::DoNotOptimize(is_power);
    }

    benchmark::ClobberMemory();
    state.counters["base"] = 2.0;
    state.counters["is_power"] = is_power ? 1.0 : 0.0;
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "quadratic power operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!all_power) {
        silex::bench_contract::fail(
                state, "quadratic square fixture was not a power",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    if (!square_matches(root, value, field)) {
        silex::bench_contract::fail(
                state, "quadratic power root did not re-power to input",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_fac_elt_evaluate_cubic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_cubic_field(field)) {
        silex::bench_contract::fail(
                state, "cubic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::FactoredElement value(field);
    silex::Element evaluated(field);
    silex::Element expected(field);
    if (!set_factored_unit_like(value, field) ||
        !set_factored_unit_like_value(expected, field)) {
        silex::bench_contract::fail(
                state, "factored element fixture setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool call_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = value.evaluate(evaluated);
        call_ok = iteration_ok && call_ok;
    }

    benchmark::ClobberMemory();
    state.counters["factor_count"] = static_cast<double>(value.length());
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "factored element evaluation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!evaluated.equal(expected)) {
        silex::bench_contract::fail(
                state, "factored evaluation disagreed with exact fixture",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_fac_elt_is_power_reduced_cubic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_cubic_field(field)) {
        silex::bench_contract::fail(
                state, "cubic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::FactoredElement value(field);
    silex::FactoredElement root(field);
    silex::Element expected(field);
    if (!set_factored_unit_like(value, field) ||
        !set_factored_unit_like_value(expected, field)) {
        silex::bench_contract::fail(
                state, "factored power fixture setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool is_power = false;
    bool call_ok = true;
    bool all_power = true;
    for (auto _ : state) {
        const bool iteration_ok = value.is_power_si(
                is_power, root, 2, silex::FactoredRootStrategy::reduced);
        call_ok = iteration_ok && call_ok;
        all_power = is_power && all_power;
        benchmark::DoNotOptimize(is_power);
    }

    benchmark::ClobberMemory();
    state.counters["base"] = 2.0;
    state.counters["is_power"] = is_power ? 1.0 : 0.0;
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "reduced factored power operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!all_power) {
        silex::bench_contract::fail(
                state, "reduced factored fixture was not a power",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    if (!factored_square_matches(root, expected, field)) {
        silex::bench_contract::fail(
                state, "reduced factored root did not square to input",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_fac_elt_is_power_compact_cubic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_cubic_field(field)) {
        silex::bench_contract::fail(
                state, "cubic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::FactoredElement value(field);
    silex::FactoredElement root(field);
    silex::Element expected(field);
    if (!set_factored_unit_like(value, field) ||
        !set_factored_unit_like_value(expected, field)) {
        silex::bench_contract::fail(
                state, "compact power fixture setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool is_power = false;
    bool call_ok = true;
    bool all_power = true;
    for (auto _ : state) {
        const bool iteration_ok = value.is_power_si(
                is_power, root, 2, silex::FactoredRootStrategy::compact);
        call_ok = iteration_ok && call_ok;
        all_power = is_power && all_power;
        benchmark::DoNotOptimize(is_power);
    }

    benchmark::ClobberMemory();
    state.counters["base"] = 2.0;
    state.counters["is_power"] = is_power ? 1.0 : 0.0;
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "compact factored power operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!all_power) {
        silex::bench_contract::fail(
                state, "compact factored fixture was not a power",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    if (!factored_square_matches(root, expected, field)) {
        silex::bench_contract::fail(
                state, "compact factored root did not square to input",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_compact_elt_set_fac_elt_cubic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_cubic_field(field)) {
        silex::bench_contract::fail(
                state, "cubic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::FactoredElement value(field);
    silex::CompactElement compact(field);
    silex::Element expected(field);
    if (!set_factored_unit_like(value, field) ||
        !set_factored_unit_like_value(expected, field)) {
        silex::bench_contract::fail(
                state, "compact conversion fixture setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool call_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = compact.set_factored_element(value, 2);
        call_ok = iteration_ok && call_ok;
        benchmark::DoNotOptimize(compact.length());
    }

    benchmark::ClobberMemory();
    state.counters["base"] = static_cast<double>(compact.base());
    state.counters["compact_length"] =
            static_cast<double>(compact.length());
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "compact conversion operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (compact.base() != 2) {
        silex::bench_contract::fail(
                state, "compact conversion did not preserve base two",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    silex::Element evaluated(field);
    if (!compact.evaluate(evaluated)) {
        silex::bench_contract::fail(
                state, "compact conversion could not be evaluated",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!evaluated.equal(expected)) {
        silex::bench_contract::fail(
                state, "compact conversion disagreed with source value",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_compact_elt_root_base_cubic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_cubic_field(field)) {
        silex::bench_contract::fail(
                state, "cubic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::FactoredElement value(field);
    silex::FactoredElement root(field);
    silex::CompactElement compact(field);
    silex::Element expected(field);
    if (!set_factored_unit_like(value, field) ||
        !set_factored_unit_like_value(expected, field) ||
        !compact.set_factored_element(value, 2)) {
        silex::bench_contract::fail(
                state, "compact root fixture setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool call_is_power = true;
    for (auto _ : state) {
        const bool iteration_is_power = compact.root_base(root);
        call_is_power = iteration_is_power && call_is_power;
        benchmark::DoNotOptimize(root.length());
    }

    benchmark::ClobberMemory();
    state.counters["base"] = static_cast<double>(compact.base());
    state.counters["is_power"] = call_is_power ? 1.0 : 0.0;
    if (state.iterations() == 0 || !call_is_power) {
        silex::bench_contract::fail(
                state, "compact base-two fixture produced no square root",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    if (compact.base() != 2) {
        silex::bench_contract::fail(
                state, "compact root fixture lost base two",
                silex::bench_contract::FailureReason::invariant);
        return;
    }
    if (!factored_square_matches(root, expected, field)) {
        silex::bench_contract::fail(
                state, "compact base root did not square to source",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_nf_fac_elt_log_embed_cubic(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    silex::NumberField field;
    if (!set_cubic_field(field)) {
        silex::bench_contract::fail(
                state, "cubic field setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::EmbeddingContext embeddings(field);
    silex::FactoredElement value(field);
    sflint::ArbVec logs(2);
    sflint::ArbVec direct_logs(2);
    silex::Element expanded(field);
    silex::Element expected(field);
    if (!embeddings.refine(128) || !set_factored_unit_like(value, field) ||
        !value.evaluate(expanded) ||
        !set_factored_unit_like_value(expected, field) ||
        !expanded.equal(expected) ||
        !silex::logarithmic_embedding(
                sflint::ArbVecRef(direct_logs), embeddings, expanded,
                silex::LogEmbeddingMode::product, 128)) {
        silex::bench_contract::fail(
                state, "log embedding reference setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    bool call_ok = true;
    for (auto _ : state) {
        const bool iteration_ok = value.logarithmic_embedding(
                sflint::ArbVecRef(logs), embeddings,
                silex::LogEmbeddingMode::product, 128);
        call_ok = iteration_ok && call_ok;
        benchmark::DoNotOptimize(logs.data());
    }

    benchmark::ClobberMemory();
    if (state.iterations() == 0 || !call_ok) {
        silex::bench_contract::fail(
                state, "factored logarithmic embedding failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    const bool overlaps =
            sflint::arb_overlaps(
                    sflint::ArbConstRef(logs.data() + 0),
                    sflint::ArbConstRef(direct_logs.data() + 0)) &&
            sflint::arb_overlaps(
                    sflint::ArbConstRef(logs.data() + 1),
                    sflint::ArbConstRef(direct_logs.data() + 1));
    state.counters["direct_overlap"] = overlaps ? 1.0 : 0.0;
    if (!overlaps) {
        silex::bench_contract::fail(
                state, "factored logs disagreed with direct element logs",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

}  // namespace

BENCHMARK(BM_nf_elt_is_square_quadratic);
BENCHMARK(BM_nf_elt_is_power_quadratic_square);
BENCHMARK(BM_nf_fac_elt_evaluate_cubic);
BENCHMARK(BM_nf_fac_elt_is_power_reduced_cubic);
BENCHMARK(BM_nf_fac_elt_is_power_compact_cubic);
BENCHMARK(BM_nf_compact_elt_set_fac_elt_cubic);
BENCHMARK(BM_nf_compact_elt_root_base_cubic);
BENCHMARK(BM_nf_fac_elt_log_embed_cubic);
BENCHMARK_MAIN();
