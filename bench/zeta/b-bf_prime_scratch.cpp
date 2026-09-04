#include <benchmark/benchmark.h>

#include <flint/arb.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/prime_ideal.hpp>

#include "benchmark_contract.hpp"

namespace {

namespace sflint = silex::flint;

constexpr ulong kBfCutoff = 20007;
constexpr ulong kBfCutoffDiv9 = kBfCutoff / 9;
constexpr slong kBfPrecision = 192;

static_assert(kBfCutoff % 9 == 0);

enum class FieldId : std::size_t {
    cubic_blocker,
    quartic_blocker,
    cubic_control,
    quartic_control,
    quintic_control,
    sextic_control,
};

enum class ScratchVariant {
    per_prime,
    per_term,
};

struct FieldSpec {
    const char* name;
    std::vector<slong> coefficients;
};

struct PrimeSchedule {
    ulong rational_prime = 0;
    std::vector<slong> degree_coefficients;
};

struct BfSummandScratch {
    sflint::Arb log_norm;
    sflint::Arb sqrt_norm_power;
    sflint::Arb term1;
    sflint::Arb term2;
};

struct BfPrimeTermScratch {
    sflint::Arb term;
    sflint::Arb weighted_term;
    BfSummandScratch summand;
};

const std::array<FieldSpec, 6>& field_specs() {
    static const std::array<FieldSpec, 6> specs{{
            {"cubic_blocker", {-5, -2, 0, 1}},
            {"quartic_blocker", {4, 1, 3, -8, 1}},
            {"cubic_control", {-9, -4, -16, 1}},
            {"quartic_control", {13, 6, 9, 4, 1}},
            {"quintic_control", {-1, 0, -15, 14, 6, 1}},
            {"sextic_control", {1, 4, -3, -4, 4, 1, 1}},
    }};
    return specs;
}

bool power_ui_less_than(ulong& out,
                        ulong base,
                        slong exponent,
                        ulong cutoff) noexcept {
    if (base < 2 || exponent <= 0 || cutoff <= 2) {
        return false;
    }

    const ulong limit = cutoff - 1;
    ulong value = 1;
    for (slong i = 0; i < exponent; ++i) {
        if (value > limit / base) {
            return false;
        }
        value *= base;
    }

    out = value;
    return true;
}

void summand_ui(sflint::Arb& out,
                BfSummandScratch& scratch,
                ulong norm_power,
                slong exponent,
                const sflint::Arb& sqrt_cutoff_log_cutoff) noexcept {
    sflint::arb_sqrt_ui(scratch.sqrt_norm_power, norm_power, kBfPrecision);
    sflint::arb_div_ui(scratch.term1, sqrt_cutoff_log_cutoff, norm_power,
                       kBfPrecision);
    sflint::arb_div_ui(scratch.term1, scratch.term1,
                       static_cast<ulong>(exponent), kBfPrecision);
    sflint::arb_div(scratch.term2, scratch.log_norm,
                    sflint::ArbConstRef(scratch.sqrt_norm_power),
                    kBfPrecision);
    sflint::arb_sub(out, scratch.term1, scratch.term2, kBfPrecision);
}

void accumulate_weighted_prime_powers(
        sflint::Arb& sum,
        BfPrimeTermScratch& scratch,
        ulong rational_prime,
        const std::vector<slong>& degree_coefficients,
        slong coefficient_sign,
        ulong cutoff,
        const sflint::Arb& sqrt_cutoff_log_cutoff) noexcept {
    for (slong degree = 1;
         degree < static_cast<slong>(degree_coefficients.size()); ++degree) {
        const slong coefficient =
                coefficient_sign *
                degree_coefficients[static_cast<std::size_t>(degree)];
        if (coefficient == 0) {
            continue;
        }

        ulong norm = 0;
        if (!power_ui_less_than(norm, rational_prime, degree, cutoff)) {
            continue;
        }

        sflint::arb_log_ui(scratch.summand.log_norm, norm, kBfPrecision);
        const slong abs_coefficient =
                coefficient < 0 ? -coefficient : coefficient;
        ulong norm_power = norm;
        const ulong limit = cutoff - 1;
        for (slong exponent = 1;; ++exponent) {
            summand_ui(scratch.term, scratch.summand, norm_power, exponent,
                       sqrt_cutoff_log_cutoff);

            const sflint::Arb* addend = &scratch.term;
            if (abs_coefficient != 1) {
                sflint::arb_mul_ui(
                        scratch.weighted_term, scratch.term,
                        static_cast<ulong>(abs_coefficient), kBfPrecision);
                addend = &scratch.weighted_term;
            }

            if (coefficient > 0) {
                sflint::arb_add(sum, sum, *addend, kBfPrecision);
            } else {
                sflint::arb_sub(sum, sum, *addend, kBfPrecision);
            }

            if (norm_power > limit / norm) {
                break;
            }
            norm_power *= norm;
        }
    }
}

void add_prime_terms(sflint::Arb& sum,
                     BfPrimeTermScratch& scratch,
                     const PrimeSchedule& schedule,
                     const sflint::Arb& sqrt_cutoff_log_cutoff,
                     const sflint::Arb& sqrt_cutoff9_log_cutoff9) noexcept {
    accumulate_weighted_prime_powers(
            sum, scratch, schedule.rational_prime,
            schedule.degree_coefficients, 1, kBfCutoff,
            sqrt_cutoff_log_cutoff);
    if (schedule.rational_prime < kBfCutoffDiv9) {
        accumulate_weighted_prime_powers(
                sum, scratch, schedule.rational_prime,
                schedule.degree_coefficients, -1, kBfCutoffDiv9,
                sqrt_cutoff9_log_cutoff9);
    }
}

// Independent direct transcription of the BF prime-power sum implemented by
// src/zeta/zeta.cpp:bf_accumulate_weighted_prime_powers.  Keep this oracle
// separate from the benchmark helpers above so a shared scratch-kernel error
// cannot validate itself.
bool accumulate_mathematical_oracle(
        sflint::Arb& sum,
        ulong rational_prime,
        const std::vector<slong>& degree_coefficients,
        slong coefficient_sign,
        ulong cutoff) noexcept {
    if ((coefficient_sign != 1 && coefficient_sign != -1) ||
        rational_prime < 2 || cutoff <= 2) {
        return false;
    }

    sflint::Arb sqrt_cutoff;
    sflint::Arb log_cutoff;
    sflint::Arb scaled_cutoff;
    sflint::arb_set_ui(sqrt_cutoff, cutoff);
    sflint::arb_sqrt(sqrt_cutoff, sqrt_cutoff, kBfPrecision);
    sflint::arb_log_ui(log_cutoff, cutoff, kBfPrecision);
    sflint::arb_mul(scaled_cutoff, sqrt_cutoff, log_cutoff, kBfPrecision);

    const ulong limit = cutoff - 1;
    for (std::size_t degree = 1; degree < degree_coefficients.size(); ++degree) {
        const slong coefficient =
                coefficient_sign * degree_coefficients[degree];
        if (coefficient == 0) {
            continue;
        }

        ulong norm = 1;
        bool norm_is_below_cutoff = true;
        for (std::size_t factor = 0; factor < degree; ++factor) {
            if (norm > limit / rational_prime) {
                norm_is_below_cutoff = false;
                break;
            }
            norm *= rational_prime;
        }
        if (!norm_is_below_cutoff) {
            continue;
        }

        sflint::Arb log_norm;
        sflint::arb_log_ui(log_norm, norm, kBfPrecision);
        const ulong coefficient_magnitude = static_cast<ulong>(
                coefficient < 0 ? -coefficient : coefficient);
        ulong norm_power = norm;
        for (ulong exponent = 1;; ++exponent) {
            sflint::Arb sqrt_norm_power;
            sflint::Arb first_term;
            sflint::Arb second_term;
            sflint::Arb summand;
            sflint::arb_sqrt_ui(sqrt_norm_power, norm_power, kBfPrecision);
            sflint::arb_div_ui(first_term, scaled_cutoff, norm_power,
                               kBfPrecision);
            sflint::arb_div_ui(first_term, first_term, exponent, kBfPrecision);
            sflint::arb_div(second_term, log_norm,
                            sflint::ArbConstRef(sqrt_norm_power),
                            kBfPrecision);
            sflint::arb_sub(summand, first_term, second_term, kBfPrecision);
            if (coefficient_magnitude != 1) {
                sflint::arb_mul_ui(summand, summand, coefficient_magnitude,
                                   kBfPrecision);
            }
            if (coefficient > 0) {
                sflint::arb_add(sum, sum, summand, kBfPrecision);
            } else {
                sflint::arb_sub(sum, sum, summand, kBfPrecision);
            }

            if (norm_power > limit / norm) {
                break;
            }
            norm_power *= norm;
        }
    }
    return true;
}

bool evaluate_mathematical_oracle(
        sflint::Arb& sum,
        const std::vector<PrimeSchedule>& schedules) noexcept {
    sflint::arb_zero(sum);
    for (const PrimeSchedule& schedule : schedules) {
        if (!accumulate_mathematical_oracle(
                    sum, schedule.rational_prime,
                    schedule.degree_coefficients, 1, kBfCutoff)) {
            return false;
        }
        if (schedule.rational_prime < kBfCutoffDiv9 &&
            !accumulate_mathematical_oracle(
                    sum, schedule.rational_prime,
                    schedule.degree_coefficients, -1, kBfCutoffDiv9)) {
            return false;
        }
    }
    return arb_is_finite(sum.raw()) != 0;
}

struct FieldFixture {
    explicit FieldFixture(const FieldSpec& spec) {
        sflint::FmpqPoly polynomial;
        for (std::size_t i = 0; i < spec.coefficients.size(); ++i) {
            sflint::fmpq_poly_set_coeff_si(
                    polynomial, static_cast<slong>(i), spec.coefficients[i]);
        }
        if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial)) ||
            !equation_order.define_equation_order(field) ||
            !maximal_order.define(field) ||
            !maximal_order.maximal_order(equation_order)) {
            throw std::runtime_error(std::string(spec.name) +
                                     ": maximal-order setup failed");
        }

        const slong degree = maximal_order.degree();
        sflint::Fmpz rational_prime;
        for (ulong p = 2; p < kBfCutoff; p = n_nextprime(p, 1)) {
            sflint::fmpz_set_ui(sflint::FmpzRef(rational_prime), p);
            silex::PrimeIdealList primes;
            if (!silex::decompose_prime(
                        primes, maximal_order,
                        sflint::FmpzConstRef(rational_prime))) {
                throw std::runtime_error(
                        std::string(spec.name) +
                        ": prime decomposition failed at p=" +
                        std::to_string(p));
            }

            PrimeSchedule schedule;
            schedule.rational_prime = p;
            schedule.degree_coefficients.assign(
                    static_cast<std::size_t>(degree + 1), 0);
            schedule.degree_coefficients[1] = -1;
            for (slong i = 0; i < primes.size(); ++i) {
                const silex::PrimeIdeal* prime = primes.at(i);
                if (prime == nullptr || prime->residue_degree() <= 0 ||
                    prime->residue_degree() > degree) {
                    throw std::runtime_error(
                            std::string(spec.name) +
                            ": invalid residue degree at p=" +
                            std::to_string(p));
                }
                ++schedule.degree_coefficients[static_cast<std::size_t>(
                        prime->residue_degree())];
            }
            schedules.push_back(std::move(schedule));
        }

        sflint::Arb sqrt_cutoff;
        sflint::Arb log_cutoff;
        sflint::arb_set_ui(sqrt_cutoff, kBfCutoff);
        sflint::arb_sqrt(sqrt_cutoff, sqrt_cutoff, kBfPrecision);
        sflint::arb_log_ui(log_cutoff, kBfCutoff, kBfPrecision);
        sflint::arb_mul(sqrt_cutoff_log_cutoff, sqrt_cutoff, log_cutoff,
                        kBfPrecision);

        sflint::Arb sqrt_cutoff9;
        sflint::Arb log_cutoff9;
        sflint::arb_set_ui(sqrt_cutoff9, kBfCutoffDiv9);
        sflint::arb_sqrt(sqrt_cutoff9, sqrt_cutoff9, kBfPrecision);
        sflint::arb_log_ui(log_cutoff9, kBfCutoffDiv9, kBfPrecision);
        sflint::arb_mul(sqrt_cutoff9_log_cutoff9, sqrt_cutoff9, log_cutoff9,
                        kBfPrecision);
        if (!evaluate_mathematical_oracle(mathematical_oracle, schedules)) {
            throw std::runtime_error(std::string(spec.name) +
                                     ": mathematical oracle failed");
        }
    }

    silex::NumberField field;
    silex::Order equation_order;
    silex::Order maximal_order;
    std::vector<PrimeSchedule> schedules;
    sflint::Arb sqrt_cutoff_log_cutoff;
    sflint::Arb sqrt_cutoff9_log_cutoff9;
    sflint::Arb mathematical_oracle;
};

const FieldFixture& field_fixture(FieldId id) {
    static const std::array<std::unique_ptr<FieldFixture>, 6> fixtures = [] {
        std::array<std::unique_ptr<FieldFixture>, 6> out;
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = std::make_unique<FieldFixture>(field_specs()[i]);
        }
        return out;
    }();
    return *fixtures.at(static_cast<std::size_t>(id));
}

void evaluate_per_prime(sflint::Arb& sum,
                        const FieldFixture& fixture) noexcept {
    sflint::arb_zero(sum);
    for (const PrimeSchedule& schedule : fixture.schedules) {
        BfPrimeTermScratch scratch;
        add_prime_terms(sum, scratch, schedule,
                        fixture.sqrt_cutoff_log_cutoff,
                        fixture.sqrt_cutoff9_log_cutoff9);
    }
}

void evaluate_per_term(sflint::Arb& sum,
                       const FieldFixture& fixture) noexcept {
    sflint::arb_zero(sum);
    BfPrimeTermScratch scratch;
    for (const PrimeSchedule& schedule : fixture.schedules) {
        add_prime_terms(sum, scratch, schedule,
                        fixture.sqrt_cutoff_log_cutoff,
                        fixture.sqrt_cutoff9_log_cutoff9);
    }
}

void verify_equivalence() {
    for (std::size_t i = 0; i < field_specs().size(); ++i) {
        const FieldFixture& fixture =
                field_fixture(static_cast<FieldId>(i));
        sflint::Arb per_prime;
        sflint::Arb per_term;
        evaluate_per_prime(per_prime, fixture);
        evaluate_per_term(per_term, fixture);
        if (arb_equal(per_prime.raw(), per_term.raw()) == 0 ||
            arb_equal(per_prime.raw(), fixture.mathematical_oracle.raw()) == 0 ||
            arb_is_finite(per_term.raw()) == 0) {
            throw std::runtime_error(
                    std::string(field_specs()[i].name) +
                    ": scratch-lifetime or mathematical-oracle mismatch");
        }
    }
}

void BM_bf_prime_scratch(benchmark::State& state,
                         FieldId field,
                         ScratchVariant variant) {
    silex::bench_contract::initialize(state);
    const FieldFixture& fixture = field_fixture(field);
    for (auto _ : state) {
        sflint::Arb sum;
        if (variant == ScratchVariant::per_prime) {
            evaluate_per_prime(sum, fixture);
        } else {
            evaluate_per_term(sum, fixture);
        }
        benchmark::DoNotOptimize(sflint::ArbConstRef(sum).raw());
        benchmark::ClobberMemory();
    }
    state.counters["cutoff"] = static_cast<double>(kBfCutoff);
    state.counters["precision"] = static_cast<double>(kBfPrecision);
    state.counters["primes"] =
            static_cast<double>(fixture.schedules.size());
    state.SetItemsProcessed(
            state.iterations() *
            static_cast<std::int64_t>(fixture.schedules.size()));

    sflint::Arb observed;
    sflint::Arb reference;
    if (variant == ScratchVariant::per_prime) {
        evaluate_per_prime(observed, fixture);
        evaluate_per_term(reference, fixture);
    } else {
        evaluate_per_term(observed, fixture);
        evaluate_per_prime(reference, fixture);
    }
    const bool lifetime_reference_equal =
            arb_equal(observed.raw(), reference.raw()) != 0;
    const bool mathematical_oracle_equal =
            arb_equal(observed.raw(), fixture.mathematical_oracle.raw()) != 0;
    const bool result_finite = arb_is_finite(observed.raw()) != 0;
    state.counters["lifetime_reference_equal"] =
            lifetime_reference_equal ? 1.0 : 0.0;
    state.counters["mathematical_oracle_equal"] =
            mathematical_oracle_equal ? 1.0 : 0.0;
    state.counters["result_finite"] = result_finite ? 1.0 : 0.0;
    if (!lifetime_reference_equal || !result_finite) {
        silex::bench_contract::fail(
                state, "prime-scratch result differs from lifetime reference",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    if (!mathematical_oracle_equal) {
        silex::bench_contract::fail(
                state, "prime-scratch result differs from mathematical oracle",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

#define SILEX_BENCHMARK_BF_SCRATCH(field_name, field_id)                  \
    BENCHMARK_CAPTURE(BM_bf_prime_scratch, field_name##_per_prime,        \
                      field_id, ScratchVariant::per_prime);               \
    BENCHMARK_CAPTURE(BM_bf_prime_scratch, field_name##_per_term,         \
                      field_id, ScratchVariant::per_term)

SILEX_BENCHMARK_BF_SCRATCH(cubic_blocker, FieldId::cubic_blocker);
SILEX_BENCHMARK_BF_SCRATCH(quartic_blocker, FieldId::quartic_blocker);
SILEX_BENCHMARK_BF_SCRATCH(cubic_control, FieldId::cubic_control);
SILEX_BENCHMARK_BF_SCRATCH(quartic_control, FieldId::quartic_control);
SILEX_BENCHMARK_BF_SCRATCH(quintic_control, FieldId::quintic_control);
SILEX_BENCHMARK_BF_SCRATCH(sextic_control, FieldId::sextic_control);

#undef SILEX_BENCHMARK_BF_SCRATCH

}  // namespace

int main(int argc, char** argv) {
    verify_equivalence();
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
