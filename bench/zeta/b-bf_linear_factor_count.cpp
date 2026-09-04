#include <benchmark/benchmark.h>

#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "silex/flint/nmod_poly.hpp"
#include "silex/flint/nmod_poly_factor.hpp"
#include "benchmark_contract.hpp"

namespace {

using silex::flint::NmodPoly;
using silex::flint::NmodPolyFactor;

constexpr ulong kBfCutoff = 20000;
constexpr std::size_t kMaxDegree = 6;
constexpr std::size_t kMaxLength = kMaxDegree + 1;

enum class FieldId : std::size_t {
    cubic_blocker,
    quartic_blocker,
    cubic_control,
    quartic_control,
    quintic_control,
    sextic_control,
};

enum class Variant {
    current,
    array_preinverse,
    array_binexp,
    roots,
    distinct_degree,
};

struct FieldSpec {
    const char* name;
    std::vector<slong> coefficients;
};

struct ReducedInput {
    ulong modulus;
    NmodPoly polynomial;

    ReducedInput(ulong p, const std::vector<slong>& coefficients)
        : modulus(p), polynomial(p) {
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            const slong coefficient = coefficients[i];
            ulong residue = 0;
            if (coefficient >= 0) {
                residue = static_cast<ulong>(coefficient) % p;
            } else {
                const ulong magnitude =
                        static_cast<ulong>(-(coefficient + 1)) + 1;
                residue = magnitude % p;
                if (residue != 0) {
                    residue = p - residue;
                }
            }
            nmod_poly_set_coeff_ui(polynomial.raw(),
                                   static_cast<slong>(i), residue);
        }
    }
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

const FieldSpec& field_spec(FieldId id) {
    return field_specs().at(static_cast<std::size_t>(id));
}

std::vector<ReducedInput> build_inputs(const FieldSpec& spec) {
    std::vector<ReducedInput> inputs;
    for (ulong p = 2; p < kBfCutoff; p = n_nextprime(p, 1)) {
        if (p <= (kBfCutoff - 1) / p) {
            continue;
        }
        ReducedInput input(p, spec.coefficients);
        if (nmod_poly_degree(input.polynomial.raw()) !=
            static_cast<slong>(spec.coefficients.size() - 1)) {
            throw std::runtime_error(std::string(spec.name) +
                                     ": degree changed modulo p");
        }
        if (nmod_poly_is_squarefree(input.polynomial.raw()) == 0) {
            continue;
        }
        inputs.push_back(std::move(input));
    }
    return inputs;
}

const std::vector<ReducedInput>& field_inputs(FieldId id) {
    static const std::array<std::vector<ReducedInput>, 6> inputs{{
            build_inputs(field_spec(FieldId::cubic_blocker)),
            build_inputs(field_spec(FieldId::quartic_blocker)),
            build_inputs(field_spec(FieldId::cubic_control)),
            build_inputs(field_spec(FieldId::quartic_control)),
            build_inputs(field_spec(FieldId::quintic_control)),
            build_inputs(field_spec(FieldId::sextic_control)),
    }};
    return inputs.at(static_cast<std::size_t>(id));
}

slong normalized_length(const ulong* coefficients, slong length) noexcept {
    while (length > 0 && coefficients[length - 1] == 0) {
        --length;
    }
    return length;
}

slong current_count(const NmodPoly& reduced, ulong modulus) {
    const slong degree = nmod_poly_degree(reduced.raw());
    NmodPoly reverse(modulus);
    NmodPoly inverse(modulus);
    NmodPoly frobenius(modulus);
    NmodPoly x(modulus);
    NmodPoly roots(modulus);

    nmod_poly_reverse(reverse.raw(), reduced.raw(), reduced.raw()->length);
    nmod_poly_inv_series(inverse.raw(), reverse.raw(), reduced.raw()->length);
    nmod_poly_powmod_x_ui_preinv(frobenius.raw(), modulus, reduced.raw(),
                                 inverse.raw());
    nmod_poly_set_coeff_ui(x.raw(), 1, 1);
    nmod_poly_sub(frobenius.raw(), frobenius.raw(), x.raw());
    nmod_poly_gcd(roots.raw(), frobenius.raw(), reduced.raw());
    const slong count = nmod_poly_degree(roots.raw());
    if (count < 0 || count > degree) {
        throw std::runtime_error("invalid current linear-factor count");
    }
    return count;
}

slong array_preinverse_count(const NmodPoly& reduced, ulong modulus) {
    const slong length = reduced.raw()->length;
    const slong degree = length - 1;
    if (degree < 2 || degree > static_cast<slong>(kMaxDegree)) {
        throw std::runtime_error("array preinverse degree out of range");
    }

    std::array<ulong, kMaxLength> reverse{};
    std::array<ulong, kMaxLength> inverse{};
    std::array<ulong, kMaxLength> frobenius{};
    std::array<ulong, kMaxLength> roots{};
    for (slong i = 0; i < length; ++i) {
        reverse[static_cast<std::size_t>(i)] =
                reduced.raw()->coeffs[length - 1 - i];
    }

    _nmod_poly_inv_series(inverse.data(), reverse.data(), length, length,
                          reduced.raw()->mod);
    const slong inverse_length = normalized_length(inverse.data(), length);
    _nmod_poly_powmod_x_ui_preinv(
            frobenius.data(), modulus, reduced.raw()->coeffs, length,
            inverse.data(), inverse_length, reduced.raw()->mod);
    frobenius[1] = nmod_sub(frobenius[1], 1, reduced.raw()->mod);
    const slong frobenius_length =
            normalized_length(frobenius.data(), degree);
    if (frobenius_length == 0) {
        return degree;
    }

    const slong gcd_length = _nmod_poly_gcd(
            roots.data(), reduced.raw()->coeffs, length, frobenius.data(),
            frobenius_length, reduced.raw()->mod);
    return gcd_length - 1;
}

slong array_binexp_count(const NmodPoly& reduced, ulong modulus) {
    const slong length = reduced.raw()->length;
    const slong degree = length - 1;
    if (degree < 2 || degree > static_cast<slong>(kMaxDegree)) {
        throw std::runtime_error("array binexp degree out of range");
    }

    std::array<ulong, kMaxLength> x{};
    std::array<ulong, kMaxLength> frobenius{};
    std::array<ulong, kMaxLength> roots{};
    x[1] = 1;
    _nmod_poly_powmod_ui_binexp(
            frobenius.data(), x.data(), modulus, reduced.raw()->coeffs,
            length, reduced.raw()->mod);
    frobenius[1] = nmod_sub(frobenius[1], 1, reduced.raw()->mod);
    const slong frobenius_length =
            normalized_length(frobenius.data(), degree);
    if (frobenius_length == 0) {
        return degree;
    }

    const slong gcd_length = _nmod_poly_gcd(
            roots.data(), reduced.raw()->coeffs, length, frobenius.data(),
            frobenius_length, reduced.raw()->mod);
    return gcd_length - 1;
}

slong roots_count(const NmodPoly& reduced, ulong) {
    NmodPolyFactor roots;
    nmod_poly_roots(roots.raw(), reduced.raw(), 0);
    return roots.raw()->num;
}

slong distinct_degree_count(const NmodPoly& reduced, ulong) {
    NmodPolyFactor factorization;
    std::array<slong, kMaxDegree> degrees{};
    slong* degree_data = degrees.data();
    nmod_poly_factor_distinct_deg(factorization.raw(), reduced.raw(),
                                  &degree_data);
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        if (degrees[static_cast<std::size_t>(i)] == 1) {
            return nmod_poly_degree(factorization.raw()->p + i);
        }
    }
    return 0;
}

using CountFunction = slong (*)(const NmodPoly&, ulong);

CountFunction count_function(Variant variant) {
    switch (variant) {
    case Variant::current:
        return current_count;
    case Variant::array_preinverse:
        return array_preinverse_count;
    case Variant::array_binexp:
        return array_binexp_count;
    case Variant::roots:
        return roots_count;
    case Variant::distinct_degree:
        return distinct_degree_count;
    }
    throw std::runtime_error("unknown linear-factor benchmark variant");
}

void verify_equivalence() {
    for (std::size_t field_index = 0; field_index < field_specs().size();
         ++field_index) {
        const auto id = static_cast<FieldId>(field_index);
        const FieldSpec& spec = field_spec(id);
        for (const ReducedInput& input : field_inputs(id)) {
            const slong expected =
                    current_count(input.polynomial, input.modulus);
            for (Variant variant : {Variant::array_preinverse,
                                    Variant::array_binexp, Variant::roots,
                                    Variant::distinct_degree}) {
                const slong actual = count_function(variant)(
                        input.polynomial, input.modulus);
                if (actual != expected) {
                    throw std::runtime_error(
                            std::string(spec.name) + ": count mismatch at p=" +
                            std::to_string(input.modulus));
                }
            }
        }
    }
}

void BM_bf_linear_factor_count(benchmark::State& state,
                               FieldId field,
                               Variant variant) {
    silex::bench_contract::initialize(state);
    const std::vector<ReducedInput>& inputs = field_inputs(field);
    const CountFunction count = count_function(variant);
    slong observed_checksum = 0;
    for (auto _ : state) {
        observed_checksum = 0;
        for (const ReducedInput& input : inputs) {
            observed_checksum += count(input.polynomial, input.modulus);
        }
        benchmark::DoNotOptimize(observed_checksum);
    }
    state.SetItemsProcessed(
            state.iterations() * static_cast<std::int64_t>(inputs.size()));

    slong reference_checksum = 0;
    for (const ReducedInput& input : inputs) {
        reference_checksum += roots_count(input.polynomial, input.modulus);
    }
    if (observed_checksum != reference_checksum) {
        silex::bench_contract::fail(
                state, "linear-factor count differs from roots reference",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    state.counters["checksum"] = static_cast<double>(observed_checksum);
    state.counters["inputs"] = static_cast<double>(inputs.size());
    silex::bench_contract::succeed(state);
}

#define SILEX_BENCHMARK_BF_VARIANTS(field_name, field_id)                   \
    BENCHMARK_CAPTURE(BM_bf_linear_factor_count, field_name##_current,      \
                      field_id, Variant::current);                          \
    BENCHMARK_CAPTURE(BM_bf_linear_factor_count,                            \
                      field_name##_array_preinverse, field_id,              \
                      Variant::array_preinverse);                            \
    BENCHMARK_CAPTURE(BM_bf_linear_factor_count, field_name##_array_binexp, \
                      field_id, Variant::array_binexp);                      \
    BENCHMARK_CAPTURE(BM_bf_linear_factor_count, field_name##_roots,        \
                      field_id, Variant::roots);                             \
    BENCHMARK_CAPTURE(BM_bf_linear_factor_count,                            \
                      field_name##_distinct_degree, field_id,               \
                      Variant::distinct_degree)

SILEX_BENCHMARK_BF_VARIANTS(cubic_blocker, FieldId::cubic_blocker);
SILEX_BENCHMARK_BF_VARIANTS(quartic_blocker, FieldId::quartic_blocker);
SILEX_BENCHMARK_BF_VARIANTS(cubic_control, FieldId::cubic_control);
SILEX_BENCHMARK_BF_VARIANTS(quartic_control, FieldId::quartic_control);
SILEX_BENCHMARK_BF_VARIANTS(quintic_control, FieldId::quintic_control);
SILEX_BENCHMARK_BF_VARIANTS(sextic_control, FieldId::sextic_control);

#undef SILEX_BENCHMARK_BF_VARIANTS

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
