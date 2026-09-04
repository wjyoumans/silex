#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/embedding.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>

#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>

namespace {
namespace sflint = silex::flint;

struct DiagnosticCounters {
    struct ProfileFrame {
        silex::DiagnosticsModule module = silex::DiagnosticsModule::core;
        const char* function = nullptr;
        const char* label = nullptr;
        std::chrono::steady_clock::time_point start;
    };

    slong verbose = 0;
    slong logs = 0;
    slong debug_failures = 0;
    slong profile_begin = 0;
    slong profile_end = 0;
    slong profile_events = 0;
    slong profile_stack_depth = 0;
    slong profile_stack_overflow = 0;
    slong profile_unmatched_end = 0;
    std::int64_t profile_elapsed_ns = 0;
    std::array<ProfileFrame, 256> profile_stack;
};

bool is_env_on(const char* value) noexcept {
    return value != nullptr && std::strcmp(value, "0") != 0 &&
           std::strcmp(value, "off") != 0 &&
           std::strcmp(value, "false") != 0;
}

silex::VerboseLevel verbose_level_from_env(const char* value) noexcept {
    if (value == nullptr || std::strcmp(value, "off") == 0) {
        return silex::VerboseLevel::off;
    }
    if (std::strcmp(value, "trace") == 0) {
        return silex::VerboseLevel::trace;
    }
    if (std::strcmp(value, "detail") == 0) {
        return silex::VerboseLevel::detail;
    }
    return silex::VerboseLevel::progress;
}

silex::LogLevel log_level_from_env(const char* value) noexcept {
    if (value == nullptr || std::strcmp(value, "off") == 0) {
        return silex::LogLevel::off;
    }
    if (std::strcmp(value, "trace") == 0) {
        return silex::LogLevel::trace;
    }
    if (std::strcmp(value, "detail") == 0) {
        return silex::LogLevel::detail;
    }
    if (std::strcmp(value, "info") == 0) {
        return silex::LogLevel::info;
    }
    if (std::strcmp(value, "warn") == 0) {
        return silex::LogLevel::warn;
    }
    return silex::LogLevel::error;
}

silex::DebugLevel debug_level_from_env(const char* value) noexcept {
    if (value == nullptr || std::strcmp(value, "off") == 0) {
        return silex::DebugLevel::off;
    }
    if (std::strcmp(value, "exhaustive") == 0) {
        return silex::DebugLevel::exhaustive;
    }
    if (std::strcmp(value, "expensive") == 0) {
        return silex::DebugLevel::expensive;
    }
    if (std::strcmp(value, "normal") == 0) {
        return silex::DebugLevel::normal;
    }
    return silex::DebugLevel::cheap;
}

void verbose_callback(void* user,
                      silex::DiagnosticsModule module,
                      silex::VerboseLevel level,
                      const char* function,
                      const char* message,
                      const char* detail) noexcept {
    auto* counters = static_cast<DiagnosticCounters*>(user);
    ++counters->verbose;
    std::cerr << "[verbose:" << silex::diagnostics_module_name(module)
              << ":" << silex::verbose_level_name(level) << "] "
              << function << ": " << message;
    if (detail != nullptr) {
        std::cerr << " (" << detail << ")";
    }
    std::cerr << "\n";
}

void log_callback(void* user,
                  silex::DiagnosticsModule module,
                  silex::LogLevel level,
                  const char* function,
                  const char* message,
                  const char* detail) noexcept {
    auto* counters = static_cast<DiagnosticCounters*>(user);
    ++counters->logs;
    std::cerr << "[log:" << silex::diagnostics_module_name(module)
              << ":" << silex::log_level_name(level) << "] "
              << function << ": " << message;
    if (detail != nullptr) {
        std::cerr << " (" << detail << ")";
    }
    std::cerr << "\n";
}

void debug_failure_callback(void* user,
                            silex::DiagnosticsModule module,
                            silex::DebugLevel level,
                            const char* function,
                            const char* file,
                            int line,
                            const char* label,
                            const char* expression) noexcept {
    auto* counters = static_cast<DiagnosticCounters*>(user);
    ++counters->debug_failures;
    std::cerr << "[debug-failure:" << silex::diagnostics_module_name(module)
              << ":" << silex::debug_level_name(level) << "] "
              << function << " " << file << ":" << line << ": "
              << label << ": " << expression << "\n";
}

void profile_callback(void* user,
                      silex::DiagnosticsModule module,
                      silex::ProfileEvent event,
                      const char* function,
                      const char* label) noexcept {
    auto* counters = static_cast<DiagnosticCounters*>(user);
    const char* event_name = "event";
    std::int64_t elapsed_ns = -1;
    if (event == silex::ProfileEvent::begin_scope) {
        ++counters->profile_begin;
        event_name = "begin";
        if (counters->profile_stack_depth <
            static_cast<slong>(counters->profile_stack.size())) {
            DiagnosticCounters::ProfileFrame& frame =
                    counters->profile_stack[static_cast<std::size_t>(
                            counters->profile_stack_depth)];
            frame.module = module;
            frame.function = function;
            frame.label = label;
            frame.start = std::chrono::steady_clock::now();
            ++counters->profile_stack_depth;
        } else {
            ++counters->profile_stack_overflow;
        }
    } else if (event == silex::ProfileEvent::end_scope) {
        ++counters->profile_end;
        event_name = "end";
        if (counters->profile_stack_depth > 0) {
            --counters->profile_stack_depth;
            const DiagnosticCounters::ProfileFrame& frame =
                    counters->profile_stack[static_cast<std::size_t>(
                            counters->profile_stack_depth)];
            const auto elapsed = std::chrono::steady_clock::now() - frame.start;
            elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                 elapsed)
                                 .count();
            counters->profile_elapsed_ns += elapsed_ns;
        } else {
            ++counters->profile_unmatched_end;
        }
    } else {
        ++counters->profile_events;
    }
    std::cerr << "[profile:" << silex::diagnostics_module_name(module)
              << ":" << event_name << "] " << function << ": "
              << label;
    if (elapsed_ns >= 0) {
        std::cerr << " elapsed_us=" << (elapsed_ns / 1000.0);
    }
    std::cerr << "\n";
}

silex::DiagnosticsContext diagnostics_from_env(
        DiagnosticCounters& counters) noexcept {
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);

    const auto hot_modules =
            silex::diagnostics_module_bit(silex::DiagnosticsModule::element) |
            silex::diagnostics_module_bit(silex::DiagnosticsModule::prime_ideal) |
            silex::diagnostics_module_bit(silex::DiagnosticsModule::relation) |
            silex::diagnostics_module_bit(silex::DiagnosticsModule::class_group) |
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group);

    const silex::VerboseLevel verbose_level =
            verbose_level_from_env(std::getenv("SILEX_DIAGNOSTICS_VERBOSE"));
    silex::diagnostics_set_verbose(diagnostics, verbose_level, hot_modules,
                                   verbose_callback, &counters);

    const silex::LogLevel log_level =
            log_level_from_env(std::getenv("SILEX_DIAGNOSTICS_LOG"));
    silex::diagnostics_set_logging(diagnostics, log_level, hot_modules,
                                   log_callback, &counters);

    const silex::DebugLevel debug_level =
            debug_level_from_env(std::getenv("SILEX_DIAGNOSTICS_DEBUG"));
    silex::diagnostics_set_debug_checks(diagnostics, debug_level, hot_modules,
                                        debug_failure_callback, &counters);

    silex::diagnostics_set_profiling(
            diagnostics, is_env_on(std::getenv("SILEX_DIAGNOSTICS_PROFILE")),
            hot_modules, profile_callback, &counters);
    return diagnostics;
}

silex::NumberField rational_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    silex::NumberField field =
            silex::NumberField::by_polynomial(sflint::FmpqPolyConstRef(polynomial));
    assert(field.is_defined());
    return field;
}

void run_factored_root_example(const silex::DiagnosticsContext& diagnostics) {
    silex::NumberField field = rational_field();
    silex::Element two(field);
    silex::Element three(field);
    assert(two.set_si(2));
    assert(three.set_si(3));

    silex::FactoredElement product(field);
    assert(product.push(two, 3));
    assert(product.push(three, 6));

    silex::CompactElement compact(field);
    assert(compact.set_factored_element(product, 3, &diagnostics));

    silex::FactoredElement compact_root(field);
    bool is_power = false;
    assert(product.is_power_si(is_power, compact_root, 3,
                               silex::FactoredRootStrategy::compact,
                               &diagnostics));
    assert(is_power);

    silex::FactoredElement exact_root(field);
    assert(exact_root.root_si(product, 3, &diagnostics));
}

void run_class_unit_example(const silex::DiagnosticsContext& diagnostics) {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), 2);
    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());

    sflint::Fmpz factor_base_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(factor_base_bound), 2);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    options.requested_certification = silex::CertificationMode::proven;
    options.diagnostics = &diagnostics;

    silex::ClassGroupContext class_group;
    class_group.set_diagnostics(&diagnostics);
    silex::OrderUnitGroup units;
    units.set_diagnostics(&diagnostics);

    assert(units.compute_with_class_group(
            class_group, order, sflint::FmpzConstRef(factor_base_bound),
            options, 160));
}

void run_unit_saturation_example(const silex::DiagnosticsContext& diagnostics) {
    sflint::Fmpz radicand;
    sflint::fmpz_set_si(sflint::FmpzRef(radicand), 2);
    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Order order = silex::Order::equation_order(field);
    assert(order.is_defined());
    silex::EmbeddingContext embeddings(field);

    silex::Element theta(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(theta.gen());
    assert(epsilon.add_si(theta, 1));
    assert(epsilon2.multiply(epsilon, epsilon));

    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group(order);
    group.set_diagnostics(&diagnostics);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    sflint::fmpz_set_si(sflint::FmpzRef(p), 5);
    sflint::fmpz_set_si(sflint::FmpzRef(ell), 2);

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    silex::OrderUnitGroup saturated(order);
    saturated.set_diagnostics(&diagnostics);
    bool changed = false;
    assert(saturated.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(changed);
}

}  // namespace

int main() {
    DiagnosticCounters counters;
    silex::DiagnosticsContext diagnostics = diagnostics_from_env(counters);

    run_factored_root_example(diagnostics);
    run_class_unit_example(diagnostics);
    run_unit_saturation_example(diagnostics);

    std::cout << "diagnostics hot-path example completed\n";
    std::cout << "verbose messages = " << counters.verbose << "\n";
    std::cout << "log messages = " << counters.logs << "\n";
    std::cout << "debug failures = " << counters.debug_failures << "\n";
    std::cout << "profile scopes = " << counters.profile_begin
              << " begin / " << counters.profile_end << " end\n";
    std::cout << "profile events = " << counters.profile_events << "\n";
    std::cout << "profile inclusive scope time = "
              << (counters.profile_elapsed_ns / 1000000.0) << " ms\n";
    std::cout << "profile stack issues = " << counters.profile_stack_overflow
              << " overflow / " << counters.profile_unmatched_end
              << " unmatched end\n";
    return 0;
}
