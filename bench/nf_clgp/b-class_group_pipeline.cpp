#include <benchmark/benchmark.h>

#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/fmpz_poly_factor.hpp>
#include <silex/order.hpp>
#include <silex/order_unit.hpp>
#include <silex/unit.hpp>

#include "benchmark_contract.hpp"

#include <algorithm>
#include <cstdint>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct ProfileAggregate {
    std::uint64_t count = 0;
    std::uint64_t events = 0;
    std::chrono::nanoseconds inclusive{0};
    std::chrono::nanoseconds exclusive{0};
};

struct ActiveProfileScope {
    std::string key;
    std::chrono::steady_clock::time_point start;
    std::chrono::nanoseconds child_time{0};
};

struct BenchmarkProfileCollector {
    std::unordered_map<std::string, ProfileAggregate> aggregates;
    std::vector<ActiveProfileScope> stack;

    static bool enabled_from_env() noexcept {
        const char* value = std::getenv("SILEX_BENCH_PROFILE");
        return value != nullptr && std::strcmp(value, "0") != 0 &&
               std::strcmp(value, "off") != 0 &&
               std::strcmp(value, "false") != 0;
    }

    static std::string key(silex::DiagnosticsModule module,
                           const char* label) {
        std::string out = silex::diagnostics_module_name(module);
        out += ":";
        out += label == nullptr ? "<null>" : label;
        return out;
    }

    static bool has_prefix(const std::string& value,
                           const char* prefix) noexcept {
        const std::size_t prefix_len = std::strlen(prefix);
        return value.size() >= prefix_len &&
               value.compare(0, prefix_len, prefix) == 0;
    }

    static std::string counter_name_for_profile_key(
            const std::string& key,
            const char* prefix,
            const char* counter_prefix = "") {
        std::string out = "profile_";
        out += counter_prefix;
        for (std::size_t i = std::strlen(prefix); i < key.size(); ++i) {
            const char ch = key[i];
            const bool alnum = (ch >= 'a' && ch <= 'z') ||
                               (ch >= 'A' && ch <= 'Z') ||
                               (ch >= '0' && ch <= '9');
            out += alnum ? ch : '_';
        }
        return out;
    }

    static void callback(void* user,
                         silex::DiagnosticsModule module,
                         silex::ProfileEvent event,
                         const char*,
                         const char* label) noexcept {
        auto* collector = static_cast<BenchmarkProfileCollector*>(user);
        if (collector == nullptr) {
            return;
        }

        const std::string profile_key = key(module, label);
        if (event == silex::ProfileEvent::begin_scope) {
            collector->stack.push_back(ActiveProfileScope{
                    profile_key, std::chrono::steady_clock::now(), {}});
            return;
        }
        if (event == silex::ProfileEvent::end_scope) {
            if (collector->stack.empty()) {
                return;
            }
            ActiveProfileScope scope = std::move(collector->stack.back());
            collector->stack.pop_back();
            const auto elapsed = std::chrono::steady_clock::now() - scope.start;
            auto& aggregate = collector->aggregates[scope.key];
            ++aggregate.count;
            aggregate.inclusive += elapsed;
            aggregate.exclusive += elapsed - scope.child_time;
            if (!collector->stack.empty()) {
                collector->stack.back().child_time += elapsed;
            }
            return;
        }

        auto& aggregate = collector->aggregates[profile_key];
        ++aggregate.events;
    }

    void configure(silex::DiagnosticsContext& diagnostics) noexcept {
        silex::diagnostics_context_init(diagnostics);
        const auto modules =
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::element) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::ideal) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::fractional_ideal) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::prime_ideal) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::relation) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::class_group) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group);
        silex::diagnostics_set_profiling(
                diagnostics, true, modules, callback, this);
    }

    void print(const char* benchmark_name) const {
        std::cerr << "\nSILEX_BENCH_PROFILE " << benchmark_name << "\n";
        std::vector<const std::pair<const std::string, ProfileAggregate>*> rows;
        rows.reserve(aggregates.size());
        for (const auto& row : aggregates) {
            rows.push_back(&row);
        }
        std::sort(rows.begin(), rows.end(), [](const auto* left, const auto* right) {
            return left->second.inclusive > right->second.inclusive;
        });
        for (const auto* row : rows) {
            const auto& key = row->first;
            const auto& aggregate = row->second;
            const double inclusive_ms =
                    static_cast<double>(aggregate.inclusive.count()) / 1.0e6;
            const double exclusive_ms =
                    static_cast<double>(aggregate.exclusive.count()) / 1.0e6;
            std::cerr << "  " << key
                      << " count=" << aggregate.count
                      << " events=" << aggregate.events
                      << " inclusive_ms=" << inclusive_ms
                      << " exclusive_ms=" << exclusive_ms << "\n";
        }
    }

    void set_stable_scan_counters(benchmark::State& state) const {
        constexpr const char* kPrefix =
                "unit_group:unit_group.stable_scan.";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0 || !has_prefix(row.first, kPrefix)) {
                continue;
            }
            state.counters[counter_name_for_profile_key(row.first, kPrefix)] =
                    static_cast<double>(aggregate.events);
        }
    }

    void set_proof_selector_scan_counters(benchmark::State& state) const {
        constexpr const char* kPrefix =
                "unit_group:unit_group.proof_selector_scan.";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0 || !has_prefix(row.first, kPrefix)) {
                continue;
            }
            state.counters[counter_name_for_profile_key(
                    row.first, kPrefix, "proof_selector_")] =
                    static_cast<double>(aggregate.events);
        }
    }

    void set_relation_path_counters(benchmark::State& state) const {
        constexpr const char* kIntegralPrefix =
                "relation:relation.factor_over_base_integral.";
        constexpr const char* kNonintegralPrefix =
                "relation:relation.factor_over_base_nonintegral.";
        constexpr const char* kSetGeneratorPrefix =
                "relation:relation.set_generator.";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0) {
                continue;
            }
            if (has_prefix(row.first, kIntegralPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kIntegralPrefix)] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kNonintegralPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kNonintegralPrefix, "nonintegral_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kSetGeneratorPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kSetGeneratorPrefix)] =
                        static_cast<double>(aggregate.events);
            }
        }
    }

    void set_partial_relation_counters(benchmark::State& state) const {
        constexpr const char* kDirectLargePrimePrefix =
                "class_group:class_group.partial_relation_direct_large_prime_";
        constexpr const char* kFullFactorPrefix =
                "class_group:class_group.partial_relation_full_factor_";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0) {
                continue;
            }
            if (has_prefix(row.first, kDirectLargePrimePrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kDirectLargePrimePrefix,
                        "partial_direct_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kFullFactorPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kFullFactorPrefix, "partial_full_")] =
                        static_cast<double>(aggregate.events);
            }
        }
    }

    void set_validation_recompute_counters(benchmark::State& state) const {
        constexpr const char* kPrefix =
                "unit_group:unit_group.validation_recompute.";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0 || !has_prefix(row.first, kPrefix)) {
                continue;
            }
            state.counters[counter_name_for_profile_key(
                    row.first, kPrefix, "validation_recompute_")] =
                    static_cast<double>(aggregate.events);
        }
    }

    void set_class_unit_control_flow_counters(benchmark::State& state) const {
        constexpr const char* kRetryPrefix =
                "unit_group:unit_group.candidate_retry.";
        constexpr const char* kOutcomePrefix =
                "unit_group:unit_group.validation_outcome.";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0) {
                continue;
            }
            if (has_prefix(row.first, kRetryPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kRetryPrefix, "candidate_retry_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kOutcomePrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kOutcomePrefix, "validation_outcome_")] =
                        static_cast<double>(aggregate.events);
            }
        }
    }

    void set_saturation_counters(benchmark::State& state) const {
        constexpr const char* kUnitRowPrefix =
                "unit_group:unit_group.saturation_row_";
        constexpr const char* kUnitKernelPrefix =
                "unit_group:unit_group.kernel_row_root_";
        constexpr const char* kUnitLocalPrefix =
                "unit_group:unit_group.saturate_local.";
        constexpr const char* kUnitBoundedPrefix =
                "unit_group:unit_group.saturate_bounded.";
        constexpr const char* kUnitStablePrefix =
                "unit_group:unit_group.stable_relation_saturation.";
        constexpr const char* kUnitCompactPrefix =
                "unit_group:unit_group.compact_unit_";
        constexpr const char* kClassSatPrefix =
                "class_group:class_group.relation_saturation.";
        constexpr const char* kClassProcessPrefix =
                "class_group:class_group.saturation_";
        constexpr const char* kClassLoopPrefix =
                "class_group:class_group.saturate_relation_lattice.";
        for (const auto& row : aggregates) {
            const ProfileAggregate& aggregate = row.second;
            if (aggregate.events == 0) {
                continue;
            }
            if (has_prefix(row.first, kUnitRowPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kUnitRowPrefix, "saturation_row_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kUnitKernelPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kUnitKernelPrefix, "kernel_root_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kUnitLocalPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kUnitLocalPrefix, "saturate_local_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kUnitBoundedPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kUnitBoundedPrefix, "saturate_bounded_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kUnitStablePrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kUnitStablePrefix,
                        "stable_saturation_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kUnitCompactPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kUnitCompactPrefix,
                        "compact_unit_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kClassSatPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kClassSatPrefix, "class_saturation_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kClassProcessPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kClassProcessPrefix,
                        "class_saturation_")] =
                        static_cast<double>(aggregate.events);
            } else if (has_prefix(row.first, kClassLoopPrefix)) {
                state.counters[counter_name_for_profile_key(
                        row.first, kClassLoopPrefix, "class_saturate_")] =
                        static_cast<double>(aggregate.events);
            }
        }
    }

    void set_scope_time_counter(benchmark::State& state,
                                const char* profile_key,
                                const char* counter_key) const {
        const auto pos = aggregates.find(profile_key);
        if (pos == aggregates.end()) {
            return;
        }
        const ProfileAggregate& aggregate = pos->second;
        state.counters[std::string("profile_time_") + counter_key +
                       "_inclusive_ms"] =
                static_cast<double>(aggregate.inclusive.count()) / 1.0e6;
        state.counters[std::string("profile_time_") + counter_key +
                       "_exclusive_ms"] =
                static_cast<double>(aggregate.exclusive.count()) / 1.0e6;
        state.counters[std::string("profile_time_") + counter_key +
                       "_count"] =
                static_cast<double>(aggregate.count);
    }

    void set_class_unit_stage_time_counters(benchmark::State& state) const {
        set_scope_time_counter(
                state, "unit_group:unit_group.compute_with_class_group",
                "class_unit_compute");
        set_scope_time_counter(
                state, "class_group:class_group.compute_candidate",
                "class_group_compute_candidate");
        set_scope_time_counter(
                state, "class_group:class_group.coordinate_candidate",
                "relation_coordinate_candidates");
        set_scope_time_counter(
                state, "class_group:class_group.coordinate_candidate_try_relation",
                "relation_try_generator");
        set_scope_time_counter(
                state, "class_group:class_group.append_relation",
                "relation_append");
        set_scope_time_counter(
                state, "class_group:class_group.append_relation_exact_decision",
                "relation_exact_decision");
        set_scope_time_counter(
                state, "class_group:fmpz_smat.hnf_context.add_row.refine_dependent",
                "hnf_refine_dependent");
        set_scope_time_counter(
                state, "class_group:fmpz_smat.hnf_context.refresh.hnf_transform",
                "hnf_refresh_transform");
        set_scope_time_counter(
                state, "class_group:class_group.try_auto_relation_saturation",
                "relation_saturation");
        set_scope_time_counter(
                state, "class_group:class_group.relation_dlog_column",
                "relation_saturation_dlog");
        set_scope_time_counter(
                state, "unit_group:unit_group.validation_analytic_index_bound",
                "validation_analytic_index_bound");
        set_scope_time_counter(
                state, "unit_group:unit_group.validation_zeta_bf_audit",
                "validation_zeta_bf_audit");
        set_scope_time_counter(
                state, "unit_group:unit_group.zeta_bf.term.prime_decomposition_type",
                "zeta_prime_decomposition");
        set_scope_time_counter(
                state, "unit_group:unit_group.validation_recompute_units",
                "validation_recompute_units");
        set_scope_time_counter(
                state,
                "unit_group:unit_group.set_relation_kernel_units_index_bounded_saturated",
                "unit_kernel_saturation");
        set_scope_time_counter(
                state, "class_group:class_group.partial_relation_full_factor_ideal",
                "partial_relation_factor_ideal");
    }
};

bool proof_diag_enabled_from_env() noexcept {
    const char* value = std::getenv("SILEX_BENCH_PROOF_DIAG");
    return value != nullptr && std::strcmp(value, "0") != 0 &&
           std::strcmp(value, "off") != 0 &&
           std::strcmp(value, "false") != 0;
}

struct ProofRestartLogEvent {
    long ell = -1;
    long index_bound = -1;
    long restart = -1;
    long carried_records = -1;
};

struct BenchmarkProofLogCollector {
    std::vector<long> index_bounds;
    std::vector<ProofRestartLogEvent> restarts;

    static bool parse_long_value(const char* detail,
                                 const char* key,
                                 long& out) noexcept {
        if (detail == nullptr || key == nullptr) {
            return false;
        }
        const std::size_t key_len = std::strlen(key);
        for (const char* pos = std::strstr(detail, key);
             pos != nullptr;
             pos = std::strstr(pos + 1, key)) {
            const bool token_boundary =
                    pos == detail || pos[-1] == ' ' || pos[-1] == '\t';
            const char* value = pos + key_len;
            if (!token_boundary || *value != '=') {
                continue;
            }
            ++value;
            char* end = nullptr;
            const long parsed = std::strtol(value, &end, 10);
            if (end == value) {
                return false;
            }
            out = parsed;
            return true;
        }
        return false;
    }

    static void callback(void* user,
                         silex::DiagnosticsModule module,
                         silex::LogLevel level,
                         const char* function,
                         const char* message,
                         const char* detail) noexcept {
        auto* collector = static_cast<BenchmarkProofLogCollector*>(user);
        if (collector == nullptr || message == nullptr) {
            return;
        }

        if (std::strcmp(message, "unit proof index bound") == 0) {
            long index_bound = -1;
            if (parse_long_value(detail, "index_bound", index_bound)) {
                collector->index_bounds.push_back(index_bound);
            }
        } else if (std::strcmp(message,
                              "unit proof accepted saturation root") == 0) {
            ProofRestartLogEvent event;
            parse_long_value(detail, "ell", event.ell);
            parse_long_value(detail, "index_bound", event.index_bound);
            parse_long_value(detail, "restart", event.restart);
            if (!parse_long_value(detail, "carried_records",
                                  event.carried_records)) {
                parse_long_value(detail, "records_before_reset",
                                 event.carried_records);
            }
            collector->restarts.push_back(event);
        } else if (std::strncmp(message, "validation summary target=",
                                std::strlen("validation summary target=")) ==
                   0) {
            // Fall through to the live diagnostic print below.
        } else if (std::strcmp(message, "raising relation-kernel target") ==
                   0 ||
                   std::strcmp(message,
                               "raising relation-kernel target for missing "
                               "unit rank") == 0 ||
                   std::strcmp(message,
                               "relation-kernel unit proof refinement "
                               "unavailable") == 0 ||
                   std::strcmp(message,
                               "current class-group candidate extension did "
                               "not make progress; restarting candidate") ==
                           0) {
            // Fall through to the live diagnostic print below.
        } else if (std::strcmp(
                           message,
                           "class/unit analytic certification missing "
                           "factor-base generation proof") == 0) {
            // Fall through to the live diagnostic print below.
        } else {
            return;
        }

        std::cerr << "SILEX_BENCH_PROOF_LOG "
                  << silex::diagnostics_module_name(module)
                  << ":" << silex::log_level_name(level) << " "
                  << (function == nullptr ? "<unknown>" : function)
                  << ": " << message;
        if (detail != nullptr) {
            std::cerr << " " << detail;
        }
        std::cerr << "\n";
    }

    void configure(silex::DiagnosticsContext& diagnostics) noexcept {
        silex::diagnostics_set_logging(
                diagnostics, silex::LogLevel::detail,
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group) |
                        silex::diagnostics_module_bit(
                                silex::DiagnosticsModule::class_group),
                callback, this);
    }

    void set_counters(benchmark::State& state) const {
        state.counters["unit_proof_log_index_bound_events"] =
                static_cast<double>(index_bounds.size());
        state.counters["unit_proof_log_restart_events"] =
                static_cast<double>(restarts.size());

        constexpr std::size_t kMaxSlots = 8;
        const std::size_t index_bound_count =
                std::min(kMaxSlots, index_bounds.size());
        for (std::size_t i = 0; i < index_bound_count; ++i) {
            state.counters[
                    "unit_proof_log_index_bound_" + std::to_string(i)] =
                    static_cast<double>(index_bounds[i]);
        }

        const std::size_t restart_count = std::min(kMaxSlots, restarts.size());
        for (std::size_t i = 0; i < restart_count; ++i) {
            const ProofRestartLogEvent& event = restarts[i];
            state.counters[
                    "unit_proof_log_restart_ell_" + std::to_string(i)] =
                    static_cast<double>(event.ell);
            state.counters[
                    "unit_proof_log_restart_bound_" + std::to_string(i)] =
                    static_cast<double>(event.index_bound);
            state.counters[
                    "unit_proof_log_restart_records_" + std::to_string(i)] =
                    static_cast<double>(event.carried_records);
        }
    }
};

inline constexpr slong kBenchComputeSatAuxTarget = 1;
inline constexpr slong kBenchComputeSatAuxStart = 2;
inline constexpr slong kBenchComputeSatAuxMax = 31;
inline constexpr slong kBenchComputeProofAuxMax = 1000;
inline constexpr slong kBenchComputeSatMaxPasses = 2;
void set_degree_one_polynomial(sflint::FmpqPoly& f) noexcept {
    sflint::fmpq_poly_zero(f);
    sflint::fmpq_poly_set_coeff_si(f, 1, 1);
}

void set_cubic_polynomial(sflint::FmpqPoly& f,
                          slong c1,
                          slong c0) noexcept {
    sflint::fmpq_poly_zero(f);
    sflint::fmpq_poly_set_coeff_si(f, 3, 1);
    sflint::fmpq_poly_set_coeff_si(f, 1, c1);
    sflint::fmpq_poly_set_coeff_si(f, 0, c0);
}

void set_cubic_polynomial(sflint::FmpqPoly& f,
                          slong c2,
                          slong c1,
                          slong c0) noexcept {
    sflint::fmpq_poly_zero(f);
    sflint::fmpq_poly_set_coeff_si(f, 3, 1);
    sflint::fmpq_poly_set_coeff_si(f, 2, c2);
    sflint::fmpq_poly_set_coeff_si(f, 1, c1);
    sflint::fmpq_poly_set_coeff_si(f, 0, c0);
}

void set_quartic_polynomial(sflint::FmpqPoly& f,
                            slong c3,
                            slong c2,
                            slong c1,
                            slong c0) noexcept {
    sflint::fmpq_poly_zero(f);
    sflint::fmpq_poly_set_coeff_si(f, 4, 1);
    sflint::fmpq_poly_set_coeff_si(f, 3, c3);
    sflint::fmpq_poly_set_coeff_si(f, 2, c2);
    sflint::fmpq_poly_set_coeff_si(f, 1, c1);
    sflint::fmpq_poly_set_coeff_si(f, 0, c0);
}

void set_quintic_polynomial(sflint::FmpqPoly& f,
                            slong c4,
                            slong c3,
                            slong c2,
                            slong c1,
                            slong c0) noexcept {
    sflint::fmpq_poly_zero(f);
    sflint::fmpq_poly_set_coeff_si(f, 5, 1);
    sflint::fmpq_poly_set_coeff_si(f, 4, c4);
    sflint::fmpq_poly_set_coeff_si(f, 3, c3);
    sflint::fmpq_poly_set_coeff_si(f, 2, c2);
    sflint::fmpq_poly_set_coeff_si(f, 1, c1);
    sflint::fmpq_poly_set_coeff_si(f, 0, c0);
}

void set_sextic_polynomial(sflint::FmpqPoly& f,
                           slong c5,
                           slong c4,
                           slong c3,
                           slong c2,
                           slong c1,
                           slong c0) noexcept {
    sflint::fmpq_poly_zero(f);
    sflint::fmpq_poly_set_coeff_si(f, 6, 1);
    sflint::fmpq_poly_set_coeff_si(f, 5, c5);
    sflint::fmpq_poly_set_coeff_si(f, 4, c4);
    sflint::fmpq_poly_set_coeff_si(f, 3, c3);
    sflint::fmpq_poly_set_coeff_si(f, 2, c2);
    sflint::fmpq_poly_set_coeff_si(f, 1, c1);
    sflint::fmpq_poly_set_coeff_si(f, 0, c0);
}

void set_generic_field_options(silex::ClassGroupCandidateOptions& options,
                               const silex::Order& order) noexcept {
    options.max_candidates = 5000;
    options.max_relations = 500;
    if (order.degree() == 4) {
        options.max_relations = 500;
    }
}

double certification_counter(silex::CertificationMode mode) noexcept {
    return static_cast<double>(static_cast<std::int32_t>(mode));
}

double proof_counter(silex::ProofState state) noexcept {
    return static_cast<double>(static_cast<std::int32_t>(state));
}

void set_class_group_counters(benchmark::State& state,
                              const silex::ClassGroupContext& context,
                              sflint::FmpzConstRef factor_base_bound,
                              bool ok) noexcept {
    state.counters["success"] = ok ? 1.0 : 0.0;
    state.counters["fb_bound"] = sflint::fmpz_get_d(factor_base_bound);
    if (!ok || !context.has_presentation()) {
        return;
    }

    const silex::FactorBase* base = context.factor_base();
    state.counters["fb_len"] = base == nullptr ? 0.0 : base->length();
    state.counters["accepted"] = context.relation_count();
    state.counters["rank"] = context.relation_rank();
    state.counters["kernel_rows"] = context.relation_kernel_unit_count();
    state.counters["skipped_rows"] = context.skipped_dependent_relation_count();
    state.counters["search_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::Search);
    state.counters["random_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::RandomProduct);
    state.counters["large_prime_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::LargePrimeMatch);
    state.counters["classgen_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::ClassGenerator);

    sflint::Fmpz order;
    if (context.order(sflint::FmpzRef(order))) {
        state.counters["class_order"] =
                sflint::fmpz_get_d(sflint::FmpzConstRef(order));
    }
}

bool valid_class_group_result(
        const silex::ClassGroupContext& context) noexcept {
    if (!context.has_presentation()) {
        return false;
    }
    sflint::Fmpz class_order;
    return context.order(sflint::FmpzRef(class_order)) &&
            sflint::fmpz_cmp_ui(sflint::FmpzConstRef(class_order), 0) > 0;
}

struct ExactClassGroupExpectation {
    slong class_order;
    const slong* invariants;
    slong invariant_count;
};

struct ExactClassUnitExpectation {
    ExactClassGroupExpectation class_group;
    slong unit_rank;
    slong torsion_order;
    double regulator_lower;
    double regulator_upper;
};

inline constexpr slong kC2Invariants[] = {2};
inline constexpr slong kC4Invariants[] = {4};
inline constexpr slong kC5Invariants[] = {5};
inline constexpr ExactClassGroupExpectation kTrivialClassGroup{
        1, nullptr, 0};
inline constexpr ExactClassGroupExpectation kC2ClassGroup{
        2, kC2Invariants, 1};
inline constexpr ExactClassGroupExpectation kC4ClassGroup{
        4, kC4Invariants, 1};
inline constexpr ExactClassGroupExpectation kC5ClassGroup{
        5, kC5Invariants, 1};

// The torsion orders and regulator intervals were computed independently with
// GP 2.17.4 bnf data, with the degree-greater-than-one fields certified by its
// exact certification routine.  Every positive-rank interval is much narrower
// than a factor of two, so it rejects the regulator of any proper finite-index
// subgroup.  Rank-zero regulator contracts are exact.
inline constexpr ExactClassUnitExpectation kDegreeOneClassUnit{
        kTrivialClassGroup, 0, 2, 1.0, 1.0};
inline constexpr ExactClassUnitExpectation kRealQuadraticClassUnit{
        kTrivialClassGroup, 1, 2, 0.48121182505, 0.48121182507};
inline constexpr ExactClassUnitExpectation kImaginaryQuadraticClassUnit{
        kC5ClassGroup, 0, 2, 1.0, 1.0};
inline constexpr ExactClassUnitExpectation kCubicTrivialClassUnit{
        kTrivialClassGroup, 1, 2, 1.34737734832, 1.34737734834};
inline constexpr ExactClassUnitExpectation kCubicNontrivialClassUnit{
        kC2ClassGroup, 1, 2, 2.35861081458, 2.35861081461};
inline constexpr ExactClassUnitExpectation kQuarticCyclotomicClassUnit{
        kTrivialClassGroup, 1, 8, 1.76274717403, 1.76274717405};
inline constexpr ExactClassUnitExpectation kQuarticNoncyclotomicClassUnit{
        kTrivialClassGroup, 1, 2, 1.29274426911, 1.29274426914};
inline constexpr ExactClassUnitExpectation kQuarticDisc70640ClassUnit{
        kTrivialClassGroup, 2, 2, 14.4609875438, 14.4609875440};
inline constexpr ExactClassUnitExpectation kQuarticDisc223479ClassUnit{
        kTrivialClassGroup, 2, 2, 49.8132070222, 49.8132070225};
inline constexpr ExactClassUnitExpectation kQuarticDisc35019ClassUnit{
        kTrivialClassGroup, 2, 2, 2.58652290111, 2.58652290114};
inline constexpr ExactClassUnitExpectation kQuarticDisc1412343ClassUnit{
        kC4ClassGroup, 2, 2, 28.9150416061, 28.9150416063};
inline constexpr ExactClassUnitExpectation kQuarticSmallClassUnit{
        kTrivialClassGroup, 2, 2, 0.37819933245, 0.37819933247};
inline constexpr ExactClassUnitExpectation kQuinticSmallClassUnit{
        kTrivialClassGroup, 2, 2, 0.43234387881, 0.43234387884};
inline constexpr ExactClassUnitExpectation kQuinticDisc11119ClassUnit{
        kTrivialClassGroup, 3, 2, 1.19655267603, 1.19655267606};
inline constexpr ExactClassUnitExpectation kQuinticDisc401370255ClassUnit{
        kTrivialClassGroup, 3, 2, 734.65403299, 734.65403302};
inline constexpr ExactClassUnitExpectation kQuinticDisc57895ClassUnit{
        kTrivialClassGroup, 3, 2, 4.43153392698, 4.43153392702};
inline constexpr ExactClassUnitExpectation kSexticSmallClassUnit{
        kTrivialClassGroup, 3, 2, 0.74063147262, 0.74063147264};
inline constexpr ExactClassUnitExpectation kRandomCubicClassUnit{
        kTrivialClassGroup, 1, 2, 8.29429929651, 8.29429929655};

void initialize_exact_class_group_counters(
        benchmark::State& state,
        ExactClassGroupExpectation expected) noexcept {
    state.counters["expected_class_order"] = expected.class_order;
    state.counters["class_order_exact"] = 0.0;
    state.counters["expected_class_invariant_count"] =
            expected.invariant_count;
    state.counters["class_invariant_count_exact"] = 0.0;
    state.counters["class_invariants_exact"] = 0.0;
    for (slong i = 0; i < expected.invariant_count; ++i) {
        if (expected.invariants != nullptr) {
            state.counters[
                    "expected_class_invariant_" + std::to_string(i)] =
                    expected.invariants[i];
        }
    }
}

bool validate_exact_class_group_result(
        benchmark::State& state,
        const silex::ClassGroupContext& context,
        ExactClassGroupExpectation expected) noexcept {
    initialize_exact_class_group_counters(state, expected);
    if (!context.has_presentation()) {
        return false;
    }

    sflint::Fmpz class_order;
    if (!context.order(sflint::FmpzRef(class_order))) {
        return false;
    }
    state.counters["validation_class_order"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(class_order));
    const bool class_order_exact = sflint::fmpz_equal_si(
            sflint::FmpzConstRef(class_order), expected.class_order);
    state.counters["class_order_exact"] = class_order_exact ? 1.0 : 0.0;

    const slong actual_invariant_count = context.invariant_count();
    state.counters["validation_class_invariant_count"] =
            actual_invariant_count;
    const bool invariant_count_exact =
            actual_invariant_count == expected.invariant_count;
    state.counters["class_invariant_count_exact"] =
            invariant_count_exact ? 1.0 : 0.0;

    bool invariants_exact =
            invariant_count_exact && expected.invariant_count >= 0 &&
            (expected.invariant_count == 0 || expected.invariants != nullptr);
    for (slong i = 0; i < expected.invariant_count; ++i) {
        sflint::Fmpz invariant;
        const bool available = context.invariant(sflint::FmpzRef(invariant), i);
        if (available) {
            state.counters[
                    "validation_class_invariant_" + std::to_string(i)] =
                    sflint::fmpz_get_d(sflint::FmpzConstRef(invariant));
        }
        if (!available || expected.invariants == nullptr ||
            !sflint::fmpz_equal_si(
                    sflint::FmpzConstRef(invariant),
                    expected.invariants[i])) {
            invariants_exact = false;
        }
    }
    state.counters["class_invariants_exact"] =
            invariants_exact ? 1.0 : 0.0;
    return class_order_exact && invariant_count_exact && invariants_exact;
}

bool valid_class_unit_result(
        const silex::ClassGroupContext& class_group,
        const silex::OrderUnitGroup& units,
        silex::CertificationMode requested_certification) noexcept {
    if (!valid_class_group_result(class_group) || !units.is_set() ||
        units.free_rank() < 0 || class_group.parent() == nullptr ||
        units.parent() == nullptr ||
        !silex::same_order_parent(class_group.parent(), units.parent())) {
        return false;
    }
    if (requested_certification == silex::CertificationMode::proven) {
        return class_group.certification_status() ==
                       silex::CertificationMode::proven &&
                units.certification_status() ==
                       silex::CertificationMode::proven;
    }
    return true;
}

bool element_has_exact_order(const silex::Element& element,
                             slong expected_order) noexcept {
    const silex::NumberField* parent = element.parent();
    if (parent == nullptr || expected_order <= 0) {
        return false;
    }

    silex::Element power(*parent);
    if (!power.is_defined()) {
        return false;
    }
    for (slong exponent = 1; exponent <= expected_order; ++exponent) {
        sflint::Fmpz exponent_value;
        sflint::fmpz_set_si(sflint::FmpzRef(exponent_value), exponent);
        if (!power.pow_fmpz(
                    element, sflint::FmpzConstRef(exponent_value)) ||
            (exponent < expected_order && power.equal_si(1))) {
            return false;
        }
    }
    return power.equal_si(1);
}

bool element_has_unit_norm(const silex::Element& element) noexcept {
    sflint::Fmpq norm;
    return element.norm(sflint::FmpqRef(norm)) &&
            (sflint::fmpq_equal_si(sflint::FmpqConstRef(norm), 1) ||
             sflint::fmpq_equal_si(sflint::FmpqConstRef(norm), -1));
}

bool regulator_in_expected_interval(
        const sflint::Arb& regulator,
        ExactClassUnitExpectation expected) noexcept {
    if (expected.unit_rank == 0) {
        return expected.regulator_lower == 1.0 &&
                expected.regulator_upper == 1.0 &&
                sflint::arb_is_one(regulator);
    }
    if (!(expected.regulator_lower > 0.0) ||
        !(expected.regulator_lower < expected.regulator_upper) ||
        !(expected.regulator_upper < 2.0 * expected.regulator_lower)) {
        return false;
    }

    sflint::Arb lower;
    sflint::Arb upper;
    sflint::arb_set_d(lower, expected.regulator_lower);
    sflint::arb_set_d(upper, expected.regulator_upper);
    return sflint::arb_gt(regulator, lower) &&
            sflint::arb_lt(regulator, upper);
}

void set_regulator_bounds_counters(benchmark::State& state,
                                   const char* prefix,
                                   const sflint::Arb& regulator) {
    sflint::Arf lower;
    sflint::Arf upper;
    sflint::arb_get_lbound_arf(lower, regulator, 128);
    sflint::arb_get_ubound_arf(upper, regulator, 128);
    state.counters[std::string(prefix) + "_lower"] =
            sflint::arf_get_d(lower, ARF_RND_DOWN);
    state.counters[std::string(prefix) + "_upper"] =
            sflint::arf_get_d(upper, ARF_RND_UP);
}

void initialize_exact_class_unit_counters(
        benchmark::State& state,
        ExactClassUnitExpectation expected) noexcept {
    initialize_exact_class_group_counters(state, expected.class_group);
    state.counters["expected_unit_rank"] = expected.unit_rank;
    state.counters["unit_rank_exact"] = 0.0;
    state.counters["expected_unit_rank_matches_signature"] = 0.0;
    state.counters["expected_torsion_order"] = expected.torsion_order;
    state.counters["torsion_order_exact"] = 0.0;
    state.counters["torsion_generator_available"] = 0.0;
    state.counters["torsion_generator_in_order"] = 0.0;
    state.counters["torsion_generator_unit_norm"] = 0.0;
    state.counters["torsion_generator_exact_order"] = 0.0;
    state.counters["free_generators_available"] = 0.0;
    state.counters["free_generators_in_order"] = 0.0;
    state.counters["free_generators_unit_norm"] = 0.0;
    state.counters["stored_regulator_available"] = 0.0;
    state.counters["recomputed_regulator_available"] = 0.0;
    state.counters["regulators_overlap"] = 0.0;
    state.counters["stored_regulator_in_expected_interval"] = 0.0;
    state.counters["recomputed_regulator_in_expected_interval"] = 0.0;
    state.counters["expected_regulator_lower"] = expected.regulator_lower;
    state.counters["expected_regulator_upper"] = expected.regulator_upper;
    state.counters["unit_output_exact"] = 0.0;
}

bool validate_exact_unit_output(
        benchmark::State& state,
        const silex::Order& order,
        const silex::OrderUnitGroup& units,
        ExactClassUnitExpectation expected) {
    const silex::NumberField* const field = order.parent();
    if (field == nullptr || expected.unit_rank < 0 ||
        expected.torsion_order <= 0) {
        return false;
    }

    sflint::Fmpz torsion_order;
    const bool torsion_order_available =
            units.torsion_order(sflint::FmpzRef(torsion_order));
    if (torsion_order_available) {
        state.counters["validation_torsion_order"] =
                sflint::fmpz_get_d(sflint::FmpzConstRef(torsion_order));
    }
    const bool torsion_order_exact = torsion_order_available &&
            sflint::fmpz_equal_si(
                    sflint::FmpzConstRef(torsion_order),
                    expected.torsion_order);
    state.counters["torsion_order_exact"] =
            torsion_order_exact ? 1.0 : 0.0;

    silex::OrderElement torsion_generator(order);
    const bool torsion_generator_available =
            units.torsion_generator(torsion_generator);
    state.counters["torsion_generator_available"] =
            torsion_generator_available ? 1.0 : 0.0;
    silex::Element torsion_value(*field);
    const bool torsion_value_available = torsion_generator_available &&
            torsion_generator.get_element(torsion_value);
    const bool torsion_in_order =
            torsion_value_available && order.contains(torsion_value);
    const bool torsion_unit_norm =
            torsion_value_available && element_has_unit_norm(torsion_value);
    const bool torsion_exact_order = torsion_value_available &&
            element_has_exact_order(torsion_value, expected.torsion_order);
    state.counters["torsion_generator_in_order"] =
            torsion_in_order ? 1.0 : 0.0;
    state.counters["torsion_generator_unit_norm"] =
            torsion_unit_norm ? 1.0 : 0.0;
    state.counters["torsion_generator_exact_order"] =
            torsion_exact_order ? 1.0 : 0.0;

    std::vector<silex::Element> evaluated_generators;
    evaluated_generators.reserve(
            static_cast<std::size_t>(expected.unit_rank));
    bool free_generators_available = true;
    bool free_generators_in_order = true;
    bool free_generators_unit_norm = true;
    for (slong i = 0; i < expected.unit_rank; ++i) {
        silex::FactoredElement generator(*field);
        silex::Element value(*field);
        if (!units.free_generator(generator, i) ||
            !generator.evaluate(value)) {
            free_generators_available = false;
            break;
        }
        free_generators_in_order =
                free_generators_in_order && order.contains(value);
        free_generators_unit_norm =
                free_generators_unit_norm && element_has_unit_norm(value);
        evaluated_generators.push_back(std::move(value));
    }
    state.counters["validation_free_generator_count"] =
            static_cast<double>(evaluated_generators.size());
    state.counters["free_generators_available"] =
            free_generators_available ? 1.0 : 0.0;
    state.counters["free_generators_in_order"] =
            free_generators_in_order ? 1.0 : 0.0;
    state.counters["free_generators_unit_norm"] =
            free_generators_unit_norm ? 1.0 : 0.0;

    sflint::Arb stored_regulator;
    const bool stored_regulator_available =
            units.regulator(sflint::ArbRef(stored_regulator));
    state.counters["stored_regulator_available"] =
            stored_regulator_available ? 1.0 : 0.0;
    if (stored_regulator_available) {
        set_regulator_bounds_counters(
                state, "validation_stored_regulator", stored_regulator);
    }

    sflint::Arb recomputed_regulator;
    bool recomputed_regulator_available = false;
    if (free_generators_available &&
        evaluated_generators.size() ==
                static_cast<std::size_t>(expected.unit_rank)) {
        silex::EmbeddingContext embeddings(*field);
        recomputed_regulator_available = silex::unit_regulator(
                sflint::ArbRef(recomputed_regulator), embeddings,
                silex::ElementSpan(evaluated_generators.data(),
                                   evaluated_generators.size()),
                128);
    }
    state.counters["recomputed_regulator_available"] =
            recomputed_regulator_available ? 1.0 : 0.0;
    if (recomputed_regulator_available) {
        set_regulator_bounds_counters(
                state, "validation_recomputed_regulator",
                recomputed_regulator);
    }

    const bool regulators_overlap = stored_regulator_available &&
            recomputed_regulator_available &&
            sflint::arb_overlaps(stored_regulator, recomputed_regulator);
    const bool stored_regulator_expected = stored_regulator_available &&
            regulator_in_expected_interval(stored_regulator, expected);
    const bool recomputed_regulator_expected =
            recomputed_regulator_available &&
            regulator_in_expected_interval(recomputed_regulator, expected);
    state.counters["regulators_overlap"] =
            regulators_overlap ? 1.0 : 0.0;
    state.counters["stored_regulator_in_expected_interval"] =
            stored_regulator_expected ? 1.0 : 0.0;
    state.counters["recomputed_regulator_in_expected_interval"] =
            recomputed_regulator_expected ? 1.0 : 0.0;

    const bool exact = torsion_order_exact &&
            torsion_generator_available && torsion_in_order &&
            torsion_unit_norm && torsion_exact_order &&
            free_generators_available && free_generators_in_order &&
            free_generators_unit_norm && stored_regulator_available &&
            recomputed_regulator_available && regulators_overlap &&
            stored_regulator_expected && recomputed_regulator_expected;
    state.counters["unit_output_exact"] = exact ? 1.0 : 0.0;
    return exact;
}

bool validate_exact_class_unit_result(
        benchmark::State& state,
        const silex::Order& order,
        const silex::ClassGroupContext& class_group,
        const silex::OrderUnitGroup& units,
        silex::CertificationMode requested_certification,
        ExactClassUnitExpectation expected) noexcept {
    initialize_exact_class_unit_counters(state, expected);

    const bool class_order_exact = validate_exact_class_group_result(
            state, class_group, expected.class_group);
    state.counters["validation_unit_rank"] = units.free_rank();
    const bool unit_rank_exact = units.free_rank() == expected.unit_rank;
    state.counters["unit_rank_exact"] = unit_rank_exact ? 1.0 : 0.0;

    slong signature_rank = -1;
    const bool signature_available = order.parent() != nullptr &&
            silex::unit_rank(signature_rank, *order.parent());
    if (signature_available) {
        state.counters["signature_unit_rank"] = signature_rank;
    }
    const bool expected_matches_signature =
            signature_available && signature_rank == expected.unit_rank;
    state.counters["expected_unit_rank_matches_signature"] =
            expected_matches_signature ? 1.0 : 0.0;

    const bool unit_output_exact = validate_exact_unit_output(
            state, order, units, expected);

    return valid_class_unit_result(
                   class_group, units, requested_certification) &&
            class_order_exact && unit_rank_exact &&
            expected_matches_signature && unit_output_exact;
}

void set_partial_class_group_counters(
        benchmark::State& state,
        const silex::ClassGroupContext& context,
        sflint::FmpzConstRef factor_base_bound) noexcept {
    state.counters["fb_bound"] = sflint::fmpz_get_d(factor_base_bound);
    state.counters["has_factor_base"] = context.has_factor_base() ? 1.0 : 0.0;
    state.counters["has_presentation"] = context.has_presentation() ? 1.0
                                                                    : 0.0;

    const silex::FactorBase* base = context.factor_base();
    state.counters["fb_len"] = base == nullptr ? 0.0 : base->length();
    state.counters["accepted"] = context.relation_count();
    state.counters["rank"] = context.relation_rank();
    state.counters["skipped_rows"] = context.skipped_dependent_relation_count();
    state.counters["search_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::Search);
    state.counters["random_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::RandomProduct);
    state.counters["large_prime_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::LargePrimeMatch);
    state.counters["classgen_sources"] = context.relation_source_count(
            silex::ClassGroupRelationSource::ClassGenerator);
}

enum class ClassUnitProofSource {
    none = 0,
    trivial_or_special = 1,
    belabas_friedman_analytic = 2,
    unavailable = 3,
};

enum class ClassUnitFailureReason {
    none = 0,
    reducible_polynomial = 1,
    field_or_order_unavailable = 2,
    options_unavailable = 3,
    compute_unavailable = 4,
};

enum class ClassUnitMatrixMode {
    class_candidate = 0,
    pair_proven = 1,
};

double class_unit_proof_source(
        const silex::ClassGroupContext& class_group,
        const silex::OrderUnitGroup& units,
        silex::CertificationMode requested_certification) noexcept {
    if (requested_certification != silex::CertificationMode::proven) {
        return static_cast<double>(ClassUnitProofSource::none);
    }
    if (class_group.certification_status() !=
                silex::CertificationMode::proven ||
        units.certification_status() != silex::CertificationMode::proven) {
        return static_cast<double>(ClassUnitProofSource::unavailable);
    }
    if (class_group.zeta_bf_proof_status() == silex::ProofState::verified) {
        return static_cast<double>(
                ClassUnitProofSource::belabas_friedman_analytic);
    }
    return static_cast<double>(ClassUnitProofSource::trivial_or_special);
}

void set_unit_proof_selector_replay_counters(
        benchmark::State& state,
        const silex::OrderUnitGroup& units,
        sflint::FmpzConstRef ell,
        sflint::FmpzConstRef aux_prime_bound) {
    state.counters["unit_proof_selector_ell"] = sflint::fmpz_get_d(ell);
    state.counters["unit_proof_selector_bound"] =
            sflint::fmpz_get_d(aux_prime_bound);
    state.PauseTiming();
    silex::PrimeIdealList primes;
    bool certified = false;
    sflint::FmpzMat kernel(0, units.free_rank());
    const bool ok = units.select_saturation_proof_primes(
            primes, certified, kernel, ell, aux_prime_bound);
    const slong target_rank = sflint::fmpz_mat_ncols(kernel);
    const slong kernel_rank = sflint::fmpz_mat_nrows(kernel);
    state.ResumeTiming();

    state.counters["unit_proof_selector_ok"] = ok ? 1.0 : 0.0;
    if (!ok) {
        return;
    }

    state.counters["unit_proof_selector_target_rank"] = target_rank;
    state.counters["unit_proof_selector_kernel_rank"] = kernel_rank;
    state.counters["unit_proof_selector_achieved_rank"] =
            target_rank >= kernel_rank ? target_rank - kernel_rank : 0;
    state.counters["unit_proof_selector_primes"] = primes.size();
    state.counters["unit_proof_selector_certified"] =
            certified ? 1.0 : 0.0;
    state.counters["unit_proof_selector_cap_hit"] =
            certified ? 0.0 : 1.0;
}

void set_unit_proof_record_counters(benchmark::State& state,
                                    const silex::OrderUnitGroup& units,
                                    bool replay_selector) {
    const slong count = units.unit_proof_record_count();
    slong verified = 0;
    slong unavailable = 0;
    bool have_last_unavailable = false;
    sflint::Fmpz last_unavailable_ell;
    sflint::Fmpz last_unavailable_aux;
    slong last_unavailable_primes = 0;
    bool last_unavailable_changed = false;

    for (slong i = 0; i < count; ++i) {
        sflint::Fmpz ell;
        sflint::Fmpz aux_prime_bound;
        silex::ProofState status = silex::ProofState::not_checked;
        slong local_primes = 0;
        bool changed = false;
        if (!units.unit_proof_record(
                    sflint::FmpzRef(ell), status,
                    sflint::FmpzRef(aux_prime_bound), local_primes, changed,
                    i)) {
            continue;
        }

        if (status == silex::ProofState::verified) {
            ++verified;
        } else if (status == silex::ProofState::unavailable) {
            ++unavailable;
            have_last_unavailable = true;
            sflint::fmpz_set(sflint::FmpzRef(last_unavailable_ell),
                             sflint::FmpzConstRef(ell));
            sflint::fmpz_set(sflint::FmpzRef(last_unavailable_aux),
                             sflint::FmpzConstRef(aux_prime_bound));
            last_unavailable_primes = local_primes;
            last_unavailable_changed = changed;
        }
    }

    state.counters["unit_proof_records"] = count;
    state.counters["unit_proof_records_verified"] = verified;
    state.counters["unit_proof_records_unavailable"] = unavailable;
    if (!have_last_unavailable) {
        return;
    }

    state.counters["unit_proof_last_unavail_ell"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(last_unavailable_ell));
    state.counters["unit_proof_last_unavail_aux"] =
            sflint::fmpz_get_d(sflint::FmpzConstRef(last_unavailable_aux));
    state.counters["unit_proof_last_unavail_primes"] =
            last_unavailable_primes;
    state.counters["unit_proof_last_unavail_changed"] =
            last_unavailable_changed ? 1.0 : 0.0;

    if (replay_selector) {
        set_unit_proof_selector_replay_counters(
                state, units, sflint::FmpzConstRef(last_unavailable_ell),
                sflint::FmpzConstRef(last_unavailable_aux));
    }
}

void set_failed_unit_proof_replay_counters(
        benchmark::State& state,
        const silex::Order& order,
        sflint::FmpzConstRef factor_base_bound,
        const silex::ClassGroupComputeOptions& options,
        slong precision) {
    if (!proof_diag_enabled_from_env() ||
        options.requested_certification != silex::CertificationMode::proven ||
        !order.is_defined() || order.parent() == nullptr) {
        return;
    }

    state.PauseTiming();
    state.counters["unit_proof_diag_replay_attempted"] = 1.0;

    silex::ClassGroupCandidateOptions local_options;
    local_options.max_candidates = options.max_candidates;
    local_options.max_relations = options.max_relations;
    local_options.diagnostics = options.diagnostics;

    silex::ClassGroupContext replay_class_group;
    replay_class_group.set_diagnostics(options.diagnostics);
    const bool class_ok = replay_class_group.compute_candidate(
            order, factor_base_bound, local_options);
    state.counters["unit_proof_diag_class_ok"] = class_ok ? 1.0 : 0.0;
    if (!class_ok) {
        state.ResumeTiming();
        return;
    }
    state.counters["unit_proof_diag_relations"] =
            replay_class_group.relation_count();
    state.counters["unit_proof_diag_kernel_units"] =
            replay_class_group.relation_kernel_unit_count();

    silex::EmbeddingContext embeddings(*order.parent());
    silex::OrderUnitGroup refined(order);
    refined.set_diagnostics(options.diagnostics);
    sflint::Fmpz aux_start;
    sflint::Fmpz aux_max;
    sflint::fmpz_set_si(sflint::FmpzRef(aux_start),
                        kBenchComputeSatAuxStart);
    sflint::fmpz_set_si(sflint::FmpzRef(aux_max), kBenchComputeSatAuxMax);
    bool unit_changed = false;
    bool unit_stable = false;
    const bool units_ok =
            refined.set_relation_kernel_units_index_bounded_saturated(
                    unit_changed, unit_stable, order, replay_class_group,
                    embeddings, precision, precision,
                    kBenchComputeSatAuxTarget, sflint::FmpzConstRef(aux_start),
                    sflint::FmpzConstRef(aux_max),
                    kBenchComputeSatMaxPasses);
    state.counters["unit_proof_diag_units_ok"] = units_ok ? 1.0 : 0.0;
    state.counters["unit_proof_diag_units_changed"] =
            unit_changed ? 1.0 : 0.0;
    state.counters["unit_proof_diag_units_stable"] =
            unit_stable ? 1.0 : 0.0;
    if (!units_ok) {
        state.ResumeTiming();
        return;
    }
    state.counters["unit_proof_diag_unit_rank"] = refined.free_rank();

    silex::OrderUnitGroup proven(order);
    proven.set_diagnostics(options.diagnostics);
    sflint::Fmpz proof_aux_bound;
    sflint::fmpz_set_si(sflint::FmpzRef(proof_aux_bound),
                        kBenchComputeProofAuxMax);
    silex::ProofState status = silex::ProofState::not_checked;
    bool proof_changed = false;
    const bool proof_ok = proven.prove_index_bound(
            status, proof_changed, refined, kBenchComputeSatAuxTarget,
            sflint::FmpzConstRef(proof_aux_bound),
            kBenchComputeSatMaxPasses, embeddings, precision);
    state.counters["unit_proof_diag_proof_ok"] = proof_ok ? 1.0 : 0.0;
    state.counters["unit_proof_diag_status"] = proof_counter(status);
    state.counters["unit_proof_diag_proof_changed"] =
            proof_changed ? 1.0 : 0.0;
    state.ResumeTiming();

    if (proof_ok) {
        set_unit_proof_record_counters(state, proven, true);
    }
}

void set_class_unit_counters(
        benchmark::State& state,
        const silex::Order& order,
        const silex::ClassGroupContext& class_group,
        const silex::OrderUnitGroup& units,
        const silex::ClassGroupComputeOptions& options,
        bool ok) noexcept {
    state.counters["success"] = ok ? 1.0 : 0.0;
    if (state.counters.find("failure_reason") == state.counters.end()) {
        state.counters["failure_reason"] =
                static_cast<double>(ClassUnitFailureReason::none);
    }
    state.counters["degree"] = order.degree();
    state.counters["requested_cert"] =
            certification_counter(options.requested_certification);
    state.counters["bf_max_cutoff"] = options.zeta_bf_max_cutoff;

    if (!class_group.has_presentation() || !units.is_set()) {
        return;
    }

    sflint::Fmpz class_order;
    if (class_group.order(sflint::FmpzRef(class_order))) {
        state.counters["class_order"] =
                sflint::fmpz_get_d(sflint::FmpzConstRef(class_order));
    }
    state.counters["class_cert"] =
            certification_counter(class_group.certification_status());
    state.counters["unit_cert"] =
            certification_counter(units.certification_status());
    state.counters["unit_rank"] = units.free_rank();
    state.counters["kernel_units"] =
            class_group.relation_kernel_unit_count();
    state.counters["fb_checked"] = proof_counter(
            class_group.factor_base_generation_checked_status());
    state.counters["relation_saturation"] =
            proof_counter(class_group.relation_saturation_status());
    state.counters["unit_proof"] =
            proof_counter(class_group.unit_proof_status());
    state.counters["regulator_proof"] =
            proof_counter(class_group.regulator_proof_status());
    state.counters["analytic_hR"] =
            proof_counter(class_group.analytic_class_regulator_status());
    state.counters["zeta_bf"] =
            proof_counter(class_group.zeta_bf_proof_status());
    state.counters["proof_source"] = class_unit_proof_source(
            class_group, units, options.requested_certification);
    if (proof_diag_enabled_from_env()) {
        set_unit_proof_record_counters(state, units, true);
    }

    ulong cutoff = 0;
    ulong max_cutoff = 0;
    slong requested_precision = 0;
    slong work_precision = 0;
    sflint::Arb error_bound;
    if (class_group.zeta_bf_proof_record(
                cutoff, max_cutoff, requested_precision, work_precision,
                sflint::ArbRef(error_bound))) {
        state.counters["bf_cutoff"] = cutoff;
        state.counters["bf_record_max"] = max_cutoff;
        state.counters["bf_requested_prec"] = requested_precision;
        state.counters["bf_work_prec"] = work_precision;
    }
}

struct FieldOrderSetup {
    silex::NumberField field;
    silex::Order equation_order;
    silex::Order maximal_order;

    bool define(const sflint::FmpqPoly& polynomial) noexcept {
        if (!field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial))) {
            return false;
        }
        if (!equation_order.define_equation_order(field)) {
            return false;
        }
        if (!maximal_order.define(field)) {
            return false;
        }
        return maximal_order.maximal_order(equation_order);
    }

    bool define_quadratic(sflint::FmpzConstRef radicand) noexcept {
        if (!field.define_quadratic(radicand)) {
            return false;
        }
        if (!equation_order.define_equation_order(field)) {
            return false;
        }
        if (!maximal_order.define(field)) {
            return false;
        }
        return maximal_order.maximal_order(equation_order);
    }
};

void BM_class_group_candidate_from_poly_with_bound_scale(
        benchmark::State& state,
        const sflint::FmpqPoly& polynomial,
        ulong bound_scale,
        ExactClassGroupExpectation expected) {
    silex::bench_contract::initialize(state);
    initialize_exact_class_group_counters(state, expected);
    FieldOrderSetup setup;
    if (!setup.define(polynomial)) {
        silex::bench_contract::fail(
                state, "field/order setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    sflint::Fmpz factor_base_bound;
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), setup.maximal_order)) {
        silex::bench_contract::fail(
                state, "class-group factor-base bound unavailable",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    if (bound_scale > 1) {
        sflint::fmpz_mul_ui(sflint::FmpzRef(factor_base_bound),
                            sflint::FmpzConstRef(factor_base_bound),
                            bound_scale);
    }

    silex::ClassGroupCandidateOptions options;
    set_generic_field_options(options, setup.maximal_order);

    BenchmarkProfileCollector profile_collector;
    silex::DiagnosticsContext diagnostics;
    const bool profile_enabled =
            BenchmarkProfileCollector::enabled_from_env();
    if (profile_enabled) {
        profile_collector.configure(diagnostics);
        options.diagnostics = &diagnostics;
    }

    bool completed = false;
    bool failed = false;
    for (auto _ : state) {
        silex::ClassGroupContext context;
        const bool ok = context.compute_candidate(
                setup.maximal_order, sflint::FmpzConstRef(factor_base_bound),
                options);
        int ok_value = ok ? 1 : 0;
        benchmark::DoNotOptimize(ok_value);
        set_class_group_counters(
                state, context, sflint::FmpzConstRef(factor_base_bound), ok);
        if (!ok) {
            failed = true;
            silex::bench_contract::fail(
                    state, "class-group candidate compute failed",
                    silex::bench_contract::FailureReason::operation);
            break;
        }
        completed = true;
    }

    if (completed && !failed) {
        silex::ClassGroupCandidateOptions validation_options = options;
        validation_options.diagnostics = nullptr;
        silex::ClassGroupContext validation_context;
        const bool validation_ok = validation_context.compute_candidate(
                setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), validation_options);
        const bool invariant_ok = validation_ok &&
                validate_exact_class_group_result(
                        state, validation_context, expected);
        if (!invariant_ok) {
            failed = true;
            silex::bench_contract::fail(
                    state,
                    validation_ok
                            ? "class-group candidate invariant failed"
                            : "class-group candidate validation rerun failed",
                    validation_ok
                            ? silex::bench_contract::FailureReason::invariant
                            : silex::bench_contract::FailureReason::operation);
        }
    }

    benchmark::ClobberMemory();
    if (profile_enabled) {
        profile_collector.set_stable_scan_counters(state);
        profile_collector.set_proof_selector_scan_counters(state);
        profile_collector.set_relation_path_counters(state);
        profile_collector.set_partial_relation_counters(state);
        profile_collector.set_validation_recompute_counters(state);
        profile_collector.set_class_unit_control_flow_counters(state);
        profile_collector.set_saturation_counters(state);
        profile_collector.set_class_unit_stage_time_counters(state);
        profile_collector.print(state.name().c_str());
    }
    if (completed && !failed) {
        silex::bench_contract::succeed(state);
    }
}

void BM_class_group_candidate_from_poly(
        benchmark::State& state,
        const sflint::FmpqPoly& polynomial,
        ExactClassGroupExpectation expected) {
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 1, expected);
}

// Exact fixture contracts below are mirrored from
// test/data/class_unit_fields.json and test/t-class-unit-matrix.cpp where
// those fixtures occur.  Benchmark-only class numbers and cyclic factors were
// independently certified with an external computer algebra system when
// these contracts were introduced.  Paired unit ranks are also checked
// against Silex's signature-based unit_rank contract during the untimed
// validation rerun.

void BM_class_group_trivial_cubic_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, 0, -2);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_nontrivial_cubic_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -2, -5);
    BM_class_group_candidate_from_poly(state, polynomial, kC2ClassGroup);
}

void BM_class_group_nontrivial_cubic_candidate_fb2x(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -2, -5);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 2, kC2ClassGroup);
}

void BM_class_group_nontrivial_cubic_candidate_fb4x(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -2, -5);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 4, kC2ClassGroup);
}

void BM_class_group_mixed_cubic_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -1, -2, 1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_cubic_3_1_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -3, 1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_cubic_minus3_minus1_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -3, -1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_small_quartic_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 0, 0, -1, -1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_small_quartic_candidate_fb2x(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 0, 0, -1, -1);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 2, kTrivialClassGroup);
}

void BM_class_group_small_quartic_candidate_fb4x(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 0, 0, -1, -1);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 4, kTrivialClassGroup);
}

void BM_class_group_quartic_1_0_minus2_1_candidate(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 1, 0, -2, 1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_quartic_minus1_minus1_0_2_candidate(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -1, -1, 0, 2);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_quartic_minus1_minus1_0_2_candidate_fb2x(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -1, -1, 0, 2);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 2, kTrivialClassGroup);
}

void BM_class_group_quartic_minus1_minus1_0_2_candidate_fb4x(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -1, -1, 0, 2);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 4, kTrivialClassGroup);
}

void BM_class_group_hard_cubic_0_minus4_minus7_candidate(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, 0, -4, -7);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_hard_quartic_minus1_minus1_0_2_candidate(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -1, -1, 0, 2);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_hard_quartic_minus1_minus1_0_2_fb2x_candidate(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -1, -1, 0, 2);
    BM_class_group_candidate_from_poly_with_bound_scale(
            state, polynomial, 2, kTrivialClassGroup);
}

void BM_class_group_quintic_minus1_minus1_candidate(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quintic_polynomial(polynomial, 0, 0, 0, -1, -1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

void BM_class_group_sextic_minus1_minus1_candidate(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_sextic_polynomial(polynomial, 0, 0, 0, 0, -1, -1);
    BM_class_group_candidate_from_poly(state, polynomial, kTrivialClassGroup);
}

ulong random_sweep_next(ulong& x) noexcept {
    x = (x * 6364136223846793005UL) + 1442695040888963407UL;
    return x;
}

void set_random_sweep_polynomial(sflint::FmpqPoly& polynomial,
                                 slong degree,
                                 slong height,
                                 ulong seed,
                                 slong sample) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, degree, 1);

    ulong x = seed ^ (static_cast<ulong>(degree) * 0x9e3779b97f4a7c15UL) ^
            (static_cast<ulong>(height) * 0xbf58476d1ce4e5b9UL) ^
            (static_cast<ulong>(sample) * 0x94d049bb133111ebUL);

    for (slong i = 0; i < degree; ++i) {
        const slong width = 2 * height + 1;
        slong coefficient = static_cast<slong>(
                random_sweep_next(x) % static_cast<ulong>(width)) -
                height;
        if (i == 0 && coefficient == 0) {
            coefficient = static_cast<slong>(
                    1 + (random_sweep_next(x) % static_cast<ulong>(height)));
        }
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficient);
    }
}

bool random_sweep_polynomial_is_irreducible(
        const sflint::FmpqPoly& polynomial) noexcept {
    sflint::FmpzPoly numerator;
    sflint::FmpzPolyFactor factorization;
    sflint::fmpq_poly_get_numerator(numerator, polynomial);
    sflint::fmpz_poly_factor(sflint::FmpzPolyFactorRef(factorization),
                             sflint::FmpzPolyConstRef(numerator));
    const sflint::FmpzPolyFactorConstRef factors(factorization);
    return sflint::fmpz_poly_factor_num(factors) == 1 &&
            sflint::fmpz_poly_factor_exp(factors, 0) == 1 &&
            sflint::fmpz_poly_factor_poly_degree(factors, 0) ==
                    sflint::fmpz_poly_degree(numerator);
}

bool configure_class_unit_random_sweep_options(
        silex::ClassGroupComputeOptions& options,
        sflint::Fmpz& factor_base_bound,
        const silex::Order& order,
        slong degree,
        silex::CertificationMode requested_certification) noexcept {
    options = silex::ClassGroupComputeOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }

    if (degree >= 3) {
        options.max_candidates = 5000;
        options.max_relations = degree >= 5 ? 1000 : 500;
    }

    options.requested_certification = requested_certification;
    options.zeta_bf_max_cutoff = 20000;
    return true;
}

void set_random_sweep_metadata(
        benchmark::State& state,
        slong degree,
        slong height,
        ulong seed,
        slong sample,
        silex::CertificationMode requested_certification,
        ClassUnitFailureReason reason) {
    state.counters["degree"] = degree;
    state.counters["height"] = height;
    state.counters["seed"] = seed;
    state.counters["sample"] = sample;
    state.counters["requested_cert"] =
            certification_counter(requested_certification);
    state.counters["failure_reason"] = static_cast<double>(reason);
    state.counters["success"] =
            reason == ClassUnitFailureReason::none ? 1.0 : 0.0;
}

void set_polynomial_coefficient_counters(benchmark::State& state,
                                         const sflint::FmpqPoly& polynomial) {
    const slong degree = sflint::fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        sflint::Fmpz coefficient;
        sflint::fmpq_poly_get_coeff_fmpz(
                sflint::FmpzRef(coefficient), polynomial, i);
        const std::string name =
                "poly_c" + std::to_string(static_cast<long>(i));
        state.counters[name] =
                sflint::fmpz_get_d(sflint::FmpzConstRef(coefficient));
    }
}

void consume_failed_random_sweep_row(
        benchmark::State& state,
        ClassUnitFailureReason reason) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(reason);
    }
}

bool define_random_sweep_setup(FieldOrderSetup& setup,
                               benchmark::State& state,
                               const sflint::FmpqPoly& polynomial,
                               slong degree,
                               slong height,
                               ulong seed,
                               slong sample,
                               silex::CertificationMode requested_cert) {
    if (!random_sweep_polynomial_is_irreducible(polynomial)) {
        set_random_sweep_metadata(
                state, degree, height, seed, sample, requested_cert,
                ClassUnitFailureReason::reducible_polynomial);
        consume_failed_random_sweep_row(
                state, ClassUnitFailureReason::reducible_polynomial);
        return false;
    }
    if (!setup.define(polynomial)) {
        set_random_sweep_metadata(
                state, degree, height, seed, sample, requested_cert,
                ClassUnitFailureReason::field_or_order_unavailable);
        consume_failed_random_sweep_row(
                state, ClassUnitFailureReason::field_or_order_unavailable);
        return false;
    }
    return true;
}

void BM_class_unit_random_sweep(benchmark::State& state) {
    const slong degree = state.range(0);
    const slong height = state.range(1);
    const slong sample = state.range(2);
    const auto requested_certification = silex::CertificationMode::proven;
    const ulong seed = 0x5eed0100UL + static_cast<ulong>(state.range(4));

    sflint::FmpqPoly polynomial;
    set_random_sweep_polynomial(polynomial, degree, height, seed, sample);

    FieldOrderSetup setup;
    if (!define_random_sweep_setup(setup, state, polynomial, degree, height,
                                   seed, sample, requested_certification)) {
        return;
    }

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_class_unit_random_sweep_options(
                options, factor_base_bound, setup.maximal_order, degree,
                requested_certification)) {
        set_random_sweep_metadata(
                state, degree, height, seed, sample, requested_certification,
                ClassUnitFailureReason::options_unavailable);
        consume_failed_random_sweep_row(
                state, ClassUnitFailureReason::options_unavailable);
        return;
    }

    bool terminal_failure = false;
    bool ok = false;
    for (auto _ : state) {
        if (terminal_failure) {
            benchmark::DoNotOptimize(ok);
            continue;
        }

        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        ok = units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128);
        benchmark::DoNotOptimize(ok);
        if (!ok) {
            terminal_failure = true;
            state.counters["failure_reason"] = static_cast<double>(
                    ClassUnitFailureReason::compute_unavailable);
            set_random_sweep_metadata(
                    state, degree, height, seed, sample,
                    requested_certification,
                    ClassUnitFailureReason::compute_unavailable);
        } else {
            set_random_sweep_metadata(
                    state, degree, height, seed, sample,
                    requested_certification, ClassUnitFailureReason::none);
        }
        set_class_unit_counters(state, setup.maximal_order, class_group,
                                units, options, ok);
    }

    benchmark::ClobberMemory();
}

void BM_class_unit_random_matrix(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    const slong degree = state.range(0);
    const slong height = state.range(1);
    const slong sample = state.range(2);
    const auto mode = static_cast<ClassUnitMatrixMode>(state.range(3));
    const ulong seed = 0x5eed0100UL + static_cast<ulong>(state.range(4));
    const auto requested_certification =
            mode == ClassUnitMatrixMode::pair_proven
            ? silex::CertificationMode::proven
            : silex::CertificationMode::unknown;
    const bool release_semantic_row =
            degree == 3 && height == 4 && sample == 0 &&
            mode == ClassUnitMatrixMode::pair_proven && state.range(4) == 0;
    if (release_semantic_row) {
        initialize_exact_class_unit_counters(state, kRandomCubicClassUnit);
    }

    sflint::FmpqPoly polynomial;
    set_random_sweep_polynomial(polynomial, degree, height, seed, sample);
    set_polynomial_coefficient_counters(state, polynomial);

    FieldOrderSetup setup;
    if (!define_random_sweep_setup(setup, state, polynomial, degree, height,
                                   seed, sample, requested_certification)) {
        state.counters["mode"] = static_cast<double>(mode);
        if (release_semantic_row) {
            silex::bench_contract::fail(
                    state, "release random-matrix setup failed",
                    silex::bench_contract::FailureReason::setup);
        }
        return;
    }

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_class_unit_random_sweep_options(
                options, factor_base_bound, setup.maximal_order, degree,
                requested_certification)) {
        state.counters["mode"] = static_cast<double>(mode);
        set_random_sweep_metadata(
                state, degree, height, seed, sample, requested_certification,
                ClassUnitFailureReason::options_unavailable);
        consume_failed_random_sweep_row(
                state, ClassUnitFailureReason::options_unavailable);
        if (release_semantic_row) {
            silex::bench_contract::fail(
                    state, "release random-matrix options unavailable",
                    silex::bench_contract::FailureReason::setup);
        }
        return;
    }

    BenchmarkProfileCollector profile_collector;
    BenchmarkProofLogCollector proof_log_collector;
    silex::DiagnosticsContext diagnostics;
    const bool profile_enabled =
            BenchmarkProfileCollector::enabled_from_env();
    const bool proof_diag_enabled = proof_diag_enabled_from_env();
    if (profile_enabled || proof_diag_enabled) {
        if (profile_enabled) {
            profile_collector.configure(diagnostics);
        } else {
            silex::diagnostics_context_init(diagnostics);
        }
        if (proof_diag_enabled) {
            proof_log_collector.configure(diagnostics);
        }
        options.diagnostics = &diagnostics;
    }

    bool terminal = false;
    bool ok = false;
    bool completed = false;
    bool failed = false;
    state.counters["mode"] = static_cast<double>(mode);
    for (auto _ : state) {
        if (terminal) {
            benchmark::DoNotOptimize(ok);
            continue;
        }

        if (mode == ClassUnitMatrixMode::class_candidate) {
            silex::ClassGroupCandidateOptions candidate_options;
            candidate_options.max_candidates = options.max_candidates;
            candidate_options.max_relations = options.max_relations;
            candidate_options.diagnostics = options.diagnostics;
            silex::ClassGroupContext context;
            ok = context.compute_candidate(
                    setup.maximal_order,
                    sflint::FmpzConstRef(factor_base_bound),
                    candidate_options);
            benchmark::DoNotOptimize(ok);
            set_random_sweep_metadata(
                    state, degree, height, seed, sample,
                    requested_certification,
                    ok ? ClassUnitFailureReason::none
                       : ClassUnitFailureReason::compute_unavailable);
            set_class_group_counters(state, context,
                                     sflint::FmpzConstRef(factor_base_bound),
                                     ok);
            completed = ok;
            failed = !ok;
            if (release_semantic_row && failed) {
                silex::bench_contract::fail(
                        state,
                        "release random-matrix class-group compute failed",
                        silex::bench_contract::FailureReason::operation);
            }
        } else {
            silex::ClassGroupContext class_group;
            silex::OrderUnitGroup units;
            ok = units.compute_with_class_group(
                    class_group, setup.maximal_order,
                    sflint::FmpzConstRef(factor_base_bound), options, 128);
            benchmark::DoNotOptimize(ok);
            set_random_sweep_metadata(
                    state, degree, height, seed, sample,
                    requested_certification,
                    ok ? ClassUnitFailureReason::none
                       : ClassUnitFailureReason::compute_unavailable);
            set_class_unit_counters(state, setup.maximal_order, class_group,
                                    units, options, ok);
            if (!ok && mode == ClassUnitMatrixMode::pair_proven) {
                set_failed_unit_proof_replay_counters(
                        state, setup.maximal_order,
                        sflint::FmpzConstRef(factor_base_bound), options, 128);
            }
            completed = ok;
            failed = !ok;
            if (release_semantic_row && failed) {
                silex::bench_contract::fail(
                        state,
                        "release random-matrix class/unit compute failed",
                        silex::bench_contract::FailureReason::operation);
            }
        }
        terminal = true;
    }

    if (completed && !failed) {
        bool validation_ok = false;
        bool invariant_ok = false;
        if (mode == ClassUnitMatrixMode::class_candidate) {
            silex::ClassGroupCandidateOptions validation_options;
            validation_options.max_candidates = options.max_candidates;
            validation_options.max_relations = options.max_relations;
            silex::ClassGroupContext validation_context;
            validation_ok = validation_context.compute_candidate(
                    setup.maximal_order,
                    sflint::FmpzConstRef(factor_base_bound),
                    validation_options);
            invariant_ok = validation_ok &&
                           valid_class_group_result(validation_context);
        } else {
            silex::ClassGroupComputeOptions validation_options = options;
            validation_options.diagnostics = nullptr;
            silex::ClassGroupContext validation_class_group;
            silex::OrderUnitGroup validation_units;
            validation_ok = validation_units.compute_with_class_group(
                    validation_class_group, setup.maximal_order,
                    sflint::FmpzConstRef(factor_base_bound),
                    validation_options, 128);
            invariant_ok = validation_ok &&
                    (release_semantic_row
                             ? validate_exact_class_unit_result(
                                       state, setup.maximal_order,
                                       validation_class_group,
                                       validation_units,
                                       requested_certification,
                                       kRandomCubicClassUnit)
                             : valid_class_unit_result(
                                       validation_class_group,
                                       validation_units,
                                       requested_certification));
        }
        if (!invariant_ok) {
            failed = true;
            if (release_semantic_row) {
                silex::bench_contract::fail(
                        state,
                        validation_ok
                                ? "release random-matrix invariant failed"
                                : "release random-matrix validation rerun failed",
                        validation_ok
                                ? silex::bench_contract::FailureReason::invariant
                                : silex::bench_contract::FailureReason::operation);
            }
        }
    }

    benchmark::ClobberMemory();
    if (proof_diag_enabled) {
        proof_log_collector.set_counters(state);
    }
    if (profile_enabled) {
        profile_collector.set_stable_scan_counters(state);
        profile_collector.set_proof_selector_scan_counters(state);
        profile_collector.set_relation_path_counters(state);
        profile_collector.set_partial_relation_counters(state);
        profile_collector.set_validation_recompute_counters(state);
        profile_collector.set_class_unit_control_flow_counters(state);
        profile_collector.set_saturation_counters(state);
        profile_collector.set_class_unit_stage_time_counters(state);
        profile_collector.print(state.name().c_str());
    }
    if (completed && !failed) {
        silex::bench_contract::succeed(state);
    }
}

void set_random_boundary_polynomial(sflint::FmpqPoly& polynomial) noexcept {
    set_random_sweep_polynomial(polynomial, 3, 4, 0x5eed0100UL, 0);
}

bool configure_class_unit_boundary_options(
        silex::ClassGroupComputeOptions& options,
        sflint::Fmpz& factor_base_bound,
        const silex::Order& order,
        slong max_candidates,
        silex::CertificationMode requested_certification) noexcept {
    options = silex::ClassGroupComputeOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }
    options.max_candidates = max_candidates;
    options.max_relations = 500;
    options.requested_certification = requested_certification;
    options.zeta_bf_max_cutoff = 20000;
    return true;
}

bool define_random_boundary_setup(FieldOrderSetup& setup,
                                  benchmark::State& state,
                                  sflint::FmpqPoly& polynomial) {
    set_random_boundary_polynomial(polynomial);
    set_polynomial_coefficient_counters(state, polynomial);
    return define_random_sweep_setup(setup, state, polynomial, 3, 4,
                                     0x5eed0100UL, 0,
                                     silex::CertificationMode::unknown);
}

void BM_class_unit_random_boundary_clgp(benchmark::State& state) {
    silex::bench_contract::initialize(state);
    const slong max_candidates = state.range(0);
    sflint::FmpqPoly polynomial;
    FieldOrderSetup setup;
    if (!define_random_boundary_setup(setup, state, polynomial)) {
        silex::bench_contract::fail(
                state, "random-boundary field/order setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_class_unit_boundary_options(
                options, factor_base_bound, setup.maximal_order,
                max_candidates, silex::CertificationMode::unknown)) {
        consume_failed_random_sweep_row(
                state, ClassUnitFailureReason::options_unavailable);
        silex::bench_contract::fail(
                state, "random-boundary options unavailable",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    state.counters["max_candidates"] = max_candidates;
    state.counters["mode"] = static_cast<double>(
            ClassUnitMatrixMode::class_candidate);
    silex::ClassGroupCandidateOptions candidate_options;
    candidate_options.max_candidates = options.max_candidates;
    candidate_options.max_relations = options.max_relations;
    candidate_options.diagnostics = options.diagnostics;
    bool completed = false;
    bool failed = false;
    for (auto _ : state) {
        silex::ClassGroupContext context;
        const bool ok = context.compute_candidate(
                setup.maximal_order, sflint::FmpzConstRef(factor_base_bound),
                candidate_options);
        int ok_value = ok ? 1 : 0;
        benchmark::DoNotOptimize(ok_value);
        completed = completed || ok;
        failed = failed || !ok;
        set_class_group_counters(state, context,
                                 sflint::FmpzConstRef(factor_base_bound), ok);
    }

    benchmark::ClobberMemory();
    if (!completed || failed) {
        silex::bench_contract::fail(
                state, "random-boundary class-group compute failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    silex::bench_contract::succeed(state);
}

bool configure_class_unit_0_1_0_options(
        silex::ClassGroupComputeOptions& options,
        sflint::Fmpz& factor_base_bound,
        const silex::Order& order,
        bool is_cubic_nontrivial,
        bool is_quartic,
        bool is_quintic,
        bool is_sextic,
        silex::CertificationMode requested_certification,
        ulong zeta_bf_max_cutoff) noexcept {
    options = silex::ClassGroupComputeOptions{};
    options.max_candidates = WORD_MAX;
    options.max_relations = WORD_MAX;

    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }

    if (is_cubic_nontrivial) {
        options.max_candidates = 5000;
        options.max_relations = 500;
    } else if (is_quartic) {
        options.max_candidates = 5000;
        options.max_relations = 500;
    } else if (is_quintic || is_sextic) {
        options.max_candidates = 5000;
        options.max_relations = 1000;
    }

    options.requested_certification = requested_certification;
    options.zeta_bf_max_cutoff = zeta_bf_max_cutoff;
    return true;
}

void BM_class_unit_0_1_0_from_setup(
        benchmark::State& state,
        FieldOrderSetup& setup,
        bool is_cubic_nontrivial,
        bool is_quartic,
        bool is_quintic,
        bool is_sextic,
        silex::CertificationMode requested_certification,
        ulong zeta_bf_max_cutoff,
        ExactClassUnitExpectation expected) {
    initialize_exact_class_unit_counters(state, expected);
    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_class_unit_0_1_0_options(
                options, factor_base_bound, setup.maximal_order,
                is_cubic_nontrivial, is_quartic, is_quintic, is_sextic,
                requested_certification, zeta_bf_max_cutoff)) {
        silex::bench_contract::fail(
                state, "class/unit benchmark options unavailable",
                silex::bench_contract::FailureReason::setup);
        return;
    }

    BenchmarkProfileCollector profile_collector;
    BenchmarkProofLogCollector proof_log_collector;
    silex::DiagnosticsContext diagnostics;
    const bool profile_enabled =
            BenchmarkProfileCollector::enabled_from_env();
    const bool proof_diag_enabled = proof_diag_enabled_from_env();
    if (profile_enabled || proof_diag_enabled) {
        if (profile_enabled) {
            profile_collector.configure(diagnostics);
        } else {
            silex::diagnostics_context_init(diagnostics);
        }
        if (proof_diag_enabled) {
            proof_log_collector.configure(diagnostics);
        }
        options.diagnostics = &diagnostics;
    }

    bool completed = false;
    bool failed = false;
    for (auto _ : state) {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        const bool ok = units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128);
        int ok_value = ok ? 1 : 0;
        benchmark::DoNotOptimize(ok_value);
        set_class_unit_counters(state, setup.maximal_order, class_group, units,
                                options, ok);
        if (!ok) {
            failed = true;
            silex::bench_contract::fail(
                    state, "class/unit compute failed",
                    silex::bench_contract::FailureReason::operation);
            break;
        }
        completed = true;
    }

    if (completed && !failed) {
        silex::ClassGroupComputeOptions validation_options = options;
        validation_options.diagnostics = nullptr;
        silex::ClassGroupContext validation_class_group;
        silex::OrderUnitGroup validation_units;
        const bool validation_ok =
                validation_units.compute_with_class_group(
                        validation_class_group, setup.maximal_order,
                        sflint::FmpzConstRef(factor_base_bound),
                        validation_options, 128);
        const bool invariant_ok = validation_ok &&
                validate_exact_class_unit_result(
                        state, setup.maximal_order, validation_class_group,
                        validation_units, requested_certification, expected);
        if (!invariant_ok) {
            failed = true;
            silex::bench_contract::fail(
                    state,
                    validation_ok ? "class/unit result invariant failed"
                                  : "class/unit validation rerun failed",
                    validation_ok
                            ? silex::bench_contract::FailureReason::invariant
                            : silex::bench_contract::FailureReason::operation);
        }
    }

    benchmark::ClobberMemory();
    if (proof_diag_enabled) {
        proof_log_collector.set_counters(state);
    }
    if (profile_enabled) {
        profile_collector.set_stable_scan_counters(state);
        profile_collector.set_proof_selector_scan_counters(state);
        profile_collector.set_relation_path_counters(state);
        profile_collector.set_partial_relation_counters(state);
        profile_collector.set_validation_recompute_counters(state);
        profile_collector.set_class_unit_control_flow_counters(state);
        profile_collector.set_saturation_counters(state);
        profile_collector.set_class_unit_stage_time_counters(state);
        profile_collector.print(state.name().c_str());
    }
    if (completed && !failed) {
        silex::bench_contract::succeed(state);
    }
}

void BM_class_unit_0_1_0_from_poly(
        benchmark::State& state,
        const sflint::FmpqPoly& polynomial,
        silex::CertificationMode requested_certification,
        ulong zeta_bf_max_cutoff,
        ExactClassUnitExpectation expected,
        bool is_cubic_nontrivial = false,
        bool is_quartic = false,
        bool is_quintic = false,
        bool is_sextic = false) {
    silex::bench_contract::initialize(state);
    initialize_exact_class_unit_counters(state, expected);
    FieldOrderSetup setup;
    if (!setup.define(polynomial)) {
        silex::bench_contract::fail(
                state, "field/order setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    BM_class_unit_0_1_0_from_setup(
            state, setup, is_cubic_nontrivial, is_quartic, is_quintic,
            is_sextic, requested_certification, zeta_bf_max_cutoff,
            expected);
}

void BM_class_unit_0_1_0_from_quadratic(
        benchmark::State& state,
        slong radicand,
        silex::CertificationMode requested_certification,
        ulong zeta_bf_max_cutoff,
        ExactClassUnitExpectation expected) {
    silex::bench_contract::initialize(state);
    initialize_exact_class_unit_counters(state, expected);
    FieldOrderSetup setup;
    sflint::Fmpz d;
    sflint::fmpz_set_si(sflint::FmpzRef(d), radicand);
    if (!setup.define_quadratic(sflint::FmpzConstRef(d))) {
        silex::bench_contract::fail(
                state, "quadratic field/order setup failed",
                silex::bench_contract::FailureReason::setup);
        return;
    }
    BM_class_unit_0_1_0_from_setup(
            state, setup, false, false, false, false, requested_certification,
            zeta_bf_max_cutoff, expected);
}

void BM_class_unit_0_1_0_degree_one_proven(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_degree_one_polynomial(polynomial);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kDegreeOneClassUnit);
}

void BM_class_unit_0_1_0_real_quadratic_proven(
        benchmark::State& state) {
    BM_class_unit_0_1_0_from_quadratic(
            state, 5, silex::CertificationMode::proven, 20000,
            kRealQuadraticClassUnit);
}

void BM_class_unit_0_1_0_imag_quadratic_proven(
        benchmark::State& state) {
    BM_class_unit_0_1_0_from_quadratic(
            state, -47, silex::CertificationMode::proven, 20000,
            kImaginaryQuadraticClassUnit);
}

void BM_class_unit_0_1_0_cubic_trivial_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, 0, -2);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kCubicTrivialClassUnit);
}

void BM_class_unit_0_1_0_cubic_nontrivial_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_cubic_polynomial(polynomial, -2, -5);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kCubicNontrivialClassUnit, true);
}

void BM_class_unit_0_1_0_quartic_cyclotomic_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 0, 0, 0, 1);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticCyclotomicClassUnit, false, true);
}

void BM_class_unit_0_1_0_quartic_noncyclotomic_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -1, -1, 0, 2);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticNoncyclotomicClassUnit, false, true);
}

void BM_class_unit_0_1_0_quartic_disc70640_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 4, 1, 2, -3);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticDisc70640ClassUnit, false, true);
}

void BM_class_unit_0_1_0_quartic_disc223479_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -4, -3, -3, 3);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticDisc223479ClassUnit, false, true);
}

void BM_class_unit_0_1_0_quartic_disc35019_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 0, -2, 3, -5);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticDisc35019ClassUnit, false, true);
}

void BM_class_unit_0_1_0_quartic_disc1412343_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, -8, 3, 1, 4);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticDisc1412343ClassUnit, false, true);
}

void BM_class_unit_0_1_0_quartic_x4_minus_x_minus_1_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quartic_polynomial(polynomial, 0, 0, -1, -1);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuarticSmallClassUnit, false, true);
}

void BM_class_unit_0_1_0_quintic_proven(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quintic_polynomial(polynomial, 0, 0, 0, -1, -1);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuinticSmallClassUnit, false, false, true);
}

void BM_class_unit_0_1_0_quintic_disc11119_proven(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quintic_polynomial(polynomial, -2, 1, 1, -3, 1);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuinticDisc11119ClassUnit, false, false, true);
}

void BM_class_unit_0_1_0_quintic_disc401370255_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quintic_polynomial(polynomial, -7, -6, -3, -3, 3);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuinticDisc401370255ClassUnit, false, false, true);
}

void BM_class_unit_0_1_0_quintic_disc57895_proven(
        benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_quintic_polynomial(polynomial, 2, -1, -2, -2, 1);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kQuinticDisc57895ClassUnit, false, false, true);
}

void BM_diagnostic_class_unit_0_1_0_quintic_disc57895_cap(
        benchmark::State& state,
        silex::CertificationMode requested_certification) {
    const slong max_candidates = state.range(0);
    const slong max_relations = state.range(1);
    sflint::FmpqPoly polynomial;
    set_quintic_polynomial(polynomial, 2, -1, -2, -2, 1);

    FieldOrderSetup setup;
    if (!setup.define(polynomial)) {
        state.counters["failure_reason"] = static_cast<double>(
                ClassUnitFailureReason::field_or_order_unavailable);
        state.counters["success"] = 0.0;
        return;
    }

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions options;
    if (!configure_class_unit_0_1_0_options(
                options, factor_base_bound, setup.maximal_order, false, false,
                true, false, requested_certification, 20000)) {
        state.counters["failure_reason"] = static_cast<double>(
                ClassUnitFailureReason::options_unavailable);
        state.counters["success"] = 0.0;
        return;
    }
    options.max_candidates = max_candidates;
    options.max_relations = max_relations;

    BenchmarkProfileCollector profile_collector;
    silex::DiagnosticsContext diagnostics;
    const bool profile_enabled =
            BenchmarkProfileCollector::enabled_from_env();
    if (profile_enabled) {
        profile_collector.configure(diagnostics);
        options.diagnostics = &diagnostics;
    }

    state.counters["max_candidates"] = max_candidates;
    state.counters["max_relations"] = max_relations;
    state.counters["mode"] = 5.0;
    for (auto _ : state) {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        const bool ok = units.compute_with_class_group(
                class_group, setup.maximal_order,
                sflint::FmpzConstRef(factor_base_bound), options, 128);
        int ok_value = ok ? 1 : 0;
        benchmark::DoNotOptimize(ok_value);
        if (!ok) {
            state.counters["failure_reason"] = static_cast<double>(
                    ClassUnitFailureReason::compute_unavailable);
        }
        set_partial_class_group_counters(
                state, class_group, sflint::FmpzConstRef(factor_base_bound));
        set_class_unit_counters(state, setup.maximal_order, class_group,
                                units, options, ok);
    }

    benchmark::ClobberMemory();
    if (profile_enabled) {
        profile_collector.set_stable_scan_counters(state);
        profile_collector.set_proof_selector_scan_counters(state);
        profile_collector.set_relation_path_counters(state);
        profile_collector.set_partial_relation_counters(state);
        profile_collector.set_validation_recompute_counters(state);
        profile_collector.set_class_unit_control_flow_counters(state);
        profile_collector.set_saturation_counters(state);
        profile_collector.set_class_unit_stage_time_counters(state);
        profile_collector.print(state.name().c_str());
    }
}

void BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap(
        benchmark::State& state) {
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_cap(
            state, silex::CertificationMode::proven);
}

void BM_class_unit_0_1_0_sextic_proven(benchmark::State& state) {
    sflint::FmpqPoly polynomial;
    set_sextic_polynomial(polynomial, 0, 0, 0, 0, -1, -1);
    BM_class_unit_0_1_0_from_poly(
            state, polynomial, silex::CertificationMode::proven, 20000,
            kSexticSmallClassUnit, false, false, false, true);
}

}  // namespace

BENCHMARK(BM_class_group_trivial_cubic_candidate);
BENCHMARK(BM_class_group_nontrivial_cubic_candidate);
BENCHMARK(BM_class_group_nontrivial_cubic_candidate_fb2x);
BENCHMARK(BM_class_group_nontrivial_cubic_candidate_fb4x);
BENCHMARK(BM_class_group_mixed_cubic_candidate);
BENCHMARK(BM_class_group_cubic_3_1_candidate);
BENCHMARK(BM_class_group_cubic_minus3_minus1_candidate);
BENCHMARK(BM_class_group_small_quartic_candidate);
BENCHMARK(BM_class_group_small_quartic_candidate_fb2x);
BENCHMARK(BM_class_group_small_quartic_candidate_fb4x);
BENCHMARK(BM_class_group_quartic_1_0_minus2_1_candidate);
BENCHMARK(BM_class_group_quartic_minus1_minus1_0_2_candidate);
BENCHMARK(BM_class_group_quartic_minus1_minus1_0_2_candidate_fb2x);
BENCHMARK(BM_class_group_quartic_minus1_minus1_0_2_candidate_fb4x);
BENCHMARK(BM_class_group_hard_cubic_0_minus4_minus7_candidate);
BENCHMARK(BM_class_group_hard_quartic_minus1_minus1_0_2_candidate);
BENCHMARK(BM_class_group_hard_quartic_minus1_minus1_0_2_fb2x_candidate);
BENCHMARK(BM_class_group_quintic_minus1_minus1_candidate);
BENCHMARK(BM_class_group_sextic_minus1_minus1_candidate);

BENCHMARK(BM_class_unit_0_1_0_degree_one_proven);
BENCHMARK(BM_class_unit_0_1_0_real_quadratic_proven);
BENCHMARK(BM_class_unit_0_1_0_imag_quadratic_proven);
BENCHMARK(BM_class_unit_0_1_0_cubic_trivial_proven);
BENCHMARK(BM_class_unit_0_1_0_cubic_nontrivial_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_cyclotomic_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_noncyclotomic_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_disc70640_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_disc223479_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_disc35019_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_disc1412343_proven);
BENCHMARK(BM_class_unit_0_1_0_quartic_x4_minus_x_minus_1_proven);
BENCHMARK(BM_class_unit_0_1_0_quintic_proven);
BENCHMARK(BM_class_unit_0_1_0_quintic_disc11119_proven);
BENCHMARK(BM_class_unit_0_1_0_quintic_disc401370255_proven);
BENCHMARK(BM_class_unit_0_1_0_quintic_disc57895_proven);
BENCHMARK(BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap)
        ->Args({100, 100})
        ->Args({250, 250})
        ->Args({500, 500})
        ->Args({750, 750})
        ->Args({1000, 1000})
        ->Args({1500, 1500});
BENCHMARK(BM_class_unit_0_1_0_sextic_proven);
BENCHMARK(BM_class_unit_random_sweep)
        ->Args({2, 4, 0, 1, 0})
        ->Args({2, 8, 0, 1, 0})
        ->Args({3, 2, 0, 1, 0});
BENCHMARK(BM_class_unit_random_matrix)
        ->Args({2, 2, 0, 0, 0})
        ->Args({2, 4, 0, 0, 0})
        ->Args({2, 8, 0, 0, 0})
        ->Args({3, 2, 0, 0, 0})
        ->Args({3, 4, 0, 0, 0})
        ->Args({3, 8, 0, 0, 0})
        ->Args({4, 2, 0, 0, 0})
        ->Args({4, 4, 0, 0, 0})
        ->Args({5, 2, 0, 0, 0})
        ->Args({6, 2, 2, 0, 0})
        ->Args({2, 4, 0, 1, 0})
        ->Args({2, 8, 0, 1, 0})
        ->Args({3, 2, 0, 1, 0})
        ->Args({3, 4, 0, 1, 0})
        ->Args({4, 2, 0, 1, 0})
        ->Args({5, 2, 0, 1, 0})
        ->Iterations(1);
BENCHMARK(BM_class_unit_random_boundary_clgp)
        ->Args({250})
        ->Args({1000})
        ->Args({5000});
BENCHMARK_MAIN();
