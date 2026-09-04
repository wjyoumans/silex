#include <silex/diagnostics.hpp>

#include <cassert>
#include <cstring>

namespace {

struct Counters {
    int verbose = 0;
    int logs = 0;
    int debug_failures = 0;
    int profiles = 0;
    const char* last_verbose_message = nullptr;
    const char* last_verbose_detail = nullptr;
    const char* last_log_message = nullptr;
    const char* last_log_detail = nullptr;
};

void verbose_callback(void* user,
                      silex::DiagnosticsModule module,
                      silex::VerboseLevel level,
                      const char*,
                      const char* message,
                      const char* detail) noexcept {
    auto* counters = static_cast<Counters*>(user);
    assert(module == silex::DiagnosticsModule::core);
    assert(level == silex::VerboseLevel::progress);
    (void)module;
    (void)level;
    ++counters->verbose;
    counters->last_verbose_message = message;
    counters->last_verbose_detail = detail;
}

void log_callback(void* user,
                  silex::DiagnosticsModule module,
                  silex::LogLevel level,
                  const char*,
                  const char* message,
                  const char* detail) noexcept {
    auto* counters = static_cast<Counters*>(user);
    assert(module == silex::DiagnosticsModule::core);
    assert(level == silex::LogLevel::info);
    (void)module;
    (void)level;
    ++counters->logs;
    counters->last_log_message = message;
    counters->last_log_detail = detail;
}

void debug_failure_callback(void* user,
                            silex::DiagnosticsModule,
                            silex::DebugLevel,
                            const char*,
                            const char*,
                            int,
                            const char*,
                            const char*) noexcept {
    auto* counters = static_cast<Counters*>(user);
    ++counters->debug_failures;
}

void profile_callback(void* user,
                      silex::DiagnosticsModule,
                      silex::ProfileEvent,
                      const char*,
                      const char*) noexcept {
    auto* counters = static_cast<Counters*>(user);
    ++counters->profiles;
}

}  // namespace

int main() {
    silex::DiagnosticsContext ctx;
    silex::diagnostics_context_init(ctx);

    assert(!silex::verbose_enabled(&ctx, silex::DiagnosticsModule::core,
                                   silex::VerboseLevel::progress));
    assert(!silex::log_enabled(&ctx, silex::DiagnosticsModule::core, silex::LogLevel::info));
    assert(!silex::debug_check_enabled(&ctx, silex::DiagnosticsModule::core,
                                       silex::DebugLevel::cheap));
    assert(!silex::profiling_enabled(&ctx, silex::DiagnosticsModule::core));
    assert(std::strcmp(silex::diagnostics_module_name(silex::DiagnosticsModule::core),
                       "core") == 0);

    Counters counters;
    const auto core_mask = silex::diagnostics_module_bit(silex::DiagnosticsModule::core);
    silex::diagnostics_set_verbose(ctx, silex::VerboseLevel::progress, core_mask,
                                   verbose_callback, &counters);
    silex::diagnostics_set_logging(ctx, silex::LogLevel::info,
                                   core_mask, log_callback, &counters);
    silex::diagnostics_set_debug_checks(ctx, silex::DebugLevel::off,
                                        core_mask, debug_failure_callback, &counters);
    silex::diagnostics_set_profiling(ctx, false,
                                     core_mask, profile_callback, &counters);

    int debug_expression_evaluations = 0;
    SILEX_VERBOSE(&ctx, silex::DiagnosticsModule::core, silex::VerboseLevel::progress,
                  "progress");
    SILEX_LOG(&ctx, silex::DiagnosticsModule::core, silex::LogLevel::info, "status");
    SILEX_DEBUG_CHECK(&ctx, silex::DiagnosticsModule::core, silex::DebugLevel::cheap,
                      "disabled-debug",
                      (++debug_expression_evaluations, true));
    SILEX_PROFILE_EVENT(&ctx, silex::DiagnosticsModule::core, "disabled-profile");

    assert(counters.verbose == 1);
    assert(std::strcmp(counters.last_verbose_message, "progress") == 0);
    assert(counters.last_verbose_detail == nullptr);
    assert(counters.logs == 1);
    assert(std::strcmp(counters.last_log_message, "status") == 0);
    assert(counters.last_log_detail == nullptr);
    assert(counters.debug_failures == 0);
    assert(counters.profiles == 0);
    assert(debug_expression_evaluations == 0);
    (void)debug_expression_evaluations;

    return 0;
}
