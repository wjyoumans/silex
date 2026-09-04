#include <silex/diagnostics.hpp>

#include <cassert>
#include <cstring>

namespace {

struct Counters {
    int verbose = 0;
    int logs = 0;
    int debug_failures = 0;
    int profile_begin = 0;
    int profile_end = 0;
    int profile_events = 0;
    const char* last_verbose_message = nullptr;
    const char* last_verbose_detail = nullptr;
    const char* last_log_message = nullptr;
    const char* last_log_detail = nullptr;
    const char* last_failure_label = nullptr;
    const char* last_failure_expression = nullptr;
};

void verbose_callback(void* user,
                      silex::DiagnosticsModule module,
                      silex::VerboseLevel level,
                      const char*,
                      const char* message,
                      const char* detail) noexcept {
    auto* counters = static_cast<Counters*>(user);
    assert(module == silex::DiagnosticsModule::core);
    assert(level == silex::VerboseLevel::detail);
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
    assert(level == silex::LogLevel::detail);
    (void)module;
    (void)level;
    ++counters->logs;
    counters->last_log_message = message;
    counters->last_log_detail = detail;
}

void debug_failure_callback(void* user,
                            silex::DiagnosticsModule module,
                            silex::DebugLevel level,
                            const char*,
                            const char*,
                            int,
                            const char* label,
                            const char* expression) noexcept {
    auto* counters = static_cast<Counters*>(user);
    assert(module == silex::DiagnosticsModule::core);
    assert(level == silex::DebugLevel::expensive);
    (void)module;
    (void)level;
    ++counters->debug_failures;
    counters->last_failure_label = label;
    counters->last_failure_expression = expression;
}

void profile_callback(void* user,
                      silex::DiagnosticsModule module,
                      silex::ProfileEvent event,
                      const char*,
                      const char* label) noexcept {
    auto* counters = static_cast<Counters*>(user);
    assert(module == silex::DiagnosticsModule::core);
    assert(std::strcmp(label, "scope") == 0 || std::strcmp(label, "event") == 0);
    (void)module;
    (void)label;
    if (event == silex::ProfileEvent::begin_scope) {
        ++counters->profile_begin;
    } else if (event == silex::ProfileEvent::end_scope) {
        ++counters->profile_end;
    } else {
        ++counters->profile_events;
    }
}

}  // namespace

int main() {
    silex::DiagnosticsContext ctx;
    silex::diagnostics_context_init(ctx);

    Counters counters;
    const auto core_mask = silex::diagnostics_module_bit(silex::DiagnosticsModule::core);

    silex::diagnostics_set_verbose(ctx, silex::VerboseLevel::detail, core_mask,
                                   verbose_callback, &counters);
    silex::diagnostics_set_logging(ctx, silex::LogLevel::detail, core_mask,
                                   log_callback, &counters);
    silex::diagnostics_set_debug_checks(ctx, silex::DebugLevel::expensive, core_mask,
                                        debug_failure_callback, &counters);

    int debug_expression_evaluations = 0;
    SILEX_DEBUG_CHECK(&ctx, silex::DiagnosticsModule::core, silex::DebugLevel::expensive,
                      "expensive-ok",
                      (++debug_expression_evaluations, true));
    assert(debug_expression_evaluations == 1);
    assert(counters.verbose == 0);
    assert(counters.logs == 1);
    assert(std::strcmp(counters.last_log_message, "running expensive debug check") == 0);
    assert(std::strcmp(counters.last_log_detail, "expensive-ok") == 0);
    assert(counters.debug_failures == 0);

    silex::diagnostics_set_logging(ctx, silex::LogLevel::off, core_mask,
                                   log_callback, &counters);
    SILEX_DEBUG_CHECK(&ctx, silex::DiagnosticsModule::core, silex::DebugLevel::expensive,
                      "expensive-fail",
                      (++debug_expression_evaluations, false));
    assert(debug_expression_evaluations == 2);
    assert(counters.logs == 1);
    assert(counters.debug_failures == 1);
    assert(std::strcmp(counters.last_failure_label, "expensive-fail") == 0);
    assert(std::strstr(counters.last_failure_expression, "debug_expression_evaluations") !=
           nullptr);

    silex::diagnostics_set_debug_checks(ctx, silex::DebugLevel::off, core_mask,
                                        debug_failure_callback, &counters);
    SILEX_VERBOSE(&ctx, silex::DiagnosticsModule::core, silex::VerboseLevel::detail,
                  "verbose-detail");
    assert(counters.verbose == 1);
    assert(std::strcmp(counters.last_verbose_message, "verbose-detail") == 0);
    assert(counters.last_verbose_detail == nullptr);
    assert(counters.logs == 1);
    assert(counters.debug_failures == 1);

    SILEX_DEBUG_CHECK(&ctx, silex::DiagnosticsModule::core, silex::DebugLevel::cheap,
                      "debug-off",
                      (++debug_expression_evaluations, true));
    assert(debug_expression_evaluations == 2);

    silex::diagnostics_set_profiling(ctx, true, core_mask, profile_callback, &counters);
    {
        SILEX_PROFILE_SCOPE(&ctx, silex::DiagnosticsModule::core, "scope");
        SILEX_PROFILE_EVENT(&ctx, silex::DiagnosticsModule::core, "event");
    }
    assert(counters.profile_begin == 1);
    assert(counters.profile_events == 1);
    assert(counters.profile_end == 1);
    assert(counters.verbose == 1);
    assert(counters.logs == 1);
    assert(counters.debug_failures == 1);
    (void)debug_expression_evaluations;

    return 0;
}
