#include <silex/diagnostics.hpp>

namespace silex {

void diagnostics_context_init(DiagnosticsContext& ctx) noexcept {
    ctx = DiagnosticsContext{};
}

const char* diagnostics_module_name(DiagnosticsModule module) noexcept {
    switch (module) {
    case DiagnosticsModule::core: return "core";
    case DiagnosticsModule::flint: return "flint";
    case DiagnosticsModule::field: return "field";
    case DiagnosticsModule::element: return "element";
    case DiagnosticsModule::order: return "order";
    case DiagnosticsModule::ideal: return "ideal";
    case DiagnosticsModule::fractional_ideal: return "fractional_ideal";
    case DiagnosticsModule::prime_ideal: return "prime_ideal";
    case DiagnosticsModule::lattice: return "lattice";
    case DiagnosticsModule::residue: return "residue";
    case DiagnosticsModule::relation: return "relation";
    case DiagnosticsModule::class_group: return "class_group";
    case DiagnosticsModule::unit_group: return "unit_group";
    case DiagnosticsModule::backend: return "backend";
    }
    return "unknown";
}

const char* log_level_name(LogLevel level) noexcept {
    switch (level) {
    case LogLevel::off: return "off";
    case LogLevel::error: return "error";
    case LogLevel::warn: return "warn";
    case LogLevel::info: return "info";
    case LogLevel::detail: return "detail";
    case LogLevel::trace: return "trace";
    }
    return "unknown";
}

const char* verbose_level_name(VerboseLevel level) noexcept {
    switch (level) {
    case VerboseLevel::off: return "off";
    case VerboseLevel::progress: return "progress";
    case VerboseLevel::detail: return "detail";
    case VerboseLevel::trace: return "trace";
    }
    return "unknown";
}

const char* debug_level_name(DebugLevel level) noexcept {
    switch (level) {
    case DebugLevel::off: return "off";
    case DebugLevel::cheap: return "cheap";
    case DebugLevel::normal: return "normal";
    case DebugLevel::expensive: return "expensive";
    case DebugLevel::exhaustive: return "exhaustive";
    }
    return "unknown";
}

void diagnostics_set_verbose(DiagnosticsContext& ctx,
                             VerboseLevel level,
                             DiagnosticsModuleMask modules,
                             VerboseCallback callback,
                             void* user) noexcept {
    ctx.verbose_level = level;
    ctx.verbose_modules = modules;
    ctx.verbose_callback = callback;
    ctx.verbose_user = user;
}

void diagnostics_set_logging(DiagnosticsContext& ctx,
                             LogLevel level,
                             DiagnosticsModuleMask modules,
                             LogCallback callback,
                             void* user) noexcept {
    ctx.log_level = level;
    ctx.log_modules = modules;
    ctx.log_callback = callback;
    ctx.log_user = user;
}

void diagnostics_set_debug_checks(DiagnosticsContext& ctx,
                                  DebugLevel level,
                                  DiagnosticsModuleMask modules,
                                  DebugFailureCallback callback,
                                  void* user) noexcept {
    ctx.debug_level = level;
    ctx.debug_modules = modules;
    ctx.debug_failure_callback = callback;
    ctx.debug_failure_user = user;
}

void diagnostics_set_profiling(DiagnosticsContext& ctx,
                               bool enabled,
                               DiagnosticsModuleMask modules,
                               ProfileCallback callback,
                               void* user) noexcept {
    ctx.profiling_enabled = enabled;
    ctx.profiling_modules = modules;
    ctx.profile_callback = callback;
    ctx.profile_user = user;
}

bool verbose_enabled(const DiagnosticsContext* ctx,
                     DiagnosticsModule module,
                     VerboseLevel level) noexcept {
    if (ctx == nullptr || level == VerboseLevel::off ||
        ctx->verbose_level == VerboseLevel::off) {
        return false;
    }
    if ((ctx->verbose_modules & diagnostics_module_bit(module)) == 0) {
        return false;
    }
    return static_cast<std::int32_t>(level) <=
           static_cast<std::int32_t>(ctx->verbose_level);
}

void verbose_emit(const DiagnosticsContext* ctx,
                  DiagnosticsModule module,
                  VerboseLevel level,
                  const char* function,
                  const char* message,
                  const char* detail) noexcept {
    if (!verbose_enabled(ctx, module, level) || ctx->verbose_callback == nullptr) {
        return;
    }
    ctx->verbose_callback(ctx->verbose_user, module, level, function, message, detail);
}

bool log_enabled(const DiagnosticsContext* ctx,
                 DiagnosticsModule module,
                 LogLevel level) noexcept {
    if (ctx == nullptr || level == LogLevel::off || ctx->log_level == LogLevel::off) {
        return false;
    }
    if ((ctx->log_modules & diagnostics_module_bit(module)) == 0) {
        return false;
    }
    return static_cast<std::int32_t>(level) <= static_cast<std::int32_t>(ctx->log_level);
}

void log_emit(const DiagnosticsContext* ctx,
              DiagnosticsModule module,
              LogLevel level,
              const char* function,
              const char* message,
              const char* detail) noexcept {
    if (!log_enabled(ctx, module, level) || ctx->log_callback == nullptr) {
        return;
    }
    ctx->log_callback(ctx->log_user, module, level, function, message, detail);
}

bool debug_check_enabled(const DiagnosticsContext* ctx,
                         DiagnosticsModule module,
                         DebugLevel level) noexcept {
    if (ctx == nullptr || level == DebugLevel::off || ctx->debug_level == DebugLevel::off) {
        return false;
    }
    if ((ctx->debug_modules & diagnostics_module_bit(module)) == 0) {
        return false;
    }
    return static_cast<std::int32_t>(level) <= static_cast<std::int32_t>(ctx->debug_level);
}

void debug_check_about_to_run(const DiagnosticsContext* ctx,
                              DiagnosticsModule module,
                              DebugLevel level,
                              const char* function,
                              const char* label) noexcept {
    if (static_cast<std::int32_t>(level) < static_cast<std::int32_t>(DebugLevel::expensive)) {
        return;
    }
    log_emit(ctx, module, LogLevel::detail, function, "running expensive debug check", label);
}

void debug_check_failed(const DiagnosticsContext* ctx,
                        DiagnosticsModule module,
                        DebugLevel level,
                        const char* function,
                        const char* file,
                        int line,
                        const char* label,
                        const char* expression) noexcept {
    if (ctx == nullptr || ctx->debug_failure_callback == nullptr) {
        return;
    }
    ctx->debug_failure_callback(ctx->debug_failure_user, module, level, function, file, line,
                                label, expression);
}

bool profiling_enabled(const DiagnosticsContext* ctx,
                       DiagnosticsModule module) noexcept {
    if (ctx == nullptr || !ctx->profiling_enabled) {
        return false;
    }
    return (ctx->profiling_modules & diagnostics_module_bit(module)) != 0;
}

void profile_emit(const DiagnosticsContext* ctx,
                  DiagnosticsModule module,
                  ProfileEvent event,
                  const char* function,
                  const char* label) noexcept {
    if (!profiling_enabled(ctx, module) || ctx->profile_callback == nullptr) {
        return;
    }
    ctx->profile_callback(ctx->profile_user, module, event, function, label);
}

ProfileScope::ProfileScope(const DiagnosticsContext* ctx,
                           DiagnosticsModule module,
                           const char* function,
                           const char* label) noexcept
    : ctx_(ctx),
      module_(module),
      function_(function),
      label_(label),
      active_(profiling_enabled(ctx, module)) {
    if (active_) {
        profile_emit(ctx_, module_, ProfileEvent::begin_scope, function_, label_);
    }
}

ProfileScope::~ProfileScope() noexcept {
    if (active_) {
        profile_emit(ctx_, module_, ProfileEvent::end_scope, function_, label_);
    }
}

}  // namespace silex
