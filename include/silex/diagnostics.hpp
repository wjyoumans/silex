#pragma once

#include <cstdint>

namespace silex {

enum class DiagnosticsModule : std::uint32_t {
    core = 0,
    flint = 1,
    field = 2,
    element = 3,
    order = 4,
    ideal = 5,
    fractional_ideal = 6,
    prime_ideal = 7,
    lattice = 8,
    residue = 9,
    relation = 10,
    class_group = 11,
    unit_group = 12,
    backend = 13,
};

using DiagnosticsModuleMask = std::uint64_t;

constexpr DiagnosticsModuleMask diagnostics_module_bit(DiagnosticsModule module) noexcept {
    return DiagnosticsModuleMask{1} << static_cast<std::uint32_t>(module);
}

constexpr DiagnosticsModuleMask diagnostics_all_modules = ~DiagnosticsModuleMask{0};

enum class LogLevel : std::int32_t {
    off = 0,
    error = 1,
    warn = 2,
    info = 3,
    detail = 4,
    trace = 5,
};

enum class VerboseLevel : std::int32_t {
    off = 0,
    progress = 1,
    detail = 2,
    trace = 3,
};

enum class DebugLevel : std::int32_t {
    off = 0,
    cheap = 1,
    normal = 2,
    expensive = 3,
    exhaustive = 4,
};

enum class ProfileEvent : std::int32_t {
    begin_scope = 0,
    end_scope = 1,
    event = 2,
};

using LogCallback = void (*)(void* user,
                             DiagnosticsModule module,
                             LogLevel level,
                             const char* function,
                             const char* message,
                             const char* detail) noexcept;

using VerboseCallback = void (*)(void* user,
                                 DiagnosticsModule module,
                                 VerboseLevel level,
                                 const char* function,
                                 const char* message,
                                 const char* detail) noexcept;

using DebugFailureCallback = void (*)(void* user,
                                      DiagnosticsModule module,
                                      DebugLevel level,
                                      const char* function,
                                      const char* file,
                                      int line,
                                      const char* label,
                                      const char* expression) noexcept;

using ProfileCallback = void (*)(void* user,
                                 DiagnosticsModule module,
                                 ProfileEvent event,
                                 const char* function,
                                 const char* label) noexcept;

struct DiagnosticsContext {
    VerboseLevel verbose_level = VerboseLevel::off;
    DiagnosticsModuleMask verbose_modules = 0;
    VerboseCallback verbose_callback = nullptr;
    void* verbose_user = nullptr;

    LogLevel log_level = LogLevel::off;
    DiagnosticsModuleMask log_modules = 0;
    LogCallback log_callback = nullptr;
    void* log_user = nullptr;

    DebugLevel debug_level = DebugLevel::off;
    DiagnosticsModuleMask debug_modules = 0;
    DebugFailureCallback debug_failure_callback = nullptr;
    void* debug_failure_user = nullptr;

    bool profiling_enabled = false;
    DiagnosticsModuleMask profiling_modules = 0;
    ProfileCallback profile_callback = nullptr;
    void* profile_user = nullptr;
};

void diagnostics_context_init(DiagnosticsContext& ctx) noexcept;

const char* diagnostics_module_name(DiagnosticsModule module) noexcept;
const char* log_level_name(LogLevel level) noexcept;
const char* verbose_level_name(VerboseLevel level) noexcept;
const char* debug_level_name(DebugLevel level) noexcept;

void diagnostics_set_verbose(DiagnosticsContext& ctx,
                             VerboseLevel level,
                             DiagnosticsModuleMask modules,
                             VerboseCallback callback,
                             void* user) noexcept;

void diagnostics_set_logging(DiagnosticsContext& ctx,
                             LogLevel level,
                             DiagnosticsModuleMask modules,
                             LogCallback callback,
                             void* user) noexcept;

void diagnostics_set_debug_checks(DiagnosticsContext& ctx,
                                  DebugLevel level,
                                  DiagnosticsModuleMask modules,
                                  DebugFailureCallback callback,
                                  void* user) noexcept;

void diagnostics_set_profiling(DiagnosticsContext& ctx,
                               bool enabled,
                               DiagnosticsModuleMask modules,
                               ProfileCallback callback,
                               void* user) noexcept;

bool log_enabled(const DiagnosticsContext* ctx,
                 DiagnosticsModule module,
                 LogLevel level) noexcept;

bool verbose_enabled(const DiagnosticsContext* ctx,
                     DiagnosticsModule module,
                     VerboseLevel level) noexcept;

void verbose_emit(const DiagnosticsContext* ctx,
                  DiagnosticsModule module,
                  VerboseLevel level,
                  const char* function,
                  const char* message,
                  const char* detail = nullptr) noexcept;

void log_emit(const DiagnosticsContext* ctx,
              DiagnosticsModule module,
              LogLevel level,
              const char* function,
              const char* message,
              const char* detail = nullptr) noexcept;

bool debug_check_enabled(const DiagnosticsContext* ctx,
                         DiagnosticsModule module,
                         DebugLevel level) noexcept;

void debug_check_about_to_run(const DiagnosticsContext* ctx,
                              DiagnosticsModule module,
                              DebugLevel level,
                              const char* function,
                              const char* label) noexcept;

void debug_check_failed(const DiagnosticsContext* ctx,
                        DiagnosticsModule module,
                        DebugLevel level,
                        const char* function,
                        const char* file,
                        int line,
                        const char* label,
                        const char* expression) noexcept;

bool profiling_enabled(const DiagnosticsContext* ctx,
                       DiagnosticsModule module) noexcept;

void profile_emit(const DiagnosticsContext* ctx,
                  DiagnosticsModule module,
                  ProfileEvent event,
                  const char* function,
                  const char* label) noexcept;

class ProfileScope {
public:
    ProfileScope(const DiagnosticsContext* ctx,
                 DiagnosticsModule module,
                 const char* function,
                 const char* label) noexcept;
    ~ProfileScope() noexcept;

    ProfileScope(const ProfileScope&) = delete;
    ProfileScope& operator=(const ProfileScope&) = delete;
    ProfileScope(ProfileScope&&) = delete;
    ProfileScope& operator=(ProfileScope&&) = delete;

private:
    const DiagnosticsContext* ctx_;
    DiagnosticsModule module_;
    const char* function_;
    const char* label_;
    bool active_;
};

}  // namespace silex

#define SILEX_DETAIL_CONCAT_IMPL(a, b) a##b
#define SILEX_DETAIL_CONCAT(a, b) SILEX_DETAIL_CONCAT_IMPL(a, b)

#define SILEX_VERBOSE(ctx, module, level, message)                              \
    do {                                                                        \
        const ::silex::DiagnosticsContext* const silex_verbose_ctx = (ctx);      \
        if (::silex::verbose_enabled(silex_verbose_ctx, (module), (level))) {    \
            ::silex::verbose_emit(silex_verbose_ctx, (module), (level),          \
                                  __func__, (message));                         \
        }                                                                       \
    } while (false)

#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
#define SILEX_LOG(ctx, module, level, message)                             \
    do {                                                                   \
        const ::silex::DiagnosticsContext* const silex_log_ctx = (ctx);     \
        if (::silex::log_enabled(silex_log_ctx, (module), (level))) {       \
            ::silex::log_emit(silex_log_ctx, (module), (level), __func__,   \
                              (message));                                  \
        }                                                                  \
    } while (false)
#else
#define SILEX_LOG(ctx, module, level, message)                  \
    do {                                                        \
        (void)sizeof(ctx);                                      \
    } while (false)
#endif

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
#define SILEX_DEBUG_CHECK(ctx, module, level, label, expression)                  \
    do {                                                                         \
        const ::silex::DiagnosticsContext* const silex_debug_ctx = (ctx);         \
        if (::silex::debug_check_enabled(silex_debug_ctx, (module), (level))) {   \
            ::silex::debug_check_about_to_run(silex_debug_ctx, (module),          \
                                              (level), __func__, (label));        \
            if (!(expression)) {                                                  \
                ::silex::debug_check_failed(silex_debug_ctx, (module), (level),   \
                                            __func__, __FILE__, __LINE__,         \
                                            (label), #expression);                \
            }                                                                    \
        }                                                                        \
    } while (false)
#else
#define SILEX_DEBUG_CHECK(ctx, module, level, label, expression) \
    do {                                                        \
        (void)sizeof(ctx);                                      \
    } while (false)
#endif

#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
#define SILEX_PROFILE_SCOPE(ctx, module, label)                                  \
    ::silex::ProfileScope SILEX_DETAIL_CONCAT(silex_profile_scope_, __LINE__)(   \
        (ctx), (module), __func__, (label))

#define SILEX_PROFILE_EVENT(ctx, module, label)                                  \
    do {                                                                         \
        const ::silex::DiagnosticsContext* const silex_profile_ctx = (ctx);       \
        if (::silex::profiling_enabled(silex_profile_ctx, (module))) {           \
            ::silex::profile_emit(silex_profile_ctx, (module),                   \
                                  ::silex::ProfileEvent::event, __func__,        \
                                  (label));                                      \
        }                                                                        \
    } while (false)
#else
#define SILEX_PROFILE_SCOPE(ctx, module, label)                 \
    do {                                                        \
        (void)sizeof(ctx);                                      \
    } while (false)
#define SILEX_PROFILE_EVENT(ctx, module, label)                 \
    do {                                                        \
        (void)sizeof(ctx);                                      \
    } while (false)
#endif
