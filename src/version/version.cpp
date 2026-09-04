#include <silex/version.hpp>

#include <silex/build_config.hpp>

#include <flint/flint.h>

namespace silex {

int version_major() noexcept {
    return build_config::version_major;
}

int version_minor() noexcept {
    return build_config::version_minor;
}

int version_patch() noexcept {
    return build_config::version_patch;
}

const char* version_string() noexcept {
    return build_config::version;
}

const char* required_flint_version() noexcept {
    return build_config::required_flint_version;
}

const char* flint_pkgconfig_version() noexcept {
    return build_config::flint_pkgconfig_version;
}

const char* flint_compile_time_version() noexcept {
    return FLINT_VERSION;
}

const char* flint_runtime_version() noexcept {
    return flint_version;
}

}  // namespace silex
