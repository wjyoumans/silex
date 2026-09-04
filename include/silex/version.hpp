#pragma once

namespace silex {

int version_major() noexcept;
int version_minor() noexcept;
int version_patch() noexcept;
const char* version_string() noexcept;
const char* required_flint_version() noexcept;
const char* flint_pkgconfig_version() noexcept;
const char* flint_compile_time_version() noexcept;
const char* flint_runtime_version() noexcept;

}  // namespace silex
