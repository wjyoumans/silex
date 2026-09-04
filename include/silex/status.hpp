#pragma once

#include <cstdint>

namespace silex {

enum class Status : std::int32_t {
    ok = 0,
    invalid_argument = 1,
    domain_error = 2,
    allocation_failed = 3,
    not_implemented = 4,
    backend_unavailable = 5,
    internal_error = 6,
};

enum class CertificationMode : std::int32_t {
    unknown = 0,
    heuristic = 1,
    grh = 2,
    proven = 3,
};

enum class ProofState : std::int32_t {
    not_checked = 0,
    unavailable = 1,
    verified = 2,
};

constexpr bool ok(Status s) noexcept {
    return s == Status::ok;
}

const char* status_message(Status s) noexcept;

}  // namespace silex
