#include <silex/status.hpp>

namespace silex {

const char* status_message(Status s) noexcept {
    switch (s) {
    case Status::ok: return "ok";
    case Status::invalid_argument: return "invalid argument";
    case Status::domain_error: return "domain error";
    case Status::allocation_failed: return "allocation failed";
    case Status::not_implemented: return "not implemented";
    case Status::backend_unavailable: return "backend unavailable";
    case Status::internal_error: return "internal error";
    }
    return "unknown status";
}

}  // namespace silex
