#include <silex/silex.hpp>

#include <cassert>
#include <cstring>

int main() {
    assert(silex::ok(silex::Status::ok));
    assert(!silex::ok(silex::Status::invalid_argument));
    assert(std::strcmp(silex::status_message(silex::Status::ok), "ok") == 0);
    assert(std::strcmp(silex::status_message(silex::Status::not_implemented),
                       "not implemented") == 0);
    assert(silex::version_string() != nullptr);
    assert(std::strcmp(silex::version_string(), "0.1.0") == 0);
    assert(silex::version_major() == 0);
    assert(silex::version_minor() == 1);
    assert(silex::version_patch() == 0);
    assert(silex::required_flint_version() != nullptr);
    assert(silex::flint_compile_time_version() != nullptr);
    assert(silex::flint_runtime_version() != nullptr);
    return 0;
}
