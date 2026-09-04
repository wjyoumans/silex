#include <silex/silex.hpp>

int main() {
    const char* version = silex::flint_runtime_version();
    return (version == nullptr || version[0] == '\0') ? 1 : 0;
}
