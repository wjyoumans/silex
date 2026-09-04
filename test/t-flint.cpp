#include <silex/version.hpp>

#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>

#include <cassert>
#include <cstring>

namespace {

bool nonempty(const char* s) {
    return s != nullptr && s[0] != '\0';
}

void check_flint_lifecycle() {
    fmpz_t z;
    fmpq_poly_t f;

    fmpz_init(z);
    fmpq_poly_init(f);

    fmpz_set_si(z, 42);
    fmpq_poly_set_coeff_fmpz(f, 2, z);
    assert(fmpq_poly_degree(f) == 2);

    fmpq_poly_clear(f);
    fmpz_clear(z);
}

}  // namespace

int main() {
    assert(nonempty(silex::required_flint_version()));
    assert(nonempty(silex::flint_pkgconfig_version()));
    assert(nonempty(silex::flint_compile_time_version()));
    assert(nonempty(silex::flint_runtime_version()));

    assert(std::strcmp(silex::flint_compile_time_version(),
                       silex::flint_runtime_version()) == 0);
    assert(std::strcmp(silex::flint_pkgconfig_version(),
                       silex::flint_compile_time_version()) == 0);

    check_flint_lifecycle();

    return 0;
}
