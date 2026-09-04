# Silex — computational algebraic number theory in C++

[![CI](https://github.com/wjyoumans/silex/actions/workflows/ci.yml/badge.svg)](https://github.com/wjyoumans/silex/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/wjyoumans/silex?sort=semver)](https://github.com/wjyoumans/silex/releases/latest)
[![Documentation](https://img.shields.io/badge/docs-0.1.0-blue)](https://wjyoumans.github.io/silex/0.1.0/)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue)](LICENSE)

Silex is a native C++20 library for computational algebraic number theory. It
provides domain objects and algorithms over FLINT-backed exact arithmetic.
This project is unrelated to the website builder and other software projects
that also use the name Silex.

Version 0.1.0 is the first public development-foundation release. The 0.x API
and object layout may change without source- or binary-compatibility guarantees;
consumers should rebuild against the exact version they use. See the
[0.1.0 support matrix](https://wjyoumans.github.io/silex/0.1.0/support_matrix.html)
for the qualified platform, package configuration, mathematical scope, and
limitations.

## AI-assisted development disclosure

Silex and its companion repositories, Silex Bench and Silex Devtools, were
built almost entirely with OpenAI Codex, initially using GPT-5.5 and later
GPT-5.6, under the direction and review of William Youmans.

## Project family

- [Silex](https://github.com/wjyoumans/silex) is the native computational
  algebraic number theory library. Its in-tree benchmarks cover native
  components and release semantics.
- [Silex Bench](https://github.com/wjyoumans/silex-bench) provides
  correctness-gated comparisons and performance campaigns across Silex and
  other mathematical systems. These checks supplement rather than replace
  Silex's native tests.
- [Silex Devtools](https://github.com/wjyoumans/silex-devtools) provides the
  maintained Codex workflows used to implement, measure, review, and release
  Silex changes.

## Documentation

- [Documentation home](https://wjyoumans.github.io/silex/)
- [Stable 0.1.0 documentation](https://wjyoumans.github.io/silex/0.1.0/)
- [Development documentation](https://wjyoumans.github.io/silex/dev/)
- [Release notes](RELEASE_NOTES.md)

## Quick start

The declared minima are CMake 3.20, a C++20 compiler, `pkg-config`, FLINT
3.0.0, and Python 3 when tests are enabled. Documentation builds additionally
require Python 3.12 or newer for the pinned Sphinx 9.1.0 dependency. On a
system with the core dependencies:

```sh
git clone https://github.com/wjyoumans/silex.git
cd silex
cmake --preset default
cmake --build --preset default
ctest --preset default
```

Examples are enabled by default. To run one:

```sh
cmake --build --preset default --target example-field-order-ideal-basics
./build/default/examples/example-field-order-ideal-basics
```

To install the default FLINT-backed package:

```sh
cmake -S . -B build/install \
  -DCMAKE_BUILD_TYPE=Release \
  -DSILEX_BUILD_TESTS=OFF \
  -DSILEX_BUILD_EXAMPLES=OFF
cmake --build build/install
cmake --install build/install --prefix "$PWD/build/install-prefix"
```

The installed package exports `Silex::silex`:

```cmake
find_package(Silex 0.1.0 EXACT CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE Silex::silex)
```

The umbrella header exposes the supported native surface:

```cpp
#include <silex/silex.hpp>

#include <iostream>

int main() {
    std::cout << "Silex " << silex::version_string() << '\n';
}
```

The supported 0.1.0 package is the default FLINT-backed configuration. fplll
and flatter remain source-tree development options; their installed-package
configurations are not qualified for 0.1.0.

## Development

- [Examples](examples/README.md)
- [Native benchmarks](bench/README.md)
- [Cross-implementation benchmarks](https://github.com/wjyoumans/silex-bench)
- [Codex development workflows](https://github.com/wjyoumans/silex-devtools)
- [Contributing guide](CONTRIBUTING.md)
- [Support policy](SUPPORT.md)
- [Security policy](SECURITY.md)
- [Code of Conduct](CODE_OF_CONDUCT.md)

Public contributions are welcome. Mathematical changes must retain their
source lineage, exactness, normalization, parent, failure, and certification
contracts. The contributor guide describes the required validation and review
evidence.

## Citation, license, and provenance

Citation metadata is available in [CITATION.cff](CITATION.cff). Silex is
distributed under [GPL-3.0-or-later](LICENSE). Copyright and project notices
are in [NOTICE.md](NOTICE.md); upstream attribution and compatible license
obligations are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
