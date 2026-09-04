# Contributing to Silex

Public contributions are welcome. Silex is an early native C++ library for
computational algebraic number theory, so changes should favor mathematical
fidelity, a coherent C++ API, and small reviewable steps over speculative
breadth.

Participation in this project is governed by the
[Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).

## Before you start

Search the [issue tracker](https://github.com/wjyoumans/silex/issues) before
opening a report. For substantial work, open an issue first so the intended
scope and source basis can be agreed before implementation. Maintainer review
is required before changing:

- public mathematical behavior or certification;
- supported compilers, platforms, dependencies, or backends;
- licensing, attribution, publication, or release policy; or
- compatibility guarantees or installed package names.

Use the [bug report](https://github.com/wjyoumans/silex/issues/new?template=bug.yml)
or [feature request](https://github.com/wjyoumans/silex/issues/new?template=feature.yml)
form where appropriate. Security reports must follow [SECURITY.md](SECURITY.md),
not the public issue tracker.

## Development setup

The required tools are a C++20 compiler, CMake 3.20 or newer, `pkg-config`,
FLINT 3.0.0 or newer, and Python 3 for the test tooling. Documentation builds
additionally require Python 3.12 or newer for the pinned Sphinx dependency.
Configure and test the baseline with:

```sh
cmake --preset default
cmake --build --preset default
ctest --preset default
```

Optional documentation and benchmark dependencies are described in the
[README](README.md), the [documentation](https://wjyoumans.github.io/silex/dev/),
and the [development guide](docs/development/index.rst).

## AI-assisted contributions and companion tools

AI-assisted and AI-generated contributions are welcome under the same
requirements as other contributions. The human submitter must understand,
review, and be able to explain the submitted change, and remains accountable
for its correctness, mathematical source basis, licensing and attribution,
tests, and performance evidence. A declaration of the tool or model and prompt
transcripts are not required.

[Silex Devtools](https://github.com/wjyoumans/silex-devtools) provides the
recommended Codex workflows for implementation, performance work, review, and
release validation. It is optional, and this repository's [AGENTS.md](AGENTS.md)
remains authoritative. [Silex Bench](https://github.com/wjyoumans/silex-bench)
provides cross-implementation correctness and performance evidence where it is
useful; that evidence supplements rather than replaces focused native tests.

## Implementation expectations

- Locate definitions, callers, tests, and documented source lineage before
  editing behavior.
- Preserve exactness, normalization, parent identity, failure semantics, and
  certification contracts.
- Prefer Silex domain types and RAII-managed FLINT storage at public boundaries.
- Keep optional backend details private and retain the default FLINT path.
- Add focused invariant, regression, or differential tests for behavior changes.
- Do not replace a known algorithm with an uncited approximation or heuristic.
- Keep generated files, build trees, local benchmark output, credentials, and
  private evidence out of Git.
- Automated contributors and contributors using generated output must also
  follow [AGENTS.md](AGENTS.md).

Mathematically significant changes must identify the upstream version and
file, routine, or range; state the behavior preserved or intentionally changed;
and update the public source map and attribution when appropriate.

## Validation

Run the narrowest relevant tests while developing, then the full applicable
gates before requesting review. For ordinary native changes:

```sh
python3 tools/check-cpp-profile.py --root .
python3 tools/check-source-name-boundary.py --root .
cmake --preset default
cmake --build --preset default
ctest --preset default
```

Use the `werror` preset for compiler-facing changes and `asan-ubsan` for
ownership, lifetime, bounds, or memory-sensitive changes. Build the Sphinx
documentation with warnings as errors when public headers, examples, or docs
change. Validate the install-prefix consumer when package exports, installed
headers, or dependency discovery change.

Runtime changes require the performance-impact triage and evidence described
in the [benchmarking guide](docs/development/benchmarking.rst). Before
creating a release tag, maintainers run the complete Release semantic gate:

```sh
cmake --preset benchmarks
cmake --build --preset benchmarks
ctest --preset benchmarks
```

## Pull requests

Create a focused branch, make cohesive commits, and open a pull request against
`main`. The pull request should include:

- a plain-language description and linked issue;
- mathematical source/provenance implications, when applicable;
- changed behavior and compatibility implications;
- exact validation commands and results; and
- remaining limitations or risks.

Update user documentation for public API changes. Do not remove copyright,
license, or lineage notices when moving or rewriting code. By submitting a
contribution, you agree that it may be distributed under the repository's
[GPL-3.0-or-later license](LICENSE).

Questions that are not bug reports belong in
[GitHub Discussions](https://github.com/wjyoumans/silex/discussions); support
expectations are described in [SUPPORT.md](SUPPORT.md).
