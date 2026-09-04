# Silex 0.1.1 release notes

**Released: 2026-09-03**

Silex 0.1.1 is the first public development-foundation release of the native
C++ library. The release is identified by the public
[`0.1.1` tag](https://github.com/wjyoumans/silex/releases/tag/0.1.1), and its
versioned documentation is available at
<https://wjyoumans.github.io/silex/0.1.1/>.

## Highlights

Silex provides exact algebraic-number-theory operations over absolute simple
number fields over `QQ`, backed by FLINT. The initial public foundation includes:

- fields, elements, orders, maximal orders, ideals, and residue arithmetic;
- exact lattices, finite abelian groups, factor bases, and relations;
- compact elements, units, regulators, class groups, and explicit proof
  metadata;
- one failure-atomic paired class/unit transaction with explicit `proven` and
  `grh` requests; and
- S-class and S-unit publication from already proven ordinary class/unit
  inputs.

The exact qualified platform, package, and mathematical presentations are in
the [0.1.1 support matrix](https://wjyoumans.github.io/silex/0.1.1/support_matrix.html).
The qualified stack is Linux x86_64 with GCC 16, Clang 22, FLINT 3.6, and
CMake 4.4. CMake 3.20 and FLINT 3.0.0 are declared minima, not qualified
version ranges.

## AI-assisted development disclosure

Silex and its companion repositories, Silex Bench and Silex Devtools, were
built almost entirely with OpenAI Codex, initially using GPT-5.5 and later
GPT-5.6, under the direction and review of William Youmans.

The companion repositories are:

- [Silex Bench](https://github.com/wjyoumans/silex-bench), for
  correctness-gated cross-implementation comparisons and performance
  campaigns; and
- [Silex Devtools](https://github.com/wjyoumans/silex-devtools), for the
  maintained Codex implementation, performance, review, and release workflows.

Comparative checks supplement rather than replace Silex's native tests.

## Build and package

The supported package is the default FLINT-backed configuration. It exports
the CMake target `Silex::silex`, honors nondefault `GNUInstallDirs` layouts,
and installs `LICENSE`, `NOTICE.md`, `THIRD_PARTY_NOTICES.md`, and
`CITATION.cff` under `CMAKE_INSTALL_DOCDIR`.

Release tests cover the default installed package, a nondefault install layout,
a clean `add_subdirectory` consumer, and a clean tracked source archive.
Exported package files are checked for repository source/build path leakage.
Optional fplll and flatter integrations remain unqualified source-tree
development options.

## Certification and compatibility

`CertificationMode::proven` is published only after its route-specific proof
contract succeeds. The `grh` request/result label is preserved, but GRH
certification semantics and coverage remain work in progress; successful
`grh` results are conditional and provisional rather than unconditional
component proofs.

Silex 0.x provides neither source-compatibility nor binary-ABI guarantees.
Consumers should rebuild against the exact version they use. This release does
not claim production or long-term-support readiness. See [SUPPORT.md](SUPPORT.md)
and [SECURITY.md](SECURITY.md) for the best-effort latest-release policies.

## Deferred beyond 0.1.1

- blanket higher-degree or unlisted quadratic paired class/unit completion;
- optional-backend installed-package qualification;
- relative fields, ray class groups, class fields, and general Galois groups;
- C, Python, or Sage bindings; and
- cross-engine targets and campaigns within the native repository; those tools
  are maintained separately in Silex Bench.

## Source and licensing

Silex is distributed under [GPL-3.0-or-later](LICENSE). Project notices are in
[NOTICE.md](NOTICE.md). Algorithm lineage and source anchors are documented in
the [public manual](https://wjyoumans.github.io/silex/0.1.1/algorithms_and_sources.html),
and third-party attribution and compatible license obligations are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
