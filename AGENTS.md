# Silex Agent Instructions

These instructions govern automated work in this repository.

## Authority and scope

- Follow the user's current direction first, then the documented source
  behavior, current Silex behavior, and FLINT contracts.
- Do not invent mathematical algorithms, proof rules, normalization choices,
  failure behavior, or performance heuristics.
- Ask before changing public behavior, supported scope, licensing, release
  policy, or externally visible compatibility promises.
- Treat the 0.x line as an evolving development foundation. Do not claim
  source or binary compatibility, release readiness beyond a tagged release,
  or production support unless the public documentation explicitly does so.

## Repository boundaries

- Work only in this repository for a Silex code or documentation task.
- Keep this repository focused on the native library, public documentation,
  tests, examples, native microbenchmarks, profiling and replay tools, and
  noninstalled Silex adapters.
- Do not add cross-project campaign orchestration, private corpora, reusable
  agent infrastructure, or general-purpose repository-management tooling.
- Do not commit build trees, benchmark scratch output, private source traces,
  credentials, local profiles, or generated documentation.

## Worktrees and commits

- Use one writing task, one repository, and one Git worktree unless the user
  explicitly requests a different workflow. Reviewers are read-only unless
  given a separate correction task.
- Parallel agents each own a separate worktree and task branch. Never edit,
  reset, or commit in another task's worktree or on `main` unless you are the
  integrating session. Development orchestration and task tracking live outside
  this repository.
- Keep commits focused and validated. Stage explicit paths; never use broad
  staging commands such as `git add .`.
- Do not create remotes, push, tag, publish, rewrite history, or destructively
  clean without explicit approval.

## Mathematical and source fidelity

- Trace a behavior to its documented source before changing it. Record the
  relevant upstream file, routine or range, and behavioral contract when the
  change is mathematically significant.
- Preserve parent identity, normalization, certification labels, failure
  semantics, and exactness contracts unless an approved decision changes them.
- Prefer the smallest test that establishes a mathematical claim. Add
  differential or invariant tests when porting behavior.
- Keep required attribution in documentation and legal notices. Fresh Git
  history does not erase upstream authorship or license obligations.
- Experiments must be isolated, labeled, and have an acceptance or stop
  criterion. Do not silently promote an experiment into production behavior.

## C++ and API rules

- Use C++20 and the restricted profile checked by
  `tools/check-cpp-profile.py`.
- Keep exceptions and RTTI disabled in the core unless an approved policy
  change says otherwise.
- Prefer Silex domain types, explicit ownership, RAII-managed FLINT storage,
  named expensive operations, and explicit options/results.
- Keep optional backend types out of public APIs and preserve a FLINT-backed
  default path.
- Do not add a C API or Python binding speculatively. Add adapters only for a
  concrete consumer or extraction task.
- Update source-name checks when a rename creates a newly retired active-tree
  identifier.

## Verification

Classify performance impact before editing runtime implementation, installed
inline or template code, backend routing, or performance-affecting build
options. Scale the evidence to the classification:

- **No runtime effect** (tests, docs, examples, noninstalled tools): state the
  concrete reason in the handoff. No timing is needed.
- **Runtime change without a performance claim** (features, fixes,
  refactors): run the affected native benchmark rows once as a semantic gate,
  then do a quick baseline/candidate comparison of the affected rows and one
  unaffected control. Escalate to the full comparison below if the quick
  comparison suggests a slowdown.
- **Optimization or performance claim**, or a change to shared hot paths
  (arithmetic kernels, relation search, linear algebra, LLL, backend routing):
  profile first, then run the full comparison.

The full comparison uses separate fresh Release build trees with identical
compiler, dependency, and build settings, a committed baseline, correctness
checks and a one-iteration semantic gate before timing, at least nine
repetitions of the affected rows and a control, and preserved raw rows and
metadata. A repeatable slowdown outside control drift blocks the change unless
the tradeoff is explicitly accepted. An inconclusive result is not a pass, and
there is no repository-wide percentage threshold. See
`docs/development/benchmarking.rst` and `tools/bench/native-benchmark.py`.

Run checks proportional to the change. For ordinary native changes, start with:

```sh
python3 tools/check-cpp-profile.py --root .
python3 tools/check-source-name-boundary.py --root .
cmake --preset default
cmake --build --preset default
ctest --preset default
```

For warning-sensitive code, also run:

```sh
cmake --preset werror
cmake --build --preset werror
ctest --preset werror
```

For ownership, lifetime, or memory-sensitive code, also run:

```sh
cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan
```

Before creating or moving a tag, Google Benchmark is required. Run the full
Release build and semantic suite:

```sh
cmake --preset benchmarks
cmake --build --preset benchmarks
ctest --preset benchmarks
```

Do not treat an unavailable benchmark dependency, an absent benchmark CTest,
or an unclassified native benchmark row as a skipped or passing release gate.
The suite audits the complete row inventory and executes the rows classified
for release semantics. It does not impose wall-clock timing thresholds.

Build documentation with warnings as errors when public APIs, examples, or
Sphinx sources change. Validate the installed-package consumer when CMake,
headers, exports, or dependencies change.

Performance claims require a reproducible benchmark command, raw
machine-readable output, comparable baseline and candidate metadata, and the
selected controls. Do not run timing measurements that support a claim while
other CPU-intensive work, including other agents' builds and tests, is active.

## Completion evidence

Report changed paths, exact validation commands and results, source/provenance
implications when relevant, and any remaining risk. Passing a plan checklist is
not evidence by itself.
