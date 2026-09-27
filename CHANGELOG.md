# Changelog

Notable user-visible changes since the last tagged release. Entries move to
[`RELEASE_NOTES.md`](RELEASE_NOTES.md) when a release is prepared. This file
follows the spirit of [Keep a Changelog](https://keepachangelog.com/), scoped
to Silex's own conventions (see `docs/development/source_fidelity.rst` for
what a durable mathematical change record requires and where it lives).

## Unreleased

### Fixed

- `Element::is_power` and `Element::is_square` now find non-integral roots of
  norm `+-1` that were previously reported unsupported, for example
  `((3 + 4i) / 5)^3` in `Q(i)`. The Hensel reconstruction only recovers an
  integral candidate; Silex now rescales a non-integral input `a` by its
  power-basis denominator `d`, lifts a root of `y^n = a d^n`, and returns
  `y / d`. Behavior for algebraic integers is unchanged. See the upstream
  source citation and the new leading-coefficient pre-filter in
  [`docs/reference/algorithms_and_sources.rst`](docs/reference/algorithms_and_sources.rst),
  "Element powers and roots".
