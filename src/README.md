# Source layout

Change modules in vertical slices. Do not add broad abstractions before source
lineage and a parity target exist.

Current layout uses native C++ module names matching the public headers under
`include/silex/`:

```text
src/abelian_group/
src/archimedean/
src/aut/
src/class_group/
src/diagnostics/
src/element/
src/embedding/
src/factor_base/
src/factored_element/
src/fmpz_smat/
src/hom/
src/ideal/
src/ideal_factorization/
src/lat/
src/number_field/
src/order/
src/order_element/
src/order_unit/
src/prime_ideal/
src/relation/
src/residue_field/
src/residue_ring/
src/signature/
src/status/
src/sunit/
src/unit/
src/version/
src/zeta/
```

Keep FLINT wrappers local and minimal until they stabilize.

Use private implementation headers under these directories only when a module
needs to be split across multiple `.cpp` files. Installed public headers remain
under `include/silex/`.

Large implementation files may temporarily use private `.inc` fragments when a
mechanical section split is useful but a true multi-translation-unit split
would require refactoring anonymous-namespace helper dependencies. Prefer a
private internal header plus real `.cpp` files once the helper boundary is
clear. Any remaining `.inc` fragments are not public headers and should be
included only by their owning `.cpp` file.
