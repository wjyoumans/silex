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
