Factored Elements, Compact Presentations, and Zeta Audits
=========================================================

This page documents the native C++ APIs for formal products of number-field
elements, Hecke-style compact presentations, and zeta/class-regulator product
audit records.  These APIs are used directly in examples and internally by
class-group and order-unit proof drivers.

Factored Elements
-----------------

``silex::FactoredElement`` is a formal product of nonzero field elements and
stores a ``NumberField`` parent handle by value.  It uses the native general
factored-element model: owned factors with signed ``slong`` exponents.

The implemented operations are:

* ``one`` and ``set_element``;
* ``push`` and ``normalize``;
* structural multiplication, division, inversion, signed powers, and exact
  roots;
* exact square and positive signed-power detection;
* explicit evaluation back to ``Element``;
* direct logarithmic embedding accumulation without expanding the product.

The formal product is structural.  Normalization removes zero-exponent
entries, but it does not merge duplicate factors or compare factors by exact
field equality.  Expansion is explicit through ``evaluate``.  Outputs publish
only after a full operation succeeds, preserving previous output values on
documented failure paths.

``FactoredRootStrategy`` selects the implemented power-root route:

``automatic`` and ``pari_reduced``
   Use the PARI-reduced-style exponent splitting path.

``hecke_compact``
   Build a ``CompactElement`` and use the compact base-root path.

The native API intentionally uses signed ``slong`` exponents for now.
Arbitrary-size exponent/base APIs are an intentional future design question,
not part of the current native C++ surface.

Compact Elements
----------------

``silex::CompactElement`` stores a ``NumberField`` parent handle, a positive
signed-word base, and factored coefficients.  It represents a compact product
over powers of that base and supports:

* ``set_factored_element`` to split a factored element into compact
  coefficients;
* ``evaluate`` to expand the compact presentation;
* ``logarithmic_embedding`` to accumulate weighted coefficient embeddings;
* ``root_base`` to extract the base-th root when the leading coefficient and
  higher compact coefficients certify one.

The compact representation is a performance and proof-driver data structure,
not a hidden representation of ordinary ``Element`` arithmetic.  Use explicit
``evaluate`` or embedding calls when a concrete field element or analytic row
is needed.

``examples/factored_compact_roots.cpp`` is the minimal compiled example.  It
builds a formal product in ``Q``, converts it to a compact presentation,
checks a Hecke-compact power root, and verifies the structural root by
evaluation.

Zeta and Class-Regulator Product APIs
-------------------------------------

``silex/zeta.hpp`` provides zeta residue, log residue, and class-regulator
product helpers for:

* degree-one orders;
* maximal quadratic orders with conductor-one metadata through FLINT
  Dirichlet ``L(1, chi_D)`` support;
* higher-degree maximal orders through the Belabas-Friedman fallback.

The public helpers include both direct value-only calls and explicit
Belabas-Friedman audit calls:

``zeta_log_residue_bf_audit``
   Compute a log residue with audit metadata.

``zeta_residue_bf_audit``
   Compute a residue with audit metadata.

``zeta_class_regulator_product_bf_audit``
   Compute the class-regulator product with audit metadata.

Each audit API has an explicit output-buffer form and an owned-return
``std::optional<ZetaBfAuditResult>`` form.  ``ZetaBfAuditResult`` stores:

``value``
   The computed Arb ball value.

``error_bound``
   The certified error-bound Arb ball.

``cutoff``
   The final Belabas-Friedman cutoff used by the computation.

``work_precision``
   The Arb working precision selected by the computation.

``examples/maximal_order_zeta.cpp`` is the minimal compiled example.  It
constructs a cubic maximal order, reads the order index and discriminant, and
prints the owned BF audit record returned by
``zeta_class_regulator_product_bf_audit``.

Audit Records Are Not Certificates
----------------------------------

The BF audit result records an analytic computation.  It is not by itself a
class-group certificate.  Class-group certification consumes analytic values
through explicit proof gates, records the BF metadata on the
``ClassGroupContext`` when that path is used, and publishes
``CertificationMode::proven`` only if the required class/unit proof components
verify.

Use zeta/BF APIs directly when you need the analytic value and audit metadata.
Use ``ClassGroupContext`` and ``OrderUnitGroup`` certification metadata when
you need to know whether a class group or unit group result has been proven.

Diagnostics and Benchmarks
--------------------------

Factored-root and compact-representation hot paths accept optional
``DiagnosticsContext`` pointers in the operations that are used by proof
drivers.  Profiling callbacks are for attribution and add overhead.  Use the
Google Benchmark targets for timing comparisons:

* ``b-silex-nf_fac_elt`` for factored/compact microbenchmarks;
* ``b-silex-nf_clgp`` for class/unit proof-driver rows that consume zeta
  and compact witnesses.

The factored/compact benchmark target covers compact-root workloads.  Use
task-local machine-readable output and a comparable baseline before making a
performance claim or changing an algorithmic path.

There is no standalone ``nf_zeta`` Google Benchmark target today.  Zeta and
Belabas-Friedman proof-driver performance should therefore be measured through
class-unit rows in ``b-silex-nf_clgp`` until a standalone zeta workload is
needed.

Source Lineage
--------------

PARI formal-product behavior and Hecke ``FactoredElem`` and compact
representation routines provide the factored/compact baselines.  PARI's
analytic class-regulator behavior and the FLINT/Arb contracts provide the zeta
and interval baseline.  See :doc:`algorithms_and_sources` for the public
file-level map.
