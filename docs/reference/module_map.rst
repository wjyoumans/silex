Public Header and Module Map
============================

This page maps installed headers to the conceptual guide.  Include
``<silex/silex.hpp>`` for the complete supported surface, or include focused
headers to reduce coupling.  Files below ``silex/detail/`` are implementation
details even though they may be present in a source or installed tree.

.. list-table:: Public headers
   :header-rows: 1
   :widths: 31 31 38

   * - Header
     - Primary objects or operations
     - Guide
   * - ``silex/silex.hpp``
     - Umbrella include
     - This module map
   * - ``silex/version.hpp``, ``silex/status.hpp``
     - Version and coarse status queries
     - :doc:`api_boundary`
   * - ``silex/fmpz_smat.hpp``
     - Sparse integer matrices and HNF contexts
     - :doc:`fmpz_smat`
   * - ``silex/lat.hpp``
     - Exact row lattices and reduction
     - :doc:`lat`
   * - ``silex/number_field.hpp``, ``silex/element.hpp``,
       ``silex/signature.hpp``
     - Fields, elements, and signatures
     - :doc:`number_fields`
   * - ``silex/embedding.hpp``, ``silex/archimedean.hpp``,
       ``silex/unit.hpp``
     - Embeddings and field-level unit helpers
     - :doc:`number_fields`
   * - ``silex/hom.hpp``, ``silex/aut.hpp``
     - Field/order maps and automorphisms
     - :doc:`number_fields`
   * - ``silex/order.hpp``, ``silex/order_element.hpp``
     - Orders and integral elements
     - :doc:`orders_ideals`
   * - ``silex/ideal.hpp``
     - Integral and fractional ideals
     - :doc:`orders_ideals`
   * - ``silex/prime_ideal.hpp``, ``silex/factor_base.hpp``,
       ``silex/ideal_factorization.hpp``
     - Prime decomposition and factor bases
     - :doc:`local_algebra`
   * - ``silex/relation.hpp``
     - Relation rows and matrices
     - :doc:`local_algebra`
   * - ``silex/residue_ring.hpp``, ``silex/residue_field.hpp``
     - Finite quotient arithmetic
     - :doc:`local_algebra`
   * - ``silex/factored.hpp``, ``silex/factored_element.hpp``
     - Formal products and compact elements
     - :doc:`factored_zeta`
   * - ``silex/zeta.hpp``
     - Zeta and Belabas--Friedman audit records
     - :doc:`factored_zeta`
   * - ``silex/abelian_group.hpp``, ``silex/class_group.hpp``
     - Finite presentations and class groups
     - :doc:`class_units_compact`
   * - ``silex/order_unit.hpp``
     - Order-unit groups and paired publication
     - :doc:`class_units_compact`
   * - ``silex/sunit.hpp``
     - S-class and S-unit publication
     - :doc:`sunit_groups`
   * - ``silex/diagnostics.hpp``
     - Context-owned diagnostic channels
     - :doc:`../getting_started`
   * - ``silex/flint/*.hpp``
     - Silex-owned FLINT RAII values and borrowed views
     - :doc:`api_boundary`

The generated ``silex/build_config.hpp`` records the configured library and
dependency versions.  It is a package metadata header, not a feature-selection
API.  Optional backend types remain outside the public headers.
