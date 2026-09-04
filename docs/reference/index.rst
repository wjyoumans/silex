Conceptual C++ API Guide
========================

These pages explain the public native C++ object model, contracts, and module
relationships.  They are intentionally conceptual: they do not enumerate
every overload or reproduce declarations from the headers.  Consult the
headers installed with the selected Silex release for exact signatures, and
use :doc:`module_map` to locate them.

Objects own resources through RAII, parent objects are kept alive by
value-like handles where long-lived parent relationships are needed, and
public methods report failure with explicit boolean/status outputs rather than
throwing exceptions.

.. toctree::
   :maxdepth: 1

   module_map
   fmpz_smat
   lat
   api_boundary
   number_fields
   orders_ideals
   local_algebra
   factored_zeta
   class_units_compact
   sunit_groups
   algorithms_and_sources
