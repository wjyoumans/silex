Developer Guide
===============

These pages collect the contributor and maintainer policies that are distinct
from the conceptual public API guide.  They describe how to preserve
source-backed mathematical behavior, the restricted C++ profile, validation,
and performance evidence.

Install the documentation dependency, then use the strict documentation
presets:

.. code-block:: console

   python3 -m pip install -r docs/requirements.txt
   cmake --preset docs
   cmake --build --preset docs
   cmake --build --preset docs-linkcheck

For an explicit local channel/version override, invoke Sphinx directly:

.. code-block:: console

   SILEX_DOCS_RELEASE=0.1.0 SILEX_DOCS_CHANNEL=dev \
     sphinx-build -n -W --keep-going -b html docs build/docs/html

The CMake ``docs-html`` target runs the same warning-as-error build.  The
dependency pin is shared by local publication checks and the hosted-docs
workflow.

.. toctree::
   :maxdepth: 1

   source_fidelity
   architecture
   instrumentation
   testing
   benchmarking
