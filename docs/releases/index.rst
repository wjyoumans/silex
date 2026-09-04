Releases
========

Silex publishes versioned release notes alongside two documentation channels:

``stable``
   Built from the canonical release tag.  Use this channel when reproducing a
   released result or relying on a versioned support matrix.

``dev``
   Built from the development branch.  It may describe interfaces intended for
   the next release and carries no compatibility promise beyond its exact
   revision.

The page header reports both the release value and channel embedded by the
documentation build.  Local and hosted builds set them with
``SILEX_DOCS_RELEASE`` and ``SILEX_DOCS_CHANNEL``; the latter accepts only
``stable`` or ``dev``.

.. toctree::
   :maxdepth: 1
   :caption: Versioned notes

   0.1.0
