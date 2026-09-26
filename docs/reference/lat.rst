Lattice API
===========

The ``lat`` module represents integer row lattices in an ambient
:math:`\mathbf Z^n`.  It is the lattice substrate used by order, ideal, and
relation code.

The public surface is the native C++ API in ``silex/lat.hpp``.

Native C++ API
--------------

The C++ API exposes the move-only ``silex::lat::Lat`` domain object.  A
``Lat`` owns an integer row basis in a fixed ambient dimension and provides
explicit named operations for expensive lattice work:

``set_basis`` / ``get_basis``
    Copy basis data between the lattice and caller-owned
    ``silex::flint::FmpzMat`` storage or borrowed matrix refs.

``basis``
    Return an owned ``silex::flint::FmpzMat`` copy of the stored basis for
    ordinary native C++ callers that do not need to reuse an output buffer.

``hnf`` / ``hnf_transform``
    Compute trimmed row HNF, optionally with the transform matrix.

``contains``, ``sum``, ``intersection``, ``index``
    Perform exact row-lattice membership and finite-index operations.

``saturate``, ``lll_reduce``, ``enum_short_vectors_arb``
    Provide the implemented saturation, search-basis, and short-vector
    enumeration paths.

``check``
    Runs lattice diagnostics through an optional ``DiagnosticsContext``.

FLINT inputs and outputs passed through these methods are borrowed for the
duration of the call.  Stored basis data is owned by the lattice object.
Mutating operations use temporaries before publishing results, so documented
output aliasing with inputs is supported.

Native callers should use Silex FLINT RAII wrappers for caller-owned FLINT
storage.  The lattice API accepts wrapper refs/views; raw FLINT handles are
still available at direct FLINT call boundaries:

.. code-block:: cpp

   #include <silex/lat.hpp>
   #include <silex/flint/arb.hpp>
   #include <silex/flint/fmpz.hpp>
   #include <silex/flint/fmpz_mat.hpp>

   silex::flint::FmpzMat basis(2, 2);
   fmpz_set_si(fmpz_mat_entry(basis.raw(), 0, 0), 2);
   fmpz_set_si(fmpz_mat_entry(basis.raw(), 1, 1), 3);

   silex::lat::Lat lattice(2);
   const bool basis_ok = lattice.set_basis(basis);

   silex::lat::Lat hnf(2);
   const bool hnf_ok = lattice.hnf(hnf);

   silex::flint::Fmpz index;
   const bool index_ok = hnf.index(index, lattice);

   silex::flint::Fmpz p;
   fmpz_set_si(p.raw(), 2);
   silex::lat::Lat saturated(2);
   const bool sat_ok = lattice.saturate(saturated, p);

``Lat`` is move-only.  Use ``set`` for an explicit deep copy and move
construction or move assignment for ownership transfer.

Native Basis and Exact Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Lat::set_basis`` copies a caller-owned matrix into the lattice.  The matrix
must have exactly ``ambient_dim()`` columns.  The input matrix may alias the
current stored basis.  ``Lat::get_basis`` copies the stored basis into an
initialized output matrix with exactly ``nrows() x ambient_dim()`` entries.
``Lat::basis`` is the owned-copy convenience form of the same export.

``Lat::hnf`` computes the trimmed row HNF of a lattice.  ``hnf_transform`` also
writes a transform matrix ``U`` such that ``U * basis == hnf_basis``.  The
transform matrix must be initialized as a square matrix whose dimension is the
input row count.

``Lat::contains`` tests an ambient integer vector.  ``sum`` and
``intersection`` compute exact row-lattice sum and intersection in a common
ambient dimension.  ``index`` computes ``[L : M]`` when ``M`` has finite index
in ``L``; it returns ``false`` on domain failure and leaves the caller index
unchanged.

``Lat::saturate`` computes the ``p``-saturation for prime ``p``.  ``lll_reduce``
returns a noncanonical reduced search basis spanning the same lattice.
``enum_short_vectors_arb`` enumerates coefficient rows bounded by an Arb
squared norm and calls a user callback with borrowed coefficient storage.

Short-Vector Enumeration
~~~~~~~~~~~~~~~~~~~~~~~~

``enum_short_vectors_arb(bound_sq, max_coord, prec, callback, user)`` calls
``callback`` with a ``1 x nrows()`` coefficient row ``x`` for each nonzero
``x`` with :math:`\lVert xB\rVert^2 \le` ``bound_sq``, where :math:`B` is the
stored basis.  Coefficients are relative to the stored basis, not to its HNF.
Enumeration is a witness search:

* ``max_coord`` limits every coefficient to ``|x_i| <= max_coord``; a negative
  value means no cap.
* The zero row is never reported, and no row is reported twice in one call.
* A callback return of ``0`` stops the enumeration, and the call returns
  ``true``.
* Dependent basis rows make the Gram matrix singular, and the call returns
  ``false`` without calling ``callback``.  A negative ``bound_sq`` or an
  empty basis returns ``true`` without calling it.  A null bound or callback,
  or ``prec < 2``, returns ``false``.
* A ``false`` return can follow callbacks already delivered, for example after
  an Arb precision failure.  Increase ``prec`` to retry; a retry calls
  ``callback`` again for rows it already received.

The Gram matrix :math:`G = BB^T` and its Cholesky factor are computed in Arb
at ``prec``.  Bases with at most 32 rows whose Cholesky midpoints fit the
double range then use a double-precision route when
``0 <= max_coord <= 10000``, or when ``max_coord`` is negative and the
Fincke--Pohst coordinate bound :math:`\sqrt{b\,(G^{-1})_{ii}}` with
:math:`b` = ``bound_sq``, evaluated in double precision, is at most 10000 for
every coordinate.  That route enumerates in double arithmetic on the Arb
midpoints, widens each coordinate interval by ``1e-9`` and the bound by a
relative ``2^-40``, so it can report rows whose squared norm is slightly above
``bound_sq`` and its completeness is not certified.  Other inputs use the Arb
route, which encloses each coordinate interval with ball arithmetic and
returns ``false`` when a partial norm or pivot is undecided at ``prec``.  It
reports a candidate unless its squared norm is certainly above ``bound_sq``,
so it too can report rows slightly past the bound.  The double route never
falls back to the Arb route after it has delivered a callback.

A capped call stays on the double route even when its uncapped coordinate
interval would lie outside the machine-integer range, because the interval is
clamped to ``max_coord`` first.  Earlier versions sent such calls to the Arb
route before any callback; they now get the non-certified double route, whose
2^-40 relative slack can admit additional rows when ``bound_sq`` is large.

Implementation Lineage
----------------------

The exact row-lattice operations retain established Silex behavior and use
FLINT HNF/SNF, exact solve, and determinant contracts.  LLL and short-vector
paths additionally follow their documented FLINT, Arb, and optional-backend
contracts.  See :doc:`algorithms_and_sources` for the public source map.
