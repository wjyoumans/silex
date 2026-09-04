Sparse Integer Matrix API
=========================

The ``fmpz_smat`` module stores sparse integer rows and row-major sparse
integer matrices for exact row-module computations.  It also exposes the
persistent modular-rank context and the exact incremental HNF/index context
used by later lattice, relation, and class-group code.

The public surface is the native C++ API in ``silex/fmpz_smat.hpp``.

Native C++ API
--------------

The C++ API exposes move-only domain objects in the ``silex::fmpz_smat``
namespace:

``SparseRow``
    Owns one canonical sparse integer row.  Stored column indices are strictly
    increasing and zero values are not stored.

``SparseMat``
    Owns a row-major sparse integer matrix with a fixed nonnegative column
    count.

``ModRankContext``
    Maintains persistent finite-field row rank state for one ambient dimension
    and prime modulus.

``HnfContext``
    Maintains exact incremental row-module HNF and full-rank index state.

The native classes use explicit named operations such as ``set_entry``,
``append_row``, ``mul_fmpz_mat``, ``rank_mod_prime_ui``, and ``full_rank_index``.
Owned Silex objects are RAII-managed and move-only.  FLINT matrices and values
passed through public methods are borrowed for the duration of the call and
copied when stored.

Native callers should use Silex FLINT RAII wrappers for caller-owned FLINT
storage.  Sparse-matrix APIs accept wrapper refs/views; raw handles remain
available at direct FLINT call boundaries:

.. code-block:: cpp

   #include <silex/fmpz_smat.hpp>
   #include <silex/flint/fmpz.hpp>
   #include <silex/flint/fmpz_mat.hpp>
   #include <silex/flint/nmod_mat.hpp>

   silex::flint::FmpzMat dense(3, 3);
   fmpz_set_si(fmpz_mat_entry(dense.raw(), 0, 0), 2);
   fmpz_set_si(fmpz_mat_entry(dense.raw(), 1, 1), 3);
   fmpz_set_si(fmpz_mat_entry(dense.raw(), 2, 2), 5);

   silex::fmpz_smat::SparseMat matrix(0);
   matrix.set_fmpz_mat(dense);

   silex::flint::FmpzMat exported(matrix.nrows(), matrix.ncols());
   matrix.get_fmpz_mat(exported);

   silex::flint::NmodMat mod(matrix.nrows(), matrix.ncols(), 7);
   matrix.get_nmod_mat(mod);

   slong rank = 0;
   const bool ok = matrix.rank_mod_prime_ui(&rank, 7);

``SparseRow`` and ``SparseMat`` are move-only.  Use ``set`` for an explicit
deep copy and move construction or move assignment for ownership transfer.
``ModRankContext`` and ``HnfContext`` are also move-only and keep their
algorithm state inside the object.

Dense Import and Export
~~~~~~~~~~~~~~~~~~~~~~~

``SparseRow::set_fmpz_mat_row`` imports one row from a caller-owned matrix.
``SparseRow::get_fmpz_mat_row`` clears and writes one row of an initialized
dense output matrix.  The output matrix must have enough columns for every
stored sparse entry.  ``SparseRow::to_fmpz_mat_row`` is the owned-copy
convenience form when caller code does not need to reuse dense storage.

``SparseMat::set_fmpz_mat`` replaces the sparse matrix with a sparse copy of a
caller-owned dense matrix.  ``SparseMat::get_fmpz_mat`` writes exactly
``nrows() x ncols()`` entries to an initialized dense output matrix.
``SparseMat::to_fmpz_mat`` returns the same dense export as an owned
``silex::flint::FmpzMat``.

``SparseMat::mul_fmpz_mat`` computes sparse-times-dense multiplication.  The
dense right operand must have one row per sparse-matrix column.  The output
matrix must already have dimensions ``nrows() x fmpz_mat_ncols(right)``.  The
output may alias the dense right operand when those dimensions are compatible.

Modular Rank and HNF Contexts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SparseMat::rank_mod_prime_ui`` computes the rank of a sparse matrix over
``Z / p Z`` using FLINT dense modular rank.  The modulus ``p`` must be prime.
Non-prime moduli return ``false`` and leave the caller rank output unchanged.

``ModRankContext::set_prime_ui`` configures a nonnegative ambient dimension and
prime modulus.  ``add_row`` and ``add_fmpz_mat_row`` set ``*independent`` to
``true`` exactly when the row increases rank over the configured finite field.
``reset`` clears row state while preserving the configured dimension and prime.

``HnfContext::reset`` configures the exact incremental row-module context and
an optional modular prefilter.  Passing ``p = 0`` disables the modular filter;
passing a nonzero prime enables it.  ``add_row`` and ``add_fmpz_mat_row`` set
``*independent`` to ``true`` exactly when rational rank increases.  Dependent
rows may still refine the integer row module.  ``get_hnf`` copies the current
exact row HNF, and ``full_rank_index`` succeeds only when rank equals ambient
dimension.

Implementation Lineage
----------------------

The surface retains established Silex sparse-row invariants and uses FLINT
``fmpz``, ``fmpz_mat``, and ``nmod_mat`` contracts for dense and modular
conversion, transpose, multiplication, rank, and exact HNF/index contexts.
The public implementation exposes native C++ domain objects rather than
private sparse triangular optimization details.  See
:doc:`algorithms_and_sources` for the public source map.
