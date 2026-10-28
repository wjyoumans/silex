#pragma once

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/order.hpp>

#include <vector>

namespace silex {
class EmbeddingContext;
}  // namespace silex

namespace silex::detail {

// One prime ideal of a decomposition cached by
// OrderAccess::store_index_divisor_decomposition, held without an Order
// handle: the cache lives in the order's shared data, and an Ideal or
// PrimeIdeal there would own a handle to that same data and never be freed.
// `hnf` is the ideal's HNF in the order basis; `principal_generator` holds the
// 1 x n coordinates of the ideal's known principal generator when
// `has_principal_generator` is set, so a restored ideal carries the same
// generator as the one first published.
struct IndexDivisorPrimeRecord {
    flint::FmpzMat hnf{0, 0};
    flint::FmpzMat principal_generator{0, 0};
    slong ramification_index = 0;
    slong residue_degree = 0;
    bool has_principal_generator = false;
};

// The characteristic polynomial of one order basis element w_index (that of
// multiplication by w_index on the order, monic in Z[x] of the order's
// degree) and its discriminant, as stored by
// OrderAccess::store_basis_element_charpoly.
struct BasisElementCharpoly {
    slong index = 0;
    flint::FmpzPoly charpoly;
    flint::Fmpz discriminant;
};

// Conductor {x in order : x * overorder is contained in order} of `order` in
// an overorder, as an ideal of `order`.  `out` must already be defined over
// `order` (same order data); it keeps that parent.  Fails, leaving `out`
// unchanged, when the orders lie in different fields or `order` is not
// contained in `overorder`.  Port of reference v0.38.6
// src/NumFieldOrd/NfOrd/NfOrd.jl:1488 `conductor(R, S)`; see conductor.cpp.
bool order_conductor(Ideal& out,
                     const Order& order,
                     const Order& overorder) noexcept;

// Contraction of an ideal of an overorder to a suborder: `ideal` intersected
// with `suborder`, as an ideal of `suborder`.  `out` must already be defined
// over `suborder`; it keeps that parent.  Fails, leaving `out` unchanged,
// when `suborder` is not contained in the ideal's order.  Port of reference
// v0.38.6 src/NumFieldOrd/NfOrd/Ideal/Arithmetic.jl:952 `contract(A, O)`;
// see conductor.cpp.
bool ideal_contract(Ideal& out,
                    const Ideal& ideal,
                    const Order& suborder) noexcept;

// base^exponent in O/pO, as Order::p_radical computes the images of the
// Frobenius power x -> x^q on the basis of a non-equation order.  Vectors
// are coordinates in the order basis with entries in [0, p);
// `multiplication_table` is an order multiplication table (row i * n + j
// holds the coordinates of w_i * w_j) and `one` the coordinates of 1.  With
// `allow_word` set and p below 2^64, products run on word residues;
// otherwise on fmpz residues.  Both give the canonical residues in [0, p).
// Fails, leaving `out` unchanged, on mismatched lengths, a table that is not
// n^2 x n, p < 2, a negative exponent, or entries outside [0, p).  Test
// seam: tests compare the two product paths through it.
bool order_coords_pow_mod_p(flint::FmpzVec& out,
                            const flint::FmpzVec& base,
                            flint::FmpzConstRef exponent,
                            const flint::FmpzVec& one,
                            const flint::FmpzMat& multiplication_table,
                            flint::FmpzConstRef p,
                            bool allow_word) noexcept;

inline bool order_has_parented_basis(const Order& order) noexcept {
    return order.is_defined() && order.has_basis() && order.parent() != nullptr;
}

// Noninstalled access to Order internals.
class OrderAccess {
public:
    // Maximality hook.  Public code cannot assert or withdraw maximality:
    // Order records it only where it is computed (maximal_order, quadratic
    // conductor metadata, degree one).  This hook writes the shared flag
    // without any check, so every handle sharing the order's data sees the
    // change.  It exists for tests that must build states the public API
    // cannot reach, such as an order whose maximality is withdrawn after
    // dependent objects were built; library code must not use it to label an
    // order maximal.
    static void set_maximality_unchecked(Order& order,
                                         bool is_maximal) noexcept {
        order.record_maximality(is_maximal);
    }

    // Index [order : Z[theta]] of the equation order of the parent field in
    // `order`, as order_index(Order::equation_order(parent), order) computes
    // it.  The value is cached in the order's shared data on first success, in
    // the same unsynchronized way as the trace matrix and discriminant, so
    // later calls skip rebuilding the equation order.  Order::set copies the
    // cache; set_basis and the factories start a new, uncached order.  Fails
    // when `order` has no basis, the equation order is not defined (the field
    // polynomial is not monic integral), or `order` does not contain Z[theta];
    // nothing is cached on failure.  Mirrors the order index cached on first
    // use by reference v0.38.6 src/NumFieldOrd/NfOrd/NfOrd.jl:326 `index`,
    // which `is_index_divisor` (NfOrd.jl:352) reads; reference 2.17.3
    // src/basemath/base2.c:2260 `primedec_aux` likewise reads the stored
    // `nf_get_index`.
    static bool equation_order_index(flint::FmpzRef out,
                                     const Order& order) noexcept;

    // The order's basis rows (power-basis coordinates) as stored in its shared
    // data, read without the copy Order::get_basis makes.  Returns null when
    // `order` is undefined or has no basis.  The pointer stays valid while
    // the order data lives: no operation changes the basis of existing order
    // data (set_basis and the factories build fresh data).
    static const flint::FmpqMat* basis(const Order& order) noexcept;

    // The order's multiplication table as stored in its shared data, read
    // without the copy Order::multiplication_table makes: row i * n + j holds
    // the order coordinates of w_i * w_j for the basis w_0, ..., w_{n-1}.
    // Fills the cache on first use, as Order::multiplication_table does.
    // Returns null when `order` is undefined or has no basis, or the table
    // cannot be built; nothing is cached then.  The pointer stays valid while
    // the order data lives: once filled, the table of an order data is never
    // written again (set_basis, Order::set and the factories build fresh
    // data).  Callers keep a handle on that data, and must not reassign it
    // through set, set_basis, define, swap or clear, while they read the view.
    static const flint::FmpzMat* multiplication_table(Order& order) noexcept;

    // The embedding context detail::order_minkowski_embedding_rows keeps in
    // the order's shared data, or null when `order` has no basis.  It has
    // roots only once some Minkowski-row build has refined it (is_set());
    // its precision is the highest one requested so far.  The context is
    // written by those builds in the same unsynchronized way as the
    // equation order index, so a reader must not hold the reference across
    // a call that can build Minkowski rows of the same order data.  Used to
    // seed a caller's own context with the same certified roots
    // (EmbeddingContextAccess::seed_roots) instead of isolating them again.
    static const EmbeddingContext* minkowski_embedding_context(
            const Order& order) noexcept;

    // Whether the order's basis is exactly the identity matrix
    // (fmpq_mat_is_one: diagonal entries 1 with denominator 1, all other
    // entries 0), so order coordinates are power-basis coefficients.  The
    // answer is cached in the order's shared data on first use, in the same
    // unsynchronized way as the equation order index.  False when `order`
    // has no basis.
    static bool basis_is_identity(const Order& order) noexcept;

    // Full decomposition of the prime p in `order`, before any residue-degree
    // filter, as decompose_prime last stored it for this order data, in the
    // order decompose_prime published it.  Returns null when nothing is
    // stored for p.  decompose_prime stores only decompositions of maximal
    // non-equation orders taken through its index-divisor branch.  The
    // records live in the order's shared data in the same unsynchronized way
    // as the equation order index; the pointer stays valid while the order
    // data lives, since later stores for other primes never move it.
    // Order::set, set_basis and the factories start with nothing stored.
    // Mirrors the per-order index-divisor decomposition cache of reference
    // v0.38.6 src/NumFieldOrd/NfOrd/Ideal/Prime.jl:213-232 (`O.index_div`).
    static const std::vector<IndexDivisorPrimeRecord>*
    index_divisor_decomposition(const Order& order,
                                flint::FmpzConstRef p) noexcept;

    // Stores `primes` as the decomposition of p in `order`, unless one is
    // already stored for p.  Silent when `order` has no basis.
    static void store_index_divisor_decomposition(
            const Order& order,
            flint::FmpzConstRef p,
            std::vector<IndexDivisorPrimeRecord>&& primes) noexcept;

    // Whether decompose_prime's search for a monogenic generator (an element
    // whose powers 1, b, ..., b^(n-1) form a basis of `order`) has already
    // run to completion on this order data and found none.  The search reads
    // only the basis and its multiplication table, which never change for
    // one order data, so its failure is a property of the order: a later
    // search would test the same candidates and fail again.  Only failure is
    // recorded; a found generator is recomputed each time.  Kept in the
    // order's shared data in the same unsynchronized way as the equation
    // order index; Order::set, set_basis and the factories start with
    // nothing recorded.  False when `order` has no basis.
    static bool monogenic_generator_absent(const Order& order) noexcept;

    // Records that the search above ran to completion on `order` and found
    // no generator.  Silent when `order` has no basis.
    static void set_monogenic_generator_absent(const Order& order) noexcept;

    // The characteristic polynomial of the basis element w_index and its
    // discriminant, as last stored for this order data by
    // store_basis_element_charpoly, or null when none is stored for
    // `index`.  Both depend only on the basis and its multiplication table,
    // which never change for one order data.  Kept in the order's shared
    // data in the same unsynchronized way as the equation order index; the
    // pointer stays valid while the order data lives, since later stores
    // never move it.  Order::set, set_basis and the factories start with
    // nothing stored.
    static const BasisElementCharpoly* basis_element_charpoly(
            const Order& order,
            slong index) noexcept;

    // Stores `charpoly` and `discriminant` for w_index, unless an entry is
    // already stored for `index`, and returns the stored entry.  Returns
    // null, storing nothing, when `order` has no basis or `index` is not in
    // [0, degree).
    static const BasisElementCharpoly* store_basis_element_charpoly(
            const Order& order,
            slong index,
            flint::FmpzPoly&& charpoly,
            flint::Fmpz&& discriminant) noexcept;

    // Sets `out` to the residue degrees above the rational prime p last
    // stored for this order data by
    // store_defining_polynomial_residue_degrees, in the order they were
    // stored, and returns true; returns false, leaving `out` unchanged,
    // when none are stored for p.  The zeta defining-polynomial route
    // (bf_residue_degrees_from_equation_order_nmod in zeta.cpp) stores the
    // degrees it reads from the distinct-degree factorization of the
    // field's defining polynomial modulo p, so that the GRH factor-base
    // bound, the norm-bounded factor-base build and the Belabas-Friedman
    // sums over one order factor each p once.  The degrees depend only on
    // the parent field's polynomial and p, which never change for one
    // order data; that route asks only for primes its own gate admits
    // (p divides neither the order discriminant nor the equation-order
    // index), so index divisors are never stored here.  Kept in the order's
    // shared data in the same unsynchronized way as the equation order
    // index; Order::set, set_basis and the factories start with nothing
    // stored.  False when `order` has no basis.
    static bool defining_polynomial_residue_degrees(std::vector<slong>& out,
                                                    const Order& order,
                                                    ulong p) noexcept;

    // Stores `degrees` (nonempty) as the residue degrees above p, unless
    // some are already stored for p.  Silent when `order` has no basis.
    static void store_defining_polynomial_residue_degrees(
            const Order& order,
            ulong p,
            const std::vector<slong>& degrees) noexcept;
};

}  // namespace silex::detail
