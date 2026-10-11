#pragma once

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order.hpp>

#include <vector>

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
};

}  // namespace silex::detail
