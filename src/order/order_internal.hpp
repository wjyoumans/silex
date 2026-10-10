#pragma once

#include <silex/order.hpp>

namespace silex::detail {

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
};

}  // namespace silex::detail
