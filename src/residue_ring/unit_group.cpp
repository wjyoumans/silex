// Unit group (O/m)^x of a residue ring of a maximal order: generators, cyclic
// structure and discrete logarithm.  Internal routine; see
// unit_group_internal.hpp for the contract.
//
// Source lineage, reference v0.38.6,
// src/NumFieldOrd/NfOrd/ResidueRingMultGrp.jl:
// - `_multgrp` (line 64): factor m, build each (O/p^v)^x, join them with
//   `_direct_product` (line 1357, CRT lifts of the generators against the
//   product of the other prime powers) and publish with `snf(G, GtoR)`
//   (line 1248).  ResidueUnitGroup::set_modulus and finish.
// - `_multgrp_mod_pv` (line 116): (O/p^v)^x = G1 x G2 with G1 = (O/p)^x on
//   gen1^(exp G2) and G2 = (1 + p)/(1 + p^v); the discrete log is
//   `disc_log2`.  build_prime_power and prime_power_discrete_log.
// - `_multgrp_mod_p` (line 171) and `_primitive_element_mod_p` (line 209):
//   a primitive root of O/p and its discrete log by `pohlig_hellman`
//   (line 754) and `_pohlig_hellman_prime_power` (line 775).
//   find_field_generator and field_discrete_log.
// - `_1_plus_p_mod_1_plus_pv` (line 236) with `method = :quadratic`, which
//   runs `_iterative_method` (lines 263-330) with `base_method = :quadratic,
//   use_p_adic = false`, steps from (1 + p^k) to (1 + p^min(2k, v)) with
//   `_quadratic_method` (line 437) on `_pu_mod_pv` (line 374), joins the
//   steps with `_expand` (line 332), and puts the result in Smith normal
//   form.  build_layer, build_one_plus and one_plus_discrete_log.
//
// Deviations from the reference:
// - Method choice.  The reference's default for (O_K/p^v)^x (`method =
//   nothing`) steps with `_artin_hasse_method` and switches to
//   `_p_adic_method` once k >= 1 + e/(p - 1).  This port uses the quadratic
//   method at every step, which is the reference's `method = :quadratic`
//   path (the one its non-maximal-order code uses for O/p^m).  Every method
//   presents the same group (1 + p)/(1 + p^v), so the cyclic structure is the
//   same; the generators chosen can differ.
// - Primitive root.  `_primitive_element_mod_p` draws random elements of
//   O/p.  This port tries order elements in a fixed order (the base-p digits
//   of k = 1, 2, ... as coordinates) with the same primitivity test,
//   x^((q - 1)/l) != 1 for every prime l | q - 1, so results are
//   reproducible.
// - Discrete log in (O/p)^x.  The reference scans powers when q - 1 < 11 and
//   otherwise runs Pohlig-Hellman with baby-step giant-step in each subgroup
//   of prime order l.  This port always runs Pohlig-Hellman and replaces
//   baby-step giant-step by a scan of the l powers.  The logarithm is the
//   same; the cost is O(l) instead of O(sqrt(l)) per digit.
// - Representatives.  Products are reduced modulo p^v (modulo m at the top
//   level) by the HNF of the modulus instead of the reference's quotient-ring
//   elements and integer modulus `pnumv`.  `_expand` evaluates its relation
//   products modulo p^v rather than p^l; both give the same class modulo
//   p^l, which is all the layer discrete log reads.  The quadratic discrete
//   log skips the reference's `mod(x - 1, pv)`: coordinates that differ by an
//   element of p^upper give the same Smith coordinates.
// - Reduction of discrete-log coordinates.  The reference reduces layer
//   coordinates modulo the largest invariant; this port reduces each
//   coordinate modulo its own invariant (FiniteAbelianGroup::
//   invariant_coordinates).  Same group element.
// - Negative exponents of units are reduced modulo #(O/p^v)^x and division
//   is multiplication by x^(#(O/p^v)^x - 1), instead of quotient-ring
//   inversion.
// - Fail closed.  The iterative discrete log checks that the remaining
//   factor is 1 modulo p^v after the last layer, and every routine returns
//   false (leaving its output unchanged) where the reference would raise.
//   A non-unit input to discrete_log fails instead of raising.

#include "unit_group_internal.hpp"

#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mat.h>

#include <silex/flint/fmpz_factor.hpp>
#include <silex/ideal_factorization.hpp>

#include <utility>

namespace silex::detail {
namespace {

using flint::Fmpz;
using flint::FmpzMat;

bool element_row(FmpzMat& row, const OrderElement& element) noexcept {
    FmpzMat next(1, element.degree());
    if (!element.get_coordinates(flint::FmpzMatRef(next))) {
        return false;
    }
    row.swap(next);
    return true;
}

bool ideal_hnf(FmpzMat& out, const Ideal& ideal) noexcept {
    FmpzMat next(ideal.degree(), ideal.degree());
    if (!ideal.get_hnf(flint::FmpzMatRef(next))) {
        return false;
    }
    out.swap(next);
    return true;
}

// Reduces a 1 x n coordinate row modulo the lattice of an upper-triangular
// HNF basis, as ResidueRing::reduce_row does.
void reduce_row(FmpzMat& row, const FmpzMat& hnf) noexcept {
    const slong n = fmpz_mat_ncols(hnf.raw());
    Fmpz q;
    for (slong i = 0; i < n; ++i) {
        fmpz_fdiv_q(q.raw(), fmpz_mat_entry(row.raw(), 0, i),
                    fmpz_mat_entry(hnf.raw(), i, i));
        if (fmpz_is_zero(q.raw()) != 0) {
            continue;
        }
        for (slong j = i; j < n; ++j) {
            fmpz_submul(fmpz_mat_entry(row.raw(), 0, j), q.raw(),
                        fmpz_mat_entry(hnf.raw(), i, j));
        }
    }
}

// Coordinates of a row in the basis given by the rows of an upper-triangular
// HNF; fails when the row is not in that lattice.
bool hnf_coordinates(FmpzMat& out,
                     const FmpzMat& hnf,
                     const FmpzMat& row) noexcept {
    const slong n = fmpz_mat_ncols(hnf.raw());
    FmpzMat work(1, n);
    FmpzMat coords(1, n);
    Fmpz q;
    Fmpz r;
    fmpz_mat_set(work.raw(), row.raw());
    for (slong i = 0; i < n; ++i) {
        fmpz_fdiv_qr(q.raw(), r.raw(), fmpz_mat_entry(work.raw(), 0, i),
                     fmpz_mat_entry(hnf.raw(), i, i));
        if (fmpz_is_zero(r.raw()) == 0) {
            return false;
        }
        fmpz_set(fmpz_mat_entry(coords.raw(), 0, i), q.raw());
        if (fmpz_is_zero(q.raw()) != 0) {
            continue;
        }
        for (slong j = i; j < n; ++j) {
            fmpz_submul(fmpz_mat_entry(work.raw(), 0, j), q.raw(),
                        fmpz_mat_entry(hnf.raw(), i, j));
        }
    }
    out.swap(coords);
    return true;
}

bool reduce_element(OrderElement& element, const FmpzMat& hnf) noexcept {
    FmpzMat row(0, 0);
    if (!element_row(row, element)) {
        return false;
    }
    reduce_row(row, hnf);
    return element.set_coordinates(flint::FmpzMatConstRef(row));
}

bool set_one_mod(OrderElement& out,
                 const Order& order,
                 const FmpzMat& hnf) noexcept {
    OrderElement one(order);
    if (!one.set_si(1) || !reduce_element(one, hnf)) {
        return false;
    }
    out.swap(one);
    return true;
}

bool mul_mod(OrderElement& out,
             const OrderElement& left,
             const OrderElement& right,
             const FmpzMat& hnf) noexcept {
    if (left.parent() == nullptr) {
        return false;
    }
    OrderElement product(*left.parent());
    if (!product.multiply(left, right) || !reduce_element(product, hnf)) {
        return false;
    }
    out.swap(product);
    return true;
}

// base^exponent modulo the HNF lattice, exponent >= 0.
bool pow_mod(OrderElement& out,
             const OrderElement& base,
             const Fmpz& exponent,
             const FmpzMat& hnf) noexcept {
    if (base.parent() == nullptr || fmpz_sgn(exponent.raw()) < 0) {
        return false;
    }
    OrderElement result;
    OrderElement reduced_base(*base.parent());
    if (!set_one_mod(result, *base.parent(), hnf) ||
        !reduced_base.set(base) || !reduce_element(reduced_base, hnf)) {
        return false;
    }
    for (slong bit = static_cast<slong>(fmpz_bits(exponent.raw())) - 1;
         bit >= 0; --bit) {
        if (!mul_mod(result, result, result, hnf)) {
            return false;
        }
        if (fmpz_tstbit(exponent.raw(), static_cast<ulong>(bit)) != 0 &&
            !mul_mod(result, result, reduced_base, hnf)) {
            return false;
        }
    }
    out.swap(result);
    return true;
}

// base^exponent for a unit base and any integer exponent: the exponent is
// reduced modulo the order of the unit group first.
bool pow_mod_unit(OrderElement& out,
                  const OrderElement& base,
                  const Fmpz& exponent,
                  const Fmpz& unit_count,
                  const FmpzMat& hnf) noexcept {
    Fmpz reduced;
    fmpz_fdiv_r(reduced.raw(), exponent.raw(), unit_count.raw());
    return pow_mod(out, base, reduced, hnf);
}

bool ideal_power(Ideal& out, const Ideal& base, slong exponent) noexcept {
    if (base.parent() == nullptr || exponent < 1) {
        return false;
    }
    Ideal power(*base.parent());
    if (!power.set(base)) {
        return false;
    }
    for (slong i = 1; i < exponent; ++i) {
        Ideal next(*base.parent());
        if (!next.multiply(power, base)) {
            return false;
        }
        power.swap(next);
    }
    out.swap(power);
    return true;
}

bool field_is_one(const ResidueFieldElement& element) noexcept {
    if (element.parent() == nullptr) {
        return false;
    }
    ResidueFieldElement one(*element.parent());
    return one.one() && element.equal(one);
}

// Reference `_pohlig_hellman_prime_power` scans the l powers of the order-l
// element g_ in place of `baby_step_giant_step`.
bool scan_prime_order(Fmpz& out,
                      const ResidueFieldElement& generator,
                      const Fmpz& ell,
                      const ResidueFieldElement& target) noexcept {
    ResidueFieldElement current(*generator.parent());
    if (!current.one()) {
        return false;
    }
    Fmpz j;
    for (fmpz_zero(j.raw()); fmpz_cmp(j.raw(), ell.raw()) < 0;
         fmpz_add_ui(j.raw(), j.raw(), 1)) {
        if (current.equal(target)) {
            fmpz_set(out.raw(), j.raw());
            return true;
        }
        if (!current.multiply(current, generator)) {
            return false;
        }
    }
    return false;
}

// Discrete log of a nonzero residue to the primitive root of the component:
// reference `pohlig_hellman` (ResidueRingMultGrp.jl:754) and
// `_pohlig_hellman_prime_power` (line 775).
bool field_discrete_log(Fmpz& out,
                        const ResidueFieldElement& element,
                        const ResidueUnitPrimePower& c) noexcept {
    const ResidueFieldElement& g = c.field_generator_residue;
    if (element.is_zero()) {
        return false;
    }
    if (fmpz_is_one(c.field_unit_count.raw()) != 0) {
        fmpz_zero(out.raw());
        return true;
    }

    const ResidueField& field = *g.parent();
    Fmpz result;
    Fmpz modulus;
    fmpz_one(modulus.raw());
    for (std::size_t f = 0; f < c.field_unit_primes.size(); ++f) {
        const Fmpz& ell = c.field_unit_primes[f];
        const slong v = c.field_unit_exponents[f];
        Fmpz ell_v;
        Fmpz cofactor;
        Fmpz ell_v_minus_one;
        fmpz_pow_ui(ell_v.raw(), ell.raw(), static_cast<ulong>(v));
        fmpz_pow_ui(ell_v_minus_one.raw(), ell.raw(),
                    static_cast<ulong>(v - 1));
        fmpz_divexact(cofactor.raw(), c.field_unit_count.raw(), ell_v.raw());

        ResidueFieldElement g_r(field);
        ResidueFieldElement h(field);
        ResidueFieldElement g_low(field);
        ResidueFieldElement probe(field);
        ResidueFieldElement correction(field);
        if (!g_r.pow_fmpz(g, flint::FmpzConstRef(cofactor)) ||
            !h.pow_fmpz(element, flint::FmpzConstRef(cofactor)) ||
            !g_low.pow_fmpz(g_r, flint::FmpzConstRef(ell_v_minus_one)) ||
            !probe.pow_fmpz(h, flint::FmpzConstRef(ell_v_minus_one))) {
            return false;
        }

        Fmpz a;
        Fmpz digit;
        Fmpz neg;
        if (!scan_prime_order(a, g_low, ell, probe)) {
            return false;
        }
        fmpz_neg(neg.raw(), a.raw());
        if (!correction.pow_fmpz(g_r, flint::FmpzConstRef(neg)) ||
            !h.multiply(h, correction)) {
            return false;
        }

        Fmpz ell_i;
        Fmpz remaining;
        fmpz_one(ell_i.raw());
        fmpz_set(remaining.raw(), ell_v_minus_one.raw());
        for (slong i = 1; i < v; ++i) {
            fmpz_mul(ell_i.raw(), ell_i.raw(), ell.raw());
            fmpz_divexact(remaining.raw(), remaining.raw(), ell.raw());
            if (!probe.pow_fmpz(h, flint::FmpzConstRef(remaining)) ||
                !scan_prime_order(digit, g_low, ell, probe)) {
                return false;
            }
            fmpz_mul(digit.raw(), digit.raw(), ell_i.raw());
            fmpz_add(a.raw(), a.raw(), digit.raw());
            fmpz_neg(neg.raw(), digit.raw());
            if (!correction.pow_fmpz(g_r, flint::FmpzConstRef(neg)) ||
                !h.multiply(h, correction)) {
                return false;
            }
        }

        fmpz_mod(a.raw(), a.raw(), ell_v.raw());
        Fmpz combined;
        fmpz_CRT(combined.raw(), result.raw(), modulus.raw(), a.raw(),
                 ell_v.raw(), 0);
        result.swap(combined);
        fmpz_mul(modulus.raw(), modulus.raw(), ell_v.raw());
    }

    fmpz_mod(result.raw(), result.raw(), c.field_unit_count.raw());
    out.swap(result);
    return true;
}

// Deterministic replacement for the reference `_primitive_element_mod_p`
// (ResidueRingMultGrp.jl:209): the first order element, in base-p digit
// order of its coordinates, whose residue has order q - 1.
bool find_field_generator(ResidueUnitPrimePower& c) noexcept {
    const Order& order = *c.prime.parent();
    const slong n = order.degree();
    Fmpz p;
    Fmpz bound;
    Fmpz k;
    Fmpz t;
    Fmpz digit;
    Fmpz cofactor;
    if (!c.prime.rational_prime(flint::FmpzRef(p))) {
        return false;
    }
    fmpz_pow_ui(bound.raw(), p.raw(), static_cast<ulong>(n));

    FmpzMat row(1, n);
    OrderElement candidate(order);
    ResidueFieldElement residue(c.field);
    ResidueFieldElement power(c.field);
    for (fmpz_one(k.raw()); fmpz_cmp(k.raw(), bound.raw()) < 0;
         fmpz_add_ui(k.raw(), k.raw(), 1)) {
        fmpz_set(t.raw(), k.raw());
        for (slong i = 0; i < n; ++i) {
            fmpz_fdiv_qr(t.raw(), digit.raw(), t.raw(), p.raw());
            fmpz_set(fmpz_mat_entry(row.raw(), 0, i), digit.raw());
        }
        if (!candidate.set_coordinates(flint::FmpzMatConstRef(row)) ||
            !residue.set_order_element(candidate)) {
            return false;
        }
        if (residue.is_zero()) {
            continue;
        }
        bool primitive = true;
        for (const Fmpz& ell : c.field_unit_primes) {
            fmpz_divexact(cofactor.raw(), c.field_unit_count.raw(),
                          ell.raw());
            if (!power.pow_fmpz(residue, flint::FmpzConstRef(cofactor))) {
                return false;
            }
            if (field_is_one(power)) {
                primitive = false;
                break;
            }
        }
        if (primitive) {
            c.field_generator.swap(candidate);
            c.field_generator_residue.swap(residue);
            return true;
        }
    }
    return false;
}

// One quadratic step: reference `_pu_mod_pv` and `_quadratic_method` for
// (1 + p^lower)/(1 + p^upper).  Appends the generators 1 + x_i to
// one_plus_generators.
bool build_layer(ResidueUnitLayer& layer,
                 ResidueUnitPrimePower& c,
                 const Ideal& prime_ideal,
                 slong lower,
                 slong upper) noexcept {
    const Order& order = *c.prime.parent();
    const slong n = order.degree();
    Ideal lower_power;
    Ideal upper_power;
    FmpzMat upper_hnf(0, 0);
    if (!ideal_power(lower_power, prime_ideal, lower) ||
        !ideal_power(upper_power, prime_ideal, upper) ||
        !ideal_hnf(layer.lower_hnf, lower_power) ||
        !ideal_hnf(upper_hnf, upper_power)) {
        return false;
    }
    layer.lower = lower;
    layer.upper = upper;

    // Basis of p^upper in the basis of p^lower (`basis_matrix(pv) *
    // basis_mat_inv(pu)`), as the relation matrix of p^lower/p^upper.
    FmpzMat relations(n, n);
    FmpzMat row(1, n);
    FmpzMat coords(0, 0);
    for (slong i = 0; i < n; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(row.raw(), 0, j),
                     fmpz_mat_entry(upper_hnf.raw(), i, j));
        }
        if (!hnf_coordinates(coords, layer.lower_hnf, row)) {
            return false;
        }
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(relations.raw(), i, j),
                     fmpz_mat_entry(coords.raw(), 0, j));
        }
    }
    if (!layer.quotient.set_relation_matrix(
                flint::FmpzMatConstRef(relations))) {
        return false;
    }

    const slong s = layer.quotient.invariant_count();
    layer.first_generator =
            static_cast<slong>(c.one_plus_generators.size());
    layer.generator_count = s;
    if (s == 0) {
        return true;
    }
    FmpzMat smith_generators(s, n);
    Fmpz largest;
    if (!layer.quotient.invariant_generator_matrix(
                flint::FmpzMatRef(smith_generators)) ||
        !layer.quotient.invariant(flint::FmpzRef(largest), s - 1)) {
        return false;
    }

    // gens[i] = sum_j (x_j mod S.snf[end]) b_j, then 1 + gens[i].
    Fmpz coefficient;
    for (slong i = 0; i < s; ++i) {
        FmpzMat generator_row(1, n);
        for (slong j = 0; j < n; ++j) {
            fmpz_mod(coefficient.raw(),
                     fmpz_mat_entry(smith_generators.raw(), i, j),
                     largest.raw());
            for (slong col = 0; col < n; ++col) {
                fmpz_addmul(fmpz_mat_entry(generator_row.raw(), 0, col),
                            coefficient.raw(),
                            fmpz_mat_entry(layer.lower_hnf.raw(), j, col));
            }
        }
        OrderElement generator(order);
        OrderElement one(order);
        if (!generator.set_coordinates(
                    flint::FmpzMatConstRef(generator_row)) ||
            !one.set_si(1) || !generator.add(generator, one) ||
            !reduce_element(generator, c.power_hnf)) {
            return false;
        }
        c.one_plus_generators.push_back(std::move(generator));
    }
    return true;
}

// Reference `_quadratic_method` discrete log for one layer: Smith
// coordinates of x - 1 in p^lower/p^upper.  Fails when x is not in
// 1 + p^lower.
bool layer_discrete_log(FmpzMat& out,
                        const ResidueUnitLayer& layer,
                        const OrderElement& element) noexcept {
    OrderElement one(*element.parent());
    OrderElement shifted(*element.parent());
    FmpzMat row(0, 0);
    FmpzMat coords(0, 0);
    if (!one.set_si(1) || !shifted.subtract(element, one) ||
        !element_row(row, shifted) ||
        !hnf_coordinates(coords, layer.lower_hnf, row)) {
        return false;
    }
    FmpzMat smith(1, layer.quotient.invariant_count());
    if (!layer.quotient.invariant_coordinates(
                flint::FmpzMatRef(smith), flint::FmpzMatConstRef(coords))) {
        return false;
    }
    out.swap(smith);
    return true;
}

// Reference `_iterative_method` discrete log: coordinates of an element of
// 1 + p on the concatenated layer generators.
bool iterative_discrete_log(FmpzMat& out,
                            const ResidueUnitPrimePower& c,
                            const OrderElement& element) noexcept {
    const Order& order = *c.prime.parent();
    const slong total = static_cast<slong>(c.one_plus_generators.size());
    FmpzMat result(1, total);
    OrderElement current(order);
    if (!current.set(element) || !reduce_element(current, c.power_hnf)) {
        return false;
    }
    Fmpz exponent;
    Fmpz inverse_exponent;
    fmpz_sub_ui(inverse_exponent.raw(), c.unit_count.raw(), 1);
    for (const ResidueUnitLayer& layer : c.layers) {
        FmpzMat digits(0, 0);
        if (!layer_discrete_log(digits, layer, current)) {
            return false;
        }
        OrderElement product;
        if (!set_one_mod(product, order, c.power_hnf)) {
            return false;
        }
        for (slong j = 0; j < layer.generator_count; ++j) {
            const slong index = layer.first_generator + j;
            fmpz_set(exponent.raw(), fmpz_mat_entry(digits.raw(), 0, j));
            fmpz_set(fmpz_mat_entry(result.raw(), 0, index), exponent.raw());
            if (fmpz_is_zero(exponent.raw()) != 0) {
                continue;
            }
            OrderElement power;
            if (!pow_mod(power,
                         c.one_plus_generators[static_cast<std::size_t>(
                                 index)],
                         exponent, c.power_hnf) ||
                !mul_mod(product, product, power, c.power_hnf)) {
                return false;
            }
        }
        // b1 = divexact(b1, prod)
        OrderElement inverse;
        if (!pow_mod(inverse, product, inverse_exponent, c.power_hnf) ||
            !mul_mod(current, current, inverse, c.power_hnf)) {
            return false;
        }
    }
    OrderElement one;
    if (!set_one_mod(one, order, c.power_hnf) || !current.equal(one)) {
        return false;
    }
    out.swap(result);
    return true;
}

// Reference `_iterative_method` with `base_method = :quadratic`, followed by
// the Smith normal form of `_1_plus_p_mod_1_plus_pv`.
bool build_one_plus(ResidueUnitPrimePower& c,
                    const Ideal& prime_ideal) noexcept {
    const Order& order = *c.prime.parent();
    const slong v = c.exponent;
    FmpzMat relations(0, 0);
    slong lower = 1;
    while (lower != v) {
        const slong upper = lower * 2 < v ? lower * 2 : v;
        ResidueUnitLayer layer;
        if (!build_layer(layer, c, prime_ideal, lower, upper)) {
            return false;
        }
        const slong old_count = fmpz_mat_nrows(relations.raw());
        const slong added = layer.generator_count;
        FmpzMat expanded(old_count + added, old_count + added);
        flint::FmpzVec layer_invariants(added);
        if (!layer.quotient.invariants(flint::FmpzVecRef(
                    layer_invariants.data(), layer_invariants.length()))) {
            return false;
        }
        // `_expand`: copy the upper triangle of M, put the new invariants on
        // the diagonal, and write -disc_log(prod g_j^M[i, j]) beside row i.
        for (slong i = 0; i < old_count; ++i) {
            for (slong j = i; j < old_count; ++j) {
                fmpz_set(fmpz_mat_entry(expanded.raw(), i, j),
                         fmpz_mat_entry(relations.raw(), i, j));
            }
        }
        for (slong i = 0; i < added; ++i) {
            fmpz_set(fmpz_mat_entry(expanded.raw(), old_count + i,
                                    old_count + i),
                     layer_invariants.data() + i);
        }
        for (slong i = 0; i < old_count && added > 0; ++i) {
            OrderElement product;
            if (!set_one_mod(product, order, c.power_hnf)) {
                return false;
            }
            for (slong j = 0; j < old_count; ++j) {
                const fmpz* entry = fmpz_mat_entry(relations.raw(), i, j);
                if (fmpz_is_zero(entry) != 0) {
                    continue;
                }
                Fmpz exponent;
                fmpz_set(exponent.raw(), entry);
                OrderElement power;
                if (!pow_mod_unit(power,
                                  c.one_plus_generators[
                                          static_cast<std::size_t>(j)],
                                  exponent, c.unit_count, c.power_hnf) ||
                    !mul_mod(product, product, power, c.power_hnf)) {
                    return false;
                }
            }
            FmpzMat alpha(0, 0);
            if (!layer_discrete_log(alpha, layer, product)) {
                return false;
            }
            for (slong j = 0; j < added; ++j) {
                fmpz_neg(fmpz_mat_entry(expanded.raw(), i, old_count + j),
                         fmpz_mat_entry(alpha.raw(), 0, j));
            }
        }
        if (added > 0) {
            relations.swap(expanded);
        }
        c.layers.push_back(std::move(layer));
        lower = upper;
    }

    if (!c.one_plus_group.set_relation_matrix(
                flint::FmpzMatConstRef(relations))) {
        return false;
    }
    return true;
}

// Reference `_multgrp_mod_pv` (ResidueRingMultGrp.jl:116).
bool build_prime_power(ResidueUnitPrimePower& out,
                       const PrimeIdeal& prime,
                       slong exponent) noexcept {
    const Order* order = prime.parent();
    if (order == nullptr || !prime.has_prime_data() || exponent < 1 ||
        !order->is_maximal()) {
        return false;
    }

    ResidueUnitPrimePower c;
    c.prime = PrimeIdeal(*order);
    Ideal prime_ideal(*order);
    if (!c.prime.set(prime) || !prime.get_ideal(prime_ideal) ||
        !ideal_power(c.power, prime_ideal, exponent) ||
        !ideal_hnf(c.power_hnf, c.power) || !c.field.set_prime(prime)) {
        return false;
    }
    c.exponent = exponent;

    Fmpz q;
    if (!prime.norm(flint::FmpzRef(q))) {
        return false;
    }
    fmpz_sub_ui(c.field_unit_count.raw(), q.raw(), 1);
    fmpz_pow_ui(c.unit_count.raw(), q.raw(),
                static_cast<ulong>(exponent - 1));
    fmpz_mul(c.unit_count.raw(), c.unit_count.raw(),
             c.field_unit_count.raw());
    {
        flint::FmpzFactor factorization;
        fmpz_factor(factorization.raw(), c.field_unit_count.raw());
        for (slong i = 0; i < factorization.raw()->num; ++i) {
            Fmpz ell;
            fmpz_set(ell.raw(), factorization.raw()->p + i);
            c.field_unit_primes.push_back(std::move(ell));
            c.field_unit_exponents.push_back(
                    static_cast<slong>(factorization.raw()->exp[i]));
        }
    }

    c.field_generator = OrderElement(*order);
    c.field_generator_residue = ResidueFieldElement(c.field);
    if (!find_field_generator(c)) {
        return false;
    }

    OrderElement tame(*order);
    if (!tame.set(c.field_generator) || !reduce_element(tame, c.power_hnf)) {
        return false;
    }

    if (exponent == 1) {
        fmpz_one(c.one_plus_exponent.raw());
        fmpz_one(c.exponent_inverse.raw());
        Fmpz relation;
        fmpz_set(relation.raw(), c.field_unit_count.raw());
        c.generators.push_back(std::move(tame));
        c.relations.push_back(std::move(relation));
        out = std::move(c);
        return true;
    }

    if (!build_one_plus(c, prime_ideal)) {
        return false;
    }

    // `_1_plus_p_mod_1_plus_pv`: Smith generators of G2,
    // prod g_j^(x_j mod S.snf[end]).
    const slong t = c.one_plus_group.invariant_count();
    const slong total = static_cast<slong>(c.one_plus_generators.size());
    fmpz_one(c.one_plus_exponent.raw());
    FmpzMat smith_generators(t, total);
    if (t > 0 &&
        (!c.one_plus_group.invariant(flint::FmpzRef(c.one_plus_exponent),
                                     t - 1) ||
         !c.one_plus_group.invariant_generator_matrix(
                 flint::FmpzMatRef(smith_generators)))) {
        return false;
    }

    // gen1_obcs = gen1^(G2.snf[end]); obcs_inv = G2.snf[end]^-1 mod (q - 1).
    OrderElement tame_power;
    if (!pow_mod(tame_power, tame, c.one_plus_exponent, c.power_hnf)) {
        return false;
    }
    if (fmpz_is_one(c.field_unit_count.raw()) != 0) {
        fmpz_zero(c.exponent_inverse.raw());
    } else if (fmpz_invmod(c.exponent_inverse.raw(),
                           c.one_plus_exponent.raw(),
                           c.field_unit_count.raw()) == 0) {
        return false;
    }
    Fmpz relation;
    fmpz_set(relation.raw(), c.field_unit_count.raw());
    c.generators.push_back(std::move(tame_power));
    c.relations.push_back(std::move(relation));

    Fmpz exponent_value;
    for (slong i = 0; i < t; ++i) {
        OrderElement generator;
        if (!set_one_mod(generator, *order, c.power_hnf)) {
            return false;
        }
        for (slong j = 0; j < total; ++j) {
            fmpz_mod(exponent_value.raw(),
                     fmpz_mat_entry(smith_generators.raw(), i, j),
                     c.one_plus_exponent.raw());
            if (fmpz_is_zero(exponent_value.raw()) != 0) {
                continue;
            }
            OrderElement power;
            if (!pow_mod(power,
                         c.one_plus_generators[static_cast<std::size_t>(j)],
                         exponent_value, c.power_hnf) ||
                !mul_mod(generator, generator, power, c.power_hnf)) {
                return false;
            }
        }
        Fmpz invariant;
        if (!c.one_plus_group.invariant(flint::FmpzRef(invariant), i)) {
            return false;
        }
        c.generators.push_back(std::move(generator));
        c.relations.push_back(std::move(invariant));
    }

    out = std::move(c);
    return true;
}

// Reference `disc_log2` in `_multgrp_mod_pv`: coordinates on
// [gen1^(exp G2), Smith generators of G2].
bool prime_power_discrete_log(FmpzMat& out,
                              const ResidueUnitPrimePower& c,
                              const OrderElement& element) noexcept {
    const slong count = static_cast<slong>(c.generators.size());
    ResidueFieldElement residue(c.field);
    Fmpz log;
    if (!residue.set_order_element(element) || residue.is_zero() ||
        !field_discrete_log(log, residue, c)) {
        return false;
    }

    FmpzMat result(1, count);
    Fmpz r;
    if (fmpz_is_one(c.field_unit_count.raw()) != 0) {
        fmpz_zero(r.raw());
    } else {
        fmpz_mul(r.raw(), log.raw(), c.exponent_inverse.raw());
        fmpz_mod(r.raw(), r.raw(), c.field_unit_count.raw());
    }
    fmpz_set(fmpz_mat_entry(result.raw(), 0, 0), r.raw());
    if (c.exponent == 1) {
        out.swap(result);
        return true;
    }

    // y = x * gen1_obcs^(-r mod rel1)
    Fmpz negated;
    fmpz_neg(negated.raw(), r.raw());
    if (fmpz_is_one(c.field_unit_count.raw()) == 0) {
        fmpz_mod(negated.raw(), negated.raw(), c.field_unit_count.raw());
    } else {
        fmpz_zero(negated.raw());
    }
    OrderElement correction;
    OrderElement reduced;
    if (!pow_mod(correction, c.generators[0], negated, c.power_hnf) ||
        !mul_mod(reduced, element, correction, c.power_hnf)) {
        return false;
    }

    FmpzMat iterative(0, 0);
    if (!iterative_discrete_log(iterative, c, reduced)) {
        return false;
    }
    const slong t = c.one_plus_group.invariant_count();
    FmpzMat smith(1, t);
    if (!c.one_plus_group.invariant_coordinates(
                flint::FmpzMatRef(smith), flint::FmpzMatConstRef(iterative))) {
        return false;
    }
    for (slong i = 0; i < t; ++i) {
        fmpz_set(fmpz_mat_entry(result.raw(), 0, 1 + i),
                 fmpz_mat_entry(smith.raw(), 0, i));
    }
    out.swap(result);
    return true;
}

}  // namespace

ResidueUnitGroup::~ResidueUnitGroup() noexcept = default;

ResidueUnitGroup::ResidueUnitGroup(ResidueUnitGroup&& other) noexcept {
    swap(other);
}

ResidueUnitGroup& ResidueUnitGroup::operator=(
        ResidueUnitGroup&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void ResidueUnitGroup::swap(ResidueUnitGroup& other) noexcept {
    modulus_.swap(other.modulus_);
    modulus_hnf_.swap(other.modulus_hnf_);
    components_.swap(other.components_);
    group_.swap(other.group_);
    invariants_.swap(other.invariants_);
    generators_.swap(other.generators_);
    std::swap(defined_, other.defined_);
}

void ResidueUnitGroup::clear() noexcept {
    ResidueUnitGroup empty;
    swap(empty);
}

bool ResidueUnitGroup::set_modulus(const Ideal& modulus) noexcept {
    const Order* order = modulus.parent();
    if (order == nullptr || !modulus.has_hnf() || !order->is_maximal()) {
        return false;
    }

    IdealFactorization factorization(*order);
    if (!factorization.factor(modulus)) {
        return false;
    }
    std::vector<ResidueUnitPrimePower> components;
    for (slong i = 0; i < factorization.length(); ++i) {
        PrimeIdeal prime(*order);
        slong exponent = 0;
        ResidueUnitPrimePower component;
        if (!factorization.prime(prime, i) ||
            !factorization.exponent(exponent, i) ||
            !build_prime_power(component, prime, exponent)) {
            return false;
        }
        components.push_back(std::move(component));
    }
    return finish(components, modulus);
}

bool ResidueUnitGroup::set_prime_power(const PrimeIdeal& prime,
                                       slong exponent) noexcept {
    std::vector<ResidueUnitPrimePower> components;
    ResidueUnitPrimePower component;
    if (!build_prime_power(component, prime, exponent)) {
        return false;
    }
    Ideal modulus(*prime.parent());
    if (!modulus.set(component.power)) {
        return false;
    }
    components.push_back(std::move(component));
    return finish(components, modulus);
}

// Reference `_direct_product` (CRT lifts of the component generators) and
// `snf(G, GtoR)` (Smith generators prod G_j^(x_j mod S.snf[end])).
bool ResidueUnitGroup::finish(std::vector<ResidueUnitPrimePower>& components,
                              const Ideal& modulus) noexcept {
    const Order& order = *modulus.parent();
    ResidueUnitGroup candidate;
    candidate.modulus_ = Ideal(order);
    if (!candidate.modulus_.set(modulus) ||
        !ideal_hnf(candidate.modulus_hnf_, modulus)) {
        return false;
    }

    slong total = 0;
    for (const ResidueUnitPrimePower& c : components) {
        total += static_cast<slong>(c.generators.size());
    }

    std::vector<OrderElement> lifted;
    FmpzMat relations(total, total);
    slong offset = 0;
    for (std::size_t i = 0; i < components.size(); ++i) {
        const ResidueUnitPrimePower& c = components[i];
        OrderElement in_power(order);
        OrderElement in_rest(order);
        if (components.size() > 1) {
            // crt(g, p^v, 1, prod of the other prime powers)
            Ideal rest(order);
            if (!rest.one()) {
                return false;
            }
            for (std::size_t j = 0; j < components.size(); ++j) {
                if (j == i) {
                    continue;
                }
                Ideal next(order);
                if (!next.multiply(rest, components[j].power)) {
                    return false;
                }
                rest.swap(next);
            }
            if (!c.power.add_to_one(in_power, in_rest, rest)) {
                return false;
            }
        }
        for (std::size_t j = 0; j < c.generators.size(); ++j) {
            OrderElement g(order);
            if (!g.set(c.generators[j])) {
                return false;
            }
            if (components.size() > 1) {
                OrderElement scaled(order);
                if (!scaled.multiply(g, in_rest) ||
                    !g.add(scaled, in_power)) {
                    return false;
                }
            }
            if (!reduce_element(g, candidate.modulus_hnf_)) {
                return false;
            }
            lifted.push_back(std::move(g));
            fmpz_set(fmpz_mat_entry(relations.raw(), offset, offset),
                     c.relations[j].raw());
            ++offset;
        }
    }

    if (total > 0) {
        if (!candidate.group_.set_relation_matrix(
                    flint::FmpzMatConstRef(relations))) {
            return false;
        }
        const slong r = candidate.group_.invariant_count();
        candidate.invariants_ = flint::FmpzVec(r);
        FmpzMat smith_generators(r, total);
        if (!candidate.group_.invariants(flint::FmpzVecRef(
                    candidate.invariants_.data(), r)) ||
            !candidate.group_.invariant_generator_matrix(
                    flint::FmpzMatRef(smith_generators))) {
            return false;
        }
        Fmpz largest;
        Fmpz exponent;
        if (r > 0) {
            fmpz_set(largest.raw(), candidate.invariants_.data() + r - 1);
        }
        for (slong i = 0; i < r; ++i) {
            OrderElement generator;
            if (!set_one_mod(generator, order, candidate.modulus_hnf_)) {
                return false;
            }
            for (slong j = 0; j < total; ++j) {
                fmpz_mod(exponent.raw(),
                         fmpz_mat_entry(smith_generators.raw(), i, j),
                         largest.raw());
                if (fmpz_is_zero(exponent.raw()) != 0) {
                    continue;
                }
                OrderElement power;
                if (!pow_mod(power, lifted[static_cast<std::size_t>(j)],
                             exponent, candidate.modulus_hnf_) ||
                    !mul_mod(generator, generator, power,
                             candidate.modulus_hnf_)) {
                    return false;
                }
            }
            candidate.generators_.push_back(std::move(generator));
        }
    }

    candidate.components_.swap(components);
    candidate.defined_ = true;
    swap(candidate);
    return true;
}

bool ResidueUnitGroup::is_defined() const noexcept {
    return defined_;
}

const Order* ResidueUnitGroup::parent_order() const noexcept {
    return defined_ ? modulus_.parent() : nullptr;
}

bool ResidueUnitGroup::get_modulus(Ideal& out) const noexcept {
    return defined_ && same_order_parent(out.parent(), modulus_.parent()) &&
           out.set(modulus_);
}

slong ResidueUnitGroup::generator_count() const noexcept {
    return defined_ ? static_cast<slong>(generators_.size()) : 0;
}

bool ResidueUnitGroup::invariants(flint::FmpzVecRef out) const noexcept {
    if (!defined_ || out.length() != generator_count()) {
        return false;
    }
    for (slong i = 0; i < generator_count(); ++i) {
        fmpz_set(out.data() + i, invariants_.data() + i);
    }
    return true;
}

bool ResidueUnitGroup::order(flint::FmpzRef out) const noexcept {
    if (!defined_) {
        return false;
    }
    Fmpz product;
    fmpz_one(product.raw());
    for (slong i = 0; i < generator_count(); ++i) {
        fmpz_mul(product.raw(), product.raw(), invariants_.data() + i);
    }
    fmpz_set(out.raw(), product.raw());
    return true;
}

bool ResidueUnitGroup::generator(OrderElement& out,
                                 slong index) const noexcept {
    if (!defined_ || index < 0 || index >= generator_count() ||
        !same_order_parent(out.parent(), modulus_.parent())) {
        return false;
    }
    return out.set(generators_[static_cast<std::size_t>(index)]);
}

bool ResidueUnitGroup::discrete_log(flint::FmpzMatRef out,
                                    const OrderElement& element) const noexcept {
    if (!defined_ || fmpz_mat_nrows(out.raw()) != 1 ||
        fmpz_mat_ncols(out.raw()) != generator_count() ||
        !same_order_parent(element.parent(), modulus_.parent())) {
        return false;
    }

    slong total = 0;
    for (const ResidueUnitPrimePower& c : components_) {
        total += static_cast<slong>(c.generators.size());
    }
    // Every component must see a unit, including components whose group is
    // trivial, so the test runs before the shortcut for an empty group.
    FmpzMat row(1, total);
    slong offset = 0;
    for (const ResidueUnitPrimePower& c : components_) {
        FmpzMat part(0, 0);
        if (!prime_power_discrete_log(part, c, element)) {
            return false;
        }
        for (std::size_t j = 0; j < c.generators.size(); ++j) {
            fmpz_set(fmpz_mat_entry(row.raw(), 0, offset),
                     fmpz_mat_entry(part.raw(), 0, static_cast<slong>(j)));
            ++offset;
        }
    }

    FmpzMat result(1, generator_count());
    if (total > 0 &&
        !group_.invariant_coordinates(flint::FmpzMatRef(result),
                                      flint::FmpzMatConstRef(row))) {
        return false;
    }
    fmpz_mat_set(out.raw(), result.raw());
    return true;
}

bool ResidueUnitGroup::exp(OrderElement& out,
                           flint::FmpzMatConstRef exponents) const noexcept {
    if (!defined_ || fmpz_mat_nrows(exponents.raw()) != 1 ||
        fmpz_mat_ncols(exponents.raw()) != generator_count() ||
        !same_order_parent(out.parent(), modulus_.parent())) {
        return false;
    }
    const Order& order = *modulus_.parent();
    OrderElement result;
    if (!set_one_mod(result, order, modulus_hnf_)) {
        return false;
    }
    Fmpz exponent;
    for (slong i = 0; i < generator_count(); ++i) {
        fmpz_mod(exponent.raw(), fmpz_mat_entry(exponents.raw(), 0, i),
                 invariants_.data() + i);
        OrderElement power;
        if (!pow_mod(power, generators_[static_cast<std::size_t>(i)],
                     exponent, modulus_hnf_) ||
            !mul_mod(result, result, power, modulus_hnf_)) {
            return false;
        }
    }
    return out.set(result);
}

}  // namespace silex::detail
