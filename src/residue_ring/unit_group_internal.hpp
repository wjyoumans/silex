#pragma once

#include <flint/flint.h>

#include <silex/abelian_group.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/ideal.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_field.hpp>

#include <vector>

namespace silex::detail {

// Unit group (O/m)^x of a residue ring of a maximal order O, with generators,
// its cyclic structure and a discrete logarithm.  Internal routine: no public
// header exposes it.
//
// Port of reference v0.38.6 src/NumFieldOrd/NfOrd/ResidueRingMultGrp.jl:
// `_multgrp` (line 64), `_multgrp_mod_pv` (line 116), `_multgrp_mod_p`
// (line 171), `_1_plus_p_mod_1_plus_pv` (line 236), `_iterative_method`
// (lines 263-330) with `base_method = :quadratic, use_p_adic = false`,
// `_expand` (line 332), `_pu_mod_pv` (line 374), `_quadratic_method`
// (line 437), `pohlig_hellman` (line 754), `_direct_product` (line 1357) and
// `snf(G, GtoR)` (line 1248).  unit_group.cpp records the deviations.
//
// Group conventions (those of the reference's `snf`): the published group is
// Z/d_1 x ... x Z/d_r with 1 < d_1 | d_2 | ... | d_r, generator(i) has order
// d_i in (O/m)^x, and discrete_log returns exponents e with 0 <= e_i < d_i
// and x = prod generator(i)^e_i modulo m.

// One step (1 + p^lower) / (1 + p^upper), lower < upper <= 2 * lower, of the
// iterative method; x -> 1 + x identifies it with p^lower / p^upper.
struct ResidueUnitLayer {
    slong lower = 0;
    slong upper = 0;
    // HNF basis of p^lower, rows in order coordinates.
    flint::FmpzMat lower_hnf{0, 0};
    // p^lower / p^upper presented on the rows of lower_hnf.
    FiniteAbelianGroup quotient;
    // This layer's generators are one_plus_generators[first, first + count).
    slong first_generator = 0;
    slong generator_count = 0;
};

// (O/p^v)^x as G1 x G2 (`_multgrp_mod_pv`): G1 = (O/p)^x, cyclic of order
// q - 1 on the generator gen1^(exp G2), and G2 = (1 + p)/(1 + p^v) in Smith
// normal form.  Its relations are the diagonal entries `relations`.
struct ResidueUnitPrimePower {
    PrimeIdeal prime;
    slong exponent = 0;
    Ideal power;
    flint::FmpzMat power_hnf{0, 0};
    // #(O/p^v)^x = (q - 1) q^(v - 1); exponents of units reduce modulo it.
    flint::Fmpz unit_count;
    ResidueField field;
    flint::Fmpz field_unit_count;
    std::vector<flint::Fmpz> field_unit_primes;
    std::vector<slong> field_unit_exponents;
    OrderElement field_generator;
    ResidueFieldElement field_generator_residue;
    std::vector<ResidueUnitLayer> layers;
    std::vector<OrderElement> one_plus_generators;
    FiniteAbelianGroup one_plus_group;
    flint::Fmpz one_plus_exponent;
    flint::Fmpz exponent_inverse;
    std::vector<OrderElement> generators;
    std::vector<flint::Fmpz> relations;
};

class ResidueUnitGroup {
public:
    ResidueUnitGroup() noexcept = default;
    ~ResidueUnitGroup() noexcept;

    ResidueUnitGroup(const ResidueUnitGroup&) = delete;
    ResidueUnitGroup& operator=(const ResidueUnitGroup&) = delete;

    ResidueUnitGroup(ResidueUnitGroup&& other) noexcept;
    ResidueUnitGroup& operator=(ResidueUnitGroup&& other) noexcept;

    void swap(ResidueUnitGroup& other) noexcept;
    void clear() noexcept;

    // (O/m)^x for a nonzero ideal m of a maximal order, factored into prime
    // powers and joined by CRT (`_multgrp`).  On failure the object is
    // unchanged.
    bool set_modulus(const Ideal& modulus) noexcept;
    // (O/p^v)^x for a prime p of a maximal order and v >= 1
    // (`_multgrp_mod_pv`, published in Smith normal form as `_multgrp` does).
    // On failure the object is unchanged.
    bool set_prime_power(const PrimeIdeal& prime, slong exponent) noexcept;

    bool is_defined() const noexcept;
    const Order* parent_order() const noexcept;
    bool get_modulus(Ideal& out) const noexcept;

    slong generator_count() const noexcept;
    bool invariants(flint::FmpzVecRef out) const noexcept;
    bool order(flint::FmpzRef out) const noexcept;
    bool generator(OrderElement& out, slong index) const noexcept;

    // Fails, leaving `out` unchanged, when x is not a unit modulo m.
    // `out` must be 1 x generator_count().
    bool discrete_log(flint::FmpzMatRef out,
                      const OrderElement& element) const noexcept;
    // prod generator(i)^exponents_i reduced modulo m; `exponents` is
    // 1 x generator_count() and may hold any integers.
    bool exp(OrderElement& out,
             flint::FmpzMatConstRef exponents) const noexcept;

private:
    bool finish(std::vector<ResidueUnitPrimePower>& components,
                const Ideal& modulus) noexcept;

    Ideal modulus_;
    flint::FmpzMat modulus_hnf_{0, 0};
    std::vector<ResidueUnitPrimePower> components_;
    FiniteAbelianGroup group_;
    flint::FmpzVec invariants_{0};
    std::vector<OrderElement> generators_;
    bool defined_ = false;
};

inline void swap(ResidueUnitGroup& left, ResidueUnitGroup& right) noexcept {
    left.swap(right);
}

}  // namespace silex::detail
