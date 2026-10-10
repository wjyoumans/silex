#pragma once

#include "class_group_storage_internal.hpp"
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/nmod_poly.hpp>

#include <vector>

namespace silex::detail::relation_search {

struct OrderCoordinateElementConversion {
    const Order* order = nullptr;
    flint::FmpqMat basis{0, 0};
    flint::FmpqMat power_row{0, 0};
    flint::FmpqPoly polynomial;

    bool reset(const Order& next_order) noexcept;
    bool set(Element& out,
             flint::FmpzMatConstRef coordinates,
             flint::FmpzConstRef den) noexcept;
};

// Exact Res(f, a) for one fixed monic integral f of degree n and integral
// a of degree < n, computed multimodularly.  Modulo a prime p at which f
// splits into n distinct linear factors x - r_i, Res(f, a) = prod_i a(r_i)
// (mod p) because f is monic, so each residue costs n dot products of
// length n against the stored powers r_i^j.  The integer is recovered by
// symmetric CRT once the product of the primes used exceeds twice the
// Hadamard bound |Res(f, a)| <= ||f||_2^deg(a) ||a||_2^n, the bound FLINT's
// _fmpz_poly_resultant_modular uses (FLINT src/fmpz_poly/resultant_modular.c,
// citing Theorem 7 of doi:10.1016/j.laa.2009.08.012).  Split primes are searched lazily;
// when they do not cover a candidate's bound, the caller computes
// fmpz_poly_resultant instead, so the value is always the exact resultant.
struct SplitPrimeResultant {
    flint::FmpzPoly field_polynomial;
    slong degree = 0;
    flint_bitcnt_t field_norm_squared_bits = 0;
    ulong last_tested_prime = 0;
    slong tested_primes = 0;
    bool enabled = false;
    bool search_exhausted = false;
    std::vector<nmod_t> moduli;
    std::vector<ulong> root_powers;
    std::vector<ulong> residues;
    flint::Fmpz coordinate_norm_squared;
    flint::Fmpz modulus;
    flint::Fmpz crt_scratch;

    // Enables the evaluator for f when deg(f) is large enough for the
    // search to pay off; otherwise leaves it disabled.  A reset with the
    // same f while enabled keeps the split primes found so far.
    void reset(const fmpz_poly_struct* monic_field_polynomial) noexcept;

    // Sets out = Res(f, a) and handled = true when the split primes cover
    // the bound for a; handled = false leaves out unspecified.
    void resultant(bool& handled,
                   fmpz_t out,
                   const fmpz_poly_struct* coordinate_polynomial) noexcept;

private:
    bool add_next_split_prime() noexcept;
};

struct NormPrefilter {
    flint::Fmpz rational_prime_product;
    flint::Fmpz large_prime_norm_bound;
    flint::Fmpq norm_scratch;
    flint::Fmpz norm_resultant_scratch;
    flint::Fmpz norm_remainder;
    flint::Fmpz norm_gcd;
    flint::FmpzPoly field_polynomial;
    flint::FmpzPoly coordinate_polynomial;
    SplitPrimeResultant split_prime_resultant;
    slong rejected_count = 0;
    bool valid = false;
    bool allow_large_prime = false;
    bool has_equation_order_norm_polynomial = false;
};

void coordinates_from_lattice_combination(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef coefficients,
        flint::FmpzMatConstRef basis) noexcept;

bool build_norm_prefilter(NormPrefilter& out,
                          const FactorBase& base,
                          flint::FmpzConstRef factor_base_bound,
                          bool allow_large_prime) noexcept;

bool reduced_basis_first_row_is_scalar_rational(
        bool& scalar,
        const Order& order,
        flint::FmpzMatConstRef basis) noexcept;

bool element_is_scalar_rational(bool& scalar,
                                const Element& element) noexcept;

bool try_coordinate_candidate_den(
        ClassGroupContext& context,
        const Order& order,
        flint::FmpzMat& coordinates,
        flint::FmpzConstRef den,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        bool& partial_throttle_exit,
        detail::RelationAdmissionCache* admission_cache = nullptr,
        bool random_relation = false,
        slong* factor_attempts = nullptr,
        slong max_factor_attempts = 0,
        bool* factor_attempt_limit_reached = nullptr,
        Relation* admission_scratch = nullptr,
        OrderCoordinateElementConversion* coordinate_conversion = nullptr,
        bool admission_scratch_base_verified = false) noexcept;

bool try_coordinate_candidate(
        ClassGroupContext& context,
        const Order& order,
        flint::FmpzMat& coordinates,
        NormPrefilter* norm_prefilter,
        ClassGroupRelationSource source,
        slong target_relation_kernel_units,
        slong max_candidates,
        slong max_relations,
        slong& candidates_tried,
        slong& accepted_relations,
        CertificationMode requested_certification,
        bool& goal_reached,
        bool& partial_throttle_exit,
        OrderCoordinateElementConversion* coordinate_conversion = nullptr)
        noexcept;

}  // namespace silex::detail::relation_search
