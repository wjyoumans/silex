#pragma once

#include "class_group_storage_internal.hpp"
#include <silex/flint/arf.hpp>
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
// citing Theorem 7 of doi:10.1016/j.laa.2009.08.012).  When it is smaller,
// the bound used is instead the T2 bound
// |Res(f, a)|^2 = prod_i |a(r_i)|^2 <= (T2(a) / n)^n <= (lambda ||a||_2^2 / n)^n,
// the inequality of arithmetic and geometric means on the n values
// |a(r_i)|^2, with T2(a) = sum_i |a(r_i)|^2 the quadratic form that Cohen
// calls "size" (GTM 138, Section 4.4.2, Definition 4.4.10 and the form
// size(M_x), p. 169), followed by the
// Rayleigh-quotient bound x^T G x <= lambda_max(G) x^T x for the real
// symmetric T2 Gram matrix G of the power basis (Horn and Johnson, Matrix
// Analysis, Theorem 4.2.2), where lambda >= lambda_max(G) is certified once
// per f: lambda I - G is certainly positive definite by FLINT arb_mat_cho
// on ball enclosures of G built from arb_fmpz_poly_complex_roots, or else
// lambda = trace(G), its upper end, since G is positive semidefinite.  It
// is not computed when the exact power sums p_k of f show
// lambda_max(G) >= max |p_(2j)| >= n ||f||_2^2, where the real-valued T2
// bound is never below the real-valued Hadamard bound.  The rounded bit
// counts the CRT compares can differ by a few bits, so the skip may pass
// up a slightly smaller modulus; that affects speed only.  The bound is
// evaluated with mag_t upper-bound arithmetic, so every rounding is
// outward.  Split primes are searched lazily;
// when they do not cover a candidate's bound, the caller computes
// fmpz_poly_resultant instead, so the value is always the exact resultant.
struct SplitPrimeResultant {
    enum class T2BoundState { not_computed, available, unavailable };

    flint::FmpzPoly field_polynomial;
    slong degree = 0;
    flint_bitcnt_t field_norm_squared_bits = 0;
    ulong last_tested_prime = 0;
    slong tested_primes = 0;
    // Test-only diagnostic: split primes used by the last handled call.
    // A call that is not handled leaves it unchanged.
    slong last_prime_count = 0;
    bool enabled = false;
    bool search_exhausted = false;
    // Certified upper bound for lambda_max(G) / n, computed on the first
    // call after a reset with a new f.
    T2BoundState t2_bound_state = T2BoundState::not_computed;
    flint::Arf t2_ratio_bound;
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
