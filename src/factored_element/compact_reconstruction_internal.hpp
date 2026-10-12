#pragma once

#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/nmod_mat.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/order.hpp>

#include <cstddef>
#include <span>
#include <vector>

namespace silex::detail {

enum class BoundedCompactReconstructionStatus {
    success,
    invalid_input,
    exponent_overflow,
    coordinate_conversion_failed,
    prime_exhausted,
    certified_bound_violated,
    output_conversion_failed,
};

struct BoundedCompactReconstructionOptions {
    // reference init_modular_big starts at the high half of the machine-word
    // primes. Tests may lower this private setting to exercise bad primes.
    ulong prime_search_start = UWORD(1) << (FLINT_BITS - 1);
};

struct BoundedCompactReconstructionReport {
    flint::Fmpz centered_crt_modulus;
    BoundedCompactReconstructionStatus status =
            BoundedCompactReconstructionStatus::invalid_input;
    std::size_t factor_rows = 0;
    std::size_t prime_trials = 0;
    std::size_t primes_used = 0;
    std::size_t field_polynomial_prime_rejections = 0;
    std::size_t order_basis_prime_rejections = 0;
    std::size_t denominator_prime_rejections = 0;
    std::size_t noninvertible_prime_rejections = 0;
    flint_bitcnt_t modulus_bits = 0;
    ulong last_prime = 0;
};

const char* bounded_compact_reconstruction_status_name(
        BoundedCompactReconstructionStatus status) noexcept;

// Source trace: compact-reconstruction coordinate-bound algorithm.
// `vec_chinese_units` -> `chinese_unit` -> `FlxqX_chinese_unit`.
//
// The caller must certify that every input evaluates to an element of the order
// and that the absolute value of every resulting order-basis coordinate is at
// most certified_coordinate_bound. The reconstruction is exact once the
// centered CRT modulus is strictly greater than twice that bound. This
// internal checkpoint deliberately does not produce the bound.
bool bounded_compact_reconstruct(
        BoundedCompactReconstructionReport& report, std::vector<Element>& out,
        const Order& order, std::span<const FactoredElement> inputs,
        flint::FmpzConstRef certified_coordinate_bound,
        const BoundedCompactReconstructionOptions& options = {}) noexcept;

// One factor of a batch of modular power products: the base reduced modulo
// the modulus, its inverse modulo the modulus (required when some exponent is
// negative, otherwise may be null), and one exponent per output.
struct ModularPowerProductFactor {
    const flint::NmodPoly* base = nullptr;
    const flint::NmodPoly* inverse = nullptr;
    const std::vector<flint::Fmpz>* exponents = nullptr;
};

// Row i of power_residues receives the coefficients, in the power basis, of
// prod_k base_k^exponents_k[i] in (Z/pZ)[x]/(modulus), where modulus is monic
// of degree power_residues.ncols.  The product is an element of that ring, so
// the result does not depend on how it is evaluated.  Each base's odd powers
// are tabulated once and shared by every output, and the exponentiations of
// one output are interleaved so that they share their squarings: interleaved
// sliding-window exponentiation (Moller, "Algorithms for multi-exponentiation",
// SAC 2001, Section 3.1, the basic interleaving exponentiation method;
// sliding windows as in Menezes, van Oorschot and Vanstone, Handbook of
// Applied Cryptography, Algorithm 14.85 and Table 14.16).  Before
// tabulating, an exponent whose magnitude exceeds the exponent lambda of the
// unit group of the ring (modular_unit_group_exponent) is replaced by
// ((e - 1) mod lambda) + 1, which gives the same power of every element when
// e > 0 and of every unit when e < 0; when lambda is not available the
// exponents are used as given.
void modular_power_products(
        flint::NmodMat& power_residues,
        std::span<const ModularPowerProductFactor> factors,
        const flint::NmodPoly& modulus) noexcept;

// Sets out to the exponent lambda of the unit group of (Z/pZ)[x]/(modulus)
// and returns true when p = modulus.mod.n is prime and modulus is monic,
// nonconstant and squarefree modulo p.  The ring is then the product of the
// fields F_(p^d) over the irreducible factors of degree d, and lambda is the
// least common multiple of the p^d - 1.  Every element a of the ring
// satisfies a^(lambda + 1) = a, so a^e = a^(((e - 1) mod lambda) + 1) for
// e > 0.  Returns false, leaving out unspecified, otherwise.
bool modular_unit_group_exponent(flint::Fmpz& out,
                                 const flint::NmodPoly& modulus) noexcept;

struct WideFactoredUnits;

// The same reconstruction for a batch with arbitrary-size exponents: output
// i is prod_k inputs.bases[k]^inputs.exponents[i, k].  This is the shape of
// reference vec_chinese_units itself, whose exponent matrix is a ZM and whose
// FlxqV_factorback (defined at polarit2.c:1275-1298, called at
// buch2.c:1039) takes t_INT exponents (buch2.c:1026-1043, 1214-1240).
// The caller certifies the same two facts as above.
bool bounded_compact_reconstruct(
        BoundedCompactReconstructionReport& report, std::vector<Element>& out,
        const Order& order, const WideFactoredUnits& inputs,
        flint::FmpzConstRef certified_coordinate_bound,
        const BoundedCompactReconstructionOptions& options = {}) noexcept;

}  // namespace silex::detail
