#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nf.h>
#include <flint/nmod.h>

#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/prime_ideal.hpp>

namespace silex::detail {

// fmpz_get_nmod(value, mod), with the reduction done inline for a word-size
// value whose magnitude is below mod.n.
inline ulong fmpz_get_nmod_fast(const fmpz* value, nmod_t mod) noexcept {
    // A word-size value of magnitude below n is its own residue, or n minus
    // its magnitude when negative; zero stays zero.  Anything else takes
    // fmpz_get_nmod, which gives that same residue in [0, n).
    const fmpz word = *value;
    if (!COEFF_IS_MPZ(word)) {
        // |word| <= COEFF_MAX, so the negation cannot overflow.
        const ulong magnitude = word < 0 ? static_cast<ulong>(-word)
                                         : static_cast<ulong>(word);
        if (magnitude < mod.n) {
            return word < 0 && magnitude != 0 ? mod.n - magnitude : magnitude;
        }
    }
    return ::fmpz_get_nmod(value, mod);
}

// out = base^|exponent| mod modulus, as n_powmod2_ui_preinv gives it, for
// a modulus with precomputed inverse n_preinvert_limb(modulus).  Returns
// false only when exponent is the minimum slong, whose magnitude is not an
// slong.  |exponent| = 1 with base < modulus returns base directly.
bool pow_mod_abs_slong_ui(ulong& out,
                          ulong base,
                          slong exponent,
                          ulong modulus,
                          ulong modulus_inverse) noexcept;

// Returns the image of the NumberField defining generator, even when a direct
// maximal-quadratic prime stores its residue polynomial in the integral omega
// basis used by the order.
bool degree_one_prime_root_mod_p(flint::Fmpz& root,
                                 flint::Fmpz& p,
                                 const PrimeIdeal& prime) noexcept;
bool evaluate_fmpq_poly_at_degree_one_prime(
        flint::Fmpz& out,
        flint::FmpqPolyConstRef polynomial,
        flint::FmpzConstRef root,
        flint::FmpzConstRef p) noexcept;
// Smallest primitive root modulo the prime p, in fmpz arithmetic for any p.
// Returns false when p is not prime.
bool primitive_root_mod_prime(flint::Fmpz& out,
                              flint::FmpzConstRef p) noexcept;
// The same smallest primitive root for a prime p below the word size,
// factoring p - 1 once.  The caller guarantees that p is prime; returns 0
// only for p < 2.
ulong primitive_root_mod_prime_ui(ulong p) noexcept;
bool quotient_log_mod_prime_setup(flint::Fmpz& cofactor,
                                  flint::Fmpz& quotient_generator,
                                  flint::FmpzConstRef p,
                                  flint::FmpzConstRef ell) noexcept;
bool quotient_log_mod_prime_apply(flint::FmpzRef out,
                                  flint::FmpzConstRef value,
                                  flint::FmpzConstRef cofactor,
                                  flint::FmpzConstRef quotient_generator,
                                  flint::FmpzConstRef p,
                                  flint::FmpzConstRef ell) noexcept;
bool factored_value_at_degree_one_root(flint::Fmpz& out,
                                       const FactoredElement& element,
                                       flint::FmpzConstRef p,
                                       flint::FmpzConstRef root) noexcept;
// factored_value_at_degree_one_root in word arithmetic for a prime p below
// the word size: fails exactly where it does and otherwise gives the same
// residue.  Each factor is evaluated from its stored form num(x) / den in
// place (no polynomial copy), and the per-coefficient route of
// evaluate_fmpq_poly_at_degree_one_prime is taken only for a factor whose
// stored denominator is divisible by p.
bool factored_value_at_degree_one_root_nmod(ulong& out,
                                           const FactoredElement& element,
                                           ulong p,
                                           ulong root) noexcept;
// The reductions mod a word-size prime p that
// factored_value_at_degree_one_root_nmod makes for each factor of several
// FactoredElements (its stored denominator and numerator coefficients),
// made once so that the elements can be evaluated at several roots modulo
// p, such as the degree-one primes over one rational prime, with each
// distinct factor reduced and evaluated once for all the elements that
// share it.  Two factors are the same when they have the same field and the
// same stored form num(x) / den, coefficient for coefficient; both
// evaluation routes of factored_value_at_degree_one_root_nmod read only
// that stored form, so equal factors give equal residues and fail
// together.
//
// index(elements) finds the distinct factors and depends on the elements
// only, so it is made once for any number of primes; reduce(p) makes the
// reductions mod p; values_at_degree_one_root(out, root) then sets out[i]
// to factored_value_at_degree_one_root_nmod(out[i], elements[i], p, root)
// for every i and succeeds exactly when every one of those calls succeeds.
// Each product multiplies the factor powers of its element in the
// element's own factor order, as factored_value_at_degree_one_root_nmod
// does.  The index borrows `elements`, which must outlive it unchanged.
class FactoredElementsNmodReduction {
public:
    void index(std::span<const FactoredElement> elements) noexcept;
    void reduce(ulong p) noexcept;
    bool values_at_degree_one_root(std::span<ulong> out,
                                   ulong root) noexcept;
    // The number of distinct factors index() found.
    std::size_t distinct_factor_count() const noexcept {
        return factors_.size();
    }

private:
    // A distinct factor: its value (`factor`, in the field `field`) and,
    // after reduce(p), its stored denominator mod p and, when that is
    // nonzero, the `length` reduced numerator coefficients starting at
    // `offset` in coefficients_, lowest degree first.
    struct Factor {
        const Element* factor = nullptr;
        const nf_struct* field = nullptr;
        std::size_t offset = 0;
        std::size_t length = 0;
        ulong denominator = 0;
    };
    // A factor of an element with nonzero exponent: the exponent and the
    // index of its distinct factor in factors_.
    struct Entry {
        slong exponent = 0;
        std::size_t factor = 0;
    };

    std::size_t element_count_ = 0;
    bool indexed_ = false;
    bool reduced_ = false;
    ulong p_ = 0;
    ulong p_inverse_ = 0;
    nmod_t mod_ = {};
    std::vector<Factor> factors_;
    // The entries of element i are entries_[entry_starts_[i]] up to
    // entries_[entry_starts_[i + 1]], in the element's factor order.
    std::vector<Entry> entries_;
    std::vector<std::size_t> entry_starts_;
    std::vector<ulong> coefficients_;
    // Scratch for values_at_degree_one_root: each distinct factor's value
    // numerator / denominator at the root.
    std::vector<ulong> numerators_;
    std::vector<ulong> denominators_;
};
bool element_value_at_degree_one_root(flint::Fmpz& out,
                                      const Element& element,
                                      flint::FmpzConstRef p,
                                      flint::FmpzConstRef root) noexcept;

}  // namespace silex::detail
