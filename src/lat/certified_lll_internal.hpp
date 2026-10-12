#pragma once

#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>

#include <atomic>

namespace silex::detail {

class TestHookAccess;

// Test-only failure injected after the first phase runs: the phase's
// output is zeroed, so the fallback is exercised after the phase has
// changed the matrix.  first_phase: the phase status is also read as -1.
// certificate: the status is kept and the certificate runs on the zeroed
// matrix, which it does not certify.
enum class CertifiedLllTestFailure { none, first_phase, certificate };

// Read each time the first phase runs, at every call site.  Off by
// default and set only through detail::TestHookAccess
// (test/test_hook_access.hpp).  Relaxed atomics so a reader never races
// with a test setting them.
class CertifiedLllTestSeams {
public:
    static CertifiedLllTestFailure failure() noexcept {
        return failure_.load(std::memory_order_relaxed);
    }
    // Counts each first-phase run that an injected failure replaced.
    static void note_injected() noexcept {
        injected_.fetch_add(1, std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<CertifiedLllTestFailure> failure_{
            CertifiedLllTestFailure::none};
    inline static std::atomic<slong> injected_{0};
};

}  // namespace silex::detail

namespace silex::lat::detail {

// FLINT's first LLL phase with an exact certificate, returning what
// fmpz_lll returns.
//
// FLINT 3.6.0 (pinned b574c48): fmpz_lll(B, U, fl) is
// fmpz_lll_with_removal_ulll(B, U, 250, NULL, fl) (src/fmpz_lll/lll.c).
// For a Z_BASIS input with |fmpz_mat_max_bits| <= 250 that is exactly one
// fmpz_lll_wrapper_with_removal_knapsack(B, U, NULL, fl) call
// (lll_with_removal_ulll.c), which returns the
// fmpz_lll_d_with_removal_knapsack output when that phase does not return
// -1 and fmpz_lll_is_reduced(B, fl, 120) accepts it
// (wrapper_with_removal_knapsack.c).  FLINT's own check strips initial zero
// rows and ends in fmpz_mat_is_reduced; in this FLINT build its cheap
// double-precision stage is compiled out, so it pays a multiprecision
// stage first.
//
// The helper runs that phase and accepts its output only when an exact
// certificate holds.  The returned matrices are fmpz_lll's on every route
// because (1) a passing certificate excludes zero Gram-Schmidt vectors, so
// there are no zero rows to strip and FLINT's check accepts the same
// output; (2) when the phase returns -1 or the certificate fails, the
// inputs are restored and fmpz_lll runs from the start; (3) the double
// phase is deterministic.  Only the direction "certificate true implies
// FLINT accepts" is used.
//
// A GRAM context (the input is a Gram matrix) is reduced by fmpz_lll with
// one wrapper call at every size (lll_with_removal_ulll.c), and FLINT's
// check of the phase output is fmpz_lll_is_reduced's GRAM branch, which
// strips no rows and ends in fmpz_mat_is_reduced_gram (is_reduced.c).  The
// certificate is integer_gram_is_lll_reduced for small outputs: when it
// holds, every Gram-Schmidt norm is positive, fmpz_mat_is_reduced_gram
// computes the same exact mu and norms and accepts, so (1) to (3) carry
// over.  For larger outputs it is FLINT's own check, the same call on the
// same matrix as in fmpz_lll, so the decision is fmpz_lll's.  The transform
// does not steer the GRAM phase: in _fmpz_lll_d (lll_d.c) U is read only
// for its row count and in the update_b test, which needs more columns
// than rows and so never holds for a square Gram matrix; every other use
// writes the row operations into U (_fmpz_lll_check_babai in
// lll_check_babai.c, fmpz_mat_move_row in lll_d.c).

enum class CertifiedLllRoute {
    certified,   // The certified first phase is the result.
    fallback,    // The phase ran but was not certified; fmpz_lll ran.
    ineligible,  // Fewer than two rows, a non APPROX context, a Z_BASIS
                 // input with entries above 250 bits, a non square GRAM
                 // input, or a context that is neither; fmpz_lll ran.
};

// The certificate uses integer_basis_is_lll_reduced when the phase output
// has entries of at most this many bits and fmpz_mat_is_reduced otherwise.
// Both are exact; the cutoff only chooses the cheaper one.  Measured on a
// busy machine: the integer test was 1.5 to 23 times faster on every
// matrix the class/unit rows produced (at most 14 bits after the phase),
// and on dense reduced n x n matrices it stays at least as fast up to
// about 16 bits for n up to 81, while above about 20 bits it falls behind
// as n grows.
inline constexpr slong certified_lll_integer_test_max_bits = 16;

// For a GRAM context the certificate uses integer_gram_is_lll_reduced when
// the phase output's row count times its entry bit size is at most this,
// and FLINT's own check of the phase output (fmpz_lll_is_reduced at the
// wrapper's precision 120) otherwise.  On the scaled T2 Gram matrices of
// the class/unit corpus (busy machine), the integer test took at most 0.79
// of the time of FLINT's check up to a product of about 2500 (degree 24),
// about the same time near 3400 (degree 32), and 1.6 to 37 times as long
// from degree 40 up.
inline constexpr slong certified_lll_gram_integer_test_max_size = 3000;

// Exact test of FLINT's fmpz_mat_is_reduced_gram(gram, delta, eta) for a
// positive definite Gram matrix of at least two rows; false for fewer rows,
// a non square matrix, or a matrix that is not positive definite (FLINT's
// test may accept some of those, or divide by zero).  It reads only the
// lower triangle and the diagonal.  Equality with FLINT holds for
// eta^2 < delta, as below.  True implies fmpz_mat_is_reduced_gram for any
// parameters.
bool integer_gram_is_lll_reduced(flint::FmpzMatConstRef gram,
                                 double delta,
                                 double eta) noexcept;

// Exact test of FLINT's fmpz_mat_is_reduced(basis, delta, eta) for a basis
// of at least two rows; false for fewer rows.  Equality with FLINT holds
// for FLINT's parameter range eta^2 < delta; if delta <= eta^2 FLINT can
// accept a zero last Gram-Schmidt vector, which this test always rejects.
// True implies fmpz_mat_is_reduced for any parameters.
bool integer_basis_is_lll_reduced(flint::FmpzMatConstRef basis,
                                  double delta,
                                  double eta) noexcept;

// The exact certificate: fmpz_mat_is_reduced(basis, delta, eta) of
// `config`, decided by integer_basis_is_lll_reduced when the entries have
// at most certified_lll_integer_test_max_bits bits.  False for fewer than
// two rows.  For a GRAM `config`: integer_gram_is_lll_reduced up to
// certified_lll_gram_integer_test_max_size (unlike fmpz_mat_is_reduced_gram
// it never divides by a zero Gram-Schmidt norm), and above it FLINT's own
// fmpz_lll_is_reduced(basis, config, 120), the test fmpz_lll applies to
// the phase output, which is not an exact test but makes the helper accept
// exactly when fmpz_lll keeps the phase output.
bool lll_certificate_holds(flint::FmpzMatConstRef basis,
                           flint::FmpzLllConstRef config) noexcept;

// The first phase and the certificate only, for callers with their own
// fallback.  `transform` may be null; otherwise it has as many rows as
// `basis` and receives the phase's row operations, as in fmpz_lll.
// certified: `basis` and `transform` hold what fmpz_lll(basis, transform,
// config) returns for their input values.  fallback: the phase ran and was
// not certified; both hold unspecified values.  ineligible: both are
// unchanged.
CertifiedLllRoute lll_first_phase_certified(fmpz_mat_struct* basis,
                                            fmpz_mat_struct* transform,
                                            flint::FmpzLllConstRef config)
        noexcept;

// Replaces `basis` and `transform` (null or with as many rows as `basis`)
// by what fmpz_lll(basis, transform, config) returns: the certified first
// phase when it holds, otherwise fmpz_lll on the restored inputs.
CertifiedLllRoute certified_fmpz_lll(fmpz_mat_struct* basis,
                                     fmpz_mat_struct* transform,
                                     flint::FmpzLllConstRef config) noexcept;

}  // namespace silex::lat::detail
