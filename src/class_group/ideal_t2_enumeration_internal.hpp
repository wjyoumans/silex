#pragma once

#include <silex/class_group.hpp>
#include <silex/flint/arb_mat.hpp>

#include <atomic>
#include <vector>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;
class TestHookAccess;

}  // namespace silex::detail

namespace silex::detail::relation_search {

// Per-ideal T2 enumeration effort: the enumeration bound targets about
// 4 * kMaxFactorAttempts lattice points.  The factor-base honesty witness
// search may raise both caps in bounded stages (factor_base_honesty.cpp);
// relation search uses them unchanged.
inline constexpr slong kMaxFactorAttempts = 500;
inline constexpr slong kMaxElementSteps =
        4 * kMaxFactorAttempts * kMaxFactorAttempts;

// Starting working precisions (bits) of one T2 enumeration setup, retried at
// ceil(1.2 * precision) up to kT2EnumerationMaxPrecision.  The source starts
// its relation search at 200 bits and its unconditional class-group proof at
// 100, each with that retry (see "Class groups and order units" in
// docs/reference/algorithms_and_sources.rst).  The factor-base honesty witness
// search uses the proof value; every witness it finds is verified exactly, so
// the precision only orders the candidates.
inline constexpr slong kT2EnumerationRelationSearchInitialPrecision = 200;
inline constexpr slong kT2EnumerationProofInitialPrecision = 100;
inline constexpr slong kT2EnumerationMaxPrecision = 1000;

// How the T2 setup forms the reduced rows and their Gram matrix.  `arb`
// forms both in Arb (the relation-search route).  `double_rows` forms them in
// double from the exact reduced basis (order coordinates) and the order's
// Minkowski rows, and falls back to the whole `arb` setup at the same
// precision when the double data is not finite, a basis entry does not fit
// a double exactly, the double Gram-Schmidt pivot test fails, the product of
// the double Gram-Schmidt pivots differs from the exact d = |disc| N(I)^2 by
// more than 2^-20 in log2 (cancellation in the double rows), or no rounding
// scale passes the LLL-bound test.  It is taken at the start precision only;
// precision retries use `arb`.  The factor-base honesty witness search uses
// it: only Gram midpoints in double are read downstream and every witness is
// verified exactly, so the route only orders the candidates.
enum class T2EnumerationGram {
    arb,
    double_rows,
};

struct FiniteIdealT2EnumerationData {
    flint::FmpzMat basis{0, 0};
    std::vector<double> quadratic_form_data;
    double initial_bound_value = 0.0;
    // Rounding scale of the accepted reduced basis; 0 when the basis is the
    // first one with quadratic-form data (no scale was accepted).
    slong rounding_scale = 0;
    // Whether the double_rows route built this data (false after its
    // fallback to the Arb setup).
    bool double_rows = false;
};

// Test-only seam for the T2 double-row setup and its fallbacks.  `arb_only`
// makes the honesty witness search request the Arb setup; `nonfinite_rows`
// makes every double-row setup see a non-finite reduced row, so it falls
// back to the Arb setup; `exhausted` makes the honesty witness search treat
// every enumeration of double-row data as exhausted, so it re-enumerates the
// ideal from the Arb setup.  Only detail::TestHookAccess sets it; production
// code reads it with a relaxed atomic load.
class T2DoubleRowsTestSeams {
public:
    enum class Mode : int {
        off = 0,
        arb_only = 1,
        nonfinite_rows = 2,
        exhausted = 3,
    };
    static Mode mode() noexcept {
        return static_cast<Mode>(mode_.load(std::memory_order_relaxed));
    }

private:
    friend class detail::TestHookAccess;
    inline static std::atomic<int> mode_{0};
};

class FiniteQuadraticFormEnumerationContext {
public:
    bool reset(const std::vector<double>& quadratic_form_data, slong dimension) noexcept;
    bool start(double bound, slong max_element_steps,
               bool skip_first_scalar) noexcept;
    bool next() noexcept;
    bool current_row(flint::FmpzMatRef out) const noexcept;
    slong element_steps() const noexcept;
    bool element_step_limit_reached() const noexcept;

private:
    void clear() noexcept;
    double q(slong row, slong column) const noexcept;
    double v(slong index) const noexcept;
    bool step(slong index) noexcept;

    slong dimension_ = 0;
    std::vector<double> quadratic_form_data_;
    std::vector<slong> coefficients_;
    std::vector<slong> increments_;
    std::vector<double> lengths_;
    std::vector<double> tails_;
    slong level_ = 0;
    slong max_element_steps_ = 0;
    slong skip_first_ = 0;
    slong element_steps_ = 0;
    double bound_ = 0.0;
    bool started_ = false;
    bool exhausted_ = true;
    bool need_outer_step_ = false;
    bool element_step_limit_reached_ = false;
};

bool fmpz_mat_single_row_is_primitive(flint::FmpzMatConstRef row) noexcept;

bool build_finite_ideal_t2_enumeration_data_with_retry(
        FiniteIdealT2EnumerationData& out, const Ideal& ideal,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache,
        slong target_factor_attempts = kMaxFactorAttempts,
        slong initial_precision =
                kT2EnumerationRelationSearchInitialPrecision,
        T2EnumerationGram gram = T2EnumerationGram::arb) noexcept;

double arb_midpoint_double(flint::ArbConstRef value) noexcept;

// The midpoints of `rows` as doubles, row-major; false when one is not
// finite.
bool arb_mat_midpoints_double(std::vector<double>& out,
                              const flint::ArbMat& rows) noexcept;

// The T2 Gram matrix, row-major, of the rows reduced_basis * order_rows,
// formed in double: reduced_basis is the exact degree x degree basis in order
// coordinates and order_rows the row-major double Minkowski rows of the order
// basis.  False when a basis entry exceeds 2^53 in absolute value (so it is
// not exactly a double) or a row or Gram value is not finite.
bool reduced_basis_t2_gram_double(std::vector<double>& gram,
                                  flint::FmpzMatConstRef reduced_basis,
                                  const std::vector<double>& order_rows) noexcept;

}  // namespace silex::detail::relation_search
