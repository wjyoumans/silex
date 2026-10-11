#include <silex/element.hpp>

#include "element_internal.hpp"

#include <silex/archimedean.hpp>
#include <silex/diagnostics.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fq.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/flint/nmod_poly_factor.hpp>

#include <flint/fmpz_vec.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

// Private Hensel-lifting and root-reconstruction machinery for `Element`:
// exact square/power roots, pure-power Hensel lifting (`is_square`,
// `is_power`) and the cyclotomic root-of-unity finder
// (`detail::cyclotomic_root_hensel`).  Split out of `element.cpp` as one
// block because the cyclotomic search reuses the pure-power Hensel/
// reconstruction helpers, so the two are not separable without exporting
// all of it.
namespace silex {
namespace {

class FqPoly {
public:
    explicit FqPoly(const fq_ctx_t ctx) noexcept
        : ctx_(ctx) {
        fq_poly_init(value_, ctx_);
    }

    ~FqPoly() noexcept {
        fq_poly_clear(value_, ctx_);
    }

    FqPoly(const FqPoly&) = delete;
    FqPoly& operator=(const FqPoly&) = delete;

    fq_poly_t& raw() noexcept { return value_; }
    const fq_poly_t& raw() const noexcept { return value_; }

private:
    const fq_ctx_struct* ctx_;
    fq_poly_t value_;
};

class FqPolyFactor {
public:
    explicit FqPolyFactor(const fq_ctx_t ctx) noexcept
        : ctx_(ctx) {
        fq_poly_factor_init(value_, ctx_);
    }

    ~FqPolyFactor() noexcept {
        fq_poly_factor_clear(value_, ctx_);
    }

    FqPolyFactor(const FqPolyFactor&) = delete;
    FqPolyFactor& operator=(const FqPolyFactor&) = delete;

    fq_poly_factor_t& raw() noexcept { return value_; }
    const fq_poly_factor_t& raw() const noexcept { return value_; }

private:
    const fq_ctx_struct* ctx_;
    fq_poly_factor_t value_;
};

void multiplication_matrix(fmpq_mat_t matrix,
                           const nf_struct* field,
                           slong degree,
                           const nf_elem_t element) noexcept {
    flint::NfElem theta(field);
    flint::NfElem power(field);
    flint::NfElem next(field);
    flint::NfElem product(field);
    flint::FmpqPoly polynomial;
    flint::Fmpq coeff;

    nf_elem_gen(theta.raw(), field);
    nf_elem_one(power.raw(), field);

    for (slong j = 0; j < degree; ++j) {
        nf_elem_mul(product.raw(), element, power.raw(), field);
        nf_elem_get_fmpq_poly(polynomial.raw(), product.raw(), field);

        for (slong i = 0; i < degree; ++i) {
            fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial.raw(), i);
            fmpq_set(fmpq_mat_entry(matrix, i, j), coeff.raw());
        }

        nf_elem_mul(next.raw(), power.raw(), theta.raw(), field);
        nf_elem_swap(power.raw(), next.raw(), field);
    }
}

bool fmpq_sqrt(fmpq_t out, const fmpq_t value) noexcept {
    if (fmpq_sgn(value) < 0) {
        return false;
    }

    flint::Fmpz numerator;
    flint::Fmpz denominator;
    flint::Fmpz remainder;

    fmpz_sqrtrem(numerator.raw(), remainder.raw(), fmpq_numref(value));
    if (fmpz_is_zero(remainder.raw()) == 0) {
        return false;
    }

    fmpz_sqrtrem(denominator.raw(), remainder.raw(), fmpq_denref(value));
    if (fmpz_is_zero(remainder.raw()) == 0) {
        return false;
    }

    fmpq_set_fmpz_frac(out, numerator.raw(), denominator.raw());
    return true;
}

bool fmpz_root_signed(fmpz_t out, const fmpz_t value, slong exponent) noexcept {
    if (fmpz_sgn(value) < 0) {
        if ((exponent & 1) == 0) {
            return false;
        }

        flint::Fmpz positive;
        fmpz_neg(positive.raw(), value);
        const bool exact = fmpz_root(out, positive.raw(), exponent) != 0;
        if (exact) {
            fmpz_neg(out, out);
        }
        return exact;
    }

    return fmpz_root(out, value, exponent) != 0;
}

bool fmpq_root(fmpq_t out, const fmpq_t value, slong exponent) noexcept {
    if (exponent <= 0) {
        return false;
    }

    flint::Fmpz numerator;
    flint::Fmpz denominator;
    if (!fmpz_root_signed(numerator.raw(), fmpq_numref(value), exponent) ||
        fmpz_root(denominator.raw(), fmpq_denref(value), exponent) == 0) {
        return false;
    }

    fmpq_set_fmpz_frac(out, numerator.raw(), denominator.raw());
    return true;
}

bool fmpq_huge_root(fmpq_t out,
                    const fmpq_t value,
                    flint::FmpzConstRef exponent) noexcept {
    if (fmpz_is_one(fmpq_denref(value)) == 0) {
        return false;
    }
    if (fmpz_is_zero(fmpq_numref(value)) != 0) {
        fmpq_zero(out);
        return true;
    }
    if (fmpz_is_one(fmpq_numref(value)) != 0) {
        fmpq_one(out);
        return true;
    }
    if (fmpz_equal_si(fmpq_numref(value), -1) != 0 &&
        fmpz_is_odd(exponent.raw()) != 0) {
        fmpq_set_si(out, -1, 1);
        return true;
    }

    return false;
}

bool get_rational_constant(fmpq_t out, const Element& element) noexcept {
    flint::FmpqPoly polynomial;
    if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }
    if (fmpq_poly_degree(polynomial.raw()) > 0) {
        return false;
    }

    fmpq_poly_get_coeff_fmpq(out, polynomial.raw(), 0);
    return true;
}

bool set_rational(Element& element, const fmpq_t value) noexcept {
    flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 0, value);
    return element.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool set_quadratic_coeffs(Element& element,
                          const fmpq_t constant,
                          const fmpq_t linear) noexcept {
    flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 0, constant);
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, linear);
    return element.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial));
}

bool verified_square_root(Element& root,
                          const Element& input,
                          const fmpq_t constant,
                          const fmpq_t linear) noexcept {
    if (!set_quadratic_coeffs(root, constant, linear)) {
        return false;
    }

    Element check(*input.parent());
    return check.multiply(root, root) && check.equal(input);
}

bool quadratic_try_sqrt_m(Element& root,
                          const Element& input,
                          const fmpq_t a,
                          const fmpq_t b,
                          const fmpz_t d,
                          const fmpq_t m) noexcept {
    flint::Fmpq u;
    flint::Fmpq r;
    flint::Fmpq s;

    fmpq_add(u.raw(), a, m);
    fmpq_div_2exp(u.raw(), u.raw(), 1);

    if (!fmpq_sqrt(r.raw(), u.raw())) {
        return false;
    }

    if (fmpq_is_zero(r.raw()) == 0) {
        fmpq_div(s.raw(), b, r.raw());
        fmpq_div_2exp(s.raw(), s.raw(), 1);
        return verified_square_root(root, input, r.raw(), s.raw());
    }

    if (fmpq_is_zero(b) == 0) {
        return false;
    }

    fmpq_div_fmpz(u.raw(), a, d);
    if (!fmpq_sqrt(s.raw(), u.raw())) {
        return false;
    }

    return verified_square_root(root, input, r.raw(), s.raw());
}

bool fmpq_poly_reduce_modulus(flint::FmpzModPoly& out,
                              const fmpq_poly_t polynomial,
                              const flint::FmpzModCtx& ctx,
                              const fmpz_t modulus) noexcept {
    flint::Fmpq coeff;
    flint::Fmpz coeff_mod;
    flint::Fmpz inverse_denominator;

    fmpz_mod_poly_zero(out.raw(), ctx.raw());
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (fmpz_invmod(inverse_denominator.raw(),
                        fmpq_denref(coeff.raw()), modulus) == 0) {
            return false;
        }
        fmpz_mod_set_fmpz(coeff_mod.raw(), fmpq_numref(coeff.raw()),
                          ctx.raw());
        fmpz_mod_mul(coeff_mod.raw(), coeff_mod.raw(),
                     inverse_denominator.raw(), ctx.raw());
        fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, coeff_mod.raw(),
                                     ctx.raw());
    }
    return true;
}

bool fmpq_poly_is_monic_integral(const fmpq_poly_t polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree < 1) {
        return false;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (fmpz_is_one(fmpq_denref(coeff.raw())) == 0) {
            return false;
        }
    }

    fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, degree);
    return fmpq_is_one(coeff.raw()) != 0;
}

bool arb_abs_upper_bound_integer(flint::Fmpz& out,
                                 flint::ArbConstRef value,
                                 slong precision) noexcept {
    flint::Arf upper;
    flint::arb_get_abs_ubound_arf(upper, value, precision);
    if (!flint::arf_is_finite(upper)) {
        return false;
    }
    flint::arf_get_fmpz(out, upper, ARF_RND_CEIL);
    return fmpz_sgn(out.raw()) >= 0;
}

bool acb_abs_upper_bound_integer(flint::Fmpz& out,
                                 const acb_struct* value,
                                 slong precision) noexcept {
    flint::Arb absolute;
    flint::acb_abs(absolute, value, precision);
    return arb_abs_upper_bound_integer(out, flint::ArbConstRef(absolute),
                                       precision);
}

bool arb_midpoint_to_double(double& out, const arb_struct* value) noexcept {
    flint::Arf midpoint;
    flint::arf_set(midpoint, arb_midref(value));
    out = flint::arf_get_d(midpoint, ARF_RND_NEAR);
    return std::isfinite(out);
}

bool symmetric_eigenvalue_bounds(double& min_eigenvalue,
                                 double& max_eigenvalue,
                                 std::vector<double> matrix,
                                 slong size) noexcept {
    if (size <= 0 ||
        matrix.size() != static_cast<std::size_t>(size * size)) {
        return false;
    }

    double scale = 0.0;
    for (double value : matrix) {
        if (!std::isfinite(value)) {
            return false;
        }
        scale = std::max(scale, std::abs(value));
    }
    if (scale == 0.0) {
        return false;
    }

    const auto entry = [size](slong row, slong col) noexcept {
        return static_cast<std::size_t>(row * size + col);
    };

    const double tolerance = 1e-12 * std::max(1.0, scale);
    const slong max_iterations = 100 * size * size;
    bool converged = false;
    for (slong iteration = 0; iteration < max_iterations; ++iteration) {
        slong pivot_row = 0;
        slong pivot_col = 1;
        double max_offdiag = 0.0;
        for (slong row = 0; row < size; ++row) {
            for (slong col = row + 1; col < size; ++col) {
                const double value = std::abs(matrix[entry(row, col)]);
                if (value > max_offdiag) {
                    max_offdiag = value;
                    pivot_row = row;
                    pivot_col = col;
                }
            }
        }

        if (max_offdiag <= tolerance) {
            converged = true;
            break;
        }

        const double app = matrix[entry(pivot_row, pivot_row)];
        const double aqq = matrix[entry(pivot_col, pivot_col)];
        const double apq = matrix[entry(pivot_row, pivot_col)];
        if (apq == 0.0) {
            matrix[entry(pivot_row, pivot_col)] = 0.0;
            matrix[entry(pivot_col, pivot_row)] = 0.0;
            continue;
        }

        const double tau = (aqq - app) / (2.0 * apq);
        const double sign = tau < 0.0 ? -1.0 : 1.0;
        const double tangent =
                sign / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
        const double sine = tangent * cosine;

        for (slong k = 0; k < size; ++k) {
            if (k == pivot_row || k == pivot_col) {
                continue;
            }
            const double aik = matrix[entry(k, pivot_row)];
            const double akq = matrix[entry(k, pivot_col)];
            const double next_aik = cosine * aik - sine * akq;
            const double next_akq = sine * aik + cosine * akq;
            matrix[entry(k, pivot_row)] = next_aik;
            matrix[entry(pivot_row, k)] = next_aik;
            matrix[entry(k, pivot_col)] = next_akq;
            matrix[entry(pivot_col, k)] = next_akq;
        }

        const double cosine2 = cosine * cosine;
        const double sine2 = sine * sine;
        const double sine_cosine = sine * cosine;
        matrix[entry(pivot_row, pivot_row)] =
                cosine2 * app - 2.0 * sine_cosine * apq + sine2 * aqq;
        matrix[entry(pivot_col, pivot_col)] =
                sine2 * app + 2.0 * sine_cosine * apq + cosine2 * aqq;
        matrix[entry(pivot_row, pivot_col)] = 0.0;
        matrix[entry(pivot_col, pivot_row)] = 0.0;
    }

    if (!converged) {
        return false;
    }

    min_eigenvalue = std::numeric_limits<double>::infinity();
    max_eigenvalue = 0.0;
    for (slong i = 0; i < size; ++i) {
        const double eigenvalue = matrix[entry(i, i)];
        if (!std::isfinite(eigenvalue)) {
            return false;
        }
        min_eigenvalue = std::min(min_eigenvalue, eigenvalue);
        max_eigenvalue = std::max(max_eigenvalue, eigenvalue);
    }

    return min_eigenvalue > 0.0 && std::isfinite(max_eigenvalue);
}

bool compute_equation_order_norm_change_constants(
        double& c1,
        double& c2,
        const NumberField& parent,
        slong precision) noexcept {
    const slong degree = parent.degree();
    if (degree <= 0) {
        return false;
    }

    EmbeddingContext embeddings(parent);
    if (!embeddings.refine(precision)) {
        return false;
    }

    std::vector<double> minkowski(
            static_cast<std::size_t>(degree * degree), 0.0);
    flint::ArbMat image(1, degree);
    flint::FmpqPoly basis_polynomial;
    Element basis(parent);
    if (!basis.is_defined()) {
        return false;
    }

    for (slong row = 0; row < degree; ++row) {
        fmpq_poly_zero(basis_polynomial.raw());
        fmpq_poly_set_coeff_si(basis_polynomial.raw(), row, 1);
        if (!basis.set_fmpq_poly(flint::FmpqPolyConstRef(basis_polynomial)) ||
            !minkowski_embedding(flint::ArbMatRef(image), embeddings, basis,
                                 MinkowskiEmbeddingMode::weighted,
                                 precision)) {
            return false;
        }

        for (slong col = 0; col < degree; ++col) {
            double value = 0.0;
            if (!arb_midpoint_to_double(
                        value, arb_mat_entry(image.raw(), 0, col))) {
                return false;
            }
            minkowski[static_cast<std::size_t>(row * degree + col)] = value;
        }
    }

    std::vector<double> gram(static_cast<std::size_t>(degree * degree), 0.0);
    for (slong row = 0; row < degree; ++row) {
        for (slong col = 0; col < degree; ++col) {
            double value = 0.0;
            for (slong k = 0; k < degree; ++k) {
                value += minkowski[
                                 static_cast<std::size_t>(row * degree + k)] *
                         minkowski[
                                 static_cast<std::size_t>(col * degree + k)];
            }
            gram[static_cast<std::size_t>(row * degree + col)] = value;
        }
    }

    double min_eigenvalue = 0.0;
    double max_eigenvalue = 0.0;
    if (!symmetric_eigenvalue_bounds(min_eigenvalue, max_eigenvalue,
                                     std::move(gram), degree)) {
        return false;
    }

    c1 = max_eigenvalue;
    c2 = 1.0 / min_eigenvalue;
    return std::isfinite(c1) && std::isfinite(c2) &&
           c1 > 0.0 && c2 > 0.0;
}

// Memo of the radicand-independent Hensel data: the norm-change constants
// below and the reconstruction lattice of power_reconstruction_data.  Every
// is_square / is_power call on a field recomputes both from scratch although
// neither reads the radicand: the constants depend only on the field and the
// working precision, the reconstruction data only on the lifted prime factor
// and its modulus p^k.  The class-group saturation stage asks thousands of
// such questions in one field, mostly with the same prime and lift exponents.
//
// Exactness: an entry is keyed by every input its computation reads, compared
// for exact equality, and both computations are deterministic, so a hit
// returns bit for bit what recomputation would return.  Only successful
// computations are stored; a failure is recomputed each time.
//
// Lifetime and thread safety: both caches are thread_local, so a thread only
// ever sees entries it computed itself, no locking is needed, and entries are
// released when the thread exits.  Each cache holds a bounded number of
// entries (round-robin replacement once full), so memory stays bounded however
// many fields a thread visits.  A lookup copies the entry out and an insertion
// happens only after the computation has returned, so no reference into a
// cache is held across a call that could reach the cache again.
thread_local detail::HenselDataCacheCounters hensel_data_cache_counters_tls;

template <class Entry>
class BoundedEntryCache {
public:
    explicit BoundedEntryCache(std::size_t capacity) noexcept
        : capacity_(capacity) {}

    std::vector<Entry>& entries() noexcept { return entries_; }

    void clear() noexcept {
        entries_.clear();
        next_replacement_ = 0;
    }

    void insert(Entry entry) noexcept {
        if (entries_.size() < capacity_) {
            entries_.push_back(std::move(entry));
            return;
        }
        entries_[next_replacement_] = std::move(entry);
        next_replacement_ = (next_replacement_ + 1) % capacity_;
    }

private:
    std::size_t capacity_;
    std::size_t next_replacement_ = 0;
    std::vector<Entry> entries_;
};

// Key: the defining polynomial, backend and precision, which are all that
// compute_equation_order_norm_change_constants reads (EmbeddingContext and
// the power-basis elements are built from the field alone).
struct NormChangeConstantsEntry {
    flint::FmpqPoly polynomial;
    NumberFieldBackendKind backend = NumberFieldBackendKind::generic;
    slong precision = 0;
    double c1 = 0.0;
    double c2 = 0.0;
};

// A thread rarely works in more than a few fields at once.
constexpr std::size_t kNormChangeConstantsCacheCapacity = 8;

BoundedEntryCache<NormChangeConstantsEntry>&
norm_change_constants_cache() noexcept {
    thread_local BoundedEntryCache<NormChangeConstantsEntry> cache(
            kNormChangeConstantsCacheCapacity);
    return cache;
}

// c1, c2 of compute_equation_order_norm_change_constants, memoized per thread
// (see the cache comment above).
bool equation_order_norm_change_constants(
        double& c1,
        double& c2,
        const NumberField& parent,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_norm_change_constants");
    const nf_struct* raw_field = parent.raw_flint_field();
    if (raw_field == nullptr) {
        return false;
    }
    const NumberFieldBackendKind backend = parent.backend_kind();
    auto& cache = norm_change_constants_cache();
    for (const NormChangeConstantsEntry& entry : cache.entries()) {
        if (entry.precision == precision && entry.backend == backend &&
            fmpq_poly_equal(entry.polynomial.raw(), raw_field->pol) != 0) {
            c1 = entry.c1;
            c2 = entry.c2;
            ++hensel_data_cache_counters_tls.norm_change_constants_hits;
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.power_hensel_norm_change_constants.cache_hit");
            return true;
        }
    }

    ++hensel_data_cache_counters_tls.norm_change_constants_misses;
    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::element,
            "element.power_hensel_norm_change_constants.cache_miss");
    if (!compute_equation_order_norm_change_constants(c1, c2, parent,
                                                      precision)) {
        return false;
    }
    NormChangeConstantsEntry entry;
    fmpq_poly_set(entry.polynomial.raw(), raw_field->pol);
    entry.backend = backend;
    entry.precision = precision;
    entry.c1 = c1;
    entry.c2 = c2;
    cache.insert(std::move(entry));
    return true;
}

slong bit_length_slong(slong value) noexcept {
    slong bits = 0;
    while (value > 0) {
        ++bits;
        value >>= 1;
    }
    return bits;
}

bool power_leq_slong(slong base, slong exponent, slong limit) noexcept {
    slong value = 1;
    for (slong i = 0; i < exponent; ++i) {
        if (base != 0 && value > limit / base) {
            return false;
        }
        value *= base;
    }
    return value <= limit;
}

slong floor_integer_root_slong(slong value, slong exponent) noexcept {
    slong root = 0;
    while (power_leq_slong(root + 1, exponent, value)) {
        ++root;
    }
    return root;
}

slong equation_order_c3(slong degree) noexcept {
    if (degree == 1) {
        return 1;
    }

    slong exponent = bit_length_slong(degree);
    exponent += exponent % 2;
    return floor_integer_root_slong(degree, exponent) + 1;
}

// Field part of pure_power_conjugate_bounds: the embedding roots of the
// field at the working precision and the conjugates of f'(theta).  Neither
// reads the radicand, yet every lifting-exponent computation used to isolate
// the roots afresh (acb_poly_find_roots through a local EmbeddingContext).
//
// Key: the defining polynomial, backend, working precision and the
// derivative polynomial whose conjugates are stored.  The roots are those of
// EmbeddingContext::refine(precision) on a fresh context, which reads only
// the defining polynomial, and the entry is used only at exactly the
// precision it was computed at, so a hit never serves a lower-precision
// enclosure and evaluates bit for bit as a fresh context would.  The entry's
// context owns its own copy of the field handle, which keeps that field's
// (small) data alive until the entry is replaced or the thread exits; an
// element of another field object with the same polynomial is evaluated
// through its polynomial.
struct ConjugateBoundsFieldEntry {
    flint::FmpqPoly polynomial;
    flint::FmpqPoly derivative_polynomial;
    NumberFieldBackendKind backend = NumberFieldBackendKind::generic;
    slong precision = 0;
    EmbeddingContext embeddings;
    flint::AcbVec derivative_values{0};
};

constexpr std::size_t kConjugateBoundsFieldCacheCapacity = 8;

BoundedEntryCache<ConjugateBoundsFieldEntry>&
conjugate_bounds_field_cache() noexcept {
    thread_local BoundedEntryCache<ConjugateBoundsFieldEntry> cache(
            kConjugateBoundsFieldCacheCapacity);
    return cache;
}

bool compute_conjugate_bounds_field_entry(
        ConjugateBoundsFieldEntry& entry,
        const NumberField& parent,
        const flint::FmpqPoly& derivative_polynomial,
        slong precision) noexcept {
    const nf_struct* raw_field = parent.raw_flint_field();
    if (raw_field == nullptr || !entry.embeddings.define(parent) ||
        !entry.embeddings.refine(precision)) {
        return false;
    }

    Element derivative(parent);
    if (!derivative.is_defined() ||
        !derivative.set_fmpq_poly(flint::FmpqPolyConstRef(
                derivative_polynomial))) {
        return false;
    }

    entry.derivative_values = flint::AcbVec(parent.degree());
    if (!entry.embeddings.evaluate_all(
                flint::AcbVecRef(entry.derivative_values), derivative,
                precision)) {
        return false;
    }
    fmpq_poly_set(entry.polynomial.raw(), raw_field->pol);
    fmpq_poly_set(entry.derivative_polynomial.raw(),
                  derivative_polynomial.raw());
    entry.backend = parent.backend_kind();
    entry.precision = precision;
    return true;
}

// Conjugates of `input` and of f'(theta) at `precision`, from the field memo
// above.  No reference into a cache is held across a call that could reach
// a Hensel cache again: only EmbeddingContext and Element operations run
// while the entry is in use.
bool conjugate_bounds_values(flint::AcbVec& input_values,
                             flint::AcbVec& derivative_values,
                             Signature& sig,
                             const Element& input,
                             const flint::FmpqPoly& derivative_polynomial,
                             slong precision) noexcept {
    const NumberField* parent = input.parent();
    const nf_struct* raw_field = parent->raw_flint_field();
    if (raw_field == nullptr) {
        return false;
    }
    const NumberFieldBackendKind backend = parent->backend_kind();
    const auto evaluate_with = [&](ConjugateBoundsFieldEntry& entry) noexcept {
        const NumberField* entry_parent = entry.embeddings.parent();
        if (entry_parent == nullptr) {
            return false;
        }
        sig = entry.embeddings.signature();
        if (!derivative_values.set_from(
                    flint::AcbVecConstRef(entry.derivative_values))) {
            return false;
        }
        if (input.has_parent(*entry_parent)) {
            return entry.embeddings.evaluate_all(
                    flint::AcbVecRef(input_values), input, precision);
        }
        flint::FmpqPoly input_polynomial;
        Element converted(*entry_parent);
        return input.get_fmpq_poly(flint::FmpqPolyRef(input_polynomial)) &&
               converted.is_defined() &&
               converted.set_fmpq_poly(
                       flint::FmpqPolyConstRef(input_polynomial)) &&
               entry.embeddings.evaluate_all(flint::AcbVecRef(input_values),
                                             converted, precision);
    };

    auto& cache = conjugate_bounds_field_cache();
    for (ConjugateBoundsFieldEntry& entry : cache.entries()) {
        if (entry.precision == precision && entry.backend == backend &&
            fmpq_poly_equal(entry.polynomial.raw(), raw_field->pol) != 0 &&
            fmpq_poly_equal(entry.derivative_polynomial.raw(),
                            derivative_polynomial.raw()) != 0) {
            ++hensel_data_cache_counters_tls.conjugate_bounds_field_hits;
            return evaluate_with(entry);
        }
    }

    ++hensel_data_cache_counters_tls.conjugate_bounds_field_misses;
    ConjugateBoundsFieldEntry entry;
    if (!compute_conjugate_bounds_field_entry(entry, *parent,
                                              derivative_polynomial,
                                              precision) ||
        !evaluate_with(entry)) {
        return false;
    }
    cache.insert(std::move(entry));
    return true;
}

bool pure_power_conjugate_bounds(
        std::vector<flint::Arb>& bounds,
        const Element& input,
        const flint::FmpqPoly& derivative_polynomial,
        slong exponent,
        slong precision) noexcept {
    const NumberField* parent = input.parent();
    if (parent == nullptr || exponent <= 1) {
        return false;
    }

    const slong degree = parent->degree();
    flint::AcbVec input_values(degree);
    flint::AcbVec derivative_values(degree);
    Signature sig;
    if (!conjugate_bounds_values(input_values, derivative_values, sig, input,
                                 derivative_polynomial, precision)) {
        return false;
    }

    const slong r1 = sig.r1();
    const slong places = sig.r1() + sig.r2();
    bounds.clear();
    bounds.reserve(static_cast<std::size_t>(places));
    for (slong place = 0; place < places; ++place) {
        const slong root_index = place < r1 ? place : r1 + 2 * (place - r1);
        flint::Fmpz input_bound_integer;
        if (!acb_abs_upper_bound_integer(input_bound_integer,
                                         input_values.data() + root_index,
                                         precision)) {
            return false;
        }

        flint::Arb root_bound;
        flint::Arb derivative_absolute;
        flint::arb_set_fmpz(root_bound, input_bound_integer);
        flint::arb_root_ui(root_bound, root_bound,
                           static_cast<ulong>(exponent), precision);
        flint::acb_abs(derivative_absolute,
                       derivative_values.data() + root_index, precision);
        flint::arb_mul(root_bound, root_bound, derivative_absolute,
                       precision);
        if (!flint::arb_is_finite(root_bound)) {
            return false;
        }
        bounds.push_back(std::move(root_bound));
    }

    return bounds.size() == static_cast<std::size_t>(places);
}

// Friedrich-Fieker lifting exponent for reconstructing, in the equation
// order of the monic integral field `parent`, an element x whose scaled
// conjugates satisfy |f'(theta)^(i) x^(i)| <= root_bounds[i] at the r1 + r2
// places, from a prime of residue degree `residue_degree` above `prime`.
//
// Source trace: reference implementation
// `src/NumFieldOrd/NfOrd/Hensel.jl:_lifting_expo` with `any_order(K)` equal
// to the equation order for this monic-integral private root path.
bool friedrich_fieker_lifting_exponent(
        slong& out,
        const NumberField& parent_field,
        const std::vector<flint::Arb>& root_bounds,
        const flint::Fmpz& prime,
        slong residue_degree,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    const NumberField* parent = &parent_field;
    if (residue_degree <= 0 || fmpz_fits_si(prime.raw()) == 0 ||
        fmpz_sgn(prime.raw()) <= 0) {
        return false;
    }

    double c1 = 0.0;
    double c2 = 0.0;
    if (!equation_order_norm_change_constants(
                c1, c2, *parent, precision, diagnostics)) {
        return false;
    }

    Signature sig;
    if (!signature(sig, *parent)) {
        return false;
    }
    const slong r1 = sig.r1();
    const slong r2 = sig.r2();
    if (root_bounds.size() != static_cast<std::size_t>(r1 + r2)) {
        return false;
    }

    flint::Arb bd;
    flint::Arb square;
    flint::Arb term;
    flint::arb_zero(bd);
    for (slong i = 0; i < r1; ++i) {
        flint::arb_sqr(square, root_bounds[static_cast<std::size_t>(i)],
                       precision);
        flint::arb_add(bd, bd, square, precision);
    }
    for (slong i = 0; i < r2; ++i) {
        flint::arb_sqr(
                square,
                root_bounds[static_cast<std::size_t>(r1 + i)],
                precision);
        flint::arb_mul_ui(term, square, 2, precision);
        flint::arb_add(bd, bd, term, precision);
    }

    flint::Arb one;
    flint::Arb boundt2;
    flint::arb_one(one);
    flint::arb_max(boundt2, bd, one, precision);

    const slong degree = parent->degree();
    const slong c3 = equation_order_c3(degree);
    flint::Arb factor;
    flint::Arb tmp;
    flint::Arb c2_arb;
    flint::Arb c3_arb;
    flint::arb_set_d(factor, c1);
    flint::arb_set_d(c2_arb, c2);
    flint::arb_set_ui(c3_arb, static_cast<ulong>(c3));
    flint::arb_mul(factor, factor, c2_arb, precision);
    flint::arb_mul(factor, factor, c3_arb, precision);
    flint::arb_mul(factor, factor, boundt2, precision);
    flint::arb_mul_2exp_si(factor, factor, degree * (degree - 1) / 2 + 2);
    flint::arb_div_ui(factor, factor, static_cast<ulong>(degree), precision);
    if (!flint::arb_is_finite(factor)) {
        return false;
    }

    flint::Fmpz tmp_integer;
    if (!arb_abs_upper_bound_integer(tmp_integer,
                                     flint::ArbConstRef(factor),
                                     precision) ||
        fmpz_is_zero(tmp_integer.raw()) != 0) {
        return false;
    }

    flint::Arb tmp_arb;
    flint::Arb log_tmp;
    flint::Arb log_prime;
    flint::Arb numerator;
    flint::Arb denominator;
    flint::Arb boundk;
    flint::arb_set_fmpz(tmp_arb, tmp_integer);
    flint::arb_log(log_tmp, tmp_arb, precision);
    flint::arb_log_fmpz(log_prime, flint::FmpzConstRef(prime), precision);
    flint::arb_mul_ui(numerator, log_tmp, static_cast<ulong>(degree),
                      precision);
    flint::arb_mul_ui(denominator, log_prime,
                      static_cast<ulong>(2 * residue_degree), precision);
    flint::arb_div(boundk, numerator, denominator, precision);
    if (!flint::arb_is_finite(boundk)) {
        return false;
    }

    flint::Fmpz lift_exponent;
    if (!arb_abs_upper_bound_integer(lift_exponent,
                                     flint::ArbConstRef(boundk),
                                     precision) ||
        fmpz_fits_si(lift_exponent.raw()) == 0) {
        return false;
    }

    out = std::max<slong>(2, fmpz_get_si(lift_exponent.raw()));
    return out >= 2;
}

// Friedrich-Fieker lifting exponent of one radicand, memoized per thread.
// The class-group saturation stage tests the same radicand many times (the
// same square class reached through different relation combinations), and
// each test recomputed the conjugate bounds and the exponent.
//
// Key: everything pure_power_lifting_exponent reads, compared for exact
// equality: the defining polynomial and backend, the derivative polynomial,
// the radicand's power-basis coordinates, the root exponent, the prime and
// the residue degree (the working precision is the fixed constant below).
// The computation is deterministic, so a hit returns exactly what
// recomputation would.  Only successful computations are stored.
struct LiftingExponentEntry {
    flint::FmpqPoly polynomial;
    flint::FmpqPoly derivative_polynomial;
    flint::FmpqPoly element;
    NumberFieldBackendKind backend = NumberFieldBackendKind::generic;
    flint::Fmpz prime;
    slong residue_degree = 0;
    slong exponent = 0;
    slong lift_exponent = 0;
};

constexpr std::size_t kLiftingExponentCacheCapacity = 64;

BoundedEntryCache<LiftingExponentEntry>& lifting_exponent_cache() noexcept {
    thread_local BoundedEntryCache<LiftingExponentEntry> cache(
            kLiftingExponentCacheCapacity);
    return cache;
}

bool pure_power_lifting_exponent(
        slong& out,
        const Element& input,
        const flint::FmpqPoly& derivative_polynomial,
        const flint::Fmpz& prime,
        slong residue_degree,
        slong exponent,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_lifting_exponent");
    const NumberField* parent = input.parent();
    if (parent == nullptr || residue_degree <= 0 || exponent <= 1 ||
        fmpz_fits_si(prime.raw()) == 0 || fmpz_sgn(prime.raw()) <= 0) {
        return false;
    }
    const nf_struct* raw_field = parent->raw_flint_field();
    flint::FmpqPoly element_polynomial;
    if (raw_field == nullptr ||
        !input.get_fmpq_poly(flint::FmpqPolyRef(element_polynomial))) {
        return false;
    }
    const NumberFieldBackendKind backend = parent->backend_kind();
    for (const LiftingExponentEntry& entry :
         lifting_exponent_cache().entries()) {
        if (entry.exponent == exponent &&
            entry.residue_degree == residue_degree &&
            entry.backend == backend &&
            fmpz_equal(entry.prime.raw(), prime.raw()) != 0 &&
            fmpq_poly_equal(entry.element.raw(), element_polynomial.raw()) !=
                    0 &&
            fmpq_poly_equal(entry.polynomial.raw(), raw_field->pol) != 0 &&
            fmpq_poly_equal(entry.derivative_polynomial.raw(),
                            derivative_polynomial.raw()) != 0) {
            out = entry.lift_exponent;
            ++hensel_data_cache_counters_tls.lifting_exponent_hits;
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.power_hensel_lifting_exponent.cache_hit");
            return true;
        }
    }

    ++hensel_data_cache_counters_tls.lifting_exponent_misses;
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_lifting_exponent.cache_miss");
    constexpr slong precision = 128;
    std::vector<flint::Arb> root_bounds;
    if (!pure_power_conjugate_bounds(
                root_bounds, input, derivative_polynomial, exponent,
                precision) ||
        !friedrich_fieker_lifting_exponent(out, *parent, root_bounds, prime,
                                           residue_degree, precision,
                                           diagnostics)) {
        return false;
    }

    LiftingExponentEntry entry;
    fmpq_poly_set(entry.polynomial.raw(), raw_field->pol);
    fmpq_poly_set(entry.derivative_polynomial.raw(),
                  derivative_polynomial.raw());
    entry.element = std::move(element_polynomial);
    entry.backend = backend;
    fmpz_set(entry.prime.raw(), prime.raw());
    entry.residue_degree = residue_degree;
    entry.exponent = exponent;
    entry.lift_exponent = out;
    lifting_exponent_cache().insert(std::move(entry));
    return true;
}

bool lifting_chain(std::vector<slong>& out,
                         slong target_exponent) noexcept {
    if (target_exponent < 2) {
        target_exponent = 2;
    }

    out.clear();
    slong exponent = target_exponent;
    while (true) {
        out.push_back(exponent);
        if (exponent <= 1) {
            break;
        }
        exponent = (exponent + 1) / 2;
    }

    std::reverse(out.begin(), out.end());
    return !out.empty() && out.front() == 1 && out.back() == target_exponent;
}

void fmpz_mod_poly_center_coeff(fmpz_t out,
                                const fmpz_mod_poly_t polynomial,
                                slong index,
                                const fmpz_t modulus,
                                const fmpz_mod_ctx_t ctx) noexcept {
    fmpz_mod_poly_get_coeff_fmpz(out, polynomial, index, ctx);

    flint::Fmpz twice;
    fmpz_mul_2exp(twice.raw(), out, 1);
    if (fmpz_cmp(twice.raw(), modulus) > 0) {
        fmpz_sub(out, out, modulus);
    }
}

bool fmpz_mod_poly_to_centered_fmpq_poly(flint::FmpqPoly& out,
                                         const flint::FmpzModPoly& input,
                                         slong degree,
                                         const fmpz_t modulus,
                                         const flint::FmpzModCtx& ctx) noexcept {
    fmpq_poly_zero(out.raw());
    flint::Fmpz coeff;
    for (slong i = 0; i < degree; ++i) {
        fmpz_mod_poly_center_coeff(coeff.raw(), input.raw(), i, modulus,
                                   ctx.raw());
        fmpq_poly_set_coeff_fmpz(out.raw(), i, coeff.raw());
    }
    return true;
}

bool fmpz_mod_poly_factorization_is_squarefree(
        const flint::FmpzModPolyFactor& factorization) noexcept {
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        if (factorization.raw()->exp[i] != 1) {
            return false;
        }
    }
    return true;
}

enum class ResidueSquareStatus {
    unsupported,
    square,
    nonsquare
};

[[maybe_unused]] ResidueSquareStatus residue_field_square_status(
        const flint::FmpzModPoly& element,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx,
        const fmpz_t prime) noexcept {
    const slong residue_degree = fmpz_mod_poly_degree(modulus.raw(), ctx.raw());
    if (residue_degree <= 0) {
        return ResidueSquareStatus::unsupported;
    }

    flint::FmpzModPoly reduced(ctx);
    fmpz_mod_poly_rem(reduced.raw(), element.raw(), modulus.raw(), ctx.raw());
    if (fmpz_mod_poly_is_zero(reduced.raw(), ctx.raw()) != 0) {
        // reference skips nonsquarefree reductions of y^2 - a; for odd residue
        // characteristic this is exactly the a == 0 residue case.
        return ResidueSquareStatus::unsupported;
    }

    flint::Fmpz exponent;
    fmpz_pow_ui(exponent.raw(), prime, static_cast<ulong>(residue_degree));
    fmpz_sub_ui(exponent.raw(), exponent.raw(), 1);
    fmpz_fdiv_q_2exp(exponent.raw(), exponent.raw(), 1);

    flint::FmpzModPoly power(ctx);
    flint::FmpzModPoly one(ctx);
    fmpz_mod_poly_powmod_fmpz_binexp(
            power.raw(), reduced.raw(), exponent.raw(), modulus.raw(),
            ctx.raw());
    fmpz_mod_poly_one(one.raw(), ctx.raw());
    if (fmpz_mod_poly_equal(power.raw(), one.raw(), ctx.raw()) != 0) {
        return ResidueSquareStatus::square;
    }
    return ResidueSquareStatus::nonsquare;
}

bool residue_field_square_root(flint::FmpzModPoly& out,
                               const flint::FmpzModPoly& element,
                               const flint::FmpzModPoly& modulus,
                               const flint::FmpzModCtx& ctx) noexcept {
    if (!out.is_initialized() || !element.is_initialized() ||
        !modulus.is_initialized()) {
        return false;
    }

    flint::FmpzModPoly reduced(ctx);
    if (!reduced.is_initialized()) {
        return false;
    }
    fmpz_mod_poly_rem(reduced.raw(), element.raw(), modulus.raw(), ctx.raw());
    if (fmpz_mod_poly_is_zero(reduced.raw(), ctx.raw()) != 0) {
        return false;
    }

    flint::FqCtx field(modulus.raw(), ctx.raw(), "z");
    flint::Fq value(field);
    flint::Fq sqrt_value(field);
    if (!field.is_initialized() || !value.is_initialized() ||
        !sqrt_value.is_initialized()) {
        return false;
    }
    fq_set_fmpz_mod_poly(value.raw(), reduced.raw(), field.raw());
    if (fq_sqrt(sqrt_value.raw(), value.raw(), field.raw()) == 0) {
        return false;
    }

    fq_get_fmpz_mod_poly(out.raw(), sqrt_value.raw(), field.raw());
    fmpz_mod_poly_rem(out.raw(), out.raw(), modulus.raw(), ctx.raw());

    flint::FmpzModPoly check(ctx);
    if (!check.is_initialized()) {
        return false;
    }
    fmpz_mod_poly_mulmod(check.raw(), out.raw(), out.raw(), modulus.raw(),
                         ctx.raw());
    return fmpz_mod_poly_equal(check.raw(), reduced.raw(), ctx.raw()) != 0;
}

bool fmpz_mod_poly_change_modulus(flint::FmpzModPoly& out,
                                  const flint::FmpzModPoly& input,
                                  const flint::FmpzModCtx& out_ctx) noexcept {
    if (!out.is_initialized() || !input.is_initialized()) {
        return false;
    }

    flint::FmpzPoly lifted;
    fmpz_mod_poly_get_fmpz_poly(lifted.raw(), input.raw(), input.context());
    fmpz_mod_poly_set_fmpz_poly(out.raw(), lifted.raw(), out_ctx.raw());
    return true;
}

bool pure_square_hensel_step(flint::FmpzModPoly& next_inverse_root,
                                   const flint::FmpzModPoly& inverse_root,
                                   const flint::FmpzModPoly& radicand,
                                   const flint::FmpzModPoly& modulus,
                                   const flint::FmpzModCtx& ctx,
                                   const fmpz_t precision_modulus) noexcept {
    flint::Fmpz two;
    flint::Fmpz inverse_two;
    fmpz_set_ui(two.raw(), 2);
    if (fmpz_invmod(inverse_two.raw(), two.raw(), precision_modulus) == 0) {
        return false;
    }

    flint::Fmpz one_plus_inverse_two;
    fmpz_add_ui(one_plus_inverse_two.raw(), inverse_two.raw(), 1);

    flint::FmpzModPoly term1(ctx);
    flint::FmpzModPoly power(ctx);
    flint::FmpzModPoly product(ctx);
    flint::FmpzModPoly term2(ctx);
    if (!term1.is_initialized() || !power.is_initialized() ||
        !product.is_initialized() || !term2.is_initialized()) {
        return false;
    }

    // reference's pure branch applies Newton to the inverse root:
    // r <- r * (1 + 1/2) - a * (1/2) * r^3 modulo <p^e, g>.
    fmpz_mod_poly_scalar_mul_fmpz(term1.raw(), inverse_root.raw(),
                                  one_plus_inverse_two.raw(), ctx.raw());
    fmpz_mod_poly_powmod_ui_binexp(power.raw(), inverse_root.raw(), 3,
                                   modulus.raw(), ctx.raw());
    fmpz_mod_poly_mulmod(product.raw(), radicand.raw(), power.raw(),
                         modulus.raw(), ctx.raw());
    fmpz_mod_poly_scalar_mul_fmpz(term2.raw(), product.raw(),
                                  inverse_two.raw(), ctx.raw());
    fmpz_mod_poly_sub(next_inverse_root.raw(), term1.raw(), term2.raw(),
                      ctx.raw());
    fmpz_mod_poly_rem(next_inverse_root.raw(), next_inverse_root.raw(),
                      modulus.raw(), ctx.raw());
    return true;
}

bool scaled_root_candidate(Element& candidate,
                                 const Element& input,
                                 const flint::FmpzModPoly& scaled_root,
                                 const flint::FmpzModCtx& ctx,
                                 const fmpz_t precision_modulus,
                                 const flint::FmpqPoly& derivative_polynomial) noexcept {
    const NumberField* parent = input.parent();
    if (parent == nullptr || !candidate.has_parent(*parent)) {
        return false;
    }

    flint::FmpqPoly centered;
    fmpz_mod_poly_to_centered_fmpq_poly(centered, scaled_root,
                                        parent->degree(), precision_modulus,
                                        ctx);

    Element numerator(*parent);
    Element denominator(*parent);
    Element inverse_denominator(*parent);
    Element check(*parent);
    if (!numerator.is_defined() || !denominator.is_defined() ||
        !inverse_denominator.is_defined() || !check.is_defined() ||
        !numerator.set_fmpq_poly(flint::FmpqPolyConstRef(centered)) ||
        !denominator.set_fmpq_poly(flint::FmpqPolyConstRef(
                derivative_polynomial)) ||
        !inverse_denominator.invert(denominator) ||
        !candidate.multiply(numerator, inverse_denominator) ||
        !check.multiply(candidate, candidate)) {
        return false;
    }
    return check.equal(input);
}

bool try_full_degree_pure_square_root_at_prime(
        Element& root,
        const Element& input,
        const flint::FmpqPoly& element_polynomial,
        const nf_struct* raw_field,
        const flint::Fmpz& prime,
        const flint::FmpzModPoly& factor,
        const flint::FmpzModCtx& prime_ctx,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_full_degree_square");
    const NumberField* parent = input.parent();
    if (parent == nullptr || factor.context() == nullptr ||
        fmpz_mod_poly_degree(factor.raw(), prime_ctx.raw()) !=
                parent->degree() ||
        !fmpq_poly_is_monic_integral(raw_field->pol)) {
        return false;
    }

    flint::FmpzModPoly element_mod_prime(prime_ctx);
    flint::FmpzModPoly sqrt_mod_prime(prime_ctx);
    flint::FmpzModPoly inverse_root(prime_ctx);
    if (!element_mod_prime.is_initialized() ||
        !sqrt_mod_prime.is_initialized() || !inverse_root.is_initialized() ||
        !fmpq_poly_reduce_modulus(element_mod_prime, element_polynomial.raw(),
                                  prime_ctx, prime.raw()) ||
        !residue_field_square_root(sqrt_mod_prime, element_mod_prime, factor,
                                   prime_ctx) ||
        fmpz_mod_poly_invmod(inverse_root.raw(), sqrt_mod_prime.raw(),
                             factor.raw(), prime_ctx.raw()) == 0) {
        return false;
    }

    flint::FmpqPoly derivative_polynomial;
    fmpq_poly_derivative(derivative_polynomial.raw(), raw_field->pol);

    slong target_lift_exponent = 0;
    if (!pure_power_lifting_exponent(
                target_lift_exponent, input, derivative_polynomial, prime,
                fmpz_mod_poly_degree(factor.raw(), prime_ctx.raw()), 2,
                diagnostics)) {
        return false;
    }

    std::vector<slong> lift_chain;
    if (!lifting_chain(lift_chain, target_lift_exponent)) {
        return false;
    }

    for (slong exponent : lift_chain) {
        if (exponent <= 1) {
            continue;
        }
        flint::Fmpz precision_modulus;
        fmpz_pow_ui(precision_modulus.raw(), prime.raw(),
                    static_cast<ulong>(exponent));
        flint::FmpzModCtx ctx(precision_modulus.raw());
        if (ctx.raw() == nullptr) {
            return false;
        }

        flint::FmpzModPoly lifted_factor(ctx);
        flint::FmpzModPoly radicand(ctx);
        flint::FmpzModPoly inverse_root_lifted(ctx);
        flint::FmpzModPoly next_inverse_root(ctx);
        flint::FmpzModPoly derivative(ctx);
        flint::FmpzModPoly inverse_root_times_derivative(ctx);
        flint::FmpzModPoly scaled_root(ctx);
        if (!lifted_factor.is_initialized() || !radicand.is_initialized() ||
            !inverse_root_lifted.is_initialized() ||
            !next_inverse_root.is_initialized() ||
            !derivative.is_initialized() ||
            !inverse_root_times_derivative.is_initialized() ||
            !scaled_root.is_initialized() ||
            !fmpq_poly_reduce_modulus(lifted_factor, raw_field->pol, ctx,
                                      precision_modulus.raw()) ||
            !fmpq_poly_reduce_modulus(radicand, element_polynomial.raw(), ctx,
                                      precision_modulus.raw()) ||
            !fmpq_poly_reduce_modulus(derivative,
                                      derivative_polynomial.raw(), ctx,
                                      precision_modulus.raw()) ||
            !fmpz_mod_poly_change_modulus(inverse_root_lifted, inverse_root,
                                          ctx)) {
            return false;
        }
        fmpz_mod_poly_rem(radicand.raw(), radicand.raw(), lifted_factor.raw(),
                          ctx.raw());
        fmpz_mod_poly_rem(derivative.raw(), derivative.raw(),
                          lifted_factor.raw(), ctx.raw());

        if (!pure_square_hensel_step(
                    next_inverse_root, inverse_root_lifted, radicand,
                    lifted_factor, ctx, precision_modulus.raw())) {
            return false;
        }

        fmpz_mod_poly_mulmod(inverse_root_times_derivative.raw(),
                             next_inverse_root.raw(), derivative.raw(),
                             lifted_factor.raw(), ctx.raw());
        fmpz_mod_poly_mulmod(scaled_root.raw(),
                             inverse_root_times_derivative.raw(),
                             radicand.raw(), lifted_factor.raw(), ctx.raw());

        Element candidate(*parent);
        if (!candidate.is_defined()) {
            return false;
        }
        if (scaled_root_candidate(candidate, input, scaled_root, ctx,
                                        precision_modulus.raw(),
                                        derivative_polynomial)) {
            root.swap(candidate);
            return true;
        }

        inverse_root = std::move(next_inverse_root);
    }

    return false;
}

enum class ResiduePowerStatus {
    unsupported,
    power,
    nonpower
};

struct ResiduePowerData {
    flint::Fmpz group_order;
    flint::Fmpz exponent;
    flint::Fmpz gcd;
    flint::Fmpz membership_exponent;
    bool may_be_power = false;
};

// True when the rational prime p divides the exponent n > 0.  Hensel lifting
// of y^n - a must not use such primes, where y^n - a is inseparable modulo p;
// the residue-field membership disproof remains valid there.
bool prime_divides_exponent(const fmpz_t prime, slong exponent) noexcept {
    if (exponent <= 0 || fmpz_sgn(prime) <= 0 ||
        fmpz_abs_fits_ui(prime) == 0) {
        return false;
    }
    return static_cast<ulong>(exponent) % fmpz_get_ui(prime) == 0;
}

bool fq_power_data(ResiduePowerData& out,
                        const flint::Fq& value,
                        const flint::FqCtx& field,
                        slong exponent) noexcept {
    out.may_be_power = false;
    if (exponent <= 1) {
        return false;
    }

    flint::Fmpz field_order;
    flint::fq_ctx_order(flint::FmpzRef(field_order), field);
    flint::fmpz_sub_ui(flint::FmpzRef(out.group_order),
                       flint::FmpzConstRef(field_order), 1);
    flint::fmpz_set_si(flint::FmpzRef(out.exponent), exponent);
    flint::fmpz_gcd(flint::FmpzRef(out.gcd),
                    flint::FmpzConstRef(out.group_order),
                    flint::FmpzConstRef(out.exponent));
    flint::fmpz_divexact(flint::FmpzRef(out.membership_exponent),
                         flint::FmpzConstRef(out.group_order),
                         flint::FmpzConstRef(out.gcd));

    flint::Fq test_value(field);
    if (!test_value.is_initialized()) {
        return false;
    }
    flint::fq_pow(test_value, value,
                  flint::FmpzConstRef(out.membership_exponent), field);
    out.may_be_power = flint::fq_is_one(test_value, field);
    return true;
}

bool append_lifted_inverse_root(std::vector<flint::FmpzPoly>& inverse_roots,
                                const flint::Fq& inverse_root,
                                const flint::FqCtx& field,
                                const flint::FmpzModPoly& modulus,
                                const flint::FmpzModCtx& ctx) noexcept {
    flint::FmpzModPoly inverse_root_mod(ctx);
    if (!inverse_root_mod.is_initialized()) {
        return false;
    }
    flint::fq_get_fmpz_mod_poly(inverse_root_mod, inverse_root, field);
    fmpz_mod_poly_rem(inverse_root_mod.raw(), inverse_root_mod.raw(),
                      modulus.raw(), ctx.raw());

    flint::FmpzPoly lifted_inverse;
    fmpz_mod_poly_get_fmpz_poly(lifted_inverse.raw(), inverse_root_mod.raw(),
                                ctx.raw());
    inverse_roots.push_back(std::move(lifted_inverse));
    return true;
}

bool fq_unique_power_inverse_root(
        std::vector<flint::FmpzPoly>& inverse_roots,
        const flint::Fq& value,
        const flint::FqCtx& field,
        const ResiduePowerData& power_data,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!flint::fmpz_is_one(flint::FmpzConstRef(power_data.gcd))) {
        return false;
    }

    flint::Fmpz root_exponent;
    if (!flint::fmpz_invmod(flint::FmpzRef(root_exponent),
                            flint::FmpzConstRef(power_data.exponent),
                            flint::FmpzConstRef(power_data.group_order))) {
        return false;
    }

    flint::Fq root(field);
    flint::Fq check(field);
    flint::Fq root_power(field);
    flint::Fq inverse_root(field);
    if (!root.is_initialized() || !check.is_initialized() ||
        !root_power.is_initialized() || !inverse_root.is_initialized()) {
        return false;
    }

    // reference gen_Shanks_sqrtn uses the Bezout coefficient of n modulo
    // #F_q - 1 directly when gcd(n, #F_q - 1) = 1.
    flint::fq_pow(root, value, flint::FmpzConstRef(root_exponent), field);
    flint::fq_pow(check, root, flint::FmpzConstRef(power_data.exponent),
                  field);
    if (!flint::fq_equal(check, value, field)) {
        return false;
    }

    flint::Fmpz exponent_minus_one;
    flint::fmpz_sub_ui(flint::FmpzRef(exponent_minus_one),
                       flint::FmpzConstRef(power_data.exponent), 1);
    flint::fq_pow(root_power, root, flint::FmpzConstRef(exponent_minus_one),
                  field);
    if (flint::fq_is_zero(root_power, field)) {
        return false;
    }
    flint::fq_inv(inverse_root, root_power, field);
    if (!append_lifted_inverse_root(inverse_roots, inverse_root, field,
                                    modulus, ctx)) {
        return false;
    }

    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::element,
            "element.power_residue_inverse_roots.unique_root");
    return true;
}

bool fmpq_poly_to_integral_fmpz_poly(flint::FmpzPoly& out,
                                     const fmpq_poly_t input) noexcept {
    if (!fmpq_poly_is_monic_integral(input)) {
        return false;
    }

    fmpz_poly_zero(out.raw());
    flint::Fmpq coeff;
    const slong degree = fmpq_poly_degree(input);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), input, i);
        fmpz_poly_set_coeff_fmpz(out.raw(), i, fmpq_numref(coeff.raw()));
    }
    return true;
}

ResiduePowerStatus residue_field_power_inverse_roots(
        std::vector<flint::FmpzPoly>& inverse_roots,
        const flint::FmpzModPoly& element,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx,
        const fmpz_t prime,
        slong exponent,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_residue_inverse_roots");
    inverse_roots.clear();
    if (exponent <= 1 ||
        fmpz_mod_poly_degree(modulus.raw(), ctx.raw()) <= 0) {
        return ResiduePowerStatus::unsupported;
    }

    flint::FmpzModPoly reduced(ctx);
    if (!reduced.is_initialized()) {
        return ResiduePowerStatus::unsupported;
    }
    fmpz_mod_poly_rem(reduced.raw(), element.raw(), modulus.raw(), ctx.raw());
    if (fmpz_mod_poly_is_zero(reduced.raw(), ctx.raw()) != 0) {
        // reference's Hensel root path skips singular pure roots.
        return ResiduePowerStatus::unsupported;
    }

    flint::FqCtx field(modulus.raw(), ctx.raw(), "z");
    flint::Fq value(field);
    flint::Fq negative_value(field);
    flint::Fq one(field);
    if (!field.is_initialized() || !value.is_initialized() ||
        !negative_value.is_initialized() || !one.is_initialized()) {
        return ResiduePowerStatus::unsupported;
    }

    fq_set_fmpz_mod_poly(value.raw(), reduced.raw(), field.raw());
    // reference Fq_ispower first tests a^((#F_q - 1)/gcd(#F_q - 1, n)) == 1.
    // Use that as an exact disproof, then use the unique-root branch from
    // reference gen_Shanks_sqrtn when gcd(n, #F_q - 1) = 1, and otherwise
    // its Bezout reduction to a gcd(n, #F_q - 1)-th root below.
    ResiduePowerData power_data;
    if (!fq_power_data(power_data, value, field, exponent)) {
        return ResiduePowerStatus::unsupported;
    }
    if (!power_data.may_be_power) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.power_residue_inverse_roots.fq_ispower_disproof");
        return ResiduePowerStatus::nonpower;
    }
    // The disproof above holds at every good prime, including p | n.  Residue
    // roots are only needed for Hensel lifting, which requires p not dividing
    // n (y^n - a is inseparable modulo p otherwise).  Report unsupported so
    // the prime loops move on to the next factor or prime.
    if (prime_divides_exponent(prime, exponent)) {
        return ResiduePowerStatus::unsupported;
    }
    if (fq_unique_power_inverse_root(
                inverse_roots, value, field, power_data, modulus, ctx,
                diagnostics)) {
        return ResiduePowerStatus::power;
    }

    // Reduce y^n = a to y^g = a^u with g = gcd(n, #F_q - 1) and
    // n*u + (#F_q - 1)*v = g, following the Bezout step of reference
    // gen_Shanks_sqrtn.  Let m = #F_q - 1.  Since a^(m/g) = 1 was checked
    // above and F_q^* is cyclic, a = w^g for some w.  If y^n = a then
    // y^g = y^(n*u + m*v) = a^u.  Conversely, if y^g = a^u then
    // y^n = (a^u)^(n/g) = w^(n*u) = w^(g - m*v) = w^g = a.  So both
    // polynomials have the same roots in F_q, and the dense polynomial has
    // degree g <= n instead of n.
    flint::Fmpz bezout_exponent;
    flint::Fmpz bezout_group_coefficient;
    flint::Fmpz bezout_gcd;
    fmpz_xgcd(bezout_gcd.raw(), bezout_exponent.raw(),
              bezout_group_coefficient.raw(), power_data.exponent.raw(),
              power_data.group_order.raw());
    if (!fmpz_equal(bezout_gcd.raw(), power_data.gcd.raw()) ||
        fmpz_fits_si(bezout_gcd.raw()) == 0) {
        return ResiduePowerStatus::unsupported;
    }
    // Resource bound only: y^g - a^u has exactly g roots in F_q, and every
    // root is lifted by the caller.  Exceeding it leaves the caller in the
    // unsupported/fail-closed state.
    constexpr slong max_reduced_root_degree = slong{1} << 16;
    const slong reduced_degree = fmpz_get_si(bezout_gcd.raw());
    if (reduced_degree < 1 || reduced_degree > max_reduced_root_degree) {
        return ResiduePowerStatus::unsupported;
    }
    fmpz_mod(bezout_exponent.raw(), bezout_exponent.raw(),
             power_data.group_order.raw());

    flint::Fq reduced_value(field);
    if (!reduced_value.is_initialized()) {
        return ResiduePowerStatus::unsupported;
    }
    flint::fq_pow(reduced_value, value, flint::FmpzConstRef(bezout_exponent),
                  field);

    fq_neg(negative_value.raw(), reduced_value.raw(), field.raw());
    fq_one(one.raw(), field.raw());

    FqPoly polynomial(field.raw());
    fq_poly_set_coeff(polynomial.raw(), 0, negative_value.raw(), field.raw());
    fq_poly_set_coeff(polynomial.raw(), reduced_degree, one.raw(),
                      field.raw());

    FqPolyFactor roots(field.raw());
    fq_poly_roots(roots.raw(), polynomial.raw(), 0, field.raw());
    if (roots.raw()->num == 0) {
        // Unreachable: y^g - a^u has g distinct roots in F_q once the
        // membership test passed.  Never turn an anomaly into a disproof.
        return ResiduePowerStatus::unsupported;
    }

    for (slong i = 0; i < roots.raw()->num; ++i) {
        const fq_poly_struct* linear = roots.raw()->poly + i;
        if (fq_poly_degree(linear, field.raw()) != 1) {
            return ResiduePowerStatus::unsupported;
        }

        flint::Fq constant(field);
        flint::Fq leading(field);
        flint::Fq inverse_leading(field);
        flint::Fq root(field);
        flint::Fq root_power(field);
        flint::Fq inverse_root(field);
        if (!constant.is_initialized() || !leading.is_initialized() ||
            !inverse_leading.is_initialized() || !root.is_initialized() ||
            !root_power.is_initialized() || !inverse_root.is_initialized()) {
            return ResiduePowerStatus::unsupported;
        }

        fq_poly_get_coeff(constant.raw(), linear, 0, field.raw());
        fq_poly_get_coeff(leading.raw(), linear, 1, field.raw());
        if (fq_is_zero(leading.raw(), field.raw()) != 0) {
            return ResiduePowerStatus::unsupported;
        }
        fq_inv(inverse_leading.raw(), leading.raw(), field.raw());
        fq_neg(root.raw(), constant.raw(), field.raw());
        fq_mul(root.raw(), root.raw(), inverse_leading.raw(), field.raw());
        if (fq_is_zero(root.raw(), field.raw()) != 0) {
            return ResiduePowerStatus::unsupported;
        }
        // Check the reduction above: every root must satisfy y^n = a.
        flint::fq_pow(root_power, root,
                      flint::FmpzConstRef(power_data.exponent), field);
        if (!flint::fq_equal(root_power, value, field)) {
            return ResiduePowerStatus::unsupported;
        }

        fq_pow_ui(root_power.raw(), root.raw(),
                  static_cast<ulong>(exponent - 1), field.raw());
        if (fq_is_zero(root_power.raw(), field.raw()) != 0) {
            return ResiduePowerStatus::unsupported;
        }
        fq_inv(inverse_root.raw(), root_power.raw(), field.raw());

        if (!append_lifted_inverse_root(inverse_roots, inverse_root, field,
                                        modulus, ctx)) {
            return ResiduePowerStatus::unsupported;
        }
    }

    return inverse_roots.empty() ? ResiduePowerStatus::unsupported
                                 : ResiduePowerStatus::power;
}

bool hensel_lift_factor_to_precision(
        flint::FmpzModPoly& out,
        const fmpq_poly_t field_polynomial,
        const flint::FmpzModPoly& field_mod_prime,
        const flint::FmpzModPoly& factor_mod_prime,
        const flint::FmpzModCtx& prime_ctx,
        const fmpz_t prime,
        slong target_exponent,
        const fmpz_t target_modulus,
        const flint::FmpzModCtx& target_ctx) noexcept {
    const slong degree =
            fmpz_mod_poly_degree(field_mod_prime.raw(), prime_ctx.raw());
    if (fmpz_mod_poly_degree(factor_mod_prime.raw(), prime_ctx.raw()) ==
        degree) {
        return fmpq_poly_reduce_modulus(out, field_polynomial, target_ctx,
                                        target_modulus);
    }
    if (target_exponent <= 0) {
        return false;
    }

    flint::FmpzPoly field_integral;
    if (!fmpq_poly_to_integral_fmpz_poly(field_integral, field_polynomial)) {
        return false;
    }

    flint::FmpzModPoly quotient(prime_ctx);
    flint::FmpzModPoly remainder(prime_ctx);
    flint::FmpzModPoly gcd(prime_ctx);
    flint::FmpzModPoly a_mod(prime_ctx);
    flint::FmpzModPoly b_mod(prime_ctx);
    if (!quotient.is_initialized() || !remainder.is_initialized() ||
        !gcd.is_initialized() || !a_mod.is_initialized() ||
        !b_mod.is_initialized()) {
        return false;
    }

    fmpz_mod_poly_divrem(quotient.raw(), remainder.raw(),
                         field_mod_prime.raw(), factor_mod_prime.raw(),
                         prime_ctx.raw());
    if (fmpz_mod_poly_is_zero(remainder.raw(), prime_ctx.raw()) == 0) {
        return false;
    }
    fmpz_mod_poly_xgcd(gcd.raw(), a_mod.raw(), b_mod.raw(),
                       factor_mod_prime.raw(), quotient.raw(),
                       prime_ctx.raw());
    if (fmpz_mod_poly_is_one(gcd.raw(), prime_ctx.raw()) == 0) {
        return false;
    }

    flint::FmpzPoly g;
    flint::FmpzPoly h;
    flint::FmpzPoly a;
    flint::FmpzPoly b;
    fmpz_mod_poly_get_fmpz_poly(g.raw(), factor_mod_prime.raw(),
                                prime_ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(h.raw(), quotient.raw(), prime_ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(a.raw(), a_mod.raw(), prime_ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(b.raw(), b_mod.raw(), prime_ctx.raw());

    flint::Fmpz current_modulus;
    fmpz_set(current_modulus.raw(), prime);
    slong current_exponent = 1;
    while (current_exponent < target_exponent) {
        const slong step =
                std::min(current_exponent,
                         target_exponent - current_exponent);
        flint::Fmpz next_factor;
        fmpz_pow_ui(next_factor.raw(), prime, static_cast<ulong>(step));

        flint::FmpzPoly g_out;
        flint::FmpzPoly h_out;
        flint::FmpzPoly a_out;
        flint::FmpzPoly b_out;
        fmpz_poly_hensel_lift(g_out.raw(), h_out.raw(), a_out.raw(),
                              b_out.raw(), field_integral.raw(), g.raw(),
                              h.raw(), a.raw(), b.raw(),
                              current_modulus.raw(), next_factor.raw());
        g = std::move(g_out);
        h = std::move(h_out);
        a = std::move(a_out);
        b = std::move(b_out);
        fmpz_mul(current_modulus.raw(), current_modulus.raw(),
                 next_factor.raw());
        current_exponent += step;
    }

    fmpz_poly_scalar_mod_fmpz(g.raw(), g.raw(), target_modulus);
    fmpz_mod_poly_set_fmpz_poly(out.raw(), g.raw(), target_ctx.raw());
    return true;
}

bool get_basis(flint::FmpzMat& out,
                     const fmpz_t prime_power,
                     slong degree,
                     const flint::FmpzModPoly& lifted_factor,
                     const flint::FmpzModCtx& ctx) noexcept {
    if (flint::fmpz_mat_nrows(out) != degree ||
        flint::fmpz_mat_ncols(out) != degree) {
        return false;
    }

    const slong factor_degree =
            fmpz_mod_poly_degree(lifted_factor.raw(), ctx.raw());
    if (factor_degree <= 0 || factor_degree > degree) {
        return false;
    }

    fmpz_mat_zero(out.raw());
    for (slong j = 0; j < factor_degree; ++j) {
        fmpz_set(fmpz_mat_entry(out.raw(), j, j), prime_power);
    }

    if (factor_degree == degree) {
        return true;
    }

    flint::FmpzModPoly x(ctx);
    flint::FmpzModPoly power(ctx);
    flint::FmpzModPoly next(ctx);
    if (!x.is_initialized() || !power.is_initialized() ||
        !next.is_initialized()) {
        return false;
    }
    fmpz_mod_poly_set_coeff_ui(x.raw(), 1, 1, ctx.raw());
    fmpz_mod_poly_set_coeff_ui(power.raw(), factor_degree - 1, 1,
                               ctx.raw());

    flint::Fmpz coeff;
    for (slong row = factor_degree; row < degree; ++row) {
        fmpz_mod_poly_mulmod(next.raw(), power.raw(), x.raw(),
                             lifted_factor.raw(), ctx.raw());
        power = std::move(next);
        next = flint::FmpzModPoly(ctx);
        if (!next.is_initialized()) {
            return false;
        }

        fmpz_one(fmpz_mat_entry(out.raw(), row, row));
        const slong power_degree =
                fmpz_mod_poly_degree(power.raw(), ctx.raw());
        for (slong col = 0; col <= power_degree; ++col) {
            fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), power.raw(), col,
                                         ctx.raw());
            fmpz_neg(fmpz_mat_entry(out.raw(), row, col), coeff.raw());
        }
    }

    return true;
}

void round_nearest_ties_away(fmpz_t out,
                             const fmpz_t numerator,
                             const fmpz_t denominator) noexcept {
    flint::Fmpz den;
    flint::Fmpz num_abs;
    fmpz_set(den.raw(), denominator);
    fmpz_abs(num_abs.raw(), numerator);
    slong sign = fmpz_sgn(numerator);
    if (fmpz_sgn(den.raw()) < 0) {
        fmpz_neg(den.raw(), den.raw());
        sign = -sign;
    }

    flint::Fmpz doubled_num;
    flint::Fmpz doubled_den;
    fmpz_mul_2exp(doubled_num.raw(), num_abs.raw(), 1);
    fmpz_add(doubled_num.raw(), doubled_num.raw(), den.raw());
    fmpz_mul_2exp(doubled_den.raw(), den.raw(), 1);
    fmpz_fdiv_q(out, doubled_num.raw(), doubled_den.raw());
    if (sign < 0) {
        fmpz_neg(out, out);
    }
}

bool reconstruct_power_root_candidate(
        Element& candidate,
        const Element& input,
        const flint::FmpzModPoly& scaled_root,
        const flint::FmpzModCtx& ctx,
        const flint::FmpzMat& basis,
        const flint::FmpzMat& inverse_num,
        flint::FmpzConstRef inverse_den,
        const Element& inverse_denominator,
        flint::FmpzConstRef exponent,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_reconstruct");
    const NumberField* parent = input.parent();
    if (parent == nullptr || !candidate.has_parent(*parent) ||
        !inverse_denominator.has_parent(*parent) ||
        fmpz_sgn(exponent.raw()) <= 0 ||
        flint::fmpz_mat_nrows(basis) != parent->degree() ||
        flint::fmpz_mat_ncols(basis) != parent->degree() ||
        flint::fmpz_mat_nrows(inverse_num) != parent->degree() ||
        flint::fmpz_mat_ncols(inverse_num) != parent->degree()) {
        return false;
    }

    const slong degree = parent->degree();
    flint::FmpzMat ve(1, degree);
    flint::Fmpz coeff;
    for (slong col = 0; col < degree; ++col) {
        fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), scaled_root.raw(), col,
                                     ctx.raw());
        fmpz_set(fmpz_mat_entry(ve.raw(), 0, col), coeff.raw());
    }

    flint::FmpzMat scaled_coordinates(1, degree);
    flint::FmpzMat mu(1, degree);
    flint::FmpzMat correction(1, degree);
    fmpz_mat_mul(scaled_coordinates.raw(), ve.raw(), inverse_num.raw());
    for (slong col = 0; col < degree; ++col) {
        round_nearest_ties_away(
                fmpz_mat_entry(mu.raw(), 0, col),
                fmpz_mat_entry(scaled_coordinates.raw(), 0, col),
                inverse_den.raw());
    }
    fmpz_mat_mul(correction.raw(), mu.raw(), basis.raw());
    for (slong col = 0; col < degree; ++col) {
        fmpz_sub(fmpz_mat_entry(ve.raw(), 0, col),
                 fmpz_mat_entry(ve.raw(), 0, col),
                 fmpz_mat_entry(correction.raw(), 0, col));
    }

    flint::FmpqPoly numerator_polynomial;
    fmpq_poly_zero(numerator_polynomial.raw());
    for (slong col = 0; col < degree; ++col) {
        fmpq_poly_set_coeff_fmpz(
                numerator_polynomial.raw(), col,
                fmpz_mat_entry(ve.raw(), 0, col));
    }

    Element numerator(*parent);
    if (!numerator.is_defined() ||
        !numerator.set_fmpq_poly(
                flint::FmpqPolyConstRef(numerator_polynomial)) ||
        !candidate.multiply(numerator, inverse_denominator)) {
        return false;
    }
    return true;
}

// Necessary condition for candidate^n == input from multiplicativity of the
// norm: N(candidate)^n == N(input).  It only rejects candidates, so a
// rejection leaves the caller in its unsupported state.  When
// N(candidate) = r/s in lowest terms has max(|r|, s) >= 2, r^n/s^n is in
// lowest terms and max(|r|^n, s^n) >= 2^n, so N(input) must have a numerator
// or denominator of at least n + 1 bits.  This avoids forming candidate^n for
// huge n unless the input is itself that large.  A candidate of norm 0 or
// +-1 passes after the cheap sign/zero comparison; unit_norm reports the
// +-1 case, which the norm cannot bound.
bool norm_admits_power_root(bool& unit_norm,
                            const Element& candidate,
                            const Element& input,
                            flint::FmpzConstRef exponent) noexcept {
    unit_norm = false;
    flint::Fmpq candidate_norm;
    flint::Fmpq input_norm;
    if (!candidate.norm(flint::FmpqRef(candidate_norm)) ||
        !input.norm(flint::FmpqRef(input_norm))) {
        return false;
    }

    const fmpz* candidate_num = fmpq_numref(candidate_norm.raw());
    const fmpz* candidate_den = fmpq_denref(candidate_norm.raw());
    const fmpz* input_num = fmpq_numref(input_norm.raw());
    const fmpz* input_den = fmpq_denref(input_norm.raw());
    if (fmpz_is_zero(candidate_num) != 0) {
        return fmpz_is_zero(input_num) != 0;
    }
    if (fmpz_is_one(candidate_den) != 0 && fmpz_is_pm1(candidate_num) != 0) {
        unit_norm = true;
        if (fmpz_is_one(input_den) == 0 || fmpz_is_pm1(input_num) == 0) {
            return false;
        }
        const bool negative =
                fmpz_sgn(candidate_num) < 0 && fmpz_is_odd(exponent.raw()) != 0;
        return (fmpz_sgn(input_num) < 0) == negative;
    }

    const flint_bitcnt_t input_bits =
            std::max(fmpz_bits(input_num), fmpz_bits(input_den));
    if (fmpz_cmp_ui(exponent.raw(), input_bits) >= 0) {
        return false;
    }
    const ulong power = fmpz_get_ui(exponent.raw());
    flint::Fmpz powered;
    fmpz_pow_ui(powered.raw(), candidate_den, power);
    if (!fmpz_equal(powered.raw(), input_den)) {
        return false;
    }
    fmpz_pow_ui(powered.raw(), candidate_num, power);
    return fmpz_equal(powered.raw(), input_num) != 0;
}

flint_bitcnt_t fmpq_poly_coefficient_bits(
        const fmpq_poly_struct* polynomial) noexcept {
    flint_bitcnt_t bits = fmpz_bits(polynomial->den);
    const slong length = fmpq_poly_length(polynomial);
    for (slong i = 0; i < length; ++i) {
        bits = std::max(bits, fmpz_bits(polynomial->coeffs + i));
    }
    return bits;
}

flint_bitcnt_t element_coefficient_bits(const Element& element) noexcept {
    flint::FmpqPoly polynomial;
    if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return 0;
    }
    return fmpq_poly_coefficient_bits(polynomial.raw());
}

// Leading coefficient of the primitive integral characteristic polynomial
// P of element over Q.  It is 1 exactly when element is an algebraic integer,
// and log lc(P) is the finite-place part of log M(P) = d h(element): it is
// the sum over finite places v of d_v log max(1, |element|_v) (Bombieri and
// Gubler, Heights in Diophantine Geometry, sections 1.5-1.6).
bool primitive_characteristic_leading_coefficient(
        flint::Fmpz& out,
        const Element& element) noexcept {
    const NumberField* parent = element.parent();
    const nf_struct* raw_field =
            parent == nullptr ? nullptr : parent->raw_flint_field();
    const slong degree = parent == nullptr ? 0 : parent->degree();
    if (raw_field == nullptr || degree <= 0) {
        return false;
    }

    // The characteristic polynomial chi is monic in Q[x]; FLINT stores it as
    // an integral numerator over den with gcd(content(numerator), den) = 1,
    // so lc of the primitive part is den / content(numerator).
    flint::FmpqMat matrix(degree, degree);
    flint::FmpqPoly characteristic;
    multiplication_matrix(matrix.raw(), raw_field, degree,
                          element.raw_flint_element());
    fmpq_mat_charpoly(characteristic.raw(), matrix.raw());
    if (fmpq_poly_degree(characteristic.raw()) != degree) {
        return false;
    }
    flint::Fmpz content;
    _fmpz_vec_content(content.raw(), characteristic.raw()->coeffs,
                      degree + 1);
    if (fmpz_is_zero(content.raw()) != 0) {
        return false;
    }
    fmpz_divexact(out.raw(), characteristic.raw()->den, content.raw());
    return true;
}

// Enclosure of log M(P), where P is the primitive integral characteristic
// polynomial of element over Q and M is the Mahler measure:
//   log M(P) = log |lc(P)| + sum_{i=1}^{d} log max(1, |sigma_i(element)|),
// with sigma_1, ..., sigma_d the complex embeddings of the parent field
// Since P is a power of the primitive minimal polynomial, log M(P) =
// d * h(element) for the absolute logarithmic Weil height h (Bombieri and
// Gubler, Heights in Diophantine Geometry, sections 1.5-1.6).
bool log_mahler_measure_enclosure(flint::Arb& out,
                                  EmbeddingContext& embeddings,
                                  const Element& element,
                                  slong precision) noexcept {
    const NumberField* parent = element.parent();
    if (parent == nullptr) {
        return false;
    }
    const slong degree = parent->degree();
    flint::Fmpz leading;
    if (!primitive_characteristic_leading_coefficient(leading, element)) {
        return false;
    }

    flint::AcbVec values(degree);
    if (!embeddings.evaluate_all(flint::AcbVecRef(values), element,
                                 precision)) {
        return false;
    }

    flint::Arb one;
    flint::Arb absolute;
    flint::Arb term;
    flint::arb_one(one);
    flint::arb_log_fmpz(out, flint::FmpzConstRef(leading), precision);
    for (slong i = 0; i < degree; ++i) {
        flint::acb_abs(absolute, values.data() + i, precision);
        arb_max(term.raw(), absolute.raw(), one.raw(), precision);
        flint::arb_log(term, term, precision);
        flint::arb_add(out, out, term, precision);
    }
    return flint::arb_is_finite(out);
}

// Necessary condition for candidate^n == input from the height: if
// c^n = a then h(a) = n h(c) (ibid., section 1.5), so with both
// characteristic polynomials of degree d = [K : Q],
// n log M(P_c) = log M(P_a).  The candidate is rejected only when certified
// enclosures prove n log M(P_c) > log M(P_a) + 1, so a true root is never
// rejected.  When the test passes, h(c^n) <= h(a) + (1 + r)/d, with r the
// total width of the two enclosures, so forming c^n costs about as much as
// the input's own size; a root of unity c (M(P_c) = 1) always passes and
// binary powering keeps its coefficients bounded.  A failure to compute an
// enclosure, or a comparison the enclosures cannot decide, leaves the
// candidate to exact verification.  The working precision grows with the
// coefficient sizes of c, a and the defining polynomial; it affects only
// whether the comparison is decisive, never whether a rejection is sound.
bool height_admits_power_root(const Element& candidate,
                              const Element& input,
                              flint::FmpzConstRef exponent) noexcept {
    const NumberField* parent = input.parent();
    const nf_struct* raw_field =
            parent == nullptr ? nullptr : parent->raw_flint_field();
    if (raw_field == nullptr) {
        return true;
    }
    const slong precision =
            128 + static_cast<slong>(std::max(
                          fmpq_poly_coefficient_bits(raw_field->pol),
                          std::max(element_coefficient_bits(candidate),
                                   element_coefficient_bits(input))));
    EmbeddingContext embeddings(*parent);
    flint::Arb candidate_log;
    flint::Arb input_log;
    if (!embeddings.refine(precision) ||
        !log_mahler_measure_enclosure(candidate_log, embeddings, candidate,
                                      precision) ||
        !log_mahler_measure_enclosure(input_log, embeddings, input,
                                      precision)) {
        return true;
    }
    flint::arb_mul_fmpz(candidate_log, candidate_log, exponent, precision);
    flint::arb_add_ui(input_log, input_log, 1, precision);
    return arb_gt(candidate_log.raw(), input_log.raw()) == 0;
}

bool verify_power_root(const Element& candidate,
                             const Element& input,
                             flint::FmpzConstRef exponent,
                             const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_verify");
    const NumberField* parent = input.parent();
    if (parent == nullptr || !candidate.has_parent(*parent) ||
        fmpz_sgn(exponent.raw()) <= 0) {
        return false;
    }

    bool unit_norm = false;
    if (!norm_admits_power_root(unit_norm, candidate, input, exponent)) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_verify.norm_reject");
        return false;
    }
    if (unit_norm) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_verify_height");
        if (!height_admits_power_root(candidate, input, exponent)) {
            SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                                "element.power_hensel_verify.height_reject");
            return false;
        }
    }

    Element check(*parent);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_verify_pow");
        if (!check.is_defined() || !check.pow_fmpz(candidate, exponent)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_verify_equal");
        return check.equal(input);
    }
}

// Key: the degree, the modulus p^k (which is also the modulus of `ctx`, see
// power_reconstruction_data) and the lifted factor's coefficients reduced
// mod p^k, which are all that get_basis, fmpz_lll and fmpz_mat_inv read.
// The lifted factor is the Hensel lift of a factor of f mod p to p^k, so it
// depends on (field, prime, factor, k) and never on the radicand.
struct ReconstructionDataEntry {
    slong degree = 0;
    flint::Fmpz modulus;
    flint::FmpzPoly lifted_factor;
    flint::FmpzMat basis{0, 0};
    flint::FmpzMat inverse_num{0, 0};
    flint::Fmpz inverse_den;
};

// One field, one prime, its factors and the lift exponents of the chains in
// use; see the cache comment above equation_order_norm_change_constants.
constexpr std::size_t kReconstructionDataCacheCapacity = 256;

BoundedEntryCache<ReconstructionDataEntry>&
reconstruction_data_cache() noexcept {
    thread_local BoundedEntryCache<ReconstructionDataEntry> cache(
            kReconstructionDataCacheCapacity);
    return cache;
}

bool compute_power_reconstruction_data(
        flint::FmpzMat& basis,
        flint::FmpzMat& inverse_num,
        flint::Fmpz& inverse_den,
        slong degree,
        const flint::FmpzModPoly& lifted_factor,
        const flint::FmpzModCtx& ctx,
        const fmpz_t precision_modulus) noexcept {
    if (!get_basis(basis, precision_modulus, degree, lifted_factor, ctx)) {
        return false;
    }

    flint::FmpzMat transform(degree, degree);
    fmpz_mat_one(transform.raw());
    flint::FmpzLll lll;
    fmpz_lll(basis.raw(), transform.raw(), lll.raw());

    return flint::fmpz_mat_inv(flint::FmpzMatRef(inverse_num),
                               flint::FmpzRef(inverse_den),
                               flint::FmpzMatConstRef(basis));
}

// LLL-reduced basis of the reconstruction lattice and its inverse, memoized
// per thread (see the cache comment above equation_order_norm_change_
// constants).  The memo is used only when `ctx` has modulus
// `precision_modulus`, as at every caller.
bool power_reconstruction_data(
        flint::FmpzMat& basis,
        flint::FmpzMat& inverse_num,
        flint::Fmpz& inverse_den,
        slong degree,
        const flint::FmpzModPoly& lifted_factor,
        const flint::FmpzModCtx& ctx,
        const fmpz_t precision_modulus,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_reconstruction_data");
    if (degree <= 0 ||
        flint::fmpz_mat_nrows(basis) != degree ||
        flint::fmpz_mat_ncols(basis) != degree ||
        flint::fmpz_mat_nrows(inverse_num) != degree ||
        flint::fmpz_mat_ncols(inverse_num) != degree) {
        return false;
    }
    if (ctx.raw() == nullptr ||
        fmpz_equal(fmpz_mod_ctx_modulus(ctx.raw()), precision_modulus) == 0) {
        return compute_power_reconstruction_data(
                basis, inverse_num, inverse_den, degree, lifted_factor, ctx,
                precision_modulus);
    }

    flint::FmpzPoly key_factor;
    fmpz_mod_poly_get_fmpz_poly(key_factor.raw(), lifted_factor.raw(),
                                ctx.raw());
    auto& cache = reconstruction_data_cache();
    for (const ReconstructionDataEntry& entry : cache.entries()) {
        if (entry.degree == degree &&
            fmpz_equal(entry.modulus.raw(), precision_modulus) != 0 &&
            fmpz_poly_equal(entry.lifted_factor.raw(), key_factor.raw()) !=
                    0) {
            fmpz_mat_set(basis.raw(), entry.basis.raw());
            fmpz_mat_set(inverse_num.raw(), entry.inverse_num.raw());
            fmpz_set(inverse_den.raw(), entry.inverse_den.raw());
            ++hensel_data_cache_counters_tls.reconstruction_data_hits;
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.power_hensel_reconstruction_data.cache_hit");
            return true;
        }
    }

    ++hensel_data_cache_counters_tls.reconstruction_data_misses;
    SILEX_PROFILE_EVENT(
            diagnostics, DiagnosticsModule::element,
            "element.power_hensel_reconstruction_data.cache_miss");
    if (!compute_power_reconstruction_data(basis, inverse_num, inverse_den,
                                           degree, lifted_factor, ctx,
                                           precision_modulus)) {
        return false;
    }
    ReconstructionDataEntry entry;
    entry.degree = degree;
    fmpz_set(entry.modulus.raw(), precision_modulus);
    entry.lifted_factor = std::move(key_factor);
    entry.basis = flint::FmpzMat(degree, degree);
    entry.inverse_num = flint::FmpzMat(degree, degree);
    fmpz_mat_set(entry.basis.raw(), basis.raw());
    fmpz_mat_set(entry.inverse_num.raw(), inverse_num.raw());
    fmpz_set(entry.inverse_den.raw(), inverse_den.raw());
    cache.insert(std::move(entry));
    return true;
}

bool pure_power_hensel_step(
        flint::FmpzModPoly& next_inverse_root,
        const flint::FmpzModPoly& inverse_root,
        const flint::FmpzModPoly& radicand,
        slong exponent,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx,
        const fmpz_t precision_modulus,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_step");
    if (exponent <= 1) {
        return false;
    }

    flint::Fmpz exponent_fmpz;
    flint::Fmpz inverse_exponent;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_step_inverse_exponent");
        fmpz_set_si(exponent_fmpz.raw(), exponent);
        if (fmpz_invmod(inverse_exponent.raw(), exponent_fmpz.raw(),
                        precision_modulus) == 0) {
            return false;
        }
    }

    flint::Fmpz one_plus_inverse_exponent;
    fmpz_add_ui(one_plus_inverse_exponent.raw(),
                inverse_exponent.raw(), 1);

    flint::FmpzModPoly bp(ctx);
    flint::FmpzModPoly power(ctx);
    flint::FmpzModPoly product(ctx);
    flint::FmpzModPoly term1(ctx);
    flint::FmpzModPoly term2(ctx);
    if (!bp.is_initialized() || !power.is_initialized() ||
        !product.is_initialized() || !term1.is_initialized() ||
        !term2.is_initialized()) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_step_powmods");
        fmpz_mod_poly_powmod_ui_binexp(
                bp.raw(), radicand.raw(), static_cast<ulong>(exponent - 1),
                modulus.raw(), ctx.raw());
        fmpz_mod_poly_scalar_mul_fmpz(term1.raw(), inverse_root.raw(),
                                      one_plus_inverse_exponent.raw(),
                                      ctx.raw());
        fmpz_mod_poly_powmod_ui_binexp(
                power.raw(), inverse_root.raw(),
                static_cast<ulong>(exponent + 1), modulus.raw(), ctx.raw());
        fmpz_mod_poly_mulmod(product.raw(), bp.raw(), power.raw(),
                             modulus.raw(), ctx.raw());
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_step_combine");
        fmpz_mod_poly_scalar_mul_fmpz(term2.raw(), product.raw(),
                                      inverse_exponent.raw(), ctx.raw());
        fmpz_mod_poly_sub(next_inverse_root.raw(), term1.raw(), term2.raw(),
                          ctx.raw());
        fmpz_mod_poly_rem(next_inverse_root.raw(), next_inverse_root.raw(),
                          modulus.raw(), ctx.raw());
    }
    return true;
}

// 1 / f'(theta) in the field, memoized per thread: try_pure_power_root_at_
// prime divides every reconstructed candidate by it and used to invert it
// for each radicand.  Key: the defining polynomial and backend, the only
// inputs of the inversion; the stored value is the inverse's power-basis
// coordinates, which determine the element exactly.
struct InverseDerivativeEntry {
    flint::FmpqPoly polynomial;
    NumberFieldBackendKind backend = NumberFieldBackendKind::generic;
    flint::FmpqPoly inverse;
};

constexpr std::size_t kInverseDerivativeCacheCapacity = 8;

BoundedEntryCache<InverseDerivativeEntry>& inverse_derivative_cache() noexcept {
    thread_local BoundedEntryCache<InverseDerivativeEntry> cache(
            kInverseDerivativeCacheCapacity);
    return cache;
}

// Sets `out` (an element of `parent`) to 1 / f'(theta), where
// `derivative_polynomial` is f' for the defining polynomial f of `parent`.
bool inverse_field_derivative(Element& out,
                              const NumberField& parent,
                              const flint::FmpqPoly& derivative_polynomial)
        noexcept {
    const nf_struct* raw_field = parent.raw_flint_field();
    if (raw_field == nullptr) {
        return false;
    }
    const NumberFieldBackendKind backend = parent.backend_kind();
    for (const InverseDerivativeEntry& entry :
         inverse_derivative_cache().entries()) {
        if (entry.backend == backend &&
            fmpq_poly_equal(entry.polynomial.raw(), raw_field->pol) != 0) {
            ++hensel_data_cache_counters_tls.inverse_derivative_hits;
            return out.set_fmpq_poly(flint::FmpqPolyConstRef(entry.inverse));
        }
    }

    ++hensel_data_cache_counters_tls.inverse_derivative_misses;
    Element denominator(parent);
    InverseDerivativeEntry entry;
    if (!denominator.is_defined() ||
        !denominator.set_fmpq_poly(
                flint::FmpqPolyConstRef(derivative_polynomial)) ||
        !out.invert(denominator) ||
        !out.get_fmpq_poly(flint::FmpqPolyRef(entry.inverse))) {
        return false;
    }
    fmpq_poly_set(entry.polynomial.raw(), raw_field->pol);
    entry.backend = backend;
    inverse_derivative_cache().insert(std::move(entry));
    return true;
}

// Lifts a root of y^n = input, where input = power_input * d^n and
// root_scale = 1/d (nullptr when power_input is used as is, i.e. d = 1), and
// verifies the scaled candidate against power_input.  See
// pure_power_hensel_root.
bool try_pure_power_root_at_prime(
        Element& root,
        const Element& input,
        const Element& power_input,
        const Element* root_scale,
        const flint::FmpqPoly& element_polynomial,
        const nf_struct* raw_field,
        const flint::Fmpz& prime,
        const flint::FmpzModPoly& field_mod_prime,
        const flint::FmpzModPoly& factor,
        const flint::FmpzModCtx& prime_ctx,
        std::vector<flint::FmpzPoly> inverse_roots,
        slong exponent,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_prime");
    const NumberField* parent = input.parent();
    if (parent == nullptr || factor.context() == nullptr || exponent <= 1 ||
        !fmpq_poly_is_monic_integral(raw_field->pol) ||
        prime_divides_exponent(prime.raw(), exponent)) {
        return false;
    }

    if (exponent == 2 &&
        fmpz_mod_poly_degree(factor.raw(), prime_ctx.raw()) ==
                parent->degree() &&
        try_full_degree_pure_square_root_at_prime(
                root, input, element_polynomial, raw_field, prime, factor,
                prime_ctx, diagnostics)) {
        // root^2 == input == power_input * d^2 exactly, so (root / d)^2 ==
        // power_input.  Skip the multiply entirely when d = 1.
        return root_scale == nullptr || root.multiply(root, *root_scale);
    }

    if (inverse_roots.empty()) {
        return false;
    }

    flint::FmpqPoly derivative_polynomial;
    Element inverse_denominator(*parent);
    flint::Fmpz exponent_fmpz;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_prime_denominator");
        fmpq_poly_derivative(derivative_polynomial.raw(), raw_field->pol);
        fmpz_set_si(exponent_fmpz.raw(), exponent);
        if (!inverse_denominator.is_defined() ||
            !inverse_field_derivative(inverse_denominator, *parent,
                                      derivative_polynomial) ||
            (root_scale != nullptr &&
             !inverse_denominator.multiply(inverse_denominator,
                                           *root_scale))) {
            return false;
        }
    }

    slong target_lift_exponent = 0;
    if (!pure_power_lifting_exponent(
                target_lift_exponent, input, derivative_polynomial, prime,
                fmpz_mod_poly_degree(factor.raw(), prime_ctx.raw()),
                exponent, diagnostics)) {
        return false;
    }

    std::vector<slong> lift_chain;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_lift_chain");
        if (!lifting_chain(lift_chain, target_lift_exponent)) {
            return false;
        }
    }

    std::vector<Element> previous_candidates;
    previous_candidates.reserve(inverse_roots.size());
    std::vector<bool> have_previous(inverse_roots.size(), false);
    for (std::size_t i = 0; i < inverse_roots.size(); ++i) {
        previous_candidates.emplace_back(*parent);
        if (!previous_candidates.back().is_defined()) {
            return false;
        }
    }

    for (slong lift_exponent : lift_chain) {
        if (lift_exponent <= 1) {
            continue;
        }
        const bool final_lift = lift_exponent == lift_chain.back();
        flint::Fmpz precision_modulus;
        fmpz_pow_ui(precision_modulus.raw(), prime.raw(),
                    static_cast<ulong>(lift_exponent));
        flint::FmpzModCtx ctx(precision_modulus.raw());
        if (ctx.raw() == nullptr) {
            return false;
        }

        flint::FmpzModPoly lifted_factor(ctx);
        flint::FmpzModPoly radicand(ctx);
        flint::FmpzModPoly derivative(ctx);
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.power_hensel_lift_factor_reduce");
            if (!lifted_factor.is_initialized() ||
                !radicand.is_initialized() || !derivative.is_initialized() ||
                !hensel_lift_factor_to_precision(
                        lifted_factor, raw_field->pol, field_mod_prime,
                        factor, prime_ctx, prime.raw(), lift_exponent,
                        precision_modulus.raw(), ctx) ||
                !fmpq_poly_reduce_modulus(radicand,
                                          element_polynomial.raw(), ctx,
                                          precision_modulus.raw()) ||
                !fmpq_poly_reduce_modulus(derivative,
                                          derivative_polynomial.raw(), ctx,
                                          precision_modulus.raw())) {
                return false;
            }
            fmpz_mod_poly_rem(radicand.raw(), radicand.raw(),
                              lifted_factor.raw(), ctx.raw());
            fmpz_mod_poly_rem(derivative.raw(), derivative.raw(),
                              lifted_factor.raw(), ctx.raw());
        }

        flint::FmpzMat reconstruction_basis(parent->degree(),
                                            parent->degree());
        flint::FmpzMat reconstruction_inverse_num(parent->degree(),
                                                  parent->degree());
        flint::Fmpz reconstruction_inverse_den;
        if (!power_reconstruction_data(
                    reconstruction_basis, reconstruction_inverse_num,
                    reconstruction_inverse_den, parent->degree(),
                    lifted_factor, ctx, precision_modulus.raw(),
                    diagnostics)) {
            return false;
        }

        for (std::size_t root_index = 0; root_index < inverse_roots.size();
             ++root_index) {
            auto& inverse_root_z = inverse_roots[root_index];
            flint::FmpzModPoly inverse_root(ctx);
            flint::FmpzModPoly next_inverse_root(ctx);
            flint::FmpzModPoly inverse_times_derivative(ctx);
            flint::FmpzModPoly scaled_root(ctx);
            if (!inverse_root.is_initialized() ||
                !next_inverse_root.is_initialized() ||
                !inverse_times_derivative.is_initialized() ||
                !scaled_root.is_initialized()) {
                return false;
            }
            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                    "element.power_hensel_root_reduce");
                fmpz_mod_poly_set_fmpz_poly(inverse_root.raw(),
                                            inverse_root_z.raw(), ctx.raw());
                fmpz_mod_poly_rem(inverse_root.raw(), inverse_root.raw(),
                                  lifted_factor.raw(), ctx.raw());
            }

            if (!pure_power_hensel_step(
                        next_inverse_root, inverse_root, radicand,
                        exponent, lifted_factor, ctx,
                        precision_modulus.raw(), diagnostics)) {
                return false;
            }

            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                    "element.power_hensel_scaled_root");
                fmpz_mod_poly_get_fmpz_poly(inverse_root_z.raw(),
                                            next_inverse_root.raw(),
                                            ctx.raw());
                fmpz_mod_poly_mulmod(inverse_times_derivative.raw(),
                                     next_inverse_root.raw(),
                                     derivative.raw(), lifted_factor.raw(),
                                     ctx.raw());
                fmpz_mod_poly_mulmod(scaled_root.raw(),
                                     inverse_times_derivative.raw(),
                                     radicand.raw(), lifted_factor.raw(),
                                     ctx.raw());
            }

            Element candidate(*parent);
            if (!candidate.is_defined()) {
                return false;
            }
            if (!reconstruct_power_root_candidate(
                        candidate, input, scaled_root, ctx,
                        reconstruction_basis, reconstruction_inverse_num,
                        flint::FmpzConstRef(reconstruction_inverse_den),
                        inverse_denominator,
                        flint::FmpzConstRef(exponent_fmpz), diagnostics)) {
                return false;
            }

            bool stabilized = false;
            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                    "element.power_hensel_stabilized_check");
                stabilized = have_previous[root_index] &&
                             candidate.equal(previous_candidates[root_index]);
            }
            if (stabilized || final_lift) {
                if (verify_power_root(
                            candidate, power_input,
                            flint::FmpzConstRef(exponent_fmpz),
                            diagnostics)) {
                    root.swap(candidate);
                    return true;
                }
            } else if (!previous_candidates[root_index].set(candidate)) {
                return false;
            } else {
                have_previous[root_index] = true;
            }
            if (final_lift) {
                have_previous[root_index] = false;
            }
        }
    }

    return false;
}

// Pre-lift residue disproof for pure_power_hensel_prime_search.  Silex
// extension, not in the reference: `_roots_hensel` (reference v0.38.6
// `src/NumFieldOrd/NfOrd/Hensel.jl` lines 57-232) tests only the first
// factor of the first good prime before lifting when the degree is below 10.
// Lifting to the Friedrich-Fieker exponent is expensive and, when no lifted
// candidate verifies, Silex fails closed instead of answering "not a power".
// So before lifting, the residue of the radicand is tested at every factor of
// the chosen prime and of a few further admissible primes; one non-power
// residue is an exact proof that the radicand is not an n-th power in K.
// Finding none changes nothing: the lift runs as before.
//
// Admissible primes and why a non-power residue proves non-power.  Let
// K = Q(theta) with theta a root of the monic integral f of degree N, and let
// a in K be the radicand.  A rational prime p is used only when
//   (1) f mod p is squarefree (of degree N, as f is monic), so p does not
//       divide disc(f); hence p is unramified in K and does not divide
//       [O_K : Z[theta]], and by Dedekind-Kummer (Cohen, GTM 138, Theorem
//       4.8.13) p O_K is the product of the distinct primes
//       P_i = (p, g_i(theta)) for the monic irreducible factors g_i of
//       f mod p, with O_K / P_i = F_p[x] / (g_i) via theta -> x;
//   (2) every power-basis coordinate of a is p-integral, so a lies in
//       Z_(p)[theta], inside the valuation ring of every P_i, and its residue
//       mod P_i is a(x) mod (p, g_i);
//   (3) only factors g_i with a(x) mod g_i nonzero are tested, so
//       v_{P_i}(a) = 0 there.
// If a = b^n with b in K, then n v_{P_i}(b) = v_{P_i}(a) = 0, so b is a
// P_i-unit and reducing modulo P_i gives a residue b' with b'^n = a mod P_i.
// Hence a mod P_i is an n-th power in F_q^*, q = p^deg(g_i).  Since F_q^* is
// cyclic of order q - 1, the residue is an n-th power exactly when its
// ((q - 1) / gcd(n, q - 1))-th power is 1.  A residue failing that test
// therefore proves that a is not an n-th power in K.  Nothing in the argument
// needs p to be coprime to n.  The radicand of the prime search is input * d^n
// (d = 1 on the integral path), an n-th power exactly when input is.
//
// The factors are not split: for each degree f, the distinct-degree
// factorisation gives the product G of the factors of degree f, and
// G' = G / gcd(a, G) is the product of those where a is nonzero (G is
// squarefree).  By the Chinese remainder theorem F_p[x] / (G') is the product
// of the fields F_p[x] / (g_i), all of order q = p^f, so
// a^((q - 1) / gcd(n, q - 1)) mod G' is 1 exactly when it is 1 at every
// such g_i.  A value other than 1 is therefore a non-power residue at some
// P_i.
//
// Resource bound only: how many further admissible primes are examined after
// the chosen one.  Any value is sound; exhausting it leaves the call on the
// unchanged lifting path.  Disproofs can be rare: for n = 2 every element of
// F_p is a square in F_{p^f} when f is even, so c * b^2 with c in Q can only
// be disproved at the P_i of odd residue degree, and in a multiquadratic
// field of degree N (residue degrees at most 2) only at the primes that split
// completely, about 1/N of them.  The grh class saturation rows
// multiquadratic_m3_m5_m7 and lmfdb_8_8_73116160000_2 need up to 21 and 36
// further primes; 40 settles all of their undecided squares, while one prime
// costs a few microseconds (a distinct-degree factorisation and one modular
// power per residue degree).
constexpr slong kPreLiftDisproofExtraPrimes = 40;

// Resource bound only: how many admissible primes the residue disproof in
// degree 10 and higher examines (pure_power_residue_disproves).  Any value is
// sound; exhausting it leaves the call unsupported.  A true power is never
// disproved and always runs to the bound, at a distinct-degree factorisation
// per prime (a fraction of a millisecond in degree 20), so the bound is the
// same as kPreLiftDisproofExtraPrimes.  A residue can
// fail to be an n-th power only at a factor whose residue field order q has
// gcd(n, q - 1) > 1, so a prime p with gcd(n, p^f - 1) = 1 for every
// f <= N is skipped exactly, without factoring, and is not counted
// (residue_degree_can_disprove_power); the outer 256-attempt cap still counts
// skipped primes and ends the scan for a large n.  On the compact-root kernel rows of
// the proven runs on Q(zeta_25) and Q(zeta_44), 9 non-powers each needed 9
// to 218 admissible primes without that skip; with it, all 18 are disproved
// within the bound.
constexpr slong kHighDegreeDisproofPrimes = 40;

// Reduces the rational polynomial `polynomial` modulo the word prime of
// `out`; false when a coefficient is not p-integral.
bool fmpq_poly_reduce_nmod(flint::NmodPoly& out,
                           const fmpq_poly_t polynomial) noexcept {
    const ulong prime = out.raw()->mod.n;
    const ulong denominator = fmpz_fdiv_ui(polynomial->den, prime);
    if (denominator == 0) {
        return false;
    }
    const ulong inverse = n_invmod(denominator, prime);
    nmod_poly_zero(out.raw());
    for (slong i = 0; i < polynomial->length; ++i) {
        const ulong coefficient = fmpz_fdiv_ui(polynomial->coeffs + i, prime);
        nmod_poly_set_coeff_ui(
                out.raw(), i,
                n_mulmod2_preinv(coefficient, inverse, prime,
                                 out.raw()->mod.ninv));
    }
    return true;
}

// True when some residue field of order q = p^f, 1 <= f <= degree, has
// gcd(n, q - 1) > 1.  When this is false, no residue at any prime above p
// can disprove that a is an n-th power, so the caller may skip p without
// factoring f mod p.  Exactness: if p satisfies (1) above, f mod p is
// squarefree of degree N = `degree`, so every irreducible factor g_i has
// degree at most N and every residue field O_K / P_i = F_p[x] / (g_i) has
// order q = p^deg(g_i) with deg(g_i) in [1, N].  As stated above
// residue_disproves_power_at_prime, F_q^* is cyclic of order q - 1, so x ->
// x^n maps F_q^* onto its subgroup of index gcd(n, q - 1); when that gcd is 1
// every nonzero residue is an n-th power (Lidl and Niederreiter, Finite
// Fields, Theorem 2.8, cyclicity of F_q^*), and
// residue_disproves_power_at_prime skips that residue degree.  Hence
// residue_disproves_power_at_prime returns false at every p for which this
// returns false, admissible or not.
bool residue_degree_can_disprove_power(ulong p, slong degree,
                                       slong exponent) noexcept {
    const ulong n = static_cast<ulong>(exponent);
    const ulong base = p % n;
    ulong power = 1;  // p^f mod n
    for (slong f = 1; f <= degree; ++f) {
        power = n_mulmod2(power, base, n);
        // gcd(n, p^f - 1) = gcd(n, (p^f - 1) mod n), with gcd(n, 0) = n.
        const ulong reduced = power == 0 ? n - 1 : power - 1;
        if (n_gcd(n, reduced) > 1) {
            return true;
        }
    }
    return false;
}

// Residue test at one rational prime p (see above).  Sets `admissible` when
// p satisfies conditions (1) and (2), and returns true only when some factor
// of f mod p with a nonzero residue has a residue that is not an n-th power.
bool residue_disproves_power_at_prime(bool& admissible,
                                      const fmpq_poly_t field_polynomial,
                                      const fmpq_poly_t element_polynomial,
                                      slong degree,
                                      const fmpz_t prime,
                                      slong exponent) noexcept {
    admissible = false;
    if (exponent <= 1 || degree <= 0 || fmpz_sgn(prime) <= 0 ||
        fmpz_abs_fits_ui(prime) == 0) {
        return false;
    }
    const ulong p = fmpz_get_ui(prime);
    flint::NmodPoly field_mod(p);
    flint::NmodPoly element(p);
    if (!fmpq_poly_reduce_nmod(field_mod, field_polynomial) ||
        nmod_poly_degree(field_mod.raw()) != degree ||
        nmod_poly_is_squarefree(field_mod.raw()) == 0 ||
        !fmpq_poly_reduce_nmod(element, element_polynomial)) {
        return false;
    }
    admissible = true;
    nmod_poly_make_monic(field_mod.raw(), field_mod.raw());

    flint::NmodPolyFactor blocks;
    std::vector<slong> block_degrees(static_cast<std::size_t>(degree) + 1, 0);
    slong* block_degree_data = block_degrees.data();
    nmod_poly_factor_distinct_deg(blocks.raw(), field_mod.raw(),
                                  &block_degree_data);

    flint::NmodPoly common(p);
    flint::NmodPoly quotient(p);
    flint::NmodPoly residue(p);
    flint::NmodPoly power(p);
    flint::Fmpz group_order;
    flint::Fmpz gcd;
    for (slong i = 0; i < blocks.raw()->num; ++i) {
        const nmod_poly_struct* block = blocks.raw()->p + i;
        const slong residue_degree = block_degrees[static_cast<std::size_t>(i)];
        if (residue_degree <= 0) {
            continue;
        }
        nmod_poly_rem(residue.raw(), element.raw(), block);
        nmod_poly_gcd(common.raw(), residue.raw(), block);
        nmod_poly_div(quotient.raw(), block, common.raw());
        if (nmod_poly_degree(quotient.raw()) <= 0) {
            continue;
        }
        nmod_poly_rem(residue.raw(), residue.raw(), quotient.raw());

        fmpz_set_ui(group_order.raw(), p);
        fmpz_pow_ui(group_order.raw(), group_order.raw(),
                    static_cast<ulong>(residue_degree));
        fmpz_sub_ui(group_order.raw(), group_order.raw(), 1);
        fmpz_set_si(gcd.raw(), exponent);
        fmpz_gcd(gcd.raw(), gcd.raw(), group_order.raw());
        if (fmpz_is_one(gcd.raw()) != 0) {
            // Every element of F_q^* is an n-th power.
            continue;
        }
        fmpz_divexact(group_order.raw(), group_order.raw(), gcd.raw());
        nmod_poly_powmod_fmpz_binexp(power.raw(), residue.raw(),
                                     group_order.raw(), quotient.raw());
        if (nmod_poly_is_one(power.raw()) == 0) {
            return true;
        }
    }
    return false;
}

// True when a residue at a factor of the chosen prime or of one of the next
// kPreLiftDisproofExtraPrimes admissible primes disproves that
// `element_polynomial` is an n-th power.
bool pre_lift_residue_disproves(const flint::FmpqPoly& element_polynomial,
                                const nf_struct* raw_field,
                                slong degree,
                                const flint::Fmpz& chosen_prime,
                                slong exponent,
                                const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_prelift_residue_disproof");
    bool admissible = false;
    if (residue_disproves_power_at_prime(admissible, raw_field->pol,
                                         element_polynomial.raw(), degree,
                                         chosen_prime.raw(), exponent)) {
        return true;
    }

    flint::Fmpz prime;
    fmpz_set(prime.raw(), chosen_prime.raw());
    slong admissible_primes = 0;
    constexpr slong max_prime_attempts = 256;
    for (slong attempt = 0; attempt < max_prime_attempts &&
                            admissible_primes < kPreLiftDisproofExtraPrimes;
         ++attempt) {
        fmpz_nextprime(prime.raw(), prime.raw(), 1);
        if (residue_disproves_power_at_prime(admissible, raw_field->pol,
                                             element_polynomial.raw(), degree,
                                             prime.raw(), exponent)) {
            return true;
        }
        if (admissible) {
            ++admissible_primes;
        }
    }
    return false;
}

// Radicand-independent part of one prime of pure_power_hensel_prime_search:
// f mod p and its factorization.  Every search walks the same primes 5, 7,
// 11, ... and used to reduce and factor f mod p again for each radicand.
//
// Key: the defining polynomial, backend and the prime, which are all that
// the reduction and fmpz_mod_poly_factor read; the factorization is
// deterministic, so the stored factors and their order are those a
// recomputation returns.  `usable` records the two radicand-independent
// reasons the search skips a prime (f mod p loses degree, or is not
// squarefree).  Factors are stored as integer lifts in [0, p) and rebuilt in
// the caller's modulus context.
struct PrimeSearchSetupEntry {
    flint::FmpqPoly polynomial;
    NumberFieldBackendKind backend = NumberFieldBackendKind::generic;
    flint::Fmpz prime;
    bool usable = false;
    flint::FmpzPoly field_mod_prime;
    std::vector<flint::FmpzPoly> factors;
};

// A search rarely needs more than a few primes per field; the cap bounds
// memory for the 256-attempt worst case and for many fields.
constexpr std::size_t kPrimeSearchSetupCacheCapacity = 64;

BoundedEntryCache<PrimeSearchSetupEntry>& prime_search_setup_cache() noexcept {
    thread_local BoundedEntryCache<PrimeSearchSetupEntry> cache(
            kPrimeSearchSetupCacheCapacity);
    return cache;
}

// Fills `field_polynomial` with f mod p and `factors` with the factors of f
// mod p in fmpz_mod_poly_factor order, and sets `usable` as above.  Returns
// false only on an allocation failure.
bool prime_search_setup(bool& usable,
                        flint::FmpzModPoly& field_polynomial,
                        std::vector<flint::FmpzModPoly>& factors,
                        const nf_struct* raw_field,
                        NumberFieldBackendKind backend,
                        slong degree,
                        const flint::Fmpz& prime,
                        const flint::FmpzModCtx& ctx,
                        const DiagnosticsContext* diagnostics) noexcept {
    factors.clear();
    const auto rebuild = [&](const PrimeSearchSetupEntry& entry) noexcept {
        usable = entry.usable;
        fmpz_mod_poly_set_fmpz_poly(field_polynomial.raw(),
                                    entry.field_mod_prime.raw(), ctx.raw());
        factors.reserve(entry.factors.size());
        for (const flint::FmpzPoly& lifted : entry.factors) {
            factors.emplace_back(ctx);
            if (!factors.back().is_initialized()) {
                return false;
            }
            fmpz_mod_poly_set_fmpz_poly(factors.back().raw(), lifted.raw(),
                                        ctx.raw());
        }
        return true;
    };

    auto& cache = prime_search_setup_cache();
    for (const PrimeSearchSetupEntry& entry : cache.entries()) {
        if (entry.backend == backend &&
            fmpz_equal(entry.prime.raw(), prime.raw()) != 0 &&
            fmpq_poly_equal(entry.polynomial.raw(), raw_field->pol) != 0) {
            ++hensel_data_cache_counters_tls.prime_setup_hits;
            return rebuild(entry);
        }
    }

    ++hensel_data_cache_counters_tls.prime_setup_misses;
    PrimeSearchSetupEntry entry;
    fmpq_poly_set(entry.polynomial.raw(), raw_field->pol);
    entry.backend = backend;
    fmpz_set(entry.prime.raw(), prime.raw());
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_prime_project");
        if (!field_polynomial.is_initialized()) {
            return false;
        }
        entry.usable =
                fmpq_poly_reduce_modulus(field_polynomial, raw_field->pol,
                                         ctx, prime.raw()) &&
                fmpz_mod_poly_degree(field_polynomial.raw(), ctx.raw()) ==
                        degree;
    }
    if (entry.usable) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_field_factor");
        flint::FmpzModPolyFactor factorization(ctx);
        fmpz_mod_poly_factor(factorization.raw(), field_polynomial.raw(),
                             ctx.raw());
        entry.usable =
                fmpz_mod_poly_factorization_is_squarefree(factorization);
        fmpz_mod_poly_get_fmpz_poly(entry.field_mod_prime.raw(),
                                    field_polynomial.raw(), ctx.raw());
        if (entry.usable) {
            entry.factors.reserve(
                    static_cast<std::size_t>(factorization.raw()->num));
            for (slong i = 0; i < factorization.raw()->num; ++i) {
                flint::FmpzModPoly factor(ctx);
                if (!factor.is_initialized()) {
                    return false;
                }
                fmpz_mod_poly_factor_get_poly(factor.raw(),
                                              factorization.raw(), i,
                                              ctx.raw());
                entry.factors.emplace_back();
                fmpz_mod_poly_get_fmpz_poly(entry.factors.back().raw(),
                                            factor.raw(), ctx.raw());
            }
        }
    }
    if (!rebuild(entry)) {
        return false;
    }
    cache.insert(std::move(entry));
    return true;
}

// Prime search and Hensel lift of pure_power_hensel_root for the radicand
// `radicand` = input * d^n, whose power-basis coordinates are
// `element_polynomial`; root_scale = 1/d, or nullptr when radicand is input
// itself (d = 1).  A lifted root is published only after exact verification
// against input.
bool pure_power_hensel_prime_search(
        bool& is_power,
        Element& root,
        const Element& radicand,
        const Element& input,
        const Element* root_scale,
        const flint::FmpqPoly& element_polynomial,
        const NumberField* parent,
        const nf_struct* raw_field,
        slong exponent,
        detail::PreLiftResidueDisproof pre_lift,
        const DiagnosticsContext* diagnostics) noexcept {
    flint::Fmpz prime;
    fmpz_set_ui(prime.raw(), 3);
    constexpr slong max_prime_attempts = 256;
    for (slong attempt = 0; attempt < max_prime_attempts; ++attempt) {
        fmpz_nextprime(prime.raw(), prime.raw(), 1);
        flint::FmpzModCtx ctx(prime.raw());
        if (ctx.raw() == nullptr) {
            return false;
        }

        // The field part (f mod p, its degree and squarefree
        // factorization) comes from the per-(field, prime) memo above; a
        // prime is skipped for the same reasons as before, and an element
        // that does not reduce mod p skips it too.
        flint::FmpzModPoly field_polynomial(ctx);
        flint::FmpzModPoly element(ctx);
        bool usable = false;
        std::vector<flint::FmpzModPoly> factors;
        if (!field_polynomial.is_initialized() || !element.is_initialized()) {
            continue;
        }
        if (!prime_search_setup(usable, field_polynomial, factors, raw_field,
                                parent->backend_kind(), parent->degree(),
                                prime, ctx, diagnostics)) {
            return false;
        }
        if (!usable) {
            continue;
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.power_hensel_prime_project");
            if (!fmpq_poly_reduce_modulus(element,
                                          element_polynomial.raw(), ctx,
                                          prime.raw())) {
                continue;
            }
        }

        for (flint::FmpzModPoly& factor : factors) {
            std::vector<flint::FmpzPoly> inverse_roots;
            const ResiduePowerStatus status =
                    residue_field_power_inverse_roots(
                            inverse_roots, element, factor, ctx, prime.raw(),
                            exponent, diagnostics);
            if (status == ResiduePowerStatus::nonpower) {
                is_power = false;
                return true;
            }
            if (status == ResiduePowerStatus::power) {
                // Exact disproof before the lift; see
                // pre_lift_residue_disproves for the admissible primes.
                if (pre_lift == detail::PreLiftResidueDisproof::enabled &&
                    pre_lift_residue_disproves(
                            element_polynomial, raw_field, parent->degree(),
                            prime, exponent, diagnostics)) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::element,
                            "element.power_prelift_residue_disproof.nonpower");
                    is_power = false;
                    return true;
                }
                if (try_pure_power_root_at_prime(
                            root, radicand, input, root_scale,
                            element_polynomial, raw_field, prime,
                            field_polynomial, factor, ctx,
                            std::move(inverse_roots), exponent, diagnostics)) {
                    is_power = true;
                    return true;
                }
                return false;
            }
        }
    }

    return false;
}

// Root-of-unity search of reference `_torsion_units_gen`: a root in K of the
// cyclotomic polynomial f = Phi_m, found as reference
// `_roots_hensel(f, max_roots = 1, is_normal = true, root_bound = ones)`
// finds it (reference v0.38.6 `src/NumFieldOrd/NfOrd/TorsionUnits.jl` lines
// 414-456 and `Hensel.jl` lines 57-232 and 313-637), for monic integral
// fields.  The roots of Phi_m are the primitive m-th roots of unity, so the
// candidate is certified by its exact multiplicative order.

// |f'(theta)^(i)| at the r1 + r2 places: the scaled root bounds of
// `_roots_hensel` (Hensel.jl lines 202-205) with root_bound = 1, since every
// conjugate of a root of unity has absolute value 1.
bool unit_circle_scaled_root_bounds(
        std::vector<flint::Arb>& bounds,
        const NumberField& parent,
        const flint::FmpqPoly& derivative_polynomial,
        slong precision) noexcept {
    EmbeddingContext embeddings(parent);
    if (!embeddings.refine(precision)) {
        return false;
    }
    Element derivative(parent);
    if (!derivative.is_defined() ||
        !derivative.set_fmpq_poly(flint::FmpqPolyConstRef(
                derivative_polynomial))) {
        return false;
    }

    const slong degree = parent.degree();
    flint::AcbVec values(degree);
    if (!embeddings.evaluate_all(flint::AcbVecRef(values), derivative,
                                 precision)) {
        return false;
    }

    const Signature sig = embeddings.signature();
    const slong r1 = sig.r1();
    const slong places = sig.r1() + sig.r2();
    bounds.clear();
    bounds.reserve(static_cast<std::size_t>(places));
    for (slong place = 0; place < places; ++place) {
        const slong root_index = place < r1 ? place : r1 + 2 * (place - r1);
        flint::Arb bound;
        flint::acb_abs(bound, values.data() + root_index, precision);
        if (!flint::arb_is_finite(bound)) {
            return false;
        }
        bounds.push_back(std::move(bound));
    }
    return bounds.size() == static_cast<std::size_t>(places);
}

// The residue data `_roots_hensel` selects: a prime p, a factor g of the
// field polynomial modulo p (a prime P above p of residue degree deg g), a
// root r of f modulo P and 1/f'(r) modulo P, the last two as polynomials in
// F_p[t] reduced modulo g.
struct CyclotomicResidueRoot {
    ulong prime = 0;
    slong residue_degree = 0;
    flint::FmpzPoly factor;
    flint::FmpzPoly root;
    flint::FmpzPoly inverse_derivative;
};

bool fq_to_fmpz_poly(flint::FmpzPoly& out,
                     const fq_struct* value,
                     const flint::FqCtx& field,
                     const flint::FmpzModCtx& ctx) noexcept {
    flint::FmpzModPoly reduced(ctx);
    if (!reduced.is_initialized()) {
        return false;
    }
    fq_get_fmpz_mod_poly(reduced.raw(), value, field.raw());
    fmpz_mod_poly_get_fmpz_poly(out.raw(), reduced.raw(), ctx.raw());
    return true;
}

// A root of f modulo the prime of `field` (F_q = F_p[t]/(g)) and the inverse
// of f' at that root.  False when f has no root in F_q.
bool fq_cyclotomic_root(CyclotomicResidueRoot& out,
                        const FqPoly& f_mod,
                        const flint::FqCtx& field,
                        const flint::FmpzModCtx& ctx) noexcept {
    FqPolyFactor roots(field.raw());
    fq_poly_roots(roots.raw(), f_mod.raw(), 0, field.raw());
    if (roots.raw()->num == 0) {
        return false;
    }
    const fq_poly_struct* linear = roots.raw()->poly;
    if (fq_poly_degree(linear, field.raw()) != 1) {
        return false;
    }

    flint::Fq constant(field);
    flint::Fq leading(field);
    flint::Fq root(field);
    flint::Fq value(field);
    if (!constant.is_initialized() || !leading.is_initialized() ||
        !root.is_initialized() || !value.is_initialized()) {
        return false;
    }
    fq_poly_get_coeff(constant.raw(), linear, 0, field.raw());
    fq_poly_get_coeff(leading.raw(), linear, 1, field.raw());
    if (fq_is_zero(leading.raw(), field.raw()) != 0) {
        return false;
    }
    fq_inv(leading.raw(), leading.raw(), field.raw());
    fq_neg(root.raw(), constant.raw(), field.raw());
    fq_mul(root.raw(), root.raw(), leading.raw(), field.raw());

    FqPoly derivative(field.raw());
    fq_poly_derivative(derivative.raw(), f_mod.raw(), field.raw());
    fq_poly_evaluate_fq(value.raw(), derivative.raw(), root.raw(),
                        field.raw());
    if (fq_is_zero(value.raw(), field.raw()) != 0) {
        return false;
    }
    fq_inv(value.raw(), value.raw(), field.raw());
    return fq_to_fmpz_poly(out.root, root.raw(), field, ctx) &&
           fq_to_fmpz_poly(out.inverse_derivative, value.raw(), field, ctx);
}

enum class CyclotomicPrimeStatus { selected, no_root, unsupported };

// Prime selection of `_roots_hensel` (Hensel.jl lines 83-181) for
// is_normal = true.  Primes start above deg f + 1 and must keep the field
// polynomial squarefree of full degree.  For each prime P above p with f
// squarefree modulo P, no residue root proves that f has no root in K, and,
// because the roots of Phi_m generate a normal extension, fewer than deg f
// residue roots proves the same (lines 159-165).  Among the remaining
// primes the reference prefers fewer residue roots, then larger residue
// degree, and stops once STABILIZED (2 for n >= 10, else 0) primes in a row
// do not improve the choice.
CyclotomicPrimeStatus select_cyclotomic_prime(
        CyclotomicResidueRoot& out,
        const flint::FmpzPoly& f,
        const nf_struct* raw_field,
        slong degree) noexcept {
    const slong f_degree = fmpz_poly_degree(f.raw());
    const slong stabilized = degree >= 10 ? 2 : 0;
    slong current_num = f_degree + 1;
    slong not_better = 1;
    bool have = false;

    // Resource bound only; reaching it leaves the search unsupported.
    constexpr slong max_prime_attempts = 256;
    ulong p = static_cast<ulong>(f_degree) + 1;
    for (slong attempt = 0; attempt < max_prime_attempts; ++attempt) {
        p = n_nextprime(p, 1);
        flint::Fmpz prime;
        fmpz_set_ui(prime.raw(), p);
        flint::FmpzModCtx ctx(prime.raw());
        flint::FmpzModPoly field_polynomial(ctx);
        if (ctx.raw() == nullptr || !field_polynomial.is_initialized()) {
            return CyclotomicPrimeStatus::unsupported;
        }
        if (!fmpq_poly_reduce_modulus(field_polynomial, raw_field->pol, ctx,
                                      prime.raw()) ||
            fmpz_mod_poly_degree(field_polynomial.raw(), ctx.raw()) !=
                    degree ||
            fmpz_mod_poly_is_squarefree(field_polynomial.raw(),
                                        ctx.raw()) == 0) {
            continue;
        }

        flint::FmpzModPolyFactor factorization(ctx);
        fmpz_mod_poly_factor(factorization.raw(), field_polynomial.raw(),
                             ctx.raw());
        for (slong i = 0; i < factorization.raw()->num; ++i) {
            flint::FmpzModPoly factor(ctx);
            if (!factor.is_initialized()) {
                return CyclotomicPrimeStatus::unsupported;
            }
            fmpz_mod_poly_factor_get_poly(factor.raw(), factorization.raw(),
                                          i, ctx.raw());
            const slong residue_degree =
                    fmpz_mod_poly_degree(factor.raw(), ctx.raw());
            flint::FqCtx field(factor, ctx, "z");
            if (field.raw() == nullptr) {
                return CyclotomicPrimeStatus::unsupported;
            }
            FqPoly f_mod(field.raw());
            for (slong j = 0; j <= f_degree; ++j) {
                fq_poly_set_coeff_fmpz(f_mod.raw(), j, f.raw()->coeffs + j,
                                       field.raw());
            }
            if (fq_poly_is_squarefree(f_mod.raw(), field.raw()) == 0) {
                continue;
            }

            FqPolyFactor roots(field.raw());
            fq_poly_roots(roots.raw(), f_mod.raw(), 0, field.raw());
            const slong num = roots.raw()->num;
            if (num == 0 || num < f_degree) {
                return CyclotomicPrimeStatus::no_root;
            }

            if (num < current_num ||
                (num == current_num && residue_degree > out.residue_degree)) {
                if (!fq_cyclotomic_root(out, f_mod, field, ctx)) {
                    return CyclotomicPrimeStatus::unsupported;
                }
                out.prime = p;
                out.residue_degree = residue_degree;
                fmpz_mod_poly_get_fmpz_poly(out.factor.raw(), factor.raw(),
                                            ctx.raw());
                current_num = num;
                have = true;
            } else {
                ++not_better;
            }
            if (not_better >= stabilized) {
                return have ? CyclotomicPrimeStatus::selected
                            : CyclotomicPrimeStatus::unsupported;
            }
        }
    }
    return CyclotomicPrimeStatus::unsupported;
}

// f(x) modulo (p^k, modulus) for x reduced modulo modulus.
void evaluate_mod(flint::FmpzModPoly& out,
                  const flint::FmpzModPoly& polynomial,
                  const flint::FmpzModPoly& x,
                  const flint::FmpzModPoly& modulus,
                  const flint::FmpzModCtx& ctx) noexcept {
    fmpz_mod_poly_compose_mod(out.raw(), polynomial.raw(), x.raw(),
                              modulus.raw(), ctx.raw());
}

// Lifting and reconstruction of `_hensel` (Hensel.jl lines 313-637) for a
// non-pure f: Newton steps x -> x - f(x) / f'(x) with the inverse of f'(x)
// lifted by y -> y (2 - y f'(x)) (lines 593-596), scaling by the Kronecker
// denominator T'(theta) of the field polynomial T before the lattice
// reconstruction (lines 425 and 599-611), and an exact check of every
// candidate that stabilizes or reaches the Friedrich-Fieker exponent (lines
// 617-631).
//
// Deviation: the reference lifts every residue root when is_normal is set
// (line 411; its comment at lines 407-408 says only max_roots are needed).
// Silex lifts one.  This loses no root: P is unramified and f is squarefree
// modulo P, so the roots of f in K reduce to distinct residue roots, and when
// f is normal with a root in K it has deg f roots in K, whose reductions are
// then all the residue roots.
bool lift_cyclotomic_root(Element& root,
                          const NumberField& field,
                          const flint::FmpzPoly& f,
                          ulong m,
                          const CyclotomicResidueRoot& residue) noexcept {
    const nf_struct* raw_field = field.raw_flint_field();
    const slong degree = field.degree();
    flint::Fmpz prime;
    fmpz_set_ui(prime.raw(), residue.prime);
    flint::FmpzModCtx prime_ctx(prime.raw());
    flint::FmpzModPoly field_mod_prime(prime_ctx);
    flint::FmpzModPoly factor_mod_prime(prime_ctx);
    if (prime_ctx.raw() == nullptr || !field_mod_prime.is_initialized() ||
        !factor_mod_prime.is_initialized() ||
        !fmpq_poly_reduce_modulus(field_mod_prime, raw_field->pol, prime_ctx,
                                  prime.raw())) {
        return false;
    }
    fmpz_mod_poly_set_fmpz_poly(factor_mod_prime.raw(),
                                residue.factor.raw(), prime_ctx.raw());

    flint::FmpqPoly derivative_polynomial;
    fmpq_poly_derivative(derivative_polynomial.raw(), raw_field->pol);
    Element inverse_denominator(field);
    if (!inverse_denominator.is_defined() ||
        !inverse_denominator.set_fmpq_poly(
                flint::FmpqPolyConstRef(derivative_polynomial)) ||
        !inverse_denominator.invert(inverse_denominator)) {
        return false;
    }

    constexpr slong precision = 128;
    std::vector<flint::Arb> root_bounds;
    slong target_lift_exponent = 0;
    std::vector<slong> lift_chain;
    if (!unit_circle_scaled_root_bounds(root_bounds, field,
                                        derivative_polynomial, precision) ||
        !friedrich_fieker_lifting_exponent(
                target_lift_exponent, field, root_bounds, prime,
                residue.residue_degree, precision, nullptr) ||
        !lifting_chain(lift_chain, target_lift_exponent)) {
        return false;
    }

    flint::FmpzPoly f_derivative;
    fmpz_poly_derivative(f_derivative.raw(), f.raw());
    flint::FmpzPoly lifted_root;
    flint::FmpzPoly lifted_inverse;
    fmpz_poly_set(lifted_root.raw(), residue.root.raw());
    fmpz_poly_set(lifted_inverse.raw(), residue.inverse_derivative.raw());

    flint::Fmpz one;
    fmpz_one(one.raw());
    Element previous(field);
    bool have_previous = false;
    if (!previous.is_defined()) {
        return false;
    }
    for (std::size_t step = 1; step < lift_chain.size(); ++step) {
        const slong lift_exponent = lift_chain[step];
        const bool final_lift = step + 1 == lift_chain.size();
        flint::Fmpz precision_modulus;
        fmpz_pow_ui(precision_modulus.raw(), prime.raw(),
                    static_cast<ulong>(lift_exponent));
        flint::FmpzModCtx ctx(precision_modulus.raw());
        if (ctx.raw() == nullptr) {
            return false;
        }
        flint::FmpzModPoly lifted_factor(ctx);
        flint::FmpzModPoly f_mod(ctx);
        flint::FmpzModPoly fs_mod(ctx);
        flint::FmpzModPoly den(ctx);
        flint::FmpzModPoly x(ctx);
        flint::FmpzModPoly y(ctx);
        flint::FmpzModPoly f_value(ctx);
        flint::FmpzModPoly fs_value(ctx);
        flint::FmpzModPoly correction(ctx);
        flint::FmpzModPoly two(ctx);
        flint::FmpzModPoly scaled(ctx);
        if (!lifted_factor.is_initialized() || !f_mod.is_initialized() ||
            !fs_mod.is_initialized() || !den.is_initialized() ||
            !x.is_initialized() || !y.is_initialized() ||
            !f_value.is_initialized() || !fs_value.is_initialized() ||
            !correction.is_initialized() || !two.is_initialized() ||
            !scaled.is_initialized() ||
            !hensel_lift_factor_to_precision(
                    lifted_factor, raw_field->pol, field_mod_prime,
                    factor_mod_prime, prime_ctx, prime.raw(), lift_exponent,
                    precision_modulus.raw(), ctx) ||
            !fmpq_poly_reduce_modulus(den, derivative_polynomial.raw(), ctx,
                                      precision_modulus.raw())) {
            return false;
        }
        fmpz_mod_poly_rem(den.raw(), den.raw(), lifted_factor.raw(),
                          ctx.raw());
        fmpz_mod_poly_set_fmpz_poly(f_mod.raw(), f.raw(), ctx.raw());
        fmpz_mod_poly_set_fmpz_poly(fs_mod.raw(), f_derivative.raw(),
                                    ctx.raw());
        fmpz_mod_poly_set_fmpz_poly(x.raw(), lifted_root.raw(), ctx.raw());
        fmpz_mod_poly_rem(x.raw(), x.raw(), lifted_factor.raw(), ctx.raw());
        fmpz_mod_poly_set_fmpz_poly(y.raw(), lifted_inverse.raw(), ctx.raw());
        fmpz_mod_poly_rem(y.raw(), y.raw(), lifted_factor.raw(), ctx.raw());
        fmpz_mod_poly_set_ui(two.raw(), 2, ctx.raw());

        // y <- y (2 - y f'(x)); x <- x - f(x) y, both at the old x.
        evaluate_mod(f_value, f_mod, x, lifted_factor, ctx);
        evaluate_mod(fs_value, fs_mod, x, lifted_factor, ctx);
        fmpz_mod_poly_mulmod(correction.raw(), y.raw(), fs_value.raw(),
                             lifted_factor.raw(), ctx.raw());
        fmpz_mod_poly_sub(correction.raw(), two.raw(), correction.raw(),
                          ctx.raw());
        fmpz_mod_poly_mulmod(y.raw(), y.raw(), correction.raw(),
                             lifted_factor.raw(), ctx.raw());
        fmpz_mod_poly_mulmod(correction.raw(), f_value.raw(), y.raw(),
                             lifted_factor.raw(), ctx.raw());
        fmpz_mod_poly_sub(x.raw(), x.raw(), correction.raw(), ctx.raw());
        fmpz_mod_poly_get_fmpz_poly(lifted_root.raw(), x.raw(), ctx.raw());
        fmpz_mod_poly_get_fmpz_poly(lifted_inverse.raw(), y.raw(), ctx.raw());

        fmpz_mod_poly_mulmod(scaled.raw(), x.raw(), den.raw(),
                             lifted_factor.raw(), ctx.raw());

        flint::FmpzMat basis(degree, degree);
        flint::FmpzMat inverse_num(degree, degree);
        flint::Fmpz inverse_den;
        Element candidate(field);
        if (!candidate.is_defined() ||
            !power_reconstruction_data(basis, inverse_num, inverse_den,
                                       degree, lifted_factor, ctx,
                                       precision_modulus.raw(), nullptr) ||
            // `input` and `exponent` are placeholders here: this cyclotomic
            // root is reconstructed and then verified by exact order
            // (has_exact_order below), not by norm_admits_power_root, so
            // `input` only needs a defined parent for the `has_parent` guard
            // and `exponent` only needs to be positive.  `inverse_denominator`
            // already has the right parent, and 1 is trivially positive.
            !reconstruct_power_root_candidate(
                    candidate, inverse_denominator, scaled, ctx, basis,
                    inverse_num, flint::FmpzConstRef(inverse_den),
                    inverse_denominator, flint::FmpzConstRef(one),
                    nullptr)) {
            return false;
        }

        if ((have_previous && candidate.equal(previous)) || final_lift) {
            if (detail::has_exact_order(candidate, m)) {
                root.swap(candidate);
                return true;
            }
            // The reference returns no roots here when the final lift fails
            // and f is normal; Silex reports failure either way.
        } else {
            if (!previous.set(candidate)) {
                return false;
            }
            have_previous = true;
        }
    }
    return false;
}

}  // namespace

namespace detail {

HenselDataCacheCounters hensel_data_cache_counters() noexcept {
    return hensel_data_cache_counters_tls;
}

void clear_hensel_data_caches() noexcept {
    norm_change_constants_cache().clear();
    reconstruction_data_cache().clear();
    conjugate_bounds_field_cache().clear();
    lifting_exponent_cache().clear();
    prime_search_setup_cache().clear();
    inverse_derivative_cache().clear();
}

bool has_exact_order(const Element& z, ulong m) noexcept {
    const NumberField* parent = z.parent();
    if (parent == nullptr || m == 0) {
        return false;
    }
    Element power(*parent);
    flint::Fmpz exponent;
    fmpz_set_ui(exponent.raw(), m);
    if (!power.is_defined() ||
        !power.pow_fmpz(z, flint::FmpzConstRef(exponent)) ||
        !power.equal_si(1)) {
        return false;
    }
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, m, 1);
    for (int i = 0; i < factors.num; ++i) {
        fmpz_set_ui(exponent.raw(), m / factors.p[i]);
        if (!power.pow_fmpz(z, flint::FmpzConstRef(exponent)) ||
            power.equal_si(1)) {
            return false;
        }
    }
    return true;
}

void generic_trace(flint::FmpqRef out,
                   const nf_struct* field,
                   slong degree,
                   const nf_elem_t element) noexcept {
    flint::FmpqMat matrix(degree, degree);
    multiplication_matrix(matrix.raw(), field, degree, element);
    fmpq_mat_trace(out.raw(), matrix.raw());
}

void generic_norm(flint::FmpqRef out,
                  const nf_struct* field,
                  slong degree,
                  const nf_elem_t element) noexcept {
    (void)degree;
    nf_elem_norm(out.raw(), element, field);
}

bool is_square_rational_constant(bool& is_square,
                                 Element& root,
                                 const Element& input) noexcept {
    flint::Fmpq value;
    flint::Fmpq candidate;

    if (!get_rational_constant(value.raw(), input)) {
        return false;
    }

    if (fmpq_sqrt(candidate.raw(), value.raw())) {
        is_square = true;
        return set_rational(root, candidate.raw());
    }

    is_square = false;
    return true;
}

bool is_square_quadratic(bool& is_square,
                         Element& root,
                         const Element& input,
                         const fmpz_t radicand) noexcept {
    flint::FmpqPoly polynomial;
    flint::Fmpq a;
    flint::Fmpq b;
    flint::Fmpq norm;
    flint::Fmpq tmp;
    flint::Fmpq m;

    if (!input.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }
    fmpq_poly_get_coeff_fmpq(a.raw(), polynomial.raw(), 0);
    fmpq_poly_get_coeff_fmpq(b.raw(), polynomial.raw(), 1);

    if (fmpq_is_zero(a.raw()) != 0 && fmpq_is_zero(b.raw()) != 0) {
        is_square = true;
        return root.zero();
    }

    fmpq_mul(norm.raw(), a.raw(), a.raw());
    fmpq_mul(tmp.raw(), b.raw(), b.raw());
    fmpq_mul_fmpz(tmp.raw(), tmp.raw(), radicand);
    fmpq_sub(norm.raw(), norm.raw(), tmp.raw());

    if (!fmpq_sqrt(m.raw(), norm.raw())) {
        is_square = false;
        return true;
    }
    if (quadratic_try_sqrt_m(root, input, a.raw(), b.raw(), radicand, m.raw())) {
        is_square = true;
        return true;
    }

    fmpq_neg(m.raw(), m.raw());
    is_square = quadratic_try_sqrt_m(root, input, a.raw(), b.raw(), radicand, m.raw());
    return true;
}

bool is_power_rational_constant(bool& is_power,
                                Element& root,
                                const Element& input,
                                flint::FmpzConstRef exponent) noexcept {
    flint::Fmpq value;
    flint::Fmpq candidate;

    if (!get_rational_constant(value.raw(), input)) {
        return false;
    }

    if (fmpz_fits_si(exponent.raw()) != 0) {
        is_power = fmpq_root(candidate.raw(), value.raw(), fmpz_get_si(exponent.raw()));
    } else {
        is_power = fmpq_huge_root(candidate.raw(), value.raw(), exponent);
    }

    if (is_power) {
        return set_rational(root, candidate.raw());
    }

    return true;
}

bool pure_power_hensel_root(bool& is_power,
                                  Element& root,
                                  const Element& input,
                                  slong exponent,
                                  const DiagnosticsContext* diagnostics,
                                  PreLiftResidueDisproof pre_lift) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_hensel_root");
    // Private port of reference `_roots_hensel(y^n - a, ispure=true,
    // is_normal=true)` for monic integral degree < 10 fields.  The result is
    // published only after exact verification; unported cases remain
    // unsupported.
    const NumberField* parent = input.parent();
    if (parent == nullptr || parent->degree() <= 1 ||
        parent->degree() >= 10 || exponent <= 1) {
        return false;
    }

    flint::FmpqPoly element_polynomial;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_input_polynomial");
        if (!input.get_fmpq_poly(flint::FmpqPolyRef(element_polynomial))) {
            return false;
        }
    }

    const nf_struct* raw_field = parent->raw_flint_field();
    if (raw_field == nullptr ||
        !fmpq_poly_is_monic_integral(raw_field->pol)) {
        return false;
    }

    // The lifted candidate is reconstructed as f'(theta) c in Z[theta], which
    // holds for every integral c.  Following reference
    // `src/NumField/NfAbs/Elem.jl:is_power(::AbsSimpleNumFieldElem, ::Int)`,
    // a non-integral input a is replaced by a d^n, with d the denominator of
    // a in the power basis, and a root y of y^n = a d^n gives the root y / d
    // of a.  Since (y / d)^n = a exactly when y^n = a d^n, a residue disproof
    // for a d^n disproves a.  Silex rescales only when a is not an algebraic
    // integer, since otherwise every root is integral and already
    // reconstructible.  The scaled radicand and scale exist only on the
    // rescaled branch, so the integral path (the common case) searches the
    // input itself with no extra Element construction, copy, or
    // multiply-by-one (a null root_scale means "scale by 1").
    if (fmpz_is_one(element_polynomial.raw()->den) == 0) {
        Element scaled_radicand;
        Element scale;
        bool rescaled = false;
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.power_hensel_denominator");
            flint::Fmpz leading;
            if (!primitive_characteristic_leading_coefficient(leading,
                                                              input)) {
                return false;
            }
            if (fmpz_is_one(leading.raw()) == 0) {
                // Silex pre-filter, not in the reference: it only returns
                // unsupported, so it never changes a definite answer.  With
                // lc(P) as in primitive_characteristic_leading_coefficient,
                // lc(P_a) = lc(P_c)^n when c^n = a, because the finite-place
                // terms d_v log max(1, |c|_v) scale by n.  A non-integral a
                // has only non-integral roots, whose lc is at least 2.  So
                // lc(P_a) must be an exact n-th power of an integer >= 2,
                // which also bounds n by log_2 lc(P_a) and hence the size of
                // d^n.
                flint::Fmpz lc_root;
                flint::Fmpz one;
                fmpz_one(one.raw());
                if (fmpz_cmp_ui(leading.raw(),
                                static_cast<ulong>(1)) <= 0 ||
                    static_cast<ulong>(exponent) >= fmpz_bits(leading.raw()) ||
                    fmpz_root(lc_root.raw(), leading.raw(), exponent) == 0) {
                    return false;
                }

                flint::Fmpq inverse_denominator;
                fmpq_set_fmpz_frac(inverse_denominator.raw(), one.raw(),
                                   element_polynomial.raw()->den);
                flint::FmpqPoly inverse_scale;
                fmpq_poly_set_coeff_fmpq(inverse_scale.raw(), 0,
                                         inverse_denominator.raw());

                flint::Fmpz scale_by;
                fmpz_pow_ui(scale_by.raw(), element_polynomial.raw()->den,
                            static_cast<ulong>(exponent));
                fmpq_poly_scalar_mul_fmpz(element_polynomial.raw(),
                                          element_polynomial.raw(),
                                          scale_by.raw());
                if (!scale.define(*parent) ||
                    !scaled_radicand.define(*parent) ||
                    !scale.set_fmpq_poly(
                            flint::FmpqPolyConstRef(inverse_scale)) ||
                    !scaled_radicand.set_fmpq_poly(
                            flint::FmpqPolyConstRef(element_polynomial))) {
                    return false;
                }
                rescaled = true;
            }
        }
        if (rescaled) {
            return pure_power_hensel_prime_search(
                    is_power, root, scaled_radicand, input, &scale,
                    element_polynomial, parent, raw_field, exponent,
                    pre_lift, diagnostics);
        }
    }

    return pure_power_hensel_prime_search(
            is_power, root, input, input, nullptr, element_polynomial,
            parent, raw_field, exponent, pre_lift, diagnostics);
}

bool pure_square_hensel_root(bool& is_square,
                                   Element& root,
                                   const Element& input,
                                   const DiagnosticsContext* diagnostics) noexcept {
    return pure_power_hensel_root(is_square, root, input, 2,
                                        diagnostics);
}

bool pure_power_residue_disproves(bool& is_power,
                                        const Element& input,
                                        slong exponent,
                                        const DiagnosticsContext* diagnostics) noexcept;

bool pure_square_residue_disproves(bool& is_square,
                                         const Element& input,
                                         const DiagnosticsContext* diagnostics) noexcept {
    // reference `roots(y^2 - a, ispure=true, is_normal=true)` returns no roots
    // before Hensel lifting when a good residue field already has no root.
    // This helper ports only that exact disproof stage.
    return pure_power_residue_disproves(is_square, input, 2,
                                              diagnostics);
}

bool pure_power_residue_disproves(bool& is_power,
                                        const Element& input,
                                        slong exponent,
                                        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.power_residue_disproof");
    const NumberField* parent = input.parent();
    if (parent == nullptr || parent->degree() <= 1 || exponent <= 1) {
        return false;
    }

    flint::FmpqPoly element_polynomial;
    if (!input.get_fmpq_poly(flint::FmpqPolyRef(element_polynomial))) {
        return false;
    }

    const nf_struct* raw_field = parent->raw_flint_field();
    if (raw_field == nullptr) {
        return false;
    }

    flint::Fmpz prime;
    fmpz_set_ui(prime.raw(), 3);
    // This bound only limits an optional disproof search.  Exhaustion leaves
    // the caller in the existing unsupported/fail-closed state.
    constexpr slong max_prime_attempts = 256;

    if (parent->degree() >= 10) {
        // Degree 10 and higher, where pure_power_hensel_root is not ported.
        // reference v0.38.6 `_roots_hensel` (Hensel.jl lines 159-161)
        // returns no roots as soon as y^n - a has no root modulo a prime
        // factor of a good prime, and in degree >= 10 it examines further
        // factors before lifting (STABILIZED = 2, lines 102-106).  Silex
        // applies that residue check, as the exact membership test of
        // pre_lift_residue_disproves (reference Fq_ispower) under the same
        // admissibility conditions (1)-(3) stated above
        // kPreLiftDisproofExtraPrimes (monic integral f; Dedekind-Kummer,
        // Cohen GTM 138 Theorem 4.8.13), to every factor of every admissible
        // prime: a residue that is not an n-th power proves that `input` is
        // not an n-th power.  The argument has no condition on the degree.
        // Silex extension: the reference stops after a few factors and then
        // lifts; Silex has no lift here and scans up to
        // kHighDegreeDisproofPrimes admissible primes instead.  A true power
        // is never disproved, so it stays unsupported, as does any input with
        // no disproving residue within that bound.
        if (!fmpq_poly_is_monic_integral(raw_field->pol)) {
            return false;
        }
        slong admissible_primes = 0;
        for (slong attempt = 0; attempt < max_prime_attempts &&
                                admissible_primes < kHighDegreeDisproofPrimes;
             ++attempt) {
            fmpz_nextprime(prime.raw(), prime.raw(), 1);
            // Exact skip (see residue_degree_can_disprove_power): no residue
            // above this prime can disprove the input, so it is not
            // factored and does not count towards the bound.
            if (fmpz_abs_fits_ui(prime.raw()) == 0 ||
                !residue_degree_can_disprove_power(fmpz_get_ui(prime.raw()),
                                                   parent->degree(),
                                                   exponent)) {
                continue;
            }
            bool admissible = false;
            if (residue_disproves_power_at_prime(
                        admissible, raw_field->pol, element_polynomial.raw(),
                        parent->degree(), prime.raw(), exponent)) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::element,
                        "element.power_residue_disproof.high_degree_nonpower");
                is_power = false;
                return true;
            }
            if (admissible) {
                ++admissible_primes;
            }
        }
        return false;
    }

    for (slong attempt = 0; attempt < max_prime_attempts; ++attempt) {
        fmpz_nextprime(prime.raw(), prime.raw(), 1);
        flint::FmpzModCtx ctx(prime.raw());
        if (ctx.raw() == nullptr) {
            return false;
        }

        flint::FmpzModPoly field_polynomial(ctx);
        flint::FmpzModPoly element(ctx);
        if (!field_polynomial.is_initialized() || !element.is_initialized() ||
            !fmpq_poly_reduce_modulus(field_polynomial, raw_field->pol, ctx,
                                      prime.raw()) ||
            !fmpq_poly_reduce_modulus(element, element_polynomial.raw(), ctx,
                                      prime.raw()) ||
            fmpz_mod_poly_degree(field_polynomial.raw(), ctx.raw()) !=
                    parent->degree()) {
            continue;
        }

        flint::FmpzModPolyFactor factorization(ctx);
        fmpz_mod_poly_factor(factorization.raw(), field_polynomial.raw(),
                             ctx.raw());
        if (!fmpz_mod_poly_factorization_is_squarefree(factorization)) {
            continue;
        }

        for (slong i = 0; i < factorization.raw()->num; ++i) {
            flint::FmpzModPoly factor(ctx);
            if (!factor.is_initialized()) {
                return false;
            }
            fmpz_mod_poly_factor_get_poly(factor.raw(), factorization.raw(),
                                          i, ctx.raw());
            std::vector<flint::FmpzPoly> inverse_roots;
            const ResiduePowerStatus status =
                    residue_field_power_inverse_roots(
                            inverse_roots, element, factor, ctx, prime.raw(),
                            exponent, diagnostics);
            if (status == ResiduePowerStatus::nonpower) {
                is_power = false;
                return true;
            }
            if (status == ResiduePowerStatus::power) {
                return false;
            }
        }
    }

    return false;
}

bool cyclotomic_root_hensel(Element& root,
                            const NumberField& field,
                            ulong m) noexcept {
    const nf_struct* raw_field = field.raw_flint_field();
    const slong degree = field.degree();
    if (!field.is_defined() || raw_field == nullptr || m < 3 ||
        degree < 2 || !fmpq_poly_is_monic_integral(raw_field->pol) ||
        !root.has_parent(field)) {
        return false;
    }

    flint::FmpzPoly f;
    fmpz_poly_cyclotomic(f.raw(), m);
    CyclotomicResidueRoot residue;
    if (select_cyclotomic_prime(residue, f, raw_field, degree) !=
        CyclotomicPrimeStatus::selected) {
        return false;
    }
    Element candidate(field);
    if (!candidate.is_defined() ||
        !lift_cyclotomic_root(candidate, field, f, m, residue)) {
        return false;
    }
    return root.set(candidate);
}

}  // namespace detail

}  // namespace silex
