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

#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

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

bool equation_order_norm_change_constants(double& c1,
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

    EmbeddingContext embeddings(*parent);
    if (!embeddings.refine(precision)) {
        return false;
    }

    Element derivative(*parent);
    if (!derivative.is_defined() ||
        !derivative.set_fmpq_poly(flint::FmpqPolyConstRef(
                derivative_polynomial))) {
        return false;
    }

    const slong degree = parent->degree();
    flint::AcbVec input_values(degree);
    flint::AcbVec derivative_values(degree);
    if (!embeddings.evaluate_all(flint::AcbVecRef(input_values), input,
                                 precision) ||
        !embeddings.evaluate_all(flint::AcbVecRef(derivative_values),
                                 derivative, precision)) {
        return false;
    }

    const Signature sig = embeddings.signature();
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
    // Source trace: reference implementation
    // `src/NumFieldOrd/NfOrd/Hensel.jl:_lifting_expo` with
    // `any_order(K)` equal to the equation order for this monic-integral
    // private root path.
    const NumberField* parent = input.parent();
    if (parent == nullptr || residue_degree <= 0 || exponent <= 1 ||
        fmpz_fits_si(prime.raw()) == 0 || fmpz_sgn(prime.raw()) <= 0) {
        return false;
    }

    constexpr slong precision = 128;
    std::vector<flint::Arb> root_bounds;
    if (!pure_power_conjugate_bounds(
                root_bounds, input, derivative_polynomial, exponent,
                precision)) {
        return false;
    }

    double c1 = 0.0;
    double c2 = 0.0;
    if (!equation_order_norm_change_constants(
                c1, c2, *parent, precision)) {
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
        fmpz_mod_poly_degree(modulus.raw(), ctx.raw()) <= 0 ||
        fmpz_fdiv_ui(prime, static_cast<ulong>(exponent)) == 0) {
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
    // reference gen_Shanks_sqrtn only when gcd(n, #F_q - 1) = 1.
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
    if (fq_unique_power_inverse_root(
                inverse_roots, value, field, power_data, modulus, ctx,
                diagnostics)) {
        return ResiduePowerStatus::power;
    }

    fq_neg(negative_value.raw(), value.raw(), field.raw());
    fq_one(one.raw(), field.raw());

    FqPoly polynomial(field.raw());
    fq_poly_set_coeff(polynomial.raw(), 0, negative_value.raw(), field.raw());
    fq_poly_set_coeff(polynomial.raw(), exponent, one.raw(), field.raw());

    FqPolyFactor roots(field.raw());
    fq_poly_roots(roots.raw(), polynomial.raw(), 0, field.raw());
    if (roots.raw()->num == 0) {
        return ResiduePowerStatus::nonpower;
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
        flint::fmpz_mat_ncols(inverse_num) != degree ||
        !get_basis(basis, precision_modulus, degree, lifted_factor,
                         ctx)) {
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

bool try_pure_power_root_at_prime(
        Element& root,
        const Element& input,
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
        fmpz_fdiv_ui(prime.raw(), static_cast<ulong>(exponent)) == 0) {
        return false;
    }

    if (exponent == 2 &&
        fmpz_mod_poly_degree(factor.raw(), prime_ctx.raw()) ==
                parent->degree() &&
        try_full_degree_pure_square_root_at_prime(
                root, input, element_polynomial, raw_field, prime, factor,
                prime_ctx, diagnostics)) {
        return true;
    }

    if (inverse_roots.empty()) {
        return false;
    }

    flint::FmpqPoly derivative_polynomial;
    Element denominator(*parent);
    Element inverse_denominator(*parent);
    flint::Fmpz exponent_fmpz;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.power_hensel_prime_denominator");
        fmpq_poly_derivative(derivative_polynomial.raw(), raw_field->pol);
        fmpz_set_si(exponent_fmpz.raw(), exponent);
        if (!denominator.is_defined() || !inverse_denominator.is_defined() ||
            !denominator.set_fmpq_poly(
                    flint::FmpqPolyConstRef(derivative_polynomial)) ||
            !inverse_denominator.invert(denominator)) {
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
                            candidate, input,
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

bool pure_power_hensel_root(bool& is_power,
                                  Element& root,
                                  const Element& input,
                                  slong exponent,
                                  const DiagnosticsContext* diagnostics) noexcept {
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

    flint::Fmpz prime;
    fmpz_set_ui(prime.raw(), 3);
    constexpr slong max_prime_attempts = 256;
    for (slong attempt = 0; attempt < max_prime_attempts; ++attempt) {
        fmpz_nextprime(prime.raw(), prime.raw(), 1);
        if (fmpz_fdiv_ui(prime.raw(), static_cast<ulong>(exponent)) == 0) {
            continue;
        }
        flint::FmpzModCtx ctx(prime.raw());
        if (ctx.raw() == nullptr) {
            return false;
        }

        flint::FmpzModPoly field_polynomial(ctx);
        flint::FmpzModPoly element(ctx);
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.power_hensel_prime_project");
            if (!field_polynomial.is_initialized() ||
                !element.is_initialized() ||
                !fmpq_poly_reduce_modulus(field_polynomial, raw_field->pol,
                                          ctx, prime.raw()) ||
                !fmpq_poly_reduce_modulus(element,
                                          element_polynomial.raw(), ctx,
                                          prime.raw()) ||
                fmpz_mod_poly_degree(field_polynomial.raw(), ctx.raw()) !=
                        parent->degree()) {
                continue;
            }
        }

        flint::FmpzModPolyFactor factorization(ctx);
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.power_hensel_field_factor");
            fmpz_mod_poly_factor(factorization.raw(), field_polynomial.raw(),
                                 ctx.raw());
            if (!fmpz_mod_poly_factorization_is_squarefree(factorization)) {
                continue;
            }
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
                if (try_pure_power_root_at_prime(
                            root, input, element_polynomial, raw_field, prime,
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
    if (parent == nullptr || parent->degree() <= 1 ||
        parent->degree() >= 10 || exponent <= 1) {
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
    for (slong attempt = 0; attempt < max_prime_attempts; ++attempt) {
        fmpz_nextprime(prime.raw(), prime.raw(), 1);
        if (fmpz_fdiv_ui(prime.raw(), static_cast<ulong>(exponent)) == 0) {
            continue;
        }
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

}  // namespace

namespace detail {

bool ensure_parent(Element& out, const NumberField& field) noexcept {
    if (out.has_parent(field)) {
        return true;
    }
    if (out.parent() != nullptr) {
        return false;
    }
    out = Element(field);
    return out.is_defined();
}

}  // namespace detail

Element::Element(const NumberField& parent) noexcept {
    define(parent);
}

Element::~Element() noexcept = default;

Element::Element(Element&& other) noexcept {
    swap(other);
}

Element& Element::operator=(Element&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void Element::swap(Element& other) noexcept {
    parent_.swap(other.parent_);
    value_.swap(other.value_);
}

void Element::clear() noexcept {
    value_.clear();
    parent_.clear();
}

bool Element::define(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    Element next;
    next.parent_ = parent;
    next.value_ = flint::NfElem(parent.raw_flint_field());
    if (!next.value_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool Element::set(const Element& other) noexcept {
    if (!has_same_parent(other)) {
        return false;
    }
    if (this == &other) {
        return true;
    }
    nf_elem_set(value_.raw(), other.value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::is_defined() const noexcept {
    return parent_.is_defined() && value_.is_defined();
}

const NumberField* Element::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

bool Element::has_parent(const NumberField& parent) const noexcept {
    return is_defined() && parent_.has_same_data(parent);
}

bool Element::has_same_parent(const Element& other) const noexcept {
    return is_defined() && other.is_defined() &&
           parent_.has_same_data(other.parent_);
}

bool Element::zero() noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_zero(value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::one() noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_one(value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::gen() noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_gen(value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::set_si(slong value) noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_set_si(value_.raw(), value, parent_.raw_flint_field());
    return true;
}

bool Element::set_fmpz(flint::FmpzConstRef value) noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_set_fmpz(value_.raw(), value.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::set_si_over_si(slong numerator, slong denominator) noexcept {
    if (!is_defined() || denominator == 0) {
        return false;
    }

    Element tmp(parent_);
    nf_elem_set_si(tmp.value_.raw(), numerator, parent_.raw_flint_field());
    nf_elem_scalar_div_si(tmp.value_.raw(), tmp.value_.raw(), denominator,
                          parent_.raw_flint_field());
    swap(tmp);
    return true;
}

bool Element::set_fmpq_poly(flint::FmpqPolyConstRef polynomial) noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_set_fmpq_poly(value_.raw(), polynomial.raw(),
                          parent_.raw_flint_field());
    return true;
}

bool Element::get_fmpq_poly(flint::FmpqPolyRef polynomial) const noexcept {
    if (!is_defined()) {
        return false;
    }
    nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(),
                          parent_.raw_flint_field());
    return true;
}

std::optional<flint::FmpqPoly> Element::to_fmpq_poly() const noexcept {
    flint::FmpqPoly polynomial;
    if (!get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return std::nullopt;
    }
    return polynomial;
}

bool Element::equal(const Element& other) const noexcept {
    if (!has_same_parent(other)) {
        return false;
    }
    return nf_elem_equal(value_.raw(), other.value_.raw(),
                         parent_.raw_flint_field()) != 0;
}

bool Element::equal_si(slong value) const noexcept {
    if (!is_defined()) {
        return false;
    }
    return nf_elem_equal_si(value_.raw(), value,
                            parent_.raw_flint_field()) != 0;
}

bool Element::negate(const Element& input) noexcept {
    if (!has_same_parent(input)) {
        return false;
    }
    nf_elem_neg(value_.raw(), input.value_.raw(), parent_.raw_flint_field());
    return true;
}

bool Element::add(const Element& left, const Element& right) noexcept {
    if (!is_defined() || !left.has_same_parent(right) ||
        !parent_.has_same_data(left.parent_)) {
        return false;
    }
    nf_elem_add(value_.raw(), left.value_.raw(), right.value_.raw(),
                parent_.raw_flint_field());
    return true;
}

bool Element::add_si(const Element& input, slong value) noexcept {
    if (!has_same_parent(input)) {
        return false;
    }
    nf_elem_add_si(value_.raw(), input.value_.raw(), value,
                   parent_.raw_flint_field());
    return true;
}

bool Element::subtract(const Element& left, const Element& right) noexcept {
    if (!is_defined() || !left.has_same_parent(right) ||
        !parent_.has_same_data(left.parent_)) {
        return false;
    }
    nf_elem_sub(value_.raw(), left.value_.raw(), right.value_.raw(),
                parent_.raw_flint_field());
    return true;
}

bool Element::multiply(const Element& left, const Element& right) noexcept {
    if (!is_defined() || !left.has_same_parent(right) ||
        !parent_.has_same_data(left.parent_)) {
        return false;
    }
    nf_elem_mul(value_.raw(), left.value_.raw(), right.value_.raw(),
                parent_.raw_flint_field());
    return true;
}

bool Element::scalar_div_si(const Element& input,
                            slong denominator) noexcept {
    if (!has_same_parent(input) || denominator == 0) {
        return false;
    }

    Element tmp(parent_);
    nf_elem_scalar_div_si(tmp.value_.raw(), input.value_.raw(), denominator,
                          parent_.raw_flint_field());
    swap(tmp);
    return true;
}

bool Element::invert(const Element& input) noexcept {
    if (!has_same_parent(input) || input.equal_si(0)) {
        return false;
    }

    Element tmp(parent_);
    nf_elem_inv(tmp.value_.raw(), input.value_.raw(),
                parent_.raw_flint_field());
    swap(tmp);
    return true;
}

bool Element::pow_fmpz(const Element& input,
                       flint::FmpzConstRef exponent) noexcept {
    if (!has_same_parent(input)) {
        return false;
    }
    const int exponent_sign = fmpz_sgn(exponent.raw());
    if (exponent_sign < 0 && input.equal_si(0)) {
        return false;
    }

    const nf_struct* raw_field = parent_.raw_flint_field();
    if (flint::fmpz_abs_fits_ui(exponent)) {
        if (exponent_sign >= 0) {
            Element tmp(parent_);
            nf_elem_pow(tmp.value_.raw(), input.value_.raw(),
                        flint::fmpz_get_ui(exponent), raw_field);
            swap(tmp);
            return true;
        }

        Element inverse(parent_);
        Element tmp(parent_);
        flint::Fmpz abs_exponent;
        if (!inverse.invert(input)) {
            return false;
        }
        flint::fmpz_neg(flint::FmpzRef(abs_exponent), exponent);
        nf_elem_pow(tmp.value_.raw(), inverse.value_.raw(),
                    flint::fmpz_get_ui(flint::FmpzConstRef(abs_exponent)),
                    raw_field);
        swap(tmp);
        return true;
    }

    Element accumulator(parent_);
    Element base(parent_);
    Element product(parent_);
    flint::Fmpz k;

    accumulator.one();
    if (exponent_sign < 0) {
        if (!base.invert(input)) {
            return false;
        }
        fmpz_neg(k.raw(), exponent.raw());
    } else {
        if (!base.set(input)) {
            return false;
        }
        fmpz_set(k.raw(), exponent.raw());
    }

    while (fmpz_is_zero(k.raw()) == 0) {
        if (fmpz_is_odd(k.raw()) != 0) {
            if (!product.multiply(accumulator, base) ||
                !accumulator.set(product)) {
                return false;
            }
        }

        fmpz_fdiv_q_2exp(k.raw(), k.raw(), 1);
        if (fmpz_is_zero(k.raw()) == 0) {
            if (!product.multiply(base, base) || !base.set(product)) {
                return false;
            }
        }
    }

    swap(accumulator);
    return true;
}

bool Element::is_square(
        bool& is_square,
        Element& root,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.is_square");
    if (!has_same_parent(root)) {
        return false;
    }

    Element candidate(parent_);
    bool flag = false;
    bool known = is_square_rational_constant(flag, candidate, *this);

    if (known && !flag && parent_.degree() != 1) {
        known = false;
    }

    if (!known && parent_.degree() != 1 &&
        parent_.backend_kind() == NumberFieldBackendKind::quadratic) {
        flint::Fmpz radicand;
        if (parent_.quadratic_radicand(flint::FmpzRef(radicand))) {
            known = is_square_quadratic(flag, candidate, *this, radicand.raw());
        }
    }

    if (!known &&
        pure_square_hensel_root(flag, candidate, *this, diagnostics)) {
        known = true;
    }

    if (!known &&
        pure_square_residue_disproves(flag, *this, diagnostics)) {
        known = true;
    }

    if (!known) {
        return false;
    }

    is_square = flag;
    if (flag) {
        root.swap(candidate);
    }
    return true;
}

bool Element::is_power(bool& is_power,
                       Element& root,
                       flint::FmpzConstRef exponent,
                       const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.is_power");
    if (!has_same_parent(root) || fmpz_cmp_ui(exponent.raw(), 1) < 0) {
        return false;
    }

    if (fmpz_is_one(exponent.raw()) != 0) {
        Element candidate(parent_);
        if (!candidate.set(*this)) {
            return false;
        }
        is_power = true;
        root.swap(candidate);
        return true;
    }

    Element candidate(parent_);
    bool flag = false;
    if (is_power_rational_constant(flag, candidate, *this, exponent)) {
        if (flag) {
            is_power = true;
            root.swap(candidate);
            return true;
        }
        if (parent_.degree() == 1) {
            is_power = false;
            return true;
        }
    }

    if (fmpz_equal_ui(exponent.raw(), 2) != 0) {
        return is_square(is_power, root, diagnostics);
    }

    if (fmpz_fits_si(exponent.raw()) != 0) {
        const slong exponent_si = fmpz_get_si(exponent.raw());
        if (exponent_si > 1) {
            if (pure_power_hensel_root(
                        flag, candidate, *this, exponent_si, diagnostics)) {
                is_power = flag;
                if (flag) {
                    root.swap(candidate);
                }
                return true;
            }
            if (pure_power_residue_disproves(
                        flag, *this, exponent_si, diagnostics)) {
                is_power = flag;
                return true;
            }
        }
    }

    return false;
}

bool Element::trace(flint::FmpqRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }

    const nf_struct* field = parent_.raw_flint_field();
    if (parent_.backend_kind() == NumberFieldBackendKind::quadratic) {
        flint::FmpqPoly polynomial;
        flint::Fmpq constant;

        nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(), field);
        fmpq_poly_get_coeff_fmpq(constant.raw(), polynomial.raw(), 0);
        fmpq_mul_2exp(out.raw(), constant.raw(), 1);
        return true;
    }

    generic_trace(out, field, parent_.degree(), value_.raw());
    return true;
}

bool Element::norm(flint::FmpqRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }

    const nf_struct* field = parent_.raw_flint_field();
    flint::Fmpz radicand;
    if (parent_.backend_kind() == NumberFieldBackendKind::quadratic &&
        parent_.quadratic_radicand(flint::FmpzRef(radicand))) {
        flint::FmpqPoly polynomial;
        flint::Fmpq constant;
        flint::Fmpq linear;
        flint::Fmpq term;

        nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(), field);
        fmpq_poly_get_coeff_fmpq(constant.raw(), polynomial.raw(), 0);
        fmpq_poly_get_coeff_fmpq(linear.raw(), polynomial.raw(), 1);

        fmpq_mul(out.raw(), constant.raw(), constant.raw());
        fmpq_mul(term.raw(), linear.raw(), linear.raw());
        fmpq_mul_fmpz(term.raw(), term.raw(), radicand.raw());
        fmpq_sub(out.raw(), out.raw(), term.raw());
        return true;
    }

    generic_norm(out, field, parent_.degree(), value_.raw());
    return true;
}

bool Element::conjugate(Element& out) const noexcept {
    if (!has_same_parent(out) ||
        parent_.backend_kind() != NumberFieldBackendKind::quadratic) {
        return false;
    }

    const nf_struct* field = parent_.raw_flint_field();
    flint::FmpqPoly polynomial;
    flint::Fmpq linear;

    nf_elem_get_fmpq_poly(polynomial.raw(), value_.raw(), field);
    fmpq_poly_get_coeff_fmpq(linear.raw(), polynomial.raw(), 1);
    fmpq_neg(linear.raw(), linear.raw());
    fmpq_poly_set_coeff_fmpq(polynomial.raw(), 1, linear.raw());
    nf_elem_set_fmpq_poly(out.value_.raw(), polynomial.raw(), field);
    return true;
}

flint::NfElemRef Element::flint_element_ref() noexcept {
    return flint::NfElemRef(value_);
}

flint::NfElemConstRef Element::flint_element_ref() const noexcept {
    return flint::NfElemConstRef(value_);
}

nf_elem_struct* Element::raw_flint_element() noexcept {
    return value_.raw();
}

const nf_elem_struct* Element::raw_flint_element() const noexcept {
    return value_.raw();
}

}  // namespace silex
