#include "compact_reconstruction_internal.hpp"
#include "wide_factored_internal.hpp"
#include "../unit/unit_internal.hpp"

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/nmod_mat.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/flint/nmod_poly_factor.hpp>
#include <silex/order_element.hpp>

#include <flint/fmpz_mat.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace silex::detail {
namespace {

struct PreparedFactor {
    PreparedFactor(const Element& representative_value,
                   std::size_t output_count) noexcept
        : representative(&representative_value), exponents(output_count) {}

    const Element* representative = nullptr;
    // Arbitrary-size exponents: reference FlxqV_factorback takes t_INT
    // exponents (polarit2.c:1275-1298), so a wide exponent matrix
    // reconstructs too.
    std::vector<flint::Fmpz> exponents;
    flint::FmpqPoly power_polynomial;
    bool needs_inverse = false;
};

void reset_report(BoundedCompactReconstructionReport& report) noexcept {
    report.status = BoundedCompactReconstructionStatus::invalid_input;
    report.factor_rows = 0;
    report.prime_trials = 0;
    report.primes_used = 0;
    report.field_polynomial_prime_rejections = 0;
    report.order_basis_prime_rejections = 0;
    report.denominator_prime_rejections = 0;
    report.noninvertible_prime_rejections = 0;
    report.nonsplit_prime_skips = 0;
    report.modulus_bits = 0;
    report.last_prime = 0;
    fmpz_one(report.centered_crt_modulus.raw());
}

bool checked_add_slong(slong& out, slong value) noexcept {
    const slong min = std::numeric_limits<slong>::min();
    const slong max = std::numeric_limits<slong>::max();
    if ((value > 0 && out > max - value) || (value < 0 && out < min - value)) {
        return false;
    }
    out += value;
    return true;
}

bool checked_add_exponent(flint::Fmpz& out, slong value) noexcept {
    // Word-size inputs keep their checked slong sum.
    slong sum = fmpz_get_si(out.raw());
    if (!checked_add_slong(sum, value)) {
        return false;
    }
    fmpz_set_si(out.raw(), sum);
    return true;
}

bool exponent_row_is_zero(const PreparedFactor& factor) noexcept {
    return std::all_of(factor.exponents.begin(), factor.exponents.end(),
                       [](const flint::Fmpz& exponent) {
                           return fmpz_is_zero(exponent.raw()) != 0;
                       });
}

bool prepare_exponent_matrix(BoundedCompactReconstructionReport& report,
                             std::vector<PreparedFactor>& factors,
                             const Order& order,
                             std::span<const FactoredElement> inputs) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    for (std::size_t column = 0; column < inputs.size(); ++column) {
        const FactoredElement& input = inputs[column];
        if (!input.is_defined() || input.parent() == nullptr ||
            !input.parent()->has_same_data(*field)) {
            return false;
        }

        for (const auto& entry : input.factors()) {
            if (entry.exponent == 0) {
                continue;
            }
            if (!entry.factor.is_defined() ||
                !entry.factor.has_parent(*field) ||
                (entry.exponent < 0 && entry.factor.equal_si(0))) {
                return false;
            }

            auto match = std::find_if(
                    factors.begin(), factors.end(),
                    [&entry](const PreparedFactor& factor) {
                        return factor.representative != nullptr &&
                               factor.representative->equal(entry.factor);
                    });
            if (match == factors.end()) {
                factors.emplace_back(entry.factor, inputs.size());
                match = std::prev(factors.end());
            }
            if (!checked_add_exponent(match->exponents[column],
                                      entry.exponent)) {
                report.status =
                        BoundedCompactReconstructionStatus::exponent_overflow;
                return false;
            }
        }
    }

    factors.erase(std::remove_if(factors.begin(), factors.end(),
                                 exponent_row_is_zero),
                  factors.end());
    report.factor_rows = factors.size();
    return true;
}

bool prepare_wide_exponent_matrix(BoundedCompactReconstructionReport& report,
                                  std::vector<PreparedFactor>& factors,
                                  const Order& order,
                                  const WideFactoredUnits& inputs) noexcept {
    const NumberField* field = order.parent();
    const slong outputs = flint::fmpz_mat_nrows(inputs.exponents);
    if (field == nullptr || outputs < 0 ||
        flint::fmpz_mat_ncols(inputs.exponents) !=
                static_cast<slong>(inputs.bases.size())) {
        return false;
    }

    for (std::size_t column = 0; column < inputs.bases.size(); ++column) {
        const Element& base = inputs.bases[column];
        if (!base.is_defined() || !base.has_parent(*field)) {
            return false;
        }
        factors.emplace_back(base, static_cast<std::size_t>(outputs));
        PreparedFactor& factor = factors.back();
        for (slong row = 0; row < outputs; ++row) {
            const fmpz* exponent = fmpz_mat_entry(
                    inputs.exponents.raw(), row, static_cast<slong>(column));
            if (fmpz_sgn(exponent) < 0 && base.equal_si(0)) {
                return false;
            }
            fmpz_set(factor.exponents[static_cast<std::size_t>(row)].raw(),
                     exponent);
        }
    }

    factors.erase(std::remove_if(factors.begin(), factors.end(),
                                 exponent_row_is_zero),
                  factors.end());
    report.factor_rows = factors.size();
    return true;
}

bool prepare_power_polynomials(std::vector<PreparedFactor>& factors) noexcept {
    for (PreparedFactor& factor : factors) {
        if (factor.representative == nullptr ||
            !factor.representative->get_fmpq_poly(
                    flint::FmpqPolyRef(factor.power_polynomial))) {
            return false;
        }
        factor.needs_inverse = std::any_of(
                factor.exponents.begin(), factor.exponents.end(),
                [](const flint::Fmpz& exponent) {
                    return fmpz_sgn(exponent.raw()) < 0;
                });
    }
    return true;
}

bool reduce_fmpq(ulong& out, flint::FmpqConstRef value, ulong prime,
                 nmod_t modulus) noexcept {
    const ulong denominator =
            fmpz_fdiv_ui(flint::fmpq_den_ref(value).raw(), prime);
    if (denominator == 0) {
        return false;
    }
    const ulong numerator =
            fmpz_fdiv_ui(flint::fmpq_num_ref(value).raw(), prime);
    out = nmod_mul(numerator, n_invmod(denominator, prime), modulus);
    return true;
}

bool reduce_power_polynomial(flint::NmodPoly& out, const flint::FmpqPoly& input,
                             slong maximum_degree, ulong prime) noexcept {
    flint::Fmpq coefficient;
    const nmod_t modulus = out.raw()->mod;
    nmod_poly_zero(out.raw());
    for (slong i = 0; i <= maximum_degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coefficient.raw(), input.raw(), i);
        ulong reduced = 0;
        if (!reduce_fmpq(reduced, flint::FmpqConstRef(coefficient), prime,
                         modulus)) {
            return false;
        }
        nmod_poly_set_coeff_ui(out.raw(), i, reduced);
    }
    return true;
}

bool reduce_order_basis(flint::NmodMat& out, const flint::FmpqMat& basis,
                        ulong prime) noexcept {
    const slong degree = nmod_mat_nrows(out.raw());
    const nmod_t modulus = out.raw()->mod;
    for (slong row = 0; row < degree; ++row) {
        for (slong column = 0; column < degree; ++column) {
            ulong reduced = 0;
            if (!reduce_fmpq(reduced, flint::fmpq_mat_entry(basis, row, column),
                             prime, modulus)) {
                return false;
            }
            nmod_mat_entry(out.raw(), row, column) = reduced;
        }
    }
    return true;
}

bool reduce_factor_polynomials(std::vector<flint::NmodPoly>& bases,
                               const std::vector<PreparedFactor>& factors,
                               slong degree, ulong prime) noexcept {
    for (std::size_t row = 0; row < factors.size(); ++row) {
        if (!reduce_power_polynomial(bases[row], factors[row].power_polynomial,
                                     degree - 1, prime)) {
            return false;
        }
    }
    return true;
}

// Odd powers base^1, base^3, ..., base^(2^width - 1) of one base (or of its
// inverse) in the ring of the power products: the sliding-window table of
// Handbook of Applied Cryptography, Algorithm 14.85 step 1.
template <typename Value>
struct OddPowerTable {
    std::vector<Value> odd_powers;
    ulong width = 0;
};

// One window of a sliding-window recoding: multiply by odd_power at bit
// position `position` of the left-to-right evaluation.
template <typename Value>
struct WindowStep {
    flint_bitcnt_t position = 0;
    const Value* odd_power = nullptr;
};

// The ring (Z/pZ)[x]/(modulus), elements in the power basis.
class PolynomialPowerRing {
public:
    using Value = flint::NmodPoly;

    explicit PolynomialPowerRing(const flint::NmodPoly& modulus) noexcept
        : modulus_(modulus) {}

    Value make() const noexcept {
        return Value(modulus_.raw()->mod.n);
    }

    void set(Value& out, const Value& value) const noexcept {
        nmod_poly_set(out.raw(), value.raw());
    }

    void one(Value& out) const noexcept {
        nmod_poly_one(out.raw());
    }

    void mul(Value& out, const Value& left, const Value& right) const noexcept {
        nmod_poly_mulmod(out.raw(), left.raw(), right.raw(), modulus_.raw());
    }

    void store(flint::NmodMat& out, slong row, const Value& value) const
            noexcept {
        for (slong i = 0; i < nmod_mat_ncols(out.raw()); ++i) {
            nmod_mat_entry(out.raw(), row, i) =
                    nmod_poly_get_coeff_ui(value.raw(), i);
        }
    }

private:
    const flint::NmodPoly& modulus_;
};

// The ring (Z/pZ)^length with coordinate-wise operations.
class EvaluationPowerRing {
public:
    using Value = std::vector<ulong>;

    EvaluationPowerRing(nmod_t modulus, std::size_t length) noexcept
        : modulus_(modulus), length_(length) {}

    Value make() const noexcept {
        return Value(length_, 0);
    }

    void set(Value& out, std::span<const ulong> value) const noexcept {
        std::copy(value.begin(), value.end(), out.begin());
    }

    void one(Value& out) const noexcept {
        std::fill(out.begin(), out.end(), nmod_set_ui(1, modulus_));
    }

    void mul(Value& out, std::span<const ulong> left,
             std::span<const ulong> right) const noexcept {
        for (std::size_t i = 0; i < length_; ++i) {
            out[i] = nmod_mul(left[i], right[i], modulus_);
        }
    }

    void store(flint::NmodMat& out, slong row, const Value& value) const
            noexcept {
        for (std::size_t i = 0; i < length_; ++i) {
            nmod_mat_entry(out.raw(), row, static_cast<slong>(i)) = value[i];
        }
    }

private:
    nmod_t modulus_;
    std::size_t length_ = 0;
};

const flint::NmodPoly& factor_base(
        const ModularPowerProductFactor& factor) noexcept {
    return *factor.base;
}

const flint::NmodPoly& factor_inverse(
        const ModularPowerProductFactor& factor) noexcept {
    return *factor.inverse;
}

bool factor_has_inverse(const ModularPowerProductFactor& factor) noexcept {
    return factor.inverse != nullptr;
}

std::span<const ulong> factor_base(
        const EvaluatedPowerProductFactor& factor) noexcept {
    return factor.base;
}

std::span<const ulong> factor_inverse(
        const EvaluatedPowerProductFactor& factor) noexcept {
    return factor.inverse;
}

bool factor_has_inverse(const EvaluatedPowerProductFactor& factor) noexcept {
    return !factor.inverse.empty();
}

// Window width minimizing the operation count for a table shared by exponents
// of total bit length total_bits.  Cost model: the table costs one squaring
// and 2^(w-1) - 1 multiplications when w > 1 (Handbook of Applied
// Cryptography Table 14.16, row for Algorithm 14.83, whose precomputation is
// that of Algorithm 14.85 step 1), and a t-bit exponent needs about t / (w + 1)
// window multiplications (Moller, "Algorithms for multi-exponentiation",
// Section 3.1, expected number of multiplications).  Squarings belong to the
// interleaved evaluation and vary with w by less than w per exponent (Moller
// gives b - max w_i up to b - 1), so they are left out of the choice.  The
// width changes only the cost, never the product.
struct SlidingWindowChoice {
    ulong width = 1;
    // Modeled multiplications of the table and the windows at that width.
    double cost = 0.0;
};

SlidingWindowChoice choose_sliding_window(flint_bitcnt_t total_bits,
                                          flint_bitcnt_t maximum_bits) noexcept {
    constexpr ulong maximum_width = 16;
    const double bits = static_cast<double>(total_bits);
    SlidingWindowChoice best{1, bits / 2.0};
    for (ulong width = 2;
         width <= maximum_width && width <= static_cast<ulong>(maximum_bits);
         ++width) {
        const double cost = static_cast<double>(UWORD(1) << (width - 1)) +
                            bits / static_cast<double>(width + 1);
        if (cost < best.cost) {
            best.cost = cost;
            best.width = width;
        }
    }
    return best;
}

ulong sliding_window_width(flint_bitcnt_t total_bits,
                           flint_bitcnt_t maximum_bits) noexcept {
    return choose_sliding_window(total_bits, maximum_bits).width;
}

template <typename Ring, typename Base>
void build_odd_power_table(OddPowerTable<typename Ring::Value>& table,
                           const Ring& ring, const Base& base,
                           ulong width) noexcept {
    const std::size_t size = static_cast<std::size_t>(UWORD(1) << (width - 1));
    table.width = width;
    table.odd_powers.clear();
    table.odd_powers.reserve(size);
    table.odd_powers.emplace_back(ring.make());
    ring.set(table.odd_powers.back(), base);
    if (size == 1) {
        return;
    }
    typename Ring::Value square = ring.make();
    ring.mul(square, base, base);
    for (std::size_t i = 1; i < size; ++i) {
        table.odd_powers.emplace_back(ring.make());
        ring.mul(table.odd_powers.back(), table.odd_powers[i - 1], square);
    }
}

// Left-to-right sliding-window recoding of a positive magnitude
// (Handbook of Applied Cryptography, Algorithm 14.85 step 3): each maximal
// window of at most table.width bits that starts and ends with a one becomes
// one multiplication by the matching odd power at the window's lowest bit.
template <typename Value>
void append_window_steps(std::vector<WindowStep<Value>>& steps,
                         const flint::Fmpz& magnitude,
                         const OddPowerTable<Value>& table) noexcept {
    slong high = static_cast<slong>(fmpz_bits(magnitude.raw())) - 1;
    while (high >= 0) {
        if (fmpz_tstbit(magnitude.raw(), static_cast<ulong>(high)) == 0) {
            --high;
            continue;
        }
        slong low = std::max<slong>(high - static_cast<slong>(table.width) + 1,
                                    0);
        while (fmpz_tstbit(magnitude.raw(), static_cast<ulong>(low)) == 0) {
            ++low;
        }
        ulong value = 0;
        for (slong bit = high; bit >= low; --bit) {
            value = (value << 1) |
                    static_cast<ulong>(
                            fmpz_tstbit(magnitude.raw(),
                                        static_cast<ulong>(bit)));
        }
        steps.push_back(WindowStep<Value>{static_cast<flint_bitcnt_t>(low),
                                          &table.odd_powers[(value - 1) / 2]});
        high = low - 1;
    }
}

// The roots of a modulus f of degree `degree` that splits into distinct
// linear factors modulo p, or false when nmod_poly_roots does not return
// `degree` distinct roots.
bool distinct_roots(std::vector<ulong>& roots, const flint::NmodPoly& modulus,
                    slong degree) noexcept {
    flint::NmodPolyFactor linear_factors;
    nmod_poly_roots(linear_factors.raw(), modulus.raw(), 0);
    if (linear_factors.raw()->num != degree) {
        return false;
    }
    roots.resize(static_cast<std::size_t>(degree));
    for (slong j = 0; j < degree; ++j) {
        const nmod_poly_struct* factor = linear_factors.raw()->p + j;
        if (nmod_poly_degree(factor) != 1 ||
            nmod_poly_get_coeff_ui(factor, 1) != 1) {
            return false;
        }
        roots[static_cast<std::size_t>(j)] =
                nmod_neg(nmod_poly_get_coeff_ui(factor, 0), modulus.raw()->mod);
    }
    return true;
}

// build_prime_residues at a prime where the reduced defining polynomial f is
// the product of the x - r_j over its distinct roots r_1, ..., r_n.  By the
// Chinese remainder theorem for F_p[x] with the pairwise coprime moduli
// x - r_j, reduction modulo them, which is evaluation at the r_j, is a ring
// isomorphism F_p[x]/(f) -> F_p^n, and its inverse is Lagrange interpolation
// (von zur Gathen and Gerhard, Modern Computer Algebra, Sections 5.2 and
// 5.4).  The power products are therefore computed coordinate-wise: a ring
// multiplication is n products in F_p and an inverse is n inverses in F_p.
// An element is a unit exactly when no coordinate is zero, that is when it
// is coprime to f, which is when nmod_poly_invmod succeeds, so the same
// primes are rejected.  Interpolation is folded into the change to order
// coordinates: with the rows of B the order basis in power coordinates and
// W[i][j] = r_j^i, the element with order-coordinate row c has evaluation row
// c (B W), so c = (evaluation row) (B W)^-1.  W is invertible (a Vandermonde
// matrix of distinct points), so B W is invertible exactly when B is, and the
// same primes are rejected for the order basis.  The residues equal those of
// the power-basis route because the map is an isomorphism.
//
// Deviation from the reference: FlxqX_chinese_unit (buch2.c:1025-1043)
// computes the products in F_p[x]/(T) with FlxqV_factorback; this route
// computes the same element of that ring through the isomorphism.
bool build_split_prime_residues(BoundedCompactReconstructionReport& report,
                                flint::NmodMat& residues,
                                const std::vector<PreparedFactor>& factors,
                                const flint::FmpqMat& order_basis,
                                slong degree, ulong prime,
                                const std::vector<ulong>& roots) noexcept {
    const nmod_t modulus = residues.raw()->mod;
    const std::size_t length = static_cast<std::size_t>(degree);
    flint::NmodMat basis_mod(degree, degree, prime);
    flint::NmodMat powers(degree, degree, prime);
    for (slong j = 0; j < degree; ++j) {
        ulong power = nmod_set_ui(1, modulus);
        for (slong i = 0; i < degree; ++i) {
            nmod_mat_entry(powers.raw(), i, j) = power;
            power = nmod_mul(power, roots[static_cast<std::size_t>(j)],
                             modulus);
        }
    }
    flint::NmodMat basis_values(degree, degree, prime);
    flint::NmodMat basis_values_inverse(degree, degree, prime);
    if (!reduce_order_basis(basis_mod, order_basis, prime)) {
        ++report.order_basis_prime_rejections;
        return false;
    }
    nmod_mat_mul(basis_values.raw(), basis_mod.raw(), powers.raw());
    if (nmod_mat_inv(basis_values_inverse.raw(), basis_values.raw()) == 0) {
        ++report.order_basis_prime_rejections;
        return false;
    }

    flint::NmodPoly base(prime);
    std::vector<ulong> base_values(factors.size() * length);
    std::vector<ulong> inverse_values(factors.size() * length);
    for (std::size_t row = 0; row < factors.size(); ++row) {
        if (!reduce_power_polynomial(base, factors[row].power_polynomial,
                                     degree - 1, prime)) {
            ++report.denominator_prime_rejections;
            return false;
        }
        nmod_poly_evaluate_nmod_vec(base_values.data() + row * length,
                                    base.raw(), roots.data(), degree);
    }

    std::vector<EvaluatedPowerProductFactor> power_factors(factors.size());
    for (std::size_t row = 0; row < factors.size(); ++row) {
        const std::span<const ulong> values(base_values.data() + row * length,
                                            length);
        power_factors[row].base = values;
        power_factors[row].exponents = &factors[row].exponents;
        if (!factors[row].needs_inverse) {
            continue;
        }
        ulong* inverse = inverse_values.data() + row * length;
        for (std::size_t j = 0; j < length; ++j) {
            if (values[j] == 0) {
                ++report.noninvertible_prime_rejections;
                return false;
            }
            inverse[j] = n_invmod(values[j], prime);
        }
        power_factors[row].inverse = std::span<const ulong>(inverse, length);
    }

    flint::NmodMat power_values(nmod_mat_nrows(residues.raw()), degree, prime);
    modular_power_products_evaluated(power_values, power_factors);
    nmod_mat_mul(residues.raw(), power_values.raw(),
                 basis_values_inverse.raw());
    return true;
}

// `split` asserts that the defining polynomial reduced modulo prime splits
// into distinct linear factors; the products are then computed in evaluation
// form (build_split_prime_residues).  Otherwise, or when its roots are not
// found, they are computed in the power basis.
bool build_prime_residues(BoundedCompactReconstructionReport& report,
                          flint::NmodMat& residues,
                          const std::vector<PreparedFactor>& factors,
                          const flint::FmpqPoly& defining_polynomial,
                          const flint::FmpqMat& order_basis, slong degree,
                          ulong prime, bool split) noexcept {
    flint::NmodPoly modulus(prime);
    if (!reduce_power_polynomial(modulus, defining_polynomial, degree, prime) ||
        nmod_poly_degree(modulus.raw()) != degree ||
        nmod_poly_get_coeff_ui(modulus.raw(), degree) != 1) {
        ++report.field_polynomial_prime_rejections;
        return false;
    }
    std::vector<ulong> roots;
    if (split && distinct_roots(roots, modulus, degree)) {
        return build_split_prime_residues(report, residues, factors,
                                          order_basis, degree, prime, roots);
    }

    flint::NmodMat basis_mod(degree, degree, prime);
    flint::NmodMat basis_inverse(degree, degree, prime);
    if (!reduce_order_basis(basis_mod, order_basis, prime) ||
        nmod_mat_inv(basis_inverse.raw(), basis_mod.raw()) == 0) {
        // This is reference's equation-order/index-prime exclusion expressed
        // through the selected order basis rather than a separate index.
        ++report.order_basis_prime_rejections;
        return false;
    }

    std::vector<flint::NmodPoly> bases;
    std::vector<flint::NmodPoly> inverses;
    bases.reserve(factors.size());
    inverses.reserve(factors.size());
    for (std::size_t i = 0; i < factors.size(); ++i) {
        bases.emplace_back(prime);
        inverses.emplace_back(prime);
    }
    if (!reduce_factor_polynomials(bases, factors, degree, prime)) {
        ++report.denominator_prime_rejections;
        return false;
    }

    for (std::size_t row = 0; row < factors.size(); ++row) {
        if (factors[row].needs_inverse &&
            nmod_poly_invmod(inverses[row].raw(), bases[row].raw(),
                             modulus.raw()) == 0) {
            ++report.noninvertible_prime_rejections;
            return false;
        }
    }

    const std::size_t output_count =
            factors.empty()
                    ? static_cast<std::size_t>(nmod_mat_nrows(residues.raw()))
                    : factors.front().exponents.size();
    std::vector<ModularPowerProductFactor> power_factors(factors.size());
    for (std::size_t row = 0; row < factors.size(); ++row) {
        power_factors[row].base = &bases[row];
        power_factors[row].inverse =
                factors[row].needs_inverse ? &inverses[row] : nullptr;
        power_factors[row].exponents = &factors[row].exponents;
    }
    flint::NmodMat power_residues(static_cast<slong>(output_count), degree,
                                  prime);
    modular_power_products(power_residues, power_factors, modulus);

    // If rows of B are the selected order basis in power coordinates, then
    // power_row * B^-1 is the order-coordinate row. This is reference's invzk
    // conversion performed with FLINT's modular matrix primitive.
    nmod_mat_mul(residues.raw(), power_residues.raw(), basis_inverse.raw());
    return true;
}

bool reconstructed_coordinates_within_bound(
        const flint::FmpzMat& coordinates, flint::FmpzConstRef bound) noexcept {
    flint::Fmpz absolute_value;
    const slong rows = flint::fmpz_mat_nrows(coordinates);
    const slong columns = flint::fmpz_mat_ncols(coordinates);
    for (slong row = 0; row < rows; ++row) {
        for (slong column = 0; column < columns; ++column) {
            fmpz_abs(absolute_value.raw(),
                     fmpz_mat_entry(coordinates.raw(), row, column));
            if (fmpz_cmp(absolute_value.raw(), bound.raw()) > 0) {
                return false;
            }
        }
    }
    return true;
}

bool publish_reconstructed_elements(std::vector<Element>& out,
                                    const flint::FmpzMat& coordinates,
                                    const Order& order) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    const slong rows = flint::fmpz_mat_nrows(coordinates);
    const slong degree = order.degree();
    std::vector<Element> candidate;
    candidate.reserve(static_cast<std::size_t>(rows));
    flint::FmpzMat row(1, degree);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < degree; ++j) {
            fmpz_set(fmpz_mat_entry(row.raw(), 0, j),
                     fmpz_mat_entry(coordinates.raw(), i, j));
        }
        OrderElement order_element(order);
        candidate.emplace_back(*field);
        if (!order_element.is_defined() || !candidate.back().is_defined() ||
            !order_element.set_coordinates(flint::FmpzMatConstRef(row)) ||
            !order_element.get_element(candidate.back())) {
            return false;
        }
    }
    out.swap(candidate);
    return true;
}


// For p = modulus.mod.n prime, true exactly when the monic modulus f of
// positive degree is a product of distinct linear factors over F_p: since
// x^p - x is the product of x - a over the elements a of F_p (Lidl and
// Niederreiter, Finite Fields, Lemma 2.4), f divides x^p - x, which is
// x^p = x modulo f, exactly when f is such a product.  This is the
// degree-one step of distinct-degree factorization; one modular
// exponentiation of x decides it.
bool splits_into_distinct_linear_factors(const flint::NmodPoly& modulus)
        noexcept {
    const ulong prime = modulus.raw()->mod.n;
    const slong degree = nmod_poly_degree(modulus.raw());
    if (degree < 1 || nmod_poly_get_coeff_ui(modulus.raw(), degree) != 1) {
        return false;
    }
    flint::NmodPoly inverse(prime);
    nmod_poly_reverse(inverse.raw(), modulus.raw(), degree + 1);
    nmod_poly_inv_series(inverse.raw(), inverse.raw(), degree + 1);
    flint::NmodPoly power(prime);
    nmod_poly_powmod_x_ui_preinv(power.raw(), prime, modulus.raw(),
                                 inverse.raw());
    flint::NmodPoly x(prime);
    flint::NmodPoly x_reduced(prime);
    nmod_poly_set_coeff_ui(x.raw(), 1, 1);
    nmod_poly_rem(x_reduced.raw(), x.raw(), modulus.raw());
    return nmod_poly_equal(power.raw(), x_reduced.raw()) != 0;
}

// Bit counts of the exponents of one factor row, by sign.
struct ExponentRowBits {
    flint_bitcnt_t positive_bits = 0;
    flint_bitcnt_t positive_maximum = 0;
    flint_bitcnt_t negative_bits = 0;
    flint_bitcnt_t negative_maximum = 0;

    void add(const flint::Fmpz& exponent) noexcept {
        const flint_bitcnt_t bits = fmpz_bits(exponent.raw());
        if (fmpz_sgn(exponent.raw()) > 0) {
            positive_bits += bits;
            positive_maximum = std::max(positive_maximum, bits);
        } else if (fmpz_sgn(exponent.raw()) < 0) {
            negative_bits += bits;
            negative_maximum = std::max(negative_maximum, bits);
        }
    }

    // The window model of sliding_window_width summed over the tables the
    // row needs: one for the base when some exponent is positive and one for
    // its inverse when some exponent is negative.
    double modeled_cost() const noexcept {
        double cost = 0.0;
        if (positive_bits != 0) {
            cost += choose_sliding_window(positive_bits, positive_maximum).cost;
        }
        if (negative_bits != 0) {
            cost += choose_sliding_window(negative_bits, negative_maximum).cost;
        }
        return cost;
    }
};

// out = ((exponent - 1) mod lambda) + 1, which lies in [1, lambda].
void reduce_exponent(flint::Fmpz& out, const flint::Fmpz& exponent,
                     const flint::Fmpz& lambda) noexcept {
    fmpz_sub_ui(out.raw(), exponent.raw(), 1);
    fmpz_fdiv_r(out.raw(), out.raw(), lambda.raw());
    fmpz_add_ui(out.raw(), out.raw(), 1);
}

// Exponent rows of the factors after reduction modulo the exponent lambda of
// the unit group of R, which is (Z/pZ)[x]/(*modulus) when modulus is given
// and (Z/pZ)^n otherwise.  A row that changes is copied into storage and
// reduced there; every other row stays the caller's.
//
// Soundness.  modular_unit_group_exponent returns lambda only when p is prime
// and the modulus is monic, nonconstant and squarefree modulo p (a product of
// distinct linear factors, its lambda = p - 1 case, is squarefree).  Then R is
// the product of the fields F_q = F_p[x]/(f_i) over the irreducible factors
// f_i of the modulus (Chinese remainder theorem; Lidl and Niederreiter,
// Finite Fields, Theorem 1.61), and lambda = lcm(q - 1).  Without a modulus,
// lambda = p - 1 is used only when p is prime: R is then the product of n
// copies of the field F_p.
//  (a) e > 0, any base: in each F_q, a^q = a for every a, zero included
//      (Lidl and Niederreiter, Lemma 2.3), and q - 1 divides lambda, so
//      a^(lambda + 1) = a in R and a^e = a^(((e - 1) mod lambda) + 1).
//      Plain e mod lambda would turn a^lambda into a^0 = 1, which is wrong
//      for a zero divisor.
//  (b) e < 0, only when the factor carries an inverse (the caller's
//      contract; build_prime_residues rejects the prime when a base with a
//      negative exponent is not invertible modulo it), so the base a is a
//      unit of R.  Its order divides lambda, so a^e = a^(e + k lambda) for
//      every integer k, and ((e - 1) mod lambda) + 1 is such an exponent.
//      The result is positive, so the entry then uses the base's table.
// Which of the equal forms is used changes only the cost, never the residue.
// Every entry of magnitude above lambda is reduced, positive by (a) or (b),
// so no exponent grows.  For a factor with an inverse, every entry of the
// row is also made positive by (b) when the window model of
// sliding_window_width charges less for one table of the reduced positive
// exponents than for a base table and an inverse table; squarings are left
// out of that comparison as they are of the width choice, and after the
// reduction no entry exceeds lambda in either form.
template <typename Factor>
void reduce_power_product_exponents(
        std::vector<const std::vector<flint::Fmpz>*>& rows,
        std::vector<std::vector<flint::Fmpz>>& storage,
        std::span<const Factor> factors, ulong prime,
        const flint::NmodPoly* modulus) noexcept {
    rows.resize(factors.size());
    for (std::size_t row = 0; row < factors.size(); ++row) {
        rows[row] = factors[row].exponents;
    }

    // lambda >= p - 1, so nothing reduces unless some magnitude exceeds
    // p - 1; only then is lambda computed.
    const auto reducible = [](const flint::Fmpz& exponent,
                              const Factor& factor) {
        return fmpz_sgn(exponent.raw()) > 0 ||
               (fmpz_sgn(exponent.raw()) < 0 && factor_has_inverse(factor));
    };
    flint::Fmpz floor;
    fmpz_set_ui(floor.raw(), prime - 1);
    bool large = false;
    for (std::size_t row = 0; row < factors.size() && !large; ++row) {
        for (const flint::Fmpz& exponent : *factors[row].exponents) {
            if (reducible(exponent, factors[row]) &&
                fmpz_cmpabs(exponent.raw(), floor.raw()) > 0) {
                large = true;
                break;
            }
        }
    }
    flint::Fmpz lambda;
    if (!large) {
        return;
    }
    if (modulus != nullptr) {
        if (!modular_unit_group_exponent(lambda, *modulus)) {
            return;
        }
    } else if (n_is_prime(prime) != 0) {
        fmpz_set_ui(lambda.raw(), prime - 1);
    } else {
        return;
    }

    storage.clear();
    storage.reserve(factors.size());
    std::vector<flint::Fmpz> reduced;
    for (std::size_t row = 0; row < factors.size(); ++row) {
        const Factor& factor = factors[row];
        const std::vector<flint::Fmpz>& exponents = *factor.exponents;
        // reduced[i] = ((e - 1) mod lambda) + 1 for every nonzero entry e
        // that may be reduced; zero and the other entries stay as given.
        reduced.resize(exponents.size());
        ExponentRowBits large_only;
        ExponentRowBits all_positive;
        bool large_entry = false;
        bool small_negative = false;
        for (std::size_t i = 0; i < exponents.size(); ++i) {
            const flint::Fmpz& exponent = exponents[i];
            if (fmpz_is_zero(exponent.raw()) != 0 ||
                !reducible(exponent, factor)) {
                fmpz_set(reduced[i].raw(), exponent.raw());
                large_only.add(exponent);
                all_positive.add(exponent);
                continue;
            }
            reduce_exponent(reduced[i], exponent, lambda);
            all_positive.add(reduced[i]);
            if (fmpz_cmpabs(exponent.raw(), lambda.raw()) > 0) {
                large_entry = true;
                large_only.add(reduced[i]);
            } else {
                small_negative =
                        small_negative || fmpz_sgn(exponent.raw()) < 0;
                large_only.add(exponent);
            }
        }
        const bool make_positive =
                small_negative &&
                all_positive.modeled_cost() < large_only.modeled_cost();
        if (!large_entry && !make_positive) {
            continue;
        }
        storage.emplace_back(exponents.size());
        std::vector<flint::Fmpz>& kept = storage.back();
        for (std::size_t i = 0; i < exponents.size(); ++i) {
            const bool use_reduced =
                    make_positive ||
                    fmpz_cmpabs(exponents[i].raw(), lambda.raw()) > 0;
            fmpz_set(kept[i].raw(),
                     use_reduced ? reduced[i].raw() : exponents[i].raw());
        }
        rows[row] = &kept;
    }
}

// The power products of modular_power_products in the ring of `ring`; the
// factor exponents are those of exponent_rows.
template <typename Ring, typename Factor>
void interleaved_power_products(
        flint::NmodMat& out, const Ring& ring,
        std::span<const Factor> factors,
        const std::vector<const std::vector<flint::Fmpz>*>& exponent_rows)
        noexcept {
    using Value = typename Ring::Value;
    const slong output_count = nmod_mat_nrows(out.raw());

    // One table per base and exponent sign that occurs, built once and
    // shared by every output.
    std::vector<OddPowerTable<Value>> positive_tables(factors.size());
    std::vector<OddPowerTable<Value>> negative_tables(factors.size());
    for (std::size_t row = 0; row < factors.size(); ++row) {
        flint_bitcnt_t positive_bits = 0;
        flint_bitcnt_t negative_bits = 0;
        flint_bitcnt_t positive_maximum = 0;
        flint_bitcnt_t negative_maximum = 0;
        for (const flint::Fmpz& exponent : *exponent_rows[row]) {
            const flint_bitcnt_t bits = fmpz_bits(exponent.raw());
            if (fmpz_sgn(exponent.raw()) > 0) {
                positive_bits += bits;
                positive_maximum = std::max(positive_maximum, bits);
            } else if (fmpz_sgn(exponent.raw()) < 0) {
                negative_bits += bits;
                negative_maximum = std::max(negative_maximum, bits);
            }
        }
        if (positive_bits != 0) {
            build_odd_power_table(
                    positive_tables[row], ring, factor_base(factors[row]),
                    sliding_window_width(positive_bits, positive_maximum));
        }
        if (negative_bits != 0) {
            build_odd_power_table(
                    negative_tables[row], ring, factor_inverse(factors[row]),
                    sliding_window_width(negative_bits, negative_maximum));
        }
    }

    std::vector<WindowStep<Value>> steps;
    Value value = ring.make();
    Value product = ring.make();
    flint::Fmpz magnitude;
    for (slong column = 0; column < output_count; ++column) {
        const std::size_t index = static_cast<std::size_t>(column);
        steps.clear();
        for (std::size_t row = 0; row < factors.size(); ++row) {
            const flint::Fmpz& exponent = (*exponent_rows[row])[index];
            if (fmpz_is_zero(exponent.raw()) != 0) {
                continue;
            }
            fmpz_abs(magnitude.raw(), exponent.raw());
            append_window_steps(steps, magnitude,
                                fmpz_sgn(exponent.raw()) < 0
                                        ? negative_tables[row]
                                        : positive_tables[row]);
        }
        std::sort(steps.begin(), steps.end(),
                  [](const WindowStep<Value>& left,
                     const WindowStep<Value>& right) {
                      return left.position > right.position;
                  });

        // Interleaved left-to-right evaluation: one squaring per bit position
        // below the highest window, shared by all bases of this output.
        ring.one(value);
        bool started = false;
        std::size_t next = 0;
        if (!steps.empty()) {
            for (slong position = static_cast<slong>(steps.front().position);
                 position >= 0; --position) {
                if (started) {
                    ring.mul(product, value, value);
                    value.swap(product);
                }
                while (next < steps.size() &&
                       steps[next].position ==
                               static_cast<flint_bitcnt_t>(position)) {
                    if (started) {
                        ring.mul(product, value, *steps[next].odd_power);
                        value.swap(product);
                    } else {
                        ring.set(value, *steps[next].odd_power);
                        started = true;
                    }
                    ++next;
                }
            }
        }
        ring.store(out, column, value);
    }
}

}  // namespace

bool modular_unit_group_exponent(flint::Fmpz& out,
                                 const flint::NmodPoly& modulus) noexcept {
    const ulong prime = modulus.raw()->mod.n;
    const slong degree = nmod_poly_degree(modulus.raw());
    if (degree < 1 || nmod_poly_get_coeff_ui(modulus.raw(), degree) != 1 ||
        n_is_prime(prime) == 0) {
        return false;
    }
    // A modulus with distinct linear factors only: R = F_p^degree, every
    // factor degree is 1, and the distinct-degree factorization below would
    // give the same lambda = p - 1.
    if (splits_into_distinct_linear_factors(modulus)) {
        fmpz_set_ui(out.raw(), prime - 1);
        return true;
    }
    if (nmod_poly_is_squarefree(modulus.raw()) == 0) {
        return false;
    }

    // Distinct-degree factorization: factor i is the product of the
    // irreducible factors of degree degrees[i], so lambda is the least common
    // multiple of p^degrees[i] - 1.
    flint::NmodPolyFactor factorization;
    std::vector<slong> degrees(static_cast<std::size_t>(degree) + 1);
    slong* degree_data = degrees.data();
    nmod_poly_factor_distinct_deg(factorization.raw(), modulus.raw(),
                                  &degree_data);
    flint::Fmpz term;
    slong total_degree = 0;
    fmpz_one(out.raw());
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        const slong factor_degree = degrees[static_cast<std::size_t>(i)];
        const slong product_degree =
                nmod_poly_degree(factorization.raw()->p + i);
        if (factor_degree <= 0 || product_degree <= 0 ||
            product_degree % factor_degree != 0) {
            return false;
        }
        total_degree += product_degree;
        fmpz_set_ui(term.raw(), prime);
        fmpz_pow_ui(term.raw(), term.raw(), static_cast<ulong>(factor_degree));
        fmpz_sub_ui(term.raw(), term.raw(), 1);
        fmpz_lcm(out.raw(), out.raw(), term.raw());
    }
    return total_degree == degree;
}

void modular_power_products(
        flint::NmodMat& power_residues,
        std::span<const ModularPowerProductFactor> factors,
        const flint::NmodPoly& modulus) noexcept {
    std::vector<const std::vector<flint::Fmpz>*> exponent_rows;
    std::vector<std::vector<flint::Fmpz>> reduced_exponents;
    reduce_power_product_exponents(exponent_rows, reduced_exponents, factors,
                                   modulus.raw()->mod.n, &modulus);
    interleaved_power_products(power_residues, PolynomialPowerRing(modulus),
                               factors, exponent_rows);
}

void modular_power_products_evaluated(
        flint::NmodMat& values,
        std::span<const EvaluatedPowerProductFactor> factors) noexcept {
    const nmod_t modulus = values.raw()->mod;
    std::vector<const std::vector<flint::Fmpz>*> exponent_rows;
    std::vector<std::vector<flint::Fmpz>> reduced_exponents;
    reduce_power_product_exponents(exponent_rows, reduced_exponents, factors,
                                   modulus.n, nullptr);
    interleaved_power_products(
            values,
            EvaluationPowerRing(modulus, static_cast<std::size_t>(
                                                 nmod_mat_ncols(values.raw()))),
            factors, exponent_rows);
}

const char* bounded_compact_reconstruction_status_name(
        BoundedCompactReconstructionStatus status) noexcept {
    switch (status) {
        case BoundedCompactReconstructionStatus::success:
            return "success";
        case BoundedCompactReconstructionStatus::invalid_input:
            return "invalid_input";
        case BoundedCompactReconstructionStatus::exponent_overflow:
            return "exponent_overflow";
        case BoundedCompactReconstructionStatus::coordinate_conversion_failed:
            return "coordinate_conversion_failed";
        case BoundedCompactReconstructionStatus::prime_exhausted:
            return "prime_exhausted";
        case BoundedCompactReconstructionStatus::certified_bound_violated:
            return "certified_bound_violated";
        case BoundedCompactReconstructionStatus::output_conversion_failed:
            return "output_conversion_failed";
    }
    return "unknown";
}

namespace {

// The defining polynomial reduced modulo prime keeps its degree, is monic and
// splits into distinct linear factors there.  When it is exactly Phi_n or
// psi_n (`recognized`, from cyclotomic_defining_polynomial) and prime does
// not divide n, the answer is the closed form of
// cyclotomic_defining_polynomial_splits, which gives the same answer as the
// test of the reduction below: Phi_n and psi_n are monic and integral, so
// they reduce to a monic polynomial of the same degree for every prime.
bool defining_polynomial_splits(
        const CyclotomicDefiningPolynomial& recognized,
        const flint::FmpqPoly& defining_polynomial, slong degree,
        ulong prime) noexcept {
    bool splits = false;
    if (cyclotomic_defining_polynomial_splits(splits, recognized, prime)) {
        return splits;
    }
    flint::NmodPoly modulus(prime);
    return reduce_power_polynomial(modulus, defining_polynomial, degree,
                                   prime) &&
           nmod_poly_degree(modulus.raw()) == degree &&
           splits_into_distinct_linear_factors(modulus);
}

// Modeled ring multiplications of the power-basis products at one prime: the
// window model of sliding_window_width summed over the tables
// (ExponentRowBits::modeled_cost) plus the squarings of the interleaved
// evaluation, which Moller ("Algorithms for multi-exponentiation",
// Section 3.1) puts at up to b - 1 for an output whose longest exponent has
// b bits.
double modeled_power_product_multiplications(
        const std::vector<PreparedFactor>& factors,
        std::size_t output_count) noexcept {
    std::vector<flint_bitcnt_t> output_bits(output_count, 0);
    double multiplications = 0.0;
    for (const PreparedFactor& factor : factors) {
        ExponentRowBits row;
        for (std::size_t column = 0;
             column < output_count && column < factor.exponents.size();
             ++column) {
            const flint::Fmpz& exponent = factor.exponents[column];
            row.add(exponent);
            output_bits[column] = std::max<flint_bitcnt_t>(
                    output_bits[column], fmpz_bits(exponent.raw()));
        }
        multiplications += row.modeled_cost();
    }
    for (const flint_bitcnt_t bits : output_bits) {
        if (bits > 1) {
            multiplications += static_cast<double>(bits - 1);
        }
    }
    return multiplications;
}

// Whether the CRT primes of a batch are searched for split primes, at which
// the products are computed in evaluation form (build_split_prime_residues):
// a ring multiplication there is degree products in F_p instead of one
// multiplication modulo the defining polynomial.
//
// Long exponents.  When some exponent has more than FLINT_BITS bits, it
// exceeds p - 1 for every word-size prime p, and at a split prime the unit
// group exponent of the residue ring is p - 1, so modular_power_products
// reduces every exponent below p and may drop the inverse tables.
//
// Short exponents.  The search is preferred when the modeled power-basis
// multiplications M at one prime (modeled_power_product_multiplications)
// reach degree * FLINT_BITS.  A split test (splits_into_distinct_linear_
// factors) is one exponentiation x^p modulo the defining polynomial;
// left-to-right binary exponentiation (Menezes, van Oorschot and Vanstone,
// Handbook of Applied Cryptography, Algorithm 14.79) takes at most
// FLINT_BITS - 1 squarings for p < 2^FLINT_BITS, the multiplications by x
// being shifts, so a test is charged FLINT_BITS multiplications.  If a
// candidate splits with probability delta, a scan of at most L candidates
// finds a split prime with probability P = 1 - (1 - delta)^L after P / delta
// tests on average, and a split prime saves about M, so the expected gain
// P (M - FLINT_BITS / delta) is not negative exactly when
// M >= FLINT_BITS / delta, whatever L is.  The threshold takes
// delta = 1 / degree, the Chebotarev density of the completely split primes
// of a Galois field; for any other field delta is smaller and the search
// can cost up to split_prime_scan_limit tests without finding a prime.  The
// model leaves out the root finding and change of basis at a split prime
// and the inverses modulo the defining polynomial at other primes.
bool prefers_split_primes(const std::vector<PreparedFactor>& factors,
                          std::size_t output_count, slong degree) noexcept {
    for (const PreparedFactor& factor : factors) {
        for (const flint::Fmpz& exponent : factor.exponents) {
            if (fmpz_bits(exponent.raw()) > FLINT_BITS) {
                return true;
            }
        }
    }
    return modeled_power_product_multiplications(factors, output_count) >=
           static_cast<double>(degree) * static_cast<double>(FLINT_BITS);
}

bool reconstruct_prepared(BoundedCompactReconstructionReport& report,
                          std::vector<Element>& out, const Order& order,
                          std::vector<PreparedFactor>& factors,
                          std::size_t output_count,
                          flint::FmpzConstRef certified_coordinate_bound,
                          const BoundedCompactReconstructionOptions& options)
        noexcept {
    const NumberField* field = order.parent();
    if (!prepare_power_polynomials(factors)) {
        report.status = BoundedCompactReconstructionStatus::
                coordinate_conversion_failed;
        return false;
    }

    const slong degree = order.degree();
    const nf_struct* raw_field = field->raw_flint_field();
    flint::FmpqPoly defining_polynomial;
    flint::FmpqMat order_basis(degree, degree);
    if (raw_field == nullptr ||
        !order.get_basis(flint::FmpqMatRef(order_basis))) {
        report.status = BoundedCompactReconstructionStatus::
                coordinate_conversion_failed;
        return false;
    }
    fmpq_poly_set(defining_polynomial.raw(), raw_field->pol);

    flint::Fmpz twice_bound;
    fmpz_mul_2exp(twice_bound.raw(), certified_coordinate_bound.raw(), 1);
    flint::Fmpz modulus;
    fmpz_one(modulus.raw());
    flint::FmpzMat reconstructed(static_cast<slong>(output_count), degree);

    // Split-prime preference (prefers_split_primes).  Every prime gives the
    // same reconstruction, so the choice changes only the cost.  Up to
    // split_prime_scan_limit candidates are tested for each prime the CRT
    // needs: by the Chebotarev density theorem the completely split primes
    // have density 1 / [L : Q] for the Galois closure L of the field, which
    // is 1 / degree for a Galois field, so the limit of 4 * degree
    // candidates finds one with probability
    // 1 - (1 - 1 / degree)^(4 degree) > 1 - e^-4 there.  When the limit is
    // reached, that candidate is tried as before and the preference ends.
    // A prime found to split here is also passed to build_prime_residues as
    // split, which then computes the products in evaluation form; every other
    // prime, including every prime of a batch without the preference, uses
    // the power basis.
    // The split test has a closed form when the defining polynomial is
    // exactly Phi_n or psi_n; it is recognized once, only when the
    // preference applies.
    bool prefer_split_primes =
            prefers_split_primes(factors, output_count, degree);
    const CyclotomicDefiningPolynomial recognized =
            prefer_split_primes ? cyclotomic_defining_polynomial(*field)
                                : CyclotomicDefiningPolynomial{};
    const std::size_t split_prime_scan_limit =
            4 * static_cast<std::size_t>(degree);
    std::size_t split_prime_scan = 0;

    ulong prime_cursor = options.prime_search_start;
    while (fmpz_cmp(modulus.raw(), twice_bound.raw()) <= 0 ||
           report.primes_used == 0) {
        if (prime_cursor >= UWORD_MAX_PRIME) {
            report.status = BoundedCompactReconstructionStatus::prime_exhausted;
            return false;
        }
        const ulong prime = n_nextprime(prime_cursor, 1);
        prime_cursor = prime;
        bool split = false;
        if (prefer_split_primes) {
            if (defining_polynomial_splits(recognized, defining_polynomial,
                                           degree, prime)) {
                split = true;
                split_prime_scan = 0;
            } else if (++split_prime_scan < split_prime_scan_limit) {
                ++report.nonsplit_prime_skips;
                continue;
            } else {
                prefer_split_primes = false;
            }
        }
        report.last_prime = prime;
        ++report.prime_trials;

        flint::NmodMat residues(static_cast<slong>(output_count), degree,
                                prime);
        if (!build_prime_residues(report, residues, factors,
                                  defining_polynomial, order_basis, degree,
                                  prime, split)) {
            continue;
        }

        if (report.primes_used == 0) {
            fmpz_mat_set_nmod_mat(reconstructed.raw(), residues.raw());
        } else {
            fmpz_mat_CRT_ui(reconstructed.raw(), reconstructed.raw(),
                            modulus.raw(), residues.raw(), 1);
        }
        fmpz_mul_ui(modulus.raw(), modulus.raw(), prime);
        ++report.primes_used;
        fmpz_set(report.centered_crt_modulus.raw(), modulus.raw());
        report.modulus_bits = fmpz_bits(modulus.raw());
    }

    if (!reconstructed_coordinates_within_bound(reconstructed,
                                                certified_coordinate_bound)) {
        report.status =
                BoundedCompactReconstructionStatus::certified_bound_violated;
        return false;
    }
    if (!publish_reconstructed_elements(out, reconstructed, order)) {
        report.status =
                BoundedCompactReconstructionStatus::output_conversion_failed;
        return false;
    }

    report.status = BoundedCompactReconstructionStatus::success;
    return true;
}

}  // namespace

bool bounded_compact_reconstruct(
        BoundedCompactReconstructionReport& report, std::vector<Element>& out,
        const Order& order, std::span<const FactoredElement> inputs,
        flint::FmpzConstRef certified_coordinate_bound,
        const BoundedCompactReconstructionOptions& options) noexcept {
    reset_report(report);
    const NumberField* field = order.parent();
    if (!order.is_defined() || field == nullptr || !order.has_basis() ||
        order.degree() <= 0 || fmpz_sgn(certified_coordinate_bound.raw()) < 0 ||
        inputs.size() >
                static_cast<std::size_t>(std::numeric_limits<slong>::max()) ||
        options.prime_search_start >= UWORD_MAX_PRIME) {
        return false;
    }

    if (inputs.empty()) {
        std::vector<Element> empty;
        out.swap(empty);
        report.status = BoundedCompactReconstructionStatus::success;
        return true;
    }

    std::vector<PreparedFactor> factors;
    if (!prepare_exponent_matrix(report, factors, order, inputs)) {
        if (report.status !=
            BoundedCompactReconstructionStatus::exponent_overflow) {
            report.status = BoundedCompactReconstructionStatus::invalid_input;
        }
        return false;
    }
    return reconstruct_prepared(report, out, order, factors, inputs.size(),
                                certified_coordinate_bound, options);
}

bool bounded_compact_reconstruct(
        BoundedCompactReconstructionReport& report, std::vector<Element>& out,
        const Order& order, const WideFactoredUnits& inputs,
        flint::FmpzConstRef certified_coordinate_bound,
        const BoundedCompactReconstructionOptions& options) noexcept {
    reset_report(report);
    const NumberField* field = order.parent();
    const slong outputs = flint::fmpz_mat_nrows(inputs.exponents);
    if (!order.is_defined() || field == nullptr || !order.has_basis() ||
        order.degree() <= 0 || fmpz_sgn(certified_coordinate_bound.raw()) < 0 ||
        outputs < 0 || options.prime_search_start >= UWORD_MAX_PRIME) {
        return false;
    }

    if (outputs == 0) {
        std::vector<Element> empty;
        out.swap(empty);
        report.status = BoundedCompactReconstructionStatus::success;
        return true;
    }

    std::vector<PreparedFactor> factors;
    if (!prepare_wide_exponent_matrix(report, factors, order, inputs)) {
        report.status = BoundedCompactReconstructionStatus::invalid_input;
        return false;
    }
    return reconstruct_prepared(report, out, order, factors,
                                static_cast<std::size_t>(outputs),
                                certified_coordinate_bound, options);
}
}  // namespace silex::detail
