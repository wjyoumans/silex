#include "compact_reconstruction_internal.hpp"
#include "wide_factored_internal.hpp"

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
// inverse) modulo the modulus: the sliding-window table of Handbook of
// Applied Cryptography, Algorithm 14.85 step 1.
struct OddPowerTable {
    std::vector<flint::NmodPoly> odd_powers;
    ulong width = 0;
};

// One window of a sliding-window recoding: multiply by odd_power at bit
// position `position` of the left-to-right evaluation.
struct WindowStep {
    flint_bitcnt_t position = 0;
    const flint::NmodPoly* odd_power = nullptr;
};

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

void build_odd_power_table(OddPowerTable& table, const flint::NmodPoly& base,
                           ulong width, const flint::NmodPoly& modulus) noexcept {
    const ulong prime = modulus.raw()->mod.n;
    const std::size_t size = static_cast<std::size_t>(UWORD(1) << (width - 1));
    table.width = width;
    table.odd_powers.clear();
    table.odd_powers.reserve(size);
    table.odd_powers.emplace_back(prime);
    nmod_poly_set(table.odd_powers.back().raw(), base.raw());
    if (size == 1) {
        return;
    }
    flint::NmodPoly square(prime);
    nmod_poly_mulmod(square.raw(), base.raw(), base.raw(), modulus.raw());
    for (std::size_t i = 1; i < size; ++i) {
        table.odd_powers.emplace_back(prime);
        nmod_poly_mulmod(table.odd_powers.back().raw(),
                         table.odd_powers[i - 1].raw(), square.raw(),
                         modulus.raw());
    }
}

// Left-to-right sliding-window recoding of a positive magnitude
// (Handbook of Applied Cryptography, Algorithm 14.85 step 3): each maximal
// window of at most table.width bits that starts and ends with a one becomes
// one multiplication by the matching odd power at the window's lowest bit.
void append_window_steps(std::vector<WindowStep>& steps,
                         const flint::Fmpz& magnitude,
                         const OddPowerTable& table) noexcept {
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
        steps.push_back(WindowStep{static_cast<flint_bitcnt_t>(low),
                                   &table.odd_powers[(value - 1) / 2]});
        high = low - 1;
    }
}

bool build_prime_residues(BoundedCompactReconstructionReport& report,
                          flint::NmodMat& residues,
                          const std::vector<PreparedFactor>& factors,
                          const flint::FmpqPoly& defining_polynomial,
                          const flint::FmpqMat& order_basis, slong degree,
                          ulong prime) noexcept {
    flint::NmodPoly modulus(prime);
    if (!reduce_power_polynomial(modulus, defining_polynomial, degree, prime) ||
        nmod_poly_degree(modulus.raw()) != degree ||
        nmod_poly_get_coeff_ui(modulus.raw(), degree) != 1) {
        ++report.field_polynomial_prime_rejections;
        return false;
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
// the unit group of R = (Z/pZ)[x]/(modulus).  A row that changes is copied
// into storage and reduced there; every other row stays the caller's.
//
// Soundness.  modular_unit_group_exponent returns lambda only when p is prime
// and the modulus is monic, nonconstant and squarefree modulo p (a product of
// distinct linear factors, its lambda = p - 1 case, is squarefree).  Then R is
// the product of the fields F_q = F_p[x]/(f_i) over the irreducible factors
// f_i of the modulus (Chinese remainder theorem; Lidl and Niederreiter,
// Finite Fields, Theorem 1.61), and lambda = lcm(q - 1).
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
void reduce_power_product_exponents(
        std::vector<const std::vector<flint::Fmpz>*>& rows,
        std::vector<std::vector<flint::Fmpz>>& storage,
        std::span<const ModularPowerProductFactor> factors,
        const flint::NmodPoly& modulus) noexcept {
    rows.resize(factors.size());
    for (std::size_t row = 0; row < factors.size(); ++row) {
        rows[row] = factors[row].exponents;
    }

    // lambda >= p - 1, so nothing reduces unless some magnitude exceeds
    // p - 1; only then is lambda computed.
    const ulong prime = modulus.raw()->mod.n;
    const auto reducible = [](const flint::Fmpz& exponent,
                              const ModularPowerProductFactor& factor) {
        return fmpz_sgn(exponent.raw()) > 0 ||
               (fmpz_sgn(exponent.raw()) < 0 && factor.inverse != nullptr);
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
    if (!large || !modular_unit_group_exponent(lambda, modulus)) {
        return;
    }

    storage.clear();
    storage.reserve(factors.size());
    std::vector<flint::Fmpz> reduced;
    for (std::size_t row = 0; row < factors.size(); ++row) {
        const ModularPowerProductFactor& factor = factors[row];
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
    const ulong prime = modulus.raw()->mod.n;
    const slong output_count = nmod_mat_nrows(power_residues.raw());
    const slong degree = nmod_mat_ncols(power_residues.raw());

    std::vector<const std::vector<flint::Fmpz>*> exponent_rows;
    std::vector<std::vector<flint::Fmpz>> reduced_exponents;
    reduce_power_product_exponents(exponent_rows, reduced_exponents, factors,
                                   modulus);

    // One table per base and exponent sign that occurs, built once and
    // shared by every output.
    std::vector<OddPowerTable> positive_tables(factors.size());
    std::vector<OddPowerTable> negative_tables(factors.size());
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
                    positive_tables[row], *factors[row].base,
                    sliding_window_width(positive_bits, positive_maximum),
                    modulus);
        }
        if (negative_bits != 0) {
            build_odd_power_table(
                    negative_tables[row], *factors[row].inverse,
                    sliding_window_width(negative_bits, negative_maximum),
                    modulus);
        }
    }

    std::vector<WindowStep> steps;
    flint::NmodPoly value(prime);
    flint::NmodPoly product(prime);
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
                  [](const WindowStep& left, const WindowStep& right) {
                      return left.position > right.position;
                  });

        // Interleaved left-to-right evaluation: one squaring per bit position
        // below the highest window, shared by all bases of this output.
        nmod_poly_one(value.raw());
        bool started = false;
        std::size_t next = 0;
        if (!steps.empty()) {
            for (slong position = static_cast<slong>(steps.front().position);
                 position >= 0; --position) {
                if (started) {
                    nmod_poly_mulmod(product.raw(), value.raw(), value.raw(),
                                     modulus.raw());
                    value.swap(product);
                }
                while (next < steps.size() &&
                       steps[next].position ==
                               static_cast<flint_bitcnt_t>(position)) {
                    if (started) {
                        nmod_poly_mulmod(product.raw(), value.raw(),
                                         steps[next].odd_power->raw(),
                                         modulus.raw());
                        value.swap(product);
                    } else {
                        nmod_poly_set(value.raw(),
                                      steps[next].odd_power->raw());
                        started = true;
                    }
                    ++next;
                }
            }
        }
        for (slong i = 0; i < degree; ++i) {
            nmod_mat_entry(power_residues.raw(), column, i) =
                    nmod_poly_get_coeff_ui(value.raw(), i);
        }
    }
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
// splits into distinct linear factors there.
bool defining_polynomial_splits(const flint::FmpqPoly& defining_polynomial,
                                slong degree, ulong prime) noexcept {
    flint::NmodPoly modulus(prime);
    return reduce_power_polynomial(modulus, defining_polynomial, degree,
                                   prime) &&
           nmod_poly_degree(modulus.raw()) == degree &&
           splits_into_distinct_linear_factors(modulus);
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

    // Split-prime preference.  When some exponent has more than FLINT_BITS
    // bits, it exceeds p - 1 for every word-size prime p, and at a prime
    // where the reduced defining polynomial splits into distinct linear
    // factors the unit group exponent of the residue ring is p - 1, so
    // modular_power_products reduces every exponent below p and may drop the
    // inverse tables.  Every prime gives the same reconstruction, so the
    // choice changes only the cost.  Up to split_prime_scan_limit candidates
    // are tested for each prime the CRT needs: by the Chebotarev density
    // theorem the completely split primes have density 1 / [L : Q] for the
    // Galois closure L of the field, which is 1 / degree for a Galois field,
    // so the limit of 4 * degree candidates finds one with probability
    // 1 - (1 - 1 / degree)^(4 degree) > 1 - e^-4 there.  When the limit is
    // reached, that candidate is tried as before and the preference ends.
    bool prefer_split_primes = false;
    for (const PreparedFactor& factor : factors) {
        for (const flint::Fmpz& exponent : factor.exponents) {
            if (fmpz_bits(exponent.raw()) > FLINT_BITS) {
                prefer_split_primes = true;
            }
        }
    }
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
        if (prefer_split_primes) {
            if (defining_polynomial_splits(defining_polynomial, degree,
                                           prime)) {
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
                                  prime)) {
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
