#include "compact_reconstruction_internal.hpp"

#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/nmod_mat.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/order_element.hpp>

#include <flint/fmpz_mat.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
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
        : representative(&representative_value), exponents(output_count, 0) {}

    const Element* representative = nullptr;
    std::vector<slong> exponents;
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

bool exponent_row_is_zero(const PreparedFactor& factor) noexcept {
    return std::all_of(factor.exponents.begin(), factor.exponents.end(),
                       [](slong exponent) { return exponent == 0; });
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
            if (!checked_add_slong(match->exponents[column], entry.exponent)) {
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

bool prepare_power_polynomials(std::vector<PreparedFactor>& factors) noexcept {
    for (PreparedFactor& factor : factors) {
        if (factor.representative == nullptr ||
            !factor.representative->get_fmpq_poly(
                    flint::FmpqPolyRef(factor.power_polynomial))) {
            return false;
        }
        factor.needs_inverse =
                std::any_of(factor.exponents.begin(), factor.exponents.end(),
                            [](slong exponent) { return exponent < 0; });
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

ulong exponent_magnitude(slong exponent) noexcept {
    if (exponent >= 0) {
        return static_cast<ulong>(exponent);
    }
    ulong magnitude = static_cast<ulong>(-(exponent + 1));
    return magnitude + 1;
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
    flint::NmodMat power_residues(static_cast<slong>(output_count), degree,
                                  prime);
    for (std::size_t column = 0; column < output_count; ++column) {
        flint::NmodPoly value(prime);
        flint::NmodPoly power(prime);
        flint::NmodPoly product(prime);
        nmod_poly_one(value.raw());
        for (std::size_t row = 0; row < factors.size(); ++row) {
            const slong exponent = factors[row].exponents[column];
            if (exponent == 0) {
                continue;
            }
            const flint::NmodPoly& base =
                    exponent < 0 ? inverses[row] : bases[row];
            nmod_poly_powmod_ui_binexp(power.raw(), base.raw(),
                                       exponent_magnitude(exponent),
                                       modulus.raw());
            nmod_poly_mulmod(product.raw(), value.raw(), power.raw(),
                             modulus.raw());
            value.swap(product);
        }
        for (slong i = 0; i < degree; ++i) {
            nmod_mat_entry(power_residues.raw(), static_cast<slong>(column),
                           i) = nmod_poly_get_coeff_ui(value.raw(), i);
        }
    }

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

}  // namespace

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
    flint::FmpzMat reconstructed(static_cast<slong>(inputs.size()), degree);

    ulong prime_cursor = options.prime_search_start;
    while (fmpz_cmp(modulus.raw(), twice_bound.raw()) <= 0 ||
           report.primes_used == 0) {
        if (prime_cursor >= UWORD_MAX_PRIME) {
            report.status = BoundedCompactReconstructionStatus::prime_exhausted;
            return false;
        }
        const ulong prime = n_nextprime(prime_cursor, 1);
        prime_cursor = prime;
        report.last_prime = prime;
        ++report.prime_trials;

        flint::NmodMat residues(static_cast<slong>(inputs.size()), degree,
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

}  // namespace silex::detail
