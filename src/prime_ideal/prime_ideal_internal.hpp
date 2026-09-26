#pragma once

#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/prime_ideal.hpp>

namespace silex::detail {

enum class RetainedQuadraticPrimeKind {
    inert,
    ramified,
    split,
};

// True when prime ideals of `order` store their residue polynomial in the
// integral generator omega of a maximal quadratic-backend order (basis
// [1, omega]).  Otherwise the residue polynomial is in the NumberField
// generator alpha.  See PrimeIdeal::residue_polynomial.
bool residue_polynomial_uses_quadratic_integral_generator(
        const Order& order) noexcept;

// Writes `element` as numerator(v) / denominator, where v is the variable of
// `prime`'s stored residue polynomial (alpha, or omega per the predicate
// above) and numerator has integer coefficients.  On the alpha convention the
// element's alpha-polynomial is used directly: this is the image of its order
// coordinates under the order basis matrix, not the coordinates themselves.
// Fails, leaving outputs unspecified, when p divides the denominator or, on
// the alpha convention, when p divides the order basis denominator (for an
// order containing Z[alpha], exactly when p divides [O : Z[alpha]]).  Reducing
// numerator modulo (p, residue polynomial) and multiplying by
// denominator^{-1} mod p gives the residue class of `element` modulo P.
bool residue_variable_numerator(flint::FmpzPoly& numerator,
                                flint::Fmpz& denominator,
                                const PrimeIdeal& prime,
                                const Element& element) noexcept;

class MaximalQuadraticPrimeAccess {
public:
    static bool set_from_integral_generator_factor(
            PrimeIdeal& out,
            const Order& order,
            flint::FmpzConstRef p,
            flint::FmpzModPolyConstRef factor,
            slong ramification_index,
            const Element& integral_generator,
            const flint::FmpzModCtx& context) noexcept;

    static bool set_first_degree_one_prime(
            PrimeIdeal& out,
            RetainedQuadraticPrimeKind& kind,
            const Order& order,
            flint::FmpzConstRef p,
            const DiagnosticsContext* diagnostics) noexcept;
};

}  // namespace silex::detail
