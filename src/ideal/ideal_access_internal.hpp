#pragma once

#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>

namespace silex::detail {

// Noninstalled access to Ideal internals.  set_hnf_direct is for tests only;
// set_known_ideal_hnf is for library code.
class IdealAccess {
public:
    enum class HnfInput { unnormalized, normalized };
    enum class IdealTest { always, assert_only };

    // The one routine behind Ideal::set_hnf, set_known_ideal_hnf and the
    // private Ideal::set_known_ideal_rows.  `unnormalized` input is put in
    // normalized HNF (failing if singular or misshapen); `normalized` input
    // is only shape-checked and copied.  IdealTest::always rejects a lattice
    // that is not an ideal; assert_only only asserts it (NDEBUG removes it).
    // On failure `ideal` is unchanged.
    static bool store_hnf(Ideal& ideal, flint::FmpzMatConstRef hnf,
                          HnfInput input, IdealTest test) noexcept;

    // Stores the normalized HNF of the full-rank lattice spanned by the rows
    // of the n x n matrix `hnf`, as Ideal::set_hnf does, but without testing
    // that the lattice is an ideal of the parent order; builds with assertions
    // enabled (NDEBUG undefined) still assert that test.  Fails, leaving
    // `ideal` unchanged, where set_hnf fails for a reason other than that
    // test: an undefined `ideal`, a misshapen or singular matrix.  Only for a
    // lattice that is an ideal by construction, with the argument stated at
    // the call.
    static bool set_known_ideal_hnf(Ideal& ideal,
                                    flint::FmpzMatConstRef hnf) noexcept;


    // The one routine behind Ideal::set_principal and
    // detail::set_principal_ideal_known_norm.  `norm` is null, or the norm of
    // `generator` as the trust helper below states; see the definition in
    // ideal.cpp.
    static bool set_principal(Ideal& ideal, const OrderElement& generator,
                              const fmpq* norm) noexcept;

    // Implementations of detail::set_known_two_generator_ideal and
    // detail::set_integral_ideal_known_hnf, defined in ideal.cpp.
    static bool set_known_two_generator(
            Ideal& out, flint::FmpzConstRef scalar_generator,
            const OrderElement& element_generator) noexcept;
    static bool set_integral_known_hnf(FractionalIdeal& out,
                                       const Ideal& numerator) noexcept;
    // Implementation of detail::set_known_linear_two_generator_ideal, defined
    // in ideal.cpp.
    static bool set_known_linear_two_generator(
            Ideal& out, flint::FmpzConstRef scalar_generator,
            flint::FmpzConstRef root) noexcept;
    // Implementation of detail::set_known_factor_two_generator_ideal,
    // defined in ideal.cpp.
    static bool set_known_factor_two_generator(
            Ideal& out, flint::FmpzConstRef scalar_generator,
            flint::FmpzPolyConstRef factor) noexcept;

    // Implementations of detail::multiply_integral_ideal_by_two_generator
    // and detail::multiply_integral_ideal_by_element, defined in ideal.cpp.
    static bool multiply_two_generator(
            Ideal& out, const Ideal& ideal,
            flint::FmpzConstRef scalar_generator,
            const OrderElement& element_generator,
            const DiagnosticsContext* diagnostics) noexcept;
    static bool multiply_by_element(
            Ideal& out, const Ideal& ideal, const Element& multiplier,
            const DiagnosticsContext* diagnostics) noexcept;

    // Stores `hnf` as the ideal's basis matrix without normalizing it or
    // checking that its row span is an ideal; only the n x n shape is
    // checked.  Every public construction route stores a full-rank upper
    // triangular HNF.  This hook exists so tests can build the other stored
    // shapes that reach the defensive lattice fallbacks of Ideal::contains.
    // Library code must not use it.
    static bool set_hnf_direct(Ideal& ideal,
                               flint::FmpzMatConstRef hnf) noexcept {
        return ideal.set_hnf_direct(hnf);
    }
};

// Trust helpers, noninstalled.  set_known_two_generator_ideal stores the HNF
// of a O + beta O for a positive scalar a and an element beta of the same
// order, trusting that the lattice is an ideal.  set_integral_ideal_known_hnf
// sets `out` to the integral fractional ideal with the given numerator.
// Both return false, leaving `out` unchanged, on undefined or mismatched
// input.
inline bool set_known_two_generator_ideal(
        Ideal& out, flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator) noexcept {
    return IdealAccess::set_known_two_generator(out, scalar_generator,
                                                element_generator);
}

// Trust helper, noninstalled.  The closed form of
// set_known_two_generator_ideal(out, a, theta - r) on an equation order
// Z[theta] whose defining polynomial F is monic integral: stores the same HNF
// of a O + (theta - r) O without a modular HNF.  It applies only when
// F(r) = 0 mod a and r is zero or a unit modulo a; otherwise, or on undefined
// input, a non-equation order or a <= 0, it returns false and leaves `out`
// unchanged, and the caller takes the general routine.
inline bool set_known_linear_two_generator_ideal(
        Ideal& out, flint::FmpzConstRef scalar_generator,
        flint::FmpzConstRef root) noexcept {
    return IdealAccess::set_known_linear_two_generator(out, scalar_generator,
                                                       root);
}

// Trust helper, noninstalled.  The closed form of
// set_known_two_generator_ideal(out, a, g(theta)) on an equation order
// Z[theta] whose defining polynomial F is monic integral, for a monic integer
// polynomial g of degree 1 <= f <= n: stores the same HNF of
// a O + g(theta) O without a modular HNF.  It applies only when g divides F
// modulo a and g(0) is a unit modulo a; otherwise, or on undefined input, a
// non-equation order, a non-monic g or a <= 0, it returns false and leaves
// `out` unchanged, and the caller takes the general routine.
inline bool set_known_factor_two_generator_ideal(
        Ideal& out, flint::FmpzConstRef scalar_generator,
        flint::FmpzPolyConstRef factor) noexcept {
    return IdealAccess::set_known_factor_two_generator(out, scalar_generator,
                                                       factor);
}

// Trust helper, noninstalled.  Sets `out` to the principal ideal
// `generator` O exactly as Ideal::set_principal does (same stored HNF, same
// failures, `out` unchanged on failure).  The caller vouches that `norm` is
// N_{K/Q}(generator); when it is a nonzero integer, its absolute value is the
// modulus of a determinant-modular HNF instead of a full HNF.  Any other
// `norm` (noninteger or zero) takes the Ideal::set_principal route.
inline bool set_principal_ideal_known_norm(Ideal& out,
                                           const OrderElement& generator,
                                           flint::FmpqConstRef norm) noexcept {
    return IdealAccess::set_principal(out, generator, norm.raw());
}

inline bool set_integral_ideal_known_hnf(FractionalIdeal& out,
                                         const Ideal& numerator) noexcept {
    return IdealAccess::set_integral_known_hnf(out, numerator);
}

// Noninstalled integral-ideal products.  The first sets `out` to (a O + beta O)
// times `ideal`; the second to `multiplier` times `ideal`, requiring an
// integral result.  See the comment above the definitions in ideal.cpp.
inline bool multiply_integral_ideal_by_two_generator(
        Ideal& out,
        const Ideal& ideal,
        flint::FmpzConstRef scalar_generator,
        const OrderElement& element_generator,
        const DiagnosticsContext* diagnostics = nullptr) noexcept {
    return IdealAccess::multiply_two_generator(
            out, ideal, scalar_generator, element_generator, diagnostics);
}

inline bool multiply_integral_ideal_by_element(
        Ideal& out,
        const Ideal& ideal,
        const Element& multiplier,
        const DiagnosticsContext* diagnostics = nullptr) noexcept {
    return IdealAccess::multiply_by_element(out, ideal, multiplier,
                                            diagnostics);
}

}  // namespace silex::detail
