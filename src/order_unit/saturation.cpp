#include "order_unit_internal.hpp"
#include "relation_unit_internal.hpp"

#include "../class_group/class_group_internal.hpp"
#include "../element/element_internal.hpp"
#include "../prime_ideal/prime_ideal_internal.hpp"
#include "../residue_field/residue_field_internal.hpp"
#include "../residue_field/word_frobenius_internal.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/unit.hpp>

namespace silex {
namespace {

bool embedding_has_parent(EmbeddingContext& embeddings,
                          const NumberField* field) noexcept {
    return field != nullptr && embeddings.parent() != nullptr &&
           embeddings.parent()->has_same_data(*field);
}

bool append_column(flint::FmpzMat& matrix,
                   const flint::FmpzMat& column) noexcept {
    const slong rows = flint::fmpz_mat_nrows(matrix);
    const slong old_cols = flint::fmpz_mat_ncols(matrix);
    if (flint::fmpz_mat_nrows(column) != rows ||
        flint::fmpz_mat_ncols(column) != 1) {
        return false;
    }

    flint::FmpzMat candidate(rows, old_cols + 1);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < old_cols; ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(candidate, i, j),
                            flint::FmpzConstRef(
                                    flint::fmpz_mat_entry(matrix, i, j).raw()));
        }
        flint::fmpz_set(flint::fmpz_mat_entry(candidate, i, old_cols),
                        flint::fmpz_mat_entry(column, i, 0));
    }

    matrix = std::move(candidate);
    return true;
}

void set_rank_one_candidate_kernel(flint::FmpzMat& kernel) noexcept {
    flint::FmpzMat candidate(1, 1);
    flint::fmpz_one(flint::fmpz_mat_entry(candidate, 0, 0));
    kernel = std::move(candidate);
}

bool field_polynomial_is_monic_integral(
        const fmpq_poly_struct* polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree < 1) {
        return false;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (!flint::fmpz_is_one(flint::fmpq_den_ref(coeff))) {
            return false;
        }
    }

    fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, degree);
    return flint::fmpz_is_one(flint::fmpq_num_ref(coeff));
}

class NmodPoly {
public:
    explicit NmodPoly(ulong modulus) noexcept {
        nmod_poly_init(value_, modulus);
    }

    ~NmodPoly() noexcept {
        nmod_poly_clear(value_);
    }

    NmodPoly(const NmodPoly&) = delete;
    NmodPoly& operator=(const NmodPoly&) = delete;

    nmod_poly_t& raw() noexcept {
        return value_;
    }

    const nmod_poly_t& raw() const noexcept {
        return value_;
    }

private:
    nmod_poly_t value_;
};

class NmodPolyFactor {
public:
    NmodPolyFactor() noexcept {
        nmod_poly_factor_init(value_);
    }

    ~NmodPolyFactor() noexcept {
        nmod_poly_factor_clear(value_);
    }

    NmodPolyFactor(const NmodPolyFactor&) = delete;
    NmodPolyFactor& operator=(const NmodPolyFactor&) = delete;

    nmod_poly_factor_t& raw() noexcept {
        return value_;
    }

    const nmod_poly_factor_t& raw() const noexcept {
        return value_;
    }

private:
    nmod_poly_factor_t value_;
};

void fmpq_poly_get_integral_fmpz_mod_poly(
        flint::FmpzModPoly& out,
        const fmpq_poly_struct* polynomial,
        const flint::FmpzModCtx& ctx) noexcept {
    flint::Fmpz coeff;
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpz(coeff.raw(), polynomial, i);
        fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, coeff.raw(), ctx.raw());
    }
}

void fmpq_poly_get_integral_nmod_poly(NmodPoly& out,
                                      const fmpq_poly_struct* polynomial,
                                      ulong modulus) noexcept {
    flint::Fmpz coeff;
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpz(coeff.raw(), polynomial, i);
        nmod_poly_set_coeff_ui(out.raw(), i,
                               fmpz_fdiv_ui(coeff.raw(), modulus));
    }
}

bool degree_one_root_from_mod_poly(flint::Fmpz& root,
                                   const fmpz_mod_poly_struct* polynomial,
                                   flint::FmpzConstRef p,
                                   const flint::FmpzModCtx& ctx) noexcept {
    if (fmpz_mod_poly_degree(polynomial, ctx.raw()) != 1) {
        return false;
    }

    flint::Fmpz constant;
    flint::Fmpz leading;
    flint::Fmpz inverse_leading;
    fmpz_mod_poly_get_coeff_fmpz(constant.raw(), polynomial, 0, ctx.raw());
    fmpz_mod_poly_get_coeff_fmpz(leading.raw(), polynomial, 1, ctx.raw());
    fmpz_mod(leading.raw(), leading.raw(), p.raw());
    if (fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) == 0) {
        return false;
    }

    fmpz_neg(root.raw(), constant.raw());
    fmpz_mod(root.raw(), root.raw(), p.raw());
    fmpz_mul(root.raw(), root.raw(), inverse_leading.raw());
    fmpz_mod(root.raw(), root.raw(), p.raw());
    return true;
}

bool degree_one_root_from_nmod_poly(flint::Fmpz& root,
                                    const nmod_poly_struct* polynomial)
        noexcept {
    if (nmod_poly_degree(polynomial) != 1) {
        return false;
    }

    const ulong constant = nmod_poly_get_coeff_ui(polynomial, 0);
    const ulong leading = nmod_poly_get_coeff_ui(polynomial, 1);
    if (leading == 0) {
        return false;
    }

    ulong root_ui = nmod_neg(constant, polynomial->mod);
    if (leading != 1) {
        root_ui = nmod_mul(root_ui, nmod_inv(leading, polynomial->mod),
                           polynomial->mod);
    }
    flint::fmpz_set_ui(flint::FmpzRef(root), root_ui);
    return true;
}

// Appends the roots of the linear factors in `factorization`.
void push_nmod_linear_factor_roots(std::vector<flint::Fmpz>& roots,
                                   const NmodPolyFactor& factorization)
        noexcept {
    roots.reserve(roots.size() +
                  static_cast<std::size_t>(factorization.raw()->num));
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        flint::Fmpz root;
        if (degree_one_root_from_nmod_poly(root, factorization.raw()->p + i)) {
            roots.push_back(std::move(root));
        }
    }
}

// Whether every value in `roots` is a root of `polynomial`.
bool nmod_poly_vanishes_at_roots(const NmodPoly& polynomial,
                                 const std::vector<flint::Fmpz>& roots)
        noexcept {
    for (const flint::Fmpz& root : roots) {
        if (nmod_poly_evaluate_nmod(
                    polynomial.raw(),
                    flint::fmpz_get_ui(flint::FmpzConstRef(root))) != 0) {
            return false;
        }
    }
    return true;
}

// Whether g = x^n + c, n >= 1: monic with no middle terms.  The root
// finder checks this once and tries the closed forms only for binomials.
bool nmod_poly_is_monic_binomial(const NmodPoly& g) noexcept {
    const slong degree = nmod_poly_degree(g.raw());
    const ulong* coeffs = g.raw()->coeffs;
    if (degree < 1 || coeffs[degree] != 1) {
        return false;
    }
    for (slong i = 1; i < degree; ++i) {
        if (coeffs[i] != 0) {
            return false;
        }
    }
    return true;
}

// For a binomial g = x^n + c (the reduction with its zero roots shifted
// out) and an odd prime q not dividing n, appends the roots of g in closed
// form (detail::binomial_roots_mod_prime), in ascending order, and returns
// true; returns false, appending nothing, for any other g or q or if the
// closed form fails, so the caller takes the generic route.
bool push_binomial_roots(std::vector<flint::Fmpz>& roots,
                         const NmodPoly& g,
                         ulong q,
                         const DiagnosticsContext* diagnostics) noexcept {
    std::vector<ulong> binomial_roots;
    if (!detail::binomial_roots_mod_prime(binomial_roots, g.raw()->coeffs,
                                          nmod_poly_degree(g.raw()), q)) {
        return false;
    }
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.stable_scan."
                        "degree_one_roots_nmod_binomial_roots");
    // A fixed order for the closed-form roots; the generic routes keep the
    // order FLINT's root splitting gives.
    std::sort(binomial_roots.begin(), binomial_roots.end());
    roots.reserve(roots.size() + binomial_roots.size());
    for (ulong root : binomial_roots) {
        roots.emplace_back();
        flint::fmpz_set_ui(flint::FmpzRef(roots.back()), root);
    }
    return true;
}

// The distinct roots of the monic integral `polynomial` modulo the
// word-size prime q.  For odd q, after the root 0, which occurs exactly
// when q divides the constant coefficient, write f = x^k g with g(0) != 0.
// For a binomial g = x^n + c with q not dividing n, the roots of g come in
// closed form, in ascending order (push_binomial_roots).  Otherwise this is
// the first step of Rabin's root finding, as FLINT's nmod_poly_roots
// (nmod_poly_factor/roots.c, _nmod_poly_push_roots) takes it, and the roots
// come in the order it finds them: write h = x^((q - 1) / 2) mod g.  The
// roots r of g have r^((q - 1) / 2) = 1 or -1, so they are the roots of
// gcd(h - 1, g) and of gcd(h + 1, g), each a squarefree product of linear
// factors (it divides x^((q - 1) / 2) -+ 1); a factor of degree one gives
// its root directly, a factor of degree two gives its roots in closed form
// in the order nmod_poly_roots would list them
// (detail::split_quadratic_roots_mod_prime), and a larger one goes to
// nmod_poly_roots.  For a
// binomial g = x^n + c that the closed form declines (q | n, or a failed
// step), h is the single term (-c)^k x^j with
// (q - 1) / 2 = k n + j, 0 <= j < n, since x^n = -c mod g
// (detail::binomial_x_power_mod).  Otherwise h comes
// from the word-size kernel word_x_power_mod
// (residue_field/word_frobenius_internal.hpp) for monic g (as here: the
// defining polynomial is monic integral and shifting out x^k keeps the
// leading coefficient 1) with deg g from 2 to 20 and q below 2^30
// (deg g <= 8) or 2^29, and from FLINT's
// nmod_poly_powmod_x_ui_preinv otherwise; x^((q - 1) / 2) mod g is unique,
// so every route gives the same h.  Every root found is checked to vanish
// on the reduction before it is returned.  With
// `require_squarefree`, returns false, with `roots` empty, when the
// reduction is not squarefree; also returns false, with `roots` empty, if a
// root fails the check.
bool degree_one_roots_mod_q_nmod(std::vector<flint::Fmpz>& roots,
                                 const fmpq_poly_struct* polynomial,
                                 ulong q,
                                 bool require_squarefree,
                                 const DiagnosticsContext* diagnostics)
        noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.stable_scan.degree_one_roots_nmod");
    NmodPoly reduced(q);
    NmodPolyFactor factorization;

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_nmod_reduce_poly");
        fmpq_poly_get_integral_nmod_poly(reduced, polynomial, q);
    }
    if (require_squarefree) {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_squarefree");
        if (nmod_poly_is_squarefree(reduced.raw()) == 0) {
            return false;
        }
    }

    if (q == 2) {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_nmod_factor");
        nmod_poly_roots(factorization.raw(), reduced.raw(), 0);
        push_nmod_linear_factor_roots(roots, factorization);
    } else {
        NmodPoly nonzero_part(q);
        NmodPoly reverse(q);
        NmodPoly preinverse(q);
        NmodPoly half_power(q);
        NmodPoly shifted(q);
        NmodPoly residue_part(q);
        NmodPoly nonresidue_part(q);

        slong zero_order = 0;
        while (zero_order < reduced.raw()->length &&
               nmod_poly_get_coeff_ui(reduced.raw(), zero_order) == 0) {
            ++zero_order;
        }
        if (zero_order > 0) {
            roots.emplace_back();
        }
        nmod_poly_shift_right(nonzero_part.raw(), reduced.raw(), zero_order);

        const bool binomial = nmod_poly_is_monic_binomial(nonzero_part);
        if (nmod_poly_degree(nonzero_part.raw()) >= 1 &&
            !(binomial &&
              push_binomial_roots(roots, nonzero_part, q, diagnostics))) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.stable_scan."
                        "degree_one_roots_nmod_frobenius_gcd");
                // h = x^((q - 1) / 2) mod g in closed form for a binomial
                // g = x^n + c (detail::binomial_x_power_mod), by the word
                // kernel where it covers the degree and q
                // (word_frobenius_internal.hpp), otherwise by FLINT; x^e
                // mod g is unique, so all three give the same reduced h.
                const slong nonzero_degree =
                        nmod_poly_degree(nonzero_part.raw());
                std::array<ulong, static_cast<std::size_t>(
                                          detail::kWordFrobeniusMaxDegree)>
                        word_power;
                nmod_poly_fit_length(half_power.raw(), nonzero_degree);
                if (binomial &&
                    detail::binomial_x_power_mod(
                            half_power.raw()->coeffs,
                            nonzero_part.raw()->coeffs, nonzero_degree,
                            (q - 1) / 2, q)) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.stable_scan."
                            "degree_one_roots_nmod_binomial_power");
                    _nmod_poly_set_length(half_power.raw(), nonzero_degree);
                    _nmod_poly_normalise(half_power.raw());
                } else if (detail::word_x_power_mod_dispatch(
                            word_power.data(), nonzero_part.raw()->coeffs,
                            nonzero_degree, (q - 1) / 2, q)) {
                    SILEX_PROFILE_EVENT(
                            diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.stable_scan."
                            "degree_one_roots_nmod_word_power");
                    for (slong i = 0; i < nonzero_degree; ++i) {
                        half_power.raw()->coeffs[i] =
                                word_power[static_cast<std::size_t>(i)];
                    }
                    _nmod_poly_set_length(half_power.raw(), nonzero_degree);
                    _nmod_poly_normalise(half_power.raw());
                } else {
                    nmod_poly_reverse(reverse.raw(), nonzero_part.raw(),
                                      nonzero_part.raw()->length);
                    nmod_poly_inv_series_newton(preinverse.raw(),
                                                reverse.raw(),
                                                reverse.raw()->length);
                    nmod_poly_powmod_x_ui_preinv(half_power.raw(),
                                                 (q - 1) / 2,
                                                 nonzero_part.raw(),
                                                 preinverse.raw());
                }
                nmod_poly_sub_ui(shifted.raw(), half_power.raw(), 1);
                nmod_poly_gcd(residue_part.raw(), shifted.raw(),
                              nonzero_part.raw());
                nmod_poly_add_ui(shifted.raw(), half_power.raw(), 1);
                nmod_poly_gcd(nonresidue_part.raw(), shifted.raw(),
                              nonzero_part.raw());
            }

            const NmodPoly* parts[] = {&residue_part, &nonresidue_part};
            for (const NmodPoly* part : parts) {
                const slong degree = nmod_poly_degree(part->raw());
                if (degree == 1) {
                    flint::Fmpz root;
                    if (degree_one_root_from_nmod_poly(root, part->raw())) {
                        roots.push_back(std::move(root));
                    }
                } else if (degree > 1) {
                    std::array<ulong, 2> quadratic_roots = {0, 0};
                    if (degree == 2 && part->raw()->coeffs[2] == 1 &&
                        detail::split_quadratic_roots_mod_prime(
                                quadratic_roots, part->raw()->coeffs[1],
                                part->raw()->coeffs[0], q)) {
                        SILEX_PROFILE_EVENT(
                                diagnostics, DiagnosticsModule::unit_group,
                                "unit_group.stable_scan."
                                "degree_one_roots_nmod_quadratic_roots");
                        for (ulong root : quadratic_roots) {
                            roots.emplace_back();
                            flint::fmpz_set_ui(flint::FmpzRef(roots.back()),
                                               root);
                        }
                        continue;
                    }
                    SILEX_PROFILE_SCOPE(
                            diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.stable_scan."
                            "degree_one_roots_nmod_factor");
                    nmod_poly_roots(factorization.raw(), part->raw(), 0);
                    push_nmod_linear_factor_roots(roots, factorization);
                }
            }
        }
    }

    if (!nmod_poly_vanishes_at_roots(reduced, roots)) {
        roots.clear();
        return false;
    }
    return true;
}

bool degree_one_roots_mod_q(std::vector<flint::Fmpz>& roots,
                            const Order& order,
                            flint::FmpzConstRef q,
                            const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.stable_scan.degree_one_roots_mod_q");
    roots.clear();
    const NumberField* field = order.parent();
    const nf_struct* raw_field = field == nullptr ? nullptr
                                                  : field->raw_flint_field();
    if (raw_field == nullptr || (!order.is_equation_order() &&
                                 !order.is_maximal()) ||
        !flint::fmpz_is_prime(q) ||
        !field_polynomial_is_monic_integral(raw_field->pol)) {
        return false;
    }

    // For equation orders, the caller has already skipped rational primes
    // dividing the order discriminant, so the reduced defining polynomial is
    // squarefree at this point.  Maximal orders that are not the equation
    // order keep the squarefree guard, since q may still divide the index
    // [O : Z[theta]].
    if (flint::fmpz_abs_fits_ui(q)) {
        const ulong q_ui = flint::fmpz_get_ui(q);
        if (q_ui > 1) {
            return degree_one_roots_mod_q_nmod(roots, raw_field->pol, q_ui,
                                               !order.is_equation_order(),
                                               diagnostics);
        }
    }

    flint::FmpzModCtx ctx(q.raw());
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPolyFactor factorization(ctx);
    if (!reduced.is_initialized()) {
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_reduce_poly");
        fmpq_poly_get_integral_fmpz_mod_poly(reduced, raw_field->pol, ctx);
    }
    {
        if (!order.is_equation_order()) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.degree_one_roots_squarefree");
            if (fmpz_mod_poly_is_squarefree(reduced.raw(), ctx.raw()) == 0) {
                return false;
            }
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_factor");
        fmpz_mod_poly_roots(factorization.raw(), reduced.raw(), 0, ctx.raw());
    }

    roots.reserve(static_cast<std::size_t>(factorization.raw()->num));
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        flint::Fmpz root;
        if (degree_one_root_from_mod_poly(
                    root, factorization.raw()->poly + i, q, ctx)) {
            roots.push_back(std::move(root));
        }
    }
    return true;
}

// The roots r of degree_one_roots_mod_q, one for each degree-one prime
// (q, theta - r) of the order over q, in the order degree_one_roots_mod_q
// lists them.  Returns false, so the caller decomposes q itself, when
// degree_one_roots_mod_q declines and, where q has at least two degree-one
// primes, for maximal quadratic orders, whose decomposition factors the
// minimal polynomial of the integral generator omega instead.
bool stable_scan_degree_one_roots(std::vector<flint::Fmpz>& roots,
                                  const Order& order,
                                  flint::FmpzConstRef q,
                                  const DiagnosticsContext* diagnostics)
        noexcept {
    if (!degree_one_roots_mod_q(roots, order, q, diagnostics)) {
        return false;
    }
    if (roots.size() >= 2 &&
        detail::residue_polynomial_uses_quadratic_integral_generator(order)) {
        return false;
    }
    return true;
}

// The proof columns of detail::stable_proof_prime_columns at q, built one
// at a time in the same order and passed to `visit_column` as each is
// built; the scan stops building them once `visit_column` returns false.
// Returns false, before any column is built, exactly where
// detail::stable_proof_prime_columns returns false (the scan then skips q);
// otherwise returns true, also when `visit_column` stopped the walk.  No
// column that is built can make the walk give up on q: a column that is
// unavailable is left out and the walk moves to the next prime over q.
// The direct columns go through `column_state`, which a caller may share
// across the values of q of one scan over `group` and `ell`; the columns
// are those of a new state either way.
template <typename VisitColumn>
bool visit_stable_proof_prime_columns(
        const OrderUnitGroup& group,
        flint::FmpzConstRef q,
        flint::FmpzConstRef ell,
        detail::DegreeOneRootColumnState& column_state,
        VisitColumn& visit_column) noexcept {
    const Order* order = group.parent();
    slong proof_rank = 0;
    if (!group.is_set() || order == nullptr ||
        !detail::dlog_proof_rank(proof_rank, group, ell)) {
        return false;
    }

    PrimeIdealList local;
    std::vector<flint::Fmpz> roots;
    if (stable_scan_degree_one_roots(roots, *order, q,
                                     group.diagnostics())) {
        if (roots.empty()) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.no_degree_one_prime");
            return true;
        }
        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_decomposition");
        bool decomposed = false;
        // The columns share the work that does not depend on the root (the
        // quotient-log setup and the generators' factors reduced mod q) and
        // the work that depends on neither the root nor q (the generators
        // and their distinct factors) through `column_state`.  It gives the
        // column and the failures of the overload without a state.
        for (const flint::Fmpz& root : roots) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.degree_one_prime");
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.dlog_column_test");
            flint::FmpzMat proof_column(proof_rank, 1);
            const bool available =
                    detail::saturation_proof_prime_column_at_degree_one_root(
                            proof_column, group, q, flint::FmpzConstRef(root),
                            ell, column_state);
            if (available) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.direct_dlog_column_available");
                if (!visit_column(proof_column)) {
                    return true;
                }
                continue;
            }

            // The direct evaluation declines, for example when a factor of a
            // generator vanishes or has a denominator at the root; the
            // decompose_prime prime over this root takes the residue-field
            // route the scan used before.
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.direct_dlog_column_unavailable");
            if (!decomposed) {
                decomposed = true;
                if (!decompose_prime(local, *order, q, 1)) {
                    // decompose_prime leaves `local` empty on failure.
                    SILEX_PROFILE_EVENT(
                            group.diagnostics(),
                            DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.decompose_failed");
                }
            }
            const PrimeIdeal* prime = nullptr;
            for (slong i = 0; i < local.size() && prime == nullptr; ++i) {
                const PrimeIdeal* candidate = local.at(i);
                flint::Fmpz candidate_root;
                flint::Fmpz candidate_q;
                if (candidate != nullptr &&
                    detail::degree_one_prime_root_mod_p(candidate_root, candidate_q,
                                                *candidate) &&
                    flint::fmpz_equal(flint::FmpzConstRef(candidate_q), q) &&
                    flint::fmpz_equal(flint::FmpzConstRef(candidate_root),
                                      flint::FmpzConstRef(root))) {
                    prime = candidate;
                }
            }
            if (prime == nullptr ||
                !detail::saturation_proof_prime_column(proof_column, group,
                                                       *prime, ell)) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.dlog_column_unavailable");
                continue;
            }
            if (!visit_column(proof_column)) {
                return true;
            }
        }
        return true;
    }

    if (!decompose_prime(local, *order, q, 1)) {
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.decompose_failed");
        return false;
    }
    if (local.size() == 0) {
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.no_degree_one_prime");
        return true;
    }

    SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.degree_one_decomposition");
    for (slong i = 0; i < local.size(); ++i) {
        const PrimeIdeal* prime = local.at(i);
        if (prime == nullptr) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.null_prime");
            continue;
        }
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.degree_one_prime");
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.dlog_column_test");
        flint::FmpzMat proof_column(proof_rank, 1);
        if (!detail::saturation_proof_prime_column(proof_column, group,
                                                   *prime, ell)) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.dlog_column_unavailable");
            continue;
        }
        if (!visit_column(proof_column)) {
            return true;
        }
    }
    return true;
}


bool order_discriminant_for_saturation(flint::Fmpz& out,
                                       const Order& order) noexcept {
    return order.discriminant(flint::FmpzRef(out)) &&
           !flint::fmpz_is_zero(flint::FmpzConstRef(out));
}

bool rational_prime_divides_order_discriminant(
        flint::FmpzConstRef discriminant,
        flint::FmpzConstRef q) noexcept {
    return flint::fmpz_divisible(discriminant, q);
}

#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
slong log_slong_from_fmpz(flint::FmpzConstRef value) noexcept {
    return flint::fmpz_fits_si(value) ? flint::fmpz_get_si(value) : -1;
}

void log_unit_proof_index_bound(const DiagnosticsContext* diagnostics,
                                flint::FmpzConstRef index_bound,
                                slong restarts,
                                bool changed) noexcept {
    if (!log_enabled(diagnostics, DiagnosticsModule::unit_group,
                     LogLevel::detail)) {
        return;
    }

    char detail[160];
    std::snprintf(detail, sizeof(detail),
                  "index_bound=%ld restarts=%ld changed=%ld",
                  static_cast<long>(log_slong_from_fmpz(index_bound)),
                  static_cast<long>(restarts),
                  changed ? 1L : 0L);
    log_emit(diagnostics, DiagnosticsModule::unit_group, LogLevel::detail,
             "prove_index_bound", "unit proof index bound", detail);
}

void log_unit_proof_restart(const DiagnosticsContext* diagnostics,
                            flint::FmpzConstRef ell,
                            flint::FmpzConstRef index_bound,
                            slong restart,
                            slong max_restarts,
                            slong carried_records) noexcept {
    if (!log_enabled(diagnostics, DiagnosticsModule::unit_group,
                     LogLevel::detail)) {
        return;
    }

    char detail[192];
    std::snprintf(detail, sizeof(detail),
                  "ell=%ld index_bound=%ld restart=%ld max_restarts=%ld "
                  "carried_records=%ld",
                  static_cast<long>(log_slong_from_fmpz(ell)),
                  static_cast<long>(log_slong_from_fmpz(index_bound)),
                  static_cast<long>(restart),
                  static_cast<long>(max_restarts),
                  static_cast<long>(carried_records));
    log_emit(diagnostics, DiagnosticsModule::unit_group, LogLevel::detail,
             "prove_index_bound", "unit proof accepted saturation root",
             detail);
}
#endif

struct RankOneProofSource {
    FactoredElement free_generator;
    OrderElement torsion_generator;
    bool has_free_generator = false;
    bool has_torsion_generator = false;
};

bool set_rank_one_proof_source(RankOneProofSource& out,
                               const OrderUnitGroup& group,
                               const Order& order) noexcept {
    out.free_generator.clear();
    out.torsion_generator.clear();
    out.has_free_generator = false;
    out.has_torsion_generator = false;

    if (group.free_rank() == 1) {
        const NumberField* field = order.parent();
        if (field == nullptr || !out.free_generator.define(*field) ||
            !group.free_generator(out.free_generator, 0)) {
            return false;
        }
        out.has_free_generator = true;
        return true;
    }

    if (group.free_rank() == 0) {
        if (!out.torsion_generator.define(order) ||
            !group.torsion_generator(out.torsion_generator)) {
            return false;
        }
        out.has_torsion_generator = true;
        return true;
    }

    return false;
}

bool rank_one_proof_source_nonzero(bool& nonzero,
                                   const RankOneProofSource& source,
                                   const OrderUnitGroup& group,
                                   const PrimeIdeal& prime,
                                   flint::FmpzConstRef ell) noexcept {
    if (source.has_free_generator) {
        return detail::saturation_proof_prime_known_rank_one_free_generator_nonzero(
                nonzero, group, source.free_generator, prime, ell);
    }
    if (source.has_torsion_generator) {
        return detail::saturation_proof_prime_known_rank_one_torsion_nonzero(
                nonzero, group, source.torsion_generator, prime, ell);
    }
    return false;
}

enum class ProofPrimeScanResult {
    keep_scanning,
    done,
    failed,
    exhausted,
};

// Resource guard, not a proof parameter: the most rational primes
// q = 1 mod ell that one `select_stable_proof_primes` scan examines.  The
// reference `compute_candidates_for_saturate` iterates `PrimesSet` without an
// upper end; a scan that reaches this count stops and the ell-local proof
// fails closed as `unavailable`.  The guard is per scan.  A scan counts at
// most d degree-one primes per q (d = field degree), so once the doubled
// `saturate!` stability threshold exceeds d * 2^16 a scan can end only with
// an empty kernel or at the guard: one ell runs at most 17 + log2 d scans
// that adjoin no root.  Scans that adjoin a root are bounded by the ell-part
// of the finite index [O^x : U] (and by the caller's restart limit), so one
// ell examines at most about (17 + log2 d + log_ell [O^x : U]) * 2^16 primes.
// This also bounds the ell >= aux_bound scans, which previously ran to the
// end of the machine-word range of q.
inline constexpr slong kStableProofMaxAuxPrimes = slong{1} << 16;

// `range_exhausted` is set when the scan returns a nonempty kernel because no
// rational prime q = 1 mod ell is left in the machine-word range, so a later
// scan with a larger threshold would return the same kernel.
bool select_stable_proof_primes(
        slong& local_prime_count,
        bool& certified,
        bool& range_exhausted,
        flint::FmpzMat& kernel,
        const OrderUnitGroup& group,
        flint::FmpzConstRef ell,
        double stable_threshold) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_proof_primes");
    range_exhausted = false;
    const Order* order = group.parent();
    slong proof_rank = 0;
    if (!group.is_set() || order == nullptr ||
        !detail::dlog_proof_rank(proof_rank, group, ell) ||
        stable_threshold <= 0.0) {
        return false;
    }
    local_prime_count = 0;

    if (proof_rank == 0) {
        flint::FmpzMat candidate_kernel(0, 0);
        kernel = std::move(candidate_kernel);
        certified = true;
        return true;
    }

    RankOneProofSource rank_one_source;
    if (proof_rank == 1 &&
        !set_rank_one_proof_source(rank_one_source, group, *order)) {
        return false;
    }

    flint::Fmpz order_discriminant;
    if (!order_discriminant_for_saturation(order_discriminant, *order)) {
        return false;
    }

    flint::FmpzMat proof_matrix(proof_rank, 0);
    flint::FmpzMat current_kernel(0, proof_rank);
    slong current_dimension = proof_rank;
    // reference starts the unchanged-kernel counter at `i = 1` for the initial
    // identity candidate matrix in `compute_candidates_for_saturate`.
    slong stable_count = 1;
    bool have_kernel = false;

    auto process_rank_one_image =
            [&](bool nonzero) noexcept -> ProofPrimeScanResult {
        ++local_prime_count;
        if (nonzero) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.rank_one_image_nonzero");
            flint::FmpzMat candidate_kernel(0, 1);
            kernel = std::move(candidate_kernel);
            certified = true;
            return ProofPrimeScanResult::done;
        }

        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.rank_one_image_zero");
        have_kernel = true;
        set_rank_one_candidate_kernel(current_kernel);
        ++stable_count;
        if (static_cast<double>(stable_count) >
            stable_threshold * static_cast<double>(current_dimension)) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.stable_threshold");
            kernel = std::move(current_kernel);
            certified = false;
            return ProofPrimeScanResult::done;
        }
        return ProofPrimeScanResult::keep_scanning;
    };

    auto process_direct_rank_one_q =
            [&](flint::FmpzConstRef q,
                ProofPrimeScanResult& result) noexcept -> bool {
        SILEX_PROFILE_SCOPE(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.stable_scan.direct_rank_one_q");
        result = ProofPrimeScanResult::keep_scanning;
        if (!rank_one_source.has_free_generator) {
            return false;
        }

        std::vector<flint::Fmpz> roots;
        if (!degree_one_roots_mod_q(roots, *order, q, group.diagnostics())) {
            return false;
        }

        if (roots.empty()) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.no_degree_one_prime");
            return true;
        }

        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_decomposition");
        bool saw_usable_root = false;
        for (const flint::Fmpz& root : roots) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.degree_one_prime");
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.rank_one_image_test");
            bool nonzero = false;
            {
                SILEX_PROFILE_SCOPE(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.rank_one_image_eval");
                bool image_available = false;
                if (flint::fmpz_abs_fits_ui(q) &&
                    flint::fmpz_abs_fits_ui(flint::FmpzConstRef(root))) {
                    image_available =
                            detail::rank_one_factored_image_nonzero_at_degree_one_root_nmod(
                                    nonzero, rank_one_source.free_generator,
                                    flint::fmpz_get_ui(q),
                                    flint::fmpz_get_ui(
                                            flint::FmpzConstRef(root)),
                                    ell);
                }
                if (!image_available) {
                    image_available =
                            detail::rank_one_factored_image_nonzero_at_degree_one_root(
                                    nonzero, rank_one_source.free_generator, q,
                                    flint::FmpzConstRef(root), ell);
                }
                if (!image_available) {
                    SILEX_PROFILE_EVENT(
                            group.diagnostics(),
                            DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.rank_one_image_unavailable");
                    continue;
                }
            }

            saw_usable_root = true;
            result = process_rank_one_image(nonzero);
            if (result != ProofPrimeScanResult::keep_scanning) {
                return true;
            }
        }

        return saw_usable_root;
    };

    // reference RelSaturate.compute_candidates_for_saturate scans q = 1 mod p
    // and stops once the kernel dimension has stayed the same for more than
    // `stable` times the number of input elements (`saturate!` passes
    // stable = 3.5 and doubles it after a candidate that is not a p-th
    // power).  This scan multiplies by the current kernel dimension instead,
    // and starts at q = ell + 1 where the reference `PrimesSet` starts at
    // `p_start` = next_prime(2^60).  Any degree-one prime q = 1 mod ell not
    // dividing disc(O) gives a valid character, so neither choice changes the
    // proof condition (an empty kernel).  No auxiliary bound limits the scan;
    // only the kStableProofMaxAuxPrimes resource guard does.
    const slong max_aux_primes_override =
            detail::StableProofTestSeams::max_aux_primes();
    const slong max_aux_primes = max_aux_primes_override > 0
                                         ? max_aux_primes_override
                                         : kStableProofMaxAuxPrimes;
    slong examined_primes = 0;
    // One column state for the whole scan: `group` and `ell` are fixed, so
    // the generators and their distinct factors are found once and only
    // the parts that depend on q are redone at each q.
    detail::DegreeOneRootColumnState column_state;
    // The proof columns at the degree-one primes over q, in the order of
    // detail::stable_proof_prime_columns, each appended to the
    // proof matrix as soon as it is built.  The walk over the primes above q
    // stops at the column that ends the scan (full rank, the stability
    // threshold, or a failure), so the columns of the remaining primes over
    // q are never built.
    //
    // Why stopping early leaves the certificate unchanged: the scan is a
    // prefix process.  Each appended column updates the kernel, the
    // stability counter and local_prime_count from the columns before it
    // alone, and the scan returns at the first column that reaches full
    // rank or the threshold; the columns after it were built only to be
    // discarded, so the published kernel, certified flag and
    // local_prime_count depend on the same column prefix either way.
    // Building those columns cannot change that prefix either: whether q is
    // skipped is decided before its first column is built (see
    // visit_stable_proof_prime_columns), and a column that turns out to be
    // unavailable is left out without affecting the others.  Only the
    // diagnostics counters of the unbuilt columns are no longer recorded.
    auto process_dlog_columns =
            [&](flint::FmpzConstRef q) noexcept -> ProofPrimeScanResult {
        ProofPrimeScanResult result = ProofPrimeScanResult::keep_scanning;
        auto process_column =
                [&](const flint::FmpzMat& proof_column) noexcept -> bool {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.dlog_column");
            ++local_prime_count;
            if (!append_column(proof_matrix, proof_column)) {
                result = ProofPrimeScanResult::failed;
                return false;
            }
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_dlog_kernel_nullspace");
            if (!detail::dlog_kernel_from_matrix(current_kernel, proof_matrix,
                                                 ell)) {
                result = ProofPrimeScanResult::failed;
                return false;
            }

            if (flint::fmpz_mat_nrows(current_kernel) == 0) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.certified_full_rank");
                kernel = std::move(current_kernel);
                certified = true;
                result = ProofPrimeScanResult::done;
                return false;
            }

            have_kernel = true;
            const slong next_dimension = flint::fmpz_mat_nrows(current_kernel);
            if (next_dimension == current_dimension) {
                ++stable_count;
            } else {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.kernel_dimension_changed");
                current_dimension = next_dimension;
                stable_count = 0;
            }
            if (static_cast<double>(stable_count) >
                stable_threshold * static_cast<double>(current_dimension)) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.stable_threshold");
                kernel = std::move(current_kernel);
                certified = false;
                result = ProofPrimeScanResult::done;
                return false;
            }
            return true;
        };
        if (!visit_stable_proof_prime_columns(group, q, ell, column_state,
                                              process_column)) {
            return ProofPrimeScanResult::keep_scanning;
        }
        return result;
    };
    auto process_prime_q =
            [&](flint::FmpzConstRef q) noexcept -> ProofPrimeScanResult {
        if (examined_primes >= max_aux_primes) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.resource_limit");
            return ProofPrimeScanResult::exhausted;
        }
        ++examined_primes;
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.candidate_q");
        if (rational_prime_divides_order_discriminant(
                    flint::FmpzConstRef(order_discriminant), q)) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.discriminant_skip");
            return ProofPrimeScanResult::keep_scanning;
        }
        if (proof_rank == 1 && rank_one_source.has_free_generator) {
            ProofPrimeScanResult direct_result =
                    ProofPrimeScanResult::keep_scanning;
            if (process_direct_rank_one_q(q, direct_result)) {
                return direct_result;
            }
        }

        if (proof_rank != 1) {
            return process_dlog_columns(q);
        }

        // proof_rank == 1 without a usable direct rank-one image.
        PrimeIdealList local;
        if (!decompose_prime(local, *order, q, 1)) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.decompose_failed");
            return ProofPrimeScanResult::keep_scanning;
        }

        if (local.size() == 0) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.no_degree_one_prime");
            return ProofPrimeScanResult::keep_scanning;
        }

        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_decomposition");
        for (slong i = 0; i < local.size(); ++i) {
            const PrimeIdeal* prime = local.at(i);
            if (prime == nullptr) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.null_prime");
                continue;
            }
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.degree_one_prime");

            bool nonzero = false;
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.rank_one_image_test");
            if (!rank_one_proof_source_nonzero(
                        nonzero, rank_one_source, group, *prime, ell)) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.rank_one_image_unavailable");
                continue;
            }

            const ProofPrimeScanResult result =
                    process_rank_one_image(nonzero);
            if (result != ProofPrimeScanResult::keep_scanning) {
                return result;
            }
        }

        return ProofPrimeScanResult::keep_scanning;
    };

    const ulong max_q_override =
            detail::StableProofTestSeams::max_q();
    const ulong max_q_ui = max_q_override > 0
                                   ? max_q_override
                                   : std::numeric_limits<ulong>::max();
    bool used_machine_scan = false;
    if (flint::fmpz_abs_fits_ui(ell)) {
        const ulong ell_ui = flint::fmpz_get_ui(ell);
        if (ell_ui > 0 && ell_ui < max_q_ui) {
            used_machine_scan = true;
            flint::Fmpz q_word;
            SILEX_PROFILE_SCOPE(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.machine_prime_scan");
            for (ulong q_ui = ell_ui + 1;;) {
                if (n_is_prime(q_ui) != 0) {
                    flint::fmpz_set_ui(flint::FmpzRef(q_word), q_ui);
                    const ProofPrimeScanResult result =
                            process_prime_q(flint::FmpzConstRef(q_word));
                    if (result == ProofPrimeScanResult::done) {
                        return true;
                    }
                    if (result == ProofPrimeScanResult::failed ||
                        result == ProofPrimeScanResult::exhausted) {
                        return false;
                    }
                }
                if (q_ui > max_q_ui - ell_ui) {
                    break;
                }
                q_ui += ell_ui;
            }
        }
    }

    if (!used_machine_scan) {
        flint::Fmpz q;
        flint::Fmpz max_q;
        flint::fmpz_add_ui(flint::FmpzRef(q), ell, 1);
        flint::fmpz_set_ui(flint::FmpzRef(max_q), max_q_ui);

        SILEX_PROFILE_SCOPE(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.stable_scan.fmpz_prime_scan");
        while (flint::fmpz_cmp(flint::FmpzConstRef(q),
                               flint::FmpzConstRef(max_q)) <= 0) {
            if (flint::fmpz_is_prime(flint::FmpzConstRef(q))) {
                const ProofPrimeScanResult result =
                        process_prime_q(flint::FmpzConstRef(q));
                if (result == ProofPrimeScanResult::done) {
                    return true;
                }
                if (result == ProofPrimeScanResult::failed ||
                    result == ProofPrimeScanResult::exhausted) {
                    return false;
                }
            }
            flint::fmpz_add(flint::FmpzRef(q), flint::FmpzConstRef(q), ell);
        }
    }

    if (!have_kernel) {
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.no_kernel");
        return false;
    }

    kernel = std::move(current_kernel);
    certified = false;
    range_exhausted = true;
    SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.exhausted_with_kernel");
    return true;
}

enum class StableRelationSaturationStatus {
    verified,
    accepted_root,
    wasted_non_power,
    unavailable,
};

struct StableRelationSaturationStep {
    StableRelationSaturationStatus status =
            StableRelationSaturationStatus::unavailable;
    slong local_primes = 0;
};

bool stable_relation_saturation_step(OrderUnitGroup& out,
                                  StableRelationSaturationStep& step,
                                  const OrderUnitGroup& group,
                                  flint::FmpzConstRef ell,
                                  double stable_threshold,
                                  bool defer_undecidable_roots,
                                  EmbeddingContext& embeddings,
                                  slong precision,
                                  detail::CompactFieldModulusCache&
                                          field_modulus_cache) noexcept {
    step = StableRelationSaturationStep{};
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!out.is_defined() || !group.is_set() ||
        !same_order_parent(out.parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || precision <= 0 ||
        stable_threshold <= 0.0) {
        return false;
    }

    slong local_prime_count = 0;
    bool certified = false;
    bool range_exhausted = false;
    flint::FmpzMat kernel(0, rank);
    SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_prime_selection");
    if (!select_stable_proof_primes(
                local_prime_count, certified, range_exhausted, kernel, group,
                ell, stable_threshold)) {
        step.status = StableRelationSaturationStatus::unavailable;
        return true;
    }
    step.local_primes = local_prime_count;

    if (flint::fmpz_mat_nrows(kernel) == 0) {
        step.status = StableRelationSaturationStatus::verified;
        return true;
    }

    bool saw_candidate = false;
    // A kernel row whose l-th root lies in O_K but not in the order O stays
    // in the residue-character kernel: at every degree-one prime p of O not
    // dividing disc(O) the root reduces into O/p, so the row remains an l-th
    // power there.  More primes cannot remove it, and the characters cannot
    // tell it from a row whose root lies in O, so this proof cannot certify
    // l-saturation of the group in O^x.  Other rows may still give roots in
    // O, so keep scanning, but never report `verified` or retry with a
    // larger threshold once such a row was seen.
    bool saw_root_outside_order = false;
    // Until the first doubling of `stable` for this ell (the callers pass
    // defer_undecidable_roots), a row of exponent ell >= 3 in a field of
    // degree >= 10 that is not a formal ell-th power is deferred, not tested
    // (see kernel_row_root).  The deferred row is handled like an undecided one
    // below, so `stable` doubles and the rescan uses more primes q; the
    // rescan's rows are tested as usual.  A true ell-th power stays in every
    // kernel, so it still reaches the root test at a later scan; a
    // non-power that more primes remove costs no root test.  The root is
    // then accepted at a larger `stable`, possibly from another kernel row,
    // so the published generators and the proof record's local primes may
    // differ from those of a first-scan attempt.  When the range is
    // exhausted a rescan would return the same kernel, so the rows are
    // tested at once, as before.
    const detail::KernelRowRootMode root_mode =
            defer_undecidable_roots && !range_exhausted
                    ? detail::KernelRowRootMode::defer_undecidable
                    : detail::KernelRowRootMode::attempt;
    for (slong row = 0; row < flint::fmpz_mat_nrows(kernel); ++row) {
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_relation_saturation.row_inspected");
        // The kernel has the proof width: the free exponents and, when l
        // divides the torsion order, the torsion exponent.  Reference
        // `saturate!` forms every candidate as prod R[j]^e[j] * zeta^e[end];
        // a row whose free part alone is divisible by l still carries a
        // torsion factor and is tested like any other row.
        if (detail::kernel_row_divisible(
                    flint::FmpzMatConstRef(kernel), row,
                    flint::fmpz_mat_ncols(kernel), ell)) {
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.row_divisible_skipped");
            continue;
        }
        saw_candidate = true;
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_relation_saturation.row_candidate");

        FactoredElement root(*field);
        bool is_power = false;
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.saturation_row_root");
        const detail::KernelRowRootResult root_result =
                detail::kernel_row_root(is_power, root, group,
                                        flint::FmpzMatConstRef(kernel), row,
                                        ell, &field_modulus_cache, root_mode);
        if (root_result == detail::KernelRowRootResult::error) {
            // An internal error, not an undecided root: a rescan cannot
            // help, so fail closed at once.
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_error");
            return false;
        }
        if (root_result == detail::KernelRowRootResult::deferred) {
            // Not tested (see root_mode above); `range_exhausted` is false
            // here, so this rescans exactly as an undecided row does.
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_deferred");
            step.status = saw_root_outside_order
                                  ? StableRelationSaturationStatus::unavailable
                                  : StableRelationSaturationStatus::
                                            wasted_non_power;
            return true;
        }
        if (root_result == detail::KernelRowRootResult::undecided) {
            // No exactly verified root: the compact presentation and the
            // factored fallback both left the row undecided.  The pinned
            // reference `saturate!` (Saturate.jl:396-411) treats every
            // candidate for which `is_power` returns no root as a non-power:
            // it stops testing candidates, doubles `stable` and rescans with
            // more primes q = 1 mod ell.  Do the same.  This never certifies
            // the row: ell is verified only by a later scan with an empty
            // kernel, or the group grows only by an exactly verified root.
            // A row that is an ell-th power whose root cannot be computed
            // stays in every kernel, so the rescans end at the
            // kStableProofMaxAuxPrimes guard and ell fails closed as
            // `unavailable`.  When the scan has already used every prime of
            // the machine-word range, a rescan would return the same kernel,
            // so fail closed at once.
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_unavailable");
            step.status = saw_root_outside_order || range_exhausted
                                  ? StableRelationSaturationStatus::unavailable
                                  : StableRelationSaturationStatus::
                                            wasted_non_power;
            return true;
        }
        if (!is_power) {
            // As for an undecided row: once the scan has used every prime
            // of the machine-word range, a rescan with a larger threshold
            // would return the same kernel and the same non-power, so fail
            // closed at once instead of doubling `stable` forever.
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_not_power");
            step.status = saw_root_outside_order || range_exhausted
                                  ? StableRelationSaturationStatus::unavailable
                                  : StableRelationSaturationStatus::
                                            wasted_non_power;
            return true;
        }
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_relation_saturation.root_power");
        if (!detail::adjoinable_root_is_order_unit(*order, root)) {
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_outside_order");
            saw_root_outside_order = true;
            continue;
        }

        flint::FmpzMat relation(1, rank);
        for (slong i = 0; i < rank; ++i) {
            flint::fmpz_set(flint::fmpz_mat_entry(relation, 0, i),
                            flint::FmpzConstRef(
                                    flint::fmpz_mat_entry(kernel, row, i).raw()));
        }

        std::vector<FactoredElement> generators;
        if (!detail::relation_basis(generators, group, root,
                                    flint::FmpzMatConstRef(relation), ell)) {
            return false;
        }

        OrderUnitGroup candidate(*order);
        candidate.set_diagnostics(out.diagnostics());
        if (!candidate.is_defined() ||
            !detail::order_unit_group_set_units_internal(
                    candidate, *order,
                    FactoredElementSpan(generators.data(), generators.size()),
                    embeddings, precision, true)) {
            return false;
        }

        bool reduced = false;
        if (!detail::reduce_stored_relation_units(reduced, candidate, embeddings,
                                               precision)) {
            return false;
        }

        step.status = StableRelationSaturationStatus::accepted_root;
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_relation_saturation.accepted_root");
        out.swap(candidate);
        return true;
    }

    SILEX_PROFILE_EVENT(
            out.diagnostics(), DiagnosticsModule::unit_group,
            saw_candidate
                    ? "unit_group.stable_relation_saturation.exhausted_candidates"
                    : "unit_group.stable_relation_saturation.no_candidate");
    step.status = saw_candidate && !saw_root_outside_order
                          ? StableRelationSaturationStatus::wasted_non_power
                          : StableRelationSaturationStatus::unavailable;
    return true;
}

bool copy_selected_primes_to_list(PrimeIdealList& out,
                                  const Order& order,
                                  const std::vector<PrimeIdeal>& selected,
                                  slong target_len) noexcept {
    if (static_cast<slong>(selected.size()) != target_len) {
        return false;
    }

    PrimeIdealList candidate(order, target_len);
    if (!candidate.is_defined()) {
        return false;
    }
    for (slong i = 0; i < target_len; ++i) {
        PrimeIdeal* dest = candidate.at(i);
        if (dest == nullptr ||
            !dest->set(selected[static_cast<std::size_t>(i)])) {
            return false;
        }
    }

    out.swap(candidate);
    return true;
}

bool publish_saturation_selection(PrimeIdealList& out,
                                  flint::FmpzMat& kernel,
                                  const Order& order,
                                  const std::vector<PrimeIdeal>& selected,
                                  const flint::FmpzMat& dlog_matrix,
                                  slong target_len,
                                  flint::FmpzConstRef ell) noexcept {
    if (!copy_selected_primes_to_list(out, order, selected, target_len)) {
        return false;
    }
    return detail::dlog_kernel_from_matrix(kernel, dlog_matrix, ell);
}

// One scan over rational primes p = 2, 3, ... <= bound, shared by
// OrderUnitGroup::select_saturation_primes and
// detail::select_saturation_primes_with_kernel.  It stops once `selected`
// holds target_len primes and offers each candidate prime ideal to
// `append_prime`, which decides whether to keep it and returns false only on
// a hard failure.
//
// First pass (source_shape_fallback false): reference
// RelSaturate.compute_candidates_for_saturate scans rational q = 1 mod ell,
// skips q dividing the discriminant, and uses only degree-one primes above
// q.  Deliberate difference: the reference also skips q dividing
// [O_K : Z[alpha]] with q not dividing d_K for non-monogenic maximal orders;
// here only q | disc(O) is skipped.  Selection only, not soundness.
//
// Fallback pass (source_shape_fallback true): every prime ideal above every
// rational prime, with no congruence or discriminant filter.
template <typename AppendPrime>
bool scan_saturation_candidate_primes(
        std::vector<PrimeIdeal>& selected,
        const Order& order,
        const DiagnosticsContext* diagnostics,
        flint::FmpzConstRef order_discriminant,
        flint::FmpzConstRef ell,
        slong target_len,
        flint::FmpzConstRef bound,
        bool source_shape_fallback,
        AppendPrime& append_prime) noexcept {
    flint::Fmpz p;
    flint::Fmpz pminus;
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (static_cast<slong>(selected.size()) < target_len &&
           fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList local;
        if (source_shape_fallback) {
            if (!decompose_prime(local, order, flint::FmpzConstRef(p))) {
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
        } else {
            fmpz_sub_ui(pminus.raw(), p.raw(), 1);
            if (fmpz_divisible(pminus.raw(), ell.raw()) == 0) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.saturation_selector.congruence_skip");
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
            if (rational_prime_divides_order_discriminant(
                        order_discriminant, flint::FmpzConstRef(p))) {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.saturation_selector.discriminant_skip");
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::unit_group,
                    "unit_group.saturation_selector_decompose_degree_one");
            if (!decompose_prime(local, order, flint::FmpzConstRef(p), 1)) {
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
        }
        for (slong i = 0; i < local.size() &&
                          static_cast<slong>(selected.size()) < target_len;
             ++i) {
            const PrimeIdeal* prime = local.at(i);
            if (prime != nullptr && !append_prime(*prime)) {
                return false;
            }
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }
    return true;
}

}  // namespace

namespace detail {

bool select_saturation_primes_with_kernel(
        PrimeIdealList& out,
        flint::FmpzMat& kernel,
        const OrderUnitGroup& group,
        flint::FmpzConstRef ell,
        slong target_len,
        flint::FmpzConstRef bound) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.select_saturation_primes_with_kernel");
    const Order* order = group.parent();
    if (!group.is_set() || order == nullptr || !flint::fmpz_is_prime(ell) ||
        target_len < 0 || fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    // Character columns carry the torsion generator when l divides the
    // torsion order, as reference `compute_candidates_for_saturate`.
    slong width = 0;
    if (!detail::dlog_proof_rank(width, group, ell)) {
        return false;
    }

    if (target_len == 0) {
        // No characters: the common kernel is the whole exponent space, so
        // its basis is the width x width identity (rows are kernel vectors,
        // as dlog_kernel_from_matrix publishes them).
        PrimeIdealList candidate(*order, 0);
        flint::FmpzMat candidate_kernel(width, width);
        if (!candidate.is_defined()) {
            return false;
        }
        flint::fmpz_mat_one(flint::FmpzMatRef(candidate_kernel));
        out.swap(candidate);
        kernel = std::move(candidate_kernel);
        return true;
    }
    std::vector<PrimeIdeal> selected;
    selected.reserve(static_cast<std::size_t>(target_len));
    flint::FmpzMat dlog_matrix(width, 0);

    flint::Fmpz order_discriminant;
    if (!order_discriminant_for_saturation(order_discriminant, *order)) {
        return false;
    }

    auto try_append_prime = [&](const PrimeIdeal& prime) noexcept -> bool {
        flint::FmpzMat column(width, 1);
        if (!detail::saturation_prime_character_column(column, group, prime,
                                                       ell)) {
            return true;
        }
        selected.emplace_back(*order);
        if (!selected.back().set(prime) ||
            !append_column(dlog_matrix, column)) {
            return false;
        }
        return true;
    };

    if (!scan_saturation_candidate_primes(
                selected, *order, group.diagnostics(),
                flint::FmpzConstRef(order_discriminant), ell, target_len,
                bound, false, try_append_prime)) {
        return false;
    }
    if (static_cast<slong>(selected.size()) == target_len) {
        return publish_saturation_selection(out, kernel, *order, selected,
                                            dlog_matrix, target_len, ell);
    }

    SILEX_PROFILE_EVENT(
            group.diagnostics(), DiagnosticsModule::unit_group,
            "unit_group.saturation_selector.source_shape_fallback");
    selected.clear();
    flint::FmpzMat fallback_matrix(width, 0);
    dlog_matrix = std::move(fallback_matrix);
    if (!scan_saturation_candidate_primes(
                selected, *order, group.diagnostics(),
                flint::FmpzConstRef(order_discriminant), ell, target_len,
                bound, true, try_append_prime)) {
        return false;
    }
    return publish_saturation_selection(out, kernel, *order, selected,
                                        dlog_matrix, target_len, ell);
}

bool stable_proof_prime_columns(std::vector<flint::FmpzMat>& columns,
                                const OrderUnitGroup& group,
                                flint::FmpzConstRef q,
                                flint::FmpzConstRef ell) noexcept {
    columns.clear();
    auto collect = [&](flint::FmpzMat& column) noexcept -> bool {
        columns.push_back(std::move(column));
        return true;
    };
    detail::DegreeOneRootColumnState column_state;
    return visit_stable_proof_prime_columns(group, q, ell, column_state,
                                            collect);
}

bool binomial_x_power_mod(ulong* power,
                          const ulong* f,
                          slong degree,
                          ulong exponent,
                          ulong modulus) noexcept {
    if (degree < 1 || modulus == 0 || f[degree] != 1) {
        return false;
    }
    for (slong i = 1; i < degree; ++i) {
        if (f[i] != 0) {
            return false;
        }
    }
    nmod_t mod;
    nmod_init(&mod, modulus);
    const ulong udegree = static_cast<ulong>(degree);
    for (slong i = 0; i < degree; ++i) {
        power[i] = 0;
    }
    power[exponent % udegree] = nmod_pow_ui(nmod_neg(f[0], mod),
                                            exponent / udegree, mod);
    return true;
}

namespace {

// A y with y^d = a in F_q^*, for an odd prime q, d >= 2 dividing m = q - 1,
// and a with a^(m / d) = 1, after Adleman, Manders and Miller ("On taking
// roots in finite fields", FOCS 1977) with the correction term found by
// Pohlig and Hellman's prime-power discrete logarithm (IEEE Trans. Inform.
// Theory 24 (1978), 106-110); see also Bach and Shallit, Algorithmic Number
// Theory, vol. 1, chapter 7.  Write m = dp * t with dp the part of m made of
// the primes dividing d, so gcd(d, t) = 1, and let u = d^-1 mod t.  Then
// x1 = a^u has x1^d = a b with b = a^(d u - 1) in the cyclic subgroup H of
// order dp, generated by h = gamma^t for the primitive root gamma.  Its
// logarithm e (b = h^e) is found modulo each prime power r^s || dp digit by
// digit and combined by the Chinese remainder theorem; as a is a d-th power,
// d | e, and y = x1 h^(-e / d).  Returns false if a step fails (y^d = a is
// checked at the end), so the caller falls back to the generic route.
bool word_power_residue_root(ulong& y,
                             ulong a,
                             ulong d,
                             ulong gamma,
                             nmod_t mod) noexcept {
    const ulong m = mod.n - 1;
    ulong t = m;
    for (ulong common = n_gcd(t, d); common > 1; common = n_gcd(t, d)) {
        t /= common;
    }
    const ulong dp = m / t;
    const ulong u = t == 1 ? 0 : n_invmod(d % t, t);
    const ulong x1 = nmod_pow_ui(a, u, mod);
    const ulong b = nmod_mul(nmod_pow_ui(x1, d, mod), nmod_inv(a, mod), mod);
    const ulong h = nmod_pow_ui(gamma, t, mod);

    n_factor_t primes;
    n_factor_init(&primes);
    n_factor(&primes, d, 1);
    ulong e = 0;
    ulong e_modulus = 1;
    for (int i = 0; i < primes.num; ++i) {
        const ulong r = primes.p[i];
        ulong prime_power = 1;
        int s = 0;
        for (ulong rest = dp; rest % r == 0; rest /= r) {
            prime_power *= r;
            ++s;
        }
        const ulong cofactor = dp / prime_power;
        const ulong h_r = nmod_pow_ui(h, cofactor, mod);
        const ulong h_r_inverse = nmod_inv(h_r, mod);
        const ulong b_r = nmod_pow_ui(b, cofactor, mod);
        const ulong zeta_r = nmod_pow_ui(h_r, prime_power / r, mod);
        ulong e_r = 0;
        ulong digit_weight = 1;
        ulong top_power = prime_power / r;
        for (int k = 0; k < s; ++k) {
            const ulong w = nmod_pow_ui(
                    nmod_mul(b_r, nmod_pow_ui(h_r_inverse, e_r, mod), mod),
                    top_power, mod);
            ulong digit = 0;
            ulong zeta_power = 1;
            while (digit < r && zeta_power != w) {
                zeta_power = nmod_mul(zeta_power, zeta_r, mod);
                ++digit;
            }
            if (digit == r) {
                return false;
            }
            e_r += digit * digit_weight;
            digit_weight *= r;
            top_power /= r;
        }
        // e = e mod e_modulus and e = e_r mod prime_power.
        const ulong difference = n_submod(e_r, e % prime_power, prime_power);
        const ulong lift = n_mulmod2(
                difference, n_invmod(e_modulus % prime_power, prime_power),
                prime_power);
        e += e_modulus * lift;
        e_modulus *= prime_power;
    }
    if (e % d != 0) {
        return false;
    }
    const ulong c = nmod_pow_ui(h, (dp - e / d) % dp, mod);
    y = nmod_mul(x1, c, mod);
    return nmod_pow_ui(y, d, mod) == a;
}

}  // namespace

bool binomial_roots_mod_prime(std::vector<ulong>& roots,
                              const ulong* f,
                              slong degree,
                              ulong q) noexcept {
    roots.clear();
    if (degree < 1 || q < 3 || q % 2 == 0 || f[degree] != 1 || f[0] == 0 ||
        static_cast<ulong>(degree) % q == 0) {
        return false;
    }
    for (slong i = 1; i < degree; ++i) {
        if (f[i] != 0) {
            return false;
        }
    }
    nmod_t mod;
    nmod_init(&mod, q);
    const ulong n = static_cast<ulong>(degree);
    const ulong m = q - 1;
    const ulong a = nmod_neg(f[0], mod);
    const ulong d = n_gcd(n, m);
    if (d > 1 && nmod_pow_ui(a, m / d, mod) != 1) {
        return true;
    }

    ulong gamma = 0;
    ulong y = a;
    if (d > 1) {
        gamma = primitive_root_mod_prime_ui(q);
        if (gamma == 0 || !word_power_residue_root(y, a, d, gamma, mod)) {
            return false;
        }
    }
    // x0 = y^v with v (n / d) = 1 mod (m / d) has x0^n = a^(v n / d) = a,
    // as a has order dividing m / d.
    const ulong reduced_modulus = m / d;
    const ulong v = reduced_modulus == 1
            ? 0
            : n_invmod((n / d) % reduced_modulus, reduced_modulus);
    const ulong x0 = nmod_pow_ui(y, v, mod);
    if (nmod_pow_ui(x0, n, mod) != a) {
        return false;
    }
    const ulong zeta = d > 1 ? nmod_pow_ui(gamma, m / d, mod) : 1;
    roots.reserve(static_cast<std::size_t>(d));
    ulong root = x0;
    for (ulong i = 0; i < d; ++i) {
        roots.push_back(root);
        root = nmod_mul(root, zeta, mod);
    }
    return true;
}

bool split_quadratic_roots_mod_prime(std::array<ulong, 2>& roots,
                                     ulong linear,
                                     ulong constant,
                                     ulong q) noexcept {
    roots = {0, 0};
    if (q < 3 || q % 2 == 0 || linear >= q || constant >= q ||
        constant == 0) {
        return false;
    }
    // The roots (-b +- s) / 2 for s^2 = D = b^2 - 4c, which are distinct
    // and nonzero exactly when D is a nonzero square and c != 0 (Cohen,
    // GTM 138, Proposition 5.1.4); the square root is FLINT's n_sqrtmod
    // (Tonelli-Shanks, Algorithm 1.5.1), checked by squaring, as in
    // monic_quadratic_factor_word (prime_ideal.cpp).
    nmod_t mod;
    nmod_init(&mod, q);
    const ulong discriminant =
            nmod_sub(nmod_mul(linear, linear, mod),
                     nmod_mul(UWORD(4) % q, constant, mod), mod);
    if (discriminant == 0 || n_jacobi_unsigned(discriminant, q) != 1) {
        return false;
    }
    const ulong s = n_sqrtmod(discriminant, q);
    if (s == 0 || nmod_mul(s, s, mod) != discriminant) {
        return false;
    }
    // q is an odd prime below 2^64 - 1, so q + 1 does not overflow.
    const ulong half = (q + 1) / 2;
    const ulong minus_b = nmod_neg(linear, mod);
    const ulong r0 = nmod_mul(nmod_add(minus_b, s, mod), half, mod);
    const ulong r1 = nmod_mul(nmod_sub(minus_b, s, mod), half, mod);

    // The order of nmod_poly_roots (nmod_poly_factor/roots.c,
    // _nmod_poly_push_roots, and _nmod_poly_split_rabin in
    // nmod_poly/find_distinct_nonzero_roots.c).  Below 10 it tries every
    // residue in increasing order.  Otherwise it splits f by
    // a = gcd((x + delta)^((q - 1) / 2) - 1, f), first with delta = 0 and
    // then with delta = n_randint(state, q) for a state fresh from
    // flint_rand_init, until a has degree one; a holds the root r with
    // r + delta a nonzero square, the cofactor f / a the other, and it lists
    // the cofactor's root first.  Each delta separates two distinct roots
    // with probability about 1/2; past kMaxDraws draws the caller takes
    // nmod_poly_roots, which then gives the same roots.
    if (q < 10) {
        roots = {std::min(r0, r1), std::max(r0, r1)};
        return true;
    }
    auto square = [&](ulong r, ulong delta) {
        return n_jacobi_unsigned(nmod_add(r, delta, mod), q) == 1;
    };
    constexpr int kMaxDraws = 128;
    flint_rand_t state;
    flint_rand_init(state);
    ulong delta = 0;
    bool found = false;
    for (int draw = 0; draw <= kMaxDraws; ++draw) {
        if (draw > 0) {
            delta = n_randint(state, q);
        }
        const bool square0 = square(r0, delta);
        if (square0 != square(r1, delta)) {
            roots = square0 ? std::array<ulong, 2>{r1, r0}
                            : std::array<ulong, 2>{r0, r1};
            found = true;
            break;
        }
    }
    flint_rand_clear(state);
    if (!found) {
        roots = {0, 0};
    }
    return found;
}

bool degree_one_roots_mod_q_for_test(std::vector<flint::Fmpz>& roots,
                                     const Order& order,
                                     flint::FmpzConstRef q,
                                     const DiagnosticsContext* diagnostics)
        noexcept {
    return degree_one_roots_mod_q(roots, order, q, diagnostics);
}

}  // namespace detail

namespace {

bool select_saturation_proof_kernel_direct_degree_one(
        slong& local_prime_count,
        bool& certified,
        flint::FmpzMat& kernel,
        flint::FmpzMat& saturation_kernel,
        const OrderUnitGroup& group,
        flint::FmpzConstRef ell,
        flint::FmpzConstRef bound) noexcept {
    SILEX_PROFILE_SCOPE(
            group.diagnostics(), DiagnosticsModule::unit_group,
            "unit_group.select_saturation_proof_kernel_direct_degree_one");
    const Order* order = group.parent();
    const slong rank = group.free_rank();
    slong proof_rank = 0;
    if (!group.is_set() || order == nullptr ||
        !detail::dlog_proof_rank(proof_rank, group, ell) ||
        proof_rank <= 1 || fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    flint::Fmpz order_discriminant;
    if (!order_discriminant_for_saturation(order_discriminant, *order)) {
        return false;
    }

    local_prime_count = 0;
    flint::FmpzMat proof_matrix(proof_rank, 0);
    flint::FmpzMat free_matrix(rank, 0);
    flint::FmpzMat current_kernel(0, proof_rank);
    bool have_kernel = false;

    flint::Fmpz p;
    flint::Fmpz pminus;
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.proof_selector_scan.candidate_prime");
        fmpz_sub_ui(pminus.raw(), p.raw(), 1);
        if (fmpz_divisible(pminus.raw(), ell.raw()) == 0) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.congruence_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }
        if (rational_prime_divides_order_discriminant(
                    flint::FmpzConstRef(order_discriminant),
                    flint::FmpzConstRef(p))) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.discriminant_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        std::vector<flint::Fmpz> roots;
        {
            SILEX_PROFILE_SCOPE(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_degree_one_roots");
            if (!degree_one_roots_mod_q(
                        roots, *order, flint::FmpzConstRef(p),
                        group.diagnostics())) {
                return false;
            }
        }

        if (roots.empty()) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.no_degree_one_prime");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.proof_selector_scan.degree_one_decomposition");
        for (const flint::Fmpz& root : roots) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.degree_one_prime");
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.dlog_column_test");

            flint::FmpzMat proof_column(proof_rank, 1);
            if (!detail::saturation_proof_prime_column_at_degree_one_root(
                        proof_column, group, flint::FmpzConstRef(p),
                        flint::FmpzConstRef(root), ell)) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.direct_dlog_column_unavailable");
                PrimeIdeal fallback_prime(*order);
                if (!fallback_prime.is_defined() ||
                    !detail::PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
                            fallback_prime, *order, flint::FmpzConstRef(p),
                            flint::FmpzConstRef(root))) {
                    return false;
                }
                if (!detail::saturation_proof_prime_column(
                            proof_column, group, fallback_prime, ell)) {
                    SILEX_PROFILE_EVENT(
                            group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.proof_selector_scan.dlog_column_unavailable");
                    continue;
                }
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.generic_dlog_column_available");
            } else {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.direct_dlog_column_available");
            }
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.dlog_column_available");
            ++local_prime_count;

            flint::FmpzMat free_column(rank, 1);
            for (slong row = 0; row < rank; ++row) {
                flint::fmpz_set(
                        flint::fmpz_mat_entry(free_column, row, 0),
                        flint::FmpzConstRef(
                                flint::fmpz_mat_entry(
                                        proof_column, row, 0).raw()));
            }
            if (!append_column(free_matrix, free_column)) {
                return false;
            }
            if (!append_column(proof_matrix, proof_column)) {
                return false;
            }
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_dlog_kernel_nullspace");
            if (!detail::dlog_kernel_from_matrix(current_kernel, proof_matrix,
                                                 ell)) {
                return false;
            }
            have_kernel = true;

            if (flint::fmpz_mat_nrows(current_kernel) == 0) {
                flint::FmpzMat candidate_saturation_kernel(0, rank);
                saturation_kernel = std::move(candidate_saturation_kernel);
                kernel = std::move(current_kernel);
                certified = true;
                return true;
            }
        }

        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    if (!have_kernel) {
        return false;
    }

    kernel = std::move(current_kernel);
    flint::FmpzMat candidate_saturation_kernel(0, rank);
    if (!detail::dlog_kernel_from_matrix(candidate_saturation_kernel,
                                         free_matrix, ell)) {
        return false;
    }
    saturation_kernel = std::move(candidate_saturation_kernel);
    certified = false;
    return true;
}

bool saturate_row_with_cache(
        OrderUnitGroup& out,
        bool& changed,
        const OrderUnitGroup& group,
        flint::FmpzMatConstRef kernel_rows,
        slong row,
        flint::FmpzConstRef ell,
        EmbeddingContext& embeddings,
        slong precision,
        detail::CompactFieldModulusCache* field_modulus_cache) noexcept {
    SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.saturate_row");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!out.is_defined() || !group.is_set() ||
        !same_order_parent(out.parent(), order) ||
        !embedding_has_parent(embeddings, field) || precision <= 0 ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel_rows) ||
        !flint::fmpz_is_prime(ell) || !flint::fmpz_fits_si(ell)) {
        return false;
    }
    bool has_torsion_column = false;
    if (!detail::kernel_rows_torsion_column(has_torsion_column, group,
                                            kernel_rows, ell)) {
        return false;
    }

    if (detail::kernel_row_divisible(kernel_rows, row,
                                     flint::fmpz_mat_ncols(kernel_rows),
                                     ell)) {
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.saturation_row_divisible_skipped");
        changed = false;
        return out.set(group);
    }
    SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.saturation_row_candidate");

    FactoredElement root(*field);
    bool is_power = false;
    SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.saturation_row_root");
    if (detail::kernel_row_root(is_power, root, group, kernel_rows, row, ell,
                                field_modulus_cache) !=
        detail::KernelRowRootResult::decided) {
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.saturation_row_root_unavailable");
        return false;
    }
    if (!is_power) {
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.saturation_row_root_not_power");
        changed = false;
        return out.set(group);
    }
    SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.saturation_row_root_power");
    if (!detail::adjoinable_root_is_order_unit(*order, root)) {
        // The row gives no l-th root in O, so it does not enlarge the group
        // inside O^x.  This is not a proof that the row is not an l-th power
        // in O^x modulo torsion; proof callers only certify from an empty
        // residue-character kernel, which such a row never produces.
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.saturation_row_root_outside_order");
        changed = false;
        return out.set(group);
    }

    flint::FmpzMat relation(1, rank);
    for (slong i = 0; i < rank; ++i) {
        flint::fmpz_set(flint::fmpz_mat_entry(relation, 0, i),
                        flint::fmpz_mat_entry(kernel_rows, row, i));
    }

    std::vector<FactoredElement> generators;
    if (!detail::relation_basis(generators, group, root,
                                flint::FmpzMatConstRef(relation), ell)) {
        return false;
    }

    OrderUnitGroup candidate(*order);
    candidate.set_diagnostics(out.diagnostics());
    if (!candidate.is_defined() ||
        !detail::order_unit_group_set_units_internal(
                candidate, *order,
                FactoredElementSpan(generators.data(), generators.size()),
                embeddings, precision, true)) {
        return false;
    }

    bool reduced = false;
    if (!detail::reduce_stored_relation_units(reduced, candidate, embeddings,
                                           precision)) {
        return false;
    }

    changed = true;
    SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.saturation_row_appended_unit");
    out.swap(candidate);
    return true;
}

}  // namespace

bool OrderUnitGroup::saturate_row(bool& changed,
                                  const OrderUnitGroup& group,
                                  flint::FmpzMatConstRef kernel_rows,
                                  slong row,
                                  flint::FmpzConstRef ell,
                                  EmbeddingContext& embeddings,
                                  slong precision) noexcept {
    // Rows have the residue_dlog_kernel width: free_rank() columns, plus the
    // torsion exponent when ell divides the torsion order.
    slong width = 0;
    if (!detail::dlog_proof_rank(width, group, ell) ||
        flint::fmpz_mat_ncols(kernel_rows) != width) {
        return false;
    }
    return saturate_row_with_cache(*this, changed, group, kernel_rows, row,
                                   ell, embeddings, precision, nullptr);
}

bool OrderUnitGroup::residue_dlog_kernel(
        flint::FmpzMat& out,
        PrimeIdealSpan primes,
        flint::FmpzConstRef ell) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_kernel");
    return detail::residue_dlog_kernel(out, *this, primes, ell);
}

bool OrderUnitGroup::residue_dlog_proof_kernel(
        flint::FmpzMat& out,
        PrimeIdealSpan primes,
        flint::FmpzConstRef ell) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_proof_kernel");
    return detail::residue_dlog_proof_kernel(out, *this, primes, ell);
}

bool OrderUnitGroup::select_saturation_primes(
        PrimeIdealList& out,
        flint::FmpzConstRef ell,
        slong target_len,
        flint::FmpzConstRef bound) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.select_saturation_primes");
    const Order* order = parent();
    if (!is_set() || order == nullptr || !flint::fmpz_is_prime(ell) ||
        target_len < 0 || fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    if (target_len == 0) {
        PrimeIdealList candidate(*order, 0);
        if (!candidate.is_defined()) {
            return false;
        }
        out.swap(candidate);
        return true;
    }

    std::vector<PrimeIdeal> selected;
    selected.reserve(static_cast<std::size_t>(target_len));

    flint::Fmpz order_discriminant;
    if (!order_discriminant_for_saturation(order_discriminant, *order)) {
        return false;
    }

    auto append_usable_prime = [&](const PrimeIdeal& prime) noexcept -> bool {
        if (!detail::saturation_prime_usable(*this, prime, ell)) {
            return true;
        }
        selected.emplace_back(*order);
        return selected.back().set(prime);
    };

    if (!scan_saturation_candidate_primes(
                selected, *order, diagnostics_,
                flint::FmpzConstRef(order_discriminant), ell, target_len,
                bound, false, append_usable_prime)) {
        return false;
    }
    if (copy_selected_primes_to_list(out, *order, selected, target_len)) {
        return true;
    }

    SILEX_PROFILE_EVENT(
            diagnostics_, DiagnosticsModule::unit_group,
            "unit_group.saturation_selector.source_shape_fallback");
    selected.clear();
    if (!scan_saturation_candidate_primes(
                selected, *order, diagnostics_,
                flint::FmpzConstRef(order_discriminant), ell, target_len,
                bound, true, append_usable_prime)) {
        return false;
    }
    return copy_selected_primes_to_list(out, *order, selected, target_len);
}

bool OrderUnitGroup::select_saturation_proof_primes(
        PrimeIdealList& out,
        bool& certified,
        flint::FmpzMat& kernel,
        flint::FmpzConstRef ell,
        flint::FmpzConstRef bound) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.select_saturation_proof_primes");
    const Order* order = parent();
    slong proof_rank = 0;
    if (!is_set() || order == nullptr ||
        !detail::dlog_proof_rank(proof_rank, *this, ell) ||
        fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    if (proof_rank == 0) {
        PrimeIdealList candidate(*order, 0);
        flint::FmpzMat candidate_kernel(0, 0);
        if (!candidate.is_defined()) {
            return false;
        }
        out.swap(candidate);
        kernel = std::move(candidate_kernel);
        certified = true;
        return true;
    }

    RankOneProofSource rank_one_source;
    if (proof_rank == 1 &&
        !set_rank_one_proof_source(rank_one_source, *this, *order)) {
        return false;
    }

    std::vector<PrimeIdeal> selected;
    flint::FmpzMat proof_matrix(proof_rank, 0);
    flint::FmpzMat current_kernel(0, proof_rank);
    bool have_kernel = false;

    flint::Fmpz order_discriminant;
    if (!order_discriminant_for_saturation(order_discriminant, *order)) {
        return false;
    }

    auto process_direct_degree_one_roots =
            [&](flint::FmpzConstRef q,
                bool& used_direct) noexcept -> ProofPrimeScanResult {
        used_direct = false;
        std::vector<flint::Fmpz> roots;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_degree_one_roots");
            if (!degree_one_roots_mod_q(roots, *order, q, diagnostics_)) {
                return ProofPrimeScanResult::keep_scanning;
            }
        }
        used_direct = true;

        if (roots.empty()) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.no_degree_one_prime");
            return ProofPrimeScanResult::keep_scanning;
        }

        std::vector<flint::FmpzMat> local_columns;
        std::vector<PrimeIdeal> local_primes;
        local_columns.reserve(roots.size());
        local_primes.reserve(roots.size());
        for (const flint::Fmpz& root : roots) {
            flint::FmpzMat proof_column(proof_rank, 1);
            if (!detail::saturation_proof_prime_column_at_degree_one_root(
                        proof_column, *this, q, flint::FmpzConstRef(root),
                        ell)) {
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.direct_dlog_column_unavailable");
                used_direct = false;
                return ProofPrimeScanResult::keep_scanning;
            }

            local_primes.emplace_back(*order);
            if (!local_primes.back().is_defined() ||
                !detail::PrimeIdealAccess::set_degree_one_prime_ideal_from_root(
                        local_primes.back(), *order, q,
                        flint::FmpzConstRef(root))) {
                used_direct = false;
                return ProofPrimeScanResult::keep_scanning;
            }
            local_columns.push_back(std::move(proof_column));
        }

        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.proof_selector_scan.degree_one_decomposition");
        for (std::size_t i = 0; i < local_columns.size(); ++i) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.degree_one_prime");
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.dlog_column_test");
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.direct_dlog_column_available");
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.dlog_column_available");

            selected.push_back(std::move(local_primes[i]));
            if (!append_column(proof_matrix, local_columns[i])) {
                return ProofPrimeScanResult::failed;
            }
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_dlog_kernel_nullspace");
            if (!detail::dlog_kernel_from_matrix(current_kernel, proof_matrix,
                                                 ell)) {
                return ProofPrimeScanResult::failed;
            }
            have_kernel = true;

            if (flint::fmpz_mat_nrows(current_kernel) == 0) {
                PrimeIdealList candidate;
                if (!detail::copy_selected_primes(candidate, *order,
                                                  selected)) {
                    return ProofPrimeScanResult::failed;
                }
                out.swap(candidate);
                kernel = std::move(current_kernel);
                certified = true;
                return ProofPrimeScanResult::done;
            }
        }

        return ProofPrimeScanResult::keep_scanning;
    };

    flint::Fmpz p;
    flint::Fmpz pminus;
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.proof_selector_scan.candidate_prime");
        fmpz_sub_ui(pminus.raw(), p.raw(), 1);
        if (fmpz_divisible(pminus.raw(), ell.raw()) == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.congruence_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }
        if (rational_prime_divides_order_discriminant(
                    flint::FmpzConstRef(order_discriminant),
                    flint::FmpzConstRef(p))) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.discriminant_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        if (proof_rank > 1) {
            bool used_direct = false;
            const ProofPrimeScanResult direct_result =
                    process_direct_degree_one_roots(
                            flint::FmpzConstRef(p), used_direct);
            if (used_direct) {
                if (direct_result == ProofPrimeScanResult::done) {
                    return true;
                }
                if (direct_result == ProofPrimeScanResult::failed) {
                    return false;
                }
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
        }

        PrimeIdealList local;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_decompose_degree_one");
            if (!decompose_prime(local, *order, flint::FmpzConstRef(p), 1)) {
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.decompose_failed");
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
        }

        if (local.size() == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.no_degree_one_prime");
        } else {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_selector_scan.degree_one_decomposition");
            for (slong i = 0; i < local.size(); ++i) {
                const PrimeIdeal* prime = local.at(i);
                if (prime == nullptr) {
                    continue;
                }

                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.degree_one_prime");
                selected.emplace_back(*order);
                if (!selected.back().set(*prime)) {
                    return false;
                }

                if (proof_rank == 1) {
                    SILEX_PROFILE_EVENT(
                            diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_selector_scan.rank_one_image_test");
                    bool nonzero = false;
                    if (!rank_one_proof_source_nonzero(
                                nonzero, rank_one_source, *this, *prime,
                                ell)) {
                        SILEX_PROFILE_EVENT(
                                diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.proof_selector_scan.rank_one_image_unavailable");
                        selected.pop_back();
                        continue;
                    }
                    have_kernel = true;
                    if (nonzero) {
                        SILEX_PROFILE_EVENT(
                                diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.proof_selector_scan.rank_one_image_nonzero");
                        PrimeIdealList candidate;
                        flint::FmpzMat candidate_kernel(0, 1);
                        if (!detail::copy_selected_primes(candidate, *order,
                                                          selected)) {
                            return false;
                        }
                        out.swap(candidate);
                        kernel = std::move(candidate_kernel);
                        certified = true;
                        return true;
                    }
                    SILEX_PROFILE_EVENT(
                            diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_selector_scan.rank_one_image_zero");
                    set_rank_one_candidate_kernel(current_kernel);
                    continue;
                }

                flint::FmpzMat proof_column(proof_rank, 1);
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.dlog_column_test");
                if (detail::saturation_proof_prime_column_direct_degree_one(
                            proof_column, *this, *prime, ell)) {
                    SILEX_PROFILE_EVENT(
                            diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_selector_scan.direct_dlog_column_available");
                } else {
                    SILEX_PROFILE_EVENT(
                            diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_selector_scan.direct_dlog_column_unavailable");
                    if (!detail::saturation_proof_prime_column(
                                proof_column, *this, *prime, ell)) {
                        SILEX_PROFILE_EVENT(
                                diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.proof_selector_scan.dlog_column_unavailable");
                        selected.pop_back();
                        continue;
                    }
                    SILEX_PROFILE_EVENT(
                            diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_selector_scan.generic_dlog_column_available");
                }
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_selector_scan.dlog_column_available");
                if (!append_column(proof_matrix, proof_column)) {
                    return false;
                }
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_dlog_kernel_nullspace");
                if (!detail::dlog_kernel_from_matrix(current_kernel,
                                                     proof_matrix, ell)) {
                    return false;
                }
                have_kernel = true;

                if (flint::fmpz_mat_nrows(current_kernel) == 0) {
                    PrimeIdealList candidate;
                    if (!detail::copy_selected_primes(candidate, *order,
                                                      selected)) {
                        return false;
                    }
                    out.swap(candidate);
                    kernel = std::move(current_kernel);
                    certified = true;
                    return true;
                }
            }
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    if (!have_kernel) {
        return false;
    }

    PrimeIdealList candidate;
    if (!detail::copy_selected_primes(candidate, *order, selected)) {
        return false;
    }
    out.swap(candidate);
    kernel = std::move(current_kernel);
    certified = false;
    return true;
}

bool OrderUnitGroup::saturate_local_once(bool& changed,
                                         const OrderUnitGroup& group,
                                         PrimeIdealSpan primes,
                                         flint::FmpzConstRef ell,
                                         EmbeddingContext& embeddings,
                                         slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_local_once");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) || primes.empty() ||
        !flint::fmpz_is_prime(ell) || precision <= 0) {
        return false;
    }
    for (const PrimeIdeal& prime : primes) {
        if (!same_order_parent(prime.parent(), order)) {
            return false;
        }
    }

    flint::FmpzMat kernel(0, 0);
    SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturation_residue_dlog");
    if (!detail::residue_dlog_kernel(kernel, group, primes, ell)) {
        return false;
    }

    return saturate_local_with_kernel_(changed, group,
                                       flint::FmpzMatConstRef(kernel), ell,
                                       embeddings, precision);
}

bool OrderUnitGroup::saturate_local_with_kernel_(
        bool& changed,
        const OrderUnitGroup& group,
        flint::FmpzMatConstRef kernel,
        flint::FmpzConstRef ell,
        EmbeddingContext& embeddings,
        slong precision,
        detail::CompactFieldModulusCache* field_modulus_cache) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_local_with_kernel");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || precision <= 0) {
        return false;
    }
    bool has_torsion_column = false;
    if (!detail::kernel_rows_torsion_column(has_torsion_column, group,
                                            kernel, ell)) {
        return false;
    }

    detail::CompactFieldModulusCache local_field_modulus_cache;
    detail::CompactFieldModulusCache& row_field_modulus_cache =
            field_modulus_cache != nullptr ? *field_modulus_cache
                                           : local_field_modulus_cache;
    bool saw_failure = false;
    for (slong i = 0; i < flint::fmpz_mat_nrows(kernel); ++i) {
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.saturate_local.kernel_row_attempt");
        OrderUnitGroup row_candidate(*order);
        row_candidate.set_diagnostics(diagnostics_);
        bool row_changed = false;
        if (saturate_row_with_cache(row_candidate, row_changed, group,
                                    kernel, i, ell, embeddings, precision,
                                    &row_field_modulus_cache)) {
            if (row_changed) {
                SILEX_PROFILE_EVENT(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_local.kernel_row_changed");
                if (!row_candidate.copy_unit_proof_records_from_(group)) {
                    return false;
                }
                changed = true;
                swap(row_candidate);
                return true;
            }
            SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.saturate_local.kernel_row_unchanged");
        } else {
            SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.saturate_local.kernel_row_failed");
            saw_failure = true;
        }
    }

    if (saw_failure) {
        return false;
    }

    changed = false;
    return set(group);
}

bool OrderUnitGroup::saturate_bounded(bool& changed,
                                      bool& stable,
                                      const OrderUnitGroup& group,
                                      flint::FmpzConstRef ell,
                                      slong aux_target_len,
                                      flint::FmpzConstRef aux_bound,
                                      slong max_passes,
                                      EmbeddingContext& embeddings,
                                      slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_bounded");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || aux_target_len <= 0 ||
        fmpz_cmp_ui(aux_bound.raw(), 2) < 0 || max_passes < 0 ||
        precision <= 0) {
        return false;
    }

    OrderUnitGroup working(*order);
    OrderUnitGroup pass_result(*order);
    working.set_diagnostics(diagnostics_);
    pass_result.set_diagnostics(diagnostics_);
    if (!working.is_defined() || !pass_result.is_defined() ||
        !working.set(group)) {
        return false;
    }
    working.set_diagnostics(diagnostics_);

    // One field-modulus cache for every pass; its entries depend only on
    // the field, the prime and the threshold.
    detail::CompactFieldModulusCache field_modulus_cache;
    bool any_changed = false;
    if (max_passes == 0) {
        if (!set(working)) {
            return false;
        }
        changed = false;
        stable = false;
        return true;
    }

    for (slong pass = 0; pass < max_passes; ++pass) {
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.saturate_bounded.pass");
        PrimeIdealList primes;
        flint::FmpzMat kernel(0, 0);
        if (!detail::select_saturation_primes_with_kernel(
                    primes, kernel, working, ell, aux_target_len, aux_bound)) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.saturate_bounded.prime_selection_failed");
            if (!any_changed) {
                return false;
            }
            if (!set(working)) {
                return false;
            }
            changed = true;
            stable = false;
            return true;
        }

        const PrimeIdeal* first_prime = primes.at(0);
        if (first_prime == nullptr) {
            return false;
        }
        bool pass_changed = false;
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.saturation_precomputed_kernel");
        if (!pass_result.saturate_local_with_kernel_(
                    pass_changed, working, flint::FmpzMatConstRef(kernel),
                    ell, embeddings, precision, &field_modulus_cache)) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.saturate_bounded.local_saturation_failed");
            if (!any_changed) {
                return false;
            }
            if (!set(working)) {
                return false;
            }
            changed = true;
            stable = false;
            return true;
        }

        if (!pass_changed) {
            SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.saturate_bounded.pass_stable");
            if (!set(pass_result)) {
                return false;
            }
            changed = any_changed;
            stable = true;
            return true;
        }

        working.swap(pass_result);
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.saturate_bounded.pass_changed");
        any_changed = true;
    }

    if (!set(working)) {
        return false;
    }
    changed = any_changed;
    stable = false;
    return true;
}

bool detail::OrderUnitGroupAccess::saturate_index_bounded_body(
        OrderUnitGroup& out,
        bool& changed,
        bool& stable,
        const OrderUnitGroup& group,
        EmbeddingContext& embeddings,
        slong aux_target_len,
        flint::FmpzConstRef aux_bound,
        slong max_passes,
        slong precision,
        bool check_torsion) noexcept {
    const DiagnosticsContext* const diagnostics = out.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.saturate_index_bounded");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!out.is_defined() || !group.is_set() ||
        !same_order_parent(out.parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        aux_target_len <= 0 || fmpz_cmp_ui(aux_bound.raw(), 2) < 0 ||
        max_passes < 0 || precision <= 0) {
        return false;
    }

    // The result is labelled proven when the index bound is one, and the
    // passes test ell-th roots up to the torsion generator at ell | w, so
    // both rely on the input torsion.  Fail closed, before any work, unless
    // it is the torsion Silex computes for the order.
    if (check_torsion && !order_unit_torsion_is_computed(group)) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.saturate_index_bounded.torsion_not_computed");
        return false;
    }

    OrderUnitGroup working(*order);
    OrderUnitGroup pass_result(*order);
    working.set_diagnostics(diagnostics);
    pass_result.set_diagnostics(diagnostics);
    if (!working.is_defined() || !pass_result.is_defined() ||
        !working.set(group)) {
        return false;
    }
    working.set_diagnostics(diagnostics);

    flint::Fmpz index_bound;
    if (!group.regulator_index_bound(flint::FmpzRef(index_bound),
                                     precision) ||
        flint::fmpz_sgn(flint::FmpzConstRef(index_bound)) <= 0) {
        return false;
    }

    if (flint::fmpz_is_one(flint::FmpzConstRef(index_bound))) {
        working.certification_ = CertificationMode::proven;
    }

    bool any_changed = false;
    bool all_stable = true;
    flint::Fmpz ell;
    flint::fmpz_set_ui(flint::FmpzRef(ell), 1);
    fmpz_nextprime(ell.raw(), ell.raw(), 1);
    while (fmpz_cmp(ell.raw(), index_bound.raw()) <= 0) {
        bool pass_changed = false;
        bool pass_stable = false;
        if (!pass_result.saturate_bounded(
                    pass_changed, pass_stable, working,
                    flint::FmpzConstRef(ell), aux_target_len, aux_bound,
                    max_passes, embeddings, precision)) {
            if (!any_changed) {
                return false;
            }
            if (!out.set(working)) {
                return false;
            }
            out.try_certify_index_one(precision);
            changed = true;
            stable = false;
            return true;
        }

        if (pass_changed) {
            working.swap(pass_result);
            any_changed = true;
        }

        if (!pass_stable) {
            all_stable = false;
        }

        fmpz_nextprime(ell.raw(), ell.raw(), 1);
    }

    if (!out.set(working)) {
        return false;
    }
    out.try_certify_index_one(precision);
    changed = any_changed;
    stable = all_stable;
    return true;
}

bool detail::OrderUnitGroupAccess::saturate_index_bounded_adaptive_body(
        OrderUnitGroup& out,
        bool& changed,
        bool& stable,
        const OrderUnitGroup& group,
        EmbeddingContext& embeddings,
        slong aux_target_len,
        flint::FmpzConstRef aux_bound_start,
        flint::FmpzConstRef aux_bound_max,
        slong max_passes,
        slong precision,
        bool check_torsion) noexcept {
    const DiagnosticsContext* const diagnostics = out.diagnostics();
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.saturate_index_bounded_adaptive");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!out.is_defined() || !group.is_set() ||
        !same_order_parent(out.parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        aux_target_len <= 0 || fmpz_cmp_ui(aux_bound_start.raw(), 2) < 0 ||
        fmpz_cmp(aux_bound_max.raw(), aux_bound_start.raw()) < 0 ||
        max_passes <= 0 || precision <= 0) {
        return false;
    }

    // As saturate_index_bounded; the torsion is checked once here, and the
    // passes below run the unchecked body.
    if (check_torsion && !order_unit_torsion_is_computed(group)) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.saturate_index_bounded_adaptive."
                "torsion_not_computed");
        return false;
    }

    OrderUnitGroup working(*order);
    OrderUnitGroup pass_result(*order);
    working.set_diagnostics(diagnostics);
    pass_result.set_diagnostics(diagnostics);
    if (!working.is_defined() || !pass_result.is_defined() ||
        !working.set(group)) {
        return false;
    }
    working.set_diagnostics(diagnostics);

    flint::Fmpz aux_bound;
    flint::Fmpz next_bound;
    flint::fmpz_set(flint::FmpzRef(aux_bound), aux_bound_start);
    bool any_changed = false;

    for (;;) {
        bool pass_changed = false;
        bool pass_stable = false;
        // `working` keeps the torsion of `group` or recomputes it when a
        // root is adjoined, so the caller's torsion guarantee still holds.
        if (saturate_index_bounded_body(
                    pass_result, pass_changed, pass_stable, working,
                    embeddings, aux_target_len,
                    flint::FmpzConstRef(aux_bound), max_passes, precision,
                    false)) {
            if (pass_changed || pass_stable) {
                working.swap(pass_result);
                if (pass_changed) {
                    any_changed = true;
                }
            }

            if (pass_stable) {
                if (!out.set(working)) {
                    return false;
                }
                changed = any_changed;
                stable = true;
                return true;
            }
        }

        if (fmpz_cmp(aux_bound.raw(), aux_bound_max.raw()) >= 0) {
            if (!any_changed) {
                return false;
            }
            if (!out.set(working)) {
                return false;
            }
            changed = true;
            stable = false;
            return true;
        }

        fmpz_mul_2exp(next_bound.raw(), aux_bound.raw(), 1);
        if (fmpz_cmp(next_bound.raw(), aux_bound_max.raw()) > 0) {
            flint::fmpz_set(flint::FmpzRef(next_bound), aux_bound_max);
        }
        if (flint::fmpz_equal(flint::FmpzConstRef(next_bound),
                              flint::FmpzConstRef(aux_bound))) {
            if (!any_changed) {
                return false;
            }
            if (!out.set(working)) {
                return false;
            }
            changed = true;
            stable = false;
            return true;
        }
        flint::fmpz_set(flint::FmpzRef(aux_bound),
                        flint::FmpzConstRef(next_bound));
    }
}

bool OrderUnitGroup::saturate_index_bounded(bool& changed,
                                            bool& stable,
                                            const OrderUnitGroup& group,
                                            EmbeddingContext& embeddings,
                                            slong aux_target_len,
                                            flint::FmpzConstRef aux_bound,
                                            slong max_passes,
                                            slong precision) noexcept {
    return detail::OrderUnitGroupAccess::saturate_index_bounded_body(
            *this, changed, stable, group, embeddings, aux_target_len,
            aux_bound, max_passes, precision, true);
}

bool OrderUnitGroup::saturate_index_bounded_adaptive(
        bool& changed,
        bool& stable,
        const OrderUnitGroup& group,
        EmbeddingContext& embeddings,
        slong aux_target_len,
        flint::FmpzConstRef aux_bound_start,
        flint::FmpzConstRef aux_bound_max,
        slong max_passes,
        slong precision) noexcept {
    return detail::OrderUnitGroupAccess::saturate_index_bounded_adaptive_body(
            *this, changed, stable, group, embeddings, aux_target_len,
            aux_bound_start, aux_bound_max, max_passes, precision, true);
}

bool OrderUnitGroup::prove_local_saturated(ProofState& status,
                                           bool& changed,
                                           const OrderUnitGroup& group,
                                           flint::FmpzConstRef ell,
                                           slong aux_target_len,
                                           flint::FmpzConstRef aux_bound,
                                           EmbeddingContext& embeddings,
                                           slong precision) noexcept {
    return prove_local_saturated_(status, changed, group, ell, aux_target_len,
                                  aux_bound, embeddings, precision,
                                  LocalSaturationProofMode::standalone);
}

bool OrderUnitGroup::prove_local_saturated_(
        ProofState& status,
        bool& changed,
        const OrderUnitGroup& group,
        flint::FmpzConstRef ell,
        slong aux_target_len,
        flint::FmpzConstRef aux_bound,
        EmbeddingContext& embeddings,
        slong precision,
        LocalSaturationProofMode mode,
        detail::CompactFieldModulusCache* field_modulus_cache,
        slong field_unit_rank) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.prove_local_saturated");
    // The kernel-row root tests of this proof (detail::kernel_row_root, its
    // only power tests) keep the degree limit of the pure-power lift
    // (detail::kSaturationLiftDegree) for every exponent: from that degree on
    // a row's root test is decided only by the residue disproof, and an
    // undecided row is handled as before.
    detail::SaturationLiftLimit kernel_row_lift_limit(
            detail::SaturationLiftExponents::all);
    detail::CompactFieldModulusCache local_field_modulus_cache;
    detail::CompactFieldModulusCache* const shared_field_modulus_cache =
            field_modulus_cache != nullptr ? field_modulus_cache
                                           : &local_field_modulus_cache;
    const bool use_stable_proof_fallback =
            mode == LocalSaturationProofMode::index_bound_pass;
    const bool check_torsion = mode == LocalSaturationProofMode::standalone;
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || aux_target_len <= 0 ||
        fmpz_cmp_ui(aux_bound.raw(), 2) < 0 || precision <= 0) {
        return false;
    }

    slong rank = field_unit_rank;
    if ((rank < 0 && !unit_rank(rank, *field)) ||
        group.free_rank() != rank) {
        return false;
    }

    // A verified ell-local record at ell | w rests on the torsion column of
    // the saturation kernel, so it relies on the input torsion.  With
    // check_torsion this fails closed, before any work, unless that is the
    // torsion Silex computes for the order.  Only the index_bound_pass mode
    // (prove_index_bound) skips it: it has checked the same input torsion,
    // and its passes keep that torsion or recompute it, so it is not checked
    // again there.
    if (check_torsion && !detail::order_unit_torsion_is_computed(group)) {
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.prove_local_saturated.torsion_not_computed");
        return false;
    }

    OrderUnitGroup working(*order);
    OrderUnitGroup pass_result(*order);
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.prove_local_saturated_setup");
        working.set_diagnostics(diagnostics_);
        pass_result.set_diagnostics(diagnostics_);
        if (!working.is_defined() || !pass_result.is_defined()) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.prove_local_saturated_copy_input");
        if (!working.set(group)) {
            return false;
        }
    }
    working.set_diagnostics(diagnostics_);

    auto publish = [&](const OrderUnitGroup& source,
                       ProofState proof_status,
                       slong local_primes,
                       bool was_changed) noexcept -> bool {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.prove_local_saturated_publish");
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.prove_local_saturated_publish_copy");
            if (!set(source)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.prove_local_saturated_publish_mark");
            if (!mark_unit_proof(ell, proof_status, aux_bound, local_primes,
                                 was_changed)) {
                return false;
            }
        }
        status = proof_status;
        changed = was_changed;
        return true;
    };

    PrimeIdealList primes;
    flint::FmpzMat kernel(0, group.free_rank());
    flint::FmpzMat saturation_kernel(0, group.free_rank());
    bool certified = false;
    slong local_prime_count = 0;
    bool has_prime_records = false;
    bool has_saturation_kernel = false;

    // reference `saturate!` starts the stable scan at stable = 3.5.
    constexpr double kInitialStable = 3.5;
    auto run_stable_fallback = [&]() noexcept -> bool {
        double stable_threshold = kInitialStable;
        bool stable_changed = false;
        bool rescanned = false;
        for (;;) {
            StableRelationSaturationStep step;
            if (!stable_relation_saturation_step(
                        pass_result, step, working, ell, stable_threshold,
                        !rescanned, embeddings, precision,
                        *shared_field_modulus_cache)) {
                return publish(working, ProofState::unavailable,
                               step.local_primes, stable_changed);
            }

            if (step.status == StableRelationSaturationStatus::verified) {
                return publish(working, ProofState::verified,
                               step.local_primes, stable_changed);
            }

            if (step.status ==
                StableRelationSaturationStatus::accepted_root) {
                SILEX_PROFILE_SCOPE(
                        diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.prove_local_saturated_accepted_root");
                if (!pass_result.copy_unit_proof_records_from_(working)) {
                    return publish(working, ProofState::unavailable,
                                   step.local_primes, stable_changed);
                }
                working.swap(pass_result);
                stable_changed = true;
                continue;
            }

            if (step.status ==
                StableRelationSaturationStatus::wasted_non_power) {
                stable_threshold *= 2.0;
                rescanned = true;
                continue;
            }

            return publish(working, ProofState::unavailable,
                           step.local_primes, stable_changed);
        }
    };

    // In the private reference-stable proof route, avoid a bounded pre-scan when
    // it cannot possibly find a proof prime: every usable rational prime q
    // satisfies q = 1 mod ell, hence q > ell.
    if (use_stable_proof_fallback &&
        fmpz_cmp(ell.raw(), aux_bound.raw()) >= 0) {
        return run_stable_fallback();
    }

    SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.proof_prime_selection");
    const bool direct_selection_ok =
            select_saturation_proof_kernel_direct_degree_one(
                    local_prime_count, certified, kernel, saturation_kernel,
                    working, ell, aux_bound);
    if (direct_selection_ok) {
        has_prime_records = false;
        has_saturation_kernel = true;
    } else if (working.select_saturation_proof_primes(
                       primes, certified, kernel, ell, aux_bound)) {
        local_prime_count = primes.size();
        has_prime_records = true;
        has_saturation_kernel = false;
    } else {
        if (use_stable_proof_fallback) {
            return run_stable_fallback();
        }
        return publish(working, ProofState::unavailable, 0, false);
    }

    if (flint::fmpz_mat_nrows(kernel) == 0) {
        return publish(working, ProofState::verified, local_prime_count,
                       false);
    }

    if (!certified && use_stable_proof_fallback) {
        // The pre-scan left a nonempty character kernel.  Before any root
        // attempt on its rows, run the selection that opens the first step
        // of `run_stable_fallback` (the reference `saturate!` scan, which
        // has no pre-scan and tests roots only after its own scan) on the
        // same `working` and threshold.  It reads only `working`, `ell` and
        // the threshold: not `pass_result` and not the field-modulus cache.
        // If its kernel is empty, that step returns `verified` without
        // testing a root, and the group is ell-saturated in O^x, so no
        // pre-scan row has an ell-th root in O that could be adjoined: the
        // path below would end in that same step from the same `working`.
        // Publish its result now.  Otherwise (a nonempty kernel, or a
        // failed scan) nothing was published or changed: continue with the
        // pre-scan root attempts and `run_stable_fallback` as before.
        slong stable_local_primes = 0;
        bool stable_certified = false;
        bool stable_range_exhausted = false;
        flint::FmpzMat stable_kernel(0, working.free_rank());
        if (select_stable_proof_primes(stable_local_primes, stable_certified,
                                       stable_range_exhausted, stable_kernel,
                                       working, ell, kInitialStable) &&
            flint::fmpz_mat_nrows(stable_kernel) == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.proof_prescan_stable_first_verified");
            return publish(working, ProofState::verified, stable_local_primes,
                           false);
        }
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.proof_prescan_stable_first_not_verified");
    }

    const PrimeIdeal* first_prime = has_prime_records ? primes.at(0) : nullptr;
    bool pass_changed = false;
    bool pass_ok = false;
    const flint::FmpzMat& local_kernel =
            has_saturation_kernel ? saturation_kernel : kernel;
    if (flint::fmpz_mat_ncols(local_kernel) == working.free_rank()) {
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_precomputed_kernel");
        pass_ok = pass_result.saturate_local_with_kernel_(
                pass_changed, working, flint::FmpzMatConstRef(local_kernel),
                ell, embeddings, precision, shared_field_modulus_cache);
    } else if (first_prime != nullptr) {
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_precomputed_kernel_unusable");
        // Free-exponent pass over the selected primes; the torsion-twisted
        // rows follow below.
        flint::FmpzMat free_kernel(0, working.free_rank());
        pass_ok = detail::residue_dlog_free_kernel(
                          free_kernel, working,
                          PrimeIdealSpan(first_prime,
                                         static_cast<std::size_t>(
                                                 primes.size())),
                          ell) &&
                  pass_result.saturate_local_with_kernel_(
                          pass_changed, working,
                          flint::FmpzMatConstRef(free_kernel), ell,
                          embeddings, precision, shared_field_modulus_cache);
    }

    const slong free_rank = working.free_rank();
    if (pass_ok && !pass_changed &&
        flint::fmpz_mat_ncols(kernel) == free_rank + 1) {
        // l divides the torsion order, so the proof kernel carries the
        // torsion exponent t.  Reference `compute_candidates_for_saturate`
        // and `saturate!` take candidates prod u_i^{e_i} * zeta^t from this
        // torsion-extended kernel.  The free-width pass above covers t = 0;
        // try the rows with t != 0 mod l, whose l-th roots (a unit that is
        // an l-th power only up to torsion, e.g. -(3 + 2 sqrt2) at l = 2)
        // the free-width kernel cannot see.
        std::vector<slong> twisted_rows;
        for (slong row = 0; row < flint::fmpz_mat_nrows(kernel); ++row) {
            if (fmpz_divisible(flint::fmpz_mat_entry(kernel, row, free_rank)
                                       .raw(),
                               ell.raw()) == 0) {
                twisted_rows.push_back(row);
            }
        }
        if (!twisted_rows.empty()) {
            flint::FmpzMat twisted_kernel(
                    static_cast<slong>(twisted_rows.size()), free_rank + 1);
            for (slong next = 0;
                 next < static_cast<slong>(twisted_rows.size()); ++next) {
                const slong row =
                        twisted_rows[static_cast<std::size_t>(next)];
                for (slong col = 0; col <= free_rank; ++col) {
                    flint::fmpz_set(
                            flint::fmpz_mat_entry(twisted_kernel, next, col),
                            flint::FmpzConstRef(
                                    flint::fmpz_mat_entry(kernel, row, col)
                                            .raw()));
                }
            }
            SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.proof_torsion_twisted_kernel");
            pass_ok = pass_result.saturate_local_with_kernel_(
                    pass_changed, working,
                    flint::FmpzMatConstRef(twisted_kernel), ell, embeddings,
                    precision, shared_field_modulus_cache);
        }
    }

    if (pass_ok && pass_changed) {
        return publish(pass_result, ProofState::unavailable,
                       local_prime_count, true);
    }

    if (!certified && use_stable_proof_fallback) {
        // The pre-scan over q <= aux_bound ended with a nonempty character
        // kernel and adjoined no root (its kernel rows were not l-th powers,
        // or their roots could not be decided).  The bound limits only this
        // pre-scan, not the proof: continue with the reference `saturate!` scan
        // (`run_stable_fallback`), which keeps drawing degree-one primes
        // above q = 1 mod ell until the kernel is empty or stable.
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_prescan_stable_fallback");
        return run_stable_fallback();
    }

    if (pass_ok) {
        return publish(pass_result,
                       certified ? ProofState::verified
                                 : ProofState::unavailable,
                       local_prime_count, false);
    }

    return publish(working, ProofState::unavailable, local_prime_count,
                   false);
}

bool OrderUnitGroup::prove_local_saturated_stable_in_place_(
        ProofState& status,
        bool& changed,
        OrderUnitGroup& group,
        flint::FmpzConstRef ell,
        slong aux_target_len,
        flint::FmpzConstRef aux_bound,
        EmbeddingContext& embeddings,
        slong precision,
        slong stable_proof_record_prefix_len,
        detail::CompactFieldModulusCache* field_modulus_cache,
        slong field_unit_rank) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::unit_group,
            "unit_group.prove_local_saturated_stable_in_place");
    // The kernel-row root tests of this proof (detail::kernel_row_root, its
    // only power tests) keep the degree limit of the pure-power lift
    // (detail::kSaturationLiftDegree) for every exponent: from that degree on
    // a row's root test is decided only by the residue disproof, and an
    // undecided row is handled as before.
    detail::SaturationLiftLimit kernel_row_lift_limit(
            detail::SaturationLiftExponents::all);
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || aux_target_len <= 0 ||
        fmpz_cmp_ui(aux_bound.raw(), 2) < 0 || precision <= 0) {
        return false;
    }

    slong rank = field_unit_rank;
    if ((rank < 0 && !unit_rank(rank, *field)) ||
        group.free_rank() != rank) {
        return false;
    }

    auto publish_in_place = [&](ProofState proof_status,
                                slong local_primes,
                                bool was_changed) noexcept -> bool {
        SILEX_PROFILE_SCOPE(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.prove_local_saturated_in_place_publish");
        if (!group.mark_unit_proof_after_(
                    ell, proof_status, aux_bound, local_primes, was_changed,
                    stable_proof_record_prefix_len)) {
            return false;
        }
        status = proof_status;
        changed = was_changed;
        return true;
    };

    detail::CompactFieldModulusCache local_field_modulus_cache;
    detail::CompactFieldModulusCache& step_field_modulus_cache =
            field_modulus_cache != nullptr ? *field_modulus_cache
                                           : local_field_modulus_cache;
    constexpr double kInitialStable = 3.5;
    double stable_threshold = kInitialStable;
    bool stable_changed = false;
    bool rescanned = false;
    for (;;) {
        OrderUnitGroup pass_result(*order);
        pass_result.set_diagnostics(diagnostics_);
        if (!pass_result.is_defined()) {
            return false;
        }

        StableRelationSaturationStep step;
        if (!stable_relation_saturation_step(
                    pass_result, step, group, ell, stable_threshold,
                    !rescanned, embeddings, precision,
                    step_field_modulus_cache)) {
            return publish_in_place(ProofState::unavailable,
                                    step.local_primes, stable_changed);
        }

        if (step.status == StableRelationSaturationStatus::verified) {
            return publish_in_place(ProofState::verified,
                                    step.local_primes, stable_changed);
        }

        if (step.status == StableRelationSaturationStatus::accepted_root) {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.prove_local_saturated_in_place_accepted_root");
            if (!pass_result.copy_unit_proof_records_from_(group)) {
                return publish_in_place(ProofState::unavailable,
                                        step.local_primes, stable_changed);
            }
            group.swap(pass_result);
            group.set_diagnostics(diagnostics_);
            stable_changed = true;
            continue;
        }

        if (step.status == StableRelationSaturationStatus::wasted_non_power) {
            stable_threshold *= 2.0;
            rescanned = true;
            continue;
        }

        return publish_in_place(ProofState::unavailable,
                                step.local_primes, stable_changed);
    }
}

bool OrderUnitGroup::prove_index_bound(ProofState& status,
                                       bool& changed,
                                       const OrderUnitGroup& group,
                                       slong aux_target_len,
                                       flint::FmpzConstRef aux_bound,
                                       slong max_restarts,
                                       EmbeddingContext& embeddings,
                                       slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.prove_index_bound");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        aux_target_len <= 0 || fmpz_cmp_ui(aux_bound.raw(), 2) < 0 ||
        max_restarts < 0 || precision <= 0) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *field) || group.free_rank() != rank) {
        return false;
    }

    // A verified result is labelled proven, and the ell-local proofs at
    // ell | w use the torsion column, so both rely on the input torsion.
    // Fail closed, before any work, unless it is the torsion Silex computes
    // for the order.  The saturation passes below keep that torsion or
    // recompute it when they adjoin a root, so a proven result carries
    // computed torsion and callers need not check it again.
    if (!detail::order_unit_torsion_is_computed(group)) {
        SILEX_PROFILE_EVENT(
                diagnostics_, DiagnosticsModule::unit_group,
                "unit_group.prove_index_bound.torsion_not_computed");
        return false;
    }

    // `rank` is unit_rank of the field, checked above; every pass below runs
    // on a group of the same order, so they reuse it instead of recomputing
    // the signature, and each still checks its group's free rank against it.
    //
    // One field-modulus cache for every ell-local pass and restart below.
    // Its entries depend only on the field, the prime and the threshold,
    // never on the unit group or order, so roots adjoined between passes
    // leave them valid.
    detail::CompactFieldModulusCache field_modulus_cache;
    OrderUnitGroup working(*order);
    OrderUnitGroup pass_result(*order);
    working.set_diagnostics(diagnostics_);
    pass_result.set_diagnostics(diagnostics_);
    if (!working.is_defined() || !pass_result.is_defined() ||
        !working.set(group)) {
        return false;
    }
    working.set_diagnostics(diagnostics_);

    auto publish = [&](const OrderUnitGroup& source,
                       ProofState proof_status,
                       bool was_changed) noexcept -> bool {
        if (!set(source)) {
            return false;
        }
        status = proof_status;
        changed = was_changed;
        return true;
    };

    bool any_changed = false;
    slong restarts = 0;
    const slong stable_proof_record_prefix_len =
            working.unit_proof_record_count_;
    flint::Fmpz index_bound;
    flint::Fmpz ell;
    flint::fmpz_set_ui(flint::FmpzRef(ell), 2);

    for (;;) {
        if (!working.regulator_index_bound(flint::FmpzRef(index_bound),
                                           precision) ||
            flint::fmpz_sgn(flint::FmpzConstRef(index_bound)) <= 0) {
            return publish(working, ProofState::unavailable, any_changed);
        }
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
        log_unit_proof_index_bound(
                diagnostics_, flint::FmpzConstRef(index_bound),
                restarts, any_changed);
#endif

        if (flint::fmpz_is_one(flint::FmpzConstRef(index_bound))) {
            working.certification_ = CertificationMode::proven;
            SILEX_LOG(diagnostics_, DiagnosticsModule::unit_group,
                      LogLevel::detail, "unit regulator index bound is one");
            return publish(working, ProofState::verified, any_changed);
        }

        bool restart = false;
        while (fmpz_cmp(ell.raw(), index_bound.raw()) <= 0) {
            ProofState pass_status = ProofState::not_checked;
            bool pass_changed = false;
            SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                                "unit_group.prove_local_saturated_pass");
            // reference `_is_definitely_saturated` reads `U.units` without
            // copying the unit context; only `saturate!` mutates after an
            // accepted root.  Keep the same split for the stable fallback.
            const bool use_in_place_stable =
                    fmpz_cmp(ell.raw(), aux_bound.raw()) >= 0;
            bool pass_ok = false;
            if (use_in_place_stable) {
                pass_ok = prove_local_saturated_stable_in_place_(
                        pass_status, pass_changed, working,
                        flint::FmpzConstRef(ell), aux_target_len, aux_bound,
                        embeddings, precision,
                        stable_proof_record_prefix_len, &field_modulus_cache,
                        rank);
            } else {
                pass_ok = pass_result.prove_local_saturated_(
                        pass_status, pass_changed, working,
                        flint::FmpzConstRef(ell), aux_target_len, aux_bound,
                        embeddings, precision,
                        LocalSaturationProofMode::index_bound_pass,
                        &field_modulus_cache, rank);
            }
            if (!pass_ok) {
                return publish(working, ProofState::unavailable, any_changed);
            }

            if (pass_changed) {
                const bool current_prime_verified =
                        pass_status == ProofState::verified;
                if (!use_in_place_stable) {
                    working.swap(pass_result);
                }
                any_changed = true;

                if (restarts >= max_restarts) {
                    return publish(working, ProofState::unavailable, true);
                }

                ++restarts;
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
                const slong carried_records =
                        working.unit_proof_record_count();
                log_unit_proof_restart(
                        diagnostics_, flint::FmpzConstRef(ell),
                        flint::FmpzConstRef(index_bound), restarts,
                        max_restarts, carried_records);
#endif
                if (current_prime_verified) {
                    fmpz_nextprime(ell.raw(), ell.raw(), 1);
                }
                restart = true;
                break;
            }

            if (pass_status != ProofState::verified) {
                return publish(use_in_place_stable ? working : pass_result,
                               ProofState::unavailable, any_changed);
            }

            if (!use_in_place_stable) {
                working.swap(pass_result);
            }
            fmpz_nextprime(ell.raw(), ell.raw(), 1);
        }

        if (restart) {
            continue;
        }

        working.certification_ = CertificationMode::proven;
        return publish(working, ProofState::verified, any_changed);
    }
}

}  // namespace silex
