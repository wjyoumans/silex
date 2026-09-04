#include <silex/relation.hpp>

#include <silex/abelian_group.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal_factorization.hpp>
#include <silex/order_element.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly.h>

#include <optional>
#include <utility>
#include <vector>

namespace silex {
namespace {

bool valid_base(const FactorBase& base) noexcept {
    return base.is_defined() && base.parent() != nullptr &&
           base.parent()->parent() != nullptr;
}

bool sparse_row_equal_fmpz_mat_row(const fmpz_smat::SparseRow& sparse,
                                  flint::FmpzMatConstRef matrix,
                                  slong row,
                                  slong ncols) noexcept {
    if (matrix.raw() == nullptr || row < 0 ||
        row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != ncols ||
        !sparse.fits_columns(ncols)) {
        return false;
    }

    const auto columns = sparse.columns();
    const flint::FmpzVecConstRef values = sparse.values();
    slong sparse_index = 0;
    for (slong col = 0; col < ncols; ++col) {
        const fmpz* dense_entry = fmpz_mat_entry(matrix.raw(), row, col);
        if (sparse_index < sparse.length() &&
            columns[static_cast<std::size_t>(sparse_index)] == col) {
            if (fmpz_equal(dense_entry, values.data() + sparse_index) == 0) {
                return false;
            }
            ++sparse_index;
        } else if (fmpz_is_zero(dense_entry) == 0) {
            return false;
        }
    }

    return sparse_index == sparse.length();
}

enum class ElementFactorOverBaseStatus {
    unsupported,
    smooth,
    nonsmooth,
};

struct ResidueHit {
    slong index = -1;
    slong residue_degree = 0;
};

struct NormPrimeBlockFactor {
    slong block_index = -1;
    slong exponent = 0;
};

struct RelationRowWriteGuard {
    explicit RelationRowWriteGuard(flint::FmpzMatRef row) noexcept
        : row_(row.raw()) {
        fmpz_mat_zero(row_);
    }

    RelationRowWriteGuard(const RelationRowWriteGuard&) = delete;
    RelationRowWriteGuard& operator=(const RelationRowWriteGuard&) = delete;

    ~RelationRowWriteGuard() noexcept {
        if (!committed_ && row_ != nullptr) {
            fmpz_mat_zero(row_);
        }
    }

    fmpz_mat_struct* row() noexcept {
        return row_;
    }

    void commit() noexcept {
        committed_ = true;
    }

private:
    fmpz_mat_struct* row_ = nullptr;
    bool committed_ = false;
};

struct ResidueReductionBlockScratch {
    std::optional<flint::FmpzModCtx> ctx;
    std::optional<flint::FmpzModPoly> reduced;
    std::optional<flint::FmpzModPoly> modulus;
    std::optional<flint::FmpzModPoly> remainder;
    bool initialized = false;

    bool init(const fmpz_poly_struct* input,
              flint::FmpzConstRef p) noexcept {
        ctx.emplace(p.raw());
        if (ctx->raw() == nullptr) {
            return false;
        }
        reduced.emplace(*ctx);
        modulus.emplace(*ctx);
        remainder.emplace(*ctx);
        if (!reduced->is_initialized() || !modulus->is_initialized() ||
            !remainder->is_initialized()) {
            return false;
        }
        fmpz_mod_poly_set_fmpz_poly(reduced->raw(), input,
                                    ctx->raw());
        initialized = true;
        return true;
    }
};

struct LinearResidueEvaluationScratch {
    flint::Fmpz value;
};

#define SILEX_RELATION_PROFILE_EVENT(diagnostics, label) \
    SILEX_PROFILE_EVENT((diagnostics), DiagnosticsModule::relation, (label))

#define SILEX_RELATION_PROFILE_RETURN(diagnostics, label, value) \
    do {                                                        \
        SILEX_RELATION_PROFILE_EVENT((diagnostics), (label));   \
        return (value);                                         \
    } while (false)

void integral_coordinates_to_polynomial(
        flint::FmpzPoly& out,
        flint::FmpzMatConstRef coordinates) noexcept {
    fmpz_poly_zero(out.raw());
    const slong cols = flint::fmpz_mat_ncols(coordinates);
    for (slong j = 0; j < cols; ++j) {
        fmpz_poly_set_coeff_fmpz(
                out.raw(), j,
                flint::fmpz_mat_entry(coordinates, 0, j).raw());
    }
}

bool reduce_polynomial_mod_residue(
        flint::FmpzPoly& out,
        const fmpz_poly_struct* input,
        const flint::FmpzPoly& residue_polynomial,
        flint::FmpzConstRef p) noexcept {
    if (input == nullptr) {
        return false;
    }
    const slong residue_degree = fmpz_poly_degree(residue_polynomial.raw());
    if (residue_degree <= 0) {
        return false;
    }

    if (residue_degree == 1) {
        flint::Fmpz constant;
        flint::Fmpz leading;
        flint::Fmpz inverse_leading;
        flint::Fmpz root;
        flint::Fmpz value;
        fmpz_poly_get_coeff_fmpz(constant.raw(),
                                 residue_polynomial.raw(), 0);
        fmpz_poly_get_coeff_fmpz(leading.raw(),
                                 residue_polynomial.raw(), 1);
        fmpz_mod(leading.raw(), leading.raw(), p.raw());
        if (fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) != 0) {
            fmpz_neg(root.raw(), constant.raw());
            fmpz_mod(root.raw(), root.raw(), p.raw());
            fmpz_mul(root.raw(), root.raw(), inverse_leading.raw());
            fmpz_mod(root.raw(), root.raw(), p.raw());

            fmpz_zero(value.raw());
            const slong input_degree = fmpz_poly_degree(input);
            for (slong i = input_degree; i >= 0; --i) {
                fmpz_mul(value.raw(), value.raw(), root.raw());
                fmpz_poly_get_coeff_fmpz(
                        constant.raw(), input, i);
                fmpz_add(value.raw(), value.raw(), constant.raw());
                fmpz_mod(value.raw(), value.raw(), p.raw());
            }
            fmpz_poly_set_fmpz(out.raw(), value.raw());
            return true;
        }
    }

    flint::FmpzModCtx ctx(p.raw());
    if (ctx.raw() == nullptr) {
        return false;
    }
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPoly modulus(ctx);
    flint::FmpzModPoly remainder(ctx);
    if (!reduced.is_initialized() || !modulus.is_initialized() ||
        !remainder.is_initialized()) {
        return false;
    }

    fmpz_mod_poly_set_fmpz_poly(reduced.raw(), input, ctx.raw());
    fmpz_mod_poly_set_fmpz_poly(modulus.raw(),
                                residue_polynomial.raw(), ctx.raw());
    if (fmpz_mod_poly_is_zero(modulus.raw(), ctx.raw()) != 0) {
        return false;
    }
    fmpz_mod_poly_rem(remainder.raw(), reduced.raw(), modulus.raw(),
                      ctx.raw());
    fmpz_mod_poly_get_fmpz_poly(out.raw(), remainder.raw(), ctx.raw());
    return true;
}

bool reduce_polynomial_mod_residue_with_block_scratch(
        flint::FmpzPoly& out,
        ResidueReductionBlockScratch& scratch,
        const fmpz_poly_struct* input,
        const flint::FmpzPoly& residue_polynomial,
        flint::FmpzConstRef p) noexcept {
    if (input == nullptr) {
        return false;
    }
    const slong residue_degree = fmpz_poly_degree(residue_polynomial.raw());
    if (residue_degree <= 0) {
        return false;
    }
    if (residue_degree == 1 &&
        reduce_polynomial_mod_residue(out, input, residue_polynomial, p)) {
        return true;
    }
    if (!scratch.initialized && !scratch.init(input, p)) {
        return false;
    }
    if (!scratch.ctx.has_value() || !scratch.reduced.has_value() ||
        !scratch.modulus.has_value() || !scratch.remainder.has_value()) {
        return false;
    }

    fmpz_mod_poly_set_fmpz_poly(scratch.modulus->raw(),
                                residue_polynomial.raw(),
                                scratch.ctx->raw());
    if (fmpz_mod_poly_is_zero(scratch.modulus->raw(),
                              scratch.ctx->raw()) != 0) {
        return false;
    }
    fmpz_mod_poly_rem(scratch.remainder->raw(), scratch.reduced->raw(),
                      scratch.modulus->raw(), scratch.ctx->raw());
    fmpz_mod_poly_get_fmpz_poly(out.raw(), scratch.remainder->raw(),
                                scratch.ctx->raw());
    return true;
}

bool coordinate_polynomial_vanishes_at_linear_root_mod_prime(
        bool& out,
        LinearResidueEvaluationScratch& scratch,
        flint::FmpzMatConstRef coordinates,
        flint::FmpzConstRef root,
        flint::FmpzConstRef p) noexcept {
    out = false;
    if (flint::fmpz_mat_nrows(coordinates) != 1) {
        return false;
    }

    const slong cols = flint::fmpz_mat_ncols(coordinates);
    if (cols <= 0) {
        out = true;
        return true;
    }

    fmpz_zero(scratch.value.raw());
    for (slong i = cols; i-- > 0;) {
        fmpz_mul(scratch.value.raw(), scratch.value.raw(), root.raw());
        fmpz_add(scratch.value.raw(), scratch.value.raw(),
                 flint::fmpz_mat_entry(coordinates, 0, i).raw());
        fmpz_mod(scratch.value.raw(), scratch.value.raw(), p.raw());
    }
    out = fmpz_is_zero(scratch.value.raw()) != 0;
    return true;
}

bool set_element_from_fmpz_poly(Element& out,
                                const flint::FmpzPoly& polynomial) noexcept {
    flint::FmpqPoly rational_polynomial;
    fmpq_poly_set_fmpz_poly(rational_polynomial.raw(), polynomial.raw());
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(rational_polynomial));
}

bool set_element_from_integral_coordinates(Element& out,
                                           const Order& order,
                                           flint::FmpzMatConstRef coordinates)
        noexcept {
    if (out.parent() == nullptr || order.parent() == nullptr ||
        !out.parent()->has_same_data(*order.parent()) ||
        flint::fmpz_mat_nrows(coordinates) != 1 ||
        flint::fmpz_mat_ncols(coordinates) != order.degree()) {
        return false;
    }
    if (order.is_equation_order()) {
        flint::FmpzPoly polynomial;
        integral_coordinates_to_polynomial(polynomial, coordinates);
        return set_element_from_fmpz_poly(out, polynomial);
    }

    OrderElement order_element(order);
    return order_element.is_defined() &&
           order_element.set_coordinates(coordinates) &&
           order_element.get_element(out);
}

bool factor_norm_over_base_blocks(
        std::vector<NormPrimeBlockFactor>& out,
        bool& smooth_over_base,
        const FactorBase& base,
        flint::FmpzConstRef abs_norm,
        const DiagnosticsContext* diagnostics) noexcept {
    out.clear();
    smooth_over_base = true;
    if (!valid_base(base)) {
        return false;
    }
    if (flint::fmpz_cmp_ui(abs_norm, 1) <= 0) {
        return true;
    }

    flint::FmpzFactor factorization;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::relation,
                "relation.factor_over_base_integral_factor_norm");
        flint::fmpz_factor(flint::FmpzFactorRef(factorization), abs_norm);
    }
    const slong factor_count =
            flint::fmpz_factor_num(flint::FmpzFactorConstRef(factorization));
    out.reserve(static_cast<std::size_t>(factor_count));
    flint::Fmpz rational_prime;
    for (slong i = 0; i < factor_count; ++i) {
        flint::fmpz_factor_get_fmpz(
                flint::FmpzRef(rational_prime),
                flint::FmpzFactorConstRef(factorization), i);
        slong block_index = -1;
        if (!base.rational_prime_block_index_for_prime(
                    block_index, flint::FmpzConstRef(rational_prime))) {
            smooth_over_base = false;
            out.clear();
            return true;
        }
        const ulong exponent = flint::fmpz_factor_exp(
                flint::FmpzFactorConstRef(factorization), i);
        if (exponent > static_cast<ulong>(WORD_MAX)) {
            return false;
        }
        out.push_back(NormPrimeBlockFactor{
                block_index, static_cast<slong>(exponent)});
    }
    return true;
}

bool factor_integral_data_over_base_by_residue_screen(
        ElementFactorOverBaseStatus& status,
        flint::FmpzMatRef exponents,
        const Element* alpha,
        const FactorBase& base,
        const fmpq* known_norm,
        const fmpz_mat_struct* known_integral_coordinates,
        const fmpz_poly_struct* known_integral_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    // reference fb_int_doit first handles integral elements by residue-polynomial
    // membership.  If one residue hit per factor does not account for the
    // whole p-adic norm valuation, it then computes exact valuations only for
    // the residue factors that divide the element.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                        "relation.factor_over_base_integral_element");
    status = ElementFactorOverBaseStatus::unsupported;

    const Order* parent = base.parent();
    if (!valid_base(base) || parent == nullptr ||
        parent->parent() == nullptr ||
        (alpha != nullptr &&
         (!alpha->has_parent(*parent->parent()) || alpha->equal_si(0))) ||
        (alpha == nullptr &&
         (known_integral_coordinates == nullptr || known_norm == nullptr ||
          ::fmpq_is_zero(known_norm) != 0)) ||
        flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length()) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_integral.error_invalid_input",
                false);
    }
    // Current factor-base residue polynomials are in defining-generator
    // coordinates.  Non-equation maximal orders stay on the exact fallback.
    if (!parent->is_maximal()) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_integral.fallback_nonmaximal_order",
                true);
    }
    if (!parent->is_equation_order()) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_integral.fallback_non_equation_order",
                true);
    }

    OrderElement integral_alpha(*parent);
    bool integral_alpha_ready = false;
    auto ensure_integral_alpha = [&]() noexcept -> bool {
        if (integral_alpha_ready) {
            return true;
        }
        if (!integral_alpha.is_defined()) {
            return false;
        }
        if (known_integral_coordinates != nullptr) {
            if (!integral_alpha.set_coordinates(
                        flint::FmpzMatConstRef(
                                known_integral_coordinates))) {
                return false;
            }
        } else if (alpha == nullptr || !integral_alpha.set_element(*alpha)) {
            return false;
        }
        integral_alpha_ready = true;
        return true;
    };
    flint::Fmpq norm;
    const fmpq* norm_value = known_norm;
    flint::Fmpz remaining_norm;
    flint::FmpzMat coordinates(1, parent->degree());
    const fmpz_mat_struct* coordinate_rows = known_integral_coordinates;
    flint::FmpzPoly input_polynomial_storage;
    const fmpz_poly_struct* input_polynomial = known_integral_polynomial;
    flint::FmpzPoly reduced;
    std::vector<NormPrimeBlockFactor> norm_block_factors;
    bool use_factored_norm_blocks = false;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::relation,
                "relation.factor_over_base_integral_setup");
        if (known_integral_coordinates != nullptr &&
            (flint::fmpz_mat_nrows(
                     flint::FmpzMatConstRef(known_integral_coordinates)) != 1 ||
             flint::fmpz_mat_ncols(
                     flint::FmpzMatConstRef(known_integral_coordinates)) !=
                     parent->degree())) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_integral.error_coordinate_shape",
                    false);
        }
        if (coordinate_rows == nullptr && !ensure_integral_alpha()) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_integral.fallback_nonintegral_element",
                    true);
        }
        if (norm_value == nullptr) {
            if (alpha == nullptr ||
                !alpha->norm(flint::FmpqRef(norm))) {
                SILEX_RELATION_PROFILE_RETURN(
                        diagnostics,
                        "relation.factor_over_base_integral.error_norm",
                        false);
            }
            norm_value = norm.raw();
        }
        if (fmpz_is_one(fmpq_denref(norm_value)) == 0) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_integral.fallback_nonintegral_norm",
                    true);
        }
        fmpz_abs(remaining_norm.raw(), fmpq_numref(norm_value));
        if (coordinate_rows == nullptr &&
            !integral_alpha.get_coordinates(flint::FmpzMatRef(coordinates))) {
            return false;
        }
        if (coordinate_rows == nullptr) {
            coordinate_rows = coordinates.raw();
        }
        if (input_polynomial == nullptr) {
            integral_coordinates_to_polynomial(
                    input_polynomial_storage,
                    flint::FmpzMatConstRef(coordinate_rows));
            input_polynomial = input_polynomial_storage.raw();
        }

        if (known_integral_coordinates != nullptr && known_norm != nullptr) {
            bool smooth_over_base = true;
            if (!factor_norm_over_base_blocks(
                        norm_block_factors, smooth_over_base, base,
                        flint::FmpzConstRef(remaining_norm), diagnostics)) {
                return false;
            }
            if (!smooth_over_base) {
                status = ElementFactorOverBaseStatus::nonsmooth;
                flint::fmpz_mat_zero(exponents);
                return true;
            }
            use_factored_norm_blocks = true;
        }
    }

    RelationRowWriteGuard candidate(exponents);

    Ideal principal_ideal;
    bool has_principal_ideal = false;
    flint::Fmpz p;
    std::vector<ResidueHit> hits;
    LinearResidueEvaluationScratch linear_scratch;
    const slong norm_block_count = use_factored_norm_blocks
            ? static_cast<slong>(norm_block_factors.size())
            : base.rational_prime_block_count();
    for (slong norm_block_index = 0;
         norm_block_index < norm_block_count &&
         (use_factored_norm_blocks ||
          fmpz_is_one(remaining_norm.raw()) == 0);
         ++norm_block_index) {
        const slong block_index = use_factored_norm_blocks
                ? norm_block_factors[static_cast<std::size_t>(
                          norm_block_index)]
                          .block_index
                : norm_block_index;
        slong length = 0;
        if (!base.rational_prime_block_data(flint::FmpzRef(p), length,
                                            block_index)) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_integral.error_prime_block_data",
                    false);
        }

        slong remaining_vp = use_factored_norm_blocks
                ? norm_block_factors[static_cast<std::size_t>(
                      norm_block_index)]
                          .exponent
                : 0;
        if (!use_factored_norm_blocks) {
            while (fmpz_divisible(remaining_norm.raw(), p.raw()) != 0) {
                fmpz_divexact(remaining_norm.raw(), remaining_norm.raw(),
                              p.raw());
                ++remaining_vp;
            }
        }
        if (remaining_vp == 0) {
            continue;
        }

        const slong block_vp = remaining_vp;
        hits.clear();
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::relation,
                    "relation.factor_over_base_integral_residue_scan");
            ResidueReductionBlockScratch reduction_scratch;
            for (slong offset = 0; offset < length; ++offset) {
                slong index = -1;
                if (!base.rational_prime_block_index(index, block_index,
                                                     offset)) {
                    SILEX_RELATION_PROFILE_RETURN(
                            diagnostics,
                            "relation.factor_over_base_integral.error_prime_block_index",
                            false);
                }
                const PrimeIdeal* prime = base.prime_at(index);
                if (prime == nullptr) {
                    SILEX_RELATION_PROFILE_RETURN(
                            diagnostics,
                            "relation.factor_over_base_integral.error_missing_prime",
                            false);
                }
                const flint::FmpzPoly* residue_polynomial =
                        detail::residue_polynomial_ptr(*prime);
                if (residue_polynomial == nullptr) {
                    SILEX_RELATION_PROFILE_RETURN(
                            diagnostics,
                            "relation.factor_over_base_integral.fallback_missing_residue_polynomial",
                            true);
                }
                SILEX_RELATION_PROFILE_EVENT(
                        diagnostics,
                        fmpz_poly_degree(residue_polynomial->raw()) == 1
                                ? "relation.factor_over_base_integral.residue_linear"
                                : "relation.factor_over_base_integral.residue_generic");
                const flint::Fmpz* linear_root =
                        detail::linear_residue_root_ptr(*prime);
                if (linear_root != nullptr) {
                    bool residue_zero = false;
                    if (!coordinate_polynomial_vanishes_at_linear_root_mod_prime(
                                residue_zero, linear_scratch,
                                flint::FmpzMatConstRef(coordinate_rows),
                                flint::FmpzConstRef(*linear_root),
                                flint::FmpzConstRef(p))) {
                        SILEX_RELATION_PROFILE_RETURN(
                                diagnostics,
                                "relation.factor_over_base_integral.fallback_linear_residue_evaluation",
                                true);
                    }
                    if (residue_zero) {
                        hits.push_back(ResidueHit{index,
                                                  prime->residue_degree()});
                    }
                    continue;
                }
                if (!reduce_polynomial_mod_residue_with_block_scratch(
                            reduced, reduction_scratch, input_polynomial,
                            *residue_polynomial, flint::FmpzConstRef(p))) {
                    SILEX_RELATION_PROFILE_RETURN(
                            diagnostics,
                            "relation.factor_over_base_integral.fallback_residue_reduction",
                            true);
                }
                if (fmpz_poly_is_zero(reduced.raw()) != 0) {
                    hits.push_back(ResidueHit{index,
                                              prime->residue_degree()});
                }
            }
        }

        slong residue_hit_vp = 0;
        for (const ResidueHit& hit : hits) {
            if (hit.index < 0 || hit.residue_degree <= 0) {
                SILEX_RELATION_PROFILE_RETURN(
                        diagnostics,
                        "relation.factor_over_base_integral.error_bad_residue_hit",
                        false);
            }
            residue_hit_vp += hit.residue_degree;
        }

        if (residue_hit_vp == block_vp) {
            SILEX_RELATION_PROFILE_EVENT(
                    diagnostics,
                    "relation.factor_over_base_integral.path_residue_degree_exact");
            for (const ResidueHit& hit : hits) {
                fmpz_set_si(
                        fmpz_mat_entry(candidate.row(), 0, hit.index),
                        1);
            }
            continue;
        }

        if (hits.empty()) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_integral.fallback_no_residue_hits",
                    true);
        }

        slong exact_remaining_vp = block_vp;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::relation,
                    "relation.factor_over_base_integral_exact_block");
            for (const ResidueHit& hit : hits) {
                if (exact_remaining_vp == 0) {
                    break;
                }
                const PrimeIdeal* prime = base.prime_at(hit.index);
                slong valuation = -1;
                if (prime == nullptr) {
                    SILEX_RELATION_PROFILE_RETURN(
                            diagnostics,
                            "relation.factor_over_base_integral.error_missing_hit_prime",
                            false);
                }
                const bool hit_is_ramified = prime->ramification_index() != 1;
                const bool hit_is_unique_over_p =
                        prime->ramification_index() > 0 &&
                        prime->residue_degree() > 0 &&
                        parent->degree() ==
                                prime->ramification_index() *
                                        prime->residue_degree();
                if (hit_is_ramified) {
                    SILEX_RELATION_PROFILE_EVENT(
                            diagnostics,
                            "relation.factor_over_base_integral.exact_hit_ramified");
                }
                if (!hit_is_unique_over_p) {
                    SILEX_RELATION_PROFILE_EVENT(
                            diagnostics,
                            "relation.factor_over_base_integral.exact_hit_nonunique");
                }
                bool has_valuation = false;
                {
                    SILEX_PROFILE_SCOPE(
                            diagnostics, DiagnosticsModule::relation,
                            "relation.factor_over_base_integral_element_valuation");
                    if (!ensure_integral_alpha()) {
                        SILEX_RELATION_PROFILE_RETURN(
                                diagnostics,
                                "relation.factor_over_base_integral.error_integral_alpha",
                                false);
                    }
                    has_valuation =
                            detail::prime_ideal_valuation_with_norm_vp(
                                    valuation, *prime, integral_alpha,
                                    block_vp, diagnostics);
                }

                if (!has_valuation) {
                    SILEX_RELATION_PROFILE_EVENT(
                            diagnostics,
                            "relation.factor_over_base_integral.path_principal_ideal_valuation_fallback");
                    if (hit_is_ramified) {
                        SILEX_RELATION_PROFILE_EVENT(
                                diagnostics,
                                "relation.factor_over_base_integral.path_principal_ideal_valuation_fallback_ramified");
                    }
                    if (!hit_is_unique_over_p) {
                        SILEX_RELATION_PROFILE_EVENT(
                                diagnostics,
                                "relation.factor_over_base_integral.path_principal_ideal_valuation_fallback_nonunique");
                    }
                    if (!has_principal_ideal) {
                        if (!principal_ideal.is_defined() &&
                            !principal_ideal.define(*parent)) {
                            SILEX_RELATION_PROFILE_RETURN(
                                    diagnostics,
                                    "relation.factor_over_base_integral.error_principal_ideal_define",
                                    false);
                        }
                        SILEX_PROFILE_SCOPE(
                                diagnostics, DiagnosticsModule::relation,
                                "relation.factor_over_base_integral_principal_ideal");
                        if (!ensure_integral_alpha()) {
                            SILEX_RELATION_PROFILE_RETURN(
                                    diagnostics,
                                    "relation.factor_over_base_integral.error_integral_alpha",
                                    false);
                        }
                        if (!principal_ideal.set_principal(integral_alpha)) {
                            SILEX_RELATION_PROFILE_RETURN(
                                    diagnostics,
                                    "relation.factor_over_base_integral.error_principal_ideal_set",
                                    false);
                        }
                        has_principal_ideal = true;
                    }

                    SILEX_PROFILE_SCOPE(
                            diagnostics, DiagnosticsModule::relation,
                            "relation.factor_over_base_integral_valuation");
                    if (!prime->valuation(valuation, principal_ideal)) {
                        SILEX_RELATION_PROFILE_RETURN(
                                diagnostics,
                                "relation.factor_over_base_integral.error_principal_ideal_valuation",
                                false);
                    }
                } else {
                    SILEX_RELATION_PROFILE_EVENT(
                            diagnostics,
                            "relation.factor_over_base_integral.path_direct_valuation");
                    if (hit_is_ramified) {
                        SILEX_RELATION_PROFILE_EVENT(
                                diagnostics,
                                "relation.factor_over_base_integral.path_direct_valuation_ramified");
                    }
                    if (!hit_is_unique_over_p) {
                        SILEX_RELATION_PROFILE_EVENT(
                                diagnostics,
                                "relation.factor_over_base_integral.path_direct_valuation_nonunique");
                    }
                }
                if (valuation <= 0) {
                    continue;
                }
                fmpz_set_si(
                        fmpz_mat_entry(candidate.row(), 0, hit.index),
                        valuation);
                exact_remaining_vp -= valuation * hit.residue_degree;
            }
        }

        if (exact_remaining_vp != 0) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_integral.fallback_exact_vp_leftover",
                    true);
        }
    }

    if (!use_factored_norm_blocks &&
        fmpz_is_one(remaining_norm.raw()) == 0) {
        SILEX_RELATION_PROFILE_EVENT(
                diagnostics,
                "relation.factor_over_base_integral.nonsmooth_norm_leftover");
        status = ElementFactorOverBaseStatus::nonsmooth;
        flint::fmpz_mat_zero(exponents);
        return true;
    }

    SILEX_RELATION_PROFILE_EVENT(
            diagnostics, "relation.factor_over_base_integral.smooth");
    status = ElementFactorOverBaseStatus::smooth;
    candidate.commit();
    return true;
}

bool factor_integral_element_over_base_by_residue_screen(
        ElementFactorOverBaseStatus& status,
        flint::FmpzMatRef exponents,
        const Element& alpha,
        const FactorBase& base,
        const fmpq* known_norm,
        const fmpz_mat_struct* known_integral_coordinates,
        const fmpz_poly_struct* known_integral_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    return factor_integral_data_over_base_by_residue_screen(
            status, exponents, &alpha, base, known_norm,
            known_integral_coordinates, known_integral_polynomial,
            diagnostics);
}

bool factor_integral_coordinates_over_base_by_residue_screen(
        ElementFactorOverBaseStatus& status,
        flint::FmpzMatRef exponents,
        const FactorBase& base,
        flint::FmpqConstRef norm,
        flint::FmpzMatConstRef integral_coordinates,
        const flint::FmpzPoly* integral_coordinate_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    return factor_integral_data_over_base_by_residue_screen(
            status, exponents, nullptr, base, norm.raw(),
            integral_coordinates.raw(),
            integral_coordinate_polynomial == nullptr
                    ? nullptr
                    : integral_coordinate_polynomial->raw(),
            diagnostics);
}

bool subtract_rational_denominator_row(
        ElementFactorOverBaseStatus& status,
        flint::FmpzMatRef exponents,
        flint::FmpzConstRef denominator,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics) noexcept {
    status = ElementFactorOverBaseStatus::unsupported;
    if (!valid_base(base) || flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length() ||
        fmpz_sgn(denominator.raw()) <= 0) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.error_denominator_input",
                false);
    }

    const Order* parent = base.parent();
    if (parent == nullptr) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.error_missing_parent",
                false);
    }

    flint::Fmpz remaining;
    fmpz_set(remaining.raw(), denominator.raw());
    flint::Fmpz p;
    for (slong block_index = 0;
         block_index < base.rational_prime_block_count() &&
         fmpz_is_one(remaining.raw()) == 0;
         ++block_index) {
        slong length = 0;
        if (!base.rational_prime_block_data(flint::FmpzRef(p), length,
                                            block_index)) {
            SILEX_RELATION_PROFILE_RETURN(
                    diagnostics,
                    "relation.factor_over_base_nonintegral.error_denominator_block_data",
                    false);
        }

        slong denominator_vp = 0;
        while (fmpz_divisible(remaining.raw(), p.raw()) != 0) {
            fmpz_divexact(remaining.raw(), remaining.raw(), p.raw());
            ++denominator_vp;
        }
        if (denominator_vp == 0) {
            continue;
        }

        slong block_degree = 0;
        for (slong offset = 0; offset < length; ++offset) {
            slong index = -1;
            if (!base.rational_prime_block_index(index, block_index,
                                                 offset)) {
                SILEX_RELATION_PROFILE_RETURN(
                        diagnostics,
                        "relation.factor_over_base_nonintegral.error_denominator_block_index",
                        false);
            }
            const PrimeIdeal* prime = base.prime_at(index);
            if (prime == nullptr || prime->ramification_index() <= 0 ||
                prime->residue_degree() <= 0 ||
                prime->ramification_index() >
                        WORD_MAX / denominator_vp ||
                prime->ramification_index() >
                        WORD_MAX / prime->residue_degree()) {
                SILEX_RELATION_PROFILE_RETURN(
                        diagnostics,
                        "relation.factor_over_base_nonintegral.error_denominator_prime",
                        false);
            }

            const slong prime_degree =
                    prime->ramification_index() *
                    prime->residue_degree();
            if (block_degree > WORD_MAX - prime_degree) {
                SILEX_RELATION_PROFILE_RETURN(
                        diagnostics,
                        "relation.factor_over_base_nonintegral.error_denominator_degree_overflow",
                        false);
            }
            block_degree += prime_degree;

            const slong exponent =
                    prime->ramification_index() * denominator_vp;
            fmpz_sub_ui(flint::fmpz_mat_entry(exponents, 0, index).raw(),
                        flint::fmpz_mat_entry(exponents, 0, index).raw(),
                        static_cast<ulong>(exponent));
        }

        if (block_degree != parent->degree()) {
            SILEX_RELATION_PROFILE_EVENT(
                    diagnostics,
                    "relation.factor_over_base_nonintegral.nonsmooth_incomplete_denominator_block");
            status = ElementFactorOverBaseStatus::nonsmooth;
            flint::fmpz_mat_zero(exponents);
            return true;
        }

        SILEX_RELATION_PROFILE_EVENT(
                diagnostics,
                "relation.factor_over_base_nonintegral.path_denominator_block");
    }

    if (fmpz_is_one(remaining.raw()) == 0) {
        SILEX_RELATION_PROFILE_EVENT(
                diagnostics,
                "relation.factor_over_base_nonintegral.nonsmooth_denominator_leftover");
        status = ElementFactorOverBaseStatus::nonsmooth;
        flint::fmpz_mat_zero(exponents);
        return true;
    }

    status = ElementFactorOverBaseStatus::smooth;
    return true;
}

bool factor_nonintegral_element_over_base_by_denominator_split(
        ElementFactorOverBaseStatus& status,
        flint::FmpzMatRef exponents,
        const Element& alpha,
        const FactorBase& base,
        const DiagnosticsContext* diagnostics) noexcept {
    // reference `class_group_add_relation` and `_factor!` include the order
    // denominator in the nonintegral smoothness filter, then factor the
    // element valuations.  For maximal equation orders, this guarded native
    // slice uses the equivalent numerator row minus the rational denominator
    // ideal row, matching the existing C fractional-ideal factorization
    // without materializing the full fractional principal ideal.
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                        "relation.factor_over_base_nonintegral_element");
    status = ElementFactorOverBaseStatus::unsupported;

    const Order* parent = base.parent();
    if (!valid_base(base) || parent == nullptr ||
        parent->parent() == nullptr ||
        !alpha.has_parent(*parent->parent()) || alpha.equal_si(0) ||
        flint::fmpz_mat_nrows(exponents) != 1 ||
        flint::fmpz_mat_ncols(exponents) != base.length()) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.error_invalid_input",
                false);
    }
    if (!parent->is_maximal() || !parent->is_equation_order()) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.fallback_order_shape",
                true);
    }

    flint::FmpqPoly polynomial;
    flint::Fmpz denominator;
    flint::FmpzPoly numerator_polynomial;
    if (!alpha.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.error_polynomial",
                false);
    }
    fmpq_poly_get_denominator(denominator.raw(), polynomial.raw());
    if (fmpz_is_one(denominator.raw()) != 0) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.fallback_integral_element",
                true);
    }
    fmpq_poly_get_numerator(numerator_polynomial.raw(), polynomial.raw());

    Element numerator(*parent->parent());
    if (!numerator.is_defined() ||
        !set_element_from_fmpz_poly(numerator, numerator_polynomial)) {
        SILEX_RELATION_PROFILE_RETURN(
                diagnostics,
                "relation.factor_over_base_nonintegral.error_numerator_element",
                false);
    }

    flint::FmpzMat candidate(1, base.length());
    ElementFactorOverBaseStatus numerator_status =
            ElementFactorOverBaseStatus::unsupported;
    if (!factor_integral_element_over_base_by_residue_screen(
                numerator_status, flint::FmpzMatRef(candidate), numerator,
                base, nullptr, nullptr, nullptr, diagnostics)) {
        return false;
    }
    if (numerator_status != ElementFactorOverBaseStatus::smooth) {
        status = numerator_status;
        SILEX_RELATION_PROFILE_EVENT(
                diagnostics,
                numerator_status == ElementFactorOverBaseStatus::nonsmooth
                        ? "relation.factor_over_base_nonintegral.nonsmooth_numerator"
                        : "relation.factor_over_base_nonintegral.fallback_numerator_unsupported");
        return true;
    }

    ElementFactorOverBaseStatus denominator_status =
            ElementFactorOverBaseStatus::unsupported;
    if (!subtract_rational_denominator_row(
                denominator_status, flint::FmpzMatRef(candidate),
                flint::FmpzConstRef(denominator), base, diagnostics)) {
        return false;
    }
    if (denominator_status != ElementFactorOverBaseStatus::smooth) {
        status = denominator_status;
        return true;
    }

    SILEX_RELATION_PROFILE_EVENT(
            diagnostics, "relation.factor_over_base_nonintegral.smooth");
    status = ElementFactorOverBaseStatus::smooth;
    flint::fmpz_mat_set(exponents, flint::FmpzMatConstRef(candidate));
    return true;
}

}  // namespace

namespace detail {

bool set_relation_from_known_row(Relation& out,
                                 const FactorBase& base,
                                 const Element& generator,
                                 flint::FmpzMatConstRef row) noexcept {
    if (!valid_base(base) || base.parent() == nullptr ||
        base.parent()->parent() == nullptr ||
        !generator.has_parent(*base.parent()->parent()) ||
        generator.equal_si(0) || flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != base.length()) {
        return false;
    }

    // Preserve the buffers of a context-owned relation when its base is
    // unchanged; hot partial matches replace only the generator and row.
    if (out.is_defined() && out.base_.equal(base)) {
        if (!out.generator_.set(generator)) {
            return false;
        }
        flint::fmpz_mat_set(flint::FmpzMatRef(out.exponents_), row);
        out.has_relation_ = true;
        return true;
    }

    Relation candidate(base);
    if (!candidate.is_defined() || !candidate.generator_.set(generator)) {
        return false;
    }

    flint::fmpz_mat_set(flint::FmpzMatRef(candidate.exponents_), row);
    candidate.has_relation_ = true;
    out.swap(candidate);
    return true;
}

bool set_relation_from_integral_coordinates_and_norm(
        Relation& out,
        const Element& generator,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const DiagnosticsContext* diagnostics) noexcept {
    return set_relation_from_integral_coordinates_and_norm(
            out, generator, integral_coordinates, norm, nullptr, diagnostics);
}

bool set_relation_from_integral_coordinates_and_norm(
        Relation& out,
        const Element& generator,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    return out.set_generator_impl(generator, norm.raw(),
                                  integral_coordinates.raw(),
                                  integral_coordinate_polynomial == nullptr
                                          ? nullptr
                                          : integral_coordinate_polynomial->raw(),
                                  diagnostics);
}

bool set_relation_from_integral_coordinates_and_norm(
        Relation& out,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                        "relation.set_generator_integral_coordinates");
    const Order* parent = out.parent();
    if (!out.is_defined() || parent == nullptr ||
        parent->parent() == nullptr ||
        flint::fmpz_mat_nrows(integral_coordinates) != 1 ||
        flint::fmpz_mat_ncols(integral_coordinates) != parent->degree()) {
        SILEX_LOG(diagnostics, DiagnosticsModule::relation, LogLevel::detail,
                  "integral-coordinate relation generator rejected");
        return false;
    }

    flint::fmpz_mat_zero(flint::FmpzMatRef(out.scratch_exponents_));
    ElementFactorOverBaseStatus factor_status =
            ElementFactorOverBaseStatus::unsupported;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                            "relation.set_generator_integral_screen");
        if (!factor_integral_coordinates_over_base_by_residue_screen(
                    factor_status, flint::FmpzMatRef(out.scratch_exponents_),
                    out.base_, norm, integral_coordinates,
                    integral_coordinate_polynomial, diagnostics)) {
            return false;
        }
    }

    if (factor_status == ElementFactorOverBaseStatus::unsupported) {
        Element fallback_generator(*parent->parent());
        if (!fallback_generator.is_defined() ||
            !set_element_from_integral_coordinates(fallback_generator,
                                                   *parent,
                                                   integral_coordinates)) {
            return false;
        }
        return out.set_generator_impl(
                fallback_generator, norm.raw(), integral_coordinates.raw(),
                integral_coordinate_polynomial == nullptr
                        ? nullptr
                        : integral_coordinate_polynomial->raw(),
                diagnostics);
    }

    if (factor_status != ElementFactorOverBaseStatus::smooth) {
        SILEX_LOG(diagnostics, DiagnosticsModule::relation, LogLevel::detail,
                  "integral-coordinate relation generator is not smooth over factor base");
        return false;
    }

    Element next_generator(*parent->parent());
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::relation,
                "relation.set_generator_integral_coordinates_commit");
        if (!next_generator.is_defined() ||
            !set_element_from_integral_coordinates(next_generator, *parent,
                                                   integral_coordinates)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::relation,
                      LogLevel::detail,
                      "integral-coordinate relation generator copy failed");
            return false;
        }

        flint::fmpz_mat_set(flint::FmpzMatRef(out.exponents_),
                            flint::FmpzMatConstRef(out.scratch_exponents_));
        out.generator_.swap(next_generator);
        out.has_relation_ = true;
    }
    return true;
}

bool factor_relation_row_from_integral_coordinates_and_norm(
        Relation& out,
        bool& handled,
        bool& smooth,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    handled = false;
    smooth = false;
    const Order* parent = out.parent();
    if (!out.is_defined() || parent == nullptr ||
        parent->parent() == nullptr ||
        flint::fmpz_mat_nrows(integral_coordinates) != 1 ||
        flint::fmpz_mat_ncols(integral_coordinates) != parent->degree()) {
        return false;
    }

    out.has_relation_ = false;
    flint::fmpz_mat_zero(flint::FmpzMatRef(out.scratch_exponents_));
    ElementFactorOverBaseStatus factor_status =
            ElementFactorOverBaseStatus::unsupported;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                            "relation.set_generator_integral_screen");
        if (!factor_integral_coordinates_over_base_by_residue_screen(
                    factor_status, flint::FmpzMatRef(out.scratch_exponents_),
                    out.base_, norm, integral_coordinates,
                    integral_coordinate_polynomial, diagnostics)) {
            return false;
        }
    }
    if (factor_status == ElementFactorOverBaseStatus::unsupported) {
        return true;
    }
    handled = true;
    smooth = factor_status == ElementFactorOverBaseStatus::smooth;
    return true;
}

bool commit_relation_generator_from_integral_coordinates(
        Relation& out,
        flint::FmpzMatConstRef integral_coordinates,
        const DiagnosticsContext* diagnostics) noexcept {
    const Order* parent = out.parent();
    if (!out.is_defined() || parent == nullptr ||
        parent->parent() == nullptr ||
        flint::fmpz_mat_nrows(integral_coordinates) != 1 ||
        flint::fmpz_mat_ncols(integral_coordinates) != parent->degree()) {
        return false;
    }

    Element next_generator(*parent->parent());
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::relation,
            "relation.set_generator_integral_coordinates_commit");
    if (!next_generator.is_defined() ||
        !set_element_from_integral_coordinates(next_generator, *parent,
                                               integral_coordinates)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::relation,
                  LogLevel::detail,
                  "integral-coordinate relation generator copy failed");
        return false;
    }

    flint::fmpz_mat_set(flint::FmpzMatRef(out.exponents_),
                        flint::FmpzMatConstRef(out.scratch_exponents_));
    out.generator_.swap(next_generator);
    out.has_relation_ = true;
    return true;
}

flint::FmpzMatConstRef pending_relation_exponents_ref(
        const Relation& relation) noexcept {
    return flint::FmpzMatConstRef(relation.scratch_exponents_);
}

}  // namespace detail

Relation::Relation(const FactorBase& base) noexcept {
    define(base);
}

Relation::~Relation() noexcept = default;

Relation::Relation(Relation&& other) noexcept {
    swap(other);
}

Relation& Relation::operator=(Relation&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void Relation::swap(Relation& other) noexcept {
    base_.swap(other.base_);
    generator_.swap(other.generator_);
    exponents_.swap(other.exponents_);
    scratch_exponents_.swap(other.scratch_exponents_);
    std::swap(has_relation_, other.has_relation_);
}

void Relation::clear() noexcept {
    base_.clear();
    generator_.clear();
    exponents_ = flint::FmpzMat(0, 0);
    scratch_exponents_ = flint::FmpzMat(0, 0);
    has_relation_ = false;
}

bool Relation::define(const FactorBase& base) noexcept {
    if (!valid_base(base)) {
        return false;
    }

    Element next_generator(*base.parent()->parent());
    if (!next_generator.is_defined()) {
        return false;
    }
    flint::FmpzMat next_exponents(1, base.length());
    flint::FmpzMat next_scratch_exponents(1, base.length());
    FactorBase next_base;
    if (!next_base.set(base)) {
        return false;
    }

    clear();
    base_.swap(next_base);
    generator_ = std::move(next_generator);
    exponents_ = std::move(next_exponents);
    scratch_exponents_ = std::move(next_scratch_exponents);
    return true;
}

bool Relation::set(const Relation& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    Relation copy(other.base_);
    if (!copy.is_defined()) {
        return false;
    }

    if (other.has_relation_) {
        if (!copy.generator_.set(other.generator_)) {
            return false;
        }
        flint::fmpz_mat_set(flint::FmpzMatRef(copy.exponents_),
                            flint::FmpzMatConstRef(other.exponents_));
        copy.has_relation_ = true;
    }

    swap(copy);
    return true;
}

bool Relation::is_defined() const noexcept {
    return valid_base(base_) && generator_.is_defined() &&
           flint::fmpz_mat_nrows(exponents_) == 1 &&
           flint::fmpz_mat_ncols(exponents_) == base_.length() &&
           flint::fmpz_mat_nrows(scratch_exponents_) == 1 &&
           flint::fmpz_mat_ncols(scratch_exponents_) == base_.length();
}

bool Relation::is_set() const noexcept {
    return is_defined() && has_relation_;
}

const FactorBase* Relation::factor_base() const noexcept {
    return is_defined() ? &base_ : nullptr;
}

const Order* Relation::parent() const noexcept {
    return is_defined() ? base_.parent() : nullptr;
}

slong Relation::length() const noexcept {
    return is_defined() ? base_.length() : 0;
}

bool Relation::set_generator(
        const Element& alpha,
        const DiagnosticsContext* diagnostics) noexcept {
    return set_generator_impl(alpha, nullptr, nullptr, nullptr, diagnostics);
}

bool Relation::set_generator_with_norm(
        const Element& alpha,
        flint::FmpqConstRef norm,
        const DiagnosticsContext* diagnostics) noexcept {
    return set_generator_impl(alpha, norm.raw(), nullptr, nullptr,
                              diagnostics);
}

bool Relation::set_generator_impl(
        const Element& alpha,
        const fmpq* known_norm,
        const fmpz_mat_struct* known_integral_coordinates,
        const fmpz_poly_struct* known_integral_polynomial,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                        "relation.set_generator");
    if (!is_defined() || parent() == nullptr ||
        parent()->parent() == nullptr ||
        !alpha.has_parent(*parent()->parent()) || alpha.equal_si(0)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::relation, LogLevel::detail,
                  "relation generator rejected");
        return false;
    }

    Element next_generator(*parent()->parent());
    if (!next_generator.is_defined()) {
        return false;
    }
    flint::fmpz_mat_zero(flint::FmpzMatRef(scratch_exponents_));

    ElementFactorOverBaseStatus factor_status =
            ElementFactorOverBaseStatus::unsupported;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                            "relation.set_generator_integral_screen");
        if (!factor_integral_element_over_base_by_residue_screen(
                    factor_status, flint::FmpzMatRef(scratch_exponents_),
                    alpha, base_, known_norm, known_integral_coordinates,
                    known_integral_polynomial, diagnostics)) {
            return false;
        }
    }

    if (factor_status == ElementFactorOverBaseStatus::unsupported) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                            "relation.set_generator_nonintegral_split");
        if (!factor_nonintegral_element_over_base_by_denominator_split(
                    factor_status, flint::FmpzMatRef(scratch_exponents_),
                    alpha, base_, diagnostics)) {
            return false;
        }
    }

    if (factor_status == ElementFactorOverBaseStatus::unsupported &&
        known_integral_coordinates != nullptr) {
        SILEX_RELATION_PROFILE_EVENT(
                diagnostics,
                "relation.set_generator.integral_ideal_fallback");
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                            "relation.set_generator_integral_ideal_fallback");
        OrderElement integral_alpha(*parent());
        Ideal ideal(*parent());
        if (!integral_alpha.is_defined() || !ideal.is_defined()) {
            return false;
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                                "relation.principal_integral_ideal");
            if (!integral_alpha.set_coordinates(
                        flint::FmpzMatConstRef(
                                known_integral_coordinates)) ||
                !ideal.set_principal(integral_alpha)) {
                return false;
            }
        }
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                                "relation.factor_over_base");
            factor_status = ideal_factor_over_base(
                                    flint::FmpzMatRef(scratch_exponents_),
                                    ideal, base_, diagnostics)
                    ? ElementFactorOverBaseStatus::smooth
                    : ElementFactorOverBaseStatus::nonsmooth;
        }
    }

    if (factor_status == ElementFactorOverBaseStatus::unsupported) {
        SILEX_RELATION_PROFILE_EVENT(
                diagnostics,
                "relation.set_generator.fallback_fractional_ideal");
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                                "relation.set_generator_fractional_fallback");
            FractionalIdeal ideal(*parent());
            if (!ideal.is_defined()) {
                return false;
            }
            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                                    "relation.principal_ideal");
                if (!ideal.set_principal(alpha, diagnostics)) {
                    return false;
                }
            }
            {
                SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                                    "relation.factor_over_base");
                if (!ideal_factor_over_base(
                            flint::FmpzMatRef(scratch_exponents_), ideal,
                            base_, diagnostics)) {
                    factor_status = ElementFactorOverBaseStatus::nonsmooth;
                } else {
                    factor_status = ElementFactorOverBaseStatus::smooth;
                }
            }
        }
    }

    if (factor_status != ElementFactorOverBaseStatus::smooth) {
        SILEX_LOG(diagnostics, DiagnosticsModule::relation,
                  LogLevel::detail,
                  "relation generator is not smooth over factor base");
        return false;
    }

    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::relation,
                            "relation.set_generator_commit");
        if (!next_generator.set(alpha)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::relation,
                      LogLevel::detail, "relation generator copy failed");
            return false;
        }

        flint::fmpz_mat_set(flint::FmpzMatRef(exponents_),
                            flint::FmpzMatConstRef(scratch_exponents_));
        generator_.swap(next_generator);
        has_relation_ = true;
    }
    return true;
}

bool Relation::generator(Element& out) const noexcept {
    if (!is_set() || out.parent() == nullptr ||
        generator_.parent() == nullptr ||
        !out.parent()->has_same_data(*generator_.parent())) {
        return false;
    }
    return out.set(generator_);
}

bool Relation::exponents(flint::FmpzMatRef out) const noexcept {
    if (!is_set() || flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != length()) {
        return false;
    }
    flint::fmpz_mat_set(out, flint::FmpzMatConstRef(exponents_));
    return true;
}

flint::FmpzMatConstRef Relation::exponents_ref() const noexcept {
    return flint::FmpzMatConstRef(exponents_);
}

RelationMatrix::RelationMatrix(const FactorBase& base) noexcept {
    define(base);
}

RelationMatrix::~RelationMatrix() noexcept = default;

RelationMatrix::RelationMatrix(RelationMatrix&& other) noexcept {
    swap(other);
}

RelationMatrix& RelationMatrix::operator=(RelationMatrix&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void RelationMatrix::swap(RelationMatrix& other) noexcept {
    base_.swap(other.base_);
    rows_.swap(other.rows_);
    generators_.swap(other.generators_);
}

void RelationMatrix::clear() noexcept {
    base_.clear();
    rows_ = fmpz_smat::SparseMat(0);
    generators_.clear();
}

bool RelationMatrix::define(const FactorBase& base) noexcept {
    if (!valid_base(base)) {
        return false;
    }

    clear();
    if (!base_.set(base)) {
        return false;
    }
    rows_ = fmpz_smat::SparseMat(base.length());
    return true;
}

bool RelationMatrix::set(const RelationMatrix& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    RelationMatrix copy(other.base_);
    if (!copy.is_defined() || !copy.rows_.set(other.rows_)) {
        return false;
    }

    copy.generators_.reserve(other.generators_.size());
    for (const Element& generator : other.generators_) {
        Element copied(*generator.parent());
        if (!copied.is_defined() || !copied.set(generator)) {
            return false;
        }
        copy.generators_.emplace_back(std::move(copied));
    }

    swap(copy);
    return true;
}

bool RelationMatrix::is_defined() const noexcept {
    return valid_base(base_) &&
           rows_.ncols() == base_.length() &&
           rows_.nrows() == static_cast<slong>(generators_.size());
}

const FactorBase* RelationMatrix::factor_base() const noexcept {
    return is_defined() ? &base_ : nullptr;
}

const Order* RelationMatrix::parent() const noexcept {
    return is_defined() ? base_.parent() : nullptr;
}

slong RelationMatrix::length() const noexcept {
    return rows_.nrows();
}

slong RelationMatrix::ncols() const noexcept {
    return rows_.ncols();
}

bool RelationMatrix::append(const Relation& relation) noexcept {
    if (!is_defined() || !relation.is_set() ||
        relation.factor_base() == nullptr ||
        !relation.factor_base()->equal(base_)) {
        return false;
    }

    Element generator_buffer(*parent()->parent());
    if (!relation.generator(generator_buffer) ||
        !rows_.append_fmpz_mat_row(
                relation.exponents_ref(), 0)) {
        return false;
    }

    generators_.emplace_back(std::move(generator_buffer));
    return true;
}

bool RelationMatrix::row(flint::FmpzMatRef out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length() ||
        flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != ncols()) {
        return false;
    }

    rows_.row_ref(index).get_fmpz_mat_row(out, 0);
    return true;
}

bool RelationMatrix::row_equal(flint::FmpzMatConstRef row,
                               slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length() ||
        flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != ncols()) {
        return false;
    }

    return sparse_row_equal_fmpz_mat_row(rows_.row_ref(index), row, 0,
                                         ncols());
}

bool RelationMatrix::row_first_nonzero(slong& out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length()) {
        return false;
    }

    const fmpz_smat::SparseRow& row = rows_.row_ref(index);
    if (row.length() == 0) {
        out = ncols();
    } else {
        out = row.col(0);
    }
    return true;
}

bool RelationMatrix::rows(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || flint::fmpz_mat_nrows(out) != length() ||
        flint::fmpz_mat_ncols(out) != ncols()) {
        return false;
    }

    rows_.get_fmpz_mat(out);
    return true;
}

bool RelationMatrix::generator(Element& out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length() ||
        parent() == nullptr || parent()->parent() == nullptr ||
        !out.has_parent(*parent()->parent())) {
        return false;
    }

    return out.set(generators_[static_cast<std::size_t>(index)]);
}

const Element* RelationMatrix::generator_at(slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length()) {
        return nullptr;
    }

    return &generators_[static_cast<std::size_t>(index)];
}

bool RelationMatrix::to_abelian_group(FiniteAbelianGroup& out) const noexcept {
    if (!is_defined()) {
        return false;
    }

    flint::FmpzMat dense_rows(length(), ncols());
    FiniteAbelianGroup candidate;
    if (!rows(flint::FmpzMatRef(dense_rows)) ||
        !candidate.set_relation_matrix(flint::FmpzMatConstRef(dense_rows))) {
        return false;
    }

    out.swap(candidate);
    return true;
}

#undef SILEX_RELATION_PROFILE_RETURN
#undef SILEX_RELATION_PROFILE_EVENT

}  // namespace silex
