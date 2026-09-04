#include "order_unit_internal.hpp"
#include "relation_unit_internal.hpp"

#include "../class_group/class_group_internal.hpp"

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

bool degree_one_roots_mod_q_nmod(std::vector<flint::Fmpz>& roots,
                                 const fmpq_poly_struct* polynomial,
                                 ulong q,
                                 const DiagnosticsContext* diagnostics)
        noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::unit_group,
            "unit_group.stable_scan.degree_one_roots_nmod");
    NmodPoly reduced(q);
    NmodPoly x_poly(q);
    NmodPoly reduced_reverse(q);
    NmodPoly reduced_preinverse(q);
    NmodPoly degree_one_part(q);
    NmodPolyFactor factorization;

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_nmod_reduce_poly");
        fmpq_poly_get_integral_nmod_poly(reduced, polynomial, q);
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_nmod_frobenius_gcd");
        nmod_poly_set_coeff_ui(x_poly.raw(), 1, 1);
        nmod_poly_reverse(reduced_reverse.raw(), reduced.raw(),
                          reduced.raw()->length);
        nmod_poly_inv_series_newton(reduced_preinverse.raw(),
                                    reduced_reverse.raw(),
                                    reduced_reverse.raw()->length);
        nmod_poly_powmod_x_ui_preinv(degree_one_part.raw(), q, reduced.raw(),
                                     reduced_preinverse.raw());
        nmod_poly_sub(degree_one_part.raw(), degree_one_part.raw(),
                      x_poly.raw());
        nmod_poly_gcd(degree_one_part.raw(), degree_one_part.raw(),
                      reduced.raw());
    }

    const slong degree_one_degree = nmod_poly_degree(degree_one_part.raw());
    if (degree_one_degree < 1) {
        return true;
    }
    if (degree_one_degree == 1) {
        flint::Fmpz root;
        if (degree_one_root_from_nmod_poly(root, degree_one_part.raw())) {
            roots.push_back(std::move(root));
        }
        return true;
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.stable_scan.degree_one_roots_nmod_factor");
        nmod_poly_roots(factorization.raw(), degree_one_part.raw(), 0);
    }

    roots.reserve(static_cast<std::size_t>(factorization.raw()->num));
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        flint::Fmpz root;
        if (degree_one_root_from_nmod_poly(root, factorization.raw()->p + i)) {
            roots.push_back(std::move(root));
        }
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

    if (order.is_equation_order() && flint::fmpz_abs_fits_ui(q)) {
        const ulong q_ui = flint::fmpz_get_ui(q);
        if (q_ui > 1) {
            return degree_one_roots_mod_q_nmod(roots, raw_field->pol, q_ui,
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
        // For equation orders, the caller has already skipped rational primes
        // dividing the order discriminant, so the reduced defining polynomial
        // is squarefree at this point.  Keep the guard for maximal orders whose
        // defining polynomial may still have index divisors.
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
};

bool select_stable_proof_primes(
        slong& local_prime_count,
        bool& certified,
        flint::FmpzMat& kernel,
        const OrderUnitGroup& group,
        flint::FmpzConstRef ell,
        double stable_threshold) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_proof_primes");
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
    // and stops after the candidate kernel remains stable for 3.5 times its
    // current dimension. This avoids treating the old auxiliary bound as a
    // proof condition without introducing a local cap.
    auto process_prime_q =
            [&](flint::FmpzConstRef q) noexcept -> ProofPrimeScanResult {
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

            if (proof_rank == 1) {
                bool nonzero = false;
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.rank_one_image_test");
                if (!rank_one_proof_source_nonzero(
                            nonzero, rank_one_source, group, *prime, ell)) {
                    SILEX_PROFILE_EVENT(
                            group.diagnostics(),
                            DiagnosticsModule::unit_group,
                            "unit_group.stable_scan.rank_one_image_unavailable");
                    continue;
                }

                const ProofPrimeScanResult result =
                        process_rank_one_image(nonzero);
                if (result != ProofPrimeScanResult::keep_scanning) {
                    return result;
                }
                continue;
            }

            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.dlog_column_test");
            flint::FmpzMat proof_column(proof_rank, 1);
            if (!detail::saturation_proof_prime_column(
                        proof_column, group, *prime, ell)) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.dlog_column_unavailable");
                continue;
            }
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_scan.dlog_column");
            ++local_prime_count;
            if (!append_column(proof_matrix, proof_column)) {
                return ProofPrimeScanResult::failed;
            }
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.proof_dlog_kernel_nullspace");
            if (!detail::dlog_kernel_from_matrix(current_kernel, proof_matrix,
                                                 ell)) {
                return ProofPrimeScanResult::failed;
            }

            if (flint::fmpz_mat_nrows(current_kernel) == 0) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.stable_scan.certified_full_rank");
                kernel = std::move(current_kernel);
                certified = true;
                return ProofPrimeScanResult::done;
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
                return ProofPrimeScanResult::done;
            }
        }

        return ProofPrimeScanResult::keep_scanning;
    };

    const ulong max_q_ui = std::numeric_limits<ulong>::max();
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
                    if (result == ProofPrimeScanResult::failed) {
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
                if (result == ProofPrimeScanResult::failed) {
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
                                  EmbeddingContext& embeddings,
                                  slong precision) noexcept {
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
    flint::FmpzMat kernel(0, rank);
    SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_prime_selection");
    if (!select_stable_proof_primes(
                local_prime_count, certified, kernel, group, ell,
                stable_threshold)) {
        step.status = StableRelationSaturationStatus::unavailable;
        return true;
    }
    step.local_primes = local_prime_count;

    if (flint::fmpz_mat_nrows(kernel) == 0) {
        step.status = StableRelationSaturationStatus::verified;
        return true;
    }

    detail::CompactFieldModulusCache field_modulus_cache;
    bool saw_candidate = false;
    for (slong row = 0; row < flint::fmpz_mat_nrows(kernel); ++row) {
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_relation_saturation.row_inspected");
        if (detail::kernel_row_divisible(
                    flint::FmpzMatConstRef(kernel), row, rank, ell)) {
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
        if (!detail::kernel_row_root(is_power, root, group,
                                     flint::FmpzMatConstRef(kernel), row,
                                     ell, &field_modulus_cache)) {
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_unavailable");
            step.status = StableRelationSaturationStatus::unavailable;
            return true;
        }
        if (!is_power) {
            SILEX_PROFILE_EVENT(
                    out.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.stable_relation_saturation.root_not_power");
            step.status = StableRelationSaturationStatus::wasted_non_power;
            return true;
        }
        SILEX_PROFILE_EVENT(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.stable_relation_saturation.root_power");

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
    step.status = saw_candidate ? StableRelationSaturationStatus::wasted_non_power
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
    const slong rank = group.free_rank();
    if (!group.is_set() || order == nullptr || !flint::fmpz_is_prime(ell) ||
        target_len < 0 || fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    if (target_len == 0) {
        PrimeIdealList candidate(*order, 0);
        flint::FmpzMat candidate_kernel(rank, 0);
        if (!candidate.is_defined()) {
            return false;
        }
        out.swap(candidate);
        kernel = std::move(candidate_kernel);
        return true;
    }

    std::vector<PrimeIdeal> selected;
    selected.reserve(static_cast<std::size_t>(target_len));
    flint::FmpzMat dlog_matrix(rank, 0);

    auto try_append_prime = [&](const PrimeIdeal& prime) noexcept -> bool {
        flint::FmpzMat column(rank, 1);
        if (!detail::saturation_prime_column(column, group, prime, ell)) {
            return true;
        }
        selected.emplace_back(*order);
        if (!selected.back().set(prime) ||
            !append_column(dlog_matrix, column)) {
            return false;
        }
        return true;
    };

    flint::Fmpz p;
    flint::Fmpz pminus;
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (static_cast<slong>(selected.size()) < target_len &&
           fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        fmpz_sub_ui(pminus.raw(), p.raw(), 1);
        if (fmpz_divisible(pminus.raw(), ell.raw()) == 0) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.saturation_selector.congruence_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        PrimeIdealList local;
        {
            SILEX_PROFILE_SCOPE(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.saturation_selector_decompose_degree_one");
            if (!decompose_prime(local, *order, flint::FmpzConstRef(p), 1)) {
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
        }
        for (slong i = 0; i < local.size() &&
                          static_cast<slong>(selected.size()) < target_len;
             ++i) {
            const PrimeIdeal* prime = local.at(i);
            if (prime != nullptr && !try_append_prime(*prime)) {
                return false;
            }
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    if (static_cast<slong>(selected.size()) == target_len) {
        return publish_saturation_selection(out, kernel, *order, selected,
                                            dlog_matrix, target_len, ell);
    }

    SILEX_PROFILE_EVENT(
            group.diagnostics(), DiagnosticsModule::unit_group,
            "unit_group.saturation_selector.source_shape_fallback");
    selected.clear();
    flint::FmpzMat fallback_matrix(rank, 0);
    dlog_matrix = std::move(fallback_matrix);
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (static_cast<slong>(selected.size()) < target_len &&
           fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList local;
        if (decompose_prime(local, *order, flint::FmpzConstRef(p))) {
            for (slong i = 0; i < local.size() &&
                              static_cast<slong>(selected.size()) < target_len;
                 ++i) {
                const PrimeIdeal* prime = local.at(i);
                if (prime != nullptr && !try_append_prime(*prime)) {
                    return false;
                }
            }
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    return publish_saturation_selection(out, kernel, *order, selected,
                                        dlog_matrix, target_len, ell);
}

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
                    !detail::set_degree_one_prime_ideal_from_root(
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
        flint::fmpz_mat_ncols(kernel_rows) != rank ||
        !flint::fmpz_is_prime(ell) || !flint::fmpz_fits_si(ell)) {
        return false;
    }

    if (detail::kernel_row_divisible(kernel_rows, row, rank, ell)) {
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
    if (!detail::kernel_row_root(is_power, root, group, kernel_rows, row, ell,
                                 field_modulus_cache)) {
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

    auto publish_selected =
            [&](std::vector<PrimeIdeal>& selected_primes) noexcept -> bool {
        if (static_cast<slong>(selected_primes.size()) != target_len) {
            return false;
        }

        PrimeIdealList candidate(*order, target_len);
        if (!candidate.is_defined()) {
            return false;
        }
        for (slong i = 0; i < target_len; ++i) {
            PrimeIdeal* dest = candidate.at(i);
            if (dest == nullptr ||
                !dest->set(selected_primes[static_cast<std::size_t>(i)])) {
                return false;
            }
        }

        out.swap(candidate);
        return true;
    };

    std::vector<PrimeIdeal> selected;
    selected.reserve(static_cast<std::size_t>(target_len));

    flint::Fmpz order_discriminant;
    if (!order_discriminant_for_saturation(order_discriminant, *order)) {
        return false;
    }

    flint::Fmpz p;
    flint::Fmpz pminus;
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (static_cast<slong>(selected.size()) < target_len &&
           fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        // reference RelSaturate.compute_candidates_for_saturate scans rational
        // q = 1 mod ell and uses only degree-one primes above q.
        fmpz_sub_ui(pminus.raw(), p.raw(), 1);
        if (fmpz_divisible(pminus.raw(), ell.raw()) == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.saturation_selector.congruence_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }
        if (rational_prime_divides_order_discriminant(
                    flint::FmpzConstRef(order_discriminant),
                    flint::FmpzConstRef(p))) {
            SILEX_PROFILE_EVENT(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.saturation_selector.discriminant_skip");
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        PrimeIdealList local;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics_, DiagnosticsModule::unit_group,
                    "unit_group.saturation_selector_decompose_degree_one");
            if (!decompose_prime(local, *order, flint::FmpzConstRef(p), 1)) {
                fmpz_nextprime(p.raw(), p.raw(), 1);
                continue;
            }
        }
        if (local.size() > 0) {
            for (slong i = 0; i < local.size() &&
                              static_cast<slong>(selected.size()) < target_len;
                 ++i) {
                const PrimeIdeal* prime = local.at(i);
                if (prime != nullptr &&
                    detail::saturation_prime_usable(*this, *prime, ell)) {
                    selected.emplace_back(*order);
                    if (!selected.back().set(*prime)) {
                        return false;
                    }
                }
            }
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    if (publish_selected(selected)) {
        return true;
    }

    SILEX_PROFILE_EVENT(
            diagnostics_, DiagnosticsModule::unit_group,
            "unit_group.saturation_selector.source_shape_fallback");
    selected.clear();
    flint::fmpz_set_ui(flint::FmpzRef(p), 2);
    while (static_cast<slong>(selected.size()) < target_len &&
           fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList local;
        if (decompose_prime(local, *order, flint::FmpzConstRef(p))) {
            for (slong i = 0; i < local.size() &&
                              static_cast<slong>(selected.size()) < target_len;
                 ++i) {
                const PrimeIdeal* prime = local.at(i);
                if (prime != nullptr &&
                    detail::saturation_prime_usable(*this, *prime, ell)) {
                    selected.emplace_back(*order);
                    if (!selected.back().set(*prime)) {
                        return false;
                    }
                }
            }
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    return publish_selected(selected);
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
                !detail::set_degree_one_prime_ideal_from_root(
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

    flint::FmpzMat kernel(0, group.free_rank());
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
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_local_with_kernel");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        flint::fmpz_mat_ncols(kernel) != group.free_rank() ||
        !flint::fmpz_is_prime(ell) || precision <= 0) {
        return false;
    }

    detail::CompactFieldModulusCache field_modulus_cache;
    bool saw_failure = false;
    for (slong i = 0; i < flint::fmpz_mat_nrows(kernel); ++i) {
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.saturate_local.kernel_row_attempt");
        OrderUnitGroup row_candidate(*order);
        row_candidate.set_diagnostics(diagnostics_);
        bool row_changed = false;
        if (saturate_row_with_cache(row_candidate, row_changed, group,
                                    kernel, i, ell, embeddings, precision,
                                    &field_modulus_cache)) {
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
        flint::FmpzMat kernel(0, working.free_rank());
        if (!select_saturation_primes_with_kernel(
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
                    ell, embeddings, precision)) {
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

bool OrderUnitGroup::saturate_index_bounded(bool& changed,
                                            bool& stable,
                                            const OrderUnitGroup& group,
                                            EmbeddingContext& embeddings,
                                            slong aux_target_len,
                                            flint::FmpzConstRef aux_bound,
                                            slong max_passes,
                                            slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_index_bounded");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        aux_target_len <= 0 || fmpz_cmp_ui(aux_bound.raw(), 2) < 0 ||
        max_passes < 0 || precision <= 0) {
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
            if (!set(working)) {
                return false;
            }
            try_certify_index_one(precision);
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

    if (!set(working)) {
        return false;
    }
    try_certify_index_one(precision);
    changed = any_changed;
    stable = all_stable;
    return true;
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
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.saturate_index_bounded_adaptive");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        aux_target_len <= 0 || fmpz_cmp_ui(aux_bound_start.raw(), 2) < 0 ||
        fmpz_cmp(aux_bound_max.raw(), aux_bound_start.raw()) < 0 ||
        max_passes <= 0 || precision <= 0) {
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

    flint::Fmpz aux_bound;
    flint::Fmpz next_bound;
    flint::fmpz_set(flint::FmpzRef(aux_bound), aux_bound_start);
    bool any_changed = false;

    for (;;) {
        bool pass_changed = false;
        bool pass_stable = false;
        if (pass_result.saturate_index_bounded(
                    pass_changed, pass_stable, working, embeddings,
                    aux_target_len, flint::FmpzConstRef(aux_bound),
                    max_passes, precision)) {
            if (pass_changed || pass_stable) {
                working.swap(pass_result);
                if (pass_changed) {
                    any_changed = true;
                }
            }

            if (pass_stable) {
                if (!set(working)) {
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
            if (!set(working)) {
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
            if (!set(working)) {
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

bool OrderUnitGroup::prove_local_saturated(ProofState& status,
                                           bool& changed,
                                           const OrderUnitGroup& group,
                                           flint::FmpzConstRef ell,
                                           slong aux_target_len,
                                           flint::FmpzConstRef aux_bound,
                                           EmbeddingContext& embeddings,
                                           slong precision) noexcept {
    return prove_local_saturated_(status, changed, group, ell, aux_target_len,
                                  aux_bound, embeddings, precision, false);
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
        bool use_stable_proof_fallback) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.prove_local_saturated");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || aux_target_len <= 0 ||
        fmpz_cmp_ui(aux_bound.raw(), 2) < 0 || precision <= 0) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *field) || group.free_rank() != rank) {
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

    auto run_stable_fallback = [&]() noexcept -> bool {
        constexpr double kInitialStable = 3.5;
        double stable_threshold = kInitialStable;
        bool stable_changed = false;
        for (;;) {
            StableRelationSaturationStep step;
            if (!stable_relation_saturation_step(
                        pass_result, step, working, ell, stable_threshold,
                        embeddings, precision)) {
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
                ell, embeddings, precision);
    } else if (first_prime != nullptr) {
        SILEX_PROFILE_EVENT(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.proof_precomputed_kernel_unusable");
        pass_ok = pass_result.saturate_local_once(
                pass_changed, working,
                PrimeIdealSpan(first_prime,
                               static_cast<std::size_t>(primes.size())),
                ell, embeddings, precision);
    }

    if (pass_ok) {
        if (pass_changed) {
            return publish(pass_result, ProofState::unavailable,
                           local_prime_count, true);
        }
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
        slong stable_proof_record_prefix_len) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::unit_group,
            "unit_group.prove_local_saturated_stable_in_place");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!is_defined() || !group.is_set() ||
        !same_order_parent(parent(), order) ||
        !embedding_has_parent(embeddings, field) ||
        !flint::fmpz_is_prime(ell) || aux_target_len <= 0 ||
        fmpz_cmp_ui(aux_bound.raw(), 2) < 0 || precision <= 0) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *field) || group.free_rank() != rank) {
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

    constexpr double kInitialStable = 3.5;
    double stable_threshold = kInitialStable;
    bool stable_changed = false;
    for (;;) {
        OrderUnitGroup pass_result(*order);
        pass_result.set_diagnostics(diagnostics_);
        if (!pass_result.is_defined()) {
            return false;
        }

        StableRelationSaturationStep step;
        if (!stable_relation_saturation_step(
                    pass_result, step, group, ell, stable_threshold,
                    embeddings, precision)) {
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
                        stable_proof_record_prefix_len);
            } else {
                pass_ok = pass_result.prove_local_saturated_(
                        pass_status, pass_changed, working,
                        flint::FmpzConstRef(ell), aux_target_len, aux_bound,
                        embeddings, precision, true);
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
