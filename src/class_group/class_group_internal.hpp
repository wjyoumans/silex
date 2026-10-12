#pragma once

#include <silex/class_group.hpp>
#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>

#include <cstddef>
#include <vector>

#include "class_group_storage_internal.hpp"
#include "ideal_minkowski_embedding_internal.hpp"

namespace silex {
class EmbeddingContext;
}  // namespace silex

namespace silex::detail {

struct RelationUnitCandidateBatchState;

struct ClassGroupRelationOptions {
    slong coordinate_search_radius = 2;
    slong ideal_search_radius = 1;
    slong target_relation_kernel_units = 0;
    slong post_finite_refinement_phase_budget = 0;
    slong max_candidates = WORD_MAX;
    slong max_relations = WORD_MAX;
    slong relation_saturation_aux_prime_bound = 0;
    slong relation_saturation_max_appends_per_ell = 0;
    ulong zeta_bf_max_cutoff = 0;
    CertificationMode requested_certification = CertificationMode::unknown;
    const DiagnosticsContext* diagnostics = nullptr;
};

inline ClassGroupRelationOptions class_group_relation_options(
        const Order& order,
        const ClassGroupComputeOptions& requested) noexcept {
    ClassGroupRelationOptions options;
    if (order.degree() >= 3) {
        options.ideal_search_radius = 12;
    }
    options.max_candidates = requested.max_candidates;
    options.max_relations = requested.max_relations;
    options.relation_saturation_aux_prime_bound =
            requested.relation_saturation_aux_prime_bound;
    options.relation_saturation_max_appends_per_ell =
            requested.relation_saturation_max_appends_per_ell;
    options.zeta_bf_max_cutoff = requested.zeta_bf_max_cutoff;
    options.requested_certification = requested.requested_certification;
    options.diagnostics = requested.diagnostics;
    return options;
}

// The class number h(D) of the imaginary quadratic maximal order of
// fundamental discriminant D < 0 with |D| < 2^63: counted from the reduced
// forms for |D| < 2^44, from Dirichlet's class number formula with a
// rigorous L(1, chi_D) enclosure (imaginary_quadratic_class_number_analytic)
// above.  Returns false, leaving `out` unspecified, if neither gives h(D).
bool imaginary_quadratic_class_number(
        flint::FmpzRef out,
        flint::FmpzConstRef discriminant) noexcept;

// Test controls and record of the steps of
// `imaginary_quadratic_class_number_analytic`.
struct ImaginaryQuadraticClassNumberAudit {
    // When true, the first step uses `first_log2_tolerance` (at most zero)
    // in place of the tolerance sized for h(D), so a test can make that step
    // fail to isolate h(D).
    bool override_first_log2_tolerance = false;
    slong first_log2_tolerance = 0;
    // The step that gave h(D): 0 for the tolerance sized for h(D), 1 and 2
    // for FLINT's precision mapping at the bit length of |D| plus 32 and at
    // twice that; -1 when none did.
    int step = -1;
};

// h(D) for a fundamental discriminant D < -4 with |D| < 2^63 from
// h(D) = sqrt|D| L(1, chi_D) / pi, returned only when the ball enclosing it
// contains exactly one integer; false otherwise.  D is not checked to be
// fundamental.  `audit`, when given, is used and updated as described above.
bool imaginary_quadratic_class_number_analytic(
        flint::FmpzRef out,
        flint::FmpzConstRef discriminant,
        ImaginaryQuadraticClassNumberAudit* audit = nullptr) noexcept;

struct WitnessedClassRelationHnfBasis {
    flint::FmpzMat rows{0, 0};
    flint::FmpzMat relation_coefficients{0, 0};
    std::vector<FactoredElement> witnesses;
};

enum class ClassGroupIdealRelationWitnessStatus {
    not_started,
    success,
    invalid_input,
    arithmetic_failure,
    factorization_failure,
    verification_failure,
    exhausted,
};

enum class ClassGroupIdealRelationWitnessStage {
    none,
    input_validation,
    direct_factorization,
    ideal_reduction,
    reduced_factorization,
    continuation_setup,
    continuation_enumeration,
    continuation_factorization,
    random_product_construction,
    random_product_enumeration,
    exact_verification,
    search_exhaustion,
};

struct ClassGroupIdealRelationWitnessResult {
    ClassGroupIdealRelationWitnessStatus status =
            ClassGroupIdealRelationWitnessStatus::not_started;
    ClassGroupIdealRelationWitnessStage stage =
            ClassGroupIdealRelationWitnessStage::none;
    slong candidates_tried = 0;
    slong random_products_tried = 0;
    bool used_reduction = false;
    bool used_random_product = false;
};

struct ClassGroupIdealRelationWitnessOptions {
    slong reduction_precision = 200;
    slong max_candidates = 500;
    slong max_candidates_per_ideal = 64;
    slong max_random_products = 8;
    ulong random_seed = UWORD(0x6a09e667f3bcc909);
};

const char* class_group_ideal_relation_witness_status_name(
        ClassGroupIdealRelationWitnessStatus status) noexcept;
const char* class_group_ideal_relation_witness_stage_name(
        ClassGroupIdealRelationWitnessStage stage) noexcept;

bool verify_class_group_ideal_relation_witness(
        const ClassGroupContext& context,
        const Ideal& ideal,
        const FactoredElement& multiplier,
        flint::FmpzMatConstRef factor_base_row) noexcept;
bool factor_base_row_ideal(
        FractionalIdeal& out,
        const FactorBase& base,
        flint::FmpzMatConstRef factor_base_row) noexcept;
bool multiply_factored_element_power_fmpz(
        FactoredElement& accumulator,
        const FactoredElement& base,
        flint::FmpzConstRef exponent) noexcept;
// Multiplies accumulator by base^exponent. Exponents that fit in slong are
// pushed as one factor; larger ones go through
// multiply_factored_element_power_fmpz, which keeps its exponent guard.
bool multiply_element_power_fmpz(FactoredElement& accumulator,
                                 const Element& base,
                                 flint::FmpzConstRef exponent) noexcept;

bool class_group_ideal_relation_continuation_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& original_ideal,
        const Ideal& reduced_ideal,
        const Element& reduction_multiplier,
        const ClassGroupIdealRelationWitnessOptions& options = {}) noexcept;

bool class_group_ideal_relation_witness(
        ClassGroupIdealRelationWitnessResult& result,
        FactoredElement& multiplier,
        flint::FmpzMatRef factor_base_row,
        const ClassGroupContext& context,
        const Ideal& ideal,
        const ClassGroupIdealRelationWitnessOptions& options = {}) noexcept;

// One result of compact_field_modulus_mod_p: whether the defining polynomial
// of the cache's field is usable modulo `prime` for `threshold`, and if so
// its monic reduction with integer coefficients in [0, prime).  An entry
// depends only on the defining polynomial, `prime` and `threshold`.
struct CompactFieldModulusCacheEntry {
    flint::Fmpz prime;
    flint::FmpzPoly modulus;
    slong threshold = 0;
    bool usable = false;
};

// Cache of compact_field_modulus_mod_p results for one number field.
// `field` shares (and so keeps alive) the field data the entries were
// computed for; a lookup for a field without the same data clears the
// cache first.  `entries` is sorted by (threshold, prime), with at most one
// entry per pair.
struct CompactFieldModulusCache {
    NumberField field;
    std::vector<CompactFieldModulusCacheEntry> entries;
};

// Returns the entry of `cache` for (`prime`, `threshold`) over `field`, or
// nullptr on a miss.  `position` receives the entry's index, or on a miss
// the index at which compact_field_modulus_cache_store keeps the order.
// The returned pointer is invalidated by the next store.
const CompactFieldModulusCacheEntry* compact_field_modulus_cache_find(
        std::size_t& position,
        CompactFieldModulusCache& cache,
        const NumberField& field,
        flint::FmpzConstRef prime,
        slong threshold) noexcept;

// Records a result for (`prime`, `threshold`) at `position`, which must come
// from a compact_field_modulus_cache_find miss for the same key with no
// store in between.  `modulus` is read only when `usable` is true.
void compact_field_modulus_cache_store(
        CompactFieldModulusCache& cache,
        std::size_t position,
        flint::FmpzConstRef prime,
        slong threshold,
        bool usable,
        const flint::FmpzModPoly& modulus,
        const flint::FmpzModCtx& ctx) noexcept;

bool compact_first_reduction(
        Ideal& reduced,
        FactoredElement& multiplier,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool compact_presentation_reduction(
        Ideal& residual_ideal,
        FactoredElement& multiplier,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        EmbeddingContext* shared_embeddings = nullptr,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
bool compact_presentation_power_root(
        bool& is_power,
        FactoredElement& root,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef decomposition_row,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        EmbeddingContext* shared_embeddings = nullptr,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
// N(remainder_value * scale^n), where remainder_value = prod f^r is the
// evaluation of `remainder` and scale is a rational integer, in a field of
// degree `degree`: prod N(f)^r * scale^(n * degree), by multiplicativity of
// the norm and N(c) = c^degree for c in Q.
// compact_power_root_from_presentation passes it to Element::is_power for the
// expanded element, computing the norm from the factor norms by
// multiplicativity instead of from the expanded element.  Measured compact-path
// norms are small (at most 876 bits on the target rows).  Returns false when a
// factor norm is unavailable or zero.
bool compact_power_input_norm(flint::Fmpq& norm,
                              const FactoredElement& remainder,
                              flint::FmpzConstRef scale,
                              slong n,
                              slong degree) noexcept;
bool compact_unit_presentation_power_root(
        bool& is_power,
        FactoredElement& root,
        const Order& order,
        const FactoredElement& element,
        slong n,
        slong arb_precision,
        slong short_precision,
        const DiagnosticsContext* diagnostics = nullptr,
        CompactFieldModulusCache* field_modulus_cache = nullptr) noexcept;

inline bool valid_certification_request(CertificationMode mode) noexcept {
    return mode == CertificationMode::unknown ||
           mode == CertificationMode::proven;
}

inline bool valid_paired_certification_request(
        CertificationMode mode) noexcept {
    return mode == CertificationMode::grh ||
           mode == CertificationMode::proven;
}

inline bool order_has_quadratic_source_data(const Order& order) noexcept {
    flint::Fmpz radicand;
    return order.parent() != nullptr &&
           order.parent()->backend_kind() == NumberFieldBackendKind::quadratic &&
           order.parent()->quadratic_radicand(flint::FmpzRef(radicand));
}

inline bool order_supports_exact_quadratic_class_certificate(
        const Order& order) noexcept {
    if (order_has_quadratic_source_data(order)) {
        return true;
    }

    flint::Fmpz discriminant;
    return order.parent() != nullptr && order.degree() == 2 &&
           order.discriminant(flint::FmpzRef(discriminant)) &&
           flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) < 0;
}

inline bool prime_to_fractional_ideal(FractionalIdeal& out,
                                      const PrimeIdeal& prime) noexcept {
    if (!same_order_parent(out.parent(), prime.parent())) {
        return false;
    }

    Ideal integral(*out.parent());
    return integral.is_defined() && prime.get_ideal(integral) &&
           out.set_integral(integral);
}

bool unit_candidate_witnesses(std::vector<FactoredElement>& out,
                                    RelationUnitCandidateBatchState& state,
                                    const ClassGroupContext& context,
                                    slong add)
        noexcept;
HnfFinishWorkspace* class_group_finish_workspace(
        ClassGroupContext& context) noexcept;
const HnfFinishWorkspace* class_group_finish_workspace(
        const ClassGroupContext& context) noexcept;
bool class_group_finish_torsion(
        const flint::Fmpz*& torsion_order,
        const Element*& torsion_generator,
        const ClassGroupContext& context) noexcept;
bool class_group_relation_rows(
        flint::FmpzMatRef out,
        const ClassGroupContext& context,
        slong first_relation) noexcept;
bool sync_hnf_finish_workspace(
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept;
bool hnf_finish_workspace_witness_coefficients(
        flint::FmpzMat& out,
        HnfFinishWorkspace& workspace,
        const ClassGroupContext& context) noexcept;
// Kernel basis of the relation-dependency (unit) rows: a basis of the left
// kernel {v : v * relations == 0} of the n x g relation matrix.  The LLL route
// reduces [c * relations | I_n] and takes the identity part of the reduced
// rows whose relation part is zero; it is accepted only when exactly
// n - rank such rows exist and every one satisfies v * relations == 0
// (checked by an exact product).  Otherwise the rows rank..n-1 of the HNF
// transform are used, as before the LLL route existed.  See "Class groups
// and order units" in docs/reference/algorithms_and_sources.rst.
enum class RelationKernelBasisRoute : unsigned char {
    lll_extended,
    hnf_transform,
};

// Binary exponent of the scale c applied to the relation columns.  It does
// not affect correctness: a too-small c only makes the exact checks fail,
// and the HNF transform rows are then used.  The value is not derived from
// the sufficient bound on c in the cited source.
inline constexpr ulong relation_kernel_lll_scale_bits = 40;

// out = the fmpz_lll-reduced n x (g + n) matrix [c * relations | I_n].
bool lll_relation_kernel_extension(flint::FmpzMat& out,
                                   flint::FmpzMatConstRef relations) noexcept;
// Applies the two exact checks to a reduced extension and, when both hold,
// sets out to the (n - rank) x n kernel rows in their reduced order.
// Precondition: `rank` is the exact rank of `relations`; the claim that the
// rows are a basis of the full left kernel lattice depends on it.
bool checked_relation_kernel_rows(flint::FmpzMat& out,
                                  flint::FmpzMatConstRef relations,
                                  flint::FmpzMatConstRef reduced_extension,
                                  slong rank) noexcept;
// The shared selection: the checked rows of `reduced_extension` when it is
// non-null and passes, the HNF transform rows otherwise.  Fails only on
// inconsistent shapes or rank.  Precondition: `rank` is the exact rank of
// `relations`.
bool relation_kernel_basis_from_extension(
        flint::FmpzMat& out,
        RelationKernelBasisRoute* route,
        flint::FmpzMatConstRef relations,
        flint::FmpzMatConstRef hnf_transform,
        slong rank,
        const flint::FmpzMat* reduced_extension) noexcept;
bool relation_kernel_basis(flint::FmpzMat& out,
                           RelationKernelBasisRoute* route,
                           flint::FmpzMatConstRef relations,
                           flint::FmpzMatConstRef hnf_transform,
                           slong rank) noexcept;
// The unit-row coefficients of the HNF analytic finish: relation_kernel_basis
// of the context's stored relations.  hnf_unit_log_matrix takes its rows from
// this function, so the unit logs and the unit witnesses use the same rows.
bool hnf_unit_witness_coefficients(flint::FmpzMat& out,
                                        const ClassGroupContext& context)
        noexcept;
bool hnf_unit_witnesses(std::vector<FactoredElement>& out,
                             const ClassGroupContext& context) noexcept;
bool class_relation_witnessed_hnf_basis(
        WitnessedClassRelationHnfBasis& out,
        const ClassGroupContext& context) noexcept;

// Marks covered[j] = 1 for every column j that is the leading (first nonzero)
// column of some row of the module's HNF row basis, and 0 otherwise. One HNF
// row extraction serves every factor base index of a pivot scan.
bool hnf_row_pivot_column_flags(std::vector<char>& covered,
                                fmpz_smat::HnfContext& module,
                                slong ambient_dim) noexcept;

// Whether `relation` is defined over a factor base equal to the context's.
// Local saturation evaluates this once per residue-dlog kernel and fails every
// nondivisible kernel row when it is false.
bool saturation_kernel_relation_base_matches(
        const Relation& relation,
        const ClassGroupContext& context) noexcept;

}  // namespace silex::detail
