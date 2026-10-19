#pragma once

#include <cstddef>
#include <vector>

#include <silex/class_group.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex {

class EmbeddingContext;
class OrderUnitGroup;

namespace detail {

struct RelationUnitExtractionState;

struct SaturationCandidateSearchResult {
    slong auxiliary_prime_count = 0;
    slong degree_one_prime_count = 0;
    slong quotient_log_row_count = 0;
    slong skipped_quotient_log_row_count = 0;
    slong collapse_after_quotient_log_rows = 0;
    slong final_stable_count = 0;
    slong final_candidate_rows = 0;
    slong final_candidate_cols = 0;
    bool collapsed_to_zero = false;
};

struct SaturationCandidateProcessingResult {
    bool empty_candidate_space = false;
    bool wasted = false;
    bool saw_power = false;
    bool saw_non_power = false;
    bool saw_class_relation = false;
    bool saw_dependent_unit = false;
    slong inspected_columns = 0;
    slong tested_power_candidates = 0;
    slong structural_power_candidates = 0;
    slong compact_power_candidates = 0;
    slong appended_class_relations = 0;
    slong appended_dependent_units = 0;
    slong simplified_relation_count = 0;
    slong simplified_generator_count = 0;
    slong candidate_rows = 0;
    slong candidate_cols = 0;
    slong auxiliary_prime_count = 0;
    SaturationCandidateSearchResult candidate_search;
};

bool relation_basis_extra_indices(std::vector<slong>& basis_indices,
                                        std::vector<slong>& extra_indices,
                                        const ClassGroupContext& context)
        noexcept;
bool saturation_relation_data(std::vector<FactoredElement>& relations,
                                    flint::FmpzMatRef relation_rows,
                                    const ClassGroupContext& context) noexcept;
bool simplified_saturation_relation_data(
        std::vector<FactoredElement>& relations,
        flint::FmpzMat& relation_rows,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        slong n,
        slong precision) noexcept;
bool simplified_saturation_context(
        ClassGroupContext& simplified,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        const RelationUnitExtractionState& unit_state,
        slong n,
        slong precision) noexcept;
bool saturation_candidates(
        flint::FmpzMat& candidates,
        slong& auxiliary_prime_count,
        const ClassGroupContext& context,
        slong n,
        double stable,
        SaturationCandidateSearchResult* search_result = nullptr)
        noexcept;
bool saturation_process_candidates_once(
        SaturationCandidateProcessingResult& result,
        ClassGroupContext& context,
        flint::FmpzMatConstRef candidates,
        slong n,
        bool append_class_relations,
        OrderUnitGroup* units = nullptr,
        EmbeddingContext* embeddings = nullptr,
        RelationUnitExtractionState* unit_state = nullptr,
        slong precision = 0,
        slong max_class_relation_appends = WORD_MAX) noexcept;
bool saturation_process_candidates_once(
        SaturationCandidateProcessingResult& result,
        const ClassGroupContext& input_context,
        ClassGroupContext& append_context,
        flint::FmpzMatConstRef candidates,
        slong n,
        bool append_class_relations,
        OrderUnitGroup* units = nullptr,
        EmbeddingContext* embeddings = nullptr,
        RelationUnitExtractionState* unit_state = nullptr,
        slong precision = 0,
        slong max_class_relation_appends = WORD_MAX) noexcept;
// The candidate search of saturation_candidates from the given relations,
// starting at auxiliary_prime_start and using at most max_auxiliary_primes
// auxiliary primes; exposed so tests can reach a chosen auxiliary prime.
bool saturation_candidates_from_relations_from_prime(
        flint::FmpzMat& candidates,
        slong& auxiliary_prime_count,
        const std::vector<FactoredElement>& relations,
        const Order& order,
        slong n,
        double stable,
        flint::FmpzConstRef auxiliary_prime_start,
        slong max_auxiliary_primes,
        SaturationCandidateSearchResult* search_result = nullptr) noexcept;
// The degree-one quotient-log row of saturation_candidates_from_relations at
// (q, theta - root), from the exact fmpz evaluation and from the word-size
// path used for q < 2^64; exposed so tests can compare them.
bool degree_one_quotient_log_row_general(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& relations,
        flint::FmpzConstRef q,
        flint::FmpzConstRef root,
        slong n) noexcept;
bool degree_one_quotient_log_row_word_size(
        flint::FmpzMat& out,
        const std::vector<FactoredElement>& relations,
        flint::FmpzConstRef q,
        flint::FmpzConstRef root,
        slong n) noexcept;
// The factors of the exponent-reduced inputs of the degree-one residue
// characters, as rational polynomials, copied once so that each auxiliary
// prime reads them in place.
struct DegreeOneResidueFactor {
    flint::FmpqPoly polynomial;
    slong exponent = 0;
};

using DegreeOneResidueInput = std::vector<DegreeOneResidueFactor>;

// Class-relation saturation with units: the candidate filter from the
// stacked residue characters of reference
// `Saturate.jl:compute_candidates_for_saturate`, over the relation
// generators, free units and (when include_torsion) the torsion generator.
// The stack grows one degree-one prime at a time, only when a candidate that
// the current stack does not rule out is about to take its power test, and
// stops by the reference rule with `stable` starting at 1.5 and doubling
// after a failed power test.  It starts only after the round has run
// `width` power tests (`power_tests`), a Silex count gate that the reference
// does not have.  See class_group.cpp.
struct RelationSaturationCandidateFilter {
    bool started = false;
    // False when the filter does not apply or failed; every candidate is
    // then power tested.
    bool available = false;
    // The reference stopping rule holds: the kernel dimension has stayed the
    // same for more than `stable * width` degree-one primes.
    bool stopped = false;
    slong ell = 0;
    slong width = 0;
    double stable = 1.5;
    slong stable_count = 1;
    slong stacked_prime_count = 0;
    // Power tests run in this round, counted before the filter starts too.
    slong power_tests = 0;
    const Order* order = nullptr;
    const DiagnosticsContext* diagnostics = nullptr;
    std::vector<FactoredElement> inputs;
    std::vector<DegreeOneResidueInput> prepared_inputs;
    bool prepared = false;
    flint::Fmpz discriminant;
    flint::Fmpz equation_index;
    bool equation_index_known = false;
    flint::Fmpz q;
    bool q_loaded = false;
    std::vector<flint::Fmpz> roots;
    std::size_t next_root = 0;
    // The stacked characters modulo ell as independent rows, each reduced
    // against the earlier ones (pivot entry 1); the kernel dimension is
    // width - rank.
    std::vector<std::vector<ulong>> rows;
    std::vector<slong> pivots;
    // Round state outside the filter proper (kept by start, which resets
    // only the filter): the kernel rows that have already gone through
    // saturation_relation_with_units in this round without an appended
    // relation, sorted.  The per-prime kernels of a round can repeat a row,
    // and its outcome cannot change before the round ends, so it is not
    // tested again (relation_saturation_round_row_settled).
    std::vector<std::vector<ulong>> settled_rows;
};

bool relation_saturation_candidate_filter_start(
        RelationSaturationCandidateFilter& filter,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        bool include_torsion,
        flint::FmpzConstRef ell) noexcept;
bool relation_saturation_candidate_filter_grow(
        RelationSaturationCandidateFilter& filter) noexcept;
bool relation_saturation_candidate_filter_rules_out(
        bool& ruled_out,
        const RelationSaturationCandidateFilter& filter,
        flint::FmpzMatConstRef kernel,
        slong row) noexcept;
bool relation_saturation_candidate_filter_kernel_empty(
        const RelationSaturationCandidateFilter& filter) noexcept;
void relation_saturation_candidate_filter_failed_power_test(
        RelationSaturationCandidateFilter& filter) noexcept;
// Whether row `row` of `kernel` was already settled in this round (see
// RelationSaturationCandidateFilter::settled_rows); rows with an entry
// outside [0, 2^64) are never settled.
bool relation_saturation_round_row_settled(
        const RelationSaturationCandidateFilter& filter,
        flint::FmpzMatConstRef kernel,
        slong row) noexcept;
// Records row `row` of `kernel` as settled in this round.
void relation_saturation_round_settle_row(
        RelationSaturationCandidateFilter& filter,
        flint::FmpzMatConstRef kernel,
        slong row) noexcept;
void relation_saturation_candidate_filter_screen(
        bool& ruled_out,
        RelationSaturationCandidateFilter& filter,
        const ClassGroupContext& context,
        const OrderUnitGroup& units,
        bool include_torsion,
        flint::FmpzConstRef ell,
        flint::FmpzMatConstRef kernel,
        slong row) noexcept;
// Integral Gram-Schmidt data of the rows b_0, ..., b_{m-1} of `basis`
// (Cohen, GTM 138, Algorithm 2.6.7, step 2): `d` becomes 1 x (m + 1) with
// d_0 = 1 and d_{j+1} = d_j |b*_j|^2, and `lambda` becomes m x m with
// lambda(k, j) = d_{j+1} mu(k, j) for j < k (zero elsewhere).  Fails when
// the rows are linearly dependent.
bool integral_gram_schmidt(flint::FmpzMat& d,
                           flint::FmpzMat& lambda,
                           flint::FmpzMatConstRef basis) noexcept;

// Size-reduces the 1 x n `row` against the rows of `basis`, whose integral
// Gram-Schmidt data `d` and `lambda` come from integral_gram_schmidt: for
// l = m-1, ..., 0 it applies Cohen's sub-algorithm RED (Algorithm 2.6.7),
// subtracting round(mu(row, l)) b_l, so that afterwards every
// |mu(row, l)| <= 1/2 (Babai's nearest-plane step).  Only integer
// multiples of basis rows are subtracted, so `row` changes by a vector of
// the lattice the basis spans.
bool size_reduce_row_against_basis(flint::FmpzMatRef row,
                                   flint::FmpzMatConstRef basis,
                                   flint::FmpzMatConstRef d,
                                   flint::FmpzMatConstRef lambda) noexcept;

// For tests and audits: a copy of the relation-kernel basis that
// relation_saturation_beta_coefficients size-reduces against, building it
// if needed, and whether the LLL-reduced kernel basis that
// relation_kernel_unit reads is built at return.
bool relation_saturation_kernel_basis(flint::FmpzMat& out,
                                      bool& lll_basis_built,
                                      const ClassGroupContext& context) noexcept;

// For tests and audits: a copy of the integral Gram-Schmidt data (see
// integral_gram_schmidt) that relation_saturation_beta_coefficients uses
// with that basis, building the basis and the data if needed.
bool relation_saturation_kernel_gram(flint::FmpzMat& d,
                                     flint::FmpzMat& lambda,
                                     const ClassGroupContext& context) noexcept;

class TestHookAccess;

// A failure that only detail::TestHookAccess (test/test_hook_access.hpp)
// can inject into the extension route of the relation-saturation kernel
// basis (ClassGroupRelationAccess::saturation_relation_kernel); the
// certified first phase has CertifiedLllTestSeams
// (lat/certified_lll_internal.hpp).  The route's output is zeroed after the
// reduction runs, so the fallback is exercised after the reduction has
// changed the matrix.  first_phase: the reduced extension is zeroed, so its
// kernel-row count check fails.  certificate: fmpz_mat_is_reduced runs on
// the zeroed basis, which it does not certify.  beta_solve: the kernel
// route is untouched, and the lattice basis that
// ClassGroupRelationAccess::solved_invariant_generator_relation_rows solves
// against is zeroed, so the solve fails and the beta coefficients are read
// from the left transform.
enum class SaturationKernelTestFailure {
    none,
    first_phase,
    certificate,
    beta_solve
};

// The size crossover of ClassGroupRelationAccess::saturation_relation_kernel.
// A kernel of more than one row takes the extension route first only when
// the presentation has more generators (columns of its relation matrix)
// than this; otherwise it is read from the presentation's left transform
// first, as before the extension route existed.  Chosen by measurement, not
// from a source: on the proven corpus (busy-machine timings), the
// left-transform route built the
// kernel faster on every presentation with at most 31 generators
// (real_quadratic_d100000001, 23 generators: 20 ms against 78 ms, about
// 0.08 s against 0.13 s per run), the two routes were level at 37
// (real_quadratic_d100000005), and from 57 generators the HNF transform of
// the left transform dominated (real_quadratic_d10000000005, 64
// generators: 3.8 s against 0.11 s).  Both routes return a kernel basis,
// so the cutoff changes the cost only.
inline constexpr slong kSaturationKernelTransformMaxGenerators = 32;

// Read by saturation_relation_kernel each time it builds a basis.  Off by
// default and never set by production code.  Relaxed atomics so a reader
// never races with a test setting them.
class SaturationKernelTestSeams {
public:
    static SaturationKernelTestFailure failure() noexcept {
        return failure_.load(std::memory_order_relaxed);
    }
    // When set, a kernel of more than one row is read from the
    // presentation's left transform first whatever the presentation's size,
    // as a kernel of one row is, so tests can compare the two routes and
    // reach the first phase and its fallbacks on multi-row kernels.
    // Overrides force_extension_route.
    static bool force_transform_route() noexcept {
        return force_transform_route_.load(std::memory_order_relaxed);
    }
    // When set, a kernel of more than one row takes the extension route
    // first whatever the presentation's size, so tests reach that route on
    // presentations below kSaturationKernelTransformMaxGenerators.
    static bool force_extension_route() noexcept {
        return force_extension_route_.load(std::memory_order_relaxed);
    }
    // Counts each extension-route run that an injected
    // failure replaced.
    static void note_injected() noexcept {
        injected_.fetch_add(1, std::memory_order_relaxed);
    }
    // Counts each extension-route basis that saturation_relation_kernel
    // kept.
    static void note_extension_kept() noexcept {
        extension_kept_.fetch_add(1, std::memory_order_relaxed);
    }
    // Counts each set of beta coefficient rows that
    // relation_saturation_beta_coefficients solved from the extension
    // route's image rows instead of reading the left transform.
    static void note_beta_solved() noexcept {
        beta_solved_.fetch_add(1, std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<SaturationKernelTestFailure> failure_{
            SaturationKernelTestFailure::none};
    inline static std::atomic<bool> force_transform_route_{false};
    inline static std::atomic<bool> force_extension_route_{false};
    inline static std::atomic<slong> injected_{0};
    inline static std::atomic<slong> extension_kept_{0};
    inline static std::atomic<slong> beta_solved_{0};
};

// The relation coefficients of the ell-local relation-saturation beta
// rows: one row of `rows` per invariant index i with ell | d_i, in index
// order: a coefficient row c with c * A equal to the exponent row of
// g_i^{d_i}, size-reduced against an LLL-reduced basis of the relation
// kernel kept for this purpose (see
// ClassGroupRelationAccess::saturation_relation_kernel).  The row c is the
// integer solution read off the reduced [c * A | I] when the kernel basis
// came from that extension route
// (ClassGroupRelationAccess::solved_invariant_generator_relation_rows), and
// the presentation's row (invariant_generator_relation_matrix) otherwise;
// either differs from the presentation's row by a relation-kernel vector.
// Requires `units` to be proven for the same order; see the definition for
// the soundness condition.
bool relation_saturation_beta_coefficients(flint::FmpzMat& rows,
                                           const ClassGroupContext& context,
                                           const OrderUnitGroup& units,
                                           flint::FmpzConstRef ell) noexcept;

// The beta rows of the ell-local relation-saturation test (2.17.3
// src/basemath/buch3.c check_prime, cited in
// docs/reference/algorithms_and_sources.rst): for each invariant index i with ell | d_i, a compact
// witness beta_i with (beta_i) = g_i^{d_i}.  The relation coefficients of
// each beta_i are size-reduced against an LLL-reduced relation-kernel basis, so
// beta_i differs from the presentation's witness
// (invariant_generator_power_witness) by a relation-kernel unit.  Requires
// `units` to be proven for the same order; see the definition.
bool relation_saturation_beta_rows(std::vector<FactoredElement>& beta_rows,
                                   const ClassGroupContext& context,
                                   const OrderUnitGroup& units,
                                   flint::FmpzConstRef ell) noexcept;

bool saturate_class_unit_context(
        bool& success,
        SaturationCandidateProcessingResult& last_result,
        ClassGroupContext& context,
        OrderUnitGroup& units,
        EmbeddingContext& embeddings,
        RelationUnitExtractionState& unit_state,
        slong n,
        double stable,
        slong precision) noexcept;

}  // namespace detail
}  // namespace silex
