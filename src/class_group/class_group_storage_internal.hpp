#pragma once

#include <silex/class_group.hpp>
#include <silex/element.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/nmod_poly.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_field.hpp>

#include "ideal_lattice_lll_internal.hpp"
#include "../abelian_group/abelian_group_internal.hpp"
#include "cyclotomic_galois_images_internal.hpp"
#include "relation_completion_scheduler_internal.hpp"
#include "relation_search_budget_internal.hpp"
#include "../factor_base/generation_different_route_internal.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace silex::detail {

// Smallest field degree at which the relation search uses its large-field
// rules: one random-walk ideal serves every pivot of a pass
// (lll_relation_search.cpp, collect_pivot_lll_ideal), and a grh search on a
// field defined by exactly Phi_n or psi_n adds the Galois images of its
// relations (relation_candidate.cpp, the Galois-image step).  Below it the
// search draws one ideal per pivot, and only a field defined by exactly
// Phi_n that passes capped_galois_images_selected adds images, under the
// capped rule of that step.  The crossover is measured, not sourced (see
// docs/reference/algorithms_and_sources.rst).
inline constexpr slong kLargeFieldRelationSearchMinDegree = 30;

// Below kLargeFieldRelationSearchMinDegree, a grh search on a field defined
// by exactly Phi_n adds capped Galois images only when the field degree is
// at most kCappedGaloisImageMaxDegree or the factor-base bound is at least
// kCappedGaloisImageMinFactorBaseBound.  Both limits are Silex's own, set
// by measurement, not sourced: on the corpus Phi_n rows below degree 30 the
// images cost more relation work than they save on every degree 22 and 24
// field with a factor-base bound below 300; two smaller rows that remain
// slower (Q(zeta_7), Q(zeta_60)) are accepted.
inline constexpr slong kCappedGaloisImageMaxDegree = 20;
inline constexpr ulong kCappedGaloisImageMinFactorBaseBound = 300;

inline bool capped_galois_images_selected(
        slong degree,
        flint::FmpzConstRef factor_base_bound) noexcept {
    return degree <= kCappedGaloisImageMaxDegree ||
           flint::fmpz_cmp_ui(factor_base_bound,
                              kCappedGaloisImageMinFactorBaseBound) >= 0;
}

struct FactorBaseHonestyFailure;

struct PartialRelationEntry {
    explicit PartialRelationEntry(const Order& order) noexcept
        : prime(order),
          generator(*order.parent()) {
    }

    // The large prime ideal, when built.  A deferred entry has none yet:
    // it records the order coordinates of `generator` and the norm
    // valuation at rational_prime, from which
    // order_element_deferred_large_prime builds it on the first stored or
    // incoming partial above the same rational prime.
    PrimeIdeal prime;
    Element generator;
    flint::FmpzMat row{0, 0};
    // The rational prime under the large prime, for every entry.
    flint::Fmpz rational_prime;
    flint::FmpzPoly residue_key;
    flint::FmpzMat deferred_coordinates{0, 0};
    slong deferred_norm_valuation = 0;
    bool has_residue_key = false;
    bool has_row = false;
    bool has_deferred_prime = false;
};

struct DirectPartialResidueBlockCacheEntry {
    DirectPartialResidueBlockCacheEntry(
            slong block_index,
            flint::FmpzConstRef rational_prime,
            bool allow_word_size = true) noexcept
        : context(rational_prime.raw()),
          input(context),
          remainder(context),
          word_size(allow_word_size &&
                    rational_prime_is_word_size(rational_prime)),
          word_input(word_modulus(word_size, rational_prime)),
          word_remainder(word_modulus(word_size, rational_prime)),
          block_index(block_index) {
        flint::fmpz_set(flint::FmpzRef(this->rational_prime),
                        rational_prime);
    }

    static bool rational_prime_is_word_size(
            flint::FmpzConstRef rational_prime) noexcept {
        return ::fmpz_sgn(rational_prime.raw()) > 0 &&
               ::fmpz_abs_fits_ui(rational_prime.raw()) != 0;
    }

    static ulong word_modulus(bool word_size,
                              flint::FmpzConstRef rational_prime) noexcept {
        return word_size ? ::fmpz_get_ui(rational_prime.raw()) : 2;
    }

    // Residue polynomials stored for this block: the word-size (nmod)
    // copies when `word_size` is set, otherwise the `fmpz_mod` copies.
    std::size_t residue_count() const noexcept {
        return word_size ? word_residue_polynomials.size()
                         : residue_polynomials.size();
    }

    flint::Fmpz rational_prime;
    flint::FmpzModCtx context;
    flint::FmpzModPoly input;
    flint::FmpzModPoly remainder;
    std::vector<flint::FmpzModPoly> residue_polynomials;
    // Same residue test over Z/pZ with single-limb arithmetic when the
    // rational prime fits in a word; the `fmpz_mod` members are unused
    // then.  reference FacBase_Idl.jl likewise builds the word-size
    // `FactorBaseSingleP` when `fits(Int, p)` and the big-integer one
    // otherwise.
    bool word_size = false;
    flint::NmodPoly word_input;
    flint::NmodPoly word_remainder;
    std::vector<flint::NmodPoly> word_residue_polynomials;
    std::vector<slong> prime_indices;
    std::vector<slong> residue_degrees;
    // True when the block holds every prime ideal above its rational prime
    // p: the sum of e(P/p) f(P/p) over the block equals the order degree.
    // The factor base keeps distinct primes, one block per rational prime.
    bool covers_all_primes_over_p = false;
    slong block_index = -1;
};

class TestHookAccess;

// Test entry of the direct large-prime partial relation; only
// detail::TestHookAccess (test/test_hook_access.hpp) can call it.
class DirectPartialRelationTestEntry {
    friend class TestHookAccess;

    // Test hook: runs the direct large-prime partial-relation factorization
    // (the residue-polynomial branch of reference Rel_add.jl `_factor!`) on
    // `generator` over `base` with a fresh residue-block cache.  Returns true
    // and writes the factor-base row and leftover rational prime when that
    // branch finds a single large prime; returns false otherwise.  When given,
    // `known_norm` must be the exact field norm of `generator` and
    // `known_integral_coordinates` its maximal-order coordinates.  With
    // `word_size_residues` false the residue test uses the multi-limb
    // `fmpz_mod` path even for word-size rational primes.  With
    // `deduce_last_hit_valuation` false every residue-hit valuation is computed,
    // none is deduced from the norm; when `deduced_hit_valuations` is given it
    // receives the number of valuations deduced.  For tests only.
    static bool row(
            flint::FmpzMat& out_row,
            flint::Fmpz& out_rational_large_prime,
            const Element& generator,
            const FactorBase& base,
            flint::FmpzConstRef factor_base_bound,
            const fmpq* known_norm = nullptr,
            const fmpz_mat_struct* known_integral_coordinates = nullptr,
            bool word_size_residues = true,
            bool deduce_last_hit_valuation = true,
            slong* deduced_hit_valuations = nullptr) noexcept;
};

struct FactorBaseGenerationRecord {
    flint::Fmpz p;
    ProofState status = ProofState::not_checked;
};

struct RelationSaturationRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
};

// What a relation-saturation proof record at `ell` rests on.
//
// - `ell_local_dlog`: the ell-local test of
//   ClassGroupContext::prove_relation_saturation_dlog_ell_; `rank`,
//   `target_rank`, and `local_primes` describe that computation.  When the
//   local target rank is 0 the proof is vacuous and no computation ran;
//   the record then stores zero for all three.
// - `exact_class_order`: the exact imaginary-quadratic route.  With
//   generation verified and h_cand equal to the exact class number h, the
//   computed relation lattice has index h_cand / h = 1 in the full one, so it
//   is saturated at every prime; `rank`, `target_rank`, and `local_primes`
//   are zero because no ell-local computation ran.
enum class RelationSaturationProofBasis : unsigned char {
    ell_local_dlog,
    exact_class_order,
};

struct RelationSaturationProofRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
    slong rank = 0;
    slong target_rank = 0;
    slong local_primes = 0;
    RelationSaturationProofBasis basis =
            RelationSaturationProofBasis::ell_local_dlog;
};

struct RelationAdmissionCache {
    slong target_relation_count = 0;
    slong relation_surplus = 0;
    slong relation_count_before_init = 0;
    slong retained_relation_count = 0;
    slong trivial_relation_count = 0;
    slong duplicate_relation_count = 0;
    slong skipped_dependent_relation_count = 0;
    slong missing = 0;
    std::vector<ulong> modular_basis;
    std::vector<ulong> relation_row_hashes;
};

struct RelationSearchControlState {
    RelationCompletionState completion;
    slong checkpoint_relation_count = -1;
    slong done_small = 0;
    slong small_fail = 0;
    slong fail_limit = 0;
    slong candidates_tried = 0;
    slong accepted_relations = 0;
    slong checkpoint_extra_relation_requests = 0;
    slong analytic_extra_relation_requests = 0;
    slong analytic_finish_precision = 160;
    slong analytic_finish_product_precision = 0;
    slong analytic_precision_doublings = 0;
    slong analytic_precision_inconclusive_count = 0;
    slong bounded_full_factor_base_restart_requests = 0;
    slong honesty_restart_requests = 0;
    slong squash_index = 0;
    slong finish_unit_log_rotation = 0;
    bool bounded_full_factor_base_restart_pending = false;
    bool bounded_full_factor_base_restart_allow_past_half = false;
    bool finish_unit_log_rotation_active = false;
    bool finish_full_rank_relation_active = false;
    bool honesty_checked = false;
    bool analytic_finish_product_valid = false;
    // reference relation_completion_parameters retains nfrootsof1 from compute_inverse_residue through
    // fundamental-unit publication. Keep the same private finish lifetime.
    bool analytic_finish_torsion_valid = false;
    flint::Arb analytic_finish_product;
    flint::Fmpz analytic_finish_torsion_order;
    Element analytic_finish_torsion_generator;
    // The analytic finish check runs once per relation checkpoint, always
    // at analytic_finish_precision.  Its embedding roots depend only on the
    // field and that precision, so the context built for a precision is
    // kept for the next check at the same precision.  It is rebuilt when
    // the precision changes, so every check sees the roots a fresh context
    // would compute.  It is cleared, with its log cache, when the check
    // returns finished and before the route state is stored as the live
    // finish state.
    EmbeddingContext analytic_finish_embeddings;
    slong analytic_finish_embedding_precision = 0;
};

struct RelationFinishState {
    RelationAdmissionCache cache;
    SubfactorBaseSchedule subfactor_base_schedule;
    RelationSearchControlState route_state;
    std::vector<slong> small_multiplier;
    flint::Fmpz factor_base_bound;
    slong target_relation_kernel_units = 0;
};

// The finish workspace of the processed relations.  It keeps no HNF: a
// sync that changes the relations checks that fmpz_mat_rank of all stored
// relations is the context's relation rank and records it as relation_rank
// (sync_hnf_finish_workspace).
struct HnfFinishWorkspace {
    flint::ArbMat relation_logs{0, 0};
    flint::ArbMat unit_logs{0, 0};
    // Kernel rows of the processed relations, kept until the next sync
    // changes the relations: hnf_unit_finish_coefficients for the analytic
    // finish and hnf_unit_witness_coefficients for the witness extraction.
    // unit_logs holds the logs of the finish rows when
    // unit_logs_finish_kernel is set and of the witness rows otherwise.
    flint::FmpzMat finish_kernel{0, 0};
    flint::FmpzMat witness_kernel{0, 0};
    // The reduced [c * relations | I] matrix of the first
    // finish_extension_relation_count relations that the last finish kernel
    // was read from (see hnf_finish_workspace_kernel).  Unlike the kernel
    // caches it is kept across syncs that append relations: the next finish
    // kernel starts its reduction from it.  The kernel claim for rows read
    // from it rests on its identity part being unimodular, which the exact
    // checks do not detect: only the phase may write it
    // (finish_extension_kernel_rows, from the phase output of an input with
    // a unimodular identity part).  The kernel_rows test seam
    // (FinishKernelTestFailure) violates this when rank is 0; tests use rank >= 1.  Relations are only ever appended, and a
    // relation part that does not match them only makes the checks fail.
    flint::FmpzMat finish_extension{0, 0};
    slong finish_extension_relation_count = 0;
    slong generator_count = 0;
    slong processed_relation_count = 0;
    slong relation_rank = 0;
    slong relation_log_precision = 0;
    slong logged_relation_count = 0;
    slong unit_log_precision = 0;
    bool exact_valid = false;
    bool relation_logs_valid = false;
    bool unit_logs_valid = false;
    bool unit_logs_finish_kernel = false;
    bool finish_kernel_valid = false;
    bool witness_kernel_valid = false;
    bool finish_extension_valid = false;

    void clear() noexcept {
        relation_logs = flint::ArbMat(0, 0);
        unit_logs = flint::ArbMat(0, 0);
        finish_kernel = flint::FmpzMat(0, 0);
        witness_kernel = flint::FmpzMat(0, 0);
        finish_extension = flint::FmpzMat(0, 0);
        finish_extension_relation_count = 0;
        generator_count = 0;
        processed_relation_count = 0;
        relation_rank = 0;
        relation_log_precision = 0;
        logged_relation_count = 0;
        unit_log_precision = 0;
        exact_valid = false;
        relation_logs_valid = false;
        unit_logs_valid = false;
        unit_logs_finish_kernel = false;
        finish_kernel_valid = false;
        witness_kernel_valid = false;
        finish_extension_valid = false;
    }
};

struct RandomIdealSearchState {
    std::vector<slong> base_indices;
    std::vector<slong> exponents;
    std::vector<Ideal> base_ideals;
    std::vector<Ideal> inverse_base_numerators;
    std::vector<flint::Fmpz> inverse_base_denominators;
    Ideal random_ideal;
    flint::Fmpz lower_bound;
    flint::Fmpz upper_bound;
    ulong random_state = UWORD(0);
    bool defined = false;
};

struct ReducedIdealLatticeCacheEntry {
    ReducedIdealLatticeCacheEntry(flint::FmpzMatConstRef ideal_hnf,
                                  flint::FmpzMatConstRef reduced_basis)
        : ideal_hnf(flint::fmpz_mat_nrows(ideal_hnf),
                    flint::fmpz_mat_ncols(ideal_hnf)),
          reduced_basis(flint::fmpz_mat_nrows(reduced_basis),
                        flint::fmpz_mat_ncols(reduced_basis)) {
        flint::fmpz_mat_set(this->ideal_hnf, ideal_hnf);
        flint::fmpz_mat_set(this->reduced_basis, reduced_basis);
    }

    flint::FmpzMat ideal_hnf;
    flint::FmpzMat reduced_basis;
};

struct ReducedIdealLatticeCache {
    // Match the existing source-backed post-finite total candidate budget.
    static constexpr std::size_t capacity = 64;

    std::vector<ReducedIdealLatticeCacheEntry> entries;
    std::size_t next_eviction = 0;

    void clear() noexcept {
        entries.clear();
        next_eviction = 0;
    }
};

// Residue-field quotient logarithms reused by local class-relation
// saturation with units (`saturate_local_once_with_units_`).  Every value is
// the exact `ResidueFieldQuotientLog::apply` image of one number-field
// element at one (prime, ell), so it stays valid while that element is
// unchanged.  The saturation pass checks the relation generators and unit
// factorizations against `elements` by value on every use and clears the
// cache on any mismatch; the cache never decides which primes or kernel rows
// are tried.
struct LocalSaturationDlogEntry {
    enum class State : unsigned char {
        unknown,
        unit,          // nonzero residue; `logs[i]` holds its quotient log
        zero,          // the residue is zero
        unreducible,   // the element does not reduce modulo the prime
        failed,        // the quotient log of a nonzero residue failed
    };

    flint::Fmpz ell;
    ResidueFieldQuotientLog quotient_log;
    ResidueFieldElement image;
    // Indexed like LocalSaturationDlogCache::elements.
    std::vector<State> states;
    std::vector<flint::Fmpz> logs;
};

struct LocalSaturationDlogCache {
    struct UnitTerm {
        slong element = 0;
        slong exponent = 0;
    };

    // Distinct relation generators and unit-factor bases.
    std::vector<Element> elements;
    std::vector<std::uint64_t> element_hashes;
    // Relation generator i is elements[relation_elements[i]].
    std::vector<slong> relation_elements;
    // Free unit generator u is the product of elements[t.element]^t.exponent
    // over unit_terms[u], factor by factor in its stored order.
    std::vector<std::vector<UnitTerm>> unit_terms;
    bool units_recorded = false;
    std::vector<LocalSaturationDlogEntry> entries;

    void clear() noexcept {
        elements.clear();
        element_hashes.clear();
        relation_elements.clear();
        unit_terms.clear();
        units_recorded = false;
        entries.clear();
    }
};

struct ClassGroupContextStorage {
    std::vector<PartialRelationEntry> partial_relations;
    std::vector<DirectPartialResidueBlockCacheEntry>
            direct_partial_residue_blocks;
    slong partial_relations_max = 0;
    slong partial_relations_nonproductive_streak = 0;
    bool use_partial_relations = false;
    bool partial_relations_configured = false;

    RandomIdealSearchState random_ideal_search_state;
    // The cap that ended the latest LLL relation-route slice of this context
    // short of its goal (none when no cap did, or when the slice failed for
    // another reason).  Written by run_lll_relation_route_ and the
    // tentative-candidate and extension calls around it; read only to name
    // the failure of a class/unit transaction.
    RelationSearchBudgetCap lll_route_budget_cap =
            RelationSearchBudgetCap::none;
    ReducedIdealLatticeCache native_post_finite_reduced_lattices;
    // Per-order data of the rounded-embedding ideal reduction (rounded order
    // rows, order T2 Gram, and the result of its route check), shared by
    // every run_lll_relation_route_ call of this context, so that the route
    // check runs once per order and not once per relation-route slice.
    relation_search::IdealLatticeRoundedEmbeddingCache
            lll_rounded_embedding_cache;
    // Row capture of the small-LLL torsion-orbit skip.  While armed (non-null)
    // the first generator relation that factors over the base copies its
    // exponent row here and disarms the capture; the candidate loop arms it
    // around one candidate test and disarms it after.
    flint::FmpzMat* generator_relation_row_capture = nullptr;
    bool generator_relation_row_captured = false;
    std::optional<RelationFinishState> finish_state;
    HnfFinishWorkspace hnf_finish_workspace;
    slong analytic_finish_precision = 0;
    flint::Arb analytic_finish_product;
    slong analytic_finish_product_precision = 0;
    bool analytic_finish_product_valid = false;

    std::vector<FactorBaseGenerationRecord> factor_base_generation_records;
    // GRH-conditional factor-base generation (grh class/unit requests only):
    // set when the factor base was built to contain every prime ideal of
    // norm at most `grh_generation_bound`, the GRH generation bound
    // min(BDF, Bach).  Kept apart from the unconditional generation status,
    // which it never sets.  Written only when the base is built, like
    // factor_base_generation_status_, so certification transactions do not
    // restore it.
    flint::Fmpz grh_generation_bound;
    bool grh_generation_covered = false;
    // The theorem that gave grh_generation_bound (bdf or bach).
    ClassGroupGenerationBasis grh_generation_basis =
            ClassGroupGenerationBasis::none;
    // How the unconditional generation requirement (the recorded
    // factor_base_generation_bound_, from proven_generation_requirement)
    // treats the class of the different: none for the
    // factor_base_class_group_bound requirement, kor1_generation for
    // Zimmert's twin-class requirement.
    // Written with that bound when the base is built, so certification
    // transactions do not restore it.
    GenerationDifferentRoute generation_different_route =
            GenerationDifferentRoute::none;
    // Relation and unit completeness of a `grh`-requested pair, written
    // when the pair is accepted (record_grh_acceptance_analytic_check,
    // record_grh_acceptance_exact_class_number) or relabelled
    // (publish_grh_labels).  Read only while the published label is `grh`;
    // cleared with the other certification metadata and restored by
    // certification transactions.
    CertificationMode class_unit_completeness = CertificationMode::unknown;
    ClassUnitCompletenessBasis class_unit_completeness_basis =
            ClassUnitCompletenessBasis::none;
    // Set when the published imaginary quadratic class group was certified
    // `proven` by its structure certificate, which proves that the factor
    // base generates the class group without a generation bound; the
    // generation status stays as built.  Read by
    // factor_base_generation_certification() (basis structure_certificate).
    // It describes the published presentation, so it is cleared with the
    // other certification metadata and restored by certification
    // transactions.
    bool generation_structure_certified = false;
    std::vector<RelationSaturationRecord> relation_saturation_records;
    std::vector<RelationSaturationProofRecord>
            relation_saturation_proof_records;
    LocalSaturationDlogCache local_saturation_dlog_cache;
    // Relation buffer of the local saturation passes (saturate_local_once_
    // and saturate_local_once_with_units_), kept across passes so its factor
    // base copy and storage are not rebuilt on every pass.  Each pass reuses
    // it only while its factor base equals the context base and redefines
    // it otherwise; every candidate row overwrites its generator and
    // exponents before the relation is read, so no value carries over.
    Relation saturation_relation_scratch;
    // LLL-reduced basis of the relation kernel of the published
    // presentation, read by ClassGroupContext::relation_kernel_unit.
    // Rebuilt on first use after the presentation changes; every write of
    // the presentation resets the flag.
    flint::FmpzMat reduced_relation_kernel{0, 0};
    bool reduced_relation_kernel_valid = false;
    // A basis of the same relation kernel, read only by the
    // relation-saturation beta size reduction
    // (relation_saturation_beta_coefficients) and never by
    // relation_kernel_unit.  It is a copy of reduced_relation_kernel when
    // that is already built.  Otherwise, for at least two rows over a
    // presentation with more than kSaturationKernelTransformMaxGenerators
    // generators, it is the checked kernel rows of the reduced [c * A | I]
    // certified LLL-reduced (the extension route).  At or below that
    // crossover, or when that route fails, it is read from the left
    // transform and reduced by FLINT's
    // double-precision reduction pass certified exactly LLL-reduced
    // (lat::detail::lll_first_phase_certified), or a copy of the fmpz_lll
    // basis when that certificate is not obtained
    // (ClassGroupRelationAccess::saturation_relation_kernel).  Built on first
    // use and reset together with the reduced kernel.
    flint::FmpzMat saturation_relation_kernel{0, 0};
    bool saturation_relation_kernel_valid = false;
    // When saturation_relation_kernel came from the extension route, the
    // identity parts of the other rows of the same reduced [c * A | I]: a
    // g x n matrix Q whose product Q * A is a basis of the row lattice of
    // the presentation's relation matrix A.  The relation-saturation beta
    // coefficients are then solved from it
    // (ClassGroupRelationAccess::solved_invariant_generator_relation_rows)
    // instead of read from the left transform.  Set together with the
    // kernel and reset together with it.
    flint::FmpzMat saturation_relation_image_transform{0, 0};
    bool saturation_relation_image_valid = false;
    // Integral Gram-Schmidt data of saturation_relation_kernel (Cohen,
    // GTM 138, Algorithm 2.6.7): a 1 x (m + 1) row d with d_0 = 1 and
    // d_{j+1} = d_j |b*_j|^2, and an m x m matrix lambda with
    // lambda(k, j) = d_{j+1} mu(k, j) for j < k.  Read by the
    // relation-saturation beta size reduction; built on first use and reset
    // together with the kernel bases.
    flint::FmpzMat saturation_relation_kernel_gram_d{0, 0};
    flint::FmpzMat saturation_relation_kernel_gram_lambda{0, 0};
    bool saturation_relation_kernel_gram_valid = false;

    // Smith form and right transform of the HNF basis of the last
    // publish_presentation.  Relation appends that leave the HNF basis
    // unchanged republish the same basis, and then the Smith form is taken
    // from here.  Keyed only on the basis, of which S and V are a function,
    // so it stays valid across appends, factor-base rebuilds and clear().
    FiniteAbelianGroupSmithCache published_smith_form;
    // Cyclotomic automorphism data for Galois images of relations, built
    // on first use for the order the relation search runs on
    // (relation_candidate.cpp, the Galois-image step).
    CyclotomicGaloisImages cyclotomic_galois_images;
    // Galois-image relations stored in the relation matrix of this context
    // (relation_candidate.cpp, the Galois-image step), and those of them
    // stored without raising the modular rank (kernel images; those
    // stored under the capped rule below galois_image_min_degree are in
    // galois_image_capped_kernel_relations instead).  On the routes that
    // use capped_relation_count and the kernel-image count (the LLL route,
    // the goal checks of relation_search.cpp), images are left out of the
    // relation cap max_relations, and the kernel images counted here out
    // of the kernel-unit target; other readers of the plain counts include
    // them.  Images stored on the admission-cache path are counted too.  A
    // factor-base rebuild starts a fresh context, so the counts always
    // refer to the stored relations.
    slong galois_image_relations = 0;
    slong galois_image_kernel_relations = 0;
    // Kernel images stored under the capped rule of Phi_n fields below
    // galois_image_min_degree (relation_candidate.cpp, the Galois-image
    // step).  They count toward the kernel-unit target, so they are not in
    // galois_image_kernel_relations; this count bounds how many are stored.
    slong galois_image_capped_kernel_relations = 0;
    // galois_image_capped_kernel_relations at the start of the running
    // relation round plus the round's relation target (the LLL route sets
    // it per round); a capped kernel image is stored only while the count
    // is below this limit plus |G| - 1.  Negative when no round is running:
    // no cap.
    slong galois_image_round_kernel_limit = -1;
    // Smallest field degree at which the Galois-image step runs.  Only
    // internal tests change it (ClassGroupRelationSearchAccess::
    // set_galois_image_min_degree), to exercise the step on small fields;
    // clear() keeps it.
    slong galois_image_min_degree = kLargeFieldRelationSearchMinDegree;

    void invalidate_reduced_relation_kernel() noexcept {
        reduced_relation_kernel = flint::FmpzMat(0, 0);
        reduced_relation_kernel_valid = false;
        saturation_relation_kernel = flint::FmpzMat(0, 0);
        saturation_relation_kernel_valid = false;
        saturation_relation_image_transform = flint::FmpzMat(0, 0);
        saturation_relation_image_valid = false;
        saturation_relation_kernel_gram_d = flint::FmpzMat(0, 0);
        saturation_relation_kernel_gram_lambda = flint::FmpzMat(0, 0);
        saturation_relation_kernel_gram_valid = false;
    }

    void clear() noexcept {
        partial_relations.clear();
        direct_partial_residue_blocks.clear();
        partial_relations_max = 0;
        partial_relations_nonproductive_streak = 0;
        use_partial_relations = false;
        partial_relations_configured = false;

        random_ideal_search_state.base_indices.clear();
        random_ideal_search_state.exponents.clear();
        random_ideal_search_state.base_ideals.clear();
        random_ideal_search_state.inverse_base_numerators.clear();
        random_ideal_search_state.inverse_base_denominators.clear();
        random_ideal_search_state.random_ideal.clear();
        flint::fmpz_zero(
                flint::FmpzRef(random_ideal_search_state.lower_bound));
        flint::fmpz_zero(
                flint::FmpzRef(random_ideal_search_state.upper_bound));
        random_ideal_search_state.random_state = UWORD(0);
        random_ideal_search_state.defined = false;
        native_post_finite_reduced_lattices.clear();
        lll_rounded_embedding_cache =
                relation_search::IdealLatticeRoundedEmbeddingCache();
        generator_relation_row_capture = nullptr;
        generator_relation_row_captured = false;
        finish_state.reset();
        hnf_finish_workspace.clear();
        analytic_finish_precision = 0;
        flint::arb_zero(flint::ArbRef(analytic_finish_product));
        analytic_finish_product_precision = 0;
        analytic_finish_product_valid = false;

        factor_base_generation_records.clear();
        flint::fmpz_zero(flint::FmpzRef(grh_generation_bound));
        grh_generation_covered = false;
        grh_generation_basis = ClassGroupGenerationBasis::none;
        generation_different_route = GenerationDifferentRoute::none;
        class_unit_completeness = CertificationMode::unknown;
        class_unit_completeness_basis = ClassUnitCompletenessBasis::none;
        generation_structure_certified = false;
        relation_saturation_records.clear();
        relation_saturation_proof_records.clear();
        local_saturation_dlog_cache.clear();
        saturation_relation_scratch.clear();
        invalidate_reduced_relation_kernel();
        cyclotomic_galois_images.clear();
        galois_image_relations = 0;
        galois_image_kernel_relations = 0;
        galois_image_capped_kernel_relations = 0;
        galois_image_round_kernel_limit = -1;
    }
};

class ClassGroupRelationSearchAccess {
public:
    static bool saturate_relations_bounded(
            ClassGroupContext& context,
            bool& changed,
            bool& saturated,
            flint::FmpzConstRef aux_prime_bound,
            slong max_appends_per_ell,
            slong max_appends_total) noexcept {
        return context.saturate_relations_bounded_(
                changed, saturated, aux_prime_bound, max_appends_per_ell,
                max_appends_total);
    }

    static slong relation_kernel_units_target(
            const ClassGroupContext& context) noexcept {
        return context.relation_kernel_units_target_;
    }

    static void set_relation_kernel_units_target(
            ClassGroupContext& context,
            slong target) noexcept {
        context.relation_kernel_units_target_ = target;
    }

    static bool build_relation_factor_base(
            ClassGroupContext& context,
            flint::FmpzConstRef bound) noexcept {
        return context.build_relation_factor_base_(bound);
    }

    static const RelationMatrix& relations(
            const ClassGroupContext& context) noexcept {
        return context.relations_;
    }

    static void record_skipped_dependent(
            ClassGroupContext& context) noexcept {
        ++context.skipped_dependent_relations_;
    }

    // Small-LLL torsion-orbit skip hooks (lll_relation_search.cpp).
    static bool generator_relation_policy_keeps_nonduplicate(
            const ClassGroupContext& context) noexcept {
        return context.generator_relation_policy_ ==
               ClassGroupContext::DependentRelationPolicy::keep_nonduplicate;
    }

    // Sets the generator-relation policy the LLL relation route uses
    // (keep_nonduplicate).  White-box hook for the torsion-orbit skip test,
    // which runs the small-LLL candidate loop outside the route.
    static void use_route_generator_relation_policy(
            ClassGroupContext& context) noexcept {
        context.generator_relation_policy_ =
                ClassGroupContext::DependentRelationPolicy::keep_nonduplicate;
    }

    // The same stored-row test the keep_nonduplicate append makes before it
    // counts a skipped duplicate (append_relation_with_outcome_).
    static bool relations_contain_row(
            bool& contains,
            const ClassGroupContext& context,
            flint::FmpzMatConstRef row) noexcept {
        return context.relations_contain_row_(contains, row);
    }

    static void arm_generator_relation_row_capture(
            ClassGroupContext& context,
            flint::FmpzMat* row) noexcept {
        if (context.private_storage_ == nullptr) {
            return;
        }
        context.private_storage_->generator_relation_row_capture = row;
        context.private_storage_->generator_relation_row_captured = false;
    }

    // Disarms the capture; true when a row was captured while armed.
    static bool disarm_generator_relation_row_capture(
            ClassGroupContext& context) noexcept {
        if (context.private_storage_ == nullptr) {
            return false;
        }
        const bool captured =
                context.private_storage_->generator_relation_row_captured;
        context.private_storage_->generator_relation_row_capture = nullptr;
        context.private_storage_->generator_relation_row_captured = false;
        return captured;
    }

    static bool append_search_relation(
            bool& skipped_dependent,
            ClassGroupContext& context,
            const Relation& relation) noexcept {
        ClassGroupContext::RelationAppendOutcome outcome =
                ClassGroupContext::RelationAppendOutcome::none;
        if (!context.append_relation_with_outcome_(
                    outcome, relation, ClassGroupRelationSource::Search,
                    ClassGroupContext::DependentRelationPolicy::keep)) {
            return false;
        }
        skipped_dependent = outcome ==
                ClassGroupContext::RelationAppendOutcome::skipped_dependent;
        return true;
    }

    // Append `relation` as a supplied relation under the keep_nonduplicate
    // policy: an exact repeat of a stored row is skipped (counted as a
    // skipped dependent relation, `skipped` set), any other row is kept even
    // when it is dependent.  White-box hook for the duplicate-row tests.
    static bool append_keep_nonduplicate_relation(
            bool& skipped,
            ClassGroupContext& context,
            const Relation& relation) noexcept {
        ClassGroupContext::RelationAppendOutcome outcome =
                ClassGroupContext::RelationAppendOutcome::none;
        if (!context.append_relation_with_outcome_(
                    outcome, relation, ClassGroupRelationSource::Supplied,
                    ClassGroupContext::DependentRelationPolicy::
                            keep_nonduplicate)) {
            return false;
        }
        skipped = outcome ==
                ClassGroupContext::RelationAppendOutcome::skipped_dependent;
        return true;
    }

    static bool defer_native_goal_publication(
            const ClassGroupContext& context) noexcept;

    // The unconditional generation requirement recorded with the installed
    // factor base: proven_generation_requirement for the context's order,
    // evaluated once when the base was built
    // (record_factor_base_generation_).  False, leaving `out` unchanged,
    // when nothing is recorded or `order` is not the context's order; the
    // caller then computes the requirement itself.
    static bool recorded_generation_requirement(
            flint::FmpzRef out,
            const ClassGroupContext& context,
            const Order& order) noexcept {
        if (!context.has_factor_base() ||
            !same_order_parent(&context.parent_, &order) ||
            flint::fmpz_is_zero(flint::FmpzConstRef(
                    context.factor_base_generation_bound_))) {
            return false;
        }
        flint::fmpz_set(out, flint::FmpzConstRef(
                                     context.factor_base_generation_bound_));
        return true;
    }

    // Record, or clear, in the attached class/unit transaction report the
    // unwitnessed required prime that ended a candidate attempt; no-ops
    // without a transaction.  Clearing also clears the report's
    // factor-base size-limit flag, so each attempt starts clean.
    static void record_factor_base_honesty_failure(
            ClassGroupContext& context,
            FactorBaseHonestyFailure&& failure) noexcept;
    static void clear_factor_base_honesty_failure(
            ClassGroupContext& context) noexcept;

    // Whether the attached class/unit transaction keeps a passed honesty
    // scan (TransactionFactorBaseHonesty) of a factor base equal to `base`
    // at the same active and required bounds.  False without a transaction.
    static bool transaction_factor_base_honest(
            const ClassGroupContext& context,
            const FactorBase& base,
            flint::FmpzConstRef active_bound,
            flint::FmpzConstRef required_bound) noexcept;
    // Keep, in the attached class/unit transaction, a passed honesty scan of
    // `base` at these bounds, replacing any earlier one; no-op without a
    // transaction.  Nothing is kept if the base cannot be copied.  Call only
    // after a passed scan: the kept verdict is reused without a rescan.
    static void record_transaction_factor_base_honest(
            ClassGroupContext& context,
            const FactorBase& base,
            flint::FmpzConstRef active_bound,
            flint::FmpzConstRef required_bound) noexcept;

    // The class-group candidate computation and extension entry points the
    // class/unit transaction calls for each attempt.
    static bool compute_relation_candidate(
            ClassGroupContext& context,
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const ClassGroupRelationOptions& options) noexcept {
        return context.compute_relation_candidate_(
                order, factor_base_bound, options);
    }
    static bool extend_relation_kernel_units(
            ClassGroupContext& context,
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const ClassGroupRelationOptions& options) noexcept {
        return context.extend_relation_kernel_units_(
                order, factor_base_bound, options);
    }

    // The LLL relation slice the proven class/unit continuation runs on the
    // current candidate.
    static bool extend_lll_relation_slice(
            ClassGroupContext& context,
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const ClassGroupRelationOptions& options) noexcept {
        return context.extend_lll_relation_slice_(
                order, factor_base_bound, options);
    }

    static bool sync_row_module_checkpoint(
            ClassGroupContext& context) noexcept {
        return context.sync_row_module_checkpoint_();
    }

    static slong relation_kernel_row_count(
            const ClassGroupContext& context) noexcept {
        return relation_kernel_row_count_from_dimensions(
                context.generator_count(), context.relation_count(),
                context.relation_rank());
    }

    static slong relation_kernel_row_count_from_dimensions(
            slong generators,
            slong relations,
            slong rank) noexcept {
        return rank == generators && relations > generators
                ? relations - generators
                : 0;
    }

    static std::size_t direct_partial_residue_block_cache_size(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? 0
                : context.private_storage_
                          ->direct_partial_residue_blocks.size();
    }

    // The stored partial (one-large-prime) relations of `context`, or
    // nullptr before its private storage exists.  White-box hook for the
    // large-prime matching tests.
    static std::vector<PartialRelationEntry>* partial_relations(
            ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_->partial_relations;
    }

    // `norm` must be the exact field norm of `generator` and
    // `integral_coordinates` its coordinates.
    static bool try_append_integral_generator_relation(
            ClassGroupContext& context,
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            ClassGroupRelationSource source) noexcept {
        return context.try_append_integral_generator_relation_(
                partial_throttle_exit, generator, integral_coordinates, norm,
                integral_coordinate_polynomial, source);
    }

    // `norm` must be the exact field norm of `generator`.
    static bool try_append_generator_relation_with_norm(
            ClassGroupContext& context,
            bool& partial_throttle_exit,
            const Element& generator,
            flint::FmpqConstRef norm,
            ClassGroupRelationSource source =
                    ClassGroupRelationSource::Supplied) noexcept {
        return context.try_append_generator_relation_with_norm_(
                partial_throttle_exit, generator, norm, source);
    }

    // Rank of the modular relation module (reference ModuleCtx Mp).
    static slong relation_module_rank(
            const ClassGroupContext& context) noexcept {
        return context.relation_module_.rank();
    }

    using DependentRelationPolicy =
            ClassGroupContext::DependentRelationPolicy;

    // While alive, generator relations are appended under `policy`
    // (DependentRelationPolicy); the previous policy is restored on
    // destruction.
    class GeneratorRelationPolicyScope {
    public:
        GeneratorRelationPolicyScope(
                ClassGroupContext& context,
                DependentRelationPolicy policy) noexcept
            : context_(context),
              previous_(context.generator_relation_policy_) {
            context_.generator_relation_policy_ = policy;
        }
        ~GeneratorRelationPolicyScope() noexcept {
            context_.generator_relation_policy_ = previous_;
        }
        GeneratorRelationPolicyScope(
                const GeneratorRelationPolicyScope&) = delete;
        GeneratorRelationPolicyScope& operator=(
                const GeneratorRelationPolicyScope&) = delete;

    private:
        ClassGroupContext& context_;
        DependentRelationPolicy previous_;
    };

    // Galois-image relations stored in `context`, and those stored without
    // raising the modular rank; see ClassGroupContextStorage.
    static slong galois_image_relations(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? 0
                : context.private_storage_->galois_image_relations;
    }

    static slong galois_image_kernel_relations(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? 0
                : context.private_storage_->galois_image_kernel_relations;
    }

    static slong galois_image_capped_kernel_relations(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? 0
                : context.private_storage_
                          ->galois_image_capped_kernel_relations;
    }

    static slong galois_image_round_kernel_limit(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? -1
                : context.private_storage_->galois_image_round_kernel_limit;
    }

    // Starts a relation round with relation target `need` for the capped
    // image rule; a negative `need` ends it (no cap).
    static void start_galois_image_round(ClassGroupContext& context,
                                         slong need) noexcept {
        if (context.private_storage_ != nullptr) {
            detail::ClassGroupContextStorage& storage =
                    *context.private_storage_;
            storage.galois_image_round_kernel_limit =
                    need < 0 ? -1
                             : storage.galois_image_capped_kernel_relations +
                                       need;
        }
    }

    // `capped`: the image was stored under the capped rule of Phi_n fields
    // below galois_image_min_degree, where a kernel image counts toward the
    // kernel-unit target and is counted in
    // galois_image_capped_kernel_relations instead.
    static void note_galois_image_relation(ClassGroupContext& context,
                                           bool kernel,
                                           bool capped = false) noexcept {
        if (context.private_storage_ != nullptr) {
            ++context.private_storage_->galois_image_relations;
            if (kernel && capped) {
                ++context.private_storage_
                          ->galois_image_capped_kernel_relations;
            } else if (kernel) {
                ++context.private_storage_->galois_image_kernel_relations;
            }
        }
    }

    // Stored relations counted toward the relation cap max_relations: all
    // stored relations except Galois images.
    static slong capped_relation_count(
            const ClassGroupContext& context) noexcept {
        const slong count = context.relation_count();
        const slong images = galois_image_relations(context);
        return count > images ? count - images : 0;
    }

    // The bound of the current factor base (zero when the base was set
    // directly rather than built from a bound).
    static flint::FmpzConstRef factor_base_build_bound(
            const ClassGroupContext& context) noexcept {
        return flint::FmpzConstRef(context.factor_base_build_bound_);
    }

    // Smallest field degree at which the Galois-image step runs; see
    // ClassGroupContextStorage::galois_image_min_degree.
    static slong galois_image_min_degree(
            const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? kLargeFieldRelationSearchMinDegree
                : context.private_storage_->galois_image_min_degree;
    }

    // Internal tests only: runs the Galois-image step from `degree` on.
    static bool set_galois_image_min_degree(ClassGroupContext& context,
                                            slong degree) noexcept {
        if (!context.ensure_private_storage_()) {
            return false;
        }
        context.private_storage_->galois_image_min_degree = degree;
        return true;
    }

    static CyclotomicGaloisImages* cyclotomic_galois_images(
            ClassGroupContext& context) noexcept {
        return context.ensure_private_storage_()
                ? &context.private_storage_->cyclotomic_galois_images
                : nullptr;
    }

    static ReducedIdealLatticeCache* native_post_finite_reduced_lattices(
            ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_
                           ->native_post_finite_reduced_lattices;
    }

    static const relation_search::IdealLatticeRoundedEmbeddingCache*
    lll_rounded_embedding_cache(const ClassGroupContext& context) noexcept {
        return context.private_storage_ == nullptr
                ? nullptr
                : &context.private_storage_->lll_rounded_embedding_cache;
    }
};

}  // namespace silex::detail
