#include <silex/embedding.hpp>

#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/fmpq_poly.h>

#include <silex/flint/acb_poly.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>

#include "../factored_element/factor_screen_hash_internal.hpp"
#include "embedding_internal.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <new>
#include <utility>
#include <vector>

namespace silex {

using detail::RootRefinementStatus;

namespace {

void profile_root_refinement_failure(
        const DiagnosticsContext* diagnostics,
        RootRefinementStatus status) noexcept {
    switch (status) {
    case RootRefinementStatus::invalid_input:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.invalid_input");
        break;
    case RootRefinementStatus::not_isolated:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.not_isolated");
        break;
    case RootRefinementStatus::low_accuracy:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.low_accuracy");
        break;
    case RootRefinementStatus::real_root_validation:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.real_root_validation");
        break;
    case RootRefinementStatus::root_isolation:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.root_isolation");
        break;
    case RootRefinementStatus::root_pairing:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.root_pairing");
        break;
    case RootRefinementStatus::precision_cap:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.precision_cap");
        break;
    case RootRefinementStatus::success:
        break;
    }
}

bool roots_have_accuracy(acb_srcptr roots, slong degree, slong precision)
        noexcept {
    for (slong i = 0; i < degree; ++i) {
        if (acb_rel_accuracy_bits(roots + i) < precision) {
            return false;
        }
    }
    return true;
}

bool roots_are_isolated(acb_srcptr roots, slong degree) noexcept {
    for (slong i = 0; i < degree; ++i) {
        if (arf_sgn(arb_midref(acb_imagref(roots + i))) < 0) {
            continue;
        }
        for (slong j = i + 1; j < degree; ++j) {
            if (arf_sgn(arb_midref(acb_imagref(roots + j))) >= 0 &&
                acb_overlaps(roots + i, roots + j) != 0) {
                return false;
            }
        }
    }
    return true;
}

void zero_certified_real_imaginary_parts(acb_ptr roots, slong degree) noexcept {
    for (slong i = 0; i < degree; ++i) {
        if (arb_contains_zero(acb_imagref(roots + i)) != 0) {
            arb_zero(acb_imagref(roots + i));
        }
    }
}

// Writes the certified roots in `roots` back in the place order of
// `previous`, the roots stored before this refinement, and restores the
// conjugate-pair layout: real places first, then for each complex place the
// root with positive imaginary part followed by its conjugate.
//
// Both vectors hold pairwise disjoint balls, each containing exactly one root
// of the defining polynomial.  If roots[j] overlaps previous[i] and no other
// previous ball, the root inside roots[j] (which lies in some previous ball)
// must lie in previous[i], so both balls enclose the same root.  The matching
// is accepted only when it is a bijection; otherwise this returns false and
// the caller retries at higher precision.
//
// This replaces re-sorting with FLINT _acb_vec_sort_pretty after a refine.
// That comparator orders by |Im| and falls back to the real part only when
// the |Im| difference ball contains zero, so two conjugate pairs with nearly
// equal |Im| can sort differently at different precisions.  Place indices
// must stay fixed once published because cached logarithmic embeddings are
// served across refinements.
bool restore_previous_root_order(acb_ptr roots,
                                 acb_srcptr previous,
                                 acb_ptr scratch,
                                 slong degree,
                                 slong num_real) noexcept {
    if (num_real < 0 || num_real > degree || (degree - num_real) % 2 != 0) {
        return false;
    }
    for (slong i = 0; i < degree; ++i) {
        slong matches = 0;
        for (slong j = 0; j < degree; ++j) {
            if (acb_overlaps(roots + j, previous + i) != 0) {
                ++matches;
            }
        }
        if (matches != 1) {
            return false;
        }
    }
    for (slong j = 0; j < degree; ++j) {
        slong match = -1;
        for (slong i = 0; i < degree; ++i) {
            if (acb_overlaps(roots + j, previous + i) != 0) {
                if (match >= 0) {
                    return false;
                }
                match = i;
            }
        }
        if (match < 0) {
            return false;
        }
        acb_set(scratch + match, roots + j);
    }

    for (slong i = 0; i < num_real; ++i) {
        if (acb_is_real(scratch + i) == 0) {
            return false;
        }
    }
    for (slong i = num_real; i < degree; i += 2) {
        if (arb_is_positive(acb_imagref(scratch + i)) == 0) {
            return false;
        }
        acb_conj(scratch + i + 1, scratch + i);
    }
    for (slong i = 0; i < degree; ++i) {
        acb_swap(roots + i, scratch + i);
    }
    return true;
}

// Guard bits added to the requested precision for the first attempt from the
// previous roots.  Root finding at exactly the requested precision returns
// balls a few bits short of it, so that attempt always failed the accuracy
// check and was repeated at double the precision.  This is a tuning constant
// of the precision schedule, not a correctness condition: every attempt is
// checked the same way.
constexpr slong refine_guard_bits = 32;

// Largest working precision tried from the previous roots, as a multiple of
// the starting working precision, before falling back to full isolation.
constexpr slong refine_work_precision_factor = 8;

// Largest degree evaluated through a table of root powers instead of FLINT's
// product-tree multipoint evaluation.  A tuning constant, not a correctness
// condition: both give rigorous enclosures.
constexpr slong root_power_table_max_degree = 64;

}  // namespace

namespace detail {

RootRefinementStatus refine_roots_from_initial(
        acb_ptr roots,
        acb_srcptr previous,
        const fmpq_poly_t polynomial,
        slong degree,
        slong num_real,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.embedding.refine_from_initial");
    if (degree <= 0 || precision <= 0) {
        profile_root_refinement_failure(diagnostics,
                                        RootRefinementStatus::invalid_input);
        return RootRefinementStatus::invalid_input;
    }

    flint::AcbPoly complex_polynomial;
    flint::AcbVec scratch(degree);
    RootRefinementStatus last_failure = RootRefinementStatus::invalid_input;
    const slong start_precision = precision <= COEFF_MAX - refine_guard_bits
                                          ? precision + refine_guard_bits
                                          : COEFF_MAX;
    const slong max_work_precision =
            start_precision <= COEFF_MAX / refine_work_precision_factor
                    ? start_precision * refine_work_precision_factor
                    : COEFF_MAX;

    // Deviation from the reference, whose first attempt is at the requested
    // precision: start a few guard bits above it (see refine_guard_bits).
    for (slong work_precision = start_precision; work_precision > 0;) {
        acb_poly_set_fmpq_poly(complex_polynomial.raw(), polynomial,
                               work_precision);

        // Source trace: reference Misc/acb_root_ctx.jl::_roots! refines roots by
        // calling FLINT acb_poly_find_roots with the previous root vector as
        // initial approximations, then validates the result.  Unlike the
        // reference, which re-sorts with _acb_vec_sort_pretty, the refined
        // roots are matched back to the previous place order by overlap
        // instead of resorted, so a place's index stays tied to the same
        // root across refinements and does not depend on sort order.
        const slong max_iterations =
                FLINT_MIN(FLINT_MAX(degree, work_precision / 4),
                          work_precision);
        slong isolated = 0;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::element,
                    "element.embedding.refine_from_initial.find_roots");
            isolated = acb_poly_find_roots(
                    roots, complex_polynomial.raw(), roots, max_iterations,
                    work_precision);
        }
        if (isolated != degree) {
            last_failure = RootRefinementStatus::not_isolated;
            profile_root_refinement_failure(diagnostics, last_failure);
        } else if (!roots_have_accuracy(roots, degree, precision)) {
            last_failure = RootRefinementStatus::low_accuracy;
            profile_root_refinement_failure(diagnostics, last_failure);
        } else if (acb_poly_validate_real_roots(
                           roots, complex_polynomial.raw(),
                           work_precision) == 0) {
            last_failure = RootRefinementStatus::real_root_validation;
            profile_root_refinement_failure(diagnostics, last_failure);
        } else {
            zero_certified_real_imaginary_parts(roots, degree);
            if (!roots_are_isolated(roots, degree)) {
                last_failure = RootRefinementStatus::root_isolation;
                profile_root_refinement_failure(diagnostics, last_failure);
            } else if (!restore_previous_root_order(roots, previous,
                                                    scratch.data(), degree,
                                                    num_real)) {
                last_failure = RootRefinementStatus::root_pairing;
                profile_root_refinement_failure(diagnostics, last_failure);
            } else {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::element,
                        "element.embedding.refine_from_initial.success");
                return RootRefinementStatus::success;
            }
        }

        // The reference doubles up to 2^22 bits.  Cap the attempts from the
        // previous roots at refine_work_precision_factor times the starting
        // precision (the requested precision plus refine_guard_bits); the
        // caller then falls back to full isolation.
        if (work_precision > max_work_precision / 2) {
            profile_root_refinement_failure(
                    diagnostics, RootRefinementStatus::precision_cap);
            return RootRefinementStatus::precision_cap;
        }
        work_precision *= 2;
    }
    return last_failure;
}

bool isolate_roots_in_previous_order(acb_ptr roots,
                                     acb_srcptr previous,
                                     const fmpz_poly_t numerator,
                                     slong degree,
                                     slong num_real,
                                     slong precision,
                                     const DiagnosticsContext* diagnostics,
                                     slong* attempts) noexcept {
    if (attempts != nullptr) {
        *attempts = 0;
    }
    if (degree <= 0 || precision <= 0) {
        return false;
    }

    flint::AcbVec scratch(degree);
    for (slong work_precision = precision;;) {
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.embedding.refine.full_roots");
            arb_fmpz_poly_complex_roots(roots, numerator, 0, work_precision);
        }
        if (attempts != nullptr) {
            ++*attempts;
        }
        // FLINT returns its own sorted order; map it back to the published
        // place order.  An ambiguous match means some new ball is wider than
        // the gap between previous balls, so isolate again more precisely.
        if (restore_previous_root_order(roots, previous, scratch.data(),
                                        degree, num_real)) {
            return true;
        }
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.embedding.refine.full_roots_order_retry");
        // Unreachable given the precondition: for a squarefree polynomial the
        // previous balls are at positive distance from every other root, so
        // the retry succeeds at a finite precision.
        // The guard only keeps the doubling from overflowing slong.
        if (work_precision > (COEFF_MAX / 2)) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.embedding.refine.full_roots_precision_overflow");
            return false;
        }
        work_precision *= 2;
    }
}

}  // namespace detail

namespace {

// Source of EmbeddingLogCache versions.  One counter for every cache keeps a
// version from being reused by a later cache.
std::atomic<std::uint64_t> next_log_cache_version{1};

std::uint64_t fresh_log_cache_version() noexcept {
    return next_log_cache_version.fetch_add(1, std::memory_order_relaxed);
}

}  // namespace

class EmbeddingLogCache {
public:
    // The version changes whenever a lookup could start returning different
    // values: when an existing entry is replaced (a precision upgrade), when
    // a store fails part way, and when a lookup missed and the computation
    // then failed, since a later store of that element at another precision
    // would turn the failed miss into a hit.  Adding a new entry otherwise
    // keeps the version: an earlier lookup that it would now match missed
    // and then either failed, which changed the version, or stored this
    // element already, so the store is not new.
    std::uint64_t version() const noexcept {
        return version_;
    }

    void mark_changed() noexcept {
        version_ = fresh_log_cache_version();
    }

    bool lookup(flint::ArbVecRef out,
                const Element& element,
                LogEmbeddingMode mode,
                slong precision) const noexcept {
        if (precision <= 0 || out.length() < 0) {
            return false;
        }
        const std::uint64_t hash = element_hash(element);
        const std::size_t mask = table_.size() - 1;
        for (std::size_t slot = first_slot(hash); slot != no_slot &&
                                                  table_[slot] != no_slot;
             slot = (slot + 1) & mask) {
            const std::size_t i = table_[slot];
            const Entry& entry = entries_[i];
            if (hashes_[i] == hash &&
                entry.mode == static_cast<int>(mode) &&
                entry.precision >= precision &&
                entry.values.length() == out.length() &&
                entry.element.equal(element)) {
                return out.set_from(flint::ArbVecConstRef(entry.values));
            }
        }
        return false;
    }

    bool store(const NumberField& parent,
               const Element& element,
               LogEmbeddingMode mode,
               slong precision,
               flint::ArbVecConstRef values,
               bool* updated_existing = nullptr) noexcept {
        if (updated_existing != nullptr) {
            *updated_existing = false;
        }
        if (precision <= 0 || values.length() < 0) {
            return false;
        }
        const std::uint64_t hash = element_hash(element);
        const std::size_t mask = table_.size() - 1;
        for (std::size_t slot = first_slot(hash); slot != no_slot &&
                                                  table_[slot] != no_slot;
             slot = (slot + 1) & mask) {
            const std::size_t i = table_[slot];
            Entry& entry = entries_[i];
            if (hashes_[i] == hash &&
                entry.mode == static_cast<int>(mode) &&
                entry.element.equal(element)) {
                if (updated_existing != nullptr) {
                    *updated_existing = true;
                }
                if (entry.values.length() != values.length()) {
                    return false;
                }
                if (entry.precision <= precision) {
                    version_ = fresh_log_cache_version();
                    if (!entry.values.set_from(values)) {
                        return false;
                    }
                    entry.precision = precision;
                }
                return true;
            }
        }

        entries_.emplace_back(parent, values.length());
        Entry& entry = entries_.back();
        entry.mode = static_cast<int>(mode);
        entry.precision = precision;
        const bool element_set = entry.element.set(element);
        // The hash of the element the entry actually holds, so the index
        // stays exact even when the copy failed.
        hashes_.push_back(element_set ? hash : element_hash(entry.element));
        index_last_entry();
        if (!element_set || !entry.values.set_from(values)) {
            version_ = fresh_log_cache_version();
            return false;
        }
        return true;
    }

private:
    static constexpr std::size_t no_slot = SIZE_MAX;

    // The slot where the probe for `hash` starts, or no_slot while the
    // table is empty.
    std::size_t first_slot(std::uint64_t hash) const noexcept {
        return table_.empty() ? no_slot
                              : static_cast<std::size_t>(hash) &
                                        (table_.size() - 1);
    }

    // Puts entry `index` in the first free slot of its probe sequence.
    void place(std::size_t index) noexcept {
        const std::size_t mask = table_.size() - 1;
        std::size_t slot = first_slot(hashes_[index]);
        while (table_[slot] != no_slot) {
            slot = (slot + 1) & mask;
        }
        table_[slot] = index;
    }

    // Adds the last entry to the table, first doubling the table and
    // placing every entry again in index order when it would be more than
    // half full.
    void index_last_entry() noexcept {
        if (2 * entries_.size() > table_.size()) {
            table_.assign(table_.empty() ? 16 : 2 * table_.size(), no_slot);
            for (std::size_t i = 0; i + 1 < entries_.size(); ++i) {
                place(i);
            }
        }
        place(entries_.size() - 1);
    }

    // Screen hash of an element; 0 for an undefined one, which
    // Element::equal never matches.
    static std::uint64_t element_hash(const Element& element) noexcept {
        const NumberField* parent = element.parent();
        return parent == nullptr
                       ? 0
                       : detail::factor_screen_hash(
                                 element, parent->raw_flint_field());
    }

    struct Entry {
        Entry(const NumberField& parent, slong places) noexcept
            : element(parent),
              values(places) {
        }

        Element element;
        flint::ArbVec values;
        int mode = 0;
        slong precision = 0;
    };

    // Entries in insertion order.  hashes_[i] is detail::factor_screen_hash
    // of entries_[i].element.  Element::equal compares stored forms (FLINT
    // nf_elem_equal), and equal stored forms hash equally, so an entry
    // with another hash cannot match; a hash match is still confirmed with
    // Element::equal.
    //
    // table_ is an open-addressing index of entries_ by hash, with linear
    // probing over a power-of-two number of slots kept at most half full;
    // no_slot marks a free slot.  Entries are placed in index order and
    // never removed, so the probe for a hash meets the entries with that
    // hash in insertion order, and a lookup returns the same entry as a
    // scan of entries_ in order.
    std::vector<Entry> entries_;
    std::vector<std::uint64_t> hashes_;
    std::vector<std::size_t> table_;
    std::uint64_t version_ = fresh_log_cache_version();
};

EmbeddingContext::EmbeddingContext(const NumberField& parent) noexcept {
    define(parent);
}

EmbeddingContext::~EmbeddingContext() noexcept {
    clear();
}

EmbeddingContext::EmbeddingContext(EmbeddingContext&& other) noexcept {
    swap(other);
}

EmbeddingContext& EmbeddingContext::operator=(EmbeddingContext&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void EmbeddingContext::swap(EmbeddingContext& other) noexcept {
    parent_.swap(other.parent_);
    std::swap(degree_, other.degree_);
    std::swap(sig_, other.sig_);
    roots_.swap(other.roots_);
    evaluation_tree_.swap(other.evaluation_tree_);
    root_power_table_.swap(other.root_power_table_);
    std::swap(prec_, other.prec_);
    std::swap(roots_are_set_, other.roots_are_set_);
    std::swap(log_embedding_cache_, other.log_embedding_cache_);
}

void EmbeddingContext::clear() noexcept {
    delete log_embedding_cache_;
    log_embedding_cache_ = nullptr;
    roots_ = flint::AcbVec(0);
    evaluation_tree_.clear();
    root_power_table_.clear();
    parent_.clear();
    degree_ = 0;
    sig_ = Signature();
    prec_ = 0;
    roots_are_set_ = false;
}

bool EmbeddingContext::define(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    EmbeddingContext next;
    next.parent_ = parent;
    next.degree_ = parent.degree();
    next.roots_ = flint::AcbVec(next.degree_);
    // The signature depends only on the defining polynomial, so it is known
    // before any root is computed.  It stays (0, 0) when it cannot be
    // computed, which only happens for a non-squarefree polynomial installed
    // behind the NumberField's back; refine() then fails.
    Signature sig;
    if (silex::signature(sig, parent) && sig.degree() == next.degree_) {
        next.sig_ = sig;
    }

    swap(next);
    return true;
}

bool EmbeddingContext::is_defined() const noexcept {
    return parent_.is_defined();
}

const NumberField* EmbeddingContext::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong EmbeddingContext::degree() const noexcept {
    return degree_;
}

bool EmbeddingContext::is_set() const noexcept {
    return parent_.is_defined() && roots_are_set_;
}

slong EmbeddingContext::precision() const noexcept {
    return prec_;
}

Signature EmbeddingContext::signature() const noexcept {
    return parent_.is_defined() ? sig_ : Signature();
}

bool EmbeddingContext::refine(slong precision) noexcept {
    return refine(precision, nullptr);
}

bool EmbeddingContext::refine(
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.embedding.refine");
    if (!parent_.is_defined() || precision <= 0) {
        return false;
    }
    if (roots_are_set_ && prec_ >= precision) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.embedding.refine.cache_hit");
        return true;
    }

    // define() stores the signature; a context without one has a defining
    // polynomial whose roots cannot be isolated.
    if (sig_.degree() != degree_) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.embedding.refine.no_signature");
        return false;
    }

    // Stored roots already accurate to the request pass every check a
    // refinement makes: they are isolated, validated and in place order, and
    // the accuracy check depends only on the balls.  Keep them.  Refining
    // them anyway costs FLINT 3.6.0 acb_poly_find_roots its full iteration
    // budget: when every correction is zero, its early-stop test overflows
    // and never fires.
    if (roots_are_set_ && roots_have_accuracy(roots_.data(), degree_,
                                              precision)) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.embedding.refine.accurate_roots");
        prec_ = precision;
        return true;
    }

    flint::AcbVec next_roots(degree_);
    bool roots_refined = false;
    if (roots_are_set_) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine.previous_roots_available");
        next_roots.set_from(flint::AcbVecConstRef(roots_));
        const RootRefinementStatus refine_status =
                detail::refine_roots_from_initial(
                next_roots.data(), roots_.data(),
                parent_.raw_flint_field()->pol, degree_, sig_.r1(), precision,
                diagnostics);
        roots_refined = refine_status == RootRefinementStatus::success;
    }

    if (!roots_refined) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                roots_are_set_
                        ? "element.embedding.refine.full_roots_fallback"
                        : "element.embedding.refine.full_roots_initial");
        flint::FmpzPoly numerator;
        fmpq_poly_get_numerator(numerator.raw(), parent_.raw_flint_field()->pol);
        if (roots_are_set_) {
            // Keep the published place order.
            if (!detail::isolate_roots_in_previous_order(
                        next_roots.data(), roots_.data(), numerator.raw(),
                        degree_, sig_.r1(), precision, diagnostics)) {
                return false;
            }
        } else {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.embedding.refine.full_roots");
            arb_fmpz_poly_complex_roots(next_roots.data(), numerator.raw(), 0,
                                        precision);
        }
    }

    roots_.set_from(flint::AcbVecConstRef(next_roots));
    evaluation_tree_.clear();
    root_power_table_.clear();
    prec_ = precision;
    roots_are_set_ = true;
    return true;
}

bool EmbeddingContext::get_root(flint::AcbRef out, slong index) const noexcept {
    if (!parent_.is_defined() || !roots_are_set_ ||
        index < 0 || index >= degree_) {
        return false;
    }
    acb_set(out.raw(), roots_.data() + index);
    return true;
}

bool EmbeddingContext::evaluate(flint::AcbRef out,
                                const Element& element,
                                slong index,
                                slong precision) noexcept {
    if (!parent_.is_defined() || precision <= 0 ||
        !element.has_parent(parent_) ||
        index < 0 || index >= degree_ || !refine(precision)) {
        return false;
    }

    flint::FmpqPoly polynomial;
    if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
        return false;
    }

    flint::AcbPoly complex_polynomial;
    flint::Acb value;
    acb_poly_set_fmpq_poly(complex_polynomial.raw(), polynomial.raw(), precision);
    acb_poly_evaluate(value.raw(), complex_polynomial.raw(),
                      roots_.data() + index, precision);
    acb_set(out.raw(), value.raw());
    return true;
}

bool EmbeddingContext::evaluate_all(flint::AcbVecRef out,
                                    const Element& element,
                                    slong precision,
                                    const DiagnosticsContext* diagnostics)
        noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.embedding.evaluate_all");
    if (!parent_.is_defined() || precision <= 0 ||
        !element.has_parent(parent_) ||
        out.length() != degree_ || !refine(precision, diagnostics)) {
        return false;
    }

    flint::FmpqPoly polynomial;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.embedding.evaluate_all_get_poly");
        if (!element.get_fmpq_poly(flint::FmpqPolyRef(polynomial))) {
            return false;
        }
    }

    flint::AcbVec values(degree_);
    if (degree_ <= root_power_table_max_degree &&
        polynomial.raw()->length <= degree_) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.embedding.evaluate_all_root_powers");
        // One dot product per root against its stored powers.  A table built
        // at a higher precision also serves this request.
        if (!root_power_table_.covers(degree_, precision)) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::element,
                    "element.embedding.evaluate_all_root_powers_build");
            if (!root_power_table_.build(roots_.data(), degree_,
                                         precision)) {
                return false;
            }
        } else {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.embedding.evaluate_all_root_powers_hit");
        }
        if (!root_power_table_.evaluate(values.data(), polynomial.raw(),
                                        precision)) {
            return false;
        }
        out.set_from(flint::AcbVecConstRef(values));
        return true;
    }

    flint::AcbPoly complex_polynomial;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.embedding.evaluate_all_evaluate_vec");
        acb_poly_set_fmpq_poly(complex_polynomial.raw(), polynomial.raw(),
                               precision);
        // FLINT's fast evaluator builds this documented product tree on every
        // call; reuse it while the embedding roots and precision are fixed.
        if (!evaluation_tree_.matches(degree_, precision)) {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::element,
                    "element.embedding.evaluate_all_evaluation_tree_build");
            if (!evaluation_tree_.build(roots_.data(), degree_, precision)) {
                return false;
            }
        } else {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.embedding.evaluate_all_evaluation_tree_hit");
        }
        if (!evaluation_tree_.evaluate(
                    values.data(),
                    flint::AcbPolyConstRef(complex_polynomial), precision)) {
            return false;
        }
    }
    out.set_from(flint::AcbVecConstRef(values));
    return true;
}

bool EmbeddingContext::cached_log_embedding_(
        flint::ArbVecRef out,
        const Element& element,
        LogEmbeddingMode mode,
        slong precision) const noexcept {
    if (!parent_.is_defined() || precision <= 0 ||
        !element.has_parent(parent_) || out.length() < 0) {
        return false;
    }
    return log_embedding_cache_ != nullptr &&
           log_embedding_cache_->lookup(out, element, mode, precision);
}

bool EmbeddingContext::store_log_embedding_(
        const Element& element,
        LogEmbeddingMode mode,
        slong precision,
        flint::ArbVecConstRef values,
        bool* updated_existing) noexcept {
    if (updated_existing != nullptr) {
        *updated_existing = false;
    }
    if (!parent_.is_defined() || precision <= 0 ||
        !element.has_parent(parent_) || values.length() < 0) {
        return false;
    }

    if (log_embedding_cache_ == nullptr) {
        log_embedding_cache_ = new (std::nothrow) EmbeddingLogCache();
    }
    if (log_embedding_cache_ == nullptr) {
        return false;
    }
    return log_embedding_cache_->store(parent_, element, mode, precision,
                                       values, updated_existing);
}

void EmbeddingContext::failed_log_embedding_miss_() noexcept {
    if (log_embedding_cache_ != nullptr) {
        log_embedding_cache_->mark_changed();
    }
}

namespace detail {

bool embedding_log_cache_version(std::uint64_t& out,
                                 const EmbeddingContext& embeddings) noexcept {
    out = 0;
    if (embeddings.log_embedding_cache_ == nullptr) {
        return false;
    }
    out = embeddings.log_embedding_cache_->version();
    return true;
}

}  // namespace detail

}  // namespace silex
