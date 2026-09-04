#include <silex/embedding.hpp>

#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/fmpq_poly.h>

#include <silex/flint/acb_poly.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_poly.hpp>

#include <new>
#include <utility>
#include <vector>

namespace silex {

namespace {

enum class RootRefinementStatus {
    success,
    invalid_input,
    not_isolated,
    low_accuracy,
    real_root_validation,
    root_isolation,
    root_pairing,
    precision_overflow
};

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
    case RootRefinementStatus::precision_overflow:
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine_from_initial.failure.precision_overflow");
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

bool sort_and_pair_roots(acb_ptr roots, acb_ptr scratch, slong degree) noexcept {
    _acb_vec_sort_pretty(roots, degree);

    slong num_real = 0;
    for (slong i = 0; i < degree; ++i) {
        if (acb_is_real(roots + i) != 0) {
            ++num_real;
        }
    }

    if (degree == num_real) {
        return true;
    }

    slong positive_complex = 0;
    for (slong i = num_real; i < degree; ++i) {
        if (arb_is_positive(acb_imagref(roots + i)) != 0) {
            acb_swap(scratch + positive_complex, roots + i);
            ++positive_complex;
        }
    }

    const slong expected_complex_pairs = (degree - num_real) / 2;
    if (positive_complex != expected_complex_pairs) {
        return false;
    }

    for (slong i = 0; i < expected_complex_pairs; ++i) {
        acb_swap(roots + num_real + 2 * i, scratch + i);
        acb_conj(roots + num_real + 2 * i + 1,
                 roots + num_real + 2 * i);
    }
    return true;
}

RootRefinementStatus refine_roots_from_initial(
        acb_ptr roots,
        const fmpq_poly_t polynomial,
        slong degree,
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

    for (slong work_precision = precision; work_precision > 0;) {
        acb_poly_set_fmpq_poly(complex_polynomial.raw(), polynomial,
                               work_precision);

        // Source trace: reference Misc/acb_root_ctx.jl::_roots! refines roots by
        // calling FLINT acb_poly_find_roots with the previous root vector as
        // initial approximations, then validates and sorts the result.
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
            } else if (!sort_and_pair_roots(roots, scratch.data(), degree)) {
                last_failure = RootRefinementStatus::root_pairing;
                profile_root_refinement_failure(diagnostics, last_failure);
            } else {
                SILEX_PROFILE_EVENT(
                        diagnostics, DiagnosticsModule::element,
                        "element.embedding.refine_from_initial.success");
                return RootRefinementStatus::success;
            }
        }

        if (work_precision > (COEFF_MAX / 2)) {
            profile_root_refinement_failure(
                    diagnostics, RootRefinementStatus::precision_overflow);
            return RootRefinementStatus::precision_overflow;
        }
        work_precision *= 2;
    }
    return last_failure;
}

}  // namespace

class EmbeddingLogCache {
public:
    bool lookup(flint::ArbVecRef out,
                const Element& element,
                LogEmbeddingMode mode,
                slong precision) const noexcept {
        if (precision <= 0 || out.length() < 0) {
            return false;
        }
        for (const Entry& entry : entries_) {
            if (entry.mode == static_cast<int>(mode) &&
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
        for (Entry& entry : entries_) {
            if (entry.mode == static_cast<int>(mode) &&
                entry.element.equal(element)) {
                if (updated_existing != nullptr) {
                    *updated_existing = true;
                }
                if (entry.values.length() != values.length()) {
                    return false;
                }
                if (entry.precision <= precision) {
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
        return entry.element.set(element) && entry.values.set_from(values);
    }

private:
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

    std::vector<Entry> entries_;
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
    std::swap(prec_, other.prec_);
    std::swap(roots_are_set_, other.roots_are_set_);
    std::swap(log_embedding_cache_, other.log_embedding_cache_);
}

void EmbeddingContext::clear() noexcept {
    delete log_embedding_cache_;
    log_embedding_cache_ = nullptr;
    roots_ = flint::AcbVec(0);
    evaluation_tree_.clear();
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

    Signature next_sig;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "element.embedding.refine.signature");
        if (!silex::signature(next_sig, parent_)) {
            return false;
        }
    }

    flint::AcbVec next_roots(degree_);
    bool roots_refined = false;
    if (roots_are_set_) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.embedding.refine.previous_roots_available");
        next_roots.set_from(flint::AcbVecConstRef(roots_));
        const RootRefinementStatus refine_status = refine_roots_from_initial(
                next_roots.data(), parent_.raw_flint_field()->pol, degree_,
                precision, diagnostics);
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
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "element.embedding.refine.full_roots");
            arb_fmpz_poly_complex_roots(next_roots.data(), numerator.raw(), 0,
                                        precision);
        }
    }

    roots_.set_from(flint::AcbVecConstRef(next_roots));
    evaluation_tree_.clear();
    sig_ = next_sig;
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

    flint::AcbPoly complex_polynomial;
    flint::AcbVec values(degree_);
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

}  // namespace silex
