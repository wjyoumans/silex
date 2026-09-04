#include <silex/archimedean.hpp>

#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>

#include <silex/flint/acb_vec.hpp>

namespace silex {
namespace {

bool embedding_signature(Signature& out,
                         const EmbeddingContext& embeddings) noexcept {
    if (!embeddings.is_defined()) {
        return false;
    }
    if (embeddings.is_set()) {
        out = embeddings.signature();
        return true;
    }
    const NumberField* parent = embeddings.parent();
    return parent != nullptr && signature(out, *parent);
}

bool place_root_index(slong& root_index,
                      slong& r1,
                      const EmbeddingContext& embeddings,
                      slong place) noexcept {
    Signature sig;
    if (!embedding_signature(sig, embeddings)) {
        return false;
    }

    r1 = sig.r1();
    const slong places = sig.r1() + sig.r2();
    if (place < 0 || place >= places) {
        return false;
    }

    root_index = place < r1 ? place : r1 + 2 * (place - r1);
    return true;
}

bool num_places(slong& places, const EmbeddingContext& embeddings) noexcept {
    Signature sig;
    if (!embedding_signature(sig, embeddings)) {
        return false;
    }
    places = sig.r1() + sig.r2();
    return true;
}

enum class LogCacheProfileEventKind {
    hit,
    store,
};

void profile_log_cache_precision(
        const DiagnosticsContext* diagnostics,
        LogCacheProfileEventKind kind,
        slong precision) noexcept {
    if (kind == LogCacheProfileEventKind::hit) {
        if (precision <= 32) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.logarithmic_embedding.cache_hit.precision.le32");
        } else if (precision <= 64) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.logarithmic_embedding.cache_hit.precision.le64");
        } else if (precision <= 128) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.logarithmic_embedding.cache_hit.precision.le128");
        } else if (precision <= 256) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.logarithmic_embedding.cache_hit.precision.le256");
        } else if (precision <= 512) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.logarithmic_embedding.cache_hit.precision.le512");
        } else {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "element.logarithmic_embedding.cache_hit.precision.gt512");
        }
        return;
    }

    if (precision <= 32) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.logarithmic_embedding.cache_store.precision.le32");
    } else if (precision <= 64) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.logarithmic_embedding.cache_store.precision.le64");
    } else if (precision <= 128) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.logarithmic_embedding.cache_store.precision.le128");
    } else if (precision <= 256) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.logarithmic_embedding.cache_store.precision.le256");
    } else if (precision <= 512) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.logarithmic_embedding.cache_store.precision.le512");
    } else {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.logarithmic_embedding.cache_store.precision.gt512");
    }
}

}  // namespace

bool archimedean_absolute(flint::ArbRef out,
                          EmbeddingContext& embeddings,
                          const Element& element,
                          slong place,
                          ArchAbsMode mode,
                          slong precision) noexcept {
    if (precision <= 0 ||
        (mode != ArchAbsMode::plain && mode != ArchAbsMode::product)) {
        return false;
    }

    slong root_index = 0;
    slong r1 = 0;
    if (!place_root_index(root_index, r1, embeddings, place)) {
        return false;
    }

    flint::Acb value;
    if (!embeddings.evaluate(flint::AcbRef(value), element, root_index, precision)) {
        return false;
    }

    flint::Arb absolute;
    acb_abs(absolute.raw(), value.raw(), precision);
    if (mode == ArchAbsMode::product && place >= r1) {
        arb_mul(absolute.raw(), absolute.raw(), absolute.raw(), precision);
    }

    arb_set(out.raw(), absolute.raw());
    return true;
}

bool logarithmic_embedding(flint::ArbVecRef out,
                           EmbeddingContext& embeddings,
                           const Element& element,
                           LogEmbeddingMode mode,
                           slong precision,
                           const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.logarithmic_embedding");
    if (precision <= 0 ||
        (mode != LogEmbeddingMode::plain && mode != LogEmbeddingMode::product) ||
        element.equal_si(0)) {
        return false;
    }

    slong places = 0;
    if (!num_places(places, embeddings) || out.length() != places) {
        return false;
    }
    if (embeddings.cached_log_embedding_(out, element, mode, precision)) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.logarithmic_embedding.cache_hit");
        profile_log_cache_precision(diagnostics,
                                    LogCacheProfileEventKind::hit,
                                    precision);
        return true;
    }

    const slong degree = embeddings.degree();
    flint::AcbVec values(degree);
    flint::ArbVec tmp(places);
    if (!embeddings.evaluate_all(flint::AcbVecRef(values), element,
                                 precision, diagnostics)) {
        return false;
    }

    const Signature sig = embeddings.signature();
    const slong r1 = sig.r1();
    for (slong i = 0; i < places; ++i) {
        const slong root_index = i < r1 ? i : r1 + 2 * (i - r1);
        acb_abs(tmp.data() + i, values.data() + root_index, precision);

        if (mode == LogEmbeddingMode::product && i >= r1) {
            arb_mul(tmp.data() + i, tmp.data() + i, tmp.data() + i, precision);
        }

        if (arb_contains_zero(tmp.data() + i) != 0 ||
            arb_is_positive(tmp.data() + i) == 0) {
            return false;
        }

        arb_log(tmp.data() + i, tmp.data() + i, precision);
    }

    if (!out.set_from(flint::ArbVecConstRef(tmp))) {
        return false;
    }
    bool updated_existing_cache_entry = false;
    if (embeddings.store_log_embedding_(
                element, mode, precision, flint::ArbVecConstRef(tmp),
                &updated_existing_cache_entry)) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "element.logarithmic_embedding.cache_store");
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                updated_existing_cache_entry
                        ? "element.logarithmic_embedding.cache_store_upgrade"
                        : "element.logarithmic_embedding.cache_store_new");
        profile_log_cache_precision(diagnostics,
                                    LogCacheProfileEventKind::store,
                                    precision);
    }
    return true;
}

bool minkowski_embedding(flint::ArbMatRef out,
                         EmbeddingContext& embeddings,
                         const Element& element,
                         MinkowskiEmbeddingMode mode,
                         slong precision) noexcept {
    if (precision <= 0 ||
        (mode != MinkowskiEmbeddingMode::plain &&
         mode != MinkowskiEmbeddingMode::weighted)) {
        return false;
    }

    const slong degree = embeddings.degree();
    if (arb_mat_nrows(out.raw()) != 1 || arb_mat_ncols(out.raw()) != degree) {
        return false;
    }

    Signature sig;
    if (!embedding_signature(sig, embeddings)) {
        return false;
    }

    flint::AcbVec values(degree);
    flint::ArbMat tmp(1, degree);
    flint::Arb sqrt2;
    if (mode == MinkowskiEmbeddingMode::weighted) {
        arb_set_ui(sqrt2.raw(), 2);
        arb_sqrt(sqrt2.raw(), sqrt2.raw(), precision);
    }

    if (!embeddings.evaluate_all(flint::AcbVecRef(values), element, precision)) {
        return false;
    }

    const slong r1 = sig.r1();
    const slong r2 = sig.r2();
    for (slong i = 0; i < r1; ++i) {
        arb_set(arb_mat_entry(tmp.raw(), 0, i), acb_realref(values.data() + i));
    }

    for (slong i = 0; i < r2; ++i) {
        const slong root_index = r1 + 2 * i;
        const slong col_re = r1 + 2 * i;
        const slong col_im = col_re + 1;

        arb_set(arb_mat_entry(tmp.raw(), 0, col_re),
                acb_realref(values.data() + root_index));
        arb_set(arb_mat_entry(tmp.raw(), 0, col_im),
                acb_imagref(values.data() + root_index));

        if (mode == MinkowskiEmbeddingMode::weighted) {
            arb_mul(arb_mat_entry(tmp.raw(), 0, col_re),
                    arb_mat_entry(tmp.raw(), 0, col_re), sqrt2.raw(), precision);
            arb_mul(arb_mat_entry(tmp.raw(), 0, col_im),
                    arb_mat_entry(tmp.raw(), 0, col_im), sqrt2.raw(), precision);
        }
    }

    arb_mat_set(out.raw(), tmp.raw());
    return true;
}

}  // namespace silex
