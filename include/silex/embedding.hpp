#pragma once

#include <flint/flint.h>

#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/flint/acb.hpp>
#include <silex/flint/acb_poly.hpp>
#include <silex/flint/acb_vec.hpp>
#include <silex/number_field.hpp>
#include <silex/signature.hpp>

namespace silex {

namespace flint {
class ArbVecConstRef;
class ArbVecRef;
}  // namespace flint

enum class LogEmbeddingMode {
    plain = 0,
    product = 1,
};

class EmbeddingLogCache;

class EmbeddingContext {
public:
    EmbeddingContext() noexcept = default;
    explicit EmbeddingContext(const NumberField& parent) noexcept;
    ~EmbeddingContext() noexcept;

    EmbeddingContext(const EmbeddingContext&) = delete;
    EmbeddingContext& operator=(const EmbeddingContext&) = delete;

    EmbeddingContext(EmbeddingContext&& other) noexcept;
    EmbeddingContext& operator=(EmbeddingContext&& other) noexcept;

    void swap(EmbeddingContext& other) noexcept;
    void clear() noexcept;
    bool define(const NumberField& parent) noexcept;

    bool is_defined() const noexcept;
    const NumberField* parent() const noexcept;
    slong degree() const noexcept;
    bool is_set() const noexcept;
    slong precision() const noexcept;
    Signature signature() const noexcept;

    bool refine(slong precision) noexcept;
    bool refine(slong precision,
                const DiagnosticsContext* diagnostics) noexcept;
    bool get_root(flint::AcbRef out, slong index) const noexcept;
    bool evaluate(flint::AcbRef out,
                  const Element& element,
                  slong index,
                  slong precision) noexcept;
    bool evaluate_all(flint::AcbVecRef out,
                      const Element& element,
                      slong precision,
                      const DiagnosticsContext* diagnostics = nullptr)
            noexcept;

private:
    NumberField parent_;
    slong degree_ = 0;
    Signature sig_;
    flint::AcbVec roots_{0};
    flint::AcbPolyEvaluationTree evaluation_tree_;
    slong prec_ = 0;
    bool roots_are_set_ = false;
    // Opaque implementation cache keeps STL storage out of this public header.
    EmbeddingLogCache* log_embedding_cache_ = nullptr;

    bool cached_log_embedding_(flint::ArbVecRef out,
                               const Element& element,
                               LogEmbeddingMode mode,
                               slong precision) const noexcept;
    bool store_log_embedding_(const Element& element,
                              LogEmbeddingMode mode,
                              slong precision,
                              flint::ArbVecConstRef values,
                              bool* updated_existing = nullptr) noexcept;

    friend bool logarithmic_embedding(flint::ArbVecRef out,
                                      EmbeddingContext& embeddings,
                                      const Element& element,
                                      LogEmbeddingMode mode,
                                      slong precision,
                                      const DiagnosticsContext* diagnostics)
            noexcept;
};

inline void swap(EmbeddingContext& left, EmbeddingContext& right) noexcept {
    left.swap(right);
}

}  // namespace silex
