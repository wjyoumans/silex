#pragma once

#include <flint/flint.h>

#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arb_vec.hpp>

namespace silex {

enum class ArchAbsMode {
    plain = 0,
    product = 1,
};

enum class MinkowskiEmbeddingMode {
    plain = 0,
    weighted = 1,
};

bool archimedean_absolute(flint::ArbRef out,
                          EmbeddingContext& embeddings,
                          const Element& element,
                          slong place,
                          ArchAbsMode mode,
                          slong precision) noexcept;

bool logarithmic_embedding(flint::ArbVecRef out,
                           EmbeddingContext& embeddings,
                           const Element& element,
                           LogEmbeddingMode mode,
                           slong precision,
                           const DiagnosticsContext* diagnostics = nullptr)
        noexcept;

bool minkowski_embedding(flint::ArbMatRef out,
                         EmbeddingContext& embeddings,
                         const Element& element,
                         MinkowskiEmbeddingMode mode,
                         slong precision) noexcept;

}  // namespace silex
