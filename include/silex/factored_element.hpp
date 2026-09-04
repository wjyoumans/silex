#pragma once

#include <flint/flint.h>

#include <silex/archimedean.hpp>
#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/embedding.hpp>
#include <silex/factored.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/number_field.hpp>

#include <vector>

namespace silex {

class ResidueFieldElement;

enum class FactoredRootStrategy {
    automatic,
    reduced,
    compact,
};

class FactoredElement {
public:
    FactoredElement() noexcept = default;
    explicit FactoredElement(const NumberField& parent) noexcept;
    ~FactoredElement() noexcept;

    FactoredElement(const FactoredElement&) = delete;
    FactoredElement& operator=(const FactoredElement&) = delete;

    FactoredElement(FactoredElement&& other) noexcept;
    FactoredElement& operator=(FactoredElement&& other) noexcept;

    void swap(FactoredElement& other) noexcept;
    void clear() noexcept;
    bool define(const NumberField& parent) noexcept;
    bool set(const FactoredElement& other) noexcept;

    bool is_defined() const noexcept;
    const NumberField* parent() const noexcept;
    slong length() const noexcept;
    FactorSpan<Element> factors() const noexcept;

    bool one() noexcept;
    bool set_element(const Element& element) noexcept;
    bool push(const Element& factor, slong exponent) noexcept;
    void normalize() noexcept;
    bool multiply(const FactoredElement& left,
                  const FactoredElement& right) noexcept;
    bool divide(const FactoredElement& left,
                const FactoredElement& right) noexcept;
    bool invert(const FactoredElement& input) noexcept;
    bool pow_si(const FactoredElement& input, slong exponent) noexcept;
    bool root_si(const FactoredElement& input,
                 slong exponent,
                 const DiagnosticsContext* diagnostics = nullptr) noexcept;
    bool is_square(bool& is_square,
                   FactoredElement& root,
                   const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;
    bool is_power_si(bool& is_power,
                     FactoredElement& root,
                     slong exponent,
                     const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;
    bool is_power_si(bool& is_power,
                     FactoredElement& root,
                     slong exponent,
                     FactoredRootStrategy strategy,
                     const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;
    const Element* factor(slong index) const noexcept;
    bool exponent(slong& out, slong index) const noexcept;
    bool evaluate(Element& out) const noexcept;
    bool logarithmic_embedding(flint::ArbVecRef out,
                               EmbeddingContext& embeddings,
                               LogEmbeddingMode mode,
                               slong precision,
                               const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;

private:
    bool push_copy(const Element& factor, slong exponent) noexcept;

    NumberField parent_;
    Factored<Element> value_;

    friend class ResidueFieldElement;
};

inline void swap(FactoredElement& left, FactoredElement& right) noexcept {
    left.swap(right);
}

class CompactElement {
public:
    CompactElement() noexcept = default;
    explicit CompactElement(const NumberField& parent) noexcept;
    ~CompactElement() noexcept;

    CompactElement(const CompactElement&) = delete;
    CompactElement& operator=(const CompactElement&) = delete;

    CompactElement(CompactElement&& other) noexcept;
    CompactElement& operator=(CompactElement&& other) noexcept;

    void swap(CompactElement& other) noexcept;
    void clear() noexcept;
    bool define(const NumberField& parent) noexcept;

    bool is_defined() const noexcept;
    const NumberField* parent() const noexcept;
    slong base() const noexcept;
    slong length() const noexcept;

    bool one(slong base) noexcept;
    bool set_factored_element(const FactoredElement& input,
                              slong base,
                              const DiagnosticsContext* diagnostics = nullptr)
            noexcept;
    bool evaluate(Element& out) const noexcept;
    bool logarithmic_embedding(flint::ArbVecRef out,
                               EmbeddingContext& embeddings,
                               LogEmbeddingMode mode,
                               slong precision) const noexcept;
    bool root_base(FactoredElement& root,
                   const DiagnosticsContext* diagnostics = nullptr)
            const noexcept;

private:
    enum class RootStatus {
        unsupported,
        not_power,
        power,
    };

    bool ensure_length(slong length) noexcept;
    void normalize() noexcept;
    RootStatus root_base_status(FactoredElement& root,
                                const DiagnosticsContext* diagnostics)
            const noexcept;

    NumberField parent_;
    slong base_ = 0;
    std::vector<FactoredElement> coefficients_;

    friend class FactoredElement;
};

inline void swap(CompactElement& left, CompactElement& right) noexcept {
    left.swap(right);
}

}  // namespace silex
