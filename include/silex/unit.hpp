#pragma once

#include <flint/flint.h>

#include <cstddef>

#include <silex/archimedean.hpp>
#include <silex/element.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

namespace silex {

class ElementSpan {
public:
    ElementSpan() noexcept = default;
    ElementSpan(const Element* data, std::size_t size) noexcept
        : data_(data),
          size_(size) {
    }

    const Element* data() const noexcept { return data_; }
    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    const Element& operator[](std::size_t index) const noexcept {
        return data_[index];
    }
    const Element* begin() const noexcept { return data_; }
    const Element* end() const noexcept {
        return data_ == nullptr ? nullptr : data_ + size_;
    }

private:
    const Element* data_ = nullptr;
    std::size_t size_ = 0;
};

bool unit_rank(slong& rank, const NumberField& field) noexcept;

bool roots_of_unity(flint::FmpzRef order,
                    Element& generator,
                    const NumberField& field) noexcept;
bool root_of_unity_order(flint::FmpzRef order,
                         const NumberField& field) noexcept;
bool root_of_unity_generator(Element& generator,
                             const NumberField& field) noexcept;

bool unit_lower_regulator_bound(flint::ArbRef out,
                                const NumberField& field,
                                slong precision) noexcept;

// Accepts an explicit real quadratic field or a canonical polynomial-defined
// x^2-d field with positive squarefree nonsquare integer d.
bool quadratic_fundamental_unit(Element& out,
                                const NumberField& field) noexcept;

bool unit_log_matrix(flint::ArbMatRef out,
                     EmbeddingContext& embeddings,
                     ElementSpan units,
                     LogEmbeddingMode mode,
                     slong precision) noexcept;

bool unit_regulator(flint::ArbRef out,
                    EmbeddingContext& embeddings,
                    ElementSpan units,
                    slong precision) noexcept;

bool units_independent(bool& independent,
                       EmbeddingContext& embeddings,
                       ElementSpan units,
                       slong precision) noexcept;

}  // namespace silex
