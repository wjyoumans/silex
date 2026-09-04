#pragma once

#include <flint/flint.h>

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_vec.hpp>

#include <optional>

namespace silex {

class FiniteAbelianGroup {
public:
    FiniteAbelianGroup() noexcept = default;
    ~FiniteAbelianGroup() noexcept;

    FiniteAbelianGroup(const FiniteAbelianGroup&) = delete;
    FiniteAbelianGroup& operator=(const FiniteAbelianGroup&) = delete;

    FiniteAbelianGroup(FiniteAbelianGroup&& other) noexcept;
    FiniteAbelianGroup& operator=(FiniteAbelianGroup&& other) noexcept;

    void swap(FiniteAbelianGroup& other) noexcept;
    void clear() noexcept;
    bool set(const FiniteAbelianGroup& other) noexcept;

    bool is_defined() const noexcept;
    bool set_relation_matrix(flint::FmpzMatConstRef relations) noexcept;
    bool set_relation_matrix_with_hnf_basis(
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis) noexcept;

    slong relation_count() const noexcept;
    slong generator_count() const noexcept;
    slong invariant_count() const noexcept;
    slong relation_kernel_count() const noexcept;

    bool relations(flint::FmpzMatRef out) const noexcept;
    bool invariant(flint::FmpzRef out, slong index) const noexcept;
    bool invariants(flint::FmpzVecRef out) const noexcept;
    bool order(flint::FmpzRef out) const noexcept;
    std::optional<flint::FmpzMat> relations() const noexcept;
    std::optional<flint::Fmpz> invariant(slong index) const noexcept;
    std::optional<flint::Fmpz> order() const noexcept;

    bool invariant_generator_matrix(flint::FmpzMatRef out) const noexcept;
    bool invariant_generator_relation_matrix(flint::FmpzMatRef out) const noexcept;
    bool relation_kernel_row(flint::FmpzMatRef out, slong index) const noexcept;
    bool relation_kernel_matrix(flint::FmpzMatRef out) const noexcept;
    std::optional<flint::FmpzMat> invariant_generator_matrix() const noexcept;
    std::optional<flint::FmpzMat> relation_kernel_matrix() const noexcept;
    bool reduce(flint::FmpzMatRef row) const noexcept;
    bool invariant_coordinates(flint::FmpzMatRef out,
                               flint::FmpzMatConstRef row) const noexcept;
    std::optional<flint::FmpzMat> invariant_coordinates(
            flint::FmpzMatConstRef row) const noexcept;

private:
    bool ensure_left_transform() const noexcept;

    flint::FmpzMat relations_{0, 0};
    mutable flint::FmpzMat snf_{0, 0};
    mutable flint::FmpzMat left_transform_{0, 0};
    mutable flint::FmpzMat right_transform_{0, 0};
    mutable flint::FmpzMat right_transform_inv_{0, 0};
    mutable flint::FmpzVec invariants_{0};
    slong generator_count_ = 0;
    bool defined_ = false;
    mutable bool has_left_transform_ = false;
};

inline void swap(FiniteAbelianGroup& left,
                 FiniteAbelianGroup& right) noexcept {
    left.swap(right);
}

}  // namespace silex
