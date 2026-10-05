#pragma once

#include <flint/flint.h>

#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_vec.hpp>

#include <optional>

namespace silex {

namespace detail {
class FiniteAbelianGroupAccess;
}

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
    // The relation-kernel rows form a Z-basis of the integer left kernel of
    // the relation matrix, and the invariant-generator relation combinations
    // come from the same left transform.  These witness rows are valid but not
    // canonical: they depend on the transform FLINT's hnf_transform picks when
    // there are more relations than generators, so only properties that hold
    // for every basis (the product with the relations, saturation) may be
    // relied on.  They are stable within one object.
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
    // Names in silex::detail (including detail::FiniteAbelianGroupAccess) are
    // internal to Silex and its tests, not a supported entry point.
    friend class detail::FiniteAbelianGroupAccess;

    bool set_relation_matrix_with_hnf_basis(
            flint::FmpzMatConstRef relations,
            flint::FmpzMatConstRef hnf_basis) noexcept;
    bool ensure_left_transform() const noexcept;

    flint::FmpzMat relations_{0, 0};
    flint::FmpzMat hnf_basis_{0, 0};
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
