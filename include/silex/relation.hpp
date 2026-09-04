#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/diagnostics.hpp>
#include <silex/factor_base.hpp>
#include <silex/fmpz_smat.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_poly.hpp>

#include <vector>

namespace silex {

class ClassGroupContext;
class FiniteAbelianGroup;
class Relation;

namespace detail {

bool set_relation_from_known_row(Relation& out,
                                 const FactorBase& base,
                                 const Element& generator,
                                 flint::FmpzMatConstRef row) noexcept;

bool set_relation_from_integral_coordinates_and_norm(
        Relation& out,
        const Element& generator,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool set_relation_from_integral_coordinates_and_norm(
        Relation& out,
        const Element& generator,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool set_relation_from_integral_coordinates_and_norm(
        Relation& out,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial = nullptr,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool factor_relation_row_from_integral_coordinates_and_norm(
        Relation& out,
        bool& handled,
        bool& smooth,
        flint::FmpzMatConstRef integral_coordinates,
        flint::FmpqConstRef norm,
        const flint::FmpzPoly* integral_coordinate_polynomial = nullptr,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool commit_relation_generator_from_integral_coordinates(
        Relation& out,
        flint::FmpzMatConstRef integral_coordinates,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

flint::FmpzMatConstRef pending_relation_exponents_ref(
        const Relation& relation) noexcept;

}  // namespace detail

class Relation {
public:
    Relation() noexcept = default;
    explicit Relation(const FactorBase& base) noexcept;
    ~Relation() noexcept;

    Relation(const Relation&) = delete;
    Relation& operator=(const Relation&) = delete;

    Relation(Relation&& other) noexcept;
    Relation& operator=(Relation&& other) noexcept;

    void swap(Relation& other) noexcept;
    void clear() noexcept;
    bool define(const FactorBase& base) noexcept;
    bool set(const Relation& other) noexcept;

    bool is_defined() const noexcept;
    bool is_set() const noexcept;
    const FactorBase* factor_base() const noexcept;
    const Order* parent() const noexcept;
    slong length() const noexcept;

    bool set_generator(const Element& alpha,
                       const DiagnosticsContext* diagnostics = nullptr)
            noexcept;
    // `norm` must be the exact field norm of `alpha`; hot relation-search
    // loops use this to avoid recomputing a norm they already screened.
    bool set_generator_with_norm(
            const Element& alpha,
            flint::FmpqConstRef norm,
            const DiagnosticsContext* diagnostics = nullptr) noexcept;
    bool generator(Element& out) const noexcept;
    bool exponents(flint::FmpzMatRef out) const noexcept;
    // Borrowed exponent-row view; valid while this relation object is alive.
    // Callers should first check `is_set()`.
    flint::FmpzMatConstRef exponents_ref() const noexcept;

private:
    bool set_generator_impl(const Element& alpha,
                            const fmpq* known_norm,
                            const fmpz_mat_struct* known_integral_coordinates,
                            const fmpz_poly_struct* known_integral_polynomial,
                            const DiagnosticsContext* diagnostics) noexcept;

    FactorBase base_;
    Element generator_;
    flint::FmpzMat exponents_{0, 0};
    flint::FmpzMat scratch_exponents_{0, 0};
    bool has_relation_ = false;

    friend bool detail::set_relation_from_known_row(
            Relation& out,
            const FactorBase& base,
            const Element& generator,
            flint::FmpzMatConstRef row) noexcept;
    friend bool detail::set_relation_from_integral_coordinates_and_norm(
            Relation& out,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const DiagnosticsContext* diagnostics) noexcept;
    friend bool detail::set_relation_from_integral_coordinates_and_norm(
            Relation& out,
            const Element& generator,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            const DiagnosticsContext* diagnostics) noexcept;
    friend bool detail::set_relation_from_integral_coordinates_and_norm(
            Relation& out,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            const DiagnosticsContext* diagnostics) noexcept;
    friend bool detail::factor_relation_row_from_integral_coordinates_and_norm(
            Relation& out,
            bool& handled,
            bool& smooth,
            flint::FmpzMatConstRef integral_coordinates,
            flint::FmpqConstRef norm,
            const flint::FmpzPoly* integral_coordinate_polynomial,
            const DiagnosticsContext* diagnostics) noexcept;
    friend bool detail::commit_relation_generator_from_integral_coordinates(
            Relation& out,
            flint::FmpzMatConstRef integral_coordinates,
            const DiagnosticsContext* diagnostics) noexcept;
    friend flint::FmpzMatConstRef detail::pending_relation_exponents_ref(
            const Relation& relation) noexcept;
};

class RelationMatrix {
public:
    RelationMatrix() noexcept = default;
    explicit RelationMatrix(const FactorBase& base) noexcept;
    ~RelationMatrix() noexcept;

    RelationMatrix(const RelationMatrix&) = delete;
    RelationMatrix& operator=(const RelationMatrix&) = delete;

    RelationMatrix(RelationMatrix&& other) noexcept;
    RelationMatrix& operator=(RelationMatrix&& other) noexcept;

    void swap(RelationMatrix& other) noexcept;
    void clear() noexcept;
    bool define(const FactorBase& base) noexcept;
    bool set(const RelationMatrix& other) noexcept;

    bool is_defined() const noexcept;
    const FactorBase* factor_base() const noexcept;
    const Order* parent() const noexcept;
    slong length() const noexcept;
    slong ncols() const noexcept;

    bool append(const Relation& relation) noexcept;
    bool row(flint::FmpzMatRef out, slong index) const noexcept;
    bool row_equal(flint::FmpzMatConstRef row, slong index) const noexcept;
    bool row_first_nonzero(slong& out, slong index) const noexcept;
    bool rows(flint::FmpzMatRef out) const noexcept;
    bool generator(Element& out, slong index) const noexcept;
    const Element* generator_at(slong index) const noexcept;
    bool to_abelian_group(FiniteAbelianGroup& out) const noexcept;

private:
    FactorBase base_;
    fmpz_smat::SparseMat rows_{0};
    std::vector<Element> generators_;

    friend class ClassGroupContext;
};

inline void swap(Relation& left, Relation& right) noexcept {
    left.swap(right);
}

inline void swap(RelationMatrix& left, RelationMatrix& right) noexcept {
    left.swap(right);
}

}  // namespace silex
