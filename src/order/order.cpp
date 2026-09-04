#include <silex/order.hpp>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/nf_elem.h>

#include <silex/archimedean.hpp>
#include <silex/diagnostics.hpp>
#include <silex/embedding.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/flint/fmpz_mod_poly.hpp>
#include <silex/flint/fmpz_mod_poly_factor.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/flint/nf_elem.hpp>
#include <silex/ideal.hpp>
#include <silex/order_element.hpp>

#include <utility>

namespace silex {

namespace detail {

struct OrderData {
    NumberField parent_;
    slong degree_ = 0;
    flint::FmpqMat basis_{0, 0};
    flint::FmpqMat basis_inv_{0, 0};
    mutable flint::FmpzMat trace_mat_{0, 0};
    flint::FmpzMat mul_table_{0, 0};
    mutable flint::ArbMat minkowski_embedding_rows_{0, 0};
    mutable EmbeddingContext minkowski_embedding_context_;
    mutable flint::Fmpz discriminant_;
    bool has_basis_ = false;
    bool has_basis_inv_ = false;
    mutable bool has_trace_mat_ = false;
    bool has_mul_table_ = false;
    mutable slong minkowski_embedding_precision_ = 0;
    mutable bool has_minkowski_embedding_rows_ = false;
    mutable bool has_discriminant_ = false;
    bool maximality_is_known_ = false;
    bool maximal_ = false;
    OrderSpecialization specialization_;
};

}  // namespace detail

namespace {

bool entries_are_integral(const fmpq_mat_t matrix) noexcept {
    for (slong i = 0; i < fmpq_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < fmpq_mat_ncols(matrix); ++j) {
            if (fmpz_is_one(fmpq_mat_entry_den(matrix, i, j)) == 0) {
                return false;
            }
        }
    }
    return true;
}

void set_elem_from_power_row(nf_elem_t out,
                             const fmpq_mat_t basis,
                             slong row,
                             const nf_struct* field) noexcept {
    flint::FmpqPoly polynomial;
    for (slong j = 0; j < fmpq_mat_ncols(basis); ++j) {
        fmpq_poly_set_coeff_fmpq(polynomial.raw(), j,
                                 fmpq_mat_entry(basis, row, j));
    }
    nf_elem_set_fmpq_poly(out, polynomial.raw(), field);
}

void power_row_from_elem(fmpq_mat_t row,
                         const nf_elem_t element,
                         const nf_struct* field,
                         slong degree) noexcept {
    flint::Fmpq coeff;
    fmpq_mat_zero(row);
    for (slong j = 0; j < degree; ++j) {
        nf_elem_get_coeff_fmpq(coeff.raw(), element, j, field);
        fmpq_set(fmpq_mat_entry(row, 0, j), coeff.raw());
    }
}

bool coordinates_are_integral(const fmpq_mat_t power_row,
                              const fmpq_mat_t basis_inv) noexcept {
    flint::FmpqMat coords(1, fmpq_mat_ncols(basis_inv));
    fmpq_mat_mul(coords.raw(), power_row, basis_inv);
    return entries_are_integral(coords.raw());
}

bool basis_contains_one(const fmpq_mat_t basis_inv) noexcept {
    const slong n = fmpq_mat_nrows(basis_inv);
    flint::FmpqMat one(1, n);
    fmpq_one(fmpq_mat_entry(one.raw(), 0, 0));
    return coordinates_are_integral(one.raw(), basis_inv);
}

bool basis_is_closed(const fmpq_mat_t basis,
                     const fmpq_mat_t basis_inv,
                     const nf_struct* field,
                     slong degree) noexcept {
    flint::NfElemVec elems(degree, field);
    if (degree > 0 && elems.data() == nullptr) {
        return false;
    }

    for (slong i = 0; i < degree; ++i) {
        set_elem_from_power_row(elems.data() + i, basis, i, field);
    }

    flint::NfElem product(field);
    flint::FmpqMat row(1, degree);
    for (slong i = 0; i < degree; ++i) {
        for (slong j = i; j < degree; ++j) {
            nf_elem_mul(product.raw(), elems.data() + i, elems.data() + j,
                        field);
            power_row_from_elem(row.raw(), product.raw(), field, degree);
            if (!coordinates_are_integral(row.raw(), basis_inv)) {
                return false;
            }
        }
    }
    return true;
}

bool fmpq_poly_is_monic_integral(const fmpq_poly_t polynomial) noexcept {
    const slong degree = fmpq_poly_degree(polynomial);
    if (degree < 1) {
        return false;
    }

    flint::Fmpq coeff;
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, i);
        if (fmpz_is_one(fmpq_denref(coeff.raw())) == 0) {
            return false;
        }
    }

    fmpq_poly_get_coeff_fmpq(coeff.raw(), polynomial, degree);
    return fmpz_is_one(fmpq_numref(coeff.raw())) != 0;
}

void fmpq_poly_get_fmpz_mod_poly(flint::FmpzModPoly& out,
                                 const fmpq_poly_t polynomial,
                                 const flint::FmpzModCtx& ctx) noexcept {
    flint::Fmpz coeff;
    const slong degree = fmpq_poly_degree(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        fmpq_poly_get_coeff_fmpz(coeff.raw(), polynomial, i);
        fmpz_mod_poly_set_coeff_fmpz(out.raw(), i, coeff.raw(), ctx.raw());
    }
}

void squarefree_kernel(flint::FmpzModPoly& out,
                       const flint::FmpzModPolyFactor& factorization,
                       const flint::FmpzModCtx& ctx) noexcept {
    fmpz_mod_poly_one(out.raw(), ctx.raw());
    for (slong i = 0; i < factorization.raw()->num; ++i) {
        fmpz_mod_poly_mul(out.raw(), out.raw(),
                          factorization.raw()->poly + i, ctx.raw());
    }
}

bool element_set_mod_poly(Element& out,
                          const flint::FmpzModPoly& polynomial,
                          const flint::FmpzModCtx& ctx) noexcept {
    flint::FmpqPoly rational;
    flint::Fmpz coeff;
    const slong degree = fmpz_mod_poly_degree(polynomial.raw(), ctx.raw());
    for (slong i = 0; i <= degree; ++i) {
        fmpz_mod_poly_get_coeff_fmpz(coeff.raw(), polynomial.raw(), i, ctx.raw());
        fmpq_poly_set_coeff_fmpz(rational.raw(), i, coeff.raw());
    }
    return out.set_fmpq_poly(flint::FmpqPolyConstRef(rational));
}

bool order_element_set_fmpz(OrderElement& out,
                            const Order& order,
                            flint::FmpzConstRef value) noexcept {
    Element element(*order.parent());
    flint::FmpqPoly polynomial;
    fmpq_poly_set_coeff_fmpz(polynomial.raw(), 0, value.raw());
    return element.set_fmpq_poly(flint::FmpqPolyConstRef(polynomial)) &&
           out.set_element(element);
}

bool order_one_mod_vec(flint::FmpzVec& one,
                       const fmpq_mat_t basis_inv,
                       const flint::FmpzModCtx& ctx) noexcept {
    const slong n = one.length();
    flint::FmpqMat power(1, n);
    flint::FmpqMat coords(1, n);

    fmpq_one(fmpq_mat_entry(power.raw(), 0, 0));
    fmpq_mat_mul(coords.raw(), power.raw(), basis_inv);

    for (slong j = 0; j < n; ++j) {
        if (fmpz_is_one(fmpq_mat_entry_den(coords.raw(), 0, j)) == 0) {
            return false;
        }

        fmpz_mod_set_fmpz(one.data() + j,
                          fmpq_mat_entry_num(coords.raw(), 0, j), ctx.raw());
    }
    return true;
}

bool multiplication_table_from_basis(flint::FmpzMat& out,
                                     const fmpq_mat_t basis,
                                     const fmpq_mat_t basis_inv,
                                     const nf_struct* field,
                                     slong degree) noexcept {
    flint::NfElemVec elems(degree, field);
    if (degree > 0 && elems.data() == nullptr) {
        return false;
    }

    for (slong i = 0; i < degree; ++i) {
        set_elem_from_power_row(elems.data() + i, basis, i, field);
    }

    flint::NfElem product(field);
    flint::FmpqMat power(1, degree);
    flint::FmpqMat coords(1, degree);
    flint::FmpzMat table(degree * degree, degree);

    for (slong i = 0; i < degree; ++i) {
        for (slong j = i; j < degree; ++j) {
            nf_elem_mul(product.raw(), elems.data() + i, elems.data() + j,
                        field);
            power_row_from_elem(power.raw(), product.raw(), field, degree);
            fmpq_mat_mul(coords.raw(), power.raw(), basis_inv);
            if (!entries_are_integral(coords.raw())) {
                return false;
            }
            for (slong k = 0; k < degree; ++k) {
                fmpz_set(fmpz_mat_entry(table.raw(), i * degree + j, k),
                         fmpq_mat_entry_num(coords.raw(), 0, k));
                if (i != j) {
                    fmpz_set(fmpz_mat_entry(table.raw(), j * degree + i, k),
                             fmpq_mat_entry_num(coords.raw(), 0, k));
                }
            }
        }
    }

    out.swap(table);
    return true;
}

void mod_vec_set(flint::FmpzVec& dest,
                 const flint::FmpzVec& src) noexcept {
    for (slong i = 0; i < dest.length(); ++i) {
        fmpz_set(dest.data() + i, src.data() + i);
    }
}

void mod_vec_mul(flint::FmpzVec& c,
                 const flint::FmpzVec& a,
                 const flint::FmpzVec& b,
                 const flint::FmpzMat& mul_table,
                 const flint::FmpzModCtx& ctx) noexcept {
    const slong n = c.length();
    flint::FmpzVec out(n);
    flint::Fmpz ab;
    flint::Fmpz coeff;
    flint::Fmpz term;

    for (slong i = 0; i < n; ++i) {
        if (fmpz_is_zero(a.data() + i) != 0) {
            continue;
        }

        for (slong j = 0; j < n; ++j) {
            if (fmpz_is_zero(b.data() + j) != 0) {
                continue;
            }

            fmpz_mod_mul(ab.raw(), a.data() + i, b.data() + j, ctx.raw());

            for (slong k = 0; k < n; ++k) {
                fmpz_mod_set_fmpz(coeff.raw(),
                                  fmpz_mat_entry(mul_table.raw(), i * n + j, k),
                                  ctx.raw());
                fmpz_mod_mul(term.raw(), ab.raw(), coeff.raw(), ctx.raw());
                fmpz_add(out.data() + k, out.data() + k, term.raw());
                fmpz_mod_set_fmpz(out.data() + k, out.data() + k, ctx.raw());
            }
        }
    }

    mod_vec_set(c, out);
}

void mod_vec_pow_fmpz(flint::FmpzVec& res,
                      const flint::FmpzVec& base,
                      flint::FmpzConstRef exponent,
                      const flint::FmpzVec& one,
                      const flint::FmpzMat& mul_table,
                      const flint::FmpzModCtx& ctx) noexcept {
    const slong n = res.length();
    flint::FmpzVec acc(n);
    flint::FmpzVec b(n);
    flint::Fmpz e;

    fmpz_set(e.raw(), exponent.raw());
    mod_vec_set(acc, one);
    mod_vec_set(b, base);

    while (fmpz_is_zero(e.raw()) == 0) {
        if (fmpz_is_odd(e.raw()) != 0) {
            mod_vec_mul(acc, acc, b, mul_table, ctx);
        }

        fmpz_fdiv_q_2exp(e.raw(), e.raw(), 1);

        if (fmpz_is_zero(e.raw()) == 0) {
            mod_vec_mul(b, b, b, mul_table, ctx);
        }
    }

    mod_vec_set(res, acc);
}

bool discriminant_square_index_ok(Order& suborder,
                                  Order& overorder,
                                  flint::FmpzConstRef index) noexcept {
    flint::Fmpz sub_disc;
    flint::Fmpz over_disc;
    flint::Fmpz rhs;

    if (!suborder.discriminant(flint::FmpzRef(sub_disc)) ||
        !overorder.discriminant(flint::FmpzRef(over_disc))) {
        return false;
    }

    fmpz_mul(rhs.raw(), index.raw(), index.raw());
    fmpz_mul(rhs.raw(), rhs.raw(), over_disc.raw());
    return fmpz_equal(sub_disc.raw(), rhs.raw()) != 0;
}

void power_trace_matrix(fmpq_mat_t out,
                        const nf_struct* field,
                        slong degree) noexcept {
    const slong length = 2 * degree - 1;
    flint::NfElemVec powers(length, field);
    flint::Fmpq trace;

    nf_elem_one(powers.data(), field);
    for (slong i = 1; i < length; ++i) {
        nf_elem_mul_gen(powers.data() + i, powers.data() + i - 1, field);
    }

    for (slong i = 0; i < degree; ++i) {
        for (slong j = 0; j < degree; ++j) {
            nf_elem_trace(trace.raw(), powers.data() + i + j, field);
            fmpq_set(fmpq_mat_entry(out, i, j), trace.raw());
        }
    }
}

bool output_shape(flint::FmpqMatRef matrix, slong rows, slong cols) noexcept {
    return fmpq_mat_nrows(matrix.raw()) == rows &&
           fmpq_mat_ncols(matrix.raw()) == cols;
}

bool output_shape(flint::FmpzMatRef matrix, slong rows, slong cols) noexcept {
    return fmpz_mat_nrows(matrix.raw()) == rows &&
           fmpz_mat_ncols(matrix.raw()) == cols;
}

void quadratic_field_discriminant(flint::FmpzRef out,
                                  flint::FmpzConstRef radicand) noexcept {
    if (flint::fmpz_fdiv_ui(radicand, 4) == 1) {
        flint::fmpz_set(out, radicand);
    } else {
        flint::fmpz_mul_ui(out, radicand, 4);
    }
}

void quadratic_basis(flint::FmpqMatRef out,
                     flint::FmpzConstRef radicand,
                     flint::FmpzConstRef conductor) noexcept {
    fmpq_mat_zero(out.raw());
    fmpq_one(fmpq_mat_entry(out.raw(), 0, 0));

    if (flint::fmpz_fdiv_ui(radicand, 4) == 1) {
        fmpq_set_fmpz(fmpq_mat_entry(out.raw(), 1, 0), conductor.raw());
        fmpq_div_2exp(fmpq_mat_entry(out.raw(), 1, 0),
                      fmpq_mat_entry(out.raw(), 1, 0), 1);
        fmpq_set_fmpz(fmpq_mat_entry(out.raw(), 1, 1), conductor.raw());
        fmpq_div_2exp(fmpq_mat_entry(out.raw(), 1, 1),
                      fmpq_mat_entry(out.raw(), 1, 1), 1);
    } else {
        fmpq_set_fmpz(fmpq_mat_entry(out.raw(), 1, 1), conductor.raw());
    }
}

}  // namespace

namespace {

bool detect_quadratic_conductor(flint::FmpzRef conductor,
                                const NumberField& parent,
                                flint::FmpqMatConstRef basis) noexcept {
    if (parent.backend_kind() != NumberFieldBackendKind::quadratic ||
        parent.degree() != 2 ||
        fmpq_mat_nrows(basis.raw()) != 2 ||
        fmpq_mat_ncols(basis.raw()) != 2) {
        return false;
    }

    flint::Fmpz radicand;
    if (!parent.quadratic_radicand(flint::FmpzRef(radicand))) {
        return false;
    }

    if (fmpq_mat_is_one(basis.raw()) != 0) {
        if (fmpz_fdiv_ui(radicand.raw(), 4) == 1) {
            fmpz_set_ui(conductor.raw(), 2);
        } else {
            fmpz_one(conductor.raw());
        }
        return true;
    }

    if (fmpq_equal_si(fmpq_mat_entry(basis.raw(), 0, 0), 1) == 0 ||
        fmpq_is_zero(fmpq_mat_entry(basis.raw(), 0, 1)) == 0) {
        return false;
    }

    if (fmpz_fdiv_ui(radicand.raw(), 4) == 1) {
        if (fmpq_equal(fmpq_mat_entry(basis.raw(), 1, 0),
                       fmpq_mat_entry(basis.raw(), 1, 1)) == 0) {
            return false;
        }
        flint::Fmpq doubled;
        fmpq_mul_2exp(doubled.raw(), fmpq_mat_entry(basis.raw(), 1, 1), 1);
        const bool ok = fmpz_is_one(fmpq_denref(doubled.raw())) != 0 &&
                        fmpz_sgn(fmpq_numref(doubled.raw())) > 0;
        if (ok) {
            fmpz_set(conductor.raw(), fmpq_numref(doubled.raw()));
        }
        return ok;
    }

    if (fmpq_is_zero(fmpq_mat_entry(basis.raw(), 1, 0)) == 0 ||
        fmpz_is_one(fmpq_mat_entry_den(basis.raw(), 1, 1)) == 0 ||
        fmpz_sgn(fmpq_mat_entry_num(basis.raw(), 1, 1)) <= 0) {
        return false;
    }

    fmpz_set(conductor.raw(), fmpq_mat_entry_num(basis.raw(), 1, 1));
    return true;
}

}  // namespace

namespace detail {

bool order_minkowski_embedding_rows(
        flint::ArbMatRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    if (!order.is_defined() || !order.data_->has_basis_ || precision <= 0) {
        return false;
    }

    const OrderData& data = *order.data_;
    const slong degree = data.degree_;
    if (degree <= 0 ||
        flint::arb_mat_nrows_value(out) != degree ||
        flint::arb_mat_ncols_value(out) != degree) {
        return false;
    }

    if (data.has_minkowski_embedding_rows_ &&
        data.minkowski_embedding_precision_ >= precision &&
        flint::arb_mat_nrows_value(data.minkowski_embedding_rows_) == degree &&
        flint::arb_mat_ncols_value(data.minkowski_embedding_rows_) == degree) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "order.minkowski_embedding_rows.cache_hit");
        flint::arb_mat_set(out,
                           flint::ArbMatConstRef(data.minkowski_embedding_rows_));
        return true;
    }
    if (!data.has_minkowski_embedding_rows_) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "order.minkowski_embedding_rows.cache_empty");
    } else if (data.minkowski_embedding_precision_ < precision) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "order.minkowski_embedding_rows.cache_precision_miss");
    } else {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::class_group,
                "order.minkowski_embedding_rows.cache_shape_miss");
    }

    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "order.minkowski_embedding_rows.build");
    flint::ArbMat rows(degree, degree);
    EmbeddingContext& embeddings = data.minkowski_embedding_context_;
    if (!embeddings.is_defined() && !embeddings.define(data.parent_)) {
        return false;
    }
    OrderElement order_element(order);
    Element element(data.parent_);
    flint::FmpzMat coordinates(1, degree);
    flint::ArbMat row(1, degree);
    if (!embeddings.is_defined() || !order_element.is_defined() ||
        !element.is_defined()) {
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "order.minkowski_embedding_rows.refine");
        if (!embeddings.refine(precision, diagnostics)) {
            return false;
        }
    }

    // Cache the order-basis weighted embedding rows. Ideal bases are still
    // multiplied into them by the caller.
    for (slong i = 0; i < degree; ++i) {
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "order.minkowski_embedding_rows.element");
            flint::fmpz_mat_zero(flint::FmpzMatRef(coordinates));
            flint::fmpz_one(flint::fmpz_mat_entry(
                    flint::FmpzMatRef(coordinates), 0, i));
            if (!order_element.set_coordinates(
                        flint::FmpzMatConstRef(coordinates)) ||
                !order_element.get_element(element)) {
                return false;
            }
        }

        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "order.minkowski_embedding_rows.row");
            if (!minkowski_embedding(flint::ArbMatRef(row), embeddings, element,
                                     MinkowskiEmbeddingMode::weighted,
                                     precision)) {
                return false;
            }
        }

        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "order.minkowski_embedding_rows.copy");
            for (slong j = 0; j < degree; ++j) {
                flint::arb_set(flint::arb_mat_entry_ref(rows, i, j),
                               flint::arb_mat_entry_ref(row, 0, j));
            }
        }
    }

    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::class_group,
                "order.minkowski_embedding_rows.store");
        data.minkowski_embedding_rows_ = std::move(rows);
        data.minkowski_embedding_precision_ = precision;
        data.has_minkowski_embedding_rows_ = true;
        flint::arb_mat_set(
                out, flint::ArbMatConstRef(data.minkowski_embedding_rows_));
    }
    return true;
}

}  // namespace detail

bool Order::set_quadratic_metadata(flint::FmpzConstRef conductor) noexcept {
    if (!is_defined() ||
        data_->parent_.backend_kind() != NumberFieldBackendKind::quadratic ||
        data_->degree_ != 2 || flint::fmpz_sgn(conductor) <= 0) {
        return false;
    }

    flint::Fmpz radicand;
    flint::Fmpz field_discriminant;
    if (!data_->parent_.quadratic_radicand(flint::FmpzRef(radicand))) {
        return false;
    }

    quadratic_field_discriminant(flint::FmpzRef(field_discriminant),
                                 flint::FmpzConstRef(radicand));
    flint::fmpz_mul(flint::FmpzRef(data_->discriminant_),
                    flint::FmpzConstRef(field_discriminant), conductor);
    flint::fmpz_mul(flint::FmpzRef(data_->discriminant_),
                    flint::FmpzConstRef(data_->discriminant_), conductor);
    data_->has_discriminant_ = true;
    data_->maximality_is_known_ = true;
    data_->maximal_ = flint::fmpz_is_one(conductor);
    data_->specialization_ = QuadraticOrderData(conductor);
    return true;
}

bool Order::set_quadratic_order(flint::FmpzConstRef conductor) noexcept {
    if (!is_defined() ||
        data_->parent_.backend_kind() != NumberFieldBackendKind::quadratic ||
        data_->degree_ != 2 || flint::fmpz_sgn(conductor) <= 0) {
        return false;
    }

    flint::Fmpz radicand;
    if (!data_->parent_.quadratic_radicand(flint::FmpzRef(radicand))) {
        return false;
    }

    Order candidate(data_->parent_);
    flint::FmpqMat basis(2, 2);
    quadratic_basis(flint::FmpqMatRef(basis), flint::FmpzConstRef(radicand),
                    conductor);
    if (!candidate.set_basis(flint::FmpqMatConstRef(basis)) ||
        !candidate.set_quadratic_metadata(conductor)) {
        return false;
    }

    swap(candidate);
    return true;
}

Order::Order(const NumberField& parent) noexcept {
    define(parent);
}

Order::~Order() noexcept = default;

Order::Order(const Order& other) noexcept = default;

Order& Order::operator=(const Order& other) noexcept = default;

Order::Order(Order&& other) noexcept {
    swap(other);
}

Order& Order::operator=(Order&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void Order::swap(Order& other) noexcept {
    data_.swap(other.data_);
}

void Order::clear() noexcept {
    data_.reset();
}

bool Order::define(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    Order next;
    next.data_ = std::make_shared<detail::OrderData>();
    next.data_->parent_ = parent;
    next.data_->degree_ = parent.degree();
    next.data_->basis_ = flint::FmpqMat(next.data_->degree_, next.data_->degree_);
    next.data_->basis_inv_ = flint::FmpqMat(next.data_->degree_, next.data_->degree_);
    next.data_->trace_mat_ = flint::FmpzMat(next.data_->degree_, next.data_->degree_);
    next.data_->mul_table_ = flint::FmpzMat(next.data_->degree_ * next.data_->degree_, next.data_->degree_);

    swap(next);
    return true;
}

bool Order::define_equation_order(const NumberField& parent) noexcept {
    Order next = Order::equation_order(parent);
    if (!next.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

Order Order::equation_order(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return {};
    }

    Order next(parent);
    if (!next.is_defined()) {
        return {};
    }

    flint::FmpqMat identity(parent.degree(), parent.degree());
    fmpq_mat_one(identity.raw());
    if (!next.set_basis(flint::FmpqMatConstRef(identity))) {
        return {};
    }

    next.data_->specialization_ = EquationOrderData{};
    if (next.data_->degree_ == 1) {
        next.data_->maximality_is_known_ = true;
        next.data_->maximal_ = true;
    } else if (next.data_->degree_ == 2 &&
               parent.backend_kind() == NumberFieldBackendKind::quadratic) {
        flint::Fmpz radicand;
        flint::Fmpz conductor;
        if (!parent.quadratic_radicand(flint::FmpzRef(radicand))) {
            return {};
        }
        if (fmpz_fdiv_ui(radicand.raw(), 4) == 1) {
            fmpz_set_ui(conductor.raw(), 2);
        } else {
            fmpz_one(conductor.raw());
        }
        if (!next.set_quadratic_metadata(flint::FmpzConstRef(conductor))) {
            return {};
        }
    }

    return next;
}

Order Order::from_basis(const NumberField& parent,
                        flint::FmpqMatConstRef basis) noexcept {
    Order next(parent);
    if (!next.is_defined() || !next.set_basis(basis)) {
        return {};
    }
    return next;
}

Order Order::quadratic_order(const NumberField& parent,
                             flint::FmpzConstRef conductor) noexcept {
    Order next(parent);
    if (!next.is_defined() || !next.set_quadratic_order(conductor)) {
        return {};
    }
    return next;
}

bool Order::set(const Order& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    Order copy(other.data_->parent_);
    if (!copy.is_defined()) {
        return false;
    }

    fmpq_mat_set(copy.data_->basis_.raw(), other.data_->basis_.raw());
    fmpq_mat_set(copy.data_->basis_inv_.raw(), other.data_->basis_inv_.raw());
    fmpz_mat_set(copy.data_->trace_mat_.raw(), other.data_->trace_mat_.raw());
    fmpz_mat_set(copy.data_->mul_table_.raw(), other.data_->mul_table_.raw());
    fmpz_set(copy.data_->discriminant_.raw(), other.data_->discriminant_.raw());

    copy.data_->has_basis_ = other.data_->has_basis_;
    copy.data_->has_basis_inv_ = other.data_->has_basis_inv_;
    copy.data_->has_trace_mat_ = other.data_->has_trace_mat_;
    copy.data_->has_mul_table_ = other.data_->has_mul_table_;
    copy.data_->has_discriminant_ = other.data_->has_discriminant_;
    copy.data_->maximality_is_known_ = other.data_->maximality_is_known_;
    copy.data_->maximal_ = other.data_->maximal_;
    copy.data_->specialization_ = other.data_->specialization_;

    swap(copy);
    return true;
}

bool Order::is_defined() const noexcept {
    return data_ != nullptr && data_->parent_.is_defined();
}

bool Order::has_same_data(const Order& other) const noexcept {
    return data_ != nullptr && data_ == other.data_;
}

const NumberField* Order::parent() const noexcept {
    return is_defined() ? &data_->parent_ : nullptr;
}

slong Order::degree() const noexcept {
    return is_defined() ? data_->degree_ : 0;
}

bool Order::has_basis() const noexcept {
    return is_defined() && data_->has_basis_;
}

bool Order::is_equation_order() const noexcept {
    return is_defined() && data_->has_basis_ && fmpq_mat_is_one(data_->basis_.raw()) != 0;
}

bool Order::maximality_known() const noexcept {
    return is_defined() && data_->maximality_is_known_;
}

bool Order::is_maximal() const noexcept {
    return is_defined() && data_->maximality_is_known_ && data_->maximal_;
}

void Order::set_maximality(bool is_maximal_value) noexcept {
    if (!is_defined()) {
        return;
    }
    data_->maximality_is_known_ = true;
    data_->maximal_ = is_maximal_value;
}

void Order::clear_maximality() noexcept {
    if (!is_defined()) {
        return;
    }
    data_->maximality_is_known_ = false;
    data_->maximal_ = false;
}

bool Order::set_basis(flint::FmpqMatConstRef basis) noexcept {
    if (!is_defined() ||
        fmpq_mat_nrows(basis.raw()) != data_->degree_ ||
        fmpq_mat_ncols(basis.raw()) != data_->degree_) {
        return false;
    }

    flint::FmpqMat basis_inv(data_->degree_, data_->degree_);
    if (fmpq_mat_inv(basis_inv.raw(), basis.raw()) == 0) {
        return false;
    }
    if (!basis_contains_one(basis_inv.raw())) {
        return false;
    }
    if (!basis_is_closed(basis.raw(), basis_inv.raw(),
                         data_->parent_.raw_flint_field(), data_->degree_)) {
        return false;
    }

    Order next(data_->parent_);
    if (!next.is_defined()) {
        return false;
    }
    fmpq_mat_set(next.data_->basis_.raw(), basis.raw());
    fmpq_mat_set(next.data_->basis_inv_.raw(), basis_inv.raw());
    next.data_->has_basis_ = true;
    next.data_->has_basis_inv_ = true;
    if (fmpq_mat_is_one(basis.raw()) != 0) {
        next.data_->specialization_ = EquationOrderData{};
        if (next.data_->degree_ == 1) {
            next.data_->maximality_is_known_ = true;
            next.data_->maximal_ = true;
        }
    }
    flint::Fmpz conductor;
    if (detect_quadratic_conductor(flint::FmpzRef(conductor),
                                   data_->parent_, basis)) {
        if (!next.set_quadratic_metadata(flint::FmpzConstRef(conductor))) {
            return false;
        }
    }

    swap(next);
    return true;
}

bool Order::get_basis(flint::FmpqMatRef out) const noexcept {
    if (!is_defined() || !data_->has_basis_ ||
        !output_shape(out, data_->degree_, data_->degree_)) {
        return false;
    }
    fmpq_mat_set(out.raw(), data_->basis_.raw());
    return true;
}

std::optional<flint::FmpqMat> Order::basis() const noexcept {
    if (!is_defined() || !data_->has_basis_) {
        return std::nullopt;
    }
    flint::FmpqMat out(data_->degree_, data_->degree_);
    if (!get_basis(flint::FmpqMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool Order::coordinates(flint::FmpqMatRef out,
                        const Element& element) const noexcept {
    if (!is_defined() || !data_->has_basis_ || !data_->has_basis_inv_ ||
        !element.has_parent(data_->parent_) ||
        !output_shape(out, 1, data_->degree_)) {
        return false;
    }

    flint::FmpqMat power(1, data_->degree_);
    flint::FmpqMat coords(1, data_->degree_);
    power_row_from_elem(power.raw(), element.raw_flint_element(),
                        data_->parent_.raw_flint_field(), data_->degree_);
    fmpq_mat_mul(coords.raw(), power.raw(), data_->basis_inv_.raw());
    fmpq_mat_set(out.raw(), coords.raw());
    return true;
}

std::optional<flint::FmpqMat> Order::coordinates(
        const Element& element) const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpqMat out(1, data_->degree_);
    if (!coordinates(flint::FmpqMatRef(out), element)) {
        return std::nullopt;
    }
    return out;
}

bool Order::contains(const Element& element) const noexcept {
    if (!is_defined() || !data_->has_basis_ || !element.has_parent(data_->parent_)) {
        return false;
    }

    flint::FmpqMat coords(1, data_->degree_);
    return coordinates(flint::FmpqMatRef(coords), element) &&
           entries_are_integral(coords.raw());
}

bool Order::trace_matrix(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !data_->has_basis_ ||
        !output_shape(out, data_->degree_, data_->degree_)) {
        return false;
    }

    if (!data_->has_trace_mat_) {
        flint::FmpqMat power_trace(data_->degree_, data_->degree_);
        flint::FmpqMat trace_basis_left(data_->degree_, data_->degree_);
        flint::FmpqMat trace_basis(data_->degree_, data_->degree_);
        flint::FmpqMat basis_transpose(data_->degree_, data_->degree_);
        flint::FmpzMat trace_integer(data_->degree_, data_->degree_);

        power_trace_matrix(power_trace.raw(), data_->parent_.raw_flint_field(),
                           data_->degree_);
        fmpq_mat_mul(trace_basis_left.raw(), data_->basis_.raw(),
                     power_trace.raw());
        fmpq_mat_transpose(basis_transpose.raw(), data_->basis_.raw());
        fmpq_mat_mul(trace_basis.raw(), trace_basis_left.raw(),
                     basis_transpose.raw());
        if (fmpq_mat_get_fmpz_mat(trace_integer.raw(), trace_basis.raw()) == 0) {
            return false;
        }

        fmpz_mat_set(data_->trace_mat_.raw(), trace_integer.raw());
        data_->has_trace_mat_ = true;
    }

    fmpz_mat_set(out.raw(), data_->trace_mat_.raw());
    return true;
}

bool Order::discriminant(flint::FmpzRef out) const noexcept {
    if (!is_defined() || !data_->has_basis_) {
        return false;
    }

    if (!data_->has_discriminant_) {
        flint::FmpzMat trace(data_->degree_, data_->degree_);
        if (!trace_matrix(flint::FmpzMatRef(trace))) {
            return false;
        }
        fmpz_mat_det(data_->discriminant_.raw(), data_->trace_mat_.raw());
        data_->has_discriminant_ = true;
    }

    fmpz_set(out.raw(), data_->discriminant_.raw());
    return true;
}

bool Order::multiplication_table(flint::FmpzMatRef out) noexcept {
    if (!is_defined() || !data_->has_basis_ ||
        !output_shape(out, data_->degree_ * data_->degree_, data_->degree_)) {
        return false;
    }

    if (!data_->has_mul_table_) {
        flint::FmpzMat table(data_->degree_ * data_->degree_, data_->degree_);
        if (!multiplication_table_from_basis(table, data_->basis_.raw(),
                                             data_->basis_inv_.raw(),
                                             data_->parent_.raw_flint_field(),
                                             data_->degree_)) {
            return false;
        }

        fmpz_mat_set(data_->mul_table_.raw(), table.raw());
        data_->has_mul_table_ = true;
    }

    fmpz_mat_set(out.raw(), data_->mul_table_.raw());
    return true;
}

bool Order::index_in(flint::FmpzRef out, const Order& overorder) const noexcept {
    return order_index(out, *this, overorder);
}

bool Order::p_radical(Ideal& out, flint::FmpzConstRef prime) const noexcept {
    if (!is_defined() || !data_->has_basis_ || out.parent() == nullptr ||
        !out.parent()->has_same_data(*this) ||
        fmpz_is_prime(prime.raw()) == 0) {
        return false;
    }

    const slong n = data_->degree_;

    if (!is_equation_order()) {
        flint::FmpzModCtx ctx(prime.raw());
        flint::FmpzModMat phi(n, n, ctx);
        flint::FmpzModMat kernel(n, n, ctx);
        flint::FmpzMat mul_table(n * n, n);
        flint::FmpzMat rows(n, n);
        flint::Fmpz q;
        flint::Fmpz entry;
        flint::FmpzVec one(n);
        flint::FmpzVec base(n);
        flint::FmpzVec image(n);
        Ideal candidate(*this);

        if (!order_one_mod_vec(one, data_->basis_inv_.raw(), ctx) ||
            !multiplication_table_from_basis(mul_table, data_->basis_.raw(),
                                             data_->basis_inv_.raw(),
                                             data_->parent_.raw_flint_field(),
                                             n)) {
            return false;
        }

        fmpz_set(q.raw(), prime.raw());
        while (fmpz_cmp_si(q.raw(), n) <= 0) {
            fmpz_mul(q.raw(), q.raw(), prime.raw());
        }

        for (slong j = 0; j < n; ++j) {
            _fmpz_vec_zero(base.data(), n);
            fmpz_one(base.data() + j);
            mod_vec_pow_fmpz(image, base, flint::FmpzConstRef(q), one,
                             mul_table, ctx);

            for (slong i = 0; i < n; ++i) {
                fmpz_set(fmpz_mod_mat_entry(phi.raw(), i, j), image.data() + i);
            }
        }

        const slong nullity = fmpz_mod_mat_nullspace(kernel.raw(),
                                                     phi.raw(), ctx.raw());
        rows = flint::FmpzMat(n + nullity, n);

        for (slong i = 0; i < n; ++i) {
            fmpz_set(fmpz_mat_entry(rows.raw(), i, i), prime.raw());
        }

        for (slong j = 0; j < nullity; ++j) {
            for (slong i = 0; i < n; ++i) {
                fmpz_mod_mat_get_entry(entry.raw(), kernel.raw(), i, j,
                                       ctx.raw());
                fmpz_set(fmpz_mat_entry(rows.raw(), n + j, i), entry.raw());
            }
        }

        return candidate.set_rows(flint::FmpzMatConstRef(rows)) &&
               out.set(candidate);
    }

    if (!fmpq_poly_is_monic_integral(data_->parent_.raw_flint_field()->pol)) {
        return false;
    }

    Ideal p_ideal(*this);
    Ideal generator_ideal(*this);
    Ideal candidate(*this);
    OrderElement p_element(*this);
    OrderElement generator(*this);
    Element generator_element(data_->parent_);
    if (!p_ideal.is_defined() || !generator_ideal.is_defined() ||
        !candidate.is_defined() || !p_element.is_defined() ||
        !generator.is_defined() ||
        !order_element_set_fmpz(p_element, *this, prime) ||
        !p_ideal.set_principal(p_element)) {
        return false;
    }

    if (n == 1) {
        return candidate.set(p_ideal) && out.set(candidate);
    }

    flint::FmpzModCtx ctx(prime.raw());
    flint::FmpzModPoly reduced(ctx);
    flint::FmpzModPoly kernel(ctx);
    flint::FmpzModPolyFactor factorization(ctx);

    fmpq_poly_get_fmpz_mod_poly(reduced, data_->parent_.raw_flint_field()->pol, ctx);
    fmpz_mod_poly_factor(factorization.raw(), reduced.raw(), ctx.raw());
    squarefree_kernel(kernel, factorization, ctx);

    if (fmpz_mod_poly_is_zero(kernel.raw(), ctx.raw()) != 0) {
        return candidate.set(p_ideal) && out.set(candidate);
    }

    if (!element_set_mod_poly(generator_element, kernel, ctx) ||
        !generator.set_element(generator_element) ||
        !generator_ideal.set_principal(generator) ||
        !candidate.add(p_ideal, generator_ideal)) {
        return false;
    }

    return out.set(candidate);
}

bool Order::pmaximal_overorder(const Order& input,
                               flint::FmpzConstRef prime) noexcept {
    if (!is_defined() || !input.has_basis() || this == &input ||
        !data_->parent_.has_same_data(input.data_->parent_) ||
        fmpz_is_prime(prime.raw()) == 0) {
        return false;
    }

    Order candidate(data_->parent_);
    if (data_->degree_ == 1) {
        if (!candidate.set(input)) {
            return false;
        }
        candidate.set_maximality(true);
        swap(candidate);
        return true;
    }

    if (const auto* quadratic =
            std::get_if<QuadraticOrderData>(&input.data_->specialization_)) {
        flint::Fmpz conductor;
        fmpz_remove(conductor.raw(), quadratic->conductor.raw(),
                    prime.raw());
        if (!candidate.set_quadratic_order(flint::FmpzConstRef(conductor))) {
            return false;
        }

        swap(candidate);
        return true;
    }

    Order current(data_->parent_);
    Order next(data_->parent_);
    flint::Fmpz p_squared;
    flint::Fmpz disc;
    flint::Fmpz index;

    if (!current.set(input)) {
        return false;
    }

    fmpz_mul(p_squared.raw(), prime.raw(), prime.raw());

    for (;;) {
        if (!current.discriminant(flint::FmpzRef(disc))) {
            return false;
        }

        if (fmpz_divisible(disc.raw(), p_squared.raw()) == 0) {
            return candidate.set(current) && (swap(candidate), true);
        }

        Ideal radical(current);
        if (!current.p_radical(radical, prime) ||
            !radical.multiplier_ring(next) ||
            !order_index(flint::FmpzRef(index), current, next)) {
            return false;
        }

        if (fmpz_is_one(index.raw()) != 0) {
            return candidate.set(current) && (swap(candidate), true);
        }

        if (fmpz_divisible(index.raw(), prime.raw()) == 0 ||
            !discriminant_square_index_ok(current, next,
                                          flint::FmpzConstRef(index))) {
            return false;
        }

        current.swap(next);
    }
}

bool Order::maximal_order(const Order& input) noexcept {
    if (!is_defined() || !input.has_basis() || this == &input ||
        !data_->parent_.has_same_data(input.data_->parent_)) {
        return false;
    }

    Order candidate(data_->parent_);
    if (data_->degree_ == 1) {
        if (!candidate.set(input)) {
            return false;
        }
        candidate.set_maximality(true);
        swap(candidate);
        return true;
    }

    if (std::holds_alternative<QuadraticOrderData>(input.data_->specialization_)) {
        flint::Fmpz one;
        fmpz_one(one.raw());
        if (!candidate.set_quadratic_order(flint::FmpzConstRef(one))) {
            return false;
        }

        swap(candidate);
        return true;
    }

    Order current(data_->parent_);
    Order next(data_->parent_);
    flint::Fmpz disc;
    flint::FmpzFactor factorization;

    if (!current.set(input) ||
        !current.discriminant(flint::FmpzRef(disc)) ||
        fmpz_is_zero(disc.raw()) != 0) {
        return false;
    }

    fmpz_abs(disc.raw(), disc.raw());
    fmpz_factor(factorization.raw(), disc.raw());

    for (slong i = 0; i < factorization.raw()->num; ++i) {
        if (factorization.raw()->exp[i] < 2) {
            continue;
        }

        if (!next.pmaximal_overorder(
                    current,
                    flint::FmpzConstRef(factorization.raw()->p + i))) {
            return false;
        }

        current.swap(next);
    }

    if (!candidate.set(current)) {
        return false;
    }
    candidate.set_maximality(true);
    swap(candidate);
    return true;
}

bool Order::quadratic_conductor(flint::FmpzRef out) const noexcept {
    const auto* quadratic = std::get_if<QuadraticOrderData>(&data_->specialization_);
    if (!is_defined() || quadratic == nullptr) {
        return false;
    }
    fmpz_set(out.raw(), quadratic->conductor.raw());
    return true;
}

bool order_index(flint::FmpzRef out,
                 const Order& suborder,
                 const Order& overorder) noexcept {
    if (!suborder.has_basis() || !overorder.has_basis() ||
        !suborder.data_->parent_.is_defined() ||
        !suborder.data_->parent_.has_same_data(overorder.data_->parent_)) {
        return false;
    }

    const slong n = suborder.degree();
    flint::FmpqMat conversion(n, n);
    flint::FmpzMat integer_conversion(n, n);
    flint::Fmpz det;

    fmpq_mat_mul(conversion.raw(), suborder.data_->basis_.raw(),
                 overorder.data_->basis_inv_.raw());
    if (fmpq_mat_get_fmpz_mat(integer_conversion.raw(), conversion.raw()) == 0) {
        return false;
    }

    fmpz_mat_det(det.raw(), integer_conversion.raw());
    if (fmpz_is_zero(det.raw()) != 0) {
        return false;
    }

    fmpz_abs(out.raw(), det.raw());
    return true;
}

}  // namespace silex
