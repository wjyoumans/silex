#include "class_group_internal.hpp"

#include <utility>

namespace silex::detail {

namespace {

bool compact_multiply_power_fmpz(FactoredElement& accumulator,
                                 const FactoredElement& base,
                                 flint::FmpzConstRef exponent) noexcept {
    if (flint::fmpz_is_zero(exponent)) {
        return true;
    }

    FactoredElement power(*accumulator.parent());
    FactoredElement product(*accumulator.parent());
    if (flint::fmpz_fits_si(exponent)) {
        return power.pow_si(base, flint::fmpz_get_si(exponent)) &&
               product.multiply(accumulator, power) &&
               accumulator.set(product);
    }

    Element base_value(*accumulator.parent());
    Element powered_value(*accumulator.parent());
    return base.evaluate(base_value) &&
           powered_value.pow_fmpz(base_value, exponent) &&
           power.set_element(powered_value) &&
           product.multiply(accumulator, power) &&
           accumulator.set(product);
}

bool multiply_element_power_fmpz(FactoredElement& accumulator,
                                 const Element& base,
                                 flint::FmpzConstRef exponent) noexcept {
    if (flint::fmpz_is_zero(exponent)) {
        return true;
    }
    if (flint::fmpz_fits_si(exponent)) {
        return accumulator.push(base, flint::fmpz_get_si(exponent));
    }

    FactoredElement factored_base(*accumulator.parent());
    return factored_base.set_element(base) &&
           compact_multiply_power_fmpz(accumulator, factored_base, exponent);
}

slong fmpz_mat_trimmed_nonzero_rows(flint::FmpzMatConstRef matrix) noexcept {
    slong rows = flint::fmpz_mat_nrows(matrix);
    while (rows > 0 &&
           ::fmpz_mat_is_zero_row(matrix.raw(), rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (::fmpz_mat_is_zero_row(matrix.raw(), i) != 0) {
            return -1;
        }
    }
    return rows;
}

bool fmpz_mat_copy_row(flint::FmpzMatRef out,
                       slong out_row,
                       flint::FmpzMatConstRef in,
                       slong in_row) noexcept {
    if (out_row < 0 || out_row >= flint::fmpz_mat_nrows(out) ||
        in_row < 0 || in_row >= flint::fmpz_mat_nrows(in) ||
        flint::fmpz_mat_ncols(out) != flint::fmpz_mat_ncols(in)) {
        return false;
    }
    for (slong j = 0; j < flint::fmpz_mat_ncols(in); ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, out_row, j),
                        flint::fmpz_mat_entry(in, in_row, j));
    }
    return true;
}

struct ClassRelationHnfTransform {
    flint::FmpzMat stored{0, 0};
    flint::FmpzMat hnf{0, 0};
    flint::FmpzMat transform{0, 0};
    slong rank = 0;
};

bool build_class_relation_hnf_transform(
        ClassRelationHnfTransform& out,
        const ClassGroupContext& context) noexcept {
    const Order* order = context.parent();
    if (!context.has_factor_base() || order == nullptr ||
        order->parent() == nullptr ||
        context.relation_count() < context.relation_rank() ||
        context.generator_count() <= 0 ||
        context.relation_rank() > context.generator_count()) {
        return false;
    }

    const slong relation_count = context.relation_count();
    const slong generator_count = context.generator_count();
    flint::FmpzMat stored(relation_count, generator_count);
    if (!context.relations(flint::FmpzMatRef(stored))) {
        return false;
    }

    ClassRelationHnfTransform candidate;
    candidate.stored = std::move(stored);
    candidate.hnf = flint::FmpzMat(relation_count, generator_count);
    candidate.transform = flint::FmpzMat(relation_count, relation_count);
    ::fmpz_mat_hnf_transform(candidate.hnf.raw(), candidate.transform.raw(),
                             candidate.stored.raw());

    const slong rank =
            fmpz_mat_trimmed_nonzero_rows(
                    flint::FmpzMatConstRef(candidate.hnf));
    if (rank < 0 || rank != context.relation_rank()) {
        return false;
    }
    candidate.rank = rank;
    out = std::move(candidate);
    return true;
}

bool push_context_relation_witnesses(
        FactoredElement& out,
        const ClassGroupContext& context,
        flint::FmpzMatConstRef coefficients,
        slong row) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*field) || row < 0 ||
        row >= flint::fmpz_mat_nrows(coefficients) ||
        flint::fmpz_mat_ncols(coefficients) != context.relation_count()) {
        return false;
    }

    FactoredElement candidate(*field);
    Element generator(*field);
    if (!candidate.is_defined() || !generator.is_defined() ||
        !candidate.one()) {
        return false;
    }
    for (slong col = 0; col < context.relation_count(); ++col) {
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(coefficients, row, col);
        if (flint::fmpz_is_zero(exponent)) {
            continue;
        }
        if (!context.relation_generator(generator, col) ||
            !multiply_element_power_fmpz(candidate, generator, exponent)) {
            return false;
        }
    }
    candidate.normalize();
    out.swap(candidate);
    return true;
}

}  // namespace

bool hnf_unit_witness_coefficients(
        flint::FmpzMat& out,
        const ClassGroupContext& context) noexcept {
    out = flint::FmpzMat(0, 0);
    ClassRelationHnfTransform hnf_data;
    if (!build_class_relation_hnf_transform(hnf_data, context)) {
        return false;
    }

    const slong relation_count = context.relation_count();
    const slong rank = hnf_data.rank;
    flint::FmpzMat coefficients(relation_count - rank, relation_count);
    for (slong row = rank; row < relation_count; ++row) {
        for (slong col = 0; col < relation_count; ++col) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatRef(coefficients), row - rank, col),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(hnf_data.transform), row,
                            col));
        }
    }
    out = std::move(coefficients);
    return true;
}

bool hnf_unit_witnesses(std::vector<FactoredElement>& out,
                             const ClassGroupContext& context) noexcept {
    out.clear();
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr) {
        return false;
    }

    flint::FmpzMat coefficients(0, 0);
    if (!hnf_unit_witness_coefficients(coefficients, context)) {
        return false;
    }

    const slong rows = flint::fmpz_mat_nrows(coefficients);
    out.reserve(static_cast<std::size_t>(rows));
    for (slong row = 0; row < rows; ++row) {
        FactoredElement candidate(*field);
        if (!candidate.is_defined() ||
            !push_context_relation_witnesses(
                    candidate, context,
                    flint::FmpzMatConstRef(coefficients), row)) {
            out.clear();
            return false;
        }
        out.push_back(std::move(candidate));
    }
    return true;
}

bool class_relation_witnessed_hnf_basis(
        WitnessedClassRelationHnfBasis& out,
        const ClassGroupContext& context) noexcept {
    const Order* order = context.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    ClassRelationHnfTransform hnf_data;
    if (field == nullptr ||
        !build_class_relation_hnf_transform(hnf_data, context) ||
        hnf_data.rank != context.generator_count()) {
        return false;
    }

    const slong rank = hnf_data.rank;
    const slong relation_count = context.relation_count();
    const slong generator_count = context.generator_count();
    WitnessedClassRelationHnfBasis candidate;
    candidate.rows = flint::FmpzMat(rank, generator_count);
    candidate.relation_coefficients =
            flint::FmpzMat(rank, relation_count);
    for (slong row = 0; row < rank; ++row) {
        if (!fmpz_mat_copy_row(flint::FmpzMatRef(candidate.rows), row,
                               flint::FmpzMatConstRef(hnf_data.hnf), row) ||
            !fmpz_mat_copy_row(
                    flint::FmpzMatRef(candidate.relation_coefficients), row,
                    flint::FmpzMatConstRef(hnf_data.transform), row)) {
            return false;
        }
    }

    flint::FmpzMat product(rank, generator_count);
    flint::fmpz_mat_mul(
            flint::FmpzMatRef(product),
            flint::FmpzMatConstRef(candidate.relation_coefficients),
            flint::FmpzMatConstRef(hnf_data.stored));
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(product),
                               flint::FmpzMatConstRef(candidate.rows))) {
        return false;
    }

    Ideal one(*order);
    flint::FmpzMat row(1, generator_count);
    if (!one.is_defined() || !one.one()) {
        return false;
    }
    candidate.witnesses.reserve(static_cast<std::size_t>(rank));
    for (slong i = 0; i < rank; ++i) {
        FactoredElement witness(*field);
        if (!witness.is_defined() ||
            !push_context_relation_witnesses(
                    witness, context,
                    flint::FmpzMatConstRef(candidate.relation_coefficients),
                    i) ||
            !fmpz_mat_copy_row(flint::FmpzMatRef(row), 0,
                               flint::FmpzMatConstRef(candidate.rows), i) ||
            !verify_class_group_ideal_relation_witness(
                    context, one, witness, flint::FmpzMatConstRef(row))) {
            return false;
        }
        candidate.witnesses.push_back(std::move(witness));
    }

    out = std::move(candidate);
    return true;
}

}  // namespace silex::detail
