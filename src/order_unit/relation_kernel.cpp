#include "compute_internal.hpp"
#include "relation_unit_internal.hpp"
#include "order_unit_internal.hpp"

#include "../class_group/class_group_internal.hpp"
#include "../order/order_internal.hpp"
#include "../residue_field/residue_field_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/ulong_extras.h>

#include <silex/archimedean.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpz_factor.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/residue_field.hpp>
#include <silex/unit.hpp>

namespace silex {
namespace detail {

constexpr slong kUnitCompactArbPrecision = 128;
constexpr slong kUnitCompactShortPrecision = 128;
constexpr slong kDirectFactoredExponentBound = 10;

// Private replay storage matching reference's fixed S-unit basis with integer CU
// coordinates and reference's FacElem bases with ZZRingElem exponents.
struct CompactRelationUnitGroupStorage {
    explicit CompactRelationUnitGroupStorage(const Order& source_order) noexcept
        : order(source_order) {
    }

    Order order;
    std::vector<Element> relation_generators;
    flint::FmpzMat exponents{0, 0};
    flint::Arb regulator;
    slong work_precision = 0;
    slong exponents_exceeding_slong = 0;
    slong maximum_exponent_bits = 0;
    bool exact_kernel_verified = false;
};

CompactRelationUnitGroup::CompactRelationUnitGroup() noexcept = default;

CompactRelationUnitGroup::~CompactRelationUnitGroup() noexcept = default;

CompactRelationUnitGroup::CompactRelationUnitGroup(
        CompactRelationUnitGroup&& other) noexcept = default;

CompactRelationUnitGroup& CompactRelationUnitGroup::operator=(
        CompactRelationUnitGroup&& other) noexcept = default;

void CompactRelationUnitGroup::swap(
        CompactRelationUnitGroup& other) noexcept {
    storage_.swap(other.storage_);
}

void CompactRelationUnitGroup::clear() noexcept {
    storage_.reset();
}

bool CompactRelationUnitGroup::is_set() const noexcept {
    return storage_ != nullptr;
}

slong CompactRelationUnitGroup::free_rank() const noexcept {
    return storage_ == nullptr
            ? -1
            : flint::fmpz_mat_nrows(storage_->exponents);
}

slong CompactRelationUnitGroup::relation_generator_count() const noexcept {
    return storage_ == nullptr
            ? 0
            : static_cast<slong>(storage_->relation_generators.size());
}

slong CompactRelationUnitGroup::exponents_exceeding_slong() const noexcept {
    return storage_ == nullptr ? 0 : storage_->exponents_exceeding_slong;
}

slong CompactRelationUnitGroup::maximum_exponent_bits() const noexcept {
    return storage_ == nullptr ? 0 : storage_->maximum_exponent_bits;
}

bool CompactRelationUnitGroup::exact_kernel_verified() const noexcept {
    return storage_ != nullptr && storage_->exact_kernel_verified;
}

bool CompactRelationUnitGroup::regulator(flint::ArbRef out) const noexcept {
    if (storage_ == nullptr) {
        return false;
    }
    flint::arb_set(out.raw(), storage_->regulator);
    return true;
}

bool CompactRelationUnitGroup::class_regulator_index_bound(
        flint::FmpzRef out,
        const ClassGroupContext& class_group,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision) const noexcept {
    const Order* class_order_parent = class_group.parent();
    if (storage_ == nullptr || class_order_parent == nullptr ||
        !storage_->order.has_same_data(*class_order_parent)) {
        return false;
    }

    flint::Fmpz class_order;
    flint::Arb candidate_product;
    if (!class_group.order(flint::FmpzRef(class_order))) {
        return false;
    }
    flint::arb_mul_fmpz(candidate_product, storage_->regulator,
                        flint::FmpzConstRef(class_order), precision);
    return class_regulator_index_bound_from_candidate_product(
            out, flint::ArbConstRef(candidate_product),
            analytic_class_regulator_product, precision,
            class_group.diagnostics());
}

namespace {

bool rank_one_proof_prime_residue_data(ResidueField& residue_field,
                                       flint::Fmpz& cofactor,
                                       const PrimeIdeal& prime,
                                       flint::FmpzConstRef ell) noexcept {
    if (prime.residue_degree() != 1 || !flint::fmpz_is_prime(ell) ||
        !residue_field.set_prime(prime)) {
        return false;
    }

    flint::Fmpz cardinality;
    flint::Fmpz group_order;
    if (!residue_field.cardinality(flint::FmpzRef(cardinality))) {
        return false;
    }
    fmpz_sub_ui(group_order.raw(), cardinality.raw(), 1);
    if (fmpz_divisible(group_order.raw(), ell.raw()) == 0) {
        return false;
    }
    fmpz_divexact(cofactor.raw(), group_order.raw(), ell.raw());
    return true;
}

bool rank_one_proof_prime_image_nonzero(
        bool& nonzero,
        const ResidueField& residue_field,
        const ResidueFieldElement& image,
        flint::FmpzConstRef cofactor) noexcept {
    ResidueFieldElement quotient_image(residue_field);
    ResidueFieldElement one(residue_field);
    if (!one.one() || !quotient_image.pow_fmpz(image, cofactor)) {
        return false;
    }

    nonzero = !quotient_image.equal(one);
    return true;
}

bool rank_one_factored_image_nonzero_degree_one(
        bool& nonzero,
        const FactoredElement& generator,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    nonzero = false;
    flint::Fmpz p;
    flint::Fmpz root;
    if (!degree_one_prime_root_mod_p(root, p, prime)) {
        return false;
    }

    return rank_one_factored_image_nonzero_at_degree_one_root(
            nonzero, generator, flint::FmpzConstRef(p),
            flint::FmpzConstRef(root), ell);
}

}  // namespace

bool rank_one_factored_image_nonzero_at_degree_one_root(
        bool& nonzero,
        const FactoredElement& generator,
        flint::FmpzConstRef p,
        flint::FmpzConstRef root,
        flint::FmpzConstRef ell) noexcept {
    nonzero = false;
    if (!generator.is_defined() || !flint::fmpz_is_prime(p) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::Fmpz group_order;
    flint::Fmpz cofactor;
    ::fmpz_sub_ui(group_order.raw(), p.raw(), 1);
    if (::fmpz_divisible(group_order.raw(), ell.raw()) == 0) {
        return false;
    }
    ::fmpz_divexact(cofactor.raw(), group_order.raw(), ell.raw());

    flint::Fmpz product;
    if (!factored_value_at_degree_one_root(product, generator, p, root)) {
        return false;
    }

    ::fmpz_powm(product.raw(), product.raw(), cofactor.raw(), p.raw());
    nonzero = !flint::fmpz_is_one(flint::FmpzConstRef(product));
    return true;
}

bool rank_one_factored_image_nonzero_at_degree_one_root_nmod(
        bool& nonzero,
        const FactoredElement& generator,
        ulong p,
        ulong root,
        flint::FmpzConstRef ell) noexcept {
    nonzero = false;
    if (!generator.is_defined() || p < 2 ||
        !flint::fmpz_abs_fits_ui(ell)) {
        return false;
    }

    const ulong ell_ui = flint::fmpz_get_ui(ell);
    if (ell_ui < 2 || n_is_prime(ell_ui) == 0 ||
        (p - 1) % ell_ui != 0) {
        return false;
    }
    const ulong cofactor = (p - 1) / ell_ui;

    ulong product = 0;
    if (!factored_value_at_degree_one_root_nmod(product, generator, p,
                                                root % p)) {
        return false;
    }

    product = n_powmod2_ui_preinv(product, cofactor, p, n_preinvert_limb(p));
    nonzero = product != 1;
    return true;
}

bool validate_relation_kernel_inputs(const Order& order,
                                     const ClassGroupContext& class_group,
                                     EmbeddingContext& embeddings,
                                     slong precision) noexcept {
    return detail::order_has_parented_basis(order) &&
           same_order_parent(class_group.parent(), &order) &&
           class_group.has_presentation() &&
           embeddings.parent() != nullptr &&
           embeddings.parent()->has_same_data(*order.parent()) &&
           precision > 0;
}

bool relation_kernel_generator(FactoredElement& out,
                               const Order& order,
                               const ClassGroupContext& class_group,
                               slong index) noexcept {
    const NumberField* field = order.parent();
    return field != nullptr && out.parent() != nullptr &&
           out.parent()->has_same_data(*field) &&
           class_group.relation_kernel_unit(out, index);
}

bool copy_free_generators(std::vector<FactoredElement>& out,
                          const OrderUnitGroup& group) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || rank < 0) {
        return false;
    }

    out.clear();
    out.reserve(static_cast<std::size_t>(rank));
    for (slong i = 0; i < rank; ++i) {
        out.emplace_back(*field);
        if (!group.free_generator(out.back(), i)) {
            return false;
        }
    }
    return true;
}

bool compact_multiply_power(FactoredElement& accumulator,
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

bool compact_relation_product(FactoredElement& out,
                              const OrderUnitGroup& group,
                              const FactoredElement& y,
                              flint::FmpzMatConstRef relation,
                              slong row) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || y.parent() == nullptr ||
        !y.parent()->has_same_data(*field) ||
        row < 0 || row >= flint::fmpz_mat_nrows(relation) ||
        flint::fmpz_mat_ncols(relation) != rank + 1 ||
        out.parent() == nullptr || !out.parent()->has_same_data(*field)) {
        return false;
    }

    FactoredElement product(*field);
    if (!product.one()) {
        return false;
    }

    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !compact_multiply_power(
                    product, generator,
                    flint::fmpz_mat_entry(relation, row, i))) {
            return false;
        }
    }
    if (!compact_multiply_power(
                product, y, flint::fmpz_mat_entry(relation, row, rank))) {
        return false;
    }

    return out.set(product);
}

bool relation_row_normalized(flint::FmpzMat& relation,
                             flint::Fmpz& exponent,
                             const OrderUnitGroup& group,
                             flint::FmpzMatConstRef rel,
                             slong row) noexcept;

bool relation_power_rhs(FactoredElement& rhs,
                        const OrderUnitGroup& group,
                        flint::FmpzMatConstRef relation) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || rhs.parent() == nullptr ||
        !rhs.parent()->has_same_data(*field) ||
        flint::fmpz_mat_nrows(relation) != 1 ||
        flint::fmpz_mat_ncols(relation) != rank) {
        return false;
    }

    if (!rhs.one()) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !compact_multiply_power(rhs, generator,
                                    flint::fmpz_mat_entry(relation, 0, i))) {
            return false;
        }
    }
    return true;
}

bool factored_structurally_equal(const FactoredElement& left,
                                 const FactoredElement& right) noexcept {
    if (!left.is_defined() || !right.is_defined() ||
        left.parent() == nullptr || right.parent() == nullptr ||
        !left.parent()->has_same_data(*right.parent())) {
        return false;
    }

    const auto left_factors = left.factors();
    const auto right_factors = right.factors();
    if (left_factors.size() != right_factors.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left_factors.size(); ++i) {
        if (left_factors[i].exponent != right_factors[i].exponent ||
            !left_factors[i].factor.equal(right_factors[i].factor)) {
            return false;
        }
    }
    return true;
}

bool relation_power_root_matches_y(bool& matches,
                                   flint::Fmpz& quotient_exp,
                                   const OrderUnitGroup& group,
                                   const FactoredElement& root,
                                   const FactoredElement& y) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_power_quotient");
    matches = false;
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!group.is_set() || field == nullptr ||
        root.parent() == nullptr || !root.parent()->has_same_data(*field) ||
        y.parent() == nullptr || !y.parent()->has_same_data(*field)) {
        return false;
    }

    flint::Fmpz torsion_order;
    if (!group.torsion_order(flint::FmpzRef(torsion_order)) ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(torsion_order)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(torsion_order)) <= 0) {
        return false;
    }

    FactoredElement quotient(*field);
    Element quotient_value(*field);
    OrderElement torsion(*order);
    Element torsion_value(*field);
    Element torsion_power(*field);
    if (!group.torsion_generator(torsion) ||
        !torsion.get_element(torsion_value)) {
        return false;
    }

    const slong order_value =
            flint::fmpz_get_si(flint::FmpzConstRef(torsion_order));
    for (slong i = 0; i < order_value; ++i) {
        FactoredElement candidate(*field);
        if (!candidate.set(y)) {
            return false;
        }
        if (i != 0 && !candidate.push(torsion_value, i)) {
            return false;
        }
        candidate.normalize();
        if (factored_structurally_equal(root, candidate)) {
            flint::fmpz_set_si(flint::FmpzRef(quotient_exp), i);
            matches = true;
            return true;
        }
    }

    if (!quotient.divide(root, y) ||
        !quotient.evaluate(quotient_value)) {
        return false;
    }

    for (slong i = 0; i < order_value; ++i) {
        if (!compute_power(torsion_power, torsion_value, i)) {
            return false;
        }
        if (quotient_value.equal(torsion_power)) {
            flint::fmpz_set_si(flint::FmpzRef(quotient_exp), i);
            matches = true;
            return true;
        }
    }

    return true;
}

bool relation_power_root(bool& is_relation,
                         FactoredElement& root,
                         flint::Fmpz& candidate_exp,
                         const OrderUnitGroup& group,
                         const FactoredElement& y,
                         flint::FmpzMatConstRef rel,
                         slong row,
                         bool require_y_root) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_power_root");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr ||
        root.parent() == nullptr || !root.parent()->has_same_data(*field) ||
        y.parent() == nullptr || !y.parent()->has_same_data(*field) ||
        row < 0 || row >= flint::fmpz_mat_nrows(rel) ||
        flint::fmpz_mat_ncols(rel) != rank + 1) {
        return false;
    }

    is_relation = false;
    flint::Fmpz torsion_order;
    if (!group.torsion_order(flint::FmpzRef(torsion_order)) ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(torsion_order)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(torsion_order)) <= 0) {
        return false;
    }

    flint::FmpzMat relation(1, rank);
    flint::Fmpz exponent;
    if (!relation_row_normalized(relation, exponent, group, rel, row)) {
        return false;
    }
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(exponent)) ||
        flint::fmpz_cmp_ui(flint::FmpzConstRef(exponent), 2) < 0) {
        return true;
    }
    const slong exponent_si = flint::fmpz_get_si(flint::FmpzConstRef(exponent));

    FactoredElement rhs(*field);
    if (!relation_power_rhs(rhs, group, flint::FmpzMatConstRef(relation))) {
        return false;
    }

    OrderElement torsion(*order);
    Element zeta(*field);
    if (!group.torsion_generator(torsion) ||
        !torsion.get_element(zeta)) {
        return false;
    }

    const slong order_value =
            flint::fmpz_get_si(flint::FmpzConstRef(torsion_order));
    detail::CompactFieldModulusCache field_modulus_cache;
    for (slong i = 0; i < order_value; ++i) {
        FactoredElement candidate(*field);
        if (!candidate.set(rhs)) {
            return false;
        }
        if (i != 0 && !candidate.push(zeta, i)) {
            return false;
        }

        bool is_power = false;
        {
            SILEX_PROFILE_SCOPE(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.dependent_relation_power");
            // reference `_ev` evaluates residual exponents with absolute value
            // below 10 directly and recursively halves every larger one.
            bool use_full_presentation = false;
            for (const auto& entry : candidate.factors()) {
                slong remainder = entry.exponent % exponent_si;
                if (remainder < 0) {
                    remainder += exponent_si;
                }
                if (remainder >= kDirectFactoredExponentBound) {
                    use_full_presentation = true;
                    break;
                }
            }

            const bool compact_available =
                    use_full_presentation &&
                    detail::compact_unit_presentation_power_root(
                            is_power, root, *order, candidate, exponent_si,
                            kUnitCompactArbPrecision,
                            kUnitCompactShortPrecision,
                            group.diagnostics(), &field_modulus_cache);
            if ((!use_full_presentation || !compact_available) &&
                !candidate.is_power_si(
                        is_power, root, exponent_si,
                        FactoredRootStrategy::compact,
                        group.diagnostics())) {
                return false;
            }
        }
        if (is_power) {
            if (!require_y_root) {
                flint::fmpz_set_si(flint::FmpzRef(candidate_exp), i);
                is_relation = true;
                return true;
            }
            bool matches_y = false;
            flint::Fmpz quotient_exp;
            if (!relation_power_root_matches_y(
                        matches_y, quotient_exp, group, root, y)) {
                return false;
            }
            if (matches_y) {
                if (!root.set(y)) {
                    return false;
                }
                flint::fmpz_set_si(flint::FmpzRef(candidate_exp), i);
                is_relation = true;
                return true;
            }
        }

        {
            SILEX_PROFILE_SCOPE(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.dependent_relation_power");
            if (!candidate.is_power_si(is_power, root, exponent_si,
                                       FactoredRootStrategy::reduced,
                                       group.diagnostics())) {
                return false;
            }
        }
        if (is_power) {
            if (!require_y_root) {
                flint::fmpz_set_si(flint::FmpzRef(candidate_exp), i);
                is_relation = true;
                return true;
            }
            bool matches_y = false;
            flint::Fmpz quotient_exp;
            if (!relation_power_root_matches_y(
                        matches_y, quotient_exp, group, root, y)) {
                return false;
            }
            if (matches_y) {
                if (!root.set(y)) {
                    return false;
                }
                flint::fmpz_set_si(flint::FmpzRef(candidate_exp), i);
                is_relation = true;
                return true;
            }
        }
    }

    return true;
}

bool check_relation_mod_torsion(bool& is_relation,
                                flint::FmpzRef torsion_exp,
                                const OrderUnitGroup& group,
                                const FactoredElement& y,
                                flint::FmpzMatConstRef relation,
                                slong row,
                                EmbeddingContext* embeddings,
                                bool require_torsion_exponent) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_check_torsion");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!group.is_set() || field == nullptr || y.parent() == nullptr ||
        !y.parent()->has_same_data(*field)) {
        return false;
    }

    flint::Fmpz torsion_order;
    OrderElement torsion(*order);
    Element torsion_value(*field);
    Element product_value(*field);
    FactoredElement product(*field);
    if (!group.torsion_order(flint::FmpzRef(torsion_order)) ||
        !flint::fmpz_fits_si(flint::FmpzConstRef(torsion_order)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(torsion_order)) <= 0 ||
        !group.torsion_generator(torsion) ||
        !torsion.get_element(torsion_value) ||
        !compact_relation_product(product, group, y, relation, row)) {
        return false;
    }

    // reference Unit/Relation.jl checks the compact product with
    // is_torsion_unit(FacElem), whose certified log test rejects a
    // non-torsion unit without expanding it.  Silex C likewise screens logs
    // before exact expansion in group_check_relation_mod_torsion.c, while
    // retaining the exact path when it needs the torsion exponent.  Failure
    // or an inconclusive log screen deliberately falls through to that path.
    if (embeddings != nullptr) {
        SILEX_PROFILE_SCOPE(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.dependent_relation_check_torsion.log_screen");
        RelationUnitExtractionState torsion_state;
        RelationTorsionStatus torsion_status = RelationTorsionStatus::inconclusive;
        if (unit_candidate_torsion_status(
                    torsion_status, torsion_state, product, *embeddings)) {
            if (torsion_status == RelationTorsionStatus::non_torsion) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_check_torsion.log_reject");
                is_relation = false;
                return true;
            }
            if (torsion_status == RelationTorsionStatus::torsion &&
                !require_torsion_exponent) {
                SILEX_PROFILE_EVENT(
                        group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_check_torsion.log_accept");
                is_relation = true;
                flint::fmpz_zero(torsion_exp);
                return true;
            }
        }
    }

    SILEX_PROFILE_SCOPE(
            group.diagnostics(), DiagnosticsModule::unit_group,
            "unit_group.dependent_relation_check_torsion.exact_fallback");
    if (!product.evaluate(product_value)) {
        return false;
    }

    const slong order_value =
            flint::fmpz_get_si(flint::FmpzConstRef(torsion_order));
    Element torsion_power(*field);
    for (slong k = 0; k < order_value; ++k) {
        if (!compute_power(torsion_power, torsion_value, k)) {
            return false;
        }
        if (product_value.equal(torsion_power)) {
            is_relation = true;
            flint::fmpz_set_si(torsion_exp, k);
            return true;
        }
    }

    is_relation = false;
    return true;
}

bool relation_row_normalized(flint::FmpzMat& relation,
                             flint::Fmpz& exponent,
                             const OrderUnitGroup& group,
                             flint::FmpzMatConstRef rel,
                             slong row) noexcept {
    const slong rank = group.free_rank();
    if (rank < 0 || flint::fmpz_mat_nrows(relation) != 1 ||
        flint::fmpz_mat_ncols(relation) != rank ||
        row < 0 || row >= flint::fmpz_mat_nrows(rel) ||
        flint::fmpz_mat_ncols(rel) != rank + 1) {
        return false;
    }

    flint::FmpzConstRef last = flint::fmpz_mat_entry(rel, row, rank);
    if (flint::fmpz_sgn(last) > 0) {
        flint::fmpz_set(flint::FmpzRef(exponent), last);
        for (slong i = 0; i < rank; ++i) {
            flint::fmpz_neg(flint::fmpz_mat_entry(relation, 0, i),
                            flint::fmpz_mat_entry(rel, row, i));
        }
    } else {
        flint::fmpz_neg(flint::FmpzRef(exponent), last);
        for (slong i = 0; i < rank; ++i) {
            flint::fmpz_set(flint::fmpz_mat_entry(relation, 0, i),
                            flint::fmpz_mat_entry(rel, row, i));
        }
    }
    return true;
}

bool relation_basis(std::vector<FactoredElement>& out,
                    const OrderUnitGroup& group,
                    const FactoredElement& root,
                    flint::FmpzMatConstRef relation,
                    flint::FmpzConstRef exponent) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || root.parent() == nullptr ||
        !root.parent()->has_same_data(*field) ||
        rank <= 0 || flint::fmpz_mat_nrows(relation) != 1 ||
        flint::fmpz_mat_ncols(relation) != rank) {
        return false;
    }

    std::vector<FactoredElement> old_generators;
    if (!copy_free_generators(old_generators, group)) {
        return false;
    }

    // Source trace: reference `UnitGrpCtx.jl:_find_new_basis` builds the
    // relation column, computes `h, u = hnf_with_transform(m)`, uses
    // `inv(u)[:, 2:r]`, and LLL-reduces the transposed basis before applying
    // it to `vcat(U.units, y)`.
    const slong extended_rank = rank + 1;
    flint::FmpzMat relation_column(extended_rank, 1);
    for (slong i = 0; i < rank; ++i) {
        flint::fmpz_neg(
                flint::fmpz_mat_entry(relation_column, i, 0),
                flint::FmpzConstRef(flint::fmpz_mat_entry(relation, 0, i)));
    }
    flint::fmpz_set(flint::fmpz_mat_entry(relation_column, rank, 0),
                    exponent);

    flint::FmpzMat hnf(extended_rank, 1);
    flint::FmpzMat transform(extended_rank, extended_rank);
    ::fmpz_mat_hnf_transform(hnf.raw(), transform.raw(),
                             relation_column.raw());
    if (!flint::fmpz_is_one(
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(hnf), 0, 0))) {
        return false;
    }

    flint::FmpzMat inverse_num(extended_rank, extended_rank);
    flint::Fmpz inverse_den;
    if (::fmpz_mat_inv(inverse_num.raw(), inverse_den.raw(),
                       transform.raw()) == 0 ||
        !flint::fmpz_is_pm1(flint::FmpzConstRef(inverse_den))) {
        return false;
    }
    if (flint::fmpz_sgn(flint::FmpzConstRef(inverse_den)) < 0) {
        ::fmpz_mat_neg(inverse_num.raw(), inverse_num.raw());
    }

    flint::FmpzMat basis_columns(extended_rank, rank);
    for (slong i = 0; i < extended_rank; ++i) {
        for (slong j = 0; j < rank; ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(basis_columns, i, j),
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(inverse_num), i, j + 1));
        }
    }

    flint::FmpzMat lll_rows(rank, extended_rank);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(lll_rows),
                              flint::FmpzMatConstRef(basis_columns));
    if (rank > 1) {
        flint::FmpzMat lll_transform(rank, rank);
        flint::FmpzLll lll;
        ::fmpz_lll(lll_rows.raw(), lll_transform.raw(), lll.raw());
    }
    flint::fmpz_mat_transpose(flint::FmpzMatRef(basis_columns),
                              flint::FmpzMatConstRef(lll_rows));

    out.clear();
    out.reserve(static_cast<std::size_t>(rank));
    for (slong i = 0; i < rank; ++i) {
        out.emplace_back(*field);
        FactoredElement& generator = out.back();
        if (!generator.one() ||
            !compact_multiply_power(
                    generator, root,
                    flint::FmpzConstRef(
                            flint::fmpz_mat_entry(
                                    basis_columns, rank, i).raw()))) {
            return false;
        }
        for (slong j = 0; j < rank; ++j) {
            if (!compact_multiply_power(
                        generator, old_generators[static_cast<std::size_t>(j)],
                        flint::FmpzConstRef(
                                flint::fmpz_mat_entry(
                                        basis_columns, j, i).raw()))) {
                return false;
            }
        }
    }

    return static_cast<slong>(out.size()) == rank;
}

bool kernel_row_divisible(flint::FmpzMatConstRef kernel_rows,
                          slong row,
                          slong len,
                          flint::FmpzConstRef ell) noexcept {
    for (slong i = 0; i < len; ++i) {
        if (!fmpz_divisible(
                    flint::fmpz_mat_entry(kernel_rows, row, i).raw(),
                    ell.raw())) {
            return false;
        }
    }
    return true;
}

bool kernel_row_product(FactoredElement& product,
                        const OrderUnitGroup& group,
                        flint::FmpzMatConstRef kernel_rows,
                        slong row) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr || product.parent() == nullptr ||
        !product.parent()->has_same_data(*field) ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel_rows) ||
        flint::fmpz_mat_ncols(kernel_rows) != rank) {
        return false;
    }

    if (!product.one()) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !compact_multiply_power(
                    product, generator,
                    flint::fmpz_mat_entry(kernel_rows, row, i))) {
            return false;
        }
    }
    return true;
}

bool kernel_row_root(bool& is_power,
                     FactoredElement& root,
                     const OrderUnitGroup& group,
                     flint::FmpzMatConstRef kernel_rows,
                     slong row,
                     flint::FmpzConstRef ell,
                     CompactFieldModulusCache* field_modulus_cache) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!group.is_set() || field == nullptr || root.parent() == nullptr ||
        !root.parent()->has_same_data(*field) ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel_rows) ||
        flint::fmpz_mat_ncols(kernel_rows) != group.free_rank() ||
        !flint::fmpz_is_prime(ell) || !flint::fmpz_fits_si(ell)) {
        return false;
    }

    FactoredElement product(*field);
    FactoredElement candidate(*field);
    {
        SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.kernel_row_root_product");
        if (!kernel_row_product(product, group, kernel_rows, row)) {
            return false;
        }
    }

    const slong exponent = flint::fmpz_get_si(ell);
    SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.kernel_row_root_compact_attempt");
    if (!detail::compact_unit_presentation_power_root(
                is_power, candidate, *order, product, exponent,
                kUnitCompactArbPrecision,
                kUnitCompactShortPrecision, group.diagnostics(),
                field_modulus_cache)) {
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.kernel_row_root_compact_unavailable");
        if (!product.is_power_si(is_power, candidate, exponent,
                                 FactoredRootStrategy::reduced,
                                 group.diagnostics())) {
            return false;
        }
        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                is_power
                        ? "unit_group.kernel_row_root_power"
                        : "unit_group.kernel_row_root_not_power");
    } else {
        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                is_power
                        ? "unit_group.kernel_row_root_compact_power"
                        : "unit_group.kernel_row_root_compact_not_power");
    }
    if (!is_power) {
        return true;
    }

    return root.set(candidate);
}

bool residue_dlog_matrix(flint::FmpzMat& out,
                         const OrderUnitGroup& group,
                         const PrimeIdeal& prime,
                         flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_matrix");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    if (residue_dlog_matrix_direct_degree_one(out, group, prime, ell)) {
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.residue_dlog_matrix.direct_degree_one");
        return true;
    }
    SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_matrix.generic");

    ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::FmpzMat candidate(rank, 1);
    ResidueFieldElement image(residue_field);
    ResidueFieldQuotientLog quotient_log(residue_field);
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !quotient_log.apply(flint::fmpz_mat_entry(candidate, i, 0),
                                image)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

bool residue_dlog_matrix_direct_degree_one(flint::FmpzMat& out,
                                           const OrderUnitGroup& group,
                                           const PrimeIdeal& prime,
                                           flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_direct_degree_one_matrix");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    flint::Fmpz p;
    flint::Fmpz root;
    flint::Fmpz cofactor;
    flint::Fmpz quotient_generator;
    if (!group.is_set() || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !flint::fmpz_is_prime(ell) ||
        !degree_one_prime_root_mod_p(root, p, prime) ||
        !quotient_log_mod_prime_setup(cofactor, quotient_generator,
                                      flint::FmpzConstRef(p), ell)) {
        return false;
    }

    flint::FmpzMat candidate(rank, 1);
    flint::Fmpz value;
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !factored_value_at_degree_one_root(
                    value, generator, flint::FmpzConstRef(p),
                    flint::FmpzConstRef(root)) ||
            !quotient_log_mod_prime_apply(
                    flint::fmpz_mat_entry(candidate, i, 0),
                    flint::FmpzConstRef(value), flint::FmpzConstRef(cofactor),
                    flint::FmpzConstRef(quotient_generator),
                    flint::FmpzConstRef(p), ell)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

bool dlog_proof_rank(slong& out,
                     const OrderUnitGroup& group,
                     flint::FmpzConstRef ell) noexcept {
    if (!group.is_set() || !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::Fmpz torsion_order;
    if (!group.torsion_order(flint::FmpzRef(torsion_order))) {
        return false;
    }
    out = group.free_rank() +
          (fmpz_divisible(torsion_order.raw(), ell.raw()) != 0 ? 1 : 0);
    return out >= 0;
}

bool residue_dlog_proof_matrix(flint::FmpzMat& out,
                               const OrderUnitGroup& group,
                               const PrimeIdeal& prime,
                               flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_proof_matrix");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    slong proof_rank = 0;
    if (!group.is_set() || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !dlog_proof_rank(proof_rank, group, ell)) {
        return false;
    }

    if (proof_rank == 0) {
        flint::FmpzMat candidate(0, 1);
        out = std::move(candidate);
        return true;
    }

    ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::FmpzMat candidate(proof_rank, 1);
    ResidueFieldElement image(residue_field);
    ResidueFieldQuotientLog quotient_log(residue_field);
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !quotient_log.apply(flint::fmpz_mat_entry(candidate, i, 0),
                                image)) {
            return false;
        }
    }

    if (proof_rank > rank) {
        OrderElement torsion(*order);
        if (!group.torsion_generator(torsion) ||
            !image.set_order_element(torsion) ||
            !quotient_log.apply(flint::fmpz_mat_entry(candidate, rank, 0),
                                image)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

bool saturation_proof_prime_column_direct_degree_one(
        flint::FmpzMat& out,
        const OrderUnitGroup& group,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_direct_degree_one_column");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    slong proof_rank = 0;
    flint::Fmpz p;
    flint::Fmpz root;
    if (!group.is_set() || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        !dlog_proof_rank(proof_rank, group, ell) ||
        !degree_one_prime_root_mod_p(root, p, prime)) {
        return false;
    }

    return saturation_proof_prime_column_at_degree_one_root(
            out, group, flint::FmpzConstRef(p), flint::FmpzConstRef(root),
            ell);
}

bool saturation_proof_prime_column_at_degree_one_root(
        flint::FmpzMat& out,
        const OrderUnitGroup& group,
        flint::FmpzConstRef p,
        flint::FmpzConstRef root,
        flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_degree_one_root_column");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    slong proof_rank = 0;
    flint::Fmpz cofactor;
    flint::Fmpz quotient_generator;
    if (!group.is_set() || field == nullptr ||
        !dlog_proof_rank(proof_rank, group, ell) ||
        !quotient_log_mod_prime_setup(cofactor, quotient_generator, p, ell)) {
        return false;
    }

    flint::FmpzMat candidate(proof_rank, 1);
    flint::Fmpz value;
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !factored_value_at_degree_one_root(value, generator, p, root) ||
            !quotient_log_mod_prime_apply(
                    flint::fmpz_mat_entry(candidate, i, 0),
                    flint::FmpzConstRef(value), flint::FmpzConstRef(cofactor),
                    flint::FmpzConstRef(quotient_generator),
                    p, ell)) {
            return false;
        }
    }

    if (proof_rank > rank) {
        OrderElement torsion(*order);
        Element torsion_value(*field);
        if (!group.torsion_generator(torsion) ||
            !torsion.get_element(torsion_value) ||
            !element_value_at_degree_one_root(value, torsion_value, p, root) ||
            !quotient_log_mod_prime_apply(
                    flint::fmpz_mat_entry(candidate, rank, 0),
                    flint::FmpzConstRef(value), flint::FmpzConstRef(cofactor),
                    flint::FmpzConstRef(quotient_generator),
                    p, ell)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

bool dlog_kernel_from_matrix(flint::FmpzMat& out,
                             const flint::FmpzMat& matrix,
                             flint::FmpzConstRef ell) noexcept {
    if (!flint::fmpz_is_prime(ell)) {
        return false;
    }

    const slong rank = flint::fmpz_mat_nrows(matrix);
    const slong len = flint::fmpz_mat_ncols(matrix);
    flint::FmpzModCtx ctx(ell.raw());
    flint::FmpzModMat mod_matrix(rank, len, ctx);
    flint::FmpzModMat transpose(len, rank, ctx);
    flint::FmpzModMat nullspace(rank, rank, ctx);
    flint::Fmpz coefficient;

    fmpz_mod_mat_set_fmpz_mat(mod_matrix.raw(), matrix.raw(), ctx.raw());
    fmpz_mod_mat_transpose(transpose.raw(), mod_matrix.raw(), ctx.raw());
    const slong nullity =
            fmpz_mod_mat_nullspace(nullspace.raw(), transpose.raw(),
                                   ctx.raw());

    flint::FmpzMat candidate(nullity, rank);
    for (slong j = 0; j < nullity; ++j) {
        for (slong i = 0; i < rank; ++i) {
            fmpz_mod_mat_get_entry(coefficient.raw(), nullspace.raw(), i, j,
                                   ctx.raw());
            flint::fmpz_set(flint::fmpz_mat_entry(candidate, j, i),
                            flint::FmpzConstRef(coefficient));
        }
    }

    out = std::move(candidate);
    return true;
}

bool residue_dlog_kernel(flint::FmpzMat& out,
                         const OrderUnitGroup& group,
                         PrimeIdealSpan primes,
                         flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.detail_residue_dlog_kernel");
    const Order* order = group.parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || order == nullptr || primes.empty() ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    flint::FmpzMat matrix(rank, static_cast<slong>(primes.size()));
    for (std::size_t j = 0; j < primes.size(); ++j) {
        if (!same_order_parent(primes[j].parent(), order)) {
            return false;
        }
        flint::FmpzMat column(rank, 1);
        if (!residue_dlog_matrix(column, group, primes[j], ell)) {
            return false;
        }
        for (slong i = 0; i < rank; ++i) {
            flint::fmpz_set(flint::fmpz_mat_entry(
                                    matrix, i, static_cast<slong>(j)),
                            flint::FmpzConstRef(
                                    flint::fmpz_mat_entry(column, i, 0).raw()));
        }
    }

    SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dlog_kernel_nullspace");
    return dlog_kernel_from_matrix(out, matrix, ell);
}

bool residue_dlog_proof_kernel(flint::FmpzMat& out,
                               const OrderUnitGroup& group,
                               PrimeIdealSpan primes,
                               flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.detail_residue_dlog_proof_kernel");
    const Order* order = group.parent();
    slong proof_rank = 0;
    if (!group.is_set() || order == nullptr || primes.empty() ||
        !dlog_proof_rank(proof_rank, group, ell)) {
        return false;
    }

    flint::FmpzMat matrix(proof_rank, static_cast<slong>(primes.size()));
    for (std::size_t j = 0; j < primes.size(); ++j) {
        if (!same_order_parent(primes[j].parent(), order)) {
            return false;
        }
        flint::FmpzMat column(proof_rank, 1);
        if (!residue_dlog_proof_matrix(column, group, primes[j], ell)) {
            return false;
        }
        for (slong i = 0; i < proof_rank; ++i) {
            flint::fmpz_set(flint::fmpz_mat_entry(
                                    matrix, i, static_cast<slong>(j)),
                            flint::FmpzConstRef(
                                    flint::fmpz_mat_entry(column, i, 0).raw()));
        }
    }

    SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_dlog_kernel_nullspace");
    return dlog_kernel_from_matrix(out, matrix, ell);
}

bool saturation_prime_usable(const OrderUnitGroup& group,
                             const PrimeIdeal& prime,
                             flint::FmpzConstRef ell) noexcept {
    flint::FmpzMat logs(0, 0);
    return saturation_prime_column(logs, group, prime, ell);
}

bool saturation_prime_column(flint::FmpzMat& out,
                             const OrderUnitGroup& group,
                             const PrimeIdeal& prime,
                             flint::FmpzConstRef ell) noexcept {
    const Order* order = group.parent();
    if (!group.is_set() || !same_order_parent(prime.parent(), order) ||
        !flint::fmpz_is_prime(ell)) {
        return false;
    }

    ResidueField residue_field(prime);
    flint::Fmpz cardinality;
    flint::Fmpz qminus;
    if (!residue_field.is_defined() ||
        !residue_field.cardinality(flint::FmpzRef(cardinality))) {
        return false;
    }
    fmpz_sub_ui(qminus.raw(), cardinality.raw(), 1);
    if (fmpz_divisible(qminus.raw(), ell.raw()) == 0) {
        return false;
    }

    return residue_dlog_matrix(out, group, prime, ell);
}

bool saturation_proof_prime_usable(const OrderUnitGroup& group,
                                   const PrimeIdeal& prime,
                                   flint::FmpzConstRef ell) noexcept {
    flint::FmpzMat column(0, 0);
    return saturation_proof_prime_column(column, group, prime, ell);
}

bool saturation_proof_prime_column(flint::FmpzMat& out,
                                   const OrderUnitGroup& group,
                                   const PrimeIdeal& prime,
                                   flint::FmpzConstRef ell) noexcept {
    const Order* order = group.parent();
    if (!group.is_set() || !same_order_parent(prime.parent(), order) ||
        prime.residue_degree() != 1 || !flint::fmpz_is_prime(ell)) {
        return false;
    }

    ResidueField residue_field(prime);
    flint::Fmpz cardinality;
    flint::Fmpz qminus;
    flint::FmpzMat logs(0, 0);
    if (!residue_field.is_defined() ||
        !residue_field.cardinality(flint::FmpzRef(cardinality))) {
        return false;
    }
    fmpz_sub_ui(qminus.raw(), cardinality.raw(), 1);
    if (fmpz_divisible(qminus.raw(), ell.raw()) == 0) {
        return false;
    }

    return residue_dlog_proof_matrix(out, group, prime, ell);
}

bool saturation_proof_prime_rank_one_nonzero(
        bool& nonzero,
        const OrderUnitGroup& group,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    nonzero = false;
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    slong proof_rank = 0;
    if (!group.is_set() || field == nullptr ||
        !same_order_parent(prime.parent(), order) ||
        prime.residue_degree() != 1 || !flint::fmpz_is_prime(ell) ||
        !dlog_proof_rank(proof_rank, group, ell) || proof_rank != 1) {
        return false;
    }

    if (group.free_rank() == 1) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, 0)) {
            return false;
        }
        return saturation_proof_prime_known_rank_one_free_generator_nonzero(
                nonzero, group, generator, prime, ell);
    } else if (group.free_rank() == 0) {
        OrderElement torsion(*order);
        if (!group.torsion_generator(torsion)) {
            return false;
        }
        return saturation_proof_prime_known_rank_one_torsion_nonzero(
                nonzero, group, torsion, prime, ell);
    } else {
        return false;
    }
}

bool saturation_proof_prime_known_rank_one_free_generator_nonzero(
        bool& nonzero,
        const OrderUnitGroup& group,
        const FactoredElement& generator,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_rank_one_free_generator_nonzero");
    nonzero = false;
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!group.is_set() || field == nullptr || group.free_rank() != 1 ||
        generator.parent() == nullptr ||
        !generator.parent()->has_same_data(*field) ||
        !same_order_parent(prime.parent(), order)) {
        return false;
    }

    ResidueField residue_field;
    flint::Fmpz cofactor;
    if (!rank_one_proof_prime_residue_data(residue_field, cofactor, prime,
                                           ell)) {
        return false;
    }

    if (rank_one_factored_image_nonzero_degree_one(nonzero, generator, prime,
                                                   ell)) {
        return true;
    }

    ResidueFieldElement image(residue_field);
    if (!image.set_factored_element(generator)) {
        return false;
    }

    return rank_one_proof_prime_image_nonzero(
            nonzero, residue_field, image, flint::FmpzConstRef(cofactor));
}

bool saturation_proof_prime_known_rank_one_torsion_nonzero(
        bool& nonzero,
        const OrderUnitGroup& group,
        const OrderElement& torsion,
        const PrimeIdeal& prime,
        flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.proof_rank_one_torsion_nonzero");
    nonzero = false;
    const Order* order = group.parent();
    if (!group.is_set() || group.free_rank() != 0 ||
        !same_order_parent(torsion.parent(), order) ||
        !same_order_parent(prime.parent(), order)) {
        return false;
    }

    ResidueField residue_field;
    flint::Fmpz cofactor;
    if (!rank_one_proof_prime_residue_data(residue_field, cofactor, prime,
                                           ell)) {
        return false;
    }

    ResidueFieldElement image(residue_field);
    if (!image.set_order_element(torsion)) {
        return false;
    }

    return rank_one_proof_prime_image_nonzero(
            nonzero, residue_field, image, flint::FmpzConstRef(cofactor));
}

bool copy_selected_primes(PrimeIdealList& out,
                          const Order& order,
                          const std::vector<PrimeIdeal>& selected) noexcept {
    PrimeIdealList candidate(order, static_cast<slong>(selected.size()));
    if (!candidate.is_defined()) {
        return false;
    }
    for (slong i = 0; i < static_cast<slong>(selected.size()); ++i) {
        PrimeIdeal* dest = candidate.at(i);
        if (dest == nullptr ||
            !dest->set(selected[static_cast<std::size_t>(i)])) {
            return false;
        }
    }

    out.swap(candidate);
    return true;
}

bool adjoin_verified_dependent_relation(bool& changed,
                                        OrderUnitGroup& out,
                                        const OrderUnitGroup& group,
                                        const FactoredElement& root,
                                        flint::FmpzMatConstRef rel,
                                        slong row,
                                        EmbeddingContext& embeddings,
                                        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_adjoin");
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || !same_order_parent(out.parent(), order) ||
        field == nullptr ||
        root.parent() == nullptr || !root.parent()->has_same_data(*field) ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*field) ||
        precision <= 0 || row < 0 || row >= flint::fmpz_mat_nrows(rel) ||
        flint::fmpz_mat_ncols(rel) != rank + 1) {
        return false;
    }

    if (rank == 0) {
        changed = false;
        return out.set(group);
    }

    flint::FmpzMat relation(1, rank);
    flint::Fmpz exponent;
    if (!relation_row_normalized(relation, exponent, group, rel, row)) {
        return false;
    }

    if (flint::fmpz_is_zero(flint::FmpzConstRef(exponent)) ||
        flint::fmpz_is_pm1(flint::FmpzConstRef(exponent))) {
        changed = false;
        return out.set(group);
    }

    std::vector<FactoredElement> generators;
    if (!relation_basis(generators, group, root,
                        flint::FmpzMatConstRef(relation),
                        flint::FmpzConstRef(exponent))) {
        return false;
    }

    changed = true;
    return order_unit_group_set_units_internal(
            out, *order,
            FactoredElementSpan(generators.data(), generators.size()),
            embeddings, precision, true);
}

bool dependent_relation_rank_zero(bool& recovered,
                                  FactoredElement& root,
                                  flint::FmpzMat& rel,
                                  flint::Fmpz& torsion_exp,
                                  const OrderUnitGroup& group,
                                  const FactoredElement& y) noexcept {
    flint::FmpzMat candidate(1, 1);
    flint::fmpz_one(flint::fmpz_mat_entry(candidate, 0, 0));

    bool is_relation = false;
    if (!check_relation_mod_torsion(
                is_relation, flint::FmpzRef(torsion_exp), group, y,
                flint::FmpzMatConstRef(candidate), 0, nullptr, true)) {
        return false;
    }
    if (is_relation) {
        if (!root.set(y)) {
            return false;
        }
        flint::fmpz_mat_set(flint::FmpzMatRef(rel),
                            flint::FmpzMatConstRef(candidate));
        recovered = true;
    } else {
        recovered = false;
    }
    return true;
}

bool try_denominator(bool& recovered,
                     FactoredElement& root,
                     flint::FmpzMat& rel,
                     flint::Fmpz& torsion_exp,
                     const OrderUnitGroup& group,
                     const FactoredElement& y,
                     EmbeddingContext& embeddings,
                     const flint::ArbMat& coordinates,
                     slong denominator,
                     slong precision,
                     bool require_y_root,
                     bool require_torsion_exponent) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_try_denominator");
    const slong rank = group.free_rank();
    flint::ArbMat scaled(rank, 1);
    flint::FmpzMat candidate(1, rank + 1);
    arb_mat_set(scaled.raw(), coordinates.raw());
    arb_mat_scalar_mul_si(scaled.raw(), scaled.raw(), denominator, precision);

    for (slong i = 0; i < rank; ++i) {
        if (!arb_get_unique_fmpz(
                    flint::fmpz_mat_entry(candidate, 0, i).raw(),
                    arb_mat_entry(scaled.raw(), i, 0))) {
            recovered = false;
            return true;
        }
        flint::fmpz_neg(flint::fmpz_mat_entry(candidate, 0, i),
                        flint::FmpzConstRef(
                                flint::fmpz_mat_entry(candidate, 0, i).raw()));
    }
    flint::fmpz_set_si(flint::fmpz_mat_entry(candidate, 0, rank),
                       denominator);

    bool is_relation = false;
    if (denominator >= 2) {
        flint::Fmpz candidate_exp;
        const bool power_root_ok = relation_power_root(
                is_relation, root, candidate_exp, group, y,
                flint::FmpzMatConstRef(candidate), 0, require_y_root);
        if (power_root_ok && is_relation) {
            flint::fmpz_zero(flint::FmpzRef(torsion_exp));
        } else if (!check_relation_mod_torsion(
                           is_relation, flint::FmpzRef(torsion_exp), group,
                           y, flint::FmpzMatConstRef(candidate), 0,
                           &embeddings, require_torsion_exponent)) {
            return false;
        }
    } else if (!check_relation_mod_torsion(
                       is_relation, flint::FmpzRef(torsion_exp), group, y,
                       flint::FmpzMatConstRef(candidate), 0, &embeddings,
                       require_torsion_exponent)) {
        return false;
    }
    if (is_relation) {
        // reference `_add_dependent_unit!` transforms `vcat(U.units, y)`.  The
        // exact torsion-check fallback may leave an unrelated attempted power
        // root in `root`, so the strict reference path must restore `y` here.
        if ((denominator < 2 || require_y_root) && !root.set(y)) {
            return false;
        }
        flint::fmpz_mat_set(flint::FmpzMatRef(rel),
                            flint::FmpzMatConstRef(candidate));
        recovered = true;
    } else {
        recovered = false;
    }
    return true;
}

bool find_dependent_relation_at_precision(bool& recovered,
                                          FactoredElement& root,
                                          flint::FmpzMat& rel,
                                          flint::Fmpz& torsion_exp,
                                          const OrderUnitGroup& group,
                                          const FactoredElement& y,
                                          EmbeddingContext& embeddings,
                                          slong denominator_bound,
                                          slong min_denominator,
                                          slong precision,
                                          bool require_y_root,
                                          bool require_torsion_exponent) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_at_precision");
    const slong rank = group.free_rank();
    slong places = 0;
    if (rank <= 0 || !compact_places(places, embeddings) ||
        places != rank + 1) {
        return false;
    }

    std::vector<FactoredElement> generators;
    if (!copy_free_generators(generators, group)) {
        return false;
    }

    flint::ArbMat logs(rank, places);
    flint::ArbMat matrix(rank, rank);
    flint::ArbMat rhs(rank, 1);
    flint::ArbMat coordinates(rank, 1);
    flint::ArbVec ylog(places);
    {
        SILEX_PROFILE_SCOPE(group.diagnostics(),
                            DiagnosticsModule::unit_group,
                            "unit_group.dependent_relation_log_solve");
        if (!compact_log_matrix(logs, embeddings,
                                FactoredElementSpan(generators.data(),
                                                    generators.size()),
                                precision) ||
            !y.logarithmic_embedding(flint::ArbVecRef(ylog), embeddings,
                                     LogEmbeddingMode::product, precision)) {
            return false;
        }

        for (slong i = 0; i < rank; ++i) {
            arb_set(arb_mat_entry(rhs.raw(), i, 0), ylog.data() + i);
            for (slong j = 0; j < rank; ++j) {
                arb_set(arb_mat_entry(matrix.raw(), i, j),
                        arb_mat_entry(logs.raw(), j, i));
            }
        }

        if (!arb_mat_solve(coordinates.raw(), matrix.raw(), rhs.raw(),
                           precision)) {
            return false;
        }
    }

    {
        SILEX_PROFILE_SCOPE(group.diagnostics(),
                            DiagnosticsModule::unit_group,
                            "unit_group.dependent_relation_denominator_scan");
        for (slong m = min_denominator; m <= denominator_bound; ++m) {
            if (!try_denominator(recovered, root, rel, torsion_exp, group, y,
                                 embeddings, coordinates, m, precision,
                                 require_y_root, require_torsion_exponent)) {
                return false;
            }
            if (recovered) {
                return true;
            }
        }
    }

    recovered = false;
    return true;
}

bool dependent_relation_bounded_with_inverse(
        bool& recovered,
        FactoredElement& root,
        flint::FmpzMat& rel,
        flint::Fmpz& torsion_exp,
        const OrderUnitGroup& group,
        const FactoredElement& y,
        EmbeddingContext& embeddings,
        const flint::ArbMat& inverse_cutoff,
        flint::FmpzConstRef denominator_bound,
        slong min_denominator,
        slong precision,
        bool require_y_root,
        bool require_torsion_exponent) noexcept {
    const slong rank = group.free_rank();
    if (!group.is_set() || rank <= 0 || y.parent() == nullptr ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*y.parent()) ||
        flint::fmpz_sgn(denominator_bound) <= 0 ||
        !flint::fmpz_fits_si(denominator_bound) ||
        min_denominator <= 0 || precision <= 0 ||
        flint::arb_mat_nrows_value(inverse_cutoff) != rank ||
        flint::arb_mat_ncols_value(inverse_cutoff) != rank ||
        flint::fmpz_mat_nrows(rel) != 1 ||
        flint::fmpz_mat_ncols(rel) != rank + 1) {
        return false;
    }

    const slong bound = flint::fmpz_get_si(denominator_bound);
    if (min_denominator > bound) {
        recovered = false;
        return true;
    }

    slong places = 0;
    if (!compact_places(places, embeddings) || places != rank + 1) {
        return false;
    }

    flint::ArbVec ylog(places);
    if (!y.logarithmic_embedding(flint::ArbVecRef(ylog), embeddings,
                                 LogEmbeddingMode::product, precision)) {
        return false;
    }

    flint::ArbMat yrow(1, rank);
    for (slong i = 0; i < rank; ++i) {
        arb_set(arb_mat_entry(yrow.raw(), 0, i), ylog.data() + i);
    }

    flint::ArbMat coordinate_row(1, rank);
    flint::ArbMat coordinates(rank, 1);
    arb_mat_mul(coordinate_row.raw(), yrow.raw(), inverse_cutoff.raw(),
                precision);
    for (slong i = 0; i < rank; ++i) {
        arb_set(arb_mat_entry(coordinates.raw(), i, 0),
                arb_mat_entry(coordinate_row.raw(), 0, i));
    }

    for (slong m = min_denominator; m <= bound; ++m) {
        if (!try_denominator(recovered, root, rel, torsion_exp, group, y,
                             embeddings, coordinates, m, precision,
                             require_y_root, require_torsion_exponent)) {
            return false;
        }
        if (recovered) {
            return true;
        }
    }

    recovered = false;
    return true;
}

bool dependent_relation_bounded_min_denominator(
        bool& recovered,
        FactoredElement& root,
        flint::FmpzMat& rel,
        flint::Fmpz& torsion_exp,
        const OrderUnitGroup& group,
        const FactoredElement& y,
        EmbeddingContext& embeddings,
        flint::FmpzConstRef denominator_bound,
        slong min_denominator,
        slong start_precision,
        slong max_precision,
        bool require_y_root,
        bool require_torsion_exponent) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_min_denominator");
    const slong rank = group.free_rank();
    if (!group.is_set() || y.parent() == nullptr ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*y.parent()) ||
        flint::fmpz_sgn(denominator_bound) <= 0 ||
        !flint::fmpz_fits_si(denominator_bound) ||
        min_denominator <= 0 || start_precision <= 0 ||
        max_precision < start_precision ||
        flint::fmpz_mat_nrows(rel) != 1 ||
        flint::fmpz_mat_ncols(rel) != rank + 1) {
        return false;
    }

    if (rank == 0) {
        return dependent_relation_rank_zero(recovered, root, rel, torsion_exp,
                                            group, y);
    }

    const slong bound = flint::fmpz_get_si(denominator_bound);
    if (min_denominator > bound) {
        recovered = false;
        return true;
    }

    for (slong precision = start_precision; precision <= max_precision;
         precision *= 2) {
        if (find_dependent_relation_at_precision(
                    recovered, root, rel, torsion_exp, group, y, embeddings,
                    bound, min_denominator, precision, require_y_root,
                    require_torsion_exponent)) {
            return true;
        }
        if (precision > max_precision / 2) {
            break;
        }
    }

    return false;
}

namespace {

void unit_coordinate_fail(OrdinaryUnitCoordinateResult& result,
                          OrdinaryUnitCoordinateStage stage,
                          slong work_precision = 0) noexcept {
    result.success = false;
    result.outcome = OrdinaryUnitCoordinateOutcome::unknown;
    result.stage = stage;
    result.work_precision = work_precision;
}

bool exact_order_unit_membership(bool& is_unit,
                                 const Order& order,
                                 const Element& value) noexcept {
    is_unit = false;
    if (value.equal_si(0)) {
        return true;
    }

    OrderElement integral_value(order);
    if (!integral_value.is_defined()) {
        return false;
    }
    if (!integral_value.set_element(value)) {
        return true;
    }

    Ideal principal(order);
    if (!principal.is_defined() || !principal.set_principal(integral_value)) {
        return false;
    }
    is_unit = principal.is_one();
    return true;
}

bool verify_ordinary_unit_coordinates(
        const OrderUnitGroup& group,
        const FactoredElement& target,
        const OrdinaryUnitCoordinates& coordinates) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (field == nullptr || target.parent() == nullptr ||
        !target.parent()->has_same_data(*field) ||
        flint::fmpz_mat_nrows(coordinates.free_exponents) != 1 ||
        flint::fmpz_mat_ncols(coordinates.free_exponents) != rank) {
        return false;
    }

    flint::Fmpz torsion_order;
    if (!group.torsion_order(flint::FmpzRef(torsion_order)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(torsion_order)) <= 0 ||
        flint::fmpz_sgn(
                flint::FmpzConstRef(coordinates.torsion_exponent)) < 0 ||
        flint::fmpz_cmp(
                flint::FmpzConstRef(coordinates.torsion_exponent),
                flint::FmpzConstRef(torsion_order)) >= 0) {
        return false;
    }

    FactoredElement product(*field);
    OrderElement torsion(*order);
    Element torsion_value(*field);
    FactoredElement torsion_factor(*field);
    if (!product.one() || !group.torsion_generator(torsion) ||
        !torsion.get_element(torsion_value) ||
        !torsion_factor.set_element(torsion_value) ||
        !compact_multiply_power(
                product, torsion_factor,
                flint::FmpzConstRef(coordinates.torsion_exponent))) {
        return false;
    }

    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !compact_multiply_power(
                    product, generator,
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(
                                    coordinates.free_exponents),
                            0, i))) {
            return false;
        }
    }

    Element product_value(*field);
    Element target_value(*field);
    return product.evaluate(product_value) && target.evaluate(target_value) &&
           product_value.equal(target_value);
}

bool ordinary_unit_coordinates_core(
        OrdinaryUnitCoordinateResult& result,
        OrdinaryUnitCoordinates& out,
        const OrderUnitGroup& group,
        const FactoredElement& value,
        const Element& expanded_value,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision) noexcept {
    result = {};
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() ||
        group.certification_status() != CertificationMode::proven ||
        field == nullptr || rank < 0 || value.parent() == nullptr ||
        !value.parent()->has_same_data(*field) ||
        !expanded_value.has_parent(*field) || embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*field) || start_precision <= 0 ||
        max_precision < start_precision) {
        unit_coordinate_fail(result,
                             OrdinaryUnitCoordinateStage::input_validation);
        return false;
    }

    bool is_unit = false;
    if (!exact_order_unit_membership(is_unit, *order, expanded_value)) {
        unit_coordinate_fail(result,
                             OrdinaryUnitCoordinateStage::exact_unit_test);
        return false;
    }
    if (!is_unit) {
        result.success = true;
        result.outcome = OrdinaryUnitCoordinateOutcome::not_unit;
        result.stage = OrdinaryUnitCoordinateStage::none;
        return true;
    }

    flint::Fmpz denominator_bound;
    flint::fmpz_one(flint::FmpzRef(denominator_bound));
    for (slong precision = start_precision; precision <= max_precision;) {
        bool recovered = false;
        FactoredElement root(*field);
        flint::FmpzMat relation(1, rank + 1);
        flint::Fmpz torsion_exponent;
        const bool relation_completed =
                root.is_defined() &&
                dependent_relation_bounded_min_denominator(
                        recovered, root, relation, torsion_exponent, group,
                        value, embeddings,
                        flint::FmpzConstRef(denominator_bound), 1, precision,
                        precision, true, true);
        result.work_precision = precision;
        if (relation_completed && recovered) {
            if (!flint::fmpz_is_one(flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(relation), 0, rank))) {
                unit_coordinate_fail(
                        result,
                        OrdinaryUnitCoordinateStage::exact_verification,
                        precision);
                return false;
            }

            OrdinaryUnitCoordinates candidate;
            candidate.free_exponents = flint::FmpzMat(1, rank);
            flint::fmpz_set(
                    flint::FmpzRef(candidate.torsion_exponent),
                    flint::FmpzConstRef(torsion_exponent));
            for (slong i = 0; i < rank; ++i) {
                flint::fmpz_neg(
                        flint::fmpz_mat_entry(
                                flint::FmpzMatRef(candidate.free_exponents),
                                0, i),
                        flint::fmpz_mat_entry(
                                flint::FmpzMatConstRef(relation), 0, i));
            }
            candidate.work_precision = precision;
            if (!verify_ordinary_unit_coordinates(group, value, candidate)) {
                unit_coordinate_fail(
                        result,
                        OrdinaryUnitCoordinateStage::exact_verification,
                        precision);
                return false;
            }

            candidate.defined = true;
            out = std::move(candidate);
            result.success = true;
            result.outcome = OrdinaryUnitCoordinateOutcome::verified;
            result.stage = OrdinaryUnitCoordinateStage::none;
            result.work_precision = precision;
            return true;
        }

        if (precision > max_precision / 2) {
            break;
        }
        precision *= 2;
    }

    unit_coordinate_fail(result,
                         OrdinaryUnitCoordinateStage::precision_exhausted,
                         result.work_precision);
    return false;
}

}  // namespace

const char* ordinary_unit_coordinate_stage_name(
        OrdinaryUnitCoordinateStage stage) noexcept {
    switch (stage) {
        case OrdinaryUnitCoordinateStage::none:
            return "none";
        case OrdinaryUnitCoordinateStage::input_validation:
            return "input_validation";
        case OrdinaryUnitCoordinateStage::exact_unit_test:
            return "exact_unit_test";
        case OrdinaryUnitCoordinateStage::precision_exhausted:
            return "precision_exhausted";
        case OrdinaryUnitCoordinateStage::exact_verification:
            return "exact_verification";
    }
    return "unknown";
}

bool ordinary_unit_coordinates(
        OrdinaryUnitCoordinateResult& result,
        OrdinaryUnitCoordinates& out,
        const OrderUnitGroup& group,
        const Element& value,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision) noexcept {
    result = {};
    const NumberField* field = value.parent();
    if (field == nullptr) {
        unit_coordinate_fail(result,
                             OrdinaryUnitCoordinateStage::input_validation);
        return false;
    }
    FactoredElement factored(*field);
    if (!factored.is_defined() || !factored.set_element(value)) {
        unit_coordinate_fail(result,
                             OrdinaryUnitCoordinateStage::exact_unit_test);
        return false;
    }
    return ordinary_unit_coordinates_core(
            result, out, group, factored, value, embeddings, start_precision,
            max_precision);
}

bool ordinary_unit_coordinates(
        OrdinaryUnitCoordinateResult& result,
        OrdinaryUnitCoordinates& out,
        const OrderUnitGroup& group,
        const FactoredElement& value,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision) noexcept {
    result = {};
    const NumberField* field = value.parent();
    if (field == nullptr) {
        unit_coordinate_fail(result,
                             OrdinaryUnitCoordinateStage::input_validation);
        return false;
    }
    Element expanded(*field);
    if (!expanded.is_defined() || !value.evaluate(expanded)) {
        unit_coordinate_fail(result,
                             OrdinaryUnitCoordinateStage::exact_unit_test);
        return false;
    }
    return ordinary_unit_coordinates_core(
            result, out, group, value, expanded, embeddings, start_precision,
            max_precision);
}

bool set_best_rank_one_relation_kernel_unit(
        OrderUnitGroup& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.set_best_rank_one_relation_kernel_unit");
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    const slong n = class_group.relation_kernel_unit_count();
    if (n < 1) {
        return false;
    }

    OrderUnitGroup best(order);
    best.set_diagnostics(out.diagnostics());
    FactoredElement candidate(*field);
    flint::Arb regulator;
    flint::Arb best_regulator;
    bool have_best = false;
    flint::Fmpz torsion_order;
    OrderElement torsion_generator(order);
    {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.rank_one_torsion");
        // reference's relation_completion_parameters computes nfrootsof1 once for the BNF pass;
        // reference stores torsion_units on the field/order context.  The
        // rank-one scan varies only the free generator, so reuse torsion.
        if (!torsion_generator.is_defined() ||
            !rank_zero_torsion(flint::FmpzRef(torsion_order),
                               torsion_generator, order)) {
            SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "rank-one relation-kernel torsion unavailable");
            return false;
        }
    }

    for (slong i = 0; i < n; ++i) {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.rank_one_candidate");
        if (!relation_kernel_generator(candidate, order, class_group, i)) {
            SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "rank-one relation-kernel candidate unavailable");
            continue;
        }

        OrderUnitGroup trial(order);
        trial.set_diagnostics(out.diagnostics());
        {
            SILEX_PROFILE_SCOPE(out.diagnostics(),
                                DiagnosticsModule::unit_group,
                                "unit_group.rank_one_set_units");
            if (!order_unit_group_set_units_internal(
                        trial, order, FactoredElementSpan(&candidate, 1),
                        embeddings, precision, true, &torsion_order,
                        &torsion_generator)) {
                SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "rank-one relation-kernel candidate rejected");
                continue;
            }
        }
        {
            SILEX_PROFILE_SCOPE(out.diagnostics(),
                                DiagnosticsModule::unit_group,
                                "unit_group.rank_one_regulator");
            if (!trial.regulator(flint::ArbRef(regulator))) {
                SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "rank-one relation-kernel regulator unavailable");
                continue;
            }
        }

        if (!have_best ||
            regulator_ubound_lt(regulator, best_regulator, precision)) {
            if (!best.set(trial)) {
                return false;
            }
            flint::arb_set(best_regulator, regulator);
            have_best = true;
        }
    }

    if (!have_best) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "rank-one relation-kernel scan found no usable candidate");
        return false;
    }

    out.swap(best);
    return true;
}

bool select_independent_relation_kernel_units(
        std::vector<FactoredElement>& selected,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong precision,
        slong rank) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr || rank <= 0) {
        return false;
    }

    selected.clear();
    selected.reserve(static_cast<std::size_t>(rank));
    bool independent = false;
    const slong n = class_group.relation_kernel_unit_count();

    for (slong i = 0; i < n && static_cast<slong>(selected.size()) < rank;
         ++i) {
        selected.emplace_back(*field);
        FactoredElement& candidate = selected.back();
        if (!relation_kernel_generator(candidate, order, class_group, i)) {
            selected.pop_back();
            continue;
        }

        // reference's FindUnits.jl keeps the exact unit test as @hassert; normal
        // selection works with the factored kernel element and logarithms.
        SILEX_DEBUG_CHECK(
                class_group.diagnostics(), DiagnosticsModule::unit_group,
                DebugLevel::exhaustive,
                "relation-kernel generator is an exact order unit",
                [&]() noexcept {
                    Element expanded(*field);
                    return expanded.is_defined() &&
                           candidate.evaluate(expanded) &&
                           evaluated_is_order_unit(order, expanded);
                }());

        if (!compact_independent(independent, embeddings,
                                 FactoredElementSpan(selected.data(),
                                                     selected.size()),
                                 precision)) {
            selected.pop_back();
            continue;
        }

        if (!independent) {
            selected.pop_back();
        }
    }

    return true;
}

bool relation_coefficients_are_exact_unit_kernel(
        const ClassGroupContext& class_group,
        flint::FmpzMatConstRef coefficients) noexcept {
    const slong rows = flint::fmpz_mat_nrows(coefficients);
    const slong relation_count = class_group.relation_count();
    const slong generator_count = class_group.generator_count();
    if (rows < 0 || relation_count < 0 || generator_count < 0 ||
        flint::fmpz_mat_ncols(coefficients) != relation_count) {
        return false;
    }

    // Each stored relation generator has the corresponding principal-ideal
    // valuation row.  An exact zero product therefore certifies a unit without
    // expanding its potentially enormous powers.
    flint::FmpzMat relations(relation_count, generator_count);
    flint::FmpzMat product(rows, generator_count);
    if (!class_group.relations(flint::FmpzMatRef(relations))) {
        return false;
    }
    flint::fmpz_mat_mul(flint::FmpzMatRef(product), coefficients,
                        flint::FmpzMatConstRef(relations));
    return ::fmpz_mat_is_zero(product.raw()) != 0;
}

bool build_compact_relation_unit_storage(
        std::unique_ptr<CompactRelationUnitGroupStorage>& out,
        const Order& order,
        const ClassGroupContext& class_group,
        const flint::FmpzMat& coefficients,
        flint::ArbConstRef regulator,
        slong work_precision) noexcept {
    const NumberField* field = order.parent();
    const slong relation_count = class_group.relation_count();
    if (field == nullptr || work_precision <= 0 ||
        !flint::arb_is_finite(regulator) || !flint::arb_is_positive(regulator) ||
        flint::fmpz_mat_ncols(coefficients) != relation_count ||
        !relation_coefficients_are_exact_unit_kernel(
                class_group, flint::FmpzMatConstRef(coefficients))) {
        return false;
    }

    std::unique_ptr<CompactRelationUnitGroupStorage> candidate(
            new (std::nothrow) CompactRelationUnitGroupStorage(order));
    if (candidate == nullptr || !candidate->order.is_defined()) {
        return false;
    }

    candidate->relation_generators.reserve(
            static_cast<std::size_t>(relation_count));
    for (slong relation = 0; relation < relation_count; ++relation) {
        candidate->relation_generators.emplace_back(*field);
        if (!candidate->relation_generators.back().is_defined() ||
            !class_group.relation_generator(
                    candidate->relation_generators.back(), relation)) {
            return false;
        }
    }

    candidate->exponents = flint::FmpzMat(
            flint::fmpz_mat_nrows(coefficients),
            flint::fmpz_mat_ncols(coefficients));
    flint::fmpz_mat_set(flint::FmpzMatRef(candidate->exponents),
                        flint::FmpzMatConstRef(coefficients));
    for (slong row = 0; row < flint::fmpz_mat_nrows(coefficients); ++row) {
        for (slong col = 0; col < flint::fmpz_mat_ncols(coefficients); ++col) {
            const flint::FmpzConstRef exponent =
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(coefficients), row, col);
            if (!flint::fmpz_fits_si(exponent)) {
                ++candidate->exponents_exceeding_slong;
            }
            const auto bits = static_cast<slong>(::fmpz_bits(exponent.raw()));
            candidate->maximum_exponent_bits =
                    std::max(candidate->maximum_exponent_bits, bits);
        }
    }

    flint::arb_set(flint::ArbRef(candidate->regulator), regulator);
    candidate->work_precision = work_precision;
    candidate->exact_kernel_verified = true;
    out = std::move(candidate);
    return true;
}

bool relation_kernel_independent_unit_count(
        slong& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong precision) noexcept {
    out = 0;
    if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                         precision)) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }
    if (rank <= 0) {
        return true;
    }

    std::vector<FactoredElement> selected;
    if (!select_independent_relation_kernel_units(
                selected, order, class_group, embeddings, precision, rank)) {
        return false;
    }

    out = static_cast<slong>(selected.size());
    return true;
}

RegulatorPivotOutcome hnf_independent_unit_count(
        slong& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept {
    out = 0;
    if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                         precision)) {
        return RegulatorPivotOutcome::invalid;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return RegulatorPivotOutcome::invalid;
    }
    if (rank <= 0) {
        return RegulatorPivotOutcome::success;
    }

    slong places = 0;
    Signature sig;
    if (!compact_places(places, embeddings) ||
        !signature(sig, *order.parent())) {
        return RegulatorPivotOutcome::invalid;
    }

    flint::ArbMat unit_logs(
            class_group.relation_count() - class_group.relation_rank(),
            places);
    std::vector<slong> selected;
    if (!hnf_unit_log_matrix(unit_logs, class_group, embeddings,
                                  workspace, precision)) {
        // Keep the opaque upstream bool terminal for the same reason as the
        // product consumer above; selector uncertainty remains retryable.
        return RegulatorPivotOutcome::invalid;
    }
    const RegulatorPivotOutcome pivot_outcome =
            regulator_pivot_unit_indices(
                    selected, unit_logs, sig, rank, places, precision);
    if (pivot_outcome != RegulatorPivotOutcome::success) {
        return pivot_outcome;
    }

    out = static_cast<slong>(selected.size());
    return RegulatorPivotOutcome::success;
}

namespace {

struct ExactHnfRelationUnits {
    flint::Arb expected_regulator;
    flint::FmpzMat coefficients{0, 0};
    slong precision = 0;
};

bool exact_hnf_relation_units(
        ExactHnfRelationUnits& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        HnfFinishWorkspace* workspace,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    flint::Arb candidate_hR;
    flint::Arb expected_regulator;
    flint::FmpzMat witness_coefficients(0, 0);
    flint::FmpzMat unit_matrix(0, 0);
    flint::FmpzMat transformed_coefficients(0, 0);
    slong independent_count = 0;
    const RegulatorPivotOutcome product_outcome =
            hnf_class_regulator_product_and_unit_matrix(
                    flint::ArbRef(candidate_hR), &independent_count,
                    order, class_group, embeddings, analytic_hR,
                    &unit_matrix, &expected_regulator, workspace, precision);
    slong target_rank = -1;
    if (product_outcome != RegulatorPivotOutcome::success ||
        order.parent() == nullptr ||
        !unit_rank(target_rank, *order.parent()) ||
        independent_count != target_rank) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: reconstruct_regulator matrix");
        return false;
    }

    const bool have_witness_coefficients = workspace != nullptr
            ? hnf_finish_workspace_witness_coefficients(
                      witness_coefficients, *workspace, class_group)
            : hnf_unit_witness_coefficients(
                      witness_coefficients, class_group);
    if (!have_witness_coefficients) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: witness extraction");
        return false;
    }
    if (!transformed_hnf_unit_coefficients_from_regulator_matrix(
                transformed_coefficients, class_group, witness_coefficients,
                unit_matrix, diagnostics)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: transform witnesses");
        return false;
    }

    ExactHnfRelationUnits candidate;
    candidate.coefficients = flint::FmpzMat(
            flint::fmpz_mat_nrows(transformed_coefficients),
            flint::fmpz_mat_ncols(transformed_coefficients));
    flint::fmpz_mat_set(flint::FmpzMatRef(candidate.coefficients),
                        flint::FmpzMatConstRef(transformed_coefficients));
    candidate.precision = precision;
    if (!reduce_relation_coefficients_by_log_lll(
                candidate.coefficients, class_group, embeddings,
                candidate.precision, diagnostics)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: exact coefficient log LLL");
        return false;
    }
    flint::arb_set(flint::ArbRef(candidate.expected_regulator),
                   flint::ArbConstRef(expected_regulator));
    out = std::move(candidate);
    return true;
}

}  // namespace

bool set_hnf_units(OrderUnitGroup& out,
                        const Order& order,
                        const ClassGroupContext& class_group,
                        EmbeddingContext& embeddings,
                        flint::ArbConstRef analytic_hR,
                        HnfFinishWorkspace* workspace,
                        slong precision) noexcept {
    if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                         precision) ||
        !flint::arb_is_finite(analytic_hR) ||
        !flint::arb_is_positive(analytic_hR)) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: invalid inputs");
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: unit rank unavailable");
        return false;
    }
    if (rank == 0) {
        return out.compute(order);
    }

    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    ExactHnfRelationUnits exact;
    if (!exact_hnf_relation_units(
                exact, order, class_group, embeddings, analytic_hR,
                workspace, precision, out.diagnostics())) {
        return false;
    }

    std::vector<FactoredElement> transformed;
    const flint::Fmpz* cached_torsion_order = nullptr;
    const Element* cached_field_torsion_generator = nullptr;
    OrderElement cached_torsion_generator(order);
    if (cached_torsion_generator.is_defined() &&
        class_group_finish_torsion(
                cached_torsion_order, cached_field_torsion_generator,
                class_group) &&
        cached_torsion_generator.set_element(
                *cached_field_torsion_generator)) {
        SILEX_PROFILE_EVENT(
                out.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.hnf_units.torsion_from_finish");
    } else {
        cached_torsion_order = nullptr;
    }
    if (!factored_units_from_relation_coefficients(
                       transformed, class_group, exact.coefficients,
                       out.diagnostics())) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: exact coefficient units");
    } else if (!unit_log_row_sums_are_small(
                       transformed, embeddings, exact.precision)) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit setup failed: exact validate_unit_log_sums");
    } else {
        if (publish_validated_hnf_units(
                    out, order, transformed, embeddings,
                    flint::ArbConstRef(exact.expected_regulator), rank,
                    exact.precision, cached_torsion_order,
                    cached_torsion_order == nullptr
                            ? nullptr
                            : &cached_torsion_generator)) {
            return true;
        }
    }

    return false;
}

bool set_hnf_compact_relation_units(
        CompactRelationUnitGroup& out,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::ArbConstRef analytic_hR,
        HnfFinishWorkspace* workspace,
        slong precision) noexcept {
    const DiagnosticsContext* const diagnostics = class_group.diagnostics();
    if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                         precision) ||
        !flint::arb_is_finite(analytic_hR) ||
        !flint::arb_is_positive(analytic_hR)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "compact HNF unit setup failed: invalid inputs");
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "compact HNF unit setup failed: unit rank unavailable");
        return false;
    }

    if (rank == 0) {
        flint::Arb regulator;
        flint::arb_one(regulator);
        flint::FmpzMat coefficients(0, class_group.relation_count());
        std::unique_ptr<CompactRelationUnitGroupStorage> storage;
        if (!build_compact_relation_unit_storage(
                    storage, order, class_group, coefficients,
                    flint::ArbConstRef(regulator), precision)) {
            return false;
        }
        out.storage_ = std::move(storage);
        return true;
    }

    ExactHnfRelationUnits exact;
    if (!exact_hnf_relation_units(
                exact, order, class_group, embeddings, analytic_hR,
                workspace, precision, diagnostics)) {
        return false;
    }
    if (flint::fmpz_mat_nrows(exact.coefficients) != rank ||
        !relation_coefficients_are_exact_unit_kernel(
                class_group,
                flint::FmpzMatConstRef(exact.coefficients))) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "compact HNF unit setup failed: exact relation kernel");
        return false;
    }

    slong places = 0;
    if (!compact_places(places, embeddings) || places != rank + 1) {
        return false;
    }

    slong work_precision = exact.precision;
    for (slong attempt = 0; attempt < 4; ++attempt) {
        flint::ArbMat logs(rank, places);
        flint::Arb regulator;
        if (!embeddings.refine(work_precision, diagnostics) ||
            !relation_coefficients_log_matrix(
                    logs, class_group, embeddings, exact.coefficients,
                    work_precision, diagnostics)) {
            SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "compact HNF unit setup failed: logarithmic embedding");
            return false;
        }
        if (unit_log_row_sums_are_small(logs, work_precision) &&
            compact_regulator_from_log_matrix(
                    flint::ArbRef(regulator), logs, rank, places,
                    work_precision) &&
            arb_radius_lt_2exp(regulator, -64) &&
            unit_regulator_matches_reconstruction(
                    flint::ArbConstRef(regulator),
                    flint::ArbConstRef(exact.expected_regulator),
                    work_precision)) {
            std::unique_ptr<CompactRelationUnitGroupStorage> storage;
            if (!build_compact_relation_unit_storage(
                        storage, order, class_group, exact.coefficients,
                        flint::ArbConstRef(regulator), work_precision)) {
                SILEX_LOG(
                        diagnostics, DiagnosticsModule::unit_group,
                        LogLevel::detail,
                        "compact HNF unit setup failed: publication");
                return false;
            }
            out.storage_ = std::move(storage);
            return true;
        }

        if (attempt == 3 ||
            work_precision > std::numeric_limits<slong>::max() / 2) {
            break;
        }
        work_precision *= 2;
    }

    SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
              LogLevel::detail,
              "compact HNF unit setup failed: regulator validation");
    return false;
}

bool set_initial_relation_kernel_units(OrderUnitGroup& out,
                                       const Order& order,
                                       const ClassGroupContext& class_group,
                                       EmbeddingContext& embeddings,
                                       slong precision) noexcept {
    SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.set_initial_relation_kernel_units");
    if (!validate_relation_kernel_inputs(order, class_group, embeddings,
                                         precision)) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }
    if (rank == 0) {
        return out.compute(order);
    }

    const slong n = class_group.relation_kernel_unit_count();
    if (n < rank) {
        return false;
    }
    if (rank == 1) {
        return set_best_rank_one_relation_kernel_unit(
                out, order, class_group, embeddings, precision);
    }

    const NumberField* field = order.parent();
    std::vector<FactoredElement> selected;
    if (field == nullptr ||
        !select_independent_relation_kernel_units(
                selected, order, class_group, embeddings, precision, rank)) {
        return false;
    }

    if (static_cast<slong>(selected.size()) != rank) {
        return false;
    }

    return order_unit_group_set_units_internal(
            out, order, FactoredElementSpan(selected.data(), selected.size()),
            embeddings, precision, true);
}

}  // namespace detail
}  // namespace silex
