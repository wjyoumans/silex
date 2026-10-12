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
#include <flint/nmod.h>
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
    if (flint::fmpz_fits_si(exponent)) {
        return power.pow_si(base, flint::fmpz_get_si(exponent)) &&
               accumulator.multiply(accumulator, power);
    }

    Element base_value(*accumulator.parent());
    Element powered_value(*accumulator.parent());
    return base.evaluate(base_value) &&
           powered_value.pow_fmpz(base_value, exponent) &&
           power.set_element(powered_value) &&
           accumulator.multiply(accumulator, power);
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

// Exact power-root check of the relation row (-a_1, ..., -a_r, m), m >= 2:
// for i = 0, ..., w - 1 (w the torsion order) it tests whether
// prod u_i^a_i * zeta^i is an exact m-th power.  With require_y_root it
// accepts such a root only when it equals y up to torsion, and then sets
// root = y.  Without it, the first exact m-th root found is accepted and
// returned in `root`, whether or not it is y up to torsion; the relation
// then certifies root^m in <U> * mu, not y^m.  That mode differs from the
// reference: v0.38.6 Unit/UnitGrpCtx.jl `_add_dependent_unit!`
// always decides with `_check_relation_mod_torsion(U.units, y, rel, p)`
// (line 97) and adjoins y itself, `_transform(vcat(U.units, y), m)`
// (line 131).  The root is an exact field element and is adjoined only
// after adjoinable_root_is_order_unit, so the result is still a subgroup of
// the order's units.  When the row is a true relation for y, any such root
// equals y up to torsion; when it is a wrong guess, the group may fail to
// contain y.
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

bool kernel_rows_torsion_column(bool& has_torsion_column,
                                const OrderUnitGroup& group,
                                flint::FmpzMatConstRef kernel_rows,
                                flint::FmpzConstRef ell) noexcept {
    has_torsion_column = false;
    const slong rank = group.free_rank();
    const slong columns = flint::fmpz_mat_ncols(kernel_rows);
    if (!group.is_set() || rank < 0) {
        return false;
    }
    if (columns == rank) {
        return true;
    }
    slong proof_rank = 0;
    if (columns != rank + 1 ||
        !dlog_proof_rank(proof_rank, group, ell) || proof_rank != columns) {
        return false;
    }
    has_torsion_column = true;
    return true;
}

bool kernel_row_product(FactoredElement& product,
                        const OrderUnitGroup& group,
                        flint::FmpzMatConstRef kernel_rows,
                        slong row) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    const slong columns = flint::fmpz_mat_ncols(kernel_rows);
    if (!group.is_set() || field == nullptr || product.parent() == nullptr ||
        !product.parent()->has_same_data(*field) ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel_rows) ||
        (columns != rank && columns != rank + 1)) {
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
    if (columns == rank + 1) {
        // Torsion coordinate of a proof-kernel row (see
        // residue_dlog_proof_matrix): multiply by zeta^t, as reference
        // `saturate!` adds zeta with the last candidate coordinate.
        flint::FmpzConstRef exponent =
                flint::fmpz_mat_entry(kernel_rows, row, rank);
        if (!flint::fmpz_is_zero(exponent)) {
            OrderElement torsion(*order);
            Element torsion_value(*field);
            FactoredElement torsion_factor(*field);
            if (!group.torsion_generator(torsion) ||
                !torsion.get_element(torsion_value) ||
                !torsion_factor.set_element(torsion_value) ||
                !compact_multiply_power(product, torsion_factor, exponent)) {
                return false;
            }
        }
    }
    return true;
}

bool adjoinable_root_is_order_unit(const Order& order,
                                   const FactoredElement& root) noexcept {
    if (order.is_maximal()) {
        return true;
    }
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }
    Element expanded(*field);
    return expanded.is_defined() && root.evaluate(expanded) &&
           evaluated_is_order_unit(order, expanded);
}

namespace {

// Whether Element::is_power (src/element/element.cpp) can find an
// `exponent`-th root of a non-rational element of `field`.  Exponent 2 goes
// to Element::is_square.  For an exponent >= 3 the only route that finds a
// root of a non-rational input is pure_power_hensel_root
// (src/element/hensel.cpp), which returns "unsupported" for a field of
// degree >= 10 before it looks at the input; there
// pure_power_residue_disproves can only prove that the input is not a
// power.  So for exponent >= 3 and degree >= 10, Element::is_power finds a
// root only of a rational constant (is_power_rational_constant).  The bound
// must follow the pure_power_hensel_root guard; a test pins it.
constexpr slong kElementPowerRouteDegreeLimit = 10;

bool element_power_test_roots_non_rational(const NumberField& field,
                                           slong exponent) noexcept {
    return exponent <= 2 || field.degree() < kElementPowerRouteDegreeLimit;
}

// Whether every factor exponent of `product` is divisible by `exponent`, so
// compact_unit_presentation_power_root and FactoredElement::is_power_si
// both return the formal root at once without calling Element::is_power
// (their formal splits leave no leftover exactly in this case).
bool factored_is_formal_power(const FactoredElement& product,
                              slong exponent) noexcept {
    for (const auto& entry : product.factors()) {
        if (entry.exponent % exponent != 0) {
            return false;
        }
    }
    return true;
}

}  // namespace

KernelRowRootResult kernel_row_root(
        bool& is_power,
        FactoredElement& root,
        const OrderUnitGroup& group,
        flint::FmpzMatConstRef kernel_rows,
        slong row,
        flint::FmpzConstRef ell,
        CompactFieldModulusCache* field_modulus_cache,
        KernelRowRootMode mode) noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (!group.is_set() || field == nullptr || root.parent() == nullptr ||
        !root.parent()->has_same_data(*field) ||
        row < 0 || row >= flint::fmpz_mat_nrows(kernel_rows) ||
        !flint::fmpz_is_prime(ell) || !flint::fmpz_fits_si(ell)) {
        return KernelRowRootResult::error;
    }
    bool has_torsion_column = false;
    if (!kernel_rows_torsion_column(has_torsion_column, group, kernel_rows,
                                    ell)) {
        return KernelRowRootResult::error;
    }

    FactoredElement product(*field);
    FactoredElement candidate(*field);
    {
        SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.kernel_row_root_product");
        if (!kernel_row_product(product, group, kernel_rows, row)) {
            return KernelRowRootResult::error;
        }
    }

    const slong exponent = flint::fmpz_get_si(ell);
    if (mode == KernelRowRootMode::defer_undecidable &&
        !element_power_test_roots_non_rational(*field, exponent) &&
        !factored_is_formal_power(product, exponent)) {
        // Both root routes below end in Element::is_power, which here finds
        // a root only of a rational l-th power remainder.  The attempt can
        // end only in a non-power, in `undecided` (the caller rescans after
        // either, as after `deferred`), or in the root of a true l-th
        // power, which the caller still meets at a later scan.  Skip the
        // compact presentation and the factored expansion.
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.kernel_row_root_deferred");
        return KernelRowRootResult::deferred;
    }
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
            SILEX_PROFILE_EVENT(group.diagnostics(),
                                DiagnosticsModule::unit_group,
                                "unit_group.kernel_row_root_unavailable");
            return KernelRowRootResult::undecided;
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
        return KernelRowRootResult::decided;
    }

    return root.set(candidate) ? KernelRowRootResult::decided
                               : KernelRowRootResult::error;
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

namespace {

// One pending entry of a residue_dlog_proof_matrix column on the word-size
// degree-one path: `image` is the residue value raised to (p - 1) / ell.
struct WordQuotientLogTarget {
    ulong image = 0;
    slong row = 0;
};

// One baby step of word_quotient_log_baby_giant: zeta^exponent.
struct WordQuotientLogBabyStep {
    ulong value = 0;
    ulong exponent = 0;
};

// Reference `mod_p` (NumFieldOrd/NfOrd/Clgp/Saturate.jl:6-46 at the pinned
// v0.38.6 revision) fills one table of all generator powers per prime when
// `nbits(pp) < 19` (line 32) and otherwise takes each logarithm by
// `baby_step_giant_step` (line 43), whose baby-step table is shared by every
// element of the prime.  Deviation (cost only): the walk here is linear in
// ell, the order of the quotient generator, and baby-step giant-step is
// cheaper from ell of 7 bits on (busy-machine probe data on proven
// saturation rows, where ell ran up to 16 bits; not a controlled
// measurement), so baby-step giant-step is taken from 7 bits instead of 19.
// Both routes return the unique logarithm in [0, ell), so the values are
// the same.
constexpr flint_bitcnt_t kWordQuotientLogBabyGiantBits = 7;

// Quotient logs of `pending` by baby-step giant-step with one baby-step
// table, as reference `baby_step_giant_step`
// (NumFieldOrd/NfOrd/ResidueRingMultGrp.jl:716-738 at the same revision): m =
// isqrt(ell) + 1, baby steps zeta^0, ..., zeta^m, and giant steps by
// zeta^(-m) for i = 0, ..., m - 1, the logarithm being i m + j mod ell.  The
// caller has checked that ell is prime and zeta has order ell, so the
// logarithm is the unique exponent in [0, ell) and equals the one the walk
// returns.  An entry not found is left unresolved for the caller.
void word_quotient_log_baby_giant(
        flint::FmpzMat& out,
        std::vector<bool>& resolved,
        const std::vector<WordQuotientLogTarget>& pending,
        ulong zeta,
        ulong p_ui,
        ulong ell_ui,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_proof_matrix.baby_giant");
    nmod_t mod;
    nmod_init(&mod, p_ui);
    nmod_t ell_mod;
    nmod_init(&ell_mod, ell_ui);
    const ulong m = n_sqrt(ell_ui) + 1;
    std::vector<WordQuotientLogBabyStep> baby_steps;
    baby_steps.reserve(static_cast<std::size_t>(m) + 1);
    ulong current = 1 % p_ui;
    for (ulong j = 0; j <= m; ++j) {
        WordQuotientLogBabyStep step;
        step.value = current;
        step.exponent = j;
        baby_steps.push_back(step);
        current = nmod_mul(current, zeta, mod);
    }
    // zeta^0, ..., zeta^m are distinct because m < ell, the order of zeta.
    std::sort(baby_steps.begin(), baby_steps.end(),
              [](const WordQuotientLogBabyStep& left,
                 const WordQuotientLogBabyStep& right) {
                  return left.value < right.value;
              });
    // zeta^(-m) = zeta^(ell - m) since zeta has order ell.
    const ulong giant = nmod_pow_ui(zeta, ell_ui - m, mod);
    for (const WordQuotientLogTarget& target : pending) {
        ulong y = target.image;
        for (ulong i = 0; i < m; ++i) {
            const auto found = std::lower_bound(
                    baby_steps.begin(), baby_steps.end(), y,
                    [](const WordQuotientLogBabyStep& step, ulong value) {
                        return step.value < value;
                    });
            if (found != baby_steps.end() && found->value == y) {
                // i < m and found->exponent <= m are both below ell.
                const ulong exponent = nmod_add(nmod_mul(i, m, ell_mod),
                                                found->exponent, ell_mod);
                flint::fmpz_set_ui(
                        flint::fmpz_mat_entry(out, target.row, 0), exponent);
                resolved[static_cast<std::size_t>(target.row)] = true;
                break;
            }
            y = nmod_mul(y, giant, mod);
        }
    }
}

// The quotient logs of `pending` (images x^((p - 1) / ell) mod p) with
// respect to zeta, written into the rows of `out` that each target names and
// marked in `resolved`: the first exponent i in [0, ell) with zeta^i equal
// to the image, found by one walk over the powers of zeta shared by every
// target or, from kWordQuotientLogBabyGiantBits on, when ell is prime and
// zeta has order ell, by word_quotient_log_baby_giant, which gives the same
// unique exponent.  A target with no match below ell is left unresolved.
void word_quotient_log_resolve(
        flint::FmpzMat& out,
        std::vector<bool>& resolved,
        std::vector<WordQuotientLogTarget> pending,
        ulong zeta,
        ulong p_ui,
        ulong ell_ui,
        const DiagnosticsContext* diagnostics) noexcept {
    if (pending.empty()) {
        return;
    }

    const ulong p_inverse = n_preinvert_limb(p_ui);
    if (FLINT_BIT_COUNT(ell_ui) >= kWordQuotientLogBabyGiantBits &&
        zeta != 1 && n_is_prime(ell_ui) != 0 &&
        n_powmod2_ui_preinv(zeta, ell_ui, p_ui, p_inverse) == 1) {
        word_quotient_log_baby_giant(out, resolved, pending, zeta, p_ui,
                                     ell_ui, diagnostics);
        return;
    }

    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                        "unit_group.residue_dlog_proof_matrix.shared_walk");
    nmod_t mod;
    nmod_init(&mod, p_ui);
    ulong current = 1 % p_ui;
    for (ulong i = 0; i < ell_ui && !pending.empty(); ++i) {
        std::size_t k = 0;
        while (k < pending.size()) {
            if (pending[k].image != current) {
                ++k;
                continue;
            }
            flint::fmpz_set_ui(
                    flint::fmpz_mat_entry(out, pending[k].row, 0), i);
            resolved[static_cast<std::size_t>(pending[k].row)] = true;
            pending[k] = pending.back();
            pending.pop_back();
        }
        current = nmod_mul(current, zeta, mod);
    }
}

// Quotient logs of all `images` with respect to `quotient_log`, sharing work
// across the column.  Each entry this resolves is the value
// ResidueFieldQuotientLog::apply returns on its word-size degree-one path
// (`quotient_log_mod_prime_apply_ui` in residue_field.cpp), which walks
// zeta^0, zeta^1, ..., zeta^(ell - 1) for every element separately and
// returns the first exponent i with zeta^i equal to the element raised to
// (p - 1) / ell.  Below kWordQuotientLogBabyGiantBits one walk visits the
// same powers in the same order once and records, for each element, the
// first such i; the walk stops early once every element is matched.  At or
// above it, when ell is prime and zeta has order ell, the
// logarithms are taken by word_quotient_log_baby_giant instead; the
// logarithm is then unique in [0, ell), so the values are the same.
// `resolved[j]` is set for the entries written into row j of `out`; an
// entry this cannot resolve (field not of degree one, a value off the
// word-size path, a zero residue, or no match below ell) is left unset for
// the caller to pass to apply.
void word_quotient_log_shared_walk(
        flint::FmpzMat& out,
        std::vector<bool>& resolved,
        const std::vector<ResidueFieldElement>& images,
        const ResidueField& residue_field,
        const ResidueFieldQuotientLog& quotient_log,
        flint::FmpzConstRef ell,
        const DiagnosticsContext* diagnostics) noexcept {
    flint::Fmpz p;
    flint::Fmpz scalar;
    ResidueFieldElement quotient_generator(residue_field);
    if (residue_field.degree() != 1 ||
        !residue_field.characteristic(flint::FmpzRef(p)) ||
        !flint::fmpz_abs_fits_ui(flint::FmpzConstRef(p)) ||
        !flint::fmpz_abs_fits_ui(ell) ||
        !quotient_log.quotient_generator(quotient_generator) ||
        !quotient_generator.degree_one_scalar(flint::FmpzRef(scalar)) ||
        !flint::fmpz_abs_fits_ui(flint::FmpzConstRef(scalar))) {
        return;
    }

    const ulong p_ui = flint::fmpz_get_ui(flint::FmpzConstRef(p));
    const ulong ell_ui = flint::fmpz_get_ui(ell);
    if (p_ui < 2 || ell_ui < 2 || (p_ui - 1) % ell_ui != 0) {
        return;
    }
    const ulong zeta = flint::fmpz_get_ui(flint::FmpzConstRef(scalar)) % p_ui;
    if (zeta == 0) {
        return;
    }

    // The quotient-log context stores cofactor (q - 1) / ell, and q = p on a
    // degree-one residue field.
    const ulong cofactor = (p_ui - 1) / ell_ui;
    const ulong p_inverse = n_preinvert_limb(p_ui);
    std::vector<WordQuotientLogTarget> pending;
    pending.reserve(images.size());
    for (std::size_t j = 0; j < images.size(); ++j) {
        if (!images[j].degree_one_scalar(flint::FmpzRef(scalar)) ||
            !flint::fmpz_abs_fits_ui(flint::FmpzConstRef(scalar))) {
            continue;
        }
        const ulong value =
                flint::fmpz_get_ui(flint::FmpzConstRef(scalar)) % p_ui;
        if (value == 0) {
            continue;
        }
        WordQuotientLogTarget target;
        target.image = n_powmod2_ui_preinv(value, cofactor, p_ui, p_inverse);
        target.row = static_cast<slong>(j);
        pending.push_back(target);
    }
    word_quotient_log_resolve(out, resolved, std::move(pending), zeta, p_ui,
                              ell_ui, diagnostics);
}

}  // namespace

bool residue_quotient_log_column(
        flint::FmpzMat& out,
        const std::vector<ResidueFieldElement>& images,
        const ResidueField& residue_field,
        const ResidueFieldQuotientLog& quotient_log,
        flint::FmpzConstRef ell,
        const DiagnosticsContext* diagnostics) noexcept {
    if (flint::fmpz_mat_nrows(out) != static_cast<slong>(images.size()) ||
        flint::fmpz_mat_ncols(out) != 1) {
        return false;
    }
    // Every entry takes the quotient log of the same residue field and ell,
    // so one walk over the powers of the quotient generator, or one
    // baby-step table from kWordQuotientLogBabyGiantBits on, serves the
    // whole column; entries left unresolved take ResidueFieldQuotientLog::apply as
    // before.
    std::vector<bool> resolved(images.size(), false);
    word_quotient_log_shared_walk(out, resolved, images, residue_field,
                                  quotient_log, ell, diagnostics);
    for (std::size_t j = 0; j < images.size(); ++j) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::unit_group,
                            "unit_group.residue_dlog_proof_matrix.entry");
        if (resolved[j]) {
            continue;
        }
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.residue_dlog_proof_matrix.entry_apply");
        if (!quotient_log.apply(
                    flint::fmpz_mat_entry(out, static_cast<slong>(j), 0),
                    images[j])) {
            return false;
        }
    }
    return true;
}

namespace {

// residue_dlog_proof_matrix at `prime`.  `prime_field`, when not null, is
// the residue field the caller already built from `prime`; it is used in
// place of building the same field again.
bool residue_dlog_proof_matrix_with_field(
        flint::FmpzMat& out,
        const OrderUnitGroup& group,
        const PrimeIdeal& prime,
        const ResidueField* prime_field,
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

    ResidueField local_field;
    if (prime_field == nullptr) {
        ResidueField built(prime);
        local_field = std::move(built);
        prime_field = &local_field;
    }
    const ResidueField& residue_field = *prime_field;
    if (!residue_field.is_defined()) {
        return false;
    }

    flint::FmpzMat candidate(proof_rank, 1);
    ResidueFieldQuotientLog quotient_log(residue_field);
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }

    // Reduce every entry first: the free generators, then the torsion
    // generator when ell divides the torsion order.
    std::vector<ResidueFieldElement> images;
    images.reserve(static_cast<std::size_t>(proof_rank));
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        images.emplace_back(residue_field);
        if (!group.free_generator(generator, i) ||
            !images.back().set_factored_element(generator)) {
            return false;
        }
    }

    if (proof_rank > rank) {
        OrderElement torsion(*order);
        images.emplace_back(residue_field);
        if (!group.torsion_generator(torsion) ||
            !images.back().set_order_element(torsion)) {
            return false;
        }
    }

    if (!residue_quotient_log_column(candidate, images, residue_field,
                                     quotient_log, ell,
                                     group.diagnostics())) {
        return false;
    }

    out = std::move(candidate);
    return true;
}

}  // namespace

bool residue_dlog_proof_matrix(flint::FmpzMat& out,
                               const OrderUnitGroup& group,
                               const PrimeIdeal& prime,
                               flint::FmpzConstRef ell) noexcept {
    return residue_dlog_proof_matrix_with_field(out, group, prime, nullptr,
                                                ell);
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

namespace {

// saturation_proof_prime_column_at_degree_one_root in word arithmetic, for
// p, root and ell below the word size: the free generators are evaluated by
// factored_value_at_degree_one_root_nmod, which fails exactly where
// factored_value_at_degree_one_root does and otherwise gives the same
// residue, and the quotient logs of the column are resolved together by
// word_quotient_log_resolve with zeta = quotient_generator.  zeta is a
// primitive root raised to (p - 1) / ell, so it has order ell, and the first
// exponent in [0, ell) the walk finds is the unique one
// quotient_log_mod_prime_apply returns.  An entry left unresolved takes
// quotient_log_mod_prime_apply.  Returns false, leaving `candidate`
// unspecified, when the word route does not apply or an entry cannot be
// formed; the caller then runs the fmpz route.
bool degree_one_root_column_word(flint::FmpzMat& candidate,
                                 const OrderUnitGroup& group,
                                 flint::FmpzConstRef p,
                                 flint::FmpzConstRef root,
                                 flint::FmpzConstRef ell,
                                 flint::FmpzConstRef cofactor,
                                 flint::FmpzConstRef quotient_generator)
        noexcept {
    const Order* order = group.parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr || !flint::fmpz_abs_fits_ui(p) ||
        !flint::fmpz_abs_fits_ui(root) || !flint::fmpz_abs_fits_ui(ell) ||
        !flint::fmpz_abs_fits_ui(cofactor) ||
        !flint::fmpz_abs_fits_ui(quotient_generator) ||
        flint::fmpz_sgn(p) <= 0 || flint::fmpz_sgn(root) < 0 ||
        flint::fmpz_sgn(ell) <= 0) {
        return false;
    }
    const ulong p_ui = flint::fmpz_get_ui(p);
    const ulong root_ui = flint::fmpz_get_ui(root);
    const ulong ell_ui = flint::fmpz_get_ui(ell);
    const ulong cofactor_ui = flint::fmpz_get_ui(cofactor);
    const ulong zeta = flint::fmpz_get_ui(quotient_generator) % p_ui;
    const slong rank = group.free_rank();
    const slong rows = flint::fmpz_mat_nrows(candidate);
    if (p_ui < 2 || ell_ui < 2 || root_ui >= p_ui || zeta == 0 ||
        rows < rank) {
        return false;
    }

    const ulong p_inverse = n_preinvert_limb(p_ui);
    std::vector<WordQuotientLogTarget> pending;
    pending.reserve(static_cast<std::size_t>(rows));
    std::vector<ulong> values;
    values.reserve(static_cast<std::size_t>(rows));
    for (slong i = 0; i < rank; ++i) {
        FactoredElement generator(*field);
        ulong value = 0;
        if (!group.free_generator(generator, i) ||
            !factored_value_at_degree_one_root_nmod(value, generator, p_ui,
                                                    root_ui)) {
            return false;
        }
        values.push_back(value);
    }
    if (rows > rank) {
        OrderElement torsion(*order);
        Element torsion_value(*field);
        flint::Fmpz value;
        if (rows != rank + 1 || !group.torsion_generator(torsion) ||
            !torsion.get_element(torsion_value) ||
            !element_value_at_degree_one_root(value, torsion_value, p, root)) {
            return false;
        }
        values.push_back(flint::fmpz_fdiv_ui(flint::FmpzConstRef(value),
                                             p_ui));
    }
    for (std::size_t j = 0; j < values.size(); ++j) {
        if (values[j] == 0) {
            return false;
        }
        WordQuotientLogTarget target;
        target.image =
                n_powmod2_ui_preinv(values[j], cofactor_ui, p_ui, p_inverse);
        target.row = static_cast<slong>(j);
        pending.push_back(target);
    }

    std::vector<bool> resolved(values.size(), false);
    word_quotient_log_resolve(candidate, resolved, std::move(pending), zeta,
                              p_ui, ell_ui, group.diagnostics());
    flint::Fmpz value;
    for (std::size_t j = 0; j < values.size(); ++j) {
        if (resolved[j]) {
            continue;
        }
        flint::fmpz_set_ui(flint::FmpzRef(value), values[j]);
        if (!quotient_log_mod_prime_apply(
                    flint::fmpz_mat_entry(candidate, static_cast<slong>(j), 0),
                    flint::FmpzConstRef(value), cofactor, quotient_generator,
                    p, ell)) {
            return false;
        }
    }
    return true;
}

}  // namespace

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
    if (degree_one_root_column_word(candidate, group, p, root, ell,
                                    flint::FmpzConstRef(cofactor),
                                    flint::FmpzConstRef(quotient_generator))) {
        out = std::move(candidate);
        return true;
    }
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

bool residue_dlog_character_column(flint::FmpzMat& out,
                                   const OrderUnitGroup& group,
                                   const PrimeIdeal& prime,
                                   flint::FmpzConstRef ell) noexcept {
    slong proof_rank = 0;
    if (!group.is_set() || !dlog_proof_rank(proof_rank, group, ell)) {
        return false;
    }
    if (proof_rank == group.free_rank()) {
        return residue_dlog_matrix(out, group, prime, ell);
    }
    // l divides the torsion order: append the torsion generator, as
    // reference `compute_candidates_for_saturate` appends zeta to the
    // input units.  Both routes below evaluate the free generators and zeta
    // with one character, so the column's entries are comparable.
    if (saturation_proof_prime_column_direct_degree_one(out, group, prime,
                                                        ell)) {
        return true;
    }
    return residue_dlog_proof_matrix(out, group, prime, ell);
}

bool residue_dlog_kernel(flint::FmpzMat& out,
                         const OrderUnitGroup& group,
                         PrimeIdealSpan primes,
                         flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.detail_residue_dlog_kernel");
    const Order* order = group.parent();
    slong width = 0;
    if (!group.is_set() || order == nullptr || primes.empty() ||
        !dlog_proof_rank(width, group, ell)) {
        return false;
    }

    flint::FmpzMat matrix(width, static_cast<slong>(primes.size()));
    for (std::size_t j = 0; j < primes.size(); ++j) {
        if (!same_order_parent(primes[j].parent(), order)) {
            return false;
        }
        flint::FmpzMat column(width, 1);
        if (!residue_dlog_character_column(column, group, primes[j], ell) ||
            flint::fmpz_mat_nrows(column) != width) {
            return false;
        }
        for (slong i = 0; i < width; ++i) {
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

bool residue_dlog_free_kernel(flint::FmpzMat& out,
                              const OrderUnitGroup& group,
                              PrimeIdealSpan primes,
                              flint::FmpzConstRef ell) noexcept {
    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.detail_residue_dlog_free_kernel");
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

bool saturation_prime_character_column(flint::FmpzMat& out,
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

    return residue_dlog_character_column(out, group, prime, ell);
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

    return residue_dlog_proof_matrix_with_field(out, group, prime,
                                                &residue_field, ell);
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

    // The new basis is built from `root` and the old generators and is
    // installed without the per-generator unit check.  The root is an
    // e-th root in the field of a product of order units times a root of
    // unity, hence a unit of O_K; in a non-maximal order it may still lie
    // outside O, and then it must not be adjoined.
    if (!adjoinable_root_is_order_unit(*order, root)) {
        SILEX_PROFILE_EVENT(group.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.dependent_relation_root_outside_order");
        changed = false;
        return out.set(group);
    }

    std::vector<FactoredElement> generators;
    if (!relation_basis(generators, group, root,
                        flint::FmpzMatConstRef(relation),
                        flint::FmpzConstRef(exponent))) {
        return false;
    }

    // The torsion is that of the order, which the new basis does not
    // change, so the group's torsion is passed as the cache instead of
    // recomputing it; set_units still checks it with
    // cached_torsion_is_valid.
    flint::Fmpz torsion_order;
    OrderElement torsion_generator(*order);
    if (!torsion_generator.is_defined() ||
        !group.torsion_order(flint::FmpzRef(torsion_order)) ||
        !group.torsion_generator(torsion_generator)) {
        return false;
    }

    changed = true;
    return order_unit_group_set_units_internal(
            out, *order,
            FactoredElementSpan(generators.data(), generators.size()),
            embeddings, precision, true, &torsion_order, &torsion_generator);
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

// Exactly verifies the candidate row (-a_1, ..., -a_r, m), m >= 1.  A
// candidate obtained from rounded or reconstructed coordinates is only a
// guess; this check is what makes a recovered relation usable.  For m = 1,
// or with require_y_root, it verifies y^m = prod u_i^a_i modulo torsion and
// returns root = y.  For m >= 2 without require_y_root (only
// OrderUnitGroup::set_relation_kernel_units_bounded passes false) it first
// accepts any exact m-th root of prod u_i^a_i * zeta^i, returned in `root`
// and not compared with y (see relation_power_root; a Silex divergence from
// reference `_add_dependent_unit!`, which always checks and adjoins y).  When
// that search cannot complete or finds no root, it checks y modulo torsion
// and, if that verifies the relation, returns root = y.  So `root` is always
// the element whose relation was verified; its value on entry is never
// returned.
bool verify_dependent_relation_candidate(bool& recovered,
                                         FactoredElement& root,
                                         flint::FmpzMat& rel,
                                         flint::Fmpz& torsion_exp,
                                         const OrderUnitGroup& group,
                                         const FactoredElement& y,
                                         EmbeddingContext& embeddings,
                                         const flint::FmpzMat& candidate,
                                         bool require_y_root,
                                         bool require_torsion_exponent) noexcept {
    const slong rank = group.free_rank();
    const flint::FmpzConstRef last =
            flint::fmpz_mat_entry(candidate, 0, rank);
    if (flint::fmpz_sgn(last) <= 0 || !flint::fmpz_fits_si(last)) {
        return false;
    }
    const slong denominator = flint::fmpz_get_si(last);

    bool is_relation = false;
    bool root_is_verified_power_root = false;
    if (denominator >= 2 && require_y_root) {
        // reference Unit/UnitGrpCtx.jl `_add_dependent_unit!` decides the
        // candidate with `_check_relation_mod_torsion` (Unit/Relation.jl)
        // and adjoins y itself; it extracts no m-th root.  The root path
        // below would only confirm the same relation (its root must equal y
        // up to torsion, and y is adjoined either way), so it runs only when
        // the torsion check fails to complete.
        if (!check_relation_mod_torsion(
                    is_relation, flint::FmpzRef(torsion_exp), group, y,
                    flint::FmpzMatConstRef(candidate), 0, &embeddings,
                    require_torsion_exponent)) {
            flint::Fmpz candidate_exp;
            if (!relation_power_root(
                        is_relation, root, candidate_exp, group, y,
                        flint::FmpzMatConstRef(candidate), 0, true) ||
                !is_relation) {
                return false;
            }
            flint::fmpz_zero(flint::FmpzRef(torsion_exp));
        }
    } else if (denominator >= 2) {
        flint::Fmpz candidate_exp;
        const bool power_root_ok = relation_power_root(
                is_relation, root, candidate_exp, group, y,
                flint::FmpzMatConstRef(candidate), 0, require_y_root);
        if (power_root_ok && is_relation) {
            flint::fmpz_zero(flint::FmpzRef(torsion_exp));
            root_is_verified_power_root = true;
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
        // reference `_add_dependent_unit!` transforms `vcat(U.units, y)`.
        // The torsion check verifies the relation for y and does not write
        // `root`, which may still hold a value from the caller or from a
        // failed root search, so `root` is set to y here.  Only an exact
        // m-th root found by relation_power_root without require_y_root is
        // kept: its own relation is the one verified.  With require_y_root
        // that search already returns root = y.
        if (!root_is_verified_power_root && !root.set(y)) {
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

    return verify_dependent_relation_candidate(
            recovered, root, rel, torsion_exp, group, y, embeddings,
            candidate, require_y_root, require_torsion_exponent);
}

namespace {

// Exact endpoints lower / scale <= upper / scale of a finite ball, scale a
// positive power of two.  False when the ball is not finite or its exponent
// does not fit a slong.
bool ball_rational_interval(flint::Fmpz& lower,
                            flint::Fmpz& upper,
                            flint::Fmpz& scale,
                            const arb_struct* ball) noexcept {
    if (!arb_is_finite(ball)) {
        return false;
    }
    flint::Fmpz exponent;
    arb_get_interval_fmpz_2exp(lower.raw(), upper.raw(), exponent.raw(),
                               ball);
    if (!flint::fmpz_fits_si(flint::FmpzConstRef(exponent))) {
        return false;
    }
    const slong shift = flint::fmpz_get_si(flint::FmpzConstRef(exponent));
    flint::fmpz_one(flint::FmpzRef(scale));
    if (shift >= 0) {
        fmpz_mul_2exp(lower.raw(), lower.raw(),
                      static_cast<flint_bitcnt_t>(shift));
        fmpz_mul_2exp(upper.raw(), upper.raw(),
                      static_cast<flint_bitcnt_t>(shift));
    } else {
        fmpz_mul_2exp(scale.raw(), scale.raw(),
                      static_cast<flint_bitcnt_t>(-shift));
    }
    return true;
}

}  // namespace

bool dependent_relation_denominator_step(bool& possible,
                                         slong& step,
                                         const flint::ArbMat& coordinates,
                                         slong denominator_bound,
                                         slong precision) noexcept {
    possible = true;
    step = 1;
    const slong rank = flint::arb_mat_nrows_value(coordinates);
    if (rank <= 0 || flint::arb_mat_ncols_value(coordinates) != 1 ||
        denominator_bound <= 0 || precision <= 0) {
        return false;
    }

    // This only skips denominators that the scan below would skip anyway;
    // it never accepts one.  try_denominator(m) needs every ball of
    // m * coordinate (arb_mat_scalar_mul_si at `precision`) to contain an
    // integer k.  That ball lies in m * E, where E is the coordinate ball
    // x +- r widened to radius 2 r + (|x| + r) 2^(8 - precision): FLINT
    // arb_mul_arf rounds the midpoint m x by less than one ulp
    // (<= m |x| 2^(1 - precision)) and adds that ulp to the radius m r,
    // which three upward mag roundings enlarge by a factor below
    // 1 + 2^-27.  So k / m lies in E.  (FLINT src/arb/mul.c: arb_mul_si
    // calls arb_mul_arf, which rounds the midpoint with ARB_RND =
    // ARF_RND_DOWN and, when inexact, sets the radius of an exact input
    // ball to that ulp (arf_mag_set_ulp) or adds the ulp to the radius
    // (arf_mag_fast_add_ulp / arf_mag_add_ulp); arb_mat_scalar_mul_si
    // applies arb_mul_si entrywise, src/arb_mat/templates.c.  Read in FLINT
    // 3.3.1, 3.4.0, 3.5.0 and the 3.6.0-dev tree; identical in all four.)
    // When E is narrower than 1 / bound^2 it holds at most one rational with
    // denominator <= bound (two differ by at least 1 / bound^2), and then
    // that rational is the simplest one in E (the one with least
    // denominator, FLINT `_fmpq_simplest_between`).  Its denominator q must
    // divide every such m <= bound, and if q > bound no m qualifies.  The
    // same simplest-rational reconstruction is the reference
    // Unit/Relation.jl `_find_rational_relation!`; here it only narrows the
    // denominators that the unchanged exact scan tries.
    const ulong bound = static_cast<ulong>(denominator_bound);
    flint::Arb enlarged;
    flint::Arb scratch;
    flint::Fmpz lower;
    flint::Fmpz upper;
    flint::Fmpz scale;
    flint::Fmpz width;
    flint::Fmpq simplest;
    flint::Fmpz common;
    flint::fmpz_one(flint::FmpzRef(common));
    for (slong i = 0; i < rank; ++i) {
        const arb_struct* entry = arb_mat_entry(coordinates.raw(), i, 0);
        if (!arb_is_finite(entry)) {
            continue;
        }
        arb_set(enlarged.raw(), entry);
        mag_mul_2exp_si(arb_radref(enlarged.raw()),
                        arb_radref(enlarged.raw()), 1);
        arb_get_mag(arb_radref(scratch.raw()), entry);
        mag_mul_2exp_si(arb_radref(scratch.raw()), arb_radref(scratch.raw()),
                        8 - precision);
        mag_add(arb_radref(enlarged.raw()), arb_radref(enlarged.raw()),
                arb_radref(scratch.raw()));
        if (!ball_rational_interval(lower, upper, scale, enlarged.raw())) {
            continue;
        }

        // Narrow: (upper - lower) / scale < 1 / bound^2.
        fmpz_sub(width.raw(), upper.raw(), lower.raw());
        fmpz_mul_ui(width.raw(), width.raw(), bound);
        fmpz_mul_ui(width.raw(), width.raw(), bound);
        if (fmpz_cmp(width.raw(), scale.raw()) >= 0) {
            continue;
        }

        _fmpq_simplest_between(fmpq_numref(simplest.raw()),
                               fmpq_denref(simplest.raw()), lower.raw(),
                               scale.raw(), upper.raw(), scale.raw());
        fmpz_lcm(common.raw(), common.raw(), fmpq_denref(simplest.raw()));
        if (fmpz_cmp_ui(common.raw(), bound) > 0) {
            possible = false;
            return true;
        }
    }

    step = flint::fmpz_get_si(flint::FmpzConstRef(common));
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
        // The scan tries m = min_denominator, ..., denominator_bound in
        // order.  dependent_relation_denominator_step shows that only
        // multiples of `step` can pass try_denominator's integer test, so
        // the others are skipped; the denominators tried, their order and
        // the exact checks are those of the full scan.
        bool possible = true;
        slong step = 1;
        if (!dependent_relation_denominator_step(possible, step, coordinates,
                                                 denominator_bound,
                                                 precision)) {
            return false;
        }
        if (!possible) {
            SILEX_PROFILE_EVENT(
                    group.diagnostics(), DiagnosticsModule::unit_group,
                    "unit_group.dependent_relation_denominator_scan.empty");
            recovered = false;
            return true;
        }
        slong m = min_denominator;
        const slong remainder = m % step;
        if (remainder != 0) {
            if (m > denominator_bound - (step - remainder)) {
                recovered = false;
                return true;
            }
            m += step - remainder;
        }
        while (m <= denominator_bound) {
            if (!try_denominator(recovered, root, rel, torsion_exp, group, y,
                                 embeddings, coordinates, m, precision,
                                 require_y_root, require_torsion_exponent)) {
                return false;
            }
            if (recovered) {
                return true;
            }
            if (m > denominator_bound - step) {
                break;
            }
            m += step;
        }
    }

    recovered = false;
    return true;
}

bool dependent_relation_rational_candidate(bool& found,
                                           flint::FmpzMat& candidate,
                                           const flint::ArbMat& coordinates,
                                           slong denominator_bound) noexcept {
    found = false;
    const slong rank = flint::arb_mat_nrows_value(coordinates);
    if (rank <= 0 || flint::arb_mat_ncols_value(coordinates) != 1 ||
        denominator_bound <= 0 || flint::fmpz_mat_nrows(candidate) != 1 ||
        flint::fmpz_mat_ncols(candidate) != rank + 1) {
        return false;
    }

    // reference Unit/Relation.jl `_find_rational_relation!`: replace each
    // coordinate ball by the simplest rational inside it (`simplest_inside`,
    // Nemo `simplest_rational_inside` = FLINT `_fmpq_simplest_between` on the
    // exact ball endpoints), reject a ball with radius above one, and reject
    // the attempt when a denominator or their lcm exceeds the bound.  The
    // reference `simplest_inside` accepts denominators below the bound while
    // its lcm test accepts an lcm equal to it; Silex accepts the closed range
    // 1..bound.  Every denominator divides the least m with y^m in the group
    // modulo torsion, and m <= [U:G] <= bound, so m == bound can occur.
    //
    // Why the precision loop terminates: if the true coordinate is a/b with
    // b <= bound and the ball has radius rho < 1/(2 bound^2), any two points
    // of the ball are less than 1/bound^2 apart, while every other p/q with
    // q <= b satisfies |p/q - a/b| >= 1/(qb) >= 1/bound^2, so the simplest
    // rational in the ball is a/b.  Nothing here is trusted: the caller
    // verifies the candidate exactly.
    std::vector<flint::Fmpq> approximations(static_cast<std::size_t>(rank));
    flint::Fmpz lower;
    flint::Fmpz upper;
    flint::Fmpz scale;
    flint::Fmpz common;
    flint::fmpz_one(flint::FmpzRef(common));
    for (slong i = 0; i < rank; ++i) {
        const arb_struct* entry = arb_mat_entry(coordinates.raw(), i, 0);
        if (!arb_is_finite(entry) ||
            mag_cmp_2exp_si(arb_radref(entry), 0) > 0 ||
            !ball_rational_interval(lower, upper, scale, entry)) {
            return true;
        }

        flint::Fmpq& approximation =
                approximations[static_cast<std::size_t>(i)];
        _fmpq_simplest_between(fmpq_numref(approximation.raw()),
                               fmpq_denref(approximation.raw()), lower.raw(),
                               scale.raw(), upper.raw(), scale.raw());
        if (fmpz_cmp_si(fmpq_denref(approximation.raw()), denominator_bound) >
            0) {
            return true;
        }
        fmpz_lcm(common.raw(), common.raw(),
                 fmpq_denref(approximation.raw()));
        if (fmpz_cmp_si(common.raw(), denominator_bound) > 0) {
            return true;
        }
    }

    for (slong i = 0; i < rank; ++i) {
        const flint::Fmpq& approximation =
                approximations[static_cast<std::size_t>(i)];
        fmpz* out = flint::fmpz_mat_entry(candidate, 0, i).raw();
        fmpz_divexact(out, common.raw(), fmpq_denref(approximation.raw()));
        fmpz_mul(out, out, fmpq_numref(approximation.raw()));
        fmpz_neg(out, out);
    }
    fmpz_set(flint::fmpz_mat_entry(candidate, 0, rank).raw(), common.raw());
    found = true;
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
        slong precision,
        bool require_y_root,
        bool require_torsion_exponent) noexcept {
    const slong rank = group.free_rank();
    if (!group.is_set() || rank <= 0 || y.parent() == nullptr ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*y.parent()) ||
        flint::fmpz_sgn(denominator_bound) <= 0 ||
        !flint::fmpz_fits_si(denominator_bound) || precision <= 0 ||
        flint::arb_mat_nrows_value(inverse_cutoff) != rank ||
        flint::arb_mat_ncols_value(inverse_cutoff) != rank ||
        flint::fmpz_mat_nrows(rel) != 1 ||
        flint::fmpz_mat_ncols(rel) != rank + 1) {
        return false;
    }
    recovered = false;

    const slong bound = flint::fmpz_get_si(denominator_bound);
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

    flint::FmpzMat candidate(1, rank + 1);
    bool found = false;
    {
        SILEX_PROFILE_SCOPE(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.dependent_relation_rational_reconstruction");
        if (!dependent_relation_rational_candidate(found, candidate,
                                                   coordinates, bound)) {
            return false;
        }
    }
    if (!found) {
        SILEX_PROFILE_EVENT(
                group.diagnostics(), DiagnosticsModule::unit_group,
                "unit_group.dependent_relation_rational_reconstruction."
                "inconclusive");
        return true;
    }

    SILEX_PROFILE_SCOPE(group.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.dependent_relation_verify_candidate");
    return verify_dependent_relation_candidate(
            recovered, root, rel, torsion_exp, group, y, embeddings,
            candidate, require_y_root, require_torsion_exponent);
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
    // A `grh` unit group is accepted as well as a `proven` one.  The not-unit
    // answer comes from the unconditional principal-ideal test below, and a
    // verified answer is checked exactly against the given generators, so
    // either answer is exact relative to those generators.  Coordinates carry
    // no label of their own: a caller that publishes a label takes it from
    // the unit group.
    const CertificationMode label = group.certification_status();
    if (!group.is_set() ||
        (label != CertificationMode::proven &&
         label != CertificationMode::grh) ||
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
        // The pair is checked once here with the test the cached-torsion
        // hook applies to it, so a pair the hook would refuse fails the
        // scan now instead of rejecting every candidate below.  The hook
        // still repeats the check per candidate; w is capped by the degree,
        // so that costs a few small powers.
        if (!torsion_generator.is_defined() ||
            !rank_zero_torsion(flint::FmpzRef(torsion_order),
                               torsion_generator, order) ||
            !cached_torsion_is_valid(order,
                                     flint::FmpzConstRef(torsion_order),
                                     torsion_generator)) {
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
                                 precision, class_group.diagnostics())) {
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
        slong precision,
        HnfUnitPivotProbe* probe_out) noexcept {
    out = 0;
    if (probe_out != nullptr) {
        probe_out->clear();
    }
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
    if (probe_out != nullptr) {
        probe_out->unit_logs = std::move(unit_logs);
        probe_out->selected = std::move(selected);
        probe_out->precision = precision;
        probe_out->relation_count = class_group.relation_count();
        probe_out->relation_rank = class_group.relation_rank();
        probe_out->valid = true;
    }
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
