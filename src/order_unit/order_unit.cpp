#include <silex/order_unit.hpp>

#include "compute_internal.hpp"
#include "order_unit_internal.hpp"

#include <limits>
#include <vector>

#include <flint/arb_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/mag.h>

#include <silex/class_group.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/residue_field.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include <new>
#include <utility>

namespace silex {
namespace {

}  // namespace

namespace detail {

constexpr slong kRegulatorAbsTolerance = 64;

bool compute_power(Element& out,
                   const Element& base,
                   slong exponent) noexcept {
    flint::Fmpz exp;
    flint::fmpz_set_si(flint::FmpzRef(exp), exponent);
    return out.pow_fmpz(base, flint::FmpzConstRef(exp));
}

}  // namespace detail

namespace detail {

bool signature_from_embeddings(Signature& sig,
                               EmbeddingContext& embeddings) noexcept {
    if (embeddings.is_set()) {
        sig = embeddings.signature();
        return true;
    }
    const NumberField* field = embeddings.parent();
    return field != nullptr && signature(sig, *field);
}

bool compact_places(slong& places, EmbeddingContext& embeddings) noexcept {
    Signature sig;
    if (!signature_from_embeddings(sig, embeddings)) {
        return false;
    }
    places = sig.r1() + sig.r2();
    return true;
}

bool compact_log_matrix(flint::ArbMat& out,
                        EmbeddingContext& embeddings,
                        FactoredElementSpan generators,
                        slong precision) noexcept {
    if (precision <= 0) {
        return false;
    }

    slong places = 0;
    if (!compact_places(places, embeddings) ||
        flint::arb_mat_nrows_value(out) !=
                static_cast<slong>(generators.size()) ||
        flint::arb_mat_ncols_value(out) != places) {
        return false;
    }

    flint::ArbVec row(places);
    for (std::size_t i = 0; i < generators.size(); ++i) {
        if (!generators[i].logarithmic_embedding(
                    flint::ArbVecRef(row), embeddings,
                    LogEmbeddingMode::product, precision)) {
            return false;
        }
        for (slong j = 0; j < places; ++j) {
            arb_set(arb_mat_entry(out.raw(), static_cast<slong>(i), j),
                    row.data() + j);
        }
    }

    return true;
}

bool compact_independence_lower_bound(flint::Arb& out,
                                      slong degree,
                                      slong rank,
                                      slong precision) noexcept {
    if (degree <= 1 || rank <= 0 || precision <= 0) {
        return false;
    }

    flint::Arb log_degree;
    flint::Arb term;
    flint::Arb power;
    flint::Arb inverse_rank_power;
    flint::Arb inverse_rank;
    flint::arb_log_ui(log_degree, static_cast<ulong>(degree), precision);
    flint::arb_mul_ui(term, log_degree, 21, precision);
    flint::arb_div_ui(term, term, 128, precision);
    flint::arb_div_ui(term, term, static_cast<ulong>(degree), precision);
    flint::arb_div_ui(term, term, static_cast<ulong>(degree), precision);

    flint::arb_one(power);
    for (slong i = 0; i < 2 * rank; ++i) {
        flint::arb_mul(power, power, term, precision);
    }

    flint::arb_one(inverse_rank);
    flint::arb_div_ui(inverse_rank, inverse_rank, static_cast<ulong>(rank),
                      precision);
    flint::arb_one(inverse_rank_power);
    for (slong i = 0; i < rank; ++i) {
        flint::arb_mul(inverse_rank_power, inverse_rank_power, inverse_rank,
                       precision);
    }

    flint::arb_mul(out, inverse_rank_power, power, precision);
    return flint::arb_is_finite(out);
}

bool compact_independence_from_log_matrix(bool& decided,
                                          bool& independent,
                                          const flint::ArbMat& logs,
                                          slong len,
                                          slong places,
                                          slong degree,
                                          slong precision) noexcept {
    decided = false;
    independent = false;
    if (len == 0) {
        decided = true;
        independent = true;
        return true;
    }
    if (len >= places || degree <= 1 || precision <= 0) {
        return false;
    }

    flint::ArbMat transpose(places, len);
    flint::ArbMat gram(len, len);
    flint::Arb determinant;
    arb_mat_transpose(transpose.raw(), logs.raw());
    arb_mat_mul(gram.raw(), logs.raw(), transpose.raw(), precision);
    arb_mat_det(determinant.raw(), gram.raw(), precision);

    if (!flint::arb_is_finite(determinant)) {
        return true;
    }

    flint::Arb lower_bound;
    flint::Arb difference;
    if (!compact_independence_lower_bound(lower_bound, degree, places - 1,
                                          precision)) {
        return false;
    }
    flint::arb_sub(difference, lower_bound, determinant, precision);
    if (flint::arb_is_positive(difference)) {
        decided = true;
        independent = false;
        return true;
    }
    if (flint::arb_is_positive(determinant)) {
        decided = true;
        independent = true;
        return true;
    }
    if (flint::arb_is_zero(determinant)) {
        decided = true;
        independent = false;
        return true;
    }
    return true;
}

bool compact_independent(bool& independent,
                         EmbeddingContext& embeddings,
                         FactoredElementSpan generators,
                         slong precision) noexcept {
    if (generators.empty()) {
        independent = true;
        return true;
    }

    Signature sig;
    const slong len = static_cast<slong>(generators.size());
    if (precision <= 0 || !signature_from_embeddings(sig, embeddings) ||
        len >= sig.r1() + sig.r2()) {
        return false;
    }

    const slong places = sig.r1() + sig.r2();
    slong work_precision = precision;
    for (;;) {
        flint::ArbMat logs(len, places);
        bool decided = false;
        if (!compact_log_matrix(logs, embeddings, generators,
                                work_precision) ||
            !compact_independence_from_log_matrix(
                    decided, independent, logs, len, places, sig.degree(),
                    work_precision)) {
            return false;
        }
        if (decided) {
            return true;
        }

        if (work_precision > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        work_precision *= 2;
    }
}

bool compact_regulator_from_log_matrix(flint::ArbRef out,
                                       const flint::ArbMat& logs,
                                       slong len,
                                       slong places,
                                       slong precision) noexcept {
    if (len == 0) {
        arb_set_ui(out.raw(), 1);
        return true;
    }
    if (len + 1 != places || precision <= 0) {
        return false;
    }

    flint::ArbMat minor(len, len);
    flint::Arb determinant;
    for (slong i = 0; i < len; ++i) {
        for (slong j = 0; j < len; ++j) {
            arb_set(arb_mat_entry(minor.raw(), i, j),
                    arb_mat_entry(logs.raw(), i, j));
        }
    }

    arb_mat_det(determinant.raw(), minor.raw(), precision);
    arb_abs(determinant.raw(), determinant.raw());
    if (!flint::arb_is_finite(determinant) ||
        flint::arb_contains_zero(determinant)) {
        return false;
    }

    arb_set(out.raw(), determinant.raw());
    return true;
}

bool arb_radius_lt_2exp(const flint::Arb& value, slong exponent) noexcept {
    return mag_cmp_2exp_si(arb_radref(value.raw()), exponent) < 0;
}

bool class_regulator_index_bound_from_candidate_product(
        flint::FmpzRef out,
        flint::ArbConstRef candidate_class_regulator_product,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    if (precision <= 0 ||
        !flint::arb_is_finite(candidate_class_regulator_product) ||
        !flint::arb_is_positive(candidate_class_regulator_product) ||
        !flint::arb_is_finite(analytic_class_regulator_product) ||
        !flint::arb_is_positive(analytic_class_regulator_product)) {
        return false;
    }

    flint::Arb quotient;
    flint::Arf upper;
    flint::Fmpz candidate;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.class_regulator_index_bound.quotient_bound");
        ::arb_div(quotient.raw(), candidate_class_regulator_product.raw(),
                  analytic_class_regulator_product.raw(), precision);
        if (!flint::arb_is_finite(quotient) ||
            !flint::arb_is_positive(quotient)) {
            return false;
        }

        flint::arb_get_ubound_arf(upper, quotient, precision);
        if (!flint::arf_is_finite(upper)) {
            return false;
        }
    }

    flint::arf_get_fmpz(candidate, upper, ARF_RND_CEIL);
    if (flint::fmpz_sgn(flint::FmpzConstRef(candidate)) <= 0) {
        flint::fmpz_one(flint::FmpzRef(candidate));
    }

    flint::fmpz_set(out, flint::FmpzConstRef(candidate));
    return true;
}

bool class_regulator_index_is_one_from_candidate_product(
        flint::ArbConstRef candidate_class_regulator_product,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision,
        const DiagnosticsContext* diagnostics) noexcept {
    if (precision <= 0 ||
        !flint::arb_is_finite(candidate_class_regulator_product) ||
        !flint::arb_is_positive(candidate_class_regulator_product) ||
        !flint::arb_is_finite(analytic_class_regulator_product) ||
        !flint::arb_is_positive(analytic_class_regulator_product)) {
        return false;
    }

    flint::Arb quotient;
    flint::Arf upper;
    {
        SILEX_PROFILE_SCOPE(
                diagnostics, DiagnosticsModule::unit_group,
                "unit_group.class_regulator_index_is_one.quotient");
        ::arb_div(quotient.raw(), candidate_class_regulator_product.raw(),
                  analytic_class_regulator_product.raw(), precision);
        if (!flint::arb_is_finite(quotient) ||
            !flint::arb_is_positive(quotient) ||
            !flint::arb_contains_si(quotient, 1)) {
            return false;
        }

        flint::arb_get_ubound_arf(upper, quotient, precision);
    }
    return flint::arf_is_finite(upper) &&
           flint::arf_cmp_ui(upper, UWORD(2)) < 0;
}

bool compact_regulator_adaptive(flint::ArbRef out,
                                EmbeddingContext& embeddings,
                                FactoredElementSpan generators,
                                slong start_precision,
                                slong abs_tolerance) noexcept {
    if (start_precision <= 0 || abs_tolerance <= 0) {
        return false;
    }

    const slong len = static_cast<slong>(generators.size());
    if (len == 0) {
        arb_set_ui(out.raw(), 1);
        return true;
    }

    slong places = 0;
    if (!compact_places(places, embeddings) || len + 1 != places) {
        return false;
    }

    slong precision = max_slong(start_precision, 32);
    for (;;) {
        flint::ArbMat logs(len, places);
        flint::Arb regulator;
        if (!compact_log_matrix(logs, embeddings, generators, precision)) {
            return false;
        }
        if (compact_regulator_from_log_matrix(
                    flint::ArbRef(regulator), logs, len, places,
                    precision) &&
            arb_radius_lt_2exp(regulator, -abs_tolerance)) {
            arb_set(out.raw(), regulator.raw());
            return true;
        }

        if (precision > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        precision *= 2;
    }
}

bool regulator_ubound_lt(const flint::Arb& left,
                         const flint::Arb& right,
                         slong precision) noexcept {
    flint::Arf left_upper;
    flint::Arf right_upper;
    flint::arb_get_ubound_arf(left_upper, left, precision);
    flint::arb_get_ubound_arf(right_upper, right, precision);
    return flint::arf_is_finite(left_upper) &&
           flint::arf_is_finite(right_upper) &&
           flint::arf_cmp(left_upper, right_upper) < 0;
}

bool evaluated_is_order_unit(const Order& order,
                             const Element& value) noexcept {
    OrderElement order_value(order);
    Ideal principal(order);
    return order_value.is_defined() && principal.is_defined() &&
           order_value.set_element(value) &&
           principal.set_principal(order_value) &&
           principal.is_one();
}

}  // namespace detail

namespace {

bool copy_and_validate_generators(std::vector<FactoredElement>& out,
                                  const Order& order,
                                  FactoredElementSpan generators) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    out.clear();
    out.reserve(generators.size());
    for (std::size_t i = 0; i < generators.size(); ++i) {
        if (generators[i].parent() == nullptr ||
            !generators[i].parent()->has_same_data(*field)) {
            return false;
        }

        Element expanded(*field);
        if (!generators[i].evaluate(expanded) ||
            !detail::evaluated_is_order_unit(order, expanded)) {
            return false;
        }

        out.emplace_back(*field);
        if (!out.back().set(generators[i])) {
            return false;
        }
    }

    return true;
}

bool copy_trusted_generators(std::vector<FactoredElement>& out,
                             const Order& order,
                             FactoredElementSpan generators) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    out.clear();
    out.reserve(generators.size());
    for (std::size_t i = 0; i < generators.size(); ++i) {
        if (generators[i].parent() == nullptr ||
            !generators[i].parent()->has_same_data(*field)) {
            return false;
        }

        out.emplace_back(*field);
        if (!out.back().set(generators[i])) {
            return false;
        }
    }

    return true;
}

}  // namespace

OrderUnitGroup::OrderUnitGroup(const Order& parent) noexcept {
    define(parent);
}

OrderUnitGroup::~OrderUnitGroup() noexcept {
    clear();
}

OrderUnitGroup::OrderUnitGroup(OrderUnitGroup&& other) noexcept {
    swap(other);
}

OrderUnitGroup& OrderUnitGroup::operator=(OrderUnitGroup&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void OrderUnitGroup::swap(OrderUnitGroup& other) noexcept {
    parent_.swap(other.parent_);
    std::swap(torsion_order_, other.torsion_order_);
    torsion_generator_.swap(other.torsion_generator_);
    free_generators_.swap(other.free_generators_);
    std::swap(free_generator_count_, other.free_generator_count_);
    std::swap(free_generator_capacity_, other.free_generator_capacity_);
    std::swap(regulator_, other.regulator_);
    std::swap(has_regulator_, other.has_regulator_);
    std::swap(certification_, other.certification_);
    unit_proof_records_.swap(other.unit_proof_records_);
    std::swap(unit_proof_record_count_, other.unit_proof_record_count_);
    std::swap(unit_proof_record_capacity_, other.unit_proof_record_capacity_);
    std::swap(is_set_, other.is_set_);
    std::swap(diagnostics_, other.diagnostics_);
}

void OrderUnitGroup::clear() noexcept {
    parent_.clear();
    flint::fmpz_zero(flint::FmpzRef(torsion_order_));
    torsion_generator_.clear();
    clear_free_generators_();
    flint::arb_zero(regulator_);
    has_regulator_ = false;
    certification_ = CertificationMode::unknown;
    clear_unit_proof_records_();
    is_set_ = false;
    diagnostics_ = nullptr;
}

bool OrderUnitGroup::define(const Order& parent) noexcept {
    if (!parent.is_defined() || !parent.has_basis() || parent.parent() == nullptr) {
        return false;
    }

    OrderUnitGroup next;
    next.diagnostics_ = diagnostics_;
    next.parent_ = parent;
    next.torsion_generator_ = OrderElement(parent);
    if (!next.torsion_generator_.is_defined()) {
        return false;
    }

    swap(next);
    return true;
}

bool OrderUnitGroup::reserve_free_generators_(slong capacity) noexcept {
    if (capacity < 0) {
        return false;
    }
    if (capacity <= free_generator_capacity_) {
        return true;
    }

    std::unique_ptr<FactoredElement[]> next(
            new (std::nothrow) FactoredElement[capacity]);
    if (next == nullptr) {
        return false;
    }

    for (slong i = 0; i < free_generator_count_; ++i) {
        next[i] = std::move(free_generators_[i]);
    }

    free_generators_.swap(next);
    free_generator_capacity_ = capacity;
    return true;
}

bool OrderUnitGroup::append_free_generator_copy_(
        const FactoredElement& generator) noexcept {
    const NumberField* field = parent_.parent();
    if (field == nullptr || generator.parent() == nullptr ||
        !generator.parent()->has_same_data(*field)) {
        return false;
    }

    if (free_generator_count_ == free_generator_capacity_) {
        const slong next_capacity =
                free_generator_capacity_ == 0 ? 1 : 2 * free_generator_capacity_;
        if (!reserve_free_generators_(next_capacity)) {
            return false;
        }
    }

    FactoredElement copy(*field);
    if (!copy.is_defined() || !copy.set(generator)) {
        return false;
    }

    free_generators_[free_generator_count_] = std::move(copy);
    ++free_generator_count_;
    return true;
}

void OrderUnitGroup::clear_free_generators_() noexcept {
    free_generators_.reset();
    free_generator_count_ = 0;
    free_generator_capacity_ = 0;
}

bool OrderUnitGroup::reserve_unit_proof_records_(slong capacity) noexcept {
    if (capacity < 0) {
        return false;
    }
    if (capacity <= unit_proof_record_capacity_) {
        return true;
    }

    std::unique_ptr<detail::UnitProofRecordData[]> next(
            new (std::nothrow) detail::UnitProofRecordData[capacity]);
    if (next == nullptr) {
        return false;
    }

    for (slong i = 0; i < unit_proof_record_count_; ++i) {
        flint::fmpz_set(flint::FmpzRef(next[i].ell),
                        flint::FmpzConstRef(unit_proof_records_[i].ell));
        next[i].status = unit_proof_records_[i].status;
        flint::fmpz_set(
                flint::FmpzRef(next[i].aux_prime_bound),
                flint::FmpzConstRef(unit_proof_records_[i].aux_prime_bound));
        next[i].local_primes = unit_proof_records_[i].local_primes;
        next[i].changed = unit_proof_records_[i].changed;
    }

    unit_proof_records_.swap(next);
    unit_proof_record_capacity_ = capacity;
    return true;
}

detail::UnitProofRecordData*
OrderUnitGroup::append_unit_proof_record_() noexcept {
    if (unit_proof_record_count_ == unit_proof_record_capacity_) {
        const slong next_capacity = unit_proof_record_capacity_ == 0
                ? 1
                : 2 * unit_proof_record_capacity_;
        if (!reserve_unit_proof_records_(next_capacity)) {
            return nullptr;
        }
    }

    detail::UnitProofRecordData* record =
            &unit_proof_records_[unit_proof_record_count_];
    ++unit_proof_record_count_;
    return record;
}

bool OrderUnitGroup::copy_unit_proof_records_from_(
        const OrderUnitGroup& other) noexcept {
    clear_unit_proof_records_();
    if (!reserve_unit_proof_records_(other.unit_proof_record_count_)) {
        return false;
    }
    for (slong i = 0; i < other.unit_proof_record_count_; ++i) {
        detail::UnitProofRecordData* dest = append_unit_proof_record_();
        if (dest == nullptr) {
            return false;
        }
        const detail::UnitProofRecordData& src = other.unit_proof_records_[i];
        flint::fmpz_set(flint::FmpzRef(dest->ell),
                        flint::FmpzConstRef(src.ell));
        dest->status = src.status;
        flint::fmpz_set(flint::FmpzRef(dest->aux_prime_bound),
                        flint::FmpzConstRef(src.aux_prime_bound));
        dest->local_primes = src.local_primes;
        dest->changed = src.changed;
    }
    return true;
}

void OrderUnitGroup::clear_unit_proof_records_() noexcept {
    unit_proof_records_.reset();
    unit_proof_record_count_ = 0;
    unit_proof_record_capacity_ = 0;
}

bool OrderUnitGroup::set(const OrderUnitGroup& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    OrderUnitGroup copy(other.parent_);
    if (!copy.is_defined()) {
        return false;
    }

    flint::fmpz_set(flint::FmpzRef(copy.torsion_order_),
                    flint::FmpzConstRef(other.torsion_order_));
    if (!copy.torsion_generator_.set(other.torsion_generator_)) {
        return false;
    }
    if (!copy.reserve_free_generators_(other.free_generator_count_)) {
        return false;
    }
    for (slong i = 0; i < other.free_generator_count_; ++i) {
        if (!copy.append_free_generator_copy_(other.free_generators_[i])) {
            return false;
        }
    }
    flint::arb_set(copy.regulator_, other.regulator_);
    copy.has_regulator_ = other.has_regulator_;
    copy.certification_ = other.certification_;
    if (!copy.copy_unit_proof_records_from_(other)) {
        return false;
    }
    copy.is_set_ = other.is_set_;
    copy.diagnostics_ = other.diagnostics_;

    swap(copy);
    return true;
}

bool OrderUnitGroup::is_defined() const noexcept {
    return parent_.is_defined();
}

const Order* OrderUnitGroup::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

void OrderUnitGroup::set_diagnostics(
        const DiagnosticsContext* diagnostics) noexcept {
    diagnostics_ = diagnostics;
}

const DiagnosticsContext* OrderUnitGroup::diagnostics() const noexcept {
    return diagnostics_;
}

bool OrderUnitGroup::is_set() const noexcept {
    return is_defined() && is_set_;
}

slong OrderUnitGroup::free_rank() const noexcept {
    return is_set() ? free_generator_count_ : -1;
}

CertificationMode OrderUnitGroup::certification_status() const noexcept {
    return is_set() ? certification_ : CertificationMode::unknown;
}

void OrderUnitGroup::mark_certification_proven_() noexcept {
    if (is_set()) {
        certification_ = CertificationMode::proven;
    }
}

bool OrderUnitGroup::torsion_order(flint::FmpzRef out) const noexcept {
    if (!is_set()) {
        return false;
    }
    flint::fmpz_set(out, flint::FmpzConstRef(torsion_order_));
    return true;
}

std::optional<flint::Fmpz> OrderUnitGroup::torsion_order() const noexcept {
    flint::Fmpz out;
    if (!torsion_order(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool OrderUnitGroup::torsion_generator(OrderElement& out) const noexcept {
    if (!is_set() || !same_order_parent(out.parent(), &parent_)) {
        return false;
    }
    return out.set(torsion_generator_);
}

bool OrderUnitGroup::free_generator(FactoredElement& out,
                                    slong index) const noexcept {
    if (!is_set() || index < 0 ||
        index >= free_generator_count_ ||
        parent_.parent() == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*parent_.parent())) {
        return false;
    }
    return out.set(free_generators_[index]);
}

bool OrderUnitGroup::regulator(flint::ArbRef out) const noexcept {
    if (!is_set() || !has_regulator_) {
        return false;
    }
    flint::arb_set(out.raw(), regulator_);
    return true;
}

std::optional<flint::Arb> OrderUnitGroup::regulator() const noexcept {
    flint::Arb out;
    if (!regulator(flint::ArbRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool OrderUnitGroup::class_regulator_product(
        flint::ArbRef out,
        const ClassGroupContext& class_group,
        slong precision) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.class_regulator_product");
    if (!is_set() || !same_order_parent(parent(), class_group.parent()) ||
        !class_group.has_presentation() || precision <= 0) {
        return false;
    }

    flint::Fmpz class_order;
    flint::Arb regulator_value;
    flint::Arb product;
    if (!class_group.order(flint::FmpzRef(class_order)) ||
        !regulator(flint::ArbRef(regulator_value))) {
        return false;
    }

    flint::arb_mul_fmpz(product, regulator_value,
                        flint::FmpzConstRef(class_order), precision);
    flint::arb_set(out.raw(), product);
    return true;
}

bool OrderUnitGroup::class_regulator_index_bound(
        flint::FmpzRef out,
        const ClassGroupContext& class_group,
        flint::ArbConstRef analytic_class_regulator_product,
        slong precision) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.class_regulator_index_bound");
    if (precision <= 0 ||
        !flint::arb_is_finite(analytic_class_regulator_product) ||
        !flint::arb_is_positive(analytic_class_regulator_product)) {
        return false;
    }

    flint::Arb candidate_hR;
    {
        SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                            "unit_group.class_regulator_index_bound.product");
        if (!class_regulator_product(flint::ArbRef(candidate_hR),
                                     class_group, precision)) {
            return false;
        }
    }

    return detail::class_regulator_index_bound_from_candidate_product(
            out, flint::ArbConstRef(candidate_hR),
            analytic_class_regulator_product, precision, diagnostics_);
}

slong OrderUnitGroup::unit_proof_record_count() const noexcept {
    return is_set() ? unit_proof_record_count_ : 0;
}

bool OrderUnitGroup::unit_proof_record(
        flint::FmpzRef ell,
        ProofState& status,
        flint::FmpzRef aux_prime_bound,
        slong& local_primes,
        bool& changed,
        slong index) const noexcept {
    if (!is_set() || index < 0 ||
        index >= unit_proof_record_count_) {
        return false;
    }

    const detail::UnitProofRecordData& record =
            unit_proof_records_[index];
    flint::fmpz_set(ell, flint::FmpzConstRef(record.ell));
    status = record.status;
    flint::fmpz_set(aux_prime_bound,
                    flint::FmpzConstRef(record.aux_prime_bound));
    local_primes = record.local_primes;
    changed = record.changed;
    return true;
}

std::optional<OrderUnitProofRecord> OrderUnitGroup::unit_proof_record(
        slong index) const noexcept {
    OrderUnitProofRecord record;
    if (!unit_proof_record(flint::FmpzRef(record.ell), record.status,
                           flint::FmpzRef(record.aux_prime_bound),
                           record.local_primes, record.changed, index)) {
        return std::nullopt;
    }
    return record;
}

bool OrderUnitGroup::unit_proof_verified(
        flint::FmpzConstRef ell) const noexcept {
    if (!is_set()) {
        return false;
    }
    for (slong i = 0; i < unit_proof_record_count_; ++i) {
        const detail::UnitProofRecordData& record = unit_proof_records_[i];
        if (flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell)) {
            return record.status == ProofState::verified;
        }
    }
    return false;
}

bool OrderUnitGroup::mark_unit_proof(flint::FmpzConstRef ell,
                                     ProofState status,
                                     flint::FmpzConstRef aux_prime_bound,
                                     slong local_primes,
                                     bool changed) noexcept {
    if (!is_set() || !flint::fmpz_is_prime(ell) ||
        (status != ProofState::unavailable &&
         status != ProofState::verified) ||
        fmpz_sgn(aux_prime_bound.raw()) < 0 || local_primes < 0) {
        return false;
    }

    detail::UnitProofRecordData* target = nullptr;
    for (slong i = 0; i < unit_proof_record_count_; ++i) {
        detail::UnitProofRecordData& record = unit_proof_records_[i];
        if (flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell)) {
            target = &record;
            break;
        }
    }

    if (target == nullptr) {
        target = append_unit_proof_record_();
        if (target == nullptr) {
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(target->ell), ell);
    }

    target->status = status;
    flint::fmpz_set(flint::FmpzRef(target->aux_prime_bound),
                    aux_prime_bound);
    target->local_primes = local_primes;
    target->changed = changed;
    return true;
}

bool OrderUnitGroup::mark_unit_proof_after_(
        flint::FmpzConstRef ell,
        ProofState status,
        flint::FmpzConstRef aux_prime_bound,
        slong local_primes,
        bool changed,
        slong stable_prefix_len) noexcept {
    if (!is_set() || !flint::fmpz_is_prime(ell) ||
        (status != ProofState::unavailable &&
         status != ProofState::verified) ||
        fmpz_sgn(aux_prime_bound.raw()) < 0 || local_primes < 0 ||
        stable_prefix_len < 0 || stable_prefix_len > unit_proof_record_count_) {
        return false;
    }

    detail::UnitProofRecordData* target = nullptr;
    for (slong i = 0; i < stable_prefix_len; ++i) {
        detail::UnitProofRecordData& record = unit_proof_records_[i];
        if (flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell)) {
            target = &record;
            break;
        }
    }

    if (target == nullptr && unit_proof_record_count_ > stable_prefix_len) {
        detail::UnitProofRecordData& last =
                unit_proof_records_[unit_proof_record_count_ - 1];
        if (flint::fmpz_equal(flint::FmpzConstRef(last.ell), ell)) {
            target = &last;
        } else if (fmpz_cmp(last.ell.raw(), ell.raw()) > 0) {
            // The caller promises that records appended after
            // `stable_prefix_len` are monotone for this proof scan.  If the
            // suffix does not look monotone for the next ell, fall back to the
            // general update path rather than risking a duplicate.
            for (slong i = stable_prefix_len; i < unit_proof_record_count_;
                 ++i) {
                detail::UnitProofRecordData& record = unit_proof_records_[i];
                if (flint::fmpz_equal(flint::FmpzConstRef(record.ell), ell)) {
                    target = &record;
                    break;
                }
            }
        }
    }

    if (target == nullptr) {
        target = append_unit_proof_record_();
        if (target == nullptr) {
            return false;
        }
        flint::fmpz_set(flint::FmpzRef(target->ell), ell);
    }

    target->status = status;
    flint::fmpz_set(flint::FmpzRef(target->aux_prime_bound),
                    aux_prime_bound);
    target->local_primes = local_primes;
    target->changed = changed;
    return true;
}

void OrderUnitGroup::reset_unit_proof_records() noexcept {
    clear_unit_proof_records_();
}

void OrderUnitGroup::try_certify_index_one(slong precision) noexcept {
    if (!is_set() || precision <= 0) {
        return;
    }

    flint::Fmpz index_bound;
    if (regulator_index_bound(flint::FmpzRef(index_bound), precision) &&
        flint::fmpz_is_one(flint::FmpzConstRef(index_bound))) {
        certification_ = CertificationMode::proven;
    }
}

bool OrderUnitGroup::regulator_index_bound(flint::FmpzRef out,
                                           slong precision) const noexcept {
    if (!is_set() || precision <= 0) {
        return false;
    }

    if (free_rank() == 0) {
        flint::fmpz_one(out);
        return true;
    }

    if (!has_regulator_ || !flint::arb_is_positive(regulator_)) {
        return false;
    }

    const Order* order = parent();
    const NumberField* field = order == nullptr ? nullptr : order->parent();
    if (field == nullptr) {
        return false;
    }

    flint::Arb lower;
    flint::Arb quotient;
    flint::Arf upper;
    flint::Fmpz candidate;
    if (!unit_lower_regulator_bound(flint::ArbRef(lower), *field,
                                    precision) ||
        !flint::arb_is_positive(lower)) {
        return false;
    }

    flint::arb_div(quotient, regulator_, lower, precision);
    if (!flint::arb_is_finite(quotient) ||
        !flint::arb_is_positive(quotient)) {
        return false;
    }

    flint::arb_get_ubound_arf(upper, quotient, precision);
    if (!flint::arf_is_finite(upper)) {
        return false;
    }

    flint::arf_get_fmpz(candidate, upper, ARF_RND_CEIL);
    if (flint::fmpz_sgn(flint::FmpzConstRef(candidate)) <= 0) {
        flint::fmpz_one(flint::FmpzRef(candidate));
    }

    flint::fmpz_set(out, flint::FmpzConstRef(candidate));
    return true;
}

namespace detail {

bool order_unit_group_set_units_internal(
        OrderUnitGroup& out,
        const Order& order,
        FactoredElementSpan generators,
        EmbeddingContext& embeddings,
        slong precision,
        bool trusted,
        const flint::Fmpz* cached_torsion_order,
        const OrderElement* cached_torsion_generator) noexcept {
    SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                        "unit_group.set_units_internal");
    if (!order.is_defined() || order.parent() == nullptr || !order.has_basis() ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*order.parent()) ||
        precision <= 0 ||
        (generators.size() > 0 && generators.data() == nullptr) ||
        generators.size() >
                static_cast<std::size_t>(std::numeric_limits<slong>::max())) {
        return false;
    }

    const slong len = static_cast<slong>(generators.size());
    slong rank = -1;
    if (!unit_rank(rank, *order.parent()) || len != rank) {
        return false;
    }

    OrderUnitGroup candidate(order);
    candidate.set_diagnostics(out.diagnostics());
    if (!candidate.is_defined()) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "set_units failed while defining candidate");
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.set_units_torsion");
        if (cached_torsion_order != nullptr &&
            cached_torsion_generator != nullptr) {
            const Order* torsion_parent = cached_torsion_generator->parent();
            if (torsion_parent == nullptr ||
                !torsion_parent->has_same_data(order) ||
                flint::fmpz_sgn(flint::FmpzConstRef(
                        *cached_torsion_order)) <= 0 ||
                !candidate.torsion_generator_.set(
                        *cached_torsion_generator)) {
                SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "set_units failed while copying cached torsion");
                return false;
            }
            flint::fmpz_set(flint::FmpzRef(candidate.torsion_order_),
                            flint::FmpzConstRef(*cached_torsion_order));
        } else if (cached_torsion_order != nullptr ||
                   cached_torsion_generator != nullptr ||
                   !rank_zero_torsion(
                           flint::FmpzRef(candidate.torsion_order_),
                           candidate.torsion_generator_, order)) {
            SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "set_units failed while setting torsion");
            return false;
        }
    }
    std::vector<FactoredElement> copied_generators;
    {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            trusted ? "unit_group.set_units_copy_trusted"
                                    : "unit_group.set_units_copy_validate");
        if (!(trusted ? copy_trusted_generators(copied_generators,
                                                order, generators)
                      : copy_and_validate_generators(
                                copied_generators, order,
                                generators))) {
            SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "set_units failed while copying generators");
            return false;
        }
    }
    if (!candidate.reserve_free_generators_(len)) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "set_units failed while reserving generators");
        return false;
    }
    {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.set_units_append");
        for (const FactoredElement& generator : copied_generators) {
            if (!candidate.append_free_generator_copy_(generator)) {
                SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                          LogLevel::detail,
                          "set_units failed while appending generators");
                return false;
            }
        }
    }

    slong places = 0;
    if (!compact_places(places, embeddings) || len + 1 != places) {
        SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "set_units failed while checking compact places");
        return false;
    }

    bool independent = false;
    {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.set_units_log_independence");
        if (!compact_independent(
                    independent, embeddings,
                    FactoredElementSpan(copied_generators.data(),
                                        copied_generators.size()),
                    precision) ||
            !independent) {
            SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "set_units failed while checking log independence");
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(out.diagnostics(), DiagnosticsModule::unit_group,
                            "unit_group.set_units_regulator");
        if (!compact_regulator_adaptive(
                    flint::ArbRef(candidate.regulator_), embeddings,
                    FactoredElementSpan(copied_generators.data(),
                                        copied_generators.size()),
                    precision, kRegulatorAbsTolerance)) {
            SILEX_LOG(out.diagnostics(), DiagnosticsModule::unit_group,
                      LogLevel::detail,
                      "set_units failed while computing regulator");
            return false;
        }
    }

    candidate.has_regulator_ = true;
    candidate.certification_ = CertificationMode::unknown;
    candidate.is_set_ = true;

    out.swap(candidate);
    return true;
}

}  // namespace detail

bool OrderUnitGroup::set_units(const Order& order,
                               FactoredElementSpan generators,
                               EmbeddingContext& embeddings,
                               slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.set_units");
    return detail::order_unit_group_set_units_internal(
            *this, order, generators, embeddings, precision, false);
}

bool OrderUnitGroup::set_relation_kernel_units(
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.set_relation_kernel_units");
    OrderUnitGroup candidate(order);
    candidate.set_diagnostics(diagnostics_);
    if (!candidate.is_defined() ||
        !detail::set_initial_relation_kernel_units(
                candidate, order, class_group, embeddings, precision)) {
        return false;
    }

    swap(candidate);
    return true;
}

bool OrderUnitGroup::set_relation_kernel_units_bounded(
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        flint::FmpzConstRef denominator_bound,
        slong start_precision,
        slong max_precision) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics_, DiagnosticsModule::unit_group,
                        "unit_group.set_relation_kernel_units_bounded");
    if (!detail::validate_relation_kernel_inputs(order, class_group,
                                                 embeddings, max_precision) ||
        flint::fmpz_sgn(denominator_bound) <= 0 ||
        !flint::fmpz_fits_si(denominator_bound) ||
        start_precision <= 0 || max_precision < start_precision) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }
    if (rank == 0) {
        return compute(order);
    }

    const slong count = class_group.relation_kernel_unit_count();
    if (count < rank) {
        return false;
    }

    OrderUnitGroup current(order);
    current.set_diagnostics(diagnostics_);
    if (!current.is_defined() ||
        !detail::set_initial_relation_kernel_units(
                current, order, class_group, embeddings, max_precision)) {
        return false;
    }

    const NumberField* field = order.parent();
    FactoredElement y(*field);
    FactoredElement root(*field);
    flint::FmpzMat rel(1, rank + 1);
    flint::Fmpz torsion_exp;
    OrderUnitGroup refined(order);
    refined.set_diagnostics(diagnostics_);

    for (slong pass = 0; pass < count; ++pass) {
        bool changed = false;
        for (slong i = 0; i < count; ++i) {
            if (!detail::relation_kernel_generator(y, order, class_group, i)) {
                continue;
            }

            bool recovered = false;
            if (!detail::dependent_relation_bounded_min_denominator(
                        recovered, root, rel, torsion_exp, current, y,
                        embeddings, denominator_bound, 2, start_precision,
                        max_precision, false, false)) {
                swap(current);
                return true;
            }

            if (!recovered) {
                continue;
            }

            if (!detail::adjoin_verified_dependent_relation(
                        changed, refined, current, root,
                        flint::FmpzMatConstRef(rel), 0, embeddings,
                        max_precision)) {
                swap(current);
                return true;
            }
            if (changed) {
                current.swap(refined);
                break;
            }
        }

        if (!changed) {
            break;
        }
    }

    swap(current);
    return true;
}

bool OrderUnitGroup::set_relation_kernel_units_index_bounded(
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::unit_group,
            "unit_group.set_relation_kernel_units_index_bounded");
    if (!detail::validate_relation_kernel_inputs(order, class_group,
                                                 embeddings, max_precision) ||
        start_precision <= 0 || max_precision < start_precision) {
        return false;
    }

    slong rank = -1;
    if (!unit_rank(rank, *order.parent())) {
        return false;
    }
    if (rank == 0) {
        return compute(order);
    }

    OrderUnitGroup initial(order);
    OrderUnitGroup refined(order);
    initial.set_diagnostics(diagnostics_);
    refined.set_diagnostics(diagnostics_);
    if (!initial.is_defined() || !refined.is_defined() ||
        !detail::set_initial_relation_kernel_units(
                initial, order, class_group, embeddings, max_precision)) {
        return false;
    }

    flint::Fmpz denominator_bound;
    if (initial.regulator_index_bound(flint::FmpzRef(denominator_bound),
                                      max_precision) &&
        flint::fmpz_sgn(flint::FmpzConstRef(denominator_bound)) > 0 &&
        flint::fmpz_fits_si(flint::FmpzConstRef(denominator_bound)) &&
        refined.set_relation_kernel_units_bounded(
                order, class_group, embeddings,
                flint::FmpzConstRef(denominator_bound), start_precision,
                max_precision)) {
        swap(refined);
    } else {
        swap(initial);
    }

    try_certify_index_one(max_precision);
    return true;
}

bool OrderUnitGroup::set_relation_kernel_units_index_bounded_saturated(
        bool& changed,
        bool& stable,
        const Order& order,
        const ClassGroupContext& class_group,
        EmbeddingContext& embeddings,
        slong start_precision,
        slong max_precision,
        slong aux_target_len,
        flint::FmpzConstRef aux_bound_start,
        flint::FmpzConstRef aux_bound_max,
        slong max_passes) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics_, DiagnosticsModule::unit_group,
            "unit_group.set_relation_kernel_units_index_bounded_saturated");
    if (!detail::validate_relation_kernel_inputs(order, class_group,
                                                 embeddings, max_precision) ||
        start_precision <= 0 || max_precision < start_precision ||
        aux_target_len <= 0 || fmpz_cmp_ui(aux_bound_start.raw(), 2) < 0 ||
        fmpz_cmp(aux_bound_max.raw(), aux_bound_start.raw()) < 0 ||
        max_passes <= 0) {
        return false;
    }

    OrderUnitGroup initial(order);
    OrderUnitGroup saturated(order);
    initial.set_diagnostics(diagnostics_);
    saturated.set_diagnostics(diagnostics_);
    if (!initial.is_defined() || !saturated.is_defined() ||
        !initial.set_relation_kernel_units_index_bounded(
                order, class_group, embeddings, start_precision,
                max_precision)) {
        return false;
    }

    bool sat_changed = false;
    bool sat_stable = false;
    if (saturated.saturate_index_bounded_adaptive(
                sat_changed, sat_stable, initial, embeddings, aux_target_len,
                aux_bound_start, aux_bound_max, max_passes, max_precision)) {
        swap(saturated);
        changed = sat_changed;
        stable = sat_stable;
    } else {
        swap(initial);
        changed = false;
        stable = false;
    }

    return true;
}

}  // namespace silex
