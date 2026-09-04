#include <silex/sunit.hpp>

#include "sunit_internal.hpp"
#include "sunit_storage_internal.hpp"

#include "../order_unit/order_unit_internal.hpp"

#include <new>
#include <utility>
#include <vector>

#include <flint/fmpz_mat.h>

#include <silex/unit.hpp>

namespace silex {
namespace {

void compute_fail(SUnitComputeResult& result,
                  SUnitComputeStage stage,
                  slong selected_index = -1) noexcept {
    result.success = false;
    result.stage = stage;
    result.selected_index = selected_index;
}

void membership_fail(SUnitMembershipResult& result,
                     SUnitMembershipStage stage,
                     slong work_precision = 0) noexcept {
    result.success = false;
    result.outcome = SUnitMembershipOutcome::unknown;
    result.stage = stage;
    result.work_precision = work_precision;
}

void membership_reject(SUnitMembershipResult& result,
                       SUnitMembershipStage stage) noexcept {
    result.success = true;
    result.outcome = SUnitMembershipOutcome::not_sunit;
    result.stage = stage;
    result.work_precision = 0;
}

bool copy_selected_primes(std::vector<PrimeIdeal>& out,
                          const Order& order,
                          PrimeIdealSpan input,
                          slong& failed_index) noexcept {
    failed_index = -1;
    if (input.size() > static_cast<std::size_t>(WORD_MAX) ||
        (input.size() != 0 && input.data() == nullptr)) {
        return false;
    }
    out.clear();
    out.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (!input[i].is_defined() || input[i].parent() == nullptr ||
            !input[i].parent()->has_same_data(order)) {
            failed_index = static_cast<slong>(i);
            return false;
        }
        out.emplace_back(order);
        if (!out.back().is_defined() || !out.back().set(input[i])) {
            failed_index = static_cast<slong>(i);
            return false;
        }
    }
    return true;
}

bool copy_selected_primes(std::vector<PrimeIdeal>& out,
                          const Order& order,
                          const std::vector<PrimeIdeal>& input) noexcept {
    out.clear();
    out.reserve(input.size());
    for (const PrimeIdeal& prime : input) {
        out.emplace_back(order);
        if (!out.back().is_defined() || !out.back().set(prime)) {
            return false;
        }
    }
    return true;
}

bool publish_s_class_storage(
        detail::SClassGroupStorage& out,
        const detail::SUnitClassContext& context,
        const SUnitComputeOptions& options) noexcept {
    (void) options;
    out.order = context.order;
    if (!copy_selected_primes(out.selected_primes, out.order,
                              context.selected_primes) ||
        !out.group.set(context.s_class_group)) {
        return false;
    }

    out.invariant_ideals.reserve(context.s_class_invariant_ideals.size());
    for (const FractionalIdeal& source : context.s_class_invariant_ideals) {
        out.invariant_ideals.emplace_back(out.order);
        if (!out.invariant_ideals.back().is_defined() ||
            !out.invariant_ideals.back().set(source)) {
            return false;
        }
    }

    const NumberField* field = out.order.parent();
    if (field == nullptr) {
        return false;
    }
    out.power_witnesses.reserve(context.s_class_power_witnesses.size());
    for (const FactoredElement& source : context.s_class_power_witnesses) {
        out.power_witnesses.emplace_back(*field);
        if (!out.power_witnesses.back().is_defined() ||
            !out.power_witnesses.back().set(source)) {
            return false;
        }
    }

    out.power_selected_exponents = flint::FmpzMat(
            flint::fmpz_mat_nrows(context.s_class_power_selected_exponents),
            flint::fmpz_mat_ncols(context.s_class_power_selected_exponents));
    flint::fmpz_mat_set(
            flint::FmpzMatRef(out.power_selected_exponents),
            flint::FmpzMatConstRef(
                    context.s_class_power_selected_exponents));
    out.source_class_certification = context.source_class_certification;
    out.proof_status = context.s_class_proof_status;
    return true;
}

bool compute_s_regulator(flint::Arb& out,
                         const detail::SUnitClassContext& context,
                         const OrderUnitGroup& ordinary_units,
                         slong precision) noexcept {
    if (precision <= 0 || !ordinary_units.regulator(flint::ArbRef(out)) ||
        !flint::arb_is_finite(out) || !flint::arb_is_positive(out)) {
        return false;
    }

    // reference bnfsunit preserves the ordinary regulator on its empty-S branch.
    if (context.selected_primes.empty()) {
        return true;
    }

    flint::Fmpz s_class_order;
    if (!context.s_class_group.order(flint::FmpzRef(s_class_order)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(s_class_order)) <= 0) {
        return false;
    }
    flint::arb_mul_fmpz(out, out, flint::FmpzConstRef(s_class_order),
                        precision);

    for (const PrimeIdeal& prime : context.selected_primes) {
        flint::Fmpz norm;
        flint::Arb log_norm;
        if (!prime.norm(flint::FmpzRef(norm)) ||
            flint::fmpz_cmp_ui(flint::FmpzConstRef(norm), 2) < 0) {
            return false;
        }
        flint::arb_log_fmpz(log_norm, flint::FmpzConstRef(norm), precision);
        if (!flint::arb_is_finite(log_norm) ||
            !flint::arb_is_positive(log_norm)) {
            return false;
        }
        flint::arb_mul(out, out, log_norm, precision);
    }
    return flint::arb_is_finite(out) && flint::arb_is_positive(out);
}

bool publish_s_unit_storage(
        detail::SUnitGroupStorage& out,
        SUnitComputeStage& failure_stage,
        const detail::SUnitClassContext& context,
        const ClassGroupContext& class_group,
        const OrderUnitGroup& ordinary_units,
        const SUnitComputeOptions& options) noexcept {
    failure_stage = SUnitComputeStage::s_unit_publication;
    if (!copy_selected_primes(out.selected_primes, out.order,
                              context.selected_primes) ||
        !out.ordinary_units.set(ordinary_units)) {
        return false;
    }
    out.diagnostics = options.diagnostics != nullptr
            ? options.diagnostics
            : class_group.diagnostics();
    out.ordinary_units.set_diagnostics(out.diagnostics);

    const NumberField* field = out.order.parent();
    if (field == nullptr) {
        return false;
    }
    out.nonunit_generators.reserve(context.generators_mod_units.size());
    for (const FactoredElement& source : context.generators_mod_units) {
        out.nonunit_generators.emplace_back(*field);
        if (!out.nonunit_generators.back().is_defined() ||
            !out.nonunit_generators.back().set(source)) {
            return false;
        }
    }

    out.nonunit_valuations = flint::FmpzMat(
            flint::fmpz_mat_nrows(context.valuation_rows),
            flint::fmpz_mat_ncols(context.valuation_rows));
    flint::fmpz_mat_set(flint::FmpzMatRef(out.nonunit_valuations),
                        flint::FmpzMatConstRef(context.valuation_rows));
    failure_stage = SUnitComputeStage::s_regulator;
    if (!compute_s_regulator(out.regulator, context, ordinary_units,
                             options.regulator_precision)) {
        return false;
    }

    out.regulator_precision = options.regulator_precision;
    out.source_class_certification = class_group.certification_status();
    out.source_unit_certification = ordinary_units.certification_status();
    out.source_relation_saturation =
            class_group.relation_saturation_status();
    out.source_unit_proof = class_group.unit_proof_status();
    out.source_regulator_proof = class_group.regulator_proof_status();
    out.proof_status = context.s_unit_mod_units_proof_status;
    out.regulator_proof_status = ProofState::verified;
    return true;
}

bool coordinate_shapes_match(const detail::SUnitGroupStorage& storage,
                             const SUnitCoordinates& coordinates) noexcept {
    const slong ordinary_rank = storage.ordinary_units.free_rank();
    const slong nonunit_rank =
            static_cast<slong>(storage.nonunit_generators.size());
    return coordinates.defined && ordinary_rank >= 0 &&
           flint::fmpz_mat_nrows(coordinates.ordinary_free_exponents) == 1 &&
           flint::fmpz_mat_ncols(coordinates.ordinary_free_exponents) ==
                   ordinary_rank &&
           flint::fmpz_mat_nrows(coordinates.nonunit_exponents) == 1 &&
           flint::fmpz_mat_ncols(coordinates.nonunit_exponents) ==
                   nonunit_rank;
}

bool compose_sunit_image(FactoredElement& out,
                         const detail::SUnitGroupStorage& storage,
                         const SUnitCoordinates& coordinates) noexcept {
    const NumberField* field = storage.order.parent();
    if (!storage.defined || field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*field) ||
        !coordinate_shapes_match(storage, coordinates)) {
        return false;
    }

    FactoredElement product(*field);
    OrderElement torsion(storage.order);
    Element torsion_value(*field);
    FactoredElement torsion_factor(*field);
    if (!product.one() || !storage.ordinary_units.torsion_generator(torsion) ||
        !torsion.get_element(torsion_value) ||
        !torsion_factor.set_element(torsion_value) ||
        !detail::multiply_factored_element_power_fmpz(
                product, torsion_factor,
                flint::FmpzConstRef(coordinates.torsion_exponent))) {
        return false;
    }

    const slong ordinary_rank = storage.ordinary_units.free_rank();
    for (slong i = 0; i < ordinary_rank; ++i) {
        FactoredElement generator(*field);
        if (!storage.ordinary_units.free_generator(generator, i) ||
            !detail::multiply_factored_element_power_fmpz(
                    product, generator,
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(
                                    coordinates.ordinary_free_exponents),
                            0, i))) {
            return false;
        }
    }
    for (slong i = 0;
         i < static_cast<slong>(storage.nonunit_generators.size()); ++i) {
        if (!detail::multiply_factored_element_power_fmpz(
                    product,
                    storage.nonunit_generators[static_cast<std::size_t>(i)],
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(
                                    coordinates.nonunit_exponents),
                            0, i))) {
            return false;
        }
    }
    product.normalize();
    out.swap(product);
    return true;
}

enum class ValuationSolveStatus {
    failure,
    not_integral,
    success,
};

ValuationSolveStatus solve_valuation_coordinates(
        flint::FmpzMat& out,
        flint::FmpzMatConstRef valuations,
        flint::FmpzMatConstRef target) noexcept {
    const slong rank = flint::fmpz_mat_nrows(valuations);
    if (flint::fmpz_mat_ncols(valuations) != rank ||
        flint::fmpz_mat_nrows(target) != 1 ||
        flint::fmpz_mat_ncols(target) != rank) {
        return ValuationSolveStatus::failure;
    }
    out = flint::FmpzMat(1, rank);
    if (rank == 0) {
        return ValuationSolveStatus::success;
    }

    flint::FmpzMat left(rank, rank);
    flint::FmpzMat right(rank, 1);
    flint::FmpzMat solution(rank, 1);
    flint::Fmpz denominator;
    flint::fmpz_mat_transpose(flint::FmpzMatRef(left), valuations);
    flint::fmpz_mat_transpose(flint::FmpzMatRef(right), target);
    if (::fmpz_mat_can_solve(solution.raw(), denominator.raw(), left.raw(),
                             right.raw()) == 0 ||
        flint::fmpz_is_zero(flint::FmpzConstRef(denominator))) {
        return ValuationSolveStatus::failure;
    }

    for (slong i = 0; i < rank; ++i) {
        if (::fmpz_divisible(
                    flint::fmpz_mat_entry(
                            flint::FmpzMatConstRef(solution), i, 0).raw(),
                    denominator.raw()) == 0) {
            return ValuationSolveStatus::not_integral;
        }
        ::fmpz_divexact(
                flint::fmpz_mat_entry(flint::FmpzMatRef(out), 0, i).raw(),
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(solution), i, 0).raw(),
                denominator.raw());
    }

    flint::FmpzMat product(1, rank);
    flint::fmpz_mat_mul(flint::FmpzMatRef(product),
                        flint::FmpzMatConstRef(out), valuations);
    return flint::fmpz_mat_equal(flint::FmpzMatConstRef(product), target)
            ? ValuationSolveStatus::success
            : ValuationSolveStatus::failure;
}

bool preimage_core(SUnitMembershipResult& result,
                   SUnitCoordinates& out,
                   const detail::SUnitGroupStorage& storage,
                   const FactoredElement& value,
                   EmbeddingContext& embeddings,
                   slong start_precision,
                   slong max_precision) noexcept {
    result = {};
    const NumberField* field = storage.order.parent();
    if (!storage.defined || field == nullptr || value.parent() == nullptr ||
        !value.parent()->has_same_data(*field) ||
        embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(*field) || start_precision <= 0 ||
        max_precision < start_precision) {
        membership_fail(result, SUnitMembershipStage::input_validation);
        return false;
    }

    Element expanded(*field);
    if (!value.evaluate(expanded)) {
        membership_fail(result, SUnitMembershipStage::selected_valuation);
        return false;
    }
    if (expanded.equal_si(0)) {
        membership_reject(result, SUnitMembershipStage::selected_valuation);
        return true;
    }

    const slong selected_count =
            static_cast<slong>(storage.selected_primes.size());
    flint::FmpzMat target_valuations(1, selected_count);
    for (slong i = 0; i < selected_count; ++i) {
        slong valuation = 0;
        if (!storage.selected_primes[static_cast<std::size_t>(i)].valuation(
                    valuation, value, storage.diagnostics)) {
            membership_fail(result,
                            SUnitMembershipStage::selected_valuation);
            return false;
        }
        flint::fmpz_set_si(
                flint::fmpz_mat_entry(
                        flint::FmpzMatRef(target_valuations), 0, i),
                valuation);
    }

    flint::FmpzMat nonunit_exponents(1, 0);
    const ValuationSolveStatus solve_status = solve_valuation_coordinates(
            nonunit_exponents,
            flint::FmpzMatConstRef(storage.nonunit_valuations),
            flint::FmpzMatConstRef(target_valuations));
    if (solve_status == ValuationSolveStatus::failure) {
        membership_fail(result, SUnitMembershipStage::valuation_solve);
        return false;
    }
    if (solve_status == ValuationSolveStatus::not_integral) {
        membership_reject(result, SUnitMembershipStage::valuation_solve);
        return true;
    }

    SUnitCoordinates nonunit_coordinates;
    nonunit_coordinates.ordinary_free_exponents = flint::FmpzMat(
            1, storage.ordinary_units.free_rank());
    nonunit_coordinates.nonunit_exponents =
            std::move(nonunit_exponents);
    nonunit_coordinates.defined = true;
    FactoredElement nonunit_part(*field);
    FactoredElement residual(*field);
    if (!compose_sunit_image(nonunit_part, storage, nonunit_coordinates) ||
        !residual.divide(value, nonunit_part)) {
        membership_fail(result, SUnitMembershipStage::residual_unit);
        return false;
    }

    detail::OrdinaryUnitCoordinateResult unit_result;
    detail::OrdinaryUnitCoordinates unit_coordinates;
    if (!detail::ordinary_unit_coordinates(
                unit_result, unit_coordinates, storage.ordinary_units,
                residual, embeddings, start_precision, max_precision)) {
        const SUnitMembershipStage stage =
                unit_result.stage ==
                                detail::OrdinaryUnitCoordinateStage::
                                        precision_exhausted
                ? SUnitMembershipStage::precision_exhausted
                : SUnitMembershipStage::residual_unit;
        membership_fail(result, stage, unit_result.work_precision);
        return false;
    }
    if (unit_result.outcome ==
        detail::OrdinaryUnitCoordinateOutcome::not_unit) {
        membership_reject(result, SUnitMembershipStage::residual_unit);
        return true;
    }
    if (unit_result.outcome !=
                detail::OrdinaryUnitCoordinateOutcome::verified ||
        !unit_coordinates.defined) {
        membership_fail(result, SUnitMembershipStage::residual_unit,
                        unit_result.work_precision);
        return false;
    }

    SUnitCoordinates candidate;
    flint::fmpz_set(flint::FmpzRef(candidate.torsion_exponent),
                    flint::FmpzConstRef(
                            unit_coordinates.torsion_exponent));
    candidate.ordinary_free_exponents = flint::FmpzMat(
            1, storage.ordinary_units.free_rank());
    flint::fmpz_mat_set(
            flint::FmpzMatRef(candidate.ordinary_free_exponents),
            flint::FmpzMatConstRef(unit_coordinates.free_exponents));
    candidate.nonunit_exponents = flint::FmpzMat(1, selected_count);
    flint::fmpz_mat_set(
            flint::FmpzMatRef(candidate.nonunit_exponents),
            flint::FmpzMatConstRef(
                    nonunit_coordinates.nonunit_exponents));
    candidate.defined = true;

    FactoredElement round_trip(*field);
    Element round_trip_value(*field);
    if (!compose_sunit_image(round_trip, storage, candidate) ||
        !round_trip.evaluate(round_trip_value) ||
        !round_trip_value.equal(expanded)) {
        membership_fail(result, SUnitMembershipStage::exact_verification,
                        unit_result.work_precision);
        return false;
    }

    out = std::move(candidate);
    result.success = true;
    result.outcome = SUnitMembershipOutcome::verified;
    result.stage = SUnitMembershipStage::none;
    result.work_precision = unit_result.work_precision;
    return true;
}

}  // namespace

SClassGroup::SClassGroup() noexcept = default;

SClassGroup::~SClassGroup() noexcept = default;

SClassGroup::SClassGroup(SClassGroup&& other) noexcept
    : storage_(std::move(other.storage_)) {
}

SClassGroup& SClassGroup::operator=(SClassGroup&& other) noexcept {
    if (this != &other) {
        storage_ = std::move(other.storage_);
    }
    return *this;
}

void SClassGroup::swap(SClassGroup& other) noexcept {
    storage_.swap(other.storage_);
}

void SClassGroup::clear() noexcept {
    storage_.reset();
}

bool SClassGroup::is_defined() const noexcept {
    return storage_ != nullptr && storage_->defined;
}

const Order* SClassGroup::parent() const noexcept {
    return is_defined() ? &storage_->order : nullptr;
}

slong SClassGroup::selected_prime_count() const noexcept {
    return is_defined()
            ? static_cast<slong>(storage_->selected_primes.size())
            : -1;
}

bool SClassGroup::selected_prime(PrimeIdeal& out, slong index) const noexcept {
    return is_defined() && index >= 0 && index < selected_prime_count() &&
           out.set(storage_->selected_primes[static_cast<std::size_t>(index)]);
}

slong SClassGroup::invariant_count() const noexcept {
    return is_defined() ? storage_->group.invariant_count() : -1;
}

bool SClassGroup::invariant(flint::FmpzRef out, slong index) const noexcept {
    return is_defined() && storage_->group.invariant(out, index);
}

std::optional<flint::Fmpz> SClassGroup::invariant(
        slong index) const noexcept {
    flint::Fmpz out;
    if (!invariant(flint::FmpzRef(out), index)) {
        return std::nullopt;
    }
    return out;
}

bool SClassGroup::order(flint::FmpzRef out) const noexcept {
    return is_defined() && storage_->group.order(out);
}

std::optional<flint::Fmpz> SClassGroup::order() const noexcept {
    flint::Fmpz out;
    if (!order(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool SClassGroup::invariant_generator(FractionalIdeal& out,
                                      slong index) const noexcept {
    return is_defined() && index >= 0 && index < invariant_count() &&
           out.set(storage_->invariant_ideals[
                   static_cast<std::size_t>(index)]);
}

bool SClassGroup::invariant_generator_power_witness(
        FactoredElement& out,
        slong index) const noexcept {
    return is_defined() && index >= 0 && index < invariant_count() &&
           out.set(storage_->power_witnesses[
                   static_cast<std::size_t>(index)]);
}

bool SClassGroup::invariant_generator_power_selected_exponents(
        flint::FmpzMatRef out,
        slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= invariant_count() ||
        flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != selected_prime_count()) {
        return false;
    }
    for (slong j = 0; j < selected_prime_count(); ++j) {
        flint::fmpz_set(
                flint::fmpz_mat_entry(out, 0, j),
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(
                                storage_->power_selected_exponents),
                        index, j));
    }
    return true;
}

CertificationMode SClassGroup::certification_status() const noexcept {
    return is_defined() && storage_->proof_status == ProofState::verified
            ? CertificationMode::proven
            : CertificationMode::unknown;
}

CertificationMode SClassGroup::source_class_certification() const noexcept {
    return is_defined() ? storage_->source_class_certification
                        : CertificationMode::unknown;
}

ProofState SClassGroup::proof_status() const noexcept {
    return is_defined() ? storage_->proof_status : ProofState::not_checked;
}

SUnitGroup::SUnitGroup() noexcept = default;

SUnitGroup::~SUnitGroup() noexcept = default;

SUnitGroup::SUnitGroup(SUnitGroup&& other) noexcept
    : storage_(std::move(other.storage_)) {
}

SUnitGroup& SUnitGroup::operator=(SUnitGroup&& other) noexcept {
    if (this != &other) {
        storage_ = std::move(other.storage_);
    }
    return *this;
}

void SUnitGroup::swap(SUnitGroup& other) noexcept {
    storage_.swap(other.storage_);
}

void SUnitGroup::clear() noexcept {
    storage_.reset();
}

bool SUnitGroup::is_defined() const noexcept {
    return storage_ != nullptr && storage_->defined;
}

const Order* SUnitGroup::parent() const noexcept {
    return is_defined() ? &storage_->order : nullptr;
}

slong SUnitGroup::selected_prime_count() const noexcept {
    return is_defined()
            ? static_cast<slong>(storage_->selected_primes.size())
            : -1;
}

bool SUnitGroup::selected_prime(PrimeIdeal& out, slong index) const noexcept {
    return is_defined() && index >= 0 && index < selected_prime_count() &&
           out.set(storage_->selected_primes[static_cast<std::size_t>(index)]);
}

bool SUnitGroup::torsion_order(flint::FmpzRef out) const noexcept {
    return is_defined() && storage_->ordinary_units.torsion_order(out);
}

std::optional<flint::Fmpz> SUnitGroup::torsion_order() const noexcept {
    flint::Fmpz out;
    if (!torsion_order(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool SUnitGroup::torsion_generator(OrderElement& out) const noexcept {
    return is_defined() && storage_->ordinary_units.torsion_generator(out);
}

slong SUnitGroup::ordinary_free_rank() const noexcept {
    return is_defined() ? storage_->ordinary_units.free_rank() : -1;
}

slong SUnitGroup::nonunit_rank() const noexcept {
    return is_defined()
            ? static_cast<slong>(storage_->nonunit_generators.size())
            : -1;
}

slong SUnitGroup::free_rank() const noexcept {
    return is_defined() ? ordinary_free_rank() + nonunit_rank() : -1;
}

slong SUnitGroup::generator_count() const noexcept {
    return is_defined() ? 1 + free_rank() : -1;
}

bool SUnitGroup::ordinary_free_generator(FactoredElement& out,
                                         slong index) const noexcept {
    return is_defined() && storage_->ordinary_units.free_generator(out, index);
}

bool SUnitGroup::nonunit_generator(FactoredElement& out,
                                   slong index) const noexcept {
    return is_defined() && index >= 0 && index < nonunit_rank() &&
           out.set(storage_->nonunit_generators[
                   static_cast<std::size_t>(index)]);
}

bool SUnitGroup::nonunit_valuation_row(flint::FmpzMatRef out,
                                       slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= nonunit_rank() ||
        flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != selected_prime_count()) {
        return false;
    }
    for (slong j = 0; j < selected_prime_count(); ++j) {
        flint::fmpz_set(
                flint::fmpz_mat_entry(out, 0, j),
                flint::fmpz_mat_entry(
                        flint::FmpzMatConstRef(storage_->nonunit_valuations),
                        index, j));
    }
    return true;
}

bool SUnitGroup::nonunit_valuation_matrix(
        flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || flint::fmpz_mat_nrows(out) != nonunit_rank() ||
        flint::fmpz_mat_ncols(out) != selected_prime_count()) {
        return false;
    }
    flint::fmpz_mat_set(out,
                        flint::FmpzMatConstRef(storage_->nonunit_valuations));
    return true;
}

bool SUnitGroup::regulator(flint::ArbRef out) const noexcept {
    if (!is_defined() ||
        storage_->regulator_proof_status != ProofState::verified) {
        return false;
    }
    flint::arb_set(out,
                   flint::ArbConstRef(storage_->regulator.raw()));
    return true;
}

std::optional<flint::Arb> SUnitGroup::regulator() const noexcept {
    flint::Arb out;
    if (!regulator(flint::ArbRef(out))) {
        return std::nullopt;
    }
    return out;
}

slong SUnitGroup::regulator_precision() const noexcept {
    return is_defined() ? storage_->regulator_precision : 0;
}

CertificationMode SUnitGroup::certification_status() const noexcept {
    return is_defined() && storage_->proof_status == ProofState::verified &&
                    storage_->regulator_proof_status == ProofState::verified
            ? CertificationMode::proven
            : CertificationMode::unknown;
}

CertificationMode SUnitGroup::source_class_certification() const noexcept {
    return is_defined() ? storage_->source_class_certification
                        : CertificationMode::unknown;
}

CertificationMode SUnitGroup::source_unit_certification() const noexcept {
    return is_defined() ? storage_->source_unit_certification
                        : CertificationMode::unknown;
}

ProofState SUnitGroup::source_relation_saturation_status() const noexcept {
    return is_defined() ? storage_->source_relation_saturation
                        : ProofState::not_checked;
}

ProofState SUnitGroup::source_unit_proof_status() const noexcept {
    return is_defined() ? storage_->source_unit_proof
                        : ProofState::not_checked;
}

ProofState SUnitGroup::source_regulator_proof_status() const noexcept {
    return is_defined() ? storage_->source_regulator_proof
                        : ProofState::not_checked;
}

ProofState SUnitGroup::proof_status() const noexcept {
    return is_defined() ? storage_->proof_status : ProofState::not_checked;
}

ProofState SUnitGroup::regulator_proof_status() const noexcept {
    return is_defined() ? storage_->regulator_proof_status
                        : ProofState::not_checked;
}

bool SUnitGroup::image(FactoredElement& out,
                       const SUnitCoordinates& coordinates) const noexcept {
    const NumberField* field = is_defined() ? storage_->order.parent() : nullptr;
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*field)) {
        return false;
    }
    FactoredElement candidate(*field);
    if (!candidate.is_defined() ||
        !compose_sunit_image(candidate, *storage_, coordinates)) {
        return false;
    }
    out.swap(candidate);
    return true;
}

bool SUnitGroup::image(Element& out,
                       const SUnitCoordinates& coordinates) const noexcept {
    const NumberField* field = is_defined() ? storage_->order.parent() : nullptr;
    if (field == nullptr || out.parent() == nullptr ||
        !out.parent()->has_same_data(*field)) {
        return false;
    }
    FactoredElement factored(*field);
    Element candidate(*field);
    if (!image(factored, coordinates) || !factored.evaluate(candidate)) {
        return false;
    }
    out.swap(candidate);
    return true;
}

bool SUnitGroup::preimage(SUnitMembershipResult& result,
                          SUnitCoordinates& out,
                          const FactoredElement& value,
                          EmbeddingContext& embeddings,
                          slong start_precision,
                          slong max_precision) const noexcept {
    if (!is_defined()) {
        result = {};
        membership_fail(result, SUnitMembershipStage::input_validation);
        return false;
    }
    return preimage_core(result, out, *storage_, value, embeddings,
                         start_precision, max_precision);
}

bool SUnitGroup::preimage(SUnitMembershipResult& result,
                          SUnitCoordinates& out,
                          const Element& value,
                          EmbeddingContext& embeddings,
                          slong start_precision,
                          slong max_precision) const noexcept {
    result = {};
    const NumberField* field = value.parent();
    if (!is_defined() || field == nullptr) {
        membership_fail(result, SUnitMembershipStage::input_validation);
        return false;
    }
    FactoredElement factored(*field);
    if (!factored.is_defined() || !factored.set_element(value)) {
        membership_fail(result, SUnitMembershipStage::input_validation);
        return false;
    }
    return preimage_core(result, out, *storage_, factored, embeddings,
                         start_precision, max_precision);
}

const char* sunit_compute_stage_name(SUnitComputeStage stage) noexcept {
    switch (stage) {
        case SUnitComputeStage::none:
            return "none";
        case SUnitComputeStage::input_validation:
            return "input_validation";
        case SUnitComputeStage::s_class_context:
            return "s_class_context";
        case SUnitComputeStage::s_class_publication:
            return "s_class_publication";
        case SUnitComputeStage::s_unit_publication:
            return "s_unit_publication";
        case SUnitComputeStage::s_regulator:
            return "s_regulator";
        case SUnitComputeStage::exact_verification:
            return "exact_verification";
    }
    return "unknown";
}

const char* sunit_membership_stage_name(
        SUnitMembershipStage stage) noexcept {
    switch (stage) {
        case SUnitMembershipStage::none:
            return "none";
        case SUnitMembershipStage::input_validation:
            return "input_validation";
        case SUnitMembershipStage::selected_valuation:
            return "selected_valuation";
        case SUnitMembershipStage::valuation_solve:
            return "valuation_solve";
        case SUnitMembershipStage::residual_unit:
            return "residual_unit";
        case SUnitMembershipStage::precision_exhausted:
            return "precision_exhausted";
        case SUnitMembershipStage::exact_verification:
            return "exact_verification";
    }
    return "unknown";
}

bool compute_sunit_groups(
        SUnitComputeResult& result,
        SClassGroup& s_class_group,
        SUnitGroup& s_unit_group,
        const ClassGroupContext& class_group,
        const OrderUnitGroup& ordinary_units,
        PrimeIdealSpan selected_primes,
        const SUnitComputeOptions& options) noexcept {
    result = {};
    const Order* order = class_group.parent();
    if (order == nullptr || order->parent() == nullptr || !order->is_maximal() ||
        !same_order_parent(ordinary_units.parent(), order) ||
        class_group.certification_status() != CertificationMode::proven ||
        ordinary_units.certification_status() != CertificationMode::proven ||
        options.regulator_precision <= 0) {
        compute_fail(result, SUnitComputeStage::input_validation);
        return false;
    }

    slong expected_rank = -1;
    if (!unit_rank(expected_rank, *order->parent()) ||
        ordinary_units.free_rank() != expected_rank) {
        compute_fail(result, SUnitComputeStage::input_validation);
        return false;
    }

    std::vector<PrimeIdeal> selected;
    slong failed_index = -1;
    if (!copy_selected_primes(selected, *order, selected_primes,
                              failed_index)) {
        compute_fail(result, SUnitComputeStage::input_validation,
                     failed_index);
        return false;
    }

    detail::SUnitClassContext context;
    detail::SUnitClassBuildResult context_result;
    if (!detail::build_sunit_class_context(context_result, context,
                                            class_group, selected)) {
        compute_fail(result, SUnitComputeStage::s_class_context,
                     context_result.selected_index);
        return false;
    }
    if (!context.defined ||
        context.s_class_proof_status != ProofState::verified ||
        context.s_unit_mod_units_proof_status != ProofState::verified) {
        compute_fail(result, SUnitComputeStage::exact_verification);
        return false;
    }

    SClassGroup class_candidate;
    class_candidate.storage_.reset(
            new (std::nothrow) detail::SClassGroupStorage());
    if (class_candidate.storage_ == nullptr ||
        !publish_s_class_storage(*class_candidate.storage_, context,
                                 options)) {
        compute_fail(result, SUnitComputeStage::s_class_publication);
        return false;
    }

    SUnitGroup unit_candidate;
    unit_candidate.storage_.reset(
            new (std::nothrow) detail::SUnitGroupStorage(context.order));
    if (unit_candidate.storage_ == nullptr) {
        compute_fail(result, SUnitComputeStage::s_unit_publication);
        return false;
    }
    SUnitComputeStage unit_failure_stage =
            SUnitComputeStage::s_unit_publication;
    if (!publish_s_unit_storage(*unit_candidate.storage_, unit_failure_stage,
                                context, class_group, ordinary_units,
                                options)) {
        compute_fail(result, unit_failure_stage);
        return false;
    }

    if (class_candidate.storage_->group.invariant_count() !=
                static_cast<slong>(
                        class_candidate.storage_->invariant_ideals.size()) ||
        class_candidate.storage_->group.invariant_count() !=
                static_cast<slong>(
                        class_candidate.storage_->power_witnesses.size()) ||
        flint::fmpz_mat_nrows(
                class_candidate.storage_->power_selected_exponents) !=
                class_candidate.storage_->group.invariant_count() ||
        flint::fmpz_mat_ncols(
                class_candidate.storage_->power_selected_exponents) !=
                static_cast<slong>(selected.size()) ||
        class_candidate.storage_->proof_status != ProofState::verified ||
        unit_candidate.storage_->ordinary_units.free_rank() != expected_rank ||
        static_cast<slong>(
                unit_candidate.storage_->nonunit_generators.size()) !=
                static_cast<slong>(selected.size()) ||
        flint::fmpz_mat_nrows(
                unit_candidate.storage_->nonunit_valuations) !=
                static_cast<slong>(selected.size()) ||
        flint::fmpz_mat_ncols(
                unit_candidate.storage_->nonunit_valuations) !=
                static_cast<slong>(selected.size()) ||
        unit_candidate.storage_->proof_status != ProofState::verified ||
        unit_candidate.storage_->regulator_proof_status !=
                ProofState::verified) {
        compute_fail(result, SUnitComputeStage::exact_verification);
        return false;
    }

    class_candidate.storage_->defined = true;
    unit_candidate.storage_->defined = true;
    s_class_group.swap(class_candidate);
    s_unit_group.swap(unit_candidate);
    result.success = true;
    result.stage = SUnitComputeStage::none;
    result.selected_index = -1;
    return true;
}

}  // namespace silex
