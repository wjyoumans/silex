#include <benchmark/benchmark.h>

#include <silex/class_group.hpp>
#include <silex/embedding.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_unit.hpp>

#include "benchmark_contract.hpp"
#include "factored_element/compact_reconstruction_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "order_unit/compact_reconstruction_bound_internal.hpp"

#include <array>
#include <chrono>
#include <cstddef>
#include <span>
#include <vector>

namespace {
namespace sflint = silex::flint;

enum class Workload {
    random_d15_h2_s43,
    random_d9_h4_s8,
    quartic_disc1412343,
};

struct CoordinateBoundContract {
    ulong coordinate_bound = 0;
    ulong bit_bound = 0;
    slong attempts = 0;
    slong precision = 0;
};

struct CoordinateBoundContractMatches {
    bool coordinate_bound = false;
    bool bit_bound = false;
    bool attempts = false;
    bool precision = false;

    bool all() const noexcept {
        return coordinate_bound && bit_bound && attempts && precision;
    }
};

// These are golden evaluations of the deterministic reference expression
//
//   2^ceil(log2(Minkowski bound) + log2(unit bound) + basis scale).
//
// They come from the fixed polynomial, maximal-order basis, and generator
// fixtures below, not from a warm call to the implementation under test.  In
// each case the Arb interval isolates the listed ceiling on the first 128-bit
// attempt.  The reordered quadratic basis is a unimodular permutation, while
// the index-three basis exercises the nontrivial normalization transform.
constexpr CoordinateBoundContract kQuarticDisc1412343BoundContract = {
        8'388'608, 23, 1, 128};
constexpr CoordinateBoundContract kLeadingOneQuadraticBoundContract = {
        2, 1, 1, 128};
constexpr CoordinateBoundContract kReorderedQuadraticBoundContract = {
        2, 1, 1, 128};
constexpr CoordinateBoundContract kIndexThreeQuarticBoundContract = {
        128, 7, 1, 128};

CoordinateBoundContractMatches match_coordinate_bound_contract(
        const CoordinateBoundContract& contract,
        const silex::detail::CompactCoordinateBoundReport& report,
        sflint::FmpzConstRef bound) noexcept {
    const bool output_matches =
            sflint::fmpz_cmp_ui(bound, contract.coordinate_bound) == 0;
    const bool report_matches =
            sflint::fmpz_cmp_ui(
                    sflint::FmpzConstRef(report.coordinate_bound),
                    contract.coordinate_bound) == 0;
    return {
            report.status ==
                            silex::detail::CompactCoordinateBoundStatus::
                                    success &&
                    output_matches && report_matches,
            report.bit_bound == contract.bit_bound,
            report.attempts == contract.attempts,
            report.precision == contract.precision,
    };
}

struct WorkloadSpec {
    std::span<const slong> coefficients;
    slong expected_class_order = 0;
    slong expected_unit_rank = 0;
    const CoordinateBoundContract* coordinate_bound_contract = nullptr;
};

constexpr std::array<slong, 16> kRandomD15H2S43 = {-2, 0, 0, 0, 0, 2, -1, 1,
                                                   2,  0, 2, 0, 2, 1, 1,  1};
constexpr std::array<slong, 10> kRandomD9H4S8 = {2,  -2, -3, -2, 4,
                                                 -2, -4, 1,  2,  1};
constexpr std::array<slong, 5> kQuarticDisc1412343 = {4, 1, 3, -8, 1};
constexpr std::array<slong, 3> kQuadraticFive = {-5, 0, 1};
constexpr std::array<slong, 5> kIndexThreeQuartic = {-5, 3, -2, 0, 1};

WorkloadSpec workload_spec(Workload workload) noexcept {
    switch (workload) {
        case Workload::random_d15_h2_s43:
            return {kRandomD15H2S43, 1, 8};
        case Workload::random_d9_h4_s8:
            return {kRandomD9H4S8, 1, 5};
        case Workload::quartic_disc1412343:
            return {kQuarticDisc1412343, 4, 2,
                    &kQuarticDisc1412343BoundContract};
    }
    return {};
}

class CompactUnitFixture {
   public:
    explicit CompactUnitFixture(Workload workload) {
        initialize(workload_spec(workload));
    }

    bool ready() const noexcept { return ready_; }
    const char* error() const noexcept { return error_; }
    slong degree() const noexcept { return maximal_order_.degree(); }
    slong unit_rank() const noexcept {
        return static_cast<slong>(generators_.size());
    }
    const std::vector<silex::FactoredElement>& generators() const noexcept {
        return generators_;
    }
    const std::vector<silex::Element>& expected() const noexcept {
        return expected_;
    }
    silex::Order& order() noexcept { return maximal_order_; }
    sflint::FmpzConstRef coordinate_bound() const noexcept {
        return sflint::FmpzConstRef(coordinate_bound_);
    }
    sflint::FmpzConstRef exact_coordinate_maximum() const noexcept {
        return sflint::FmpzConstRef(exact_coordinate_maximum_);
    }
    silex::EmbeddingContext& embeddings() noexcept { return embeddings_; }
    std::size_t factor_rows() const noexcept { return factor_rows_; }
    std::size_t primes_used() const noexcept { return primes_used_; }
    flint_bitcnt_t modulus_bits() const noexcept { return modulus_bits_; }
    flint_bitcnt_t coordinate_bound_bits() const noexcept {
        return fmpz_bits(coordinate_bound_.raw());
    }
    flint_bitcnt_t exact_coordinate_maximum_bits() const noexcept {
        return fmpz_bits(exact_coordinate_maximum_.raw());
    }
    ulong bit_bound() const noexcept {
        return bound_report_.bit_bound;
    }
    slong bound_attempts() const noexcept {
        return bound_report_.attempts;
    }
    slong bound_precision() const noexcept {
        return bound_report_.precision;
    }
    bool has_coordinate_bound_contract() const noexcept {
        return coordinate_bound_contract_ != nullptr;
    }
    CoordinateBoundContractMatches match_coordinate_bound_contract(
            const silex::detail::CompactCoordinateBoundReport& report,
            sflint::FmpzConstRef bound) const noexcept {
        if (coordinate_bound_contract_ == nullptr) {
            return {};
        }
        return ::match_coordinate_bound_contract(*coordinate_bound_contract_,
                                                 report, bound);
    }
    double initial_bound_milliseconds() const noexcept {
        return initial_bound_milliseconds_;
    }

   private:
    void fail(const char* message) noexcept {
        error_ = message;
        ready_ = false;
    }

    bool set_polynomial(std::span<const slong> coefficients) noexcept {
        if (coefficients.size() < 2) {
            return false;
        }
        sflint::FmpqPoly polynomial;
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                           coefficients[i]);
        }
        return field_.define_by_polynomial(
                sflint::FmpqPolyConstRef(polynomial));
    }

    bool configure_route(sflint::Fmpz& factor_base_bound,
                         silex::ClassGroupComputeOptions& options) noexcept {
        if (!silex::factor_base_class_group_bound(
                    sflint::FmpzRef(factor_base_bound), maximal_order_)) {
            return false;
        }
        if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) <
            0) {
            sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
        }
        options = silex::ClassGroupComputeOptions{};
        options.max_candidates = 5000;
        options.max_relations = maximal_order_.degree() >= 5 ? 1000 : 500;
        options.zeta_bf_max_cutoff = 20000;
        options.requested_certification = silex::CertificationMode::proven;
        return true;
    }

    bool derive_expected_outputs_and_bound() noexcept {
        sflint::fmpz_zero(sflint::FmpzRef(exact_coordinate_maximum_));
        expected_.clear();
        expected_.reserve(generators_.size());
        sflint::FmpqMat rational_coordinates(1, maximal_order_.degree());
        sflint::FmpzMat integral_coordinates(1, maximal_order_.degree());
        sflint::Fmpz absolute_value;
        for (const silex::FactoredElement& generator : generators_) {
            expected_.emplace_back(field_);
            if (!generator.evaluate(expected_.back()) ||
                !maximal_order_.coordinates(
                        sflint::FmpqMatRef(rational_coordinates),
                        expected_.back()) ||
                fmpq_mat_get_fmpz_mat(integral_coordinates.raw(),
                                      rational_coordinates.raw()) == 0) {
                return false;
            }
            for (slong column = 0; column < maximal_order_.degree(); ++column) {
                fmpz_abs(absolute_value.raw(),
                         fmpz_mat_entry(integral_coordinates.raw(), 0, column));
                if (fmpz_cmp(absolute_value.raw(),
                             exact_coordinate_maximum_.raw()) > 0) {
                    fmpz_set(exact_coordinate_maximum_.raw(),
                             absolute_value.raw());
                }
            }
        }
        return true;
    }

    bool derive_coordinate_bound() noexcept {
        const auto start = std::chrono::steady_clock::now();
        const bool success = silex::detail::compact_unit_coordinate_bound(
                bound_report_, sflint::FmpzRef(coordinate_bound_),
                maximal_order_, {generators_.data(), generators_.size()},
                embeddings_);
        initial_bound_milliseconds_ =
                std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - start)
                        .count();
        if (!success ||
            bound_report_.status !=
                    silex::detail::CompactCoordinateBoundStatus::success ||
            sflint::fmpz_cmp(sflint::FmpzConstRef(coordinate_bound_),
                             sflint::FmpzConstRef(exact_coordinate_maximum_)) <
                    0 ||
            (coordinate_bound_contract_ != nullptr &&
             !match_coordinate_bound_contract(
                      bound_report_,
                      sflint::FmpzConstRef(coordinate_bound_))
                      .all())) {
            return false;
        }
        return true;
    }

    bool outputs_equal(
            const std::vector<silex::Element>& values) const noexcept {
        if (values.size() != expected_.size()) {
            return false;
        }
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (!values[i].equal(expected_[i])) {
                return false;
            }
        }
        return true;
    }

    bool validate_reconstruction_checkpoint() noexcept {
        silex::detail::BoundedCompactReconstructionReport report;
        std::vector<silex::Element> reconstructed;
        if (!silex::detail::bounded_compact_reconstruct(
                    report, reconstructed, maximal_order_, generators_,
                    sflint::FmpzConstRef(coordinate_bound_)) ||
            report.status != silex::detail::BoundedCompactReconstructionStatus::
                                     success ||
            !outputs_equal(reconstructed)) {
            return false;
        }
        sflint::Fmpz twice_bound;
        sflint::Fmpz previous_modulus;
        fmpz_mul_2exp(twice_bound.raw(), coordinate_bound_.raw(), 1);
        if (report.last_prime == 0 ||
            fmpz_cmp(report.centered_crt_modulus.raw(), twice_bound.raw()) <=
                    0) {
            return false;
        }
        fmpz_divexact_ui(previous_modulus.raw(),
                         report.centered_crt_modulus.raw(), report.last_prime);
        if (fmpz_cmp(previous_modulus.raw(), twice_bound.raw()) > 0) {
            return false;
        }
        factor_rows_ = report.factor_rows;
        primes_used_ = report.primes_used;
        modulus_bits_ = report.modulus_bits;
        return true;
    }

    void initialize(const WorkloadSpec& spec) {
        coordinate_bound_contract_ = spec.coordinate_bound_contract;
        if (!set_polynomial(spec.coefficients) ||
            !equation_order_.define_equation_order(field_) ||
            !maximal_order_.define(field_) ||
            !maximal_order_.maximal_order(equation_order_)) {
            fail("field/maximal-order setup failed");
            return;
        }

        sflint::Fmpz factor_base_bound;
        silex::ClassGroupComputeOptions options;
        if (!configure_route(factor_base_bound, options)) {
            fail("strict-native options unavailable");
            return;
        }

        silex::detail::ClassUnitTransactionReport audit;
        if (!silex::detail::compute_class_unit_transaction(
                    units_, class_group_, maximal_order_,
                    sflint::FmpzConstRef(factor_base_bound), options, 128,
                    audit) || !audit.final_result_published ||
            units_.certification_status() != silex::CertificationMode::proven ||
            units_.free_rank() != spec.expected_unit_rank) {
            fail("strict-native source fixture failed");
            return;
        }

        sflint::Fmpz class_order;
        if (!class_group_.order(sflint::FmpzRef(class_order)) ||
            !sflint::fmpz_equal_si(class_order, spec.expected_class_order)) {
            fail("strict-native class order mismatch");
            return;
        }

        generators_.reserve(static_cast<std::size_t>(units_.free_rank()));
        for (slong i = 0; i < units_.free_rank(); ++i) {
            generators_.emplace_back(field_);
            if (!units_.free_generator(generators_.back(), i)) {
                fail("compact unit extraction failed");
                return;
            }
        }
        if (!derive_expected_outputs_and_bound()) {
            fail("direct exact bound derivation failed");
            return;
        }
        if (!derive_coordinate_bound()) {
            fail("coordinate bound differential failed");
            return;
        }
        if (!validate_reconstruction_checkpoint()) {
            fail("bounded reconstruction differential failed");
            return;
        }
        ready_ = true;
        error_ = nullptr;
    }

    silex::NumberField field_;
    silex::Order equation_order_;
    silex::Order maximal_order_;
    silex::ClassGroupContext class_group_;
    silex::OrderUnitGroup units_;
    std::vector<silex::FactoredElement> generators_;
    std::vector<silex::Element> expected_;
    silex::EmbeddingContext embeddings_;
    sflint::Fmpz exact_coordinate_maximum_;
    sflint::Fmpz coordinate_bound_;
    silex::detail::CompactCoordinateBoundReport bound_report_;
    const CoordinateBoundContract* coordinate_bound_contract_ = nullptr;
    std::size_t factor_rows_ = 0;
    std::size_t primes_used_ = 0;
    flint_bitcnt_t modulus_bits_ = 0;
    double initial_bound_milliseconds_ = 0.0;
    const char* error_ = "fixture not initialized";
    bool ready_ = false;
};

enum class BoundBasisWorkload {
    leading_one_quadratic,
    reordered_quadratic,
    index_three_quartic,
};

bool define_field(silex::NumberField& field,
                  std::span<const slong> coefficients) noexcept {
    if (coefficients.size() < 2) {
        return false;
    }
    sflint::FmpqPoly polynomial;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, static_cast<slong>(i),
                                       coefficients[i]);
    }
    return field.define_by_polynomial(sflint::FmpqPolyConstRef(polynomial));
}

bool basis_begins_with_one(const silex::Order& order) noexcept {
    sflint::FmpqMat basis(order.degree(), order.degree());
    if (!order.get_basis(sflint::FmpqMatRef(basis))) {
        return false;
    }
    for (slong column = 0; column < order.degree(); ++column) {
        const slong expected = column == 0 ? 1 : 0;
        if (!sflint::fmpq_equal_si(
                    sflint::fmpq_mat_entry(sflint::FmpqMatConstRef(basis), 0,
                                           column),
                    expected)) {
            return false;
        }
    }
    return true;
}

class BoundBasisFixture {
   public:
    explicit BoundBasisFixture(BoundBasisWorkload workload) {
        initialize(workload);
    }

    bool ready() const noexcept { return ready_; }
    const char* error() const noexcept { return error_; }
    silex::Order& order() noexcept { return order_; }
    silex::EmbeddingContext& embeddings() noexcept { return embeddings_; }
    const std::vector<silex::FactoredElement>& generators() const noexcept {
        return generators_;
    }
    slong degree() const noexcept { return order_.degree(); }
    std::size_t unit_rank() const noexcept { return generators_.size(); }
    ulong equation_index() const noexcept { return equation_index_; }
    bool leading_one_basis() const noexcept { return leading_one_basis_; }
    ulong bit_bound() const noexcept {
        return coordinate_bound_contract_->bit_bound;
    }
    slong bound_attempts() const noexcept {
        return coordinate_bound_contract_->attempts;
    }
    slong bound_precision() const noexcept {
        return coordinate_bound_contract_->precision;
    }
    flint_bitcnt_t coordinate_bound_bits() const noexcept {
        return fmpz_bits(expected_bound_.raw());
    }
    flint_bitcnt_t exact_coordinate_maximum_bits() const noexcept {
        return fmpz_bits(exact_coordinate_maximum_.raw());
    }

    CoordinateBoundContractMatches match_coordinate_bound_contract(
            const silex::detail::CompactCoordinateBoundReport& report,
            sflint::FmpzConstRef bound) const noexcept {
        if (coordinate_bound_contract_ == nullptr) {
            return {};
        }
        return ::match_coordinate_bound_contract(*coordinate_bound_contract_,
                                                 report, bound);
    }

    bool verify(
            const silex::detail::CompactCoordinateBoundReport& report,
            sflint::FmpzConstRef bound) const noexcept {
        return match_coordinate_bound_contract(report, bound).all() &&
               sflint::fmpz_equal(bound,
                                  sflint::FmpzConstRef(expected_bound_)) &&
               sflint::fmpz_cmp(
                       bound,
                       sflint::FmpzConstRef(exact_coordinate_maximum_)) >= 0;
    }

   private:
    void fail(const char* message) noexcept {
        error_ = message;
        ready_ = false;
    }

    bool define_orders(std::span<const slong> coefficients,
                       ulong expected_index) noexcept {
        if (!define_field(field_, coefficients) ||
            !equation_order_.define_equation_order(field_) ||
            !maximal_order_.define(field_) ||
            !maximal_order_.maximal_order(equation_order_)) {
            return false;
        }

        sflint::Fmpz index;
        if (!silex::order_index(sflint::FmpzRef(index), equation_order_,
                                maximal_order_) ||
            sflint::fmpz_cmp_ui(sflint::FmpzConstRef(index),
                                expected_index) != 0) {
            return false;
        }
        equation_index_ = expected_index;
        return true;
    }

    bool use_maximal_order(bool reorder_basis) noexcept {
        if (!reorder_basis) {
            return order_.set(maximal_order_);
        }

        sflint::FmpqMat basis(maximal_order_.degree(),
                              maximal_order_.degree());
        if (!maximal_order_.get_basis(sflint::FmpqMatRef(basis))) {
            return false;
        }
        for (slong column = 0; column < maximal_order_.degree(); ++column) {
            fmpq_swap(fmpq_mat_entry(basis.raw(), 0, column),
                      fmpq_mat_entry(basis.raw(), 1, column));
        }
        order_ = silex::Order::from_basis(
                field_, sflint::FmpqMatConstRef(basis));
        if (!order_.is_defined()) {
            return false;
        }
        order_.set_maximality(true);
        return true;
    }

    bool add_quadratic_unit() noexcept {
        silex::Element theta(field_);
        silex::Element one(field_);
        silex::Element numerator(field_);
        silex::Element unit(field_);
        if (!theta.gen() || !one.one() || !numerator.add(theta, one) ||
            !unit.scalar_div_si(numerator, 2) || !order_.contains(unit)) {
            return false;
        }
        generators_.emplace_back(field_);
        return generators_.back().push(unit, 1);
    }

    bool add_trivial_generators(std::size_t count) noexcept {
        silex::Element one(field_);
        if (!one.one()) {
            return false;
        }
        generators_.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            generators_.emplace_back(field_);
            if (!generators_.back().push(one, 1)) {
                return false;
            }
        }
        return true;
    }

    bool derive_exact_coordinate_maximum() noexcept {
        sflint::fmpz_zero(sflint::FmpzRef(exact_coordinate_maximum_));
        sflint::FmpqMat rational_coordinates(1, order_.degree());
        sflint::FmpzMat integral_coordinates(1, order_.degree());
        sflint::Fmpz absolute_value;
        for (const silex::FactoredElement& generator : generators_) {
            silex::Element value(field_);
            if (!generator.evaluate(value) ||
                !order_.coordinates(sflint::FmpqMatRef(rational_coordinates),
                                    value) ||
                fmpq_mat_get_fmpz_mat(integral_coordinates.raw(),
                                      rational_coordinates.raw()) == 0) {
                return false;
            }
            for (slong column = 0; column < order_.degree(); ++column) {
                fmpz_abs(absolute_value.raw(),
                         fmpz_mat_entry(integral_coordinates.raw(), 0,
                                        column));
                if (fmpz_cmp(absolute_value.raw(),
                             exact_coordinate_maximum_.raw()) > 0) {
                    fmpz_set(exact_coordinate_maximum_.raw(),
                             absolute_value.raw());
                }
            }
        }
        return fmpz_sgn(exact_coordinate_maximum_.raw()) > 0;
    }

    bool verify_coordinate_bound_contract() noexcept {
        silex::detail::CompactCoordinateBoundReport report;
        sflint::Fmpz bound;
        if (!silex::detail::compact_unit_coordinate_bound(
                    report, sflint::FmpzRef(bound), order_,
                    {generators_.data(), generators_.size()}, embeddings_) ||
            !verify(report, sflint::FmpzConstRef(bound))) {
            return false;
        }
        return true;
    }

    void initialize(BoundBasisWorkload workload) {
        bool expected_leading_one = false;
        switch (workload) {
            case BoundBasisWorkload::leading_one_quadratic:
                expected_leading_one = true;
                coordinate_bound_contract_ =
                        &kLeadingOneQuadraticBoundContract;
                if (!define_orders(kQuadraticFive, 2) ||
                    !use_maximal_order(false) || !add_quadratic_unit()) {
                    fail("leading-one fixture setup failed");
                    return;
                }
                break;
            case BoundBasisWorkload::reordered_quadratic:
                coordinate_bound_contract_ =
                        &kReorderedQuadraticBoundContract;
                if (!define_orders(kQuadraticFive, 2) ||
                    !use_maximal_order(true) || !add_quadratic_unit()) {
                    fail("reordered fixture setup failed");
                    return;
                }
                break;
            case BoundBasisWorkload::index_three_quartic:
                coordinate_bound_contract_ =
                        &kIndexThreeQuarticBoundContract;
                if (!define_orders(kIndexThreeQuartic, 3) ||
                    !use_maximal_order(false) || !add_trivial_generators(2)) {
                    fail("index-three fixture setup failed");
                    return;
                }
                break;
        }

        sflint::Fmpz benchmark_order_index;
        if (!silex::order_index(sflint::FmpzRef(benchmark_order_index),
                                equation_order_, order_) ||
            sflint::fmpz_cmp_ui(
                    sflint::FmpzConstRef(benchmark_order_index),
                    equation_index_) != 0) {
            fail("benchmark order index mismatch");
            return;
        }
        sflint::fmpz_set_ui(sflint::FmpzRef(expected_bound_),
                            coordinate_bound_contract_->coordinate_bound);
        leading_one_basis_ = basis_begins_with_one(order_);
        if (leading_one_basis_ != expected_leading_one) {
            fail("unexpected leading-one basis shape");
            return;
        }
        if (!derive_exact_coordinate_maximum() ||
            !verify_coordinate_bound_contract()) {
            fail("coordinate-bound fixture verification failed");
            return;
        }
        ready_ = true;
        error_ = nullptr;
    }

    silex::NumberField field_;
    silex::Order equation_order_;
    silex::Order maximal_order_;
    silex::Order order_;
    std::vector<silex::FactoredElement> generators_;
    silex::EmbeddingContext embeddings_;
    sflint::Fmpz exact_coordinate_maximum_;
    sflint::Fmpz expected_bound_;
    const CoordinateBoundContract* coordinate_bound_contract_ = nullptr;
    ulong equation_index_ = 0;
    bool leading_one_basis_ = false;
    const char* error_ = "fixture not initialized";
    bool ready_ = false;
};

CompactUnitFixture& fixture(Workload workload) {
    switch (workload) {
        case Workload::random_d15_h2_s43: {
            static CompactUnitFixture value(workload);
            return value;
        }
        case Workload::random_d9_h4_s8: {
            static CompactUnitFixture value(workload);
            return value;
        }
        case Workload::quartic_disc1412343: {
            static CompactUnitFixture value(workload);
            return value;
        }
    }
    static CompactUnitFixture fallback(Workload::quartic_disc1412343);
    return fallback;
}

BoundBasisFixture& bound_basis_fixture(BoundBasisWorkload workload) {
    switch (workload) {
        case BoundBasisWorkload::leading_one_quadratic: {
            static BoundBasisFixture value(workload);
            return value;
        }
        case BoundBasisWorkload::reordered_quadratic: {
            static BoundBasisFixture value(workload);
            return value;
        }
        case BoundBasisWorkload::index_three_quartic: {
            static BoundBasisFixture value(workload);
            return value;
        }
    }
    static BoundBasisFixture fallback(
            BoundBasisWorkload::leading_one_quadratic);
    return fallback;
}

bool outputs_equal(const std::vector<silex::Element>& actual,
                   const std::vector<silex::Element>& expected) noexcept {
    if (actual.size() != expected.size()) {
        return false;
    }
    for (std::size_t i = 0; i < actual.size(); ++i) {
        if (!actual[i].equal(expected[i])) {
            return false;
        }
    }
    return true;
}

void set_coordinate_bound_contract_counters(
        benchmark::State& state,
        const CoordinateBoundContractMatches& matches) {
    state.counters["coordinate_bound_exact_match"] =
            matches.coordinate_bound ? 1.0 : 0.0;
    state.counters["bit_bound_exact_match"] =
            matches.bit_bound ? 1.0 : 0.0;
    state.counters["bound_attempts_exact_match"] =
            matches.attempts ? 1.0 : 0.0;
    state.counters["bound_precision_exact_match"] =
            matches.precision ? 1.0 : 0.0;
}

void set_common_counters(benchmark::State& state,
                         const CompactUnitFixture& source) {
    state.counters["degree"] = static_cast<double>(source.degree());
    state.counters["unit_rank"] = static_cast<double>(source.unit_rank());
    state.counters["factor_rows"] = static_cast<double>(source.factor_rows());
    state.counters["coordinate_bound_bits"] =
            static_cast<double>(source.coordinate_bound_bits());
    state.counters["exact_coordinate_maximum_bits"] =
            static_cast<double>(source.exact_coordinate_maximum_bits());
    state.counters["bit_bound"] =
            static_cast<double>(source.bit_bound());
    state.counters["bound_attempts"] =
            static_cast<double>(source.bound_attempts());
    state.counters["bound_precision"] =
            static_cast<double>(source.bound_precision());
    state.counters["initial_bound_ms"] =
            source.initial_bound_milliseconds();
    state.counters["crt_primes"] = static_cast<double>(source.primes_used());
    state.counters["crt_modulus_bits"] =
            static_cast<double>(source.modulus_bits());
}

void benchmark_coordinate_bound(benchmark::State& state,
                                     Workload workload) {
    silex::bench_contract::initialize(state);
    CompactUnitFixture& source = fixture(workload);
    if (!source.ready()) {
        silex::bench_contract::fail(
                state, source.error(),
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::detail::CompactCoordinateBoundReport report;
    sflint::Fmpz bound;
    bool success = true;
    for (auto _ : state) {
        (void)_;
        const auto start = std::chrono::steady_clock::now();
        success = silex::detail::compact_unit_coordinate_bound(
                report, sflint::FmpzRef(bound), source.order(),
                {source.generators().data(), source.generators().size()},
                source.embeddings());
        const auto elapsed = std::chrono::steady_clock::now() - start;
        state.SetIterationTime(std::chrono::duration<double>(elapsed).count());
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(bound.raw());
        benchmark::ClobberMemory();
        if (!success) {
            break;
        }
    }
    set_common_counters(state, source);
    state.counters["timed_bound_attempts"] =
            static_cast<double>(report.attempts);
    state.counters["timed_bound_precision"] =
            static_cast<double>(report.precision);
    const CoordinateBoundContractMatches exact_matches =
            source.match_coordinate_bound_contract(
                    report, sflint::FmpzConstRef(bound));
    if (source.has_coordinate_bound_contract()) {
        set_coordinate_bound_contract_counters(state, exact_matches);
    }
    state.SetItemsProcessed(state.iterations() * source.unit_rank());
    if (state.iterations() == 0 || !success) {
        silex::bench_contract::fail(
                state, "coordinate bound operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!sflint::fmpz_equal(sflint::FmpzConstRef(bound),
                            source.coordinate_bound()) ||
        sflint::fmpz_cmp(sflint::FmpzConstRef(bound),
                         source.exact_coordinate_maximum()) < 0 ||
        (source.has_coordinate_bound_contract() && !exact_matches.all())) {
        silex::bench_contract::fail(
                state, "coordinate bound differential failed",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void benchmark_bound_basis(benchmark::State& state,
                           BoundBasisWorkload workload) {
    silex::bench_contract::initialize(state);
    BoundBasisFixture& source = bound_basis_fixture(workload);
    if (!source.ready()) {
        silex::bench_contract::fail(
                state, source.error(),
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::detail::CompactCoordinateBoundReport report;
    sflint::Fmpz bound;
    bool success = true;
    for (auto _ : state) {
        (void)_;
        const auto start = std::chrono::steady_clock::now();
        success = silex::detail::compact_unit_coordinate_bound(
                report, sflint::FmpzRef(bound), source.order(),
                {source.generators().data(), source.generators().size()},
                source.embeddings());
        const auto elapsed = std::chrono::steady_clock::now() - start;
        state.SetIterationTime(std::chrono::duration<double>(elapsed).count());
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(bound.raw());
        benchmark::ClobberMemory();
        if (!success) {
            break;
        }
    }

    state.counters["degree"] = static_cast<double>(source.degree());
    state.counters["unit_rank"] =
            static_cast<double>(source.unit_rank());
    state.counters["equation_index"] =
            static_cast<double>(source.equation_index());
    state.counters["basis_leading_one"] =
            source.leading_one_basis() ? 1.0 : 0.0;
    state.counters["normalization_required"] =
            source.leading_one_basis() ? 0.0 : 1.0;
    state.counters["coordinate_bound_bits"] =
            static_cast<double>(source.coordinate_bound_bits());
    state.counters["exact_coordinate_maximum_bits"] =
            static_cast<double>(source.exact_coordinate_maximum_bits());
    state.counters["bit_bound"] = static_cast<double>(source.bit_bound());
    state.counters["bound_attempts"] =
            static_cast<double>(source.bound_attempts());
    state.counters["bound_precision"] =
            static_cast<double>(source.bound_precision());
    const CoordinateBoundContractMatches exact_matches =
            source.match_coordinate_bound_contract(
                    report, sflint::FmpzConstRef(bound));
    set_coordinate_bound_contract_counters(state, exact_matches);
    state.SetItemsProcessed(state.iterations());
    if (state.iterations() == 0 || !success) {
        silex::bench_contract::fail(
                state, "coordinate-bound basis operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!exact_matches.all() ||
        !source.verify(report, sflint::FmpzConstRef(bound))) {
        silex::bench_contract::fail(
                state, "coordinate-bound basis differential failed",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void benchmark_direct_expansion(benchmark::State& state, Workload workload) {
    silex::bench_contract::initialize(state);
    CompactUnitFixture& source = fixture(workload);
    if (!source.ready()) {
        silex::bench_contract::fail(
                state, source.error(),
                silex::bench_contract::FailureReason::setup);
        return;
    }

    std::vector<silex::Element> outputs;
    bool success = true;
    for (auto _ : state) {
        (void)_;
        const auto start = std::chrono::steady_clock::now();
        outputs.clear();
        outputs.reserve(source.generators().size());
        for (const silex::FactoredElement& generator : source.generators()) {
            outputs.emplace_back(*source.order().parent());
            if (!generator.evaluate(outputs.back())) {
                success = false;
                break;
            }
        }
        const auto elapsed = std::chrono::steady_clock::now() - start;
        state.SetIterationTime(std::chrono::duration<double>(elapsed).count());
        benchmark::DoNotOptimize(outputs.data());
        benchmark::ClobberMemory();
    }
    set_common_counters(state, source);
    state.SetItemsProcessed(state.iterations() * source.unit_rank());
    if (state.iterations() == 0 || !success) {
        silex::bench_contract::fail(
                state, "direct expansion operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!outputs_equal(outputs, source.expected())) {
        silex::bench_contract::fail(
                state, "direct expansion differential failed",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void benchmark_bounded_crt(benchmark::State& state, Workload workload) {
    silex::bench_contract::initialize(state);
    CompactUnitFixture& source = fixture(workload);
    if (!source.ready()) {
        silex::bench_contract::fail(
                state, source.error(),
                silex::bench_contract::FailureReason::setup);
        return;
    }

    silex::detail::BoundedCompactReconstructionReport report;
    std::vector<silex::Element> outputs;
    bool success = true;
    for (auto _ : state) {
        (void)_;
        const auto start = std::chrono::steady_clock::now();
        success = silex::detail::bounded_compact_reconstruct(
                report, outputs, source.order(), source.generators(),
                source.coordinate_bound());
        const auto elapsed = std::chrono::steady_clock::now() - start;
        state.SetIterationTime(std::chrono::duration<double>(elapsed).count());
        benchmark::DoNotOptimize(success);
        benchmark::DoNotOptimize(outputs.data());
        benchmark::ClobberMemory();
        if (!success) {
            break;
        }
    }
    set_common_counters(state, source);
    state.counters["prime_trials"] = static_cast<double>(report.prime_trials);
    state.counters["field_polynomial_rejections"] =
            static_cast<double>(report.field_polynomial_prime_rejections);
    state.counters["order_basis_rejections"] =
            static_cast<double>(report.order_basis_prime_rejections);
    state.counters["denominator_rejections"] =
            static_cast<double>(report.denominator_prime_rejections);
    state.counters["noninvertible_rejections"] =
            static_cast<double>(report.noninvertible_prime_rejections);
    state.SetItemsProcessed(state.iterations() * source.unit_rank());
    if (state.iterations() == 0 || !success) {
        silex::bench_contract::fail(
                state, "bounded CRT operation failed",
                silex::bench_contract::FailureReason::operation);
        return;
    }
    if (!outputs_equal(outputs, source.expected())) {
        silex::bench_contract::fail(
                state, "bounded CRT differential failed",
                silex::bench_contract::FailureReason::reference_mismatch);
        return;
    }
    silex::bench_contract::succeed(state);
}

void BM_direct_random_d15_h2_s43(benchmark::State& state) {
    benchmark_direct_expansion(state, Workload::random_d15_h2_s43);
}

void BM_crt_random_d15_h2_s43(benchmark::State& state) {
    benchmark_bounded_crt(state, Workload::random_d15_h2_s43);
}

void BM_bound_random_d15_h2_s43(benchmark::State& state) {
    benchmark_coordinate_bound(state, Workload::random_d15_h2_s43);
}

void BM_direct_random_d9_h4_s8(benchmark::State& state) {
    benchmark_direct_expansion(state, Workload::random_d9_h4_s8);
}

void BM_crt_random_d9_h4_s8(benchmark::State& state) {
    benchmark_bounded_crt(state, Workload::random_d9_h4_s8);
}

void BM_bound_random_d9_h4_s8(benchmark::State& state) {
    benchmark_coordinate_bound(state, Workload::random_d9_h4_s8);
}

void BM_direct_quartic_disc1412343(benchmark::State& state) {
    benchmark_direct_expansion(state, Workload::quartic_disc1412343);
}

void BM_crt_quartic_disc1412343(benchmark::State& state) {
    benchmark_bounded_crt(state, Workload::quartic_disc1412343);
}

void BM_bound_quartic_disc1412343(benchmark::State& state) {
    benchmark_coordinate_bound(state, Workload::quartic_disc1412343);
}

void BM_bound_basis_leading_one_quadratic(benchmark::State& state) {
    benchmark_bound_basis(state,
                          BoundBasisWorkload::leading_one_quadratic);
}

void BM_bound_basis_reordered_quadratic(benchmark::State& state) {
    benchmark_bound_basis(state, BoundBasisWorkload::reordered_quadratic);
}

void BM_bound_basis_index_three_quartic(benchmark::State& state) {
    benchmark_bound_basis(state, BoundBasisWorkload::index_three_quartic);
}

BENCHMARK(BM_bound_basis_leading_one_quadratic)
        ->UseManualTime()
        ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_bound_basis_reordered_quadratic)
        ->UseManualTime()
        ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_bound_basis_index_three_quartic)
        ->UseManualTime()
        ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_direct_random_d15_h2_s43)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_crt_random_d15_h2_s43)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_bound_random_d15_h2_s43)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_direct_random_d9_h4_s8)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_crt_random_d9_h4_s8)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_bound_random_d9_h4_s8)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_direct_quartic_disc1412343)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_crt_quartic_disc1412343)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_bound_quartic_disc1412343)
        ->Iterations(1)
        ->UseManualTime()
        ->Unit(benchmark::kMillisecond);

}  // namespace

BENCHMARK_MAIN();
