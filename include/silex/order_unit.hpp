#pragma once

#include <cstdint>
#include <cstddef>
#include <memory>
#include <optional>

#include <flint/flint.h>

#include <silex/factored_element.hpp>
#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/status.hpp>

namespace silex {

class ClassGroupContext;
struct ClassGroupComputeOptions;
class OrderUnitGroup;

struct OrderUnitProofRecord {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
    flint::Fmpz aux_prime_bound;
    slong local_primes = 0;
    bool changed = false;
};

class PrimeIdealSpan {
public:
    PrimeIdealSpan() noexcept = default;
    PrimeIdealSpan(const PrimeIdeal* data, std::size_t size) noexcept
        : data_(data),
          size_(size) {
    }

    const PrimeIdeal* data() const noexcept { return data_; }
    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    const PrimeIdeal& operator[](std::size_t index) const noexcept {
        return data_[index];
    }
    const PrimeIdeal* begin() const noexcept { return data_; }
    const PrimeIdeal* end() const noexcept {
        return data_ == nullptr ? nullptr : data_ + size_;
    }

private:
    const PrimeIdeal* data_ = nullptr;
    std::size_t size_ = 0;
};

class FactoredElementSpan {
public:
    FactoredElementSpan() noexcept = default;
    FactoredElementSpan(const FactoredElement* data, std::size_t size) noexcept
        : data_(data),
          size_(size) {
    }

    const FactoredElement* data() const noexcept { return data_; }
    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    const FactoredElement& operator[](std::size_t index) const noexcept {
        return data_[index];
    }
    const FactoredElement* begin() const noexcept { return data_; }
    const FactoredElement* end() const noexcept {
        return data_ == nullptr ? nullptr : data_ + size_;
    }

private:
    const FactoredElement* data_ = nullptr;
    std::size_t size_ = 0;
};

namespace detail {

class ClassGroupCertificationAccess;

struct UnitProofRecordData {
    flint::Fmpz ell;
    ProofState status = ProofState::not_checked;
    flint::Fmpz aux_prime_bound;
    slong local_primes = 0;
    bool changed = false;
};

bool order_unit_group_set_units_internal(
        OrderUnitGroup& out,
        const Order& order,
        FactoredElementSpan generators,
        EmbeddingContext& embeddings,
        slong precision,
        bool trusted,
        const flint::Fmpz* cached_torsion_order = nullptr,
        const OrderElement* cached_torsion_generator = nullptr) noexcept;

}  // namespace detail

class OrderUnitGroup {
public:
    OrderUnitGroup() noexcept = default;
    explicit OrderUnitGroup(const Order& parent) noexcept;
    ~OrderUnitGroup() noexcept;

    OrderUnitGroup(const OrderUnitGroup&) = delete;
    OrderUnitGroup& operator=(const OrderUnitGroup&) = delete;

    OrderUnitGroup(OrderUnitGroup&& other) noexcept;
    OrderUnitGroup& operator=(OrderUnitGroup&& other) noexcept;

    void swap(OrderUnitGroup& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const OrderUnitGroup& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    void set_diagnostics(const DiagnosticsContext* diagnostics) noexcept;
    const DiagnosticsContext* diagnostics() const noexcept;
    bool is_set() const noexcept;
    slong free_rank() const noexcept;
    CertificationMode certification_status() const noexcept;

    bool torsion_order(flint::FmpzRef out) const noexcept;
    std::optional<flint::Fmpz> torsion_order() const noexcept;
    bool torsion_generator(OrderElement& out) const noexcept;
    bool free_generator(FactoredElement& out, slong index) const noexcept;
    bool regulator(flint::ArbRef out) const noexcept;
    std::optional<flint::Arb> regulator() const noexcept;
    bool class_regulator_product(flint::ArbRef out,
                                 const ClassGroupContext& class_group,
                                 slong precision) const noexcept;
    bool class_regulator_index_bound(
            flint::FmpzRef out,
            const ClassGroupContext& class_group,
            flint::ArbConstRef analytic_class_regulator_product,
            slong precision) const noexcept;
    slong unit_proof_record_count() const noexcept;
    bool unit_proof_record(flint::FmpzRef ell,
                           ProofState& status,
                           flint::FmpzRef aux_prime_bound,
                           slong& local_primes,
                           bool& changed,
                           slong index) const noexcept;
    std::optional<OrderUnitProofRecord> unit_proof_record(
            slong index) const noexcept;
    bool unit_proof_verified(flint::FmpzConstRef ell) const noexcept;
    bool regulator_index_bound(flint::FmpzRef out,
                               slong precision) const noexcept;

    // Computes a proven full unit group when an exact standalone route is
    // available; see the release support matrix. Failure preserves the
    // current object.
    bool compute(const Order& order) noexcept;
    // Publishes both groups atomically at the requested proven/GRH level.
    bool compute_with_class_group(ClassGroupContext& class_group,
                                  const Order& order,
                                  flint::FmpzConstRef factor_base_bound,
                                  const ClassGroupComputeOptions& options,
                                  slong precision) noexcept;
    // set_units and positive-rank relation-kernel operations publish exact
    // full-rank subgroups with unknown certification. Rank-zero relation-kernel
    // construction delegates to the proven compute() path. Invalid input and
    // construction failure preserve output.
    bool set_units(const Order& order,
                   FactoredElementSpan generators,
                   EmbeddingContext& embeddings,
                   slong precision) noexcept;
    bool set_relation_kernel_units(const Order& order,
                                   const ClassGroupContext& class_group,
                                   EmbeddingContext& embeddings,
                                   slong precision) noexcept;
    // Bounded refinements deterministically retain their last valid subgroup
    // when no further verified refinement is available.
    bool set_relation_kernel_units_bounded(
            const Order& order,
            const ClassGroupContext& class_group,
            EmbeddingContext& embeddings,
            flint::FmpzConstRef denominator_bound,
            slong start_precision,
            slong max_precision) noexcept;
    bool set_relation_kernel_units_index_bounded(
            const Order& order,
            const ClassGroupContext& class_group,
            EmbeddingContext& embeddings,
            slong start_precision,
            slong max_precision) noexcept;
    // A failed optional saturation phase publishes the index-bounded subgroup
    // with changed=false and stable=false; input failure preserves output.
    bool set_relation_kernel_units_index_bounded_saturated(
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
            slong max_passes) noexcept;
    bool saturate_row(bool& changed,
                      const OrderUnitGroup& group,
                      flint::FmpzMatConstRef kernel_rows,
                      slong row,
                      flint::FmpzConstRef ell,
                      EmbeddingContext& embeddings,
                      slong precision) noexcept;
    bool residue_dlog_kernel(flint::FmpzMat& out,
                             PrimeIdealSpan primes,
                             flint::FmpzConstRef ell) const noexcept;
    bool residue_dlog_proof_kernel(flint::FmpzMat& out,
                                   PrimeIdealSpan primes,
                                   flint::FmpzConstRef ell) const noexcept;
    bool select_saturation_primes(PrimeIdealList& out,
                                  flint::FmpzConstRef ell,
                                  slong target_len,
                                  flint::FmpzConstRef bound) const noexcept;
    bool select_saturation_proof_primes(PrimeIdealList& out,
                                        bool& certified,
                                        flint::FmpzMat& kernel,
                                        flint::FmpzConstRef ell,
                                        flint::FmpzConstRef bound)
            const noexcept;
    bool saturate_local_once(bool& changed,
                             const OrderUnitGroup& group,
                             PrimeIdealSpan primes,
                             flint::FmpzConstRef ell,
                             EmbeddingContext& embeddings,
                             slong precision) noexcept;
    bool saturate_bounded(bool& changed,
                          bool& stable,
                          const OrderUnitGroup& group,
                          flint::FmpzConstRef ell,
                          slong aux_target_len,
                          flint::FmpzConstRef aux_bound,
                          slong max_passes,
                          EmbeddingContext& embeddings,
                          slong precision) noexcept;
    bool saturate_index_bounded(bool& changed,
                                bool& stable,
                                const OrderUnitGroup& group,
                                EmbeddingContext& embeddings,
                                slong aux_target_len,
                                flint::FmpzConstRef aux_bound,
                                slong max_passes,
                                slong precision) noexcept;
    bool saturate_index_bounded_adaptive(
            bool& changed,
            bool& stable,
            const OrderUnitGroup& group,
            EmbeddingContext& embeddings,
            slong aux_target_len,
            flint::FmpzConstRef aux_bound_start,
            flint::FmpzConstRef aux_bound_max,
            slong max_passes,
            slong precision) noexcept;
    bool prove_local_saturated(ProofState& status,
                               bool& changed,
                               const OrderUnitGroup& group,
                               flint::FmpzConstRef ell,
                               slong aux_target_len,
                               flint::FmpzConstRef aux_bound,
                               EmbeddingContext& embeddings,
                               slong precision) noexcept;
    bool prove_index_bound(ProofState& status,
                           bool& changed,
                           const OrderUnitGroup& group,
                           slong aux_target_len,
                           flint::FmpzConstRef aux_bound,
                           slong max_restarts,
                           EmbeddingContext& embeddings,
                           slong precision) noexcept;

private:
    friend class ClassGroupContext;
    friend class detail::ClassGroupCertificationAccess;

    void mark_certification_proven_() noexcept;

    friend bool detail::order_unit_group_set_units_internal(
            OrderUnitGroup& out,
            const Order& order,
            FactoredElementSpan generators,
            EmbeddingContext& embeddings,
            slong precision,
            bool trusted,
            const flint::Fmpz* cached_torsion_order,
            const OrderElement* cached_torsion_generator) noexcept;
    bool mark_unit_proof(flint::FmpzConstRef ell,
                         ProofState status,
                         flint::FmpzConstRef aux_prime_bound,
                         slong local_primes,
                         bool changed) noexcept;
    bool mark_unit_proof_after_(
            flint::FmpzConstRef ell,
            ProofState status,
            flint::FmpzConstRef aux_prime_bound,
            slong local_primes,
            bool changed,
            slong stable_prefix_len) noexcept;
    void reset_unit_proof_records() noexcept;
    void try_certify_index_one(slong precision) noexcept;
    bool prove_local_saturated_(
            ProofState& status,
            bool& changed,
            const OrderUnitGroup& group,
            flint::FmpzConstRef ell,
            slong aux_target_len,
            flint::FmpzConstRef aux_bound,
            EmbeddingContext& embeddings,
            slong precision,
            bool use_stable_proof_fallback) noexcept;
    bool saturate_local_with_kernel_(bool& changed,
                                     const OrderUnitGroup& group,
                                     flint::FmpzMatConstRef kernel,
                                     flint::FmpzConstRef ell,
                                     EmbeddingContext& embeddings,
                                     slong precision) noexcept;
    bool prove_local_saturated_stable_in_place_(
            ProofState& status,
            bool& changed,
            OrderUnitGroup& group,
            flint::FmpzConstRef ell,
            slong aux_target_len,
            flint::FmpzConstRef aux_bound,
            EmbeddingContext& embeddings,
            slong precision,
            slong stable_proof_record_prefix_len) noexcept;
    bool compute_with_relation_class_group_(
            ClassGroupContext& class_group,
            const Order& order,
            flint::FmpzConstRef factor_base_bound,
            const ClassGroupComputeOptions& options,
            slong precision,
            slong rank,
            const DiagnosticsContext* active_diagnostics) noexcept;

    bool reserve_free_generators_(slong capacity) noexcept;
    bool append_free_generator_copy_(const FactoredElement& generator) noexcept;
    void clear_free_generators_() noexcept;
    bool reserve_unit_proof_records_(slong capacity) noexcept;
    detail::UnitProofRecordData* append_unit_proof_record_() noexcept;
    bool copy_unit_proof_records_from_(const OrderUnitGroup& other) noexcept;
    void clear_unit_proof_records_() noexcept;

    Order parent_;
    flint::Fmpz torsion_order_;
    OrderElement torsion_generator_;
    std::unique_ptr<FactoredElement[]> free_generators_;
    slong free_generator_count_ = 0;
    slong free_generator_capacity_ = 0;
    flint::Arb regulator_;
    bool has_regulator_ = false;
    CertificationMode certification_ = CertificationMode::unknown;
    std::unique_ptr<detail::UnitProofRecordData[]> unit_proof_records_;
    slong unit_proof_record_count_ = 0;
    slong unit_proof_record_capacity_ = 0;
    bool is_set_ = false;
    const DiagnosticsContext* diagnostics_ = nullptr;
};

inline void swap(OrderUnitGroup& left, OrderUnitGroup& right) noexcept {
    left.swap(right);
}

}  // namespace silex
