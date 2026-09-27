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

class OrderUnitGroupAccess;

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
    // Residue-character kernels and their rows.  A kernel row has one
    // column per free generator of `group`, in free_generator order, and,
    // when ell divides the torsion order w, one final column for the torsion
    // generator; its width is free_rank() + 1 when ell | w and free_rank()
    // otherwise.  Entries are reduced to [0, ell).  residue_dlog_kernel
    // returns a basis of the common kernel, over the given primes, of the
    // ell-th-power residue characters on these generators.  saturate_row
    // takes one such row (e_1, ..., e_r[, t]), tests whether
    // prod u_i^(e_i) * zeta^t is an ell-th power, and, when its root is a
    // unit of the order, adjoins it; rows of any other width are rejected.
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
    // residue_dlog_proof_kernel returns the same kernel as
    // residue_dlog_kernel: the same width and column order, with one
    // residue character per prime applied to every column.  The two differ
    // only in how each prime's column is evaluated (residue_dlog_kernel
    // tries a direct degree-one evaluation first, residue_dlog_proof_kernel
    // always uses the residue-field quotient log).  The kernel of a column
    // does not depend on which nontrivial character mod ell is chosen, so
    // whenever both calls succeed they publish the same matrix.
    bool residue_dlog_proof_kernel(flint::FmpzMat& out,
                                   PrimeIdealSpan primes,
                                   flint::FmpzConstRef ell) const noexcept;
    // select_saturation_primes tests each prime on the free generators
    // only, while the bounded saturation passes (select_saturation_primes_
    // with_kernel) test primes with the torsion column included when
    // ell | w.  Both require q = 1 mod ell (q the residue field size) and a
    // residue character on the free generators, and the torsion generator
    // is a unit modulo every prime, so the extra column rejects no prime:
    // the per-prime acceptance test agrees between the two.  Both selectors
    // also fall back to a second scan over every residue degree, the same
    // way, when the first scan does not fill target_len.  The scans can
    // still pick different primes, though: the public selector's first scan
    // also skips rational primes dividing the order discriminant, which the
    // bounded selector's first scan does not, and this is the only
    // difference between the two scans (e.g. Z[sqrt(7)], ell = 3).  If the
    // acceptance test changes, keep it in agreement between the two; the
    // discriminant-skip difference is a separate, tracked question.
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
    // saturate_local_once and the bounded saturation passes below test the
    // rows of the residue_dlog_kernel shape described above, so a unit that
    // is an ell-th power only up to torsion is found when ell | w.
    bool saturate_local_once(bool& changed,
                             const OrderUnitGroup& group,
                             PrimeIdealSpan primes,
                             flint::FmpzConstRef ell,
                             EmbeddingContext& embeddings,
                             slong precision) noexcept;
    // In the bounded saturation passes, stable=true means only that a pass
    // over the selected residue-character kernel adjoined no root.  It does
    // not mean the subgroup is saturated: kernel rows that are not l-th
    // powers, or whose l-th root is not a unit of this order, are skipped.
    // Saturation is certified only by a proof route that records a verified
    // unit-proof status, such as prove_local_saturated or prove_index_bound.
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

    // Internal unit-installation hooks (trusted generators, cached torsion)
    // are reachable only through the noninstalled detail::OrderUnitGroupAccess.
    friend class detail::OrderUnitGroupAccess;
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
