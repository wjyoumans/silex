#pragma once

#include <memory>
#include <optional>

#include <flint/flint.h>

#include <silex/class_group.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/order_unit.hpp>

namespace silex {

namespace detail {
struct SClassGroupStorage;
struct SUnitGroupStorage;
}  // namespace detail

enum class SUnitComputeStage {
    none,
    input_validation,
    s_class_context,
    s_class_publication,
    s_unit_publication,
    s_regulator,
    exact_verification,
};

struct SUnitComputeOptions {
    slong regulator_precision = 128;
    const DiagnosticsContext* diagnostics = nullptr;
};

struct SUnitComputeResult {
    bool success = false;
    SUnitComputeStage stage = SUnitComputeStage::none;
    slong selected_index = -1;
};

enum class SUnitMembershipOutcome {
    unknown,
    not_sunit,
    verified,
};

enum class SUnitMembershipStage {
    none,
    input_validation,
    selected_valuation,
    valuation_solve,
    residual_unit,
    precision_exhausted,
    exact_verification,
};

struct SUnitMembershipResult {
    bool success = false;
    SUnitMembershipOutcome outcome = SUnitMembershipOutcome::unknown;
    SUnitMembershipStage stage = SUnitMembershipStage::none;
    slong work_precision = 0;
};

struct SUnitCoordinates {
    // Stable coordinate order: torsion, ordinary free units, then nonunits.
    flint::Fmpz torsion_exponent;
    flint::FmpzMat ordinary_free_exponents{1, 0};
    flint::FmpzMat nonunit_exponents{1, 0};
    bool defined = false;
};

class SClassGroup {
public:
    SClassGroup() noexcept;
    ~SClassGroup() noexcept;

    SClassGroup(const SClassGroup&) = delete;
    SClassGroup& operator=(const SClassGroup&) = delete;

    SClassGroup(SClassGroup&& other) noexcept;
    SClassGroup& operator=(SClassGroup&& other) noexcept;

    void swap(SClassGroup& other) noexcept;
    void clear() noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong selected_prime_count() const noexcept;
    bool selected_prime(PrimeIdeal& out, slong index) const noexcept;

    slong invariant_count() const noexcept;
    bool invariant(flint::FmpzRef out, slong index) const noexcept;
    std::optional<flint::Fmpz> invariant(slong index) const noexcept;
    bool order(flint::FmpzRef out) const noexcept;
    std::optional<flint::Fmpz> order() const noexcept;
    bool invariant_generator(FractionalIdeal& out, slong index) const noexcept;
    bool invariant_generator_power_witness(FactoredElement& out,
                                           slong index) const noexcept;
    bool invariant_generator_power_selected_exponents(
            flint::FmpzMatRef out,
            slong index) const noexcept;

    CertificationMode certification_status() const noexcept;
    CertificationMode source_class_certification() const noexcept;
    ProofState proof_status() const noexcept;

private:
    std::unique_ptr<detail::SClassGroupStorage> storage_;

    friend bool compute_sunit_groups(
            SUnitComputeResult& result,
            SClassGroup& s_class_group,
            class SUnitGroup& s_unit_group,
            const ClassGroupContext& class_group,
            const OrderUnitGroup& ordinary_units,
            PrimeIdealSpan selected_primes,
            const SUnitComputeOptions& options) noexcept;
};

class SUnitGroup {
public:
    SUnitGroup() noexcept;
    ~SUnitGroup() noexcept;

    SUnitGroup(const SUnitGroup&) = delete;
    SUnitGroup& operator=(const SUnitGroup&) = delete;

    SUnitGroup(SUnitGroup&& other) noexcept;
    SUnitGroup& operator=(SUnitGroup&& other) noexcept;

    void swap(SUnitGroup& other) noexcept;
    void clear() noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong selected_prime_count() const noexcept;
    bool selected_prime(PrimeIdeal& out, slong index) const noexcept;

    bool torsion_order(flint::FmpzRef out) const noexcept;
    std::optional<flint::Fmpz> torsion_order() const noexcept;
    bool torsion_generator(OrderElement& out) const noexcept;
    slong ordinary_free_rank() const noexcept;
    slong nonunit_rank() const noexcept;
    slong free_rank() const noexcept;
    slong generator_count() const noexcept;
    bool ordinary_free_generator(FactoredElement& out,
                                 slong index) const noexcept;
    bool nonunit_generator(FactoredElement& out, slong index) const noexcept;
    bool nonunit_valuation_row(flint::FmpzMatRef out,
                               slong index) const noexcept;
    bool nonunit_valuation_matrix(flint::FmpzMatRef out) const noexcept;

    bool regulator(flint::ArbRef out) const noexcept;
    std::optional<flint::Arb> regulator() const noexcept;
    slong regulator_precision() const noexcept;

    CertificationMode certification_status() const noexcept;
    CertificationMode source_class_certification() const noexcept;
    CertificationMode source_unit_certification() const noexcept;
    ProofState source_relation_saturation_status() const noexcept;
    ProofState source_unit_proof_status() const noexcept;
    ProofState source_regulator_proof_status() const noexcept;
    ProofState proof_status() const noexcept;
    ProofState regulator_proof_status() const noexcept;

    bool image(FactoredElement& out,
               const SUnitCoordinates& coordinates) const noexcept;
    bool image(Element& out,
               const SUnitCoordinates& coordinates) const noexcept;
    bool preimage(SUnitMembershipResult& result,
                  SUnitCoordinates& out,
                  const FactoredElement& value,
                  EmbeddingContext& embeddings,
                  slong start_precision,
                  slong max_precision) const noexcept;
    bool preimage(SUnitMembershipResult& result,
                  SUnitCoordinates& out,
                  const Element& value,
                  EmbeddingContext& embeddings,
                  slong start_precision,
                  slong max_precision) const noexcept;

private:
    std::unique_ptr<detail::SUnitGroupStorage> storage_;

    friend bool compute_sunit_groups(
            SUnitComputeResult& result,
            SClassGroup& s_class_group,
            SUnitGroup& s_unit_group,
            const ClassGroupContext& class_group,
            const OrderUnitGroup& ordinary_units,
            PrimeIdealSpan selected_primes,
            const SUnitComputeOptions& options) noexcept;
};

const char* sunit_compute_stage_name(SUnitComputeStage stage) noexcept;
const char* sunit_membership_stage_name(SUnitMembershipStage stage) noexcept;

bool compute_sunit_groups(
        SUnitComputeResult& result,
        SClassGroup& s_class_group,
        SUnitGroup& s_unit_group,
        const ClassGroupContext& class_group,
        const OrderUnitGroup& ordinary_units,
        PrimeIdealSpan selected_primes,
        const SUnitComputeOptions& options = {}) noexcept;

inline void swap(SClassGroup& left, SClassGroup& right) noexcept {
    left.swap(right);
}

inline void swap(SUnitGroup& left, SUnitGroup& right) noexcept {
    left.swap(right);
}

}  // namespace silex
