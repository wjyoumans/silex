#pragma once

#include <silex/diagnostics.hpp>
#include <silex/order.hpp>

#include <atomic>

#include "../zeta/zeta_internal.hpp"

namespace silex::detail {

class TestHookAccess;

struct RelationFactorBasePlan {
    flint::Fmpz working_bound;
    ZetaBfResidueDegreeCache residue_degrees;
    Order order;
    bool valid = false;
    // Set by specialize_maximal_imaginary_quadratic_factor_base_plan when
    // the base up to the generation bound would exceed the relation module
    // size limit: the working bound then stays compact, and the `proven`
    // class group is certified by its structure certificate instead of the
    // generation bound (ClassGroupGenerationBasis::structure_certificate).
    bool imaginary_quadratic_structure_certificate = false;
};

ZetaBfResidueDegreeCache* relation_factor_base_plan_residue_degrees(
        RelationFactorBasePlan& plan,
        const Order& order) noexcept;

const ZetaBfResidueDegreeCache* relation_factor_base_plan_residue_degrees(
        const RelationFactorBasePlan& plan,
        const Order& order) noexcept;

namespace relation_search {

bool build_relation_factor_base_plan(
        RelationFactorBasePlan& out,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool build_maximal_imaginary_quadratic_factor_base_plan(
        RelationFactorBasePlan& out,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool specialize_maximal_imaginary_quadratic_factor_base_plan(
        RelationFactorBasePlan& plan,
        const Order& order,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool relation_factor_base_restart_limit(
        flint::Fmpz& out,
        const Order& order) noexcept;

bool next_relation_factor_base_bound_to_limit(
        flint::Fmpz& out,
        flint::FmpzConstRef current,
        flint::FmpzConstRef limit) noexcept;

bool next_relation_factor_base_bound(
        flint::Fmpz& out,
        flint::FmpzConstRef current,
        flint::FmpzConstRef limit) noexcept;

bool quadratic_prime_is_bad(flint::FmpzConstRef discriminant,
                            flint::FmpzConstRef rational_prime) noexcept;

bool quadratic_nth_suitable_ideal_bound(
        flint::Fmpz& out,
        flint::FmpzConstRef discriminant,
        slong count) noexcept;

bool quadratic_working_bound_floor(slong& out,
                                   flint::FmpzConstRef discriminant) noexcept;

}  // namespace relation_search

// Test-only override for the structure-certificate switch of
// specialize_maximal_imaginary_quadratic_factor_base_plan: a positive value
// replaces kRelationModuleMaxColumns as the number of prime ideals up to the
// generation bound past which a `proven` imaginary quadratic plan keeps its
// compact base, so a test can take that route on a small field; zero means
// off.  Only detail::TestHookAccess (test/test_hook_access.hpp) can set it;
// production code only reads it.  Relaxed atomic so a reader on another
// thread never races with a test setting it.
class QuadraticStructureCertificateTestSeams {
public:
    static slong column_limit() noexcept {
        return column_limit_.load(std::memory_order_relaxed);
    }
    // Test-only: when set, a plan that keeps its compact base because the
    // base up to the generation bound exceeds the column limit does not
    // record the structure-certificate switch, so the quadratic backend runs
    // its factor-base honesty route on the compact base (a route no plan
    // reaches in production, where a compact base always carries the
    // switch).  Off by default.
    static bool compact_base_without_certificate() noexcept {
        return compact_base_without_certificate_.load(
                std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<slong> column_limit_{0};
    inline static std::atomic<bool> compact_base_without_certificate_{false};
};

// Test-only counters of the analytic-bound searches that stop doubling at the
// cap floor(4 log(|D|)^2) (the early return of compute_analytic_bound and
// compute_quadratic_analytic_bound), so a test can observe that the early
// return is taken.  Only detail::TestHookAccess reads or resets them;
// production code increments them with relaxed atomics and no decision
// reads them.
class AnalyticBoundCapSeams {
public:
    static void record(bool quadratic) noexcept {
        counts_[quadratic ? 1 : 0].fetch_add(1, std::memory_order_relaxed);
    }

private:
    friend class TestHookAccess;
    inline static std::atomic<slong> counts_[2]{};
};

}  // namespace silex::detail
