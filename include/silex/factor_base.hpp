#pragma once

#include <flint/flint.h>

#include <silex/flint/fmpz.hpp>
#include <silex/order.hpp>
#include <silex/prime_ideal.hpp>

#include <vector>

namespace silex {

struct DiagnosticsContext;
namespace detail {
class FactorBaseBlockAccess;
}  // namespace detail

class FactorBase {
public:
    FactorBase() noexcept = default;
    explicit FactorBase(const Order& parent) noexcept;
    ~FactorBase() noexcept;

    FactorBase(const FactorBase&) = delete;
    FactorBase& operator=(const FactorBase&) = delete;

    FactorBase(FactorBase&& other) noexcept;
    FactorBase& operator=(FactorBase&& other) noexcept;

    void swap(FactorBase& other) noexcept;
    void clear() noexcept;
    bool define(const Order& parent) noexcept;
    bool set(const FactorBase& other) noexcept;

    bool is_defined() const noexcept;
    const Order* parent() const noexcept;
    slong length() const noexcept;
    slong rational_prime_block_count() const noexcept;
    bool rational_prime_blocks_are_complete() const noexcept;

    bool build(flint::FmpzConstRef bound) noexcept;
    bool build(flint::FmpzConstRef bound,
               const DiagnosticsContext* diagnostics) noexcept;
    bool build_prime_ideal_norm_bounded(flint::FmpzConstRef bound) noexcept;
    bool build_relation_completion_base(
            flint::FmpzConstRef bound) noexcept;
    bool build_lll_relation_base(
            flint::FmpzConstRef bound) noexcept;
    bool prime(PrimeIdeal& out, slong index) const noexcept;
    const PrimeIdeal* prime_at(slong index) const noexcept;
    bool rational_prime_block_data(flint::FmpzRef rational_prime,
                                   slong& length,
                                   slong block_index) const noexcept;
    bool rational_prime_block_index(slong& index,
                                    slong block_index,
                                    slong offset) const noexcept;
    bool rational_prime_block(flint::FmpzRef rational_prime,
                              slong& start,
                              slong& length,
                              slong block_index) const noexcept;
    bool rational_prime_block_index_for_prime(
            slong& block_index,
            flint::FmpzConstRef rational_prime) const noexcept;
    bool equal(const FactorBase& other) const noexcept;
    slong index(const PrimeIdeal& prime) const noexcept;
    bool contains(const PrimeIdeal& prime) const noexcept;

private:
    friend class ClassGroupContext;
    friend class detail::FactorBaseBlockAccess;

    struct PrimeBlock {
        flint::Fmpz rational_prime;
        slong start = 0;
        slong length = 0;
        std::vector<slong> indices;
        bool contiguous = true;
        bool complete = false;
    };

    bool append_prime(const PrimeIdeal& prime) noexcept;
    bool append_prime(PrimeIdeal&& prime) noexcept;
    bool build_maximal_imaginary_quadratic_relation_base(
            flint::FmpzConstRef bound,
            const DiagnosticsContext* diagnostics) noexcept;
    bool build_prime_ideal_norm_bounded_impl(
            flint::FmpzConstRef bound,
            bool include_inert_primes,
            slong degree_limit,
            bool retain_one_split_prime = false) noexcept;
    bool build_lll_relation_base_impl(
            flint::FmpzConstRef bound,
            bool incomplete) noexcept;
    bool add_block_index(flint::FmpzConstRef rational_prime,
                         slong prime_index) noexcept;
    bool rebuild_blocks() noexcept;
    slong block_position(flint::FmpzConstRef rational_prime,
                         bool& found) const noexcept;
    static void refresh_block(PrimeBlock& block) noexcept;

    Order parent_;
    std::vector<PrimeIdeal> primes_;
    std::vector<PrimeBlock> blocks_;
    bool complete_rational_prime_blocks_ = false;
};

bool factor_base_class_group_bound(flint::FmpzRef out,
                                   const Order& order) noexcept;

inline void swap(FactorBase& left, FactorBase& right) noexcept {
    left.swap(right);
}

}  // namespace silex
