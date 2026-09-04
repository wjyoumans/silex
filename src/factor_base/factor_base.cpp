#include <silex/factor_base.hpp>

#include "factor_base_internal.hpp"
#include "../prime_ideal/prime_ideal_internal.hpp"

#include <silex/diagnostics.hpp>
#include <silex/signature.hpp>

#include <flint/fmpz.h>

#include <algorithm>
#include <utility>

namespace silex {
namespace {

void floor_sqrt_div_ui(flint::Fmpz& out,
                       flint::FmpzConstRef value,
                       ulong denominator) noexcept {
    flint::Fmpz quotient;
    fmpz_fdiv_q_ui(quotient.raw(), value.raw(), denominator);
    fmpz_sqrt(out.raw(), quotient.raw());
}

void ceil_sqrt(flint::Fmpz& out, flint::FmpzConstRef value) noexcept {
    flint::Fmpz remainder;
    fmpz_sqrtrem(out.raw(), remainder.raw(), value.raw());
    if (fmpz_is_zero(remainder.raw()) == 0) {
        fmpz_add_ui(out.raw(), out.raw(), 1);
    }
}

void generic_minkowski_bound(flint::Fmpz& out,
                             flint::FmpzConstRef abs_discriminant,
                             slong degree,
                             slong complex_pairs) noexcept {
    flint::Fmpz numerator;
    flint::Fmpz denominator;
    flint::Fmpz sqrt_discriminant;

    ceil_sqrt(sqrt_discriminant, abs_discriminant);
    fmpz_fac_ui(numerator.raw(), static_cast<ulong>(degree));
    fmpz_mul(numerator.raw(), numerator.raw(), sqrt_discriminant.raw());
    fmpz_mul_2exp(numerator.raw(), numerator.raw(),
                  static_cast<ulong>(complex_pairs));
    fmpz_set_ui(denominator.raw(), static_cast<ulong>(degree));
    fmpz_pow_ui(denominator.raw(), denominator.raw(),
                static_cast<ulong>(degree));
    fmpz_cdiv_q(out.raw(), numerator.raw(), denominator.raw());
    if (fmpz_is_zero(out.raw()) != 0) {
        fmpz_one(out.raw());
    }
}

bool prime_ideal_norm_at_most_bound(const PrimeIdeal& prime,
                                    flint::FmpzConstRef rational_prime,
                                    flint::FmpzConstRef bound) noexcept {
    const slong residue_degree = prime.residue_degree();
    if (residue_degree <= 0) {
        return false;
    }

    flint::Fmpz norm;
    fmpz_pow_ui(norm.raw(), rational_prime.raw(),
                static_cast<ulong>(residue_degree));
    return fmpz_cmp(norm.raw(), bound.raw()) <= 0;
}

bool is_inert_prime_decomposition(const PrimeIdealList& decomposed,
                                  slong degree) noexcept {
    if (degree <= 0 || decomposed.size() != 1) {
        return false;
    }

    const PrimeIdeal* prime = decomposed.at(0);
    return prime != nullptr &&
           prime->ramification_index() == 1 &&
           prime->residue_degree() == degree;
}

}  // namespace

FactorBase::FactorBase(const Order& parent) noexcept {
    define(parent);
}

FactorBase::~FactorBase() noexcept = default;

FactorBase::FactorBase(FactorBase&& other) noexcept {
    swap(other);
}

FactorBase& FactorBase::operator=(FactorBase&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FactorBase::swap(FactorBase& other) noexcept {
    parent_.swap(other.parent_);
    primes_.swap(other.primes_);
    blocks_.swap(other.blocks_);
    std::swap(complete_rational_prime_blocks_,
              other.complete_rational_prime_blocks_);
}

void FactorBase::clear() noexcept {
    primes_.clear();
    blocks_.clear();
    parent_.clear();
    complete_rational_prime_blocks_ = false;
}

bool FactorBase::define(const Order& parent) noexcept {
    if (!parent.has_basis()) {
        return false;
    }

    clear();
    parent_ = parent;
    return true;
}

bool FactorBase::set(const FactorBase& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    FactorBase copy(other.parent_);
    if (!copy.is_defined()) {
        return false;
    }

    for (const PrimeIdeal& prime : other.primes_) {
        if (!copy.append_prime(prime)) {
            return false;
        }
    }
    if (copy.blocks_.size() != other.blocks_.size()) {
        return false;
    }
    for (std::size_t i = 0; i < copy.blocks_.size(); ++i) {
        copy.blocks_[i].complete = other.blocks_[i].complete;
    }
    copy.complete_rational_prime_blocks_ =
            other.complete_rational_prime_blocks_;

    swap(copy);
    return true;
}

bool FactorBase::is_defined() const noexcept {
    return parent_.has_basis();
}

const Order* FactorBase::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong FactorBase::length() const noexcept {
    return is_defined() ? static_cast<slong>(primes_.size()) : 0;
}

slong FactorBase::rational_prime_block_count() const noexcept {
    return is_defined() ? static_cast<slong>(blocks_.size()) : 0;
}

bool FactorBase::rational_prime_blocks_are_complete() const noexcept {
    return is_defined() && complete_rational_prime_blocks_;
}

bool detail::FactorBaseBlockAccess::rational_prime_block_is_complete(
        bool& complete,
        const FactorBase& base,
        slong block_index) noexcept {
    complete = false;
    if (!base.is_defined() || block_index < 0 ||
        block_index >= base.rational_prime_block_count()) {
        return false;
    }
    complete = base.complete_rational_prime_blocks_ ||
               base.blocks_[static_cast<std::size_t>(block_index)].complete;
    return true;
}

bool FactorBase::build(flint::FmpzConstRef bound) noexcept {
    return build(bound, nullptr);
}

bool FactorBase::build(flint::FmpzConstRef bound,
                       const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                        "factor_base.build_complete");
    if (!is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0 ||
        !parent_.is_maximal()) {
        return false;
    }

    FactorBase candidate(parent_);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList decomposed;
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "factor_base.decompose_prime");
            if (!decompose_prime(decomposed, parent_,
                                 flint::FmpzConstRef(p), 0,
                                 diagnostics)) {
                return false;
            }
        }

        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                                "factor_base.append_prime_block");
            for (slong i = 0; i < decomposed.size(); ++i) {
                const PrimeIdeal* prime = decomposed.at(i);
                if (prime == nullptr || !candidate.append_prime(*prime)) {
                    return false;
                }
            }
        }

        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    candidate.complete_rational_prime_blocks_ = true;
    swap(candidate);
    return true;
}

bool FactorBase::build_prime_ideal_norm_bounded(
        flint::FmpzConstRef bound) noexcept {
    return build_prime_ideal_norm_bounded_impl(bound, true, 0);
}

bool FactorBase::build_relation_completion_base(
        flint::FmpzConstRef bound) noexcept {
    return build_prime_ideal_norm_bounded_impl(bound, false, 0);
}

bool FactorBase::build_lll_relation_base(
        flint::FmpzConstRef bound) noexcept {
    return build_lll_relation_base_impl(bound, false);
}

bool FactorBase::build_maximal_imaginary_quadratic_relation_base(
        flint::FmpzConstRef bound,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(
            diagnostics, DiagnosticsModule::class_group,
            "factor_base.maximal_imaginary_quadratic.build");
    flint::Fmpz discriminant;
    if (parent_.degree() != 2 || !parent_.is_maximal() ||
        !parent_.discriminant(flint::FmpzRef(discriminant)) ||
        flint::fmpz_sgn(flint::FmpzConstRef(discriminant)) >= 0) {
        return false;
    }
    if (!is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0) {
        return false;
    }

    FactorBase candidate(parent_);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdeal retained(parent_);
        detail::RetainedQuadraticPrimeKind kind =
                detail::RetainedQuadraticPrimeKind::inert;
        bool retained_ready = false;
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "factor_base.maximal_imaginary_quadratic.classify_prime");
            retained_ready = retained.is_defined() &&
                    detail::MaximalQuadraticPrimeAccess::
                            set_first_degree_one_prime(
                                    retained, kind, parent_,
                                    flint::FmpzConstRef(p), diagnostics);
        }
        if (!retained_ready) {
            // The direct path has not published candidate. Preserve the
            // established generic decomposition route as the atomic fallback.
            return build_prime_ideal_norm_bounded_impl(bound, false, 0, true);
        }
        if (kind == detail::RetainedQuadraticPrimeKind::inert) {
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }
        {
            SILEX_PROFILE_SCOPE(
                    diagnostics, DiagnosticsModule::class_group,
                    "factor_base.maximal_imaginary_quadratic.append_prime");
            if (!candidate.append_prime(std::move(retained))) {
                return false;
            }

            bool found = false;
            const slong block = candidate.block_position(
                    flint::FmpzConstRef(p), found);
            if (!found) {
                return false;
            }
            candidate.blocks_[static_cast<std::size_t>(block)].complete =
                    kind == detail::RetainedQuadraticPrimeKind::ramified;
        }
        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    swap(candidate);
    return true;
}

bool FactorBase::build_lll_relation_base_impl(
        flint::FmpzConstRef bound,
        bool incomplete) noexcept {
    if (incomplete
                ? !build_prime_ideal_norm_bounded_impl(
                          bound, true, 0)
                : !build(bound)) {
        return false;
    }

    struct NormIndex {
        slong index = 0;
        flint::Fmpz norm;
    };

    std::vector<NormIndex> order;
    order.reserve(primes_.size());
    for (slong i = 0; i < static_cast<slong>(primes_.size()); ++i) {
        NormIndex entry;
        entry.index = i;
        if (!primes_[static_cast<std::size_t>(i)].norm(
                    flint::FmpzRef(entry.norm))) {
            return false;
        }
        order.emplace_back(std::move(entry));
    }

    std::sort(order.begin(), order.end(),
              [](const NormIndex& left, const NormIndex& right) noexcept {
                  const int cmp = flint::fmpz_cmp(
                          flint::FmpzConstRef(left.norm),
                          flint::FmpzConstRef(right.norm));
                  if (cmp != 0) {
                      return cmp > 0;
                  }
                  return left.index < right.index;
              });

    std::vector<PrimeIdeal> sorted;
    sorted.reserve(primes_.size());
    for (const NormIndex& entry : order) {
        PrimeIdeal copy(parent_);
        if (!copy.is_defined() ||
            !copy.set(primes_[static_cast<std::size_t>(entry.index)])) {
            return false;
        }
        sorted.emplace_back(std::move(copy));
    }

    primes_.swap(sorted);
    return rebuild_blocks();
}

bool FactorBase::build_prime_ideal_norm_bounded_impl(
        flint::FmpzConstRef bound,
        bool include_inert_primes,
        slong degree_limit,
        bool retain_one_split_prime) noexcept {
    if (!is_defined() || fmpz_cmp_ui(bound.raw(), 2) < 0 ||
        !parent_.is_maximal() || degree_limit < 0) {
        return false;
    }

    FactorBase candidate(parent_);
    flint::Fmpz p;
    fmpz_set_ui(p.raw(), 2);
    while (fmpz_cmp(p.raw(), bound.raw()) <= 0) {
        PrimeIdealList decomposed;
        if (!decompose_prime(decomposed, parent_,
                             flint::FmpzConstRef(p))) {
            return false;
        }
        // The source factor-base generator skips inert rational primes for the relation
        // factor base, even when p^degree is still below the tested bound.
        if (!include_inert_primes &&
            is_inert_prime_decomposition(decomposed, parent_.degree())) {
            fmpz_nextprime(p.raw(), p.raw(), 1);
            continue;
        }

        slong retained = 0;
        for (slong i = 0; i < decomposed.size(); ++i) {
            const PrimeIdeal* prime = decomposed.at(i);
            if (prime == nullptr) {
                return false;
            }
            if (degree_limit > 0 &&
                prime->residue_degree() > degree_limit) {
                continue;
            }
            if (!prime_ideal_norm_at_most_bound(
                        *prime, flint::FmpzConstRef(p), bound)) {
                continue;
            }
            if (!candidate.append_prime(*prime)) {
                return false;
            }
            ++retained;
            // reference's quadratic FBquad stores one prime form for a split
            // rational prime; its inverse represents the conjugate class.
            if (retain_one_split_prime && decomposed.size() == 2) {
                break;
            }
        }
        if (retained > 0) {
            // The source factor-base generator marks LV[p] when every prime in the cached
            // decomposition survives the relation-base filters.  Retain that
            // fact here so init_rel and subFBgen need not decompose p again.
            bool found = false;
            const slong block = candidate.block_position(
                    flint::FmpzConstRef(p), found);
            if (!found) {
                return false;
            }
            candidate.blocks_[static_cast<std::size_t>(block)].complete =
                    retained == decomposed.size();
        }

        fmpz_nextprime(p.raw(), p.raw(), 1);
    }

    swap(candidate);
    return true;
}

bool FactorBase::prime(PrimeIdeal& out, slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length() ||
        !same_order_parent(out.parent(), &parent_)) {
        return false;
    }
    return out.set(primes_[static_cast<std::size_t>(index)]);
}

const PrimeIdeal* FactorBase::prime_at(slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= length()) {
        return nullptr;
    }
    return &primes_[static_cast<std::size_t>(index)];
}

bool FactorBase::rational_prime_block_data(
        flint::FmpzRef rational_prime,
        slong& length,
        slong block_index) const noexcept {
    length = 0;
    if (!is_defined() || block_index < 0 ||
        block_index >= rational_prime_block_count()) {
        return false;
    }

    const PrimeBlock& block =
            blocks_[static_cast<std::size_t>(block_index)];
    fmpz_set(rational_prime.raw(), block.rational_prime.raw());
    length = block.length;
    return true;
}

bool FactorBase::rational_prime_block_index(slong& index,
                                            slong block_index,
                                            slong offset) const noexcept {
    index = -1;
    if (!is_defined() || block_index < 0 ||
        block_index >= rational_prime_block_count()) {
        return false;
    }

    const PrimeBlock& block =
            blocks_[static_cast<std::size_t>(block_index)];
    if (offset < 0 || offset >= static_cast<slong>(block.indices.size())) {
        return false;
    }

    index = block.indices[static_cast<std::size_t>(offset)];
    return true;
}

bool FactorBase::rational_prime_block(flint::FmpzRef rational_prime,
                                      slong& start,
                                      slong& length,
                                      slong block_index) const noexcept {
    start = 0;
    length = 0;
    if (!is_defined() || block_index < 0 ||
        block_index >= rational_prime_block_count()) {
        return false;
    }

    const PrimeBlock& block =
            blocks_[static_cast<std::size_t>(block_index)];
    if (!block.contiguous) {
        return false;
    }
    fmpz_set(rational_prime.raw(), block.rational_prime.raw());
    start = block.start;
    length = block.length;
    return true;
}

bool FactorBase::rational_prime_block_index_for_prime(
        slong& block_index,
        flint::FmpzConstRef rational_prime) const noexcept {
    block_index = -1;
    if (!is_defined()) {
        return false;
    }

    bool found = false;
    const slong position = block_position(rational_prime, found);
    if (!found) {
        return false;
    }
    block_index = position;
    return true;
}

bool FactorBase::equal(const FactorBase& other) const noexcept {
    if (!is_defined() || !other.is_defined() ||
        !same_order_parent(parent(), other.parent()) ||
        length() != other.length() ||
        rational_prime_block_count() != other.rational_prime_block_count()) {
        return false;
    }

    for (slong i = 0; i < rational_prime_block_count(); ++i) {
        const PrimeBlock& left = blocks_[static_cast<std::size_t>(i)];
        const PrimeBlock& right = other.blocks_[static_cast<std::size_t>(i)];
        if (!flint::fmpz_equal(left.rational_prime, right.rational_prime) ||
            left.start != right.start || left.length != right.length ||
            left.contiguous != right.contiguous ||
            left.indices != right.indices) {
            return false;
        }
    }

    for (slong i = 0; i < length(); ++i) {
        if (!primes_[static_cast<std::size_t>(i)].equal(
                    other.primes_[static_cast<std::size_t>(i)])) {
            return false;
        }
    }
    return true;
}

slong FactorBase::index(const PrimeIdeal& prime) const noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data()) {
        return -1;
    }

    flint::Fmpz rational_prime;
    if (!prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return -1;
    }

    bool found = false;
    const slong block =
            block_position(flint::FmpzConstRef(rational_prime), found);
    if (!found) {
        return -1;
    }

    const PrimeBlock& prime_block = blocks_[static_cast<std::size_t>(block)];
    for (slong index : prime_block.indices) {
        if (primes_[static_cast<std::size_t>(index)].equal(prime)) {
            return index;
        }
    }
    return -1;
}

bool FactorBase::append_prime(const PrimeIdeal& prime) noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data()) {
        return false;
    }

    PrimeIdeal copy(parent_);
    flint::Fmpz rational_prime;
    if (!copy.is_defined() || !copy.set(prime) ||
        !copy.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    const slong index = static_cast<slong>(primes_.size());
    primes_.emplace_back(std::move(copy));
    return add_block_index(flint::FmpzConstRef(rational_prime), index);
}

bool FactorBase::append_prime(PrimeIdeal&& prime) noexcept {
    if (!is_defined() || !same_order_parent(prime.parent(), &parent_) ||
        !prime.has_prime_data()) {
        return false;
    }

    flint::Fmpz rational_prime;
    if (!prime.rational_prime(flint::FmpzRef(rational_prime))) {
        return false;
    }

    const slong index = static_cast<slong>(primes_.size());
    primes_.emplace_back(std::move(prime));
    return add_block_index(flint::FmpzConstRef(rational_prime), index);
}

slong FactorBase::block_position(
        flint::FmpzConstRef rational_prime,
        bool& found) const noexcept {
    found = false;
    slong lo = 0;
    slong hi = static_cast<slong>(blocks_.size());
    while (lo < hi) {
        const slong mid = lo + (hi - lo) / 2;
        if (fmpz_cmp(blocks_[static_cast<std::size_t>(mid)]
                             .rational_prime.raw(),
                     rational_prime.raw()) < 0) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    if (lo < static_cast<slong>(blocks_.size()) &&
        fmpz_equal(blocks_[static_cast<std::size_t>(lo)]
                           .rational_prime.raw(),
                   rational_prime.raw()) != 0) {
        found = true;
    }
    return lo;
}

void FactorBase::refresh_block(PrimeBlock& block) noexcept {
    block.length = static_cast<slong>(block.indices.size());
    block.start = block.indices.empty() ? 0 : block.indices.front();
    block.contiguous = true;
    for (slong i = 0; i < static_cast<slong>(block.indices.size()); ++i) {
        if (block.indices[static_cast<std::size_t>(i)] != block.start + i) {
            block.contiguous = false;
            break;
        }
    }
}

bool FactorBase::add_block_index(flint::FmpzConstRef rational_prime,
                                 slong prime_index) noexcept {
    if (prime_index < 0 || prime_index >= static_cast<slong>(primes_.size())) {
        return false;
    }

    bool found = false;
    const slong pos = block_position(rational_prime, found);
    if (found) {
        PrimeBlock& block = blocks_[static_cast<std::size_t>(pos)];
        block.indices.push_back(prime_index);
        refresh_block(block);
        return true;
    }

    PrimeBlock block;
    fmpz_set(block.rational_prime.raw(), rational_prime.raw());
    block.indices.push_back(prime_index);
    refresh_block(block);
    blocks_.insert(blocks_.begin() + pos, std::move(block));
    return true;
}

bool FactorBase::rebuild_blocks() noexcept {
    std::vector<PrimeBlock> previous_blocks = std::move(blocks_);
    blocks_.clear();
    flint::Fmpz rational_prime;
    for (slong i = 0; i < static_cast<slong>(primes_.size()); ++i) {
        if (!primes_[static_cast<std::size_t>(i)].rational_prime(
                    flint::FmpzRef(rational_prime)) ||
            !add_block_index(flint::FmpzConstRef(rational_prime), i)) {
            return false;
        }
    }
    if (blocks_.size() != previous_blocks.size()) {
        return false;
    }
    for (std::size_t i = 0; i < blocks_.size(); ++i) {
        if (!flint::fmpz_equal(blocks_[i].rational_prime,
                               previous_blocks[i].rational_prime)) {
            return false;
        }
        blocks_[i].complete = previous_blocks[i].complete;
    }
    return true;
}

bool FactorBase::contains(const PrimeIdeal& prime) const noexcept {
    return index(prime) >= 0;
}

bool factor_base_class_group_bound(flint::FmpzRef out,
                                   const Order& order) noexcept {
    if (!order.is_maximal()) {
        return false;
    }

    const slong degree = order.degree();
    flint::Fmpz bound;
    if (degree == 1) {
        fmpz_one(bound.raw());
        fmpz_set(out.raw(), bound.raw());
        return true;
    }

    flint::Fmpz discriminant;
    flint::Fmpz abs_discriminant;
    if (!order.discriminant(flint::FmpzRef(discriminant)) ||
        fmpz_is_zero(discriminant.raw()) != 0 ||
        order.parent() == nullptr) {
        return false;
    }

    Signature sig;
    if (!sig.compute(*order.parent())) {
        return false;
    }

    fmpz_abs(abs_discriminant.raw(), discriminant.raw());
    if (degree == 2 && fmpz_sgn(discriminant.raw()) < 0) {
        floor_sqrt_div_ui(bound, flint::FmpzConstRef(abs_discriminant), 3);
    } else if (degree == 2) {
        fmpz_sqrt(bound.raw(), abs_discriminant.raw());
        fmpz_fdiv_q_2exp(bound.raw(), bound.raw(), 1);
    } else {
        generic_minkowski_bound(bound,
                                flint::FmpzConstRef(abs_discriminant),
                                degree,
                                sig.r2());
    }

    fmpz_set(out.raw(), bound.raw());
    return true;
}

}  // namespace silex
