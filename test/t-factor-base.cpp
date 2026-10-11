#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/order.hpp>

#include "factor_base/factor_base_internal.hpp"
#include "prime_ideal/prime_ideal_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

void poly_x(sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
}

void poly_x3_minus(sflint::FmpqPoly& polynomial, slong value) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -value);
}

void poly_x4_minus_x3_minus_x2_plus_2(
        sflint::FmpqPoly& polynomial) noexcept {
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 2);
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x(polynomial);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x3_minus(polynomial, 2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quartic_field() noexcept {
    sflint::FmpqPoly polynomial;
    poly_x4_minus_x3_minus_x2_plus_2(polynomial);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

slong count_prime_with_ef(const silex::FactorBase& base,
                          slong rational_prime,
                          slong ramification_index,
                          slong residue_degree) noexcept {
    const silex::Order* parent = base.parent();
    if (parent == nullptr) {
        return 0;
    }

    slong count = 0;
    silex::PrimeIdeal prime(*parent);
    sflint::Fmpz p;
    for (slong i = 0; i < base.length(); ++i) {
        if (!base.prime(prime, i) ||
            !prime.rational_prime(sflint::FmpzRef(p))) {
            return -1;
        }
        if (sflint::fmpz_equal_si(p, rational_prime) &&
            prime.ramification_index() == ramification_index &&
            prime.residue_degree() == residue_degree) {
            ++count;
        }
    }
    return count;
}

bool block_matches(const silex::FactorBase& base,
                   slong block_index,
                   slong rational_prime,
                   slong start,
                   slong length) noexcept {
    sflint::Fmpz p;
    slong block_start = -1;
    slong block_length = -1;
    return base.rational_prime_block(sflint::FmpzRef(p), block_start,
                                     block_length, block_index) &&
           sflint::fmpz_equal_si(p, rational_prime) &&
           block_start == start &&
           block_length == length;
}

bool block_completeness_matches_decomposition(
        const silex::FactorBase& base,
        slong block_index) noexcept {
    const silex::Order* parent = base.parent();
    sflint::Fmpz p;
    slong length = 0;
    bool complete = false;
    if (parent == nullptr ||
        !base.rational_prime_block_data(sflint::FmpzRef(p), length,
                                        block_index) ||
        !silex::detail::FactorBaseBlockAccess::
                rational_prime_block_is_complete(
                        complete, base, block_index)) {
        return false;
    }

    silex::PrimeIdealList decomposition;
    if (!silex::decompose_prime(decomposition, *parent,
                                sflint::FmpzConstRef(p))) {
        return false;
    }
    bool expected = decomposition.size() == length;
    for (slong i = 0; expected && i < decomposition.size(); ++i) {
        const silex::PrimeIdeal* prime = decomposition.at(i);
        expected = prime != nullptr && base.contains(*prime);
    }
    return complete == expected;
}

bool all_block_completeness_matches_decomposition(
        const silex::FactorBase& base) noexcept {
    for (slong i = 0; i < base.rational_prime_block_count(); ++i) {
        if (!block_completeness_matches_decomposition(base, i)) {
            return false;
        }
    }
    return true;
}

int test_degree_one() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::FactorBase base(order);
    silex::FactorBase copy(order);
    silex::FactorBase other(order);
    assert(base.is_defined());
    assert(silex::same_order_parent(base.parent(), &order));

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 7));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.rational_prime_blocks_are_complete());
    assert(base.length() == 4);
    assert(base.rational_prime_block_count() == 4);
    assert(block_matches(base, 0, 2, 0, 1));
    assert(block_matches(base, 1, 3, 1, 1));
    assert(block_matches(base, 2, 5, 2, 1));
    assert(block_matches(base, 3, 7, 3, 1));
    slong invalid_start = -1;
    slong invalid_length = -1;
    assert(!base.rational_prime_block(sflint::FmpzRef(bound), invalid_start,
                                      invalid_length, 4));
    assert(count_prime_with_ef(base, 2, 1, 1) == 1);
    assert(count_prime_with_ef(base, 3, 1, 1) == 1);
    assert(count_prime_with_ef(base, 5, 1, 1) == 1);
    assert(count_prime_with_ef(base, 7, 1, 1) == 1);

    silex::PrimeIdeal prime(order);
    sflint::Fmpz p;
    assert(base.prime(prime, 2));
    assert(prime.rational_prime(sflint::FmpzRef(p)));
    assert(sflint::fmpz_equal_si(p, 5));
    assert(base.contains(prime));
    assert(base.index(prime) == 2);
    slong block_index = -1;
    assert(base.rational_prime_block_index_for_prime(
            block_index, sflint::FmpzConstRef(p)));
    assert(block_index == 2);
    const silex::PrimeIdeal* borrowed = base.prime_at(2);
    assert(borrowed != nullptr);
    assert(borrowed->rational_prime(sflint::FmpzRef(p)));
    assert(sflint::fmpz_equal_si(p, 5));
    assert(base.prime_at(-1) == nullptr);
    assert(base.prime_at(base.length()) == nullptr);

    assert(set_fmpz_si(p, 11));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    silex::PrimeIdeal absent(order);
    assert(absent.set(*primes.at(0)));
    assert(!base.contains(absent));
    assert(base.index(absent) == -1);
    assert(!base.rational_prime_block_index_for_prime(
            block_index, sflint::FmpzConstRef(p)));
    assert(block_index == -1);

    assert(copy.set(base));
    assert(copy.rational_prime_blocks_are_complete());
    assert(copy.length() == 4);
    assert(copy.set(copy));
    assert(copy.length() == 4);
    swap(copy, other);
    assert(copy.length() == 0);
    assert(other.length() == 4);
    assert(other.rational_prime_blocks_are_complete());
    assert(count_prime_with_ef(other, 7, 1, 1) == 1);

    return 0;
}

int test_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 7));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 5);
    assert(base.rational_prime_block_count() == 4);
    assert(block_matches(base, 0, 2, 0, 1));
    assert(block_matches(base, 1, 3, 1, 1));
    assert(block_matches(base, 2, 5, 2, 1));
    assert(block_matches(base, 3, 7, 3, 2));
    assert(count_prime_with_ef(base, 2, 2, 1) == 1);
    assert(count_prime_with_ef(base, 3, 1, 2) == 1);
    assert(count_prime_with_ef(base, 5, 1, 2) == 1);
    assert(count_prime_with_ef(base, 7, 1, 1) == 2);

    silex::PrimeIdeal prime(order);
    for (slong i = 0; i < base.length(); ++i) {
        assert(base.prime(prime, i));
        assert(base.contains(prime));
        assert(base.index(prime) == i);
    }

    return 0;
}

int test_generic_cubic() {
    silex::NumberField field = cubic_field();
    silex::Order equation = silex::test::equation_order(field);

    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());

    silex::FactorBase base(maximal);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() >= 3);
    assert(count_prime_with_ef(base, 5, 1, 1) == 1);
    assert(count_prime_with_ef(base, 5, 1, 2) == 1);

    silex::PrimeIdeal prime(maximal);
    for (slong i = 0; i < base.length(); ++i) {
        assert(base.prime(prime, i));
        assert(base.contains(prime));
        assert(base.index(prime) == i);
    }

    return 0;
}

int test_prime_ideal_norm_bounded_shape() {
    silex::NumberField cubic = cubic_field();
    silex::Order cubic_equation = silex::test::equation_order(cubic);
    silex::Order cubic_maximal(cubic);
    assert(cubic_maximal.maximal_order(cubic_equation));
    assert(cubic_maximal.is_maximal());

    sflint::Fmpz cubic_bound;
    assert(set_fmpz_si(cubic_bound, 5));
    silex::FactorBase cubic_full(cubic_maximal);
    silex::FactorBase cubic_norm(cubic_maximal);
    assert(cubic_full.build(sflint::FmpzConstRef(cubic_bound)));
    assert(cubic_norm.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(cubic_bound)));
    assert(cubic_full.rational_prime_blocks_are_complete());
    assert(!cubic_norm.rational_prime_blocks_are_complete());
    assert(all_block_completeness_matches_decomposition(cubic_full));
    assert(all_block_completeness_matches_decomposition(cubic_norm));
    assert(cubic_full.length() == 4);
    assert(cubic_norm.length() == 3);
    assert(cubic_norm.rational_prime_block_count() == 3);
    assert(count_prime_with_ef(cubic_norm, 5, 1, 1) == 1);
    assert(count_prime_with_ef(cubic_norm, 5, 1, 2) == 0);

    silex::NumberField quartic = quartic_field();
    silex::Order quartic_equation = silex::test::equation_order(quartic);
    silex::Order quartic_maximal(quartic);
    assert(quartic_maximal.maximal_order(quartic_equation));
    assert(quartic_maximal.is_maximal());

    sflint::Fmpz quartic_bound;
    assert(set_fmpz_si(quartic_bound, 12));
    silex::FactorBase quartic_full(quartic_maximal);
    silex::FactorBase quartic_norm(quartic_maximal);
    assert(quartic_full.build(sflint::FmpzConstRef(quartic_bound)));
    assert(quartic_norm.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(quartic_bound)));
    assert(quartic_full.rational_prime_blocks_are_complete());
    assert(!quartic_norm.rational_prime_blocks_are_complete());
    assert(all_block_completeness_matches_decomposition(quartic_full));
    assert(all_block_completeness_matches_decomposition(quartic_norm));
    assert(quartic_full.length() == 8);
    assert(quartic_norm.length() == 4);
    assert(quartic_norm.rational_prime_block_count() == 3);
    assert(block_matches(quartic_norm, 0, 2, 0, 2));
    assert(block_matches(quartic_norm, 1, 3, 2, 1));
    assert(block_matches(quartic_norm, 2, 11, 3, 1));

    return 0;
}

int test_lll_relation_base_order() {
    silex::NumberField quartic = quartic_field();
    silex::Order quartic_equation = silex::test::equation_order(quartic);
    silex::Order quartic_maximal(quartic);
    assert(quartic_maximal.maximal_order(quartic_equation));
    assert(quartic_maximal.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 12));
    silex::FactorBase full(quartic_maximal);
    silex::FactorBase reference(quartic_maximal);
    assert(full.build(sflint::FmpzConstRef(bound)));
    assert(reference.build_lll_relation_base(sflint::FmpzConstRef(bound)));
    assert(full.rational_prime_blocks_are_complete());
    assert(reference.rational_prime_blocks_are_complete());
    assert(reference.length() == full.length());

    silex::PrimeIdeal previous(quartic_maximal);
    silex::PrimeIdeal current(quartic_maximal);
    sflint::Fmpz previous_norm;
    sflint::Fmpz current_norm;
    for (slong i = 0; i < reference.length(); ++i) {
        assert(reference.prime(current, i));
        assert(reference.contains(current));
        assert(reference.index(current) == i);
        assert(full.contains(current));
        assert(current.norm(sflint::FmpzRef(current_norm)));
        if (i > 0) {
            assert(previous.norm(sflint::FmpzRef(previous_norm)));
            assert(sflint::fmpz_cmp(sflint::FmpzConstRef(previous_norm),
                                    sflint::FmpzConstRef(current_norm)) >= 0);
        }
        assert(previous.set(current));
    }

    for (slong block = 0; block < reference.rational_prime_block_count();
         ++block) {
        sflint::Fmpz rational_prime;
        slong length = 0;
        assert(reference.rational_prime_block_data(
                sflint::FmpzRef(rational_prime), length, block));
        assert(length > 0);
        for (slong offset = 0; offset < length; ++offset) {
            slong index = -1;
            sflint::Fmpz p;
            assert(reference.rational_prime_block_index(index, block, offset));
            assert(index >= 0 && index < reference.length());
            assert(reference.prime(current, index));
            assert(current.rational_prime(sflint::FmpzRef(p)));
            assert(sflint::fmpz_equal(sflint::FmpzConstRef(p),
                                      sflint::FmpzConstRef(rational_prime)));
        }
    }

    return 0;
}

// The norm ordering applied to a prebuilt norm-bounded base (the grh
// policy's probe base) keeps the same prime ideals, sorts them by norm,
// largest first, with ties in build order, and keeps each rational-prime
// block and its completeness.
int test_lll_relation_base_from_norm_bounded() {
    silex::NumberField quartic = quartic_field();
    silex::Order quartic_equation = silex::test::equation_order(quartic);
    silex::Order quartic_maximal(quartic);
    assert(quartic_maximal.maximal_order(quartic_equation));

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 40));
    silex::FactorBase built(quartic_maximal);
    silex::FactorBase ordered(quartic_maximal);
    assert(built.build_prime_ideal_norm_bounded(sflint::FmpzConstRef(bound)));
    assert(ordered.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(silex::detail::FactorBaseBlockAccess::
                   order_lll_relation_base_by_norm(ordered));
    assert(ordered.length() == built.length());
    assert(ordered.rational_prime_block_count() ==
           built.rational_prime_block_count());
    assert(all_block_completeness_matches_decomposition(ordered));

    silex::PrimeIdeal current(quartic_maximal);
    sflint::Fmpz previous_norm;
    sflint::Fmpz current_norm;
    slong previous_index = -1;
    bool has_tie = false;
    for (slong i = 0; i < ordered.length(); ++i) {
        assert(ordered.prime(current, i));
        assert(ordered.index(current) == i);
        const slong built_index = built.index(current);
        assert(built_index >= 0);
        assert(current.norm(sflint::FmpzRef(current_norm)));
        if (i > 0) {
            const int cmp = sflint::fmpz_cmp(
                    sflint::FmpzConstRef(previous_norm),
                    sflint::FmpzConstRef(current_norm));
            assert(cmp >= 0);
            if (cmp == 0) {
                has_tie = true;
                assert(built_index > previous_index);
            }
        }
        sflint::fmpz_set(sflint::FmpzRef(previous_norm),
                         sflint::FmpzConstRef(current_norm));
        previous_index = built_index;
    }
    assert(has_tie);

    for (slong block = 0; block < ordered.rational_prime_block_count();
         ++block) {
        sflint::Fmpz ordered_prime;
        sflint::Fmpz built_prime;
        slong ordered_length = 0;
        slong built_length = 0;
        bool ordered_complete = false;
        bool built_complete = false;
        assert(ordered.rational_prime_block_data(
                sflint::FmpzRef(ordered_prime), ordered_length, block));
        assert(built.rational_prime_block_data(
                sflint::FmpzRef(built_prime), built_length, block));
        assert(sflint::fmpz_equal(sflint::FmpzConstRef(ordered_prime),
                                  sflint::FmpzConstRef(built_prime)));
        assert(ordered_length == built_length);
        assert(silex::detail::FactorBaseBlockAccess::
                       rational_prime_block_is_complete(
                               ordered_complete, ordered, block));
        assert(silex::detail::FactorBaseBlockAccess::
                       rational_prime_block_is_complete(
                               built_complete, built, block));
        assert(ordered_complete == built_complete);
    }

    return 0;
}

bool block_length_for_prime(slong& length,
                            const silex::FactorBase& base,
                            slong rational_prime) noexcept {
    sflint::Fmpz p;
    sflint::fmpz_set_si(sflint::FmpzRef(p), rational_prime);
    slong block = -1;
    slong start = 0;
    if (!base.rational_prime_block_index_for_prime(
                block, sflint::FmpzConstRef(p))) {
        length = 0;
        return true;
    }
    return base.rational_prime_block(sflint::FmpzRef(p), start, length,
                                     block);
}

// A norm-bounded build that takes its decompositions from a cache filled by
// earlier builds of the same order equals a fresh build at every bound,
// including a bound below the cache's reach and a bound at which a cached
// rational prime gains a prime ideal of larger norm (5 in Q(2^(1/3)): a
// degree-one and a degree-two prime, norms 5 and 25).  A cache from another
// order fails the build and leaves the base unchanged.
int test_norm_bounded_build_with_decomposition_cache() {
    silex::NumberField cubic = cubic_field();
    silex::Order cubic_equation = silex::test::equation_order(cubic);
    silex::Order cubic_maximal(cubic);
    assert(cubic_maximal.maximal_order(cubic_equation));

    silex::detail::PrimeDecompositionCache cache;
    const slong bounds[] = {4, 8, 16, 32, 64, 12};
    const std::size_t primes_up_to[] = {2, 4, 6, 11, 18, 18};
    slong previous_five_length = -1;
    for (std::size_t i = 0; i < sizeof(bounds) / sizeof(bounds[0]); ++i) {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, bounds[i]));
        silex::FactorBase cached(cubic_maximal);
        silex::FactorBase fresh(cubic_maximal);
        assert(silex::detail::FactorBaseBlockAccess::
                       build_prime_ideal_norm_bounded_cached(
                               cached, sflint::FmpzConstRef(bound), cache));
        assert(fresh.build_prime_ideal_norm_bounded(
                sflint::FmpzConstRef(bound)));
        assert(cache.decompositions.size() == primes_up_to[i]);
        assert(cached.equal(fresh));
        assert(fresh.equal(cached));
        assert(all_block_completeness_matches_decomposition(cached));
        for (slong block = 0; block < cached.rational_prime_block_count();
             ++block) {
            bool cached_complete = false;
            bool fresh_complete = false;
            assert(silex::detail::FactorBaseBlockAccess::
                           rational_prime_block_is_complete(
                                   cached_complete, cached, block));
            assert(silex::detail::FactorBaseBlockAccess::
                           rational_prime_block_is_complete(
                                   fresh_complete, fresh, block));
            assert(cached_complete == fresh_complete);
        }
        slong five_length = 0;
        assert(block_length_for_prime(five_length, cached, 5));
        if (bounds[i] == 16) {
            assert(five_length == 1);
        }
        if (bounds[i] == 32) {
            assert(previous_five_length == 1 && five_length == 2);
        }
        previous_five_length = five_length;
    }

    silex::NumberField quadratic = quadratic_field(2);
    silex::Order quadratic_order = silex::test::equation_order(quadratic);
    assert(quadratic_order.is_maximal());
    silex::detail::PrimeDecompositionCache other_cache;
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 16));
    silex::FactorBase quadratic_base(quadratic_order);
    assert(silex::detail::FactorBaseBlockAccess::
                   build_prime_ideal_norm_bounded_cached(
                           quadratic_base, sflint::FmpzConstRef(bound),
                           other_cache));
    silex::FactorBase kept(cubic_maximal);
    assert(kept.build_prime_ideal_norm_bounded(sflint::FmpzConstRef(bound)));
    silex::FactorBase copy(cubic_maximal);
    assert(copy.set(kept));
    assert(!silex::detail::FactorBaseBlockAccess::
                    build_prime_ideal_norm_bounded_cached(
                            kept, sflint::FmpzConstRef(bound), other_cache));
    assert(kept.equal(copy));

    return 0;
}

silex::NumberField monic_field(const slong* coefficients,
                               slong degree) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i <= degree; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }
    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

std::size_t skipped_entry_count(
        const silex::detail::PrimeDecompositionCache& cache) noexcept {
    std::size_t count = 0;
    for (const slong degree : cache.skipped_smallest_residue_degree) {
        if (degree > 0) {
            ++count;
        }
    }
    return count;
}

// Builds over one cache at the doubling bounds 2, 4, ..., 2^steps, as the
// grh factor-base probe does, and checks each against a fresh build without
// a cache: the same prime ideals in the same order, the same blocks and the
// same completeness flags.  Returns the number of placeholders after each
// build in skipped[0..steps-1] and whether any placeholder was later filled.
bool doubling_builds_match_fresh(const silex::Order& order,
                                 slong steps,
                                 std::size_t* skipped,
                                 bool& filled_placeholder) noexcept {
    silex::detail::PrimeDecompositionCache cache;
    filled_placeholder = false;
    std::vector<slong> previous;
    for (slong step = 0; step < steps; ++step) {
        sflint::Fmpz bound;
        sflint::fmpz_one(sflint::FmpzRef(bound));
        sflint::fmpz_mul_2exp(sflint::FmpzRef(bound),
                              sflint::FmpzConstRef(bound),
                              static_cast<ulong>(step + 1));
        silex::FactorBase cached(order);
        silex::FactorBase fresh(order);
        if (!silex::detail::FactorBaseBlockAccess::
                    build_prime_ideal_norm_bounded_cached(
                            cached, sflint::FmpzConstRef(bound), cache) ||
            !fresh.build_prime_ideal_norm_bounded(
                    sflint::FmpzConstRef(bound)) ||
            !cached.equal(fresh) || !fresh.equal(cached) ||
            cached.rational_prime_block_count() !=
                    fresh.rational_prime_block_count() ||
            cache.decompositions.size() !=
                    cache.skipped_smallest_residue_degree.size()) {
            return false;
        }
        for (slong block = 0; block < cached.rational_prime_block_count();
             ++block) {
            sflint::Fmpz cached_prime;
            sflint::Fmpz fresh_prime;
            slong cached_length = 0;
            slong fresh_length = 0;
            bool cached_complete = false;
            bool fresh_complete = false;
            if (!cached.rational_prime_block_data(
                        sflint::FmpzRef(cached_prime), cached_length,
                        block) ||
                !fresh.rational_prime_block_data(
                        sflint::FmpzRef(fresh_prime), fresh_length, block) ||
                !sflint::fmpz_equal(sflint::FmpzConstRef(cached_prime),
                                    sflint::FmpzConstRef(fresh_prime)) ||
                cached_length != fresh_length ||
                !silex::detail::FactorBaseBlockAccess::
                        rational_prime_block_is_complete(
                                cached_complete, cached, block) ||
                !silex::detail::FactorBaseBlockAccess::
                        rational_prime_block_is_complete(
                                fresh_complete, fresh, block) ||
                cached_complete != fresh_complete) {
                return false;
            }
        }
        for (std::size_t i = 0; i < previous.size(); ++i) {
            if (previous[i] > 0 &&
                cache.skipped_smallest_residue_degree[i] == 0) {
                filled_placeholder = true;
            }
        }
        previous = cache.skipped_smallest_residue_degree;
        skipped[step] = skipped_entry_count(cache);
    }
    return true;
}

// A cached norm-bounded build skips, without decomposing, a rational prime
// whose residue degrees (read from the defining polynomial) put every prime
// ideal above it over the bound, and decomposes it once a doubled bound
// admits one; the base is the fresh build's at every bound.  Index
// divisors and discriminant divisors are always decomposed.
int test_norm_bounded_build_skips_primes_over_the_bound() {
    // Q(2^(1/3)): 7 and 13 are inert (norms 343 and 2197); 2 and 3 are
    // ramified.  Bounds 2, 4, ..., 4096.
    {
        silex::NumberField cubic = cubic_field();
        silex::Order equation = silex::test::equation_order(cubic);
        silex::Order maximal(cubic);
        assert(maximal.maximal_order(equation));
        std::size_t skipped[12] = {};
        bool filled = false;
        assert(doubling_builds_match_fresh(maximal, 12, skipped, filled));
        assert(filled);
        bool any_skipped = false;
        for (const std::size_t count : skipped) {
            any_skipped = any_skipped || count > 0;
        }
        assert(any_skipped);
    }

    // Q(sqrt(5)) from x^2 - 5: 2 divides the index of Z[sqrt(5)] and is
    // inert (norm 4 > 2), so it is decomposed even at bound 2; 5 divides
    // the discriminant.  3 and 7 are inert (norms 9 and 49): bound 4 skips
    // 3, bound 8 skips 3 and 7, and bound 16 decomposes 3.
    {
        const slong coefficients[] = {-5, 0, 1};
        silex::NumberField quadratic = monic_field(coefficients, 2);
        silex::Order equation = silex::test::equation_order(quadratic);
        silex::Order maximal(quadratic);
        assert(maximal.maximal_order(equation));
        assert(!equation.is_maximal());

        silex::detail::PrimeDecompositionCache cache;
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        silex::FactorBase base(maximal);
        assert(silex::detail::FactorBaseBlockAccess::
                       build_prime_ideal_norm_bounded_cached(
                               base, sflint::FmpzConstRef(bound), cache));
        assert(base.length() == 0);
        assert(cache.decompositions.size() == 1);
        assert(cache.skipped_smallest_residue_degree[0] == 0);
        assert(cache.decompositions[0].size() == 1);

        std::size_t skipped[5] = {};
        bool filled = false;
        assert(doubling_builds_match_fresh(maximal, 5, skipped, filled));
        assert(filled);
        assert(skipped[0] == 0 && skipped[1] == 1 && skipped[2] == 2);
    }

    // Dedekind's cubic x^3 - x^2 - 2x - 8: 2 splits completely in the
    // maximal order, so it divides the index of every Z[theta].
    {
        const slong coefficients[] = {-8, -2, -1, 1};
        silex::NumberField cubic = monic_field(coefficients, 3);
        silex::Order equation = silex::test::equation_order(cubic);
        silex::Order maximal(cubic);
        assert(maximal.maximal_order(equation));
        assert(!equation.is_maximal());
        std::size_t skipped[10] = {};
        bool filled = false;
        assert(doubling_builds_match_fresh(maximal, 10, skipped, filled));
    }

    return 0;
}

// equal() answers copies from a shared identity and everything else from the
// element-wise comparison; both must give the same answer.
int test_equal_copies_and_independent_builds() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 7));
    sflint::Fmpz smaller_bound;
    assert(set_fmpz_si(smaller_bound, 5));

    silex::FactorBase undefined;
    assert(!undefined.equal(undefined));

    silex::FactorBase base(order);
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.equal(base));
    assert(!base.equal(undefined));
    assert(!undefined.equal(base));

    silex::FactorBase copy;
    assert(copy.set(base));
    assert(copy.equal(base));
    assert(base.equal(copy));

    // Built separately with the same contents: equal through the deep path.
    silex::FactorBase independent(order);
    assert(independent.build(sflint::FmpzConstRef(bound)));
    assert(independent.equal(base));
    assert(base.equal(independent));

    // A copy of a copy, and copies moved or swapped, keep comparing equal.
    silex::FactorBase second;
    assert(second.set(copy));
    assert(second.equal(base));
    silex::FactorBase moved(std::move(second));
    assert(moved.equal(base));
    assert(!second.equal(base));
    silex::FactorBase other(order);
    assert(other.build(sflint::FmpzConstRef(smaller_bound)));
    assert(!other.equal(base));
    other.swap(moved);
    assert(other.equal(base));
    assert(!moved.equal(base));

    // Changing either side of a copy ends the equality.
    assert(copy.build(sflint::FmpzConstRef(smaller_bound)));
    assert(!copy.equal(base));
    assert(!base.equal(copy));
    assert(copy.equal(moved));
    assert(copy.set(base));
    assert(copy.equal(base));

    // Changing the source after the copy ends the equality both ways.
    {
        silex::FactorBase source(order);
        assert(source.build(sflint::FmpzConstRef(bound)));
        silex::FactorBase source_copy;
        assert(source_copy.set(source));
        assert(source_copy.equal(source));
        assert(source.build(sflint::FmpzConstRef(smaller_bound)));
        assert(!source_copy.equal(source));
        assert(!source.equal(source_copy));
    }

    assert(copy.define(order));
    assert(copy.length() == 0);
    assert(!copy.equal(base));
    silex::FactorBase empty(order);
    assert(copy.equal(empty));
    assert(copy.set(base));
    copy.clear();
    assert(!copy.equal(base));
    assert(!copy.equal(copy));

    // Reordering a copy in place changes its blocks and order.
    silex::NumberField quartic = quartic_field();
    silex::Order quartic_equation = silex::test::equation_order(quartic);
    silex::Order quartic_maximal(quartic);
    assert(quartic_maximal.maximal_order(quartic_equation));
    sflint::Fmpz quartic_bound;
    assert(set_fmpz_si(quartic_bound, 40));
    silex::FactorBase built(quartic_maximal);
    assert(built.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(quartic_bound)));
    silex::FactorBase ordered;
    assert(ordered.set(built));
    assert(ordered.equal(built));
    assert(silex::detail::FactorBaseBlockAccess::
                   order_lll_relation_base_by_norm(ordered));
    assert(!ordered.equal(built));
    assert(!built.equal(ordered));
    silex::FactorBase ordered_independent(quartic_maximal);
    assert(ordered_independent.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(quartic_bound)));
    assert(silex::detail::FactorBaseBlockAccess::
                   order_lll_relation_base_by_norm(ordered_independent));
    assert(ordered.equal(ordered_independent));

    return 0;
}

int test_class_group_bound() {
    sflint::Fmpz bound;

    silex::NumberField rational = degree_one_field();
    silex::Order rational_order = silex::test::equation_order(rational);
    assert(set_fmpz_si(bound, 99));
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), rational_order));
    assert(sflint::fmpz_is_one(bound));

    silex::NumberField real_quad = quadratic_field(17);
    silex::Order real_equation;
    silex::Order real_maximal(real_quad);
    real_equation = silex::test::equation_order(real_quad);
    assert(real_maximal.maximal_order(real_equation));
    assert(set_fmpz_si(bound, 99));
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), real_maximal));
    assert(sflint::fmpz_equal_si(bound, 2));

    silex::NumberField imag_quad = quadratic_field(-5);
    silex::Order imag_equation;
    silex::Order imag_maximal(imag_quad);
    imag_equation = silex::test::equation_order(imag_quad);
    assert(imag_maximal.maximal_order(imag_equation));
    assert(set_fmpz_si(bound, 99));
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), imag_maximal));
    assert(sflint::fmpz_equal_si(bound, 2));

    silex::NumberField imag_large = quadratic_field(-47);
    silex::Order imag_large_equation;
    silex::Order imag_large_maximal(imag_large);
    imag_large_equation = silex::test::equation_order(imag_large);
    assert(imag_large_maximal.maximal_order(imag_large_equation));
    assert(set_fmpz_si(bound, 99));
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), imag_large_maximal));
    assert(sflint::fmpz_equal_si(bound, 3));

    silex::NumberField nonmax_field = quadratic_field(5);
    silex::Order nonmax_order = silex::test::equation_order(nonmax_field);
    assert(!nonmax_order.is_maximal());
    assert(set_fmpz_si(bound, 99));
    assert(!silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), nonmax_order));
    assert(sflint::fmpz_equal_si(bound, 99));

    silex::NumberField cubic = cubic_field();
    silex::Order cubic_order = silex::test::equation_order(cubic);
    cubic_order = silex::test::verified_maximal_order(cubic_order);
    assert(set_fmpz_si(bound, 99));
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), cubic_order));
    assert(sflint::fmpz_equal_si(bound, 3));

    return 0;
}

int test_failure_preserves_output() {
    silex::NumberField rational = degree_one_field();
    silex::Order rational_order = silex::test::equation_order(rational);

    silex::NumberField quadratic = quadratic_field(2);
    sflint::Fmpz conductor;
    assert(set_fmpz_si(conductor, 2));
    silex::Order nonmax_order = silex::Order::quadratic_order(
            quadratic, sflint::FmpzConstRef(conductor));
    assert(nonmax_order.maximality_known() && !nonmax_order.is_maximal());

    silex::FactorBase base(rational_order);
    silex::FactorBase nonmax_base(nonmax_order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 3);

    assert(set_fmpz_si(bound, 1));
    assert(!base.build(sflint::FmpzConstRef(bound)));
    assert(base.length() == 3);
    assert(count_prime_with_ef(base, 5, 1, 1) == 1);
    assert(!base.build_prime_ideal_norm_bounded(sflint::FmpzConstRef(bound)));
    assert(base.length() == 3);
    assert(count_prime_with_ef(base, 5, 1, 1) == 1);

    assert(set_fmpz_si(bound, 5));
    assert(!nonmax_base.build(sflint::FmpzConstRef(bound)));
    assert(nonmax_base.length() == 0);
    assert(!nonmax_base.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(nonmax_base.length() == 0);

    return 0;
}

// Each prime of `base` is found at its own index, so no prime ideal occurs
// twice (index compares ideals, not their residue data).
bool primes_are_distinct(const silex::FactorBase& base) noexcept {
    for (slong i = 0; i < base.length(); ++i) {
        const silex::PrimeIdeal* prime = base.prime_at(i);
        if (prime == nullptr || base.index(*prime) != i) {
            return false;
        }
    }
    return true;
}

// Copies `prime`, a degree-one prime with residue polynomial x - r in the
// alpha convention, into `out` with the residue polynomial x - r - p: the
// same residue map onto F_p, so the same prime ideal, recorded with an
// unreduced coefficient.  PrimeIdeal::equal tells the two apart.
bool set_unreduced_residue_variant(silex::PrimeIdeal& out,
                                   const silex::PrimeIdeal& prime) noexcept {
    using silex::detail::PrimeIdealAccess;
    sflint::FmpzPoly residue;
    sflint::Fmpz p;
    sflint::Fmpz constant;
    if (prime.residue_degree() != 1 ||
        PrimeIdealAccess::residue_uses_integral_generator(prime) ||
        !prime.residue_polynomial(sflint::FmpzPolyRef(residue)) ||
        !prime.rational_prime(sflint::FmpzRef(p)) || !out.set(prime)) {
        return false;
    }
    fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
    fmpz_add(constant.raw(), constant.raw(), p.raw());
    fmpz_poly_set_coeff_fmpz(residue.raw(), 0, constant.raw());
    return PrimeIdealAccess::set_alpha_residue_polynomial_for_testing(
                   out, sflint::FmpzPolyConstRef(residue)) &&
           !out.equal(prime);
}

// x^3 - 2: the base at 5 holds degree-one primes above 2, 3 and 5 with
// residue data in the alpha convention.  A copy of each with a different
// residue polynomial is the same prime and is found at the same index.
int test_contains_ignores_residue_representation() {
    silex::NumberField field = cubic_field();
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::FactorBase base(maximal);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));

    slong checked = 0;
    for (slong i = 0; i < base.length(); ++i) {
        const silex::PrimeIdeal* prime = base.prime_at(i);
        assert(prime != nullptr);
        silex::PrimeIdeal variant(maximal);
        if (!set_unreduced_residue_variant(variant, *prime)) {
            continue;
        }
        assert(base.contains(variant));
        assert(base.index(variant) == i);
        ++checked;
    }
    assert(checked >= 1);
    return 0;
}

// Every builder lists distinct primes, and append_prime refuses to repeat
// one, leaving the base unchanged.  A repeated prime would make a block
// look complete (the sum of e f reaching the degree) while a prime above p
// is missing.
int test_append_rejects_repeated_prime() {
    using silex::detail::FactorBaseBlockAccess;
    silex::NumberField field = cubic_field();
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 30));
    silex::FactorBase complete(maximal);
    silex::FactorBase norm_bounded(maximal);
    silex::FactorBase relation(maximal);
    silex::FactorBase lll(maximal);
    assert(complete.build(sflint::FmpzConstRef(bound)));
    assert(norm_bounded.build_prime_ideal_norm_bounded(
            sflint::FmpzConstRef(bound)));
    assert(relation.build_relation_completion_base(
            sflint::FmpzConstRef(bound)));
    assert(lll.build_lll_relation_base(sflint::FmpzConstRef(bound)));
    assert(primes_are_distinct(complete));
    assert(primes_are_distinct(norm_bounded));
    assert(primes_are_distinct(relation));
    assert(primes_are_distinct(lll));

    silex::FactorBase copy;
    assert(copy.set(complete));
    assert(primes_are_distinct(copy));

    // build(bound, diagnostics) and the maximal imaginary-quadratic route
    // (built from set_first_degree_one_prime).
    silex::FactorBase with_diagnostics(maximal);
    assert(with_diagnostics.build(sflint::FmpzConstRef(bound), nullptr));
    assert(primes_are_distinct(with_diagnostics));
    for (slong radicand : {-5, -47}) {
        silex::NumberField imag = quadratic_field(radicand);
        silex::Order imag_equation = silex::test::equation_order(imag);
        silex::Order imag_maximal(imag);
        assert(imag_maximal.maximal_order(imag_equation));
        silex::FactorBase imag_base(imag_maximal);
        assert(FactorBaseBlockAccess::
                       build_maximal_imaginary_quadratic_for_testing(
                               imag_base, sflint::FmpzConstRef(bound)));
        assert(imag_base.length() > 0);
        assert(primes_are_distinct(imag_base));
    }

    // 5 = P Q with f(P) = 1 and f(Q) = 2.  A second copy of the first prime
    // is refused and the other prime above 5 is still accepted.
    silex::PrimeIdealList five;
    sflint::Fmpz p;
    assert(set_fmpz_si(p, 5));
    assert(silex::decompose_prime(five, maximal, sflint::FmpzConstRef(p)));
    assert(five.size() == 2);
    silex::FactorBase pair(maximal);
    assert(FactorBaseBlockAccess::append_prime_for_testing(pair,
                                                           *five.at(0)));
    assert(!FactorBaseBlockAccess::append_prime_for_testing(pair,
                                                            *five.at(0)));
    assert(pair.length() == 1 && pair.rational_prime_block_count() == 1);
    assert(FactorBaseBlockAccess::append_prime_for_testing(pair,
                                                           *five.at(1)));
    assert(pair.length() == 2 && primes_are_distinct(pair));

    // Repeats of every base prime, also with other residue data, are
    // refused and the contents stay those of an independent build.
    silex::FactorBase independent(maximal);
    assert(independent.build(sflint::FmpzConstRef(bound)));
    const slong length = complete.length();
    slong variants = 0;
    for (slong i = 0; i < length; ++i) {
        const silex::PrimeIdeal* prime = complete.prime_at(i);
        assert(prime != nullptr);
        silex::PrimeIdeal repeat(maximal);
        assert(repeat.set(*prime));
        assert(!FactorBaseBlockAccess::append_prime_for_testing(complete,
                                                                repeat));
        silex::PrimeIdeal variant(maximal);
        if (set_unreduced_residue_variant(variant, *prime)) {
            assert(!FactorBaseBlockAccess::append_prime_for_testing(
                    complete, variant));
            ++variants;
        }
        assert(complete.length() == length);
    }
    assert(variants >= 1);
    assert(complete.equal(independent));
    assert(primes_are_distinct(complete));
    return 0;
}

silex::FactorBase local_factor_base() {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);

    silex::FactorBase base(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    assert(base.build(sflint::FmpzConstRef(bound)));
    return base;
}

int test_keeps_parent_order_alive() {
    silex::FactorBase base = local_factor_base();
    assert(base.is_defined());
    assert(base.parent() != nullptr);
    assert(base.parent()->parent() != nullptr);
    assert(base.length() == 3);

    silex::PrimeIdeal prime(*base.parent());
    sflint::Fmpz p;
    assert(base.prime(prime, 2));
    assert(prime.rational_prime(sflint::FmpzRef(p)));
    assert(sflint::fmpz_equal_si(p, 5));
    assert(silex::same_order_parent(prime.parent(), base.parent()));
    return 0;
}

}  // namespace

int main() {
    assert(test_degree_one() == 0);
    assert(test_quadratic() == 0);
    assert(test_generic_cubic() == 0);
    assert(test_prime_ideal_norm_bounded_shape() == 0);
    assert(test_lll_relation_base_order() == 0);
    assert(test_lll_relation_base_from_norm_bounded() == 0);
    assert(test_norm_bounded_build_with_decomposition_cache() == 0);
    assert(test_norm_bounded_build_skips_primes_over_the_bound() == 0);
    assert(test_equal_copies_and_independent_builds() == 0);
    assert(test_class_group_bound() == 0);
    assert(test_failure_preserves_output() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    assert(test_contains_ignores_residue_representation() == 0);
    assert(test_append_rejects_repeated_prime() == 0);
    return 0;
}
