#include <silex/factor_base.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/order.hpp>

#include "factor_base/factor_base_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <utility>

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
    cubic_order.set_maximality(true);
    assert(set_fmpz_si(bound, 99));
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), cubic_order));
    assert(sflint::fmpz_equal_si(bound, 5));

    return 0;
}

int test_failure_preserves_output() {
    silex::NumberField rational = degree_one_field();
    silex::Order rational_order = silex::test::equation_order(rational);

    silex::NumberField quadratic = quadratic_field(2);
    silex::Order nonmax_order = silex::test::equation_order(quadratic);
    nonmax_order.set_maximality(false);

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
    assert(test_class_group_bound() == 0);
    assert(test_failure_preserves_output() == 0);
    assert(test_keeps_parent_order_alive() == 0);
    return 0;
}
