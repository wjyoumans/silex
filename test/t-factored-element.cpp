#include <silex/archimedean.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_vec.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField cubic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

void set_rational(silex::Element& element,
                  slong numerator,
                  ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, numerator, denominator);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool element_is_rational(const silex::Element& element,
                         slong numerator,
                         ulong denominator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq expected;
    sflint::Fmpq actual;
    sflint::fmpq_set_si(expected, numerator, denominator);
    assert(element.get_fmpq_poly(sflint::FmpqPolyRef(polynomial)));
    sflint::fmpq_poly_get_coeff_fmpq(sflint::FmpqRef(actual), polynomial, 0);
    return sflint::fmpq_poly_degree(polynomial) <= 0 &&
           sflint::fmpq_equal(actual, expected);
}

void set_one_plus_theta(silex::Element& element) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

bool has_radius_at_most_2exp(const arb_struct* value,
                             slong exponent) noexcept {
    return ::arb_is_finite(value) != 0 &&
           ::mag_cmp_2exp_si(arb_radref(value), exponent) <= 0;
}

void set_cubic_element(silex::Element& element,
                       slong constant,
                       slong linear,
                       slong quadratic) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, linear);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, quadratic);
    assert(element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)));
}

// The balanced product used by evaluation must give exactly the left-to-right
// product of the factor powers, for small, large, and negative exponents and
// for odd and even factor counts.
void check_many_factor_evaluation() noexcept {
    silex::NumberField field = cubic_field();
    const slong exponents[] = {1,  -1, 2,  -3, 5, 9,   -9, 10,
                               -11, 23, -24, 1, 3, -2, 7,  1,
                               -1, 12, 4,  -5, 6, 1,   -17, 1025};
    constexpr slong count = static_cast<slong>(sizeof(exponents) /
                                               sizeof(exponents[0]));

    for (slong length = 0; length <= count; ++length) {
        silex::FactoredElement product(field);
        silex::Element expected(field);
        silex::Element factor(field);
        silex::Element power(field);
        silex::Element next(field);
        silex::Element evaluated(field);
        assert(expected.one());
        for (slong i = 0; i < length; ++i) {
            set_cubic_element(factor, 2 * i + 3, i % 4 - 1, (i % 3) - 1);
            assert(product.push(factor, exponents[i]));
            sflint::Fmpz exponent;
            fmpz_set_si(exponent.raw(), exponents[i]);
            assert(power.pow_fmpz(factor, sflint::FmpzConstRef(exponent)));
            assert(next.multiply(expected, power));
            assert(expected.set(next));
        }
        assert(product.evaluate(evaluated));
        assert(evaluated.equal(expected));
    }
}

// Power tests whose reduced leftover has several factors, so the leftover is
// expanded by a balanced product rather than a single power.
void check_multi_factor_leftover_roots() noexcept {
    silex::NumberField field = cubic_field();
    silex::Element base[6] = {
            silex::Element(field), silex::Element(field),
            silex::Element(field), silex::Element(field),
            silex::Element(field), silex::Element(field)};
    for (slong i = 0; i < 6; ++i) {
        set_cubic_element(base[i], i + 2, 1, i % 2);
    }
    silex::Element partial(field);
    silex::Element next(field);
    assert(partial.one());
    for (slong i = 0; i < 4; ++i) {
        assert(next.multiply(partial, base[i]));
        assert(partial.set(next));
    }
    silex::Element partial_square(field);
    assert(partial_square.multiply(partial, partial));

    bool is_power = true;
    silex::FactoredElement input(field);
    silex::FactoredElement root(field);
    silex::FactoredElement check(field);
    silex::Element input_value(field);
    silex::Element check_value(field);

    // Cubes: b0 b1 b2 b3 (b0 b1 b2 b3)^2 b4^6 b5^-6.
    for (slong i = 0; i < 4; ++i) {
        assert(input.push(base[i], 1));
    }
    assert(input.push(partial_square, 1));
    assert(input.push(base[4], 6));
    assert(input.push(base[5], -6));
    assert(input.evaluate(input_value));
    assert(input.is_power_si(is_power, root, 3,
                             silex::FactoredRootStrategy::reduced));
    assert(is_power);
    assert(check.pow_si(root, 3));
    assert(check.evaluate(check_value));
    assert(check_value.equal(input_value));

    assert(input.one());
    for (slong i = 0; i < 4; ++i) {
        assert(input.push(base[i], 1));
    }
    assert(input.push(base[4], 3));
    // b0 b1 b2 b3 b4^3 is not a cube (an independent `nfroots` finds no root of
    // x^3 - a; its ideal valuations are not all divisible by 3).  A residue
    // field disproves it, so the answer is a definite "not a cube".
    is_power = true;
    const bool answered = input.is_power_si(
            is_power, root, 3, silex::FactoredRootStrategy::reduced);
    assert(answered && !is_power);

    // Squares: b0 b1 b2 b3 (b0 b1 b2 b3) b4^4 b5^-2.
    assert(input.one());
    for (slong i = 0; i < 4; ++i) {
        assert(input.push(base[i], 1));
    }
    assert(input.push(partial, 1));
    assert(input.push(base[4], 4));
    assert(input.push(base[5], -2));
    assert(input.evaluate(input_value));
    is_power = false;
    assert(input.is_square(is_power, root));
    assert(is_power);
    assert(check.pow_si(root, 2));
    assert(check.evaluate(check_value));
    assert(check_value.equal(input_value));

    assert(input.one());
    for (slong i = 0; i < 5; ++i) {
        assert(input.push(base[i], 1));
    }
    // b0 b1 b2 b3 b4 is not a square (an independent `nfroots` finds no root of
    // x^2 - a); a residue field disproves it.
    is_power = true;
    const bool square_answered = input.is_square(is_power, root);
    assert(square_answered && !is_power);
}

// The reduced root verifies the expanded leftover against a norm computed
// from its factors, prod N(f)^r.  With non-integral factors, factors of
// negative norm, and repeated factors, its answers must equal those of
// Element::is_power on the evaluated value.
void check_reduced_root_factored_norm() noexcept {
    silex::NumberField field = cubic_field();
    silex::Element a(field);
    silex::Element b(field);
    silex::Element c(field);
    silex::Element half(field);
    set_cubic_element(a, -3, 1, 0);  // theta - 3, norm -25
    set_rational(half, 1, 2);
    assert(a.multiply(a, half));     // (theta - 3) / 2
    set_cubic_element(b, 1, 1, 1);
    set_cubic_element(c, 5, -2, 1);

    struct Entry {
        const silex::Element* factor;
        slong exponent;
    };
    struct Case {
        std::vector<Entry> entries;
        slong n;
        bool power;
    };
    const std::vector<Case> cases = {
            {{{&a, 1}, {&b, 3}, {&a, 2}}, 3, true},
            {{{&a, 2}, {&c, 3}, {&a, 3}, {&c, 2}}, 5, true},
            {{{&a, 4}, {&b, 1}, {&c, -3}, {&b, 2}, {&a, -1}}, 3, true},
            {{{&a, 1}, {&c, 2}}, 3, false},
            {{{&a, 3}, {&b, 1}, {&c, 4}}, 5, false},
    };
    for (const Case& test : cases) {
        silex::FactoredElement input(field);
        for (const Entry& entry : test.entries) {
            assert(input.push(*entry.factor, entry.exponent));
        }
        silex::Element value(field);
        assert(input.evaluate(value));

        sflint::Fmpz exponent;
        fmpz_set_si(exponent.raw(), test.n);
        bool plain_power = !test.power;
        silex::Element plain_root(field);
        const bool plain_answered = value.is_power(
                plain_power, plain_root, sflint::FmpzConstRef(exponent));

        bool is_power = !test.power;
        silex::FactoredElement root(field);
        const bool answered = input.is_power_si(
                is_power, root, test.n, silex::FactoredRootStrategy::reduced);
        assert(answered == plain_answered);
        if (!answered) {
            continue;
        }
        assert(is_power == plain_power);
        assert(is_power == test.power);
        if (is_power) {
            silex::FactoredElement check(field);
            silex::Element check_value(field);
            assert(check.pow_si(root, test.n));
            assert(check.evaluate(check_value));
            assert(check_value.equal(value));
        }
    }
}

// In degree 10 and higher the reduced root first screens the factored
// leftover by residues, and evaluates it only when the screen decides
// nothing.  Its answers must equal those of Element::is_power on the
// evaluated value: definite non-powers (also with a factor that is not
// integral at the first usable prime, where the screen gives up), a rational
// cube found after the screen, and a true cube that no step can find, which
// stays unanswered.  Q(zeta_25) has degree 20; 1 + zeta is not a cube
// (`nfroots` finds no root of x^3 - (1 + zeta)).
void check_reduced_root_high_degree_screen() noexcept {
    sflint::FmpqPoly polynomial;
    for (slong i = 0; i <= 20; i += 5) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, 1);
    }
    const silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Element zeta(field);
    silex::Element one_plus_zeta(field);
    silex::Element one_plus_zeta_seventh(field);
    silex::Element two(field);
    silex::Element four(field);
    silex::Element seven(field);
    assert(zeta.gen());
    assert(one_plus_zeta.set_si(1));
    assert(one_plus_zeta.add(one_plus_zeta, zeta));
    assert(one_plus_zeta_seventh.scalar_div_si(one_plus_zeta, 7));
    assert(two.set_si(2));
    assert(four.set_si(4));
    assert(seven.set_si(7));

    struct Entry {
        const silex::Element* factor;
        slong exponent;
    };
    struct Case {
        std::vector<Entry> entries;
        bool answered;
        bool power;
    };
    const std::vector<Case> cases = {
            {{{&one_plus_zeta, 2}, {&zeta, 3}, {&one_plus_zeta, 5}}, true,
             false},
            {{{&two, 3}, {&one_plus_zeta, 1}}, true, false},
            {{{&one_plus_zeta_seventh, 1}, {&seven, 1}}, true, false},
            {{{&two, 1}, {&four, 1}}, true, true},
            {{{&one_plus_zeta, 1}, {&one_plus_zeta, 2}}, false, false},
    };
    for (const Case& test : cases) {
        silex::FactoredElement input(field);
        for (const Entry& entry : test.entries) {
            assert(input.push(*entry.factor, entry.exponent));
        }
        silex::Element value(field);
        assert(input.evaluate(value));

        sflint::Fmpz exponent;
        fmpz_set_si(exponent.raw(), 3);
        bool plain_power = !test.power;
        silex::Element plain_root(field);
        const bool plain_answered = value.is_power(
                plain_power, plain_root, sflint::FmpzConstRef(exponent));

        bool is_power = !test.power;
        silex::FactoredElement root(field);
        const bool answered = input.is_power_si(
                is_power, root, 3, silex::FactoredRootStrategy::reduced);
        assert(answered == test.answered);
        assert(answered == plain_answered);
        if (!answered) {
            continue;
        }
        assert(is_power == plain_power);
        assert(is_power == test.power);
        if (is_power) {
            silex::FactoredElement check(field);
            silex::Element check_value(field);
            assert(check.pow_si(root, 3));
            assert(check.evaluate(check_value));
            assert(check_value.equal(value));
        }
    }
}

// Negative cases over Q, where the answer is definite: the leftover has
// several factors and is not a power.
void check_definite_non_powers() noexcept {
    silex::NumberField field = degree_one_field();
    silex::Element base[3] = {silex::Element(field), silex::Element(field),
                              silex::Element(field)};
    assert(base[0].set_si(2));
    assert(base[1].set_si(3));
    assert(base[2].set_si(5));
    silex::FactoredElement input(field);
    silex::FactoredElement root(field);
    bool is_power = true;
    // 2 * 3 * 5^2 is not a square.
    assert(input.push(base[0], 1));
    assert(input.push(base[1], 1));
    assert(input.push(base[2], 2));
    assert(input.is_square(is_power, root));
    assert(!is_power);
    // 2^3 * 3^3 * 5 is not a cube.
    assert(input.one());
    assert(input.push(base[0], 3));
    assert(input.push(base[1], 3));
    assert(input.push(base[2], 1));
    is_power = true;
    assert(input.is_power_si(is_power, root, 3,
                             silex::FactoredRootStrategy::reduced));
    assert(!is_power);
}

}  // namespace

// The merged factor list of `value`, as the reference's factored-element
// dictionary holds it: equal factors (Element::equal) in order of first
// occurrence, each with the exact fmpz sum of its entries' exponents, zero
// sums dropped.
void merge_factors(std::vector<const silex::Element*>& bases,
                   sflint::FmpzVec& exponents,
                   const silex::FactoredElement& value) noexcept {
    const auto entries = value.factors();
    sflint::FmpzVec sums(static_cast<slong>(entries.size()));
    std::vector<const silex::Element*> distinct;
    for (const auto& entry : entries) {
        std::size_t d = 0;
        while (d < distinct.size() && !distinct[d]->equal(entry.factor)) {
            ++d;
        }
        if (d == distinct.size()) {
            distinct.push_back(&entry.factor);
        }
        ::fmpz_add_si(sums.data() + d, sums.data() + d, entry.exponent);
    }
    bases.clear();
    sflint::FmpzVec kept(static_cast<slong>(distinct.size()));
    for (std::size_t d = 0; d < distinct.size(); ++d) {
        if (::fmpz_is_zero(sums.data() + d) == 0) {
            ::fmpz_set(kept.data() + bases.size(), sums.data() + d);
            bases.push_back(distinct[d]);
        }
    }
    exponents = std::move(kept);
}

// The merged factored log pass: the factor log of every merged base, in
// order, at the bucketed precision of FactoredElement::logarithmic_embedding
// for the working precision `work`, summed with its merged exponent, then
// each coordinate rounded to the shortest form within 2^-precision.  Returns
// false when a partial sum misses the target radius (the caller doubles).
bool merged_factored_log_pass(sflint::ArbVec& out,
                              silex::EmbeddingContext& embeddings,
                              const std::vector<const silex::Element*>& bases,
                              const sflint::FmpzVec& exponents,
                              silex::LogEmbeddingMode mode,
                              slong precision,
                              slong places,
                              slong work) noexcept {
    slong factor_precision = 1;
    while (factor_precision <= work) {
        factor_precision *= 2;
    }

    sflint::ArbVec row(places);
    _arb_vec_zero(out.data(), places);
    for (std::size_t k = 0; k < bases.size(); ++k) {
        if (!silex::logarithmic_embedding(sflint::ArbVecRef(row), embeddings,
                                          *bases[k], mode,
                                          factor_precision)) {
            return false;
        }
        const fmpz* exponent = exponents.data() + k;
        for (slong j = 0; j < places; ++j) {
            if (k == 0) {
                ::arb_mul_fmpz(out.data() + j, row.data() + j, exponent,
                               factor_precision);
            } else {
                ::arb_addmul_fmpz(out.data() + j, row.data() + j, exponent,
                                  factor_precision);
            }
            if (!has_radius_at_most_2exp(out.data() + j, -precision)) {
                return false;
            }
        }
    }
    for (slong j = 0; j < places; ++j) {
        arb_struct* coordinate = out.data() + j;
        if (::arb_rel_accuracy_bits(coordinate) < 0 ||
            ::arb_is_exact(coordinate) != 0) {
            continue;
        }
        slong rounding = ::arb_bits(coordinate) / 2;
        if (rounding < 2) {
            continue;
        }
        sflint::Arb rounded;
        ::arb_set_round(rounded.raw(), coordinate, rounding);
        while (has_radius_at_most_2exp(rounded.raw(), -precision) &&
               rounding > 4) {
            ::arb_set(coordinate, rounded.raw());
            rounding /= 2;
            ::arb_set_round(rounded.raw(), rounded.raw(), rounding);
        }
    }
    return true;
}

// The merged factored log with its own precision doubling: start at the
// target plus guard bits for the largest merged exponent and the number of
// merged bases, and double the working precision until every partial sum
// meets 2^-precision.  `first_pass_only` stops after the first pass (false
// when that pass misses the target).
bool merged_factored_log(sflint::ArbVec& out,
                         silex::EmbeddingContext& embeddings,
                         const silex::FactoredElement& value,
                         silex::LogEmbeddingMode mode,
                         slong precision,
                         slong places,
                         bool first_pass_only = false) noexcept {
    std::vector<const silex::Element*> bases;
    sflint::FmpzVec exponents(0);
    merge_factors(bases, exponents, value);
    if (bases.empty()) {
        _arb_vec_zero(out.data(), places);
        return true;
    }
    slong exponent_bits = 0;
    for (std::size_t k = 0; k < bases.size(); ++k) {
        const slong b = static_cast<slong>(
                ::fmpz_bits(exponents.data() + k));
        exponent_bits = b > exponent_bits ? b : exponent_bits;
    }
    slong count_bits = 0;
    for (std::size_t n = bases.size(); n != 0; n >>= 1) {
        ++count_bits;
    }
    slong work = precision + exponent_bits + count_bits;
    for (int attempt = 0; attempt < 12; ++attempt) {
        if (merged_factored_log_pass(out, embeddings, bases, exponents, mode,
                                     precision, places, work)) {
            return true;
        }
        if (first_pass_only) {
            return false;
        }
        work *= 2;
    }
    return false;
}

// The exact log embedding of `value`, enclosed far more tightly than any
// result checked against it: the log of the evaluated element at high
// precision.
void exact_factored_log(sflint::ArbVec& out,
                        const silex::NumberField& field,
                        const silex::FactoredElement& value,
                        silex::LogEmbeddingMode mode) noexcept {
    silex::Element evaluated(field);
    assert(value.evaluate(evaluated));
    silex::EmbeddingContext embeddings(field);
    assert(silex::logarithmic_embedding(sflint::ArbVecRef(out), embeddings,
                                        evaluated, mode, 1200));
}

// A factored element that repeats equal factors (as separate copies, built
// along different paths, interleaved with other factors) has the log
// embedding of its merged factor list, bit for bit, and that enclosure
// contains the exact log embedding.
void check_factored_log_of(const silex::NumberField& field,
                           const silex::FactoredElement& value,
                           slong places) noexcept {
    for (const silex::LogEmbeddingMode mode :
         {silex::LogEmbeddingMode::plain,
          silex::LogEmbeddingMode::product}) {
        sflint::ArbVec exact(places);
        exact_factored_log(exact, field, value, mode);
        for (const slong precision : {slong(30), slong(64), slong(200)}) {
            silex::EmbeddingContext embeddings(field);
            silex::EmbeddingContext reference_embeddings(field);
            sflint::ArbVec logs(places);
            sflint::ArbVec expected(places);
            assert(value.logarithmic_embedding(sflint::ArbVecRef(logs),
                                               embeddings, mode, precision));
            assert(merged_factored_log(expected, reference_embeddings,
                                       value, mode, precision, places));
            for (slong j = 0; j < places; ++j) {
                assert(::arb_equal(logs.data() + j, expected.data() + j) !=
                       0);
                assert(has_radius_at_most_2exp(logs.data() + j, -precision));
                assert(::arb_contains(logs.data() + j, exact.data() + j) !=
                       0);
            }
        }
    }
}

void check_repeated_factor_logs() noexcept {
    // Degree one, quadratic and cubic fields: the three FLINT element layouts.
    silex::NumberField rational = degree_one_field();
    silex::Element two(rational);
    silex::Element other_two(rational);
    silex::Element three(rational);
    assert(two.set_si(2));
    assert(three.set_si_over_si(3, 7));
    assert(other_two.set_si_over_si(6, 3));
    silex::FactoredElement rational_value(rational);
    assert(rational_value.push(two, 3));
    assert(rational_value.push(three, -2));
    assert(rational_value.push(other_two, -5));
    assert(rational_value.push(three, 1));
    assert(rational_value.push(two, 4));
    check_factored_log_of(rational, rational_value, 1);

    silex::NumberField quadratic = quadratic_field(2);
    silex::Element unit_like(quadratic);
    silex::Element theta(quadratic);
    silex::Element same_unit(quadratic);
    silex::Element shifted(quadratic);
    set_one_plus_theta(unit_like);
    assert(theta.gen());
    assert(same_unit.add_si(theta, 1));
    assert(shifted.add_si(theta, 5));
    silex::FactoredElement quadratic_value(quadratic);
    assert(quadratic_value.push(unit_like, 7));
    assert(quadratic_value.push(shifted, 2));
    assert(quadratic_value.push(same_unit, -3));
    assert(quadratic_value.push(shifted, -1));
    assert(quadratic_value.push(unit_like, 1));
    check_factored_log_of(quadratic, quadratic_value, 2);

    silex::NumberField cubic = cubic_field();
    silex::Element cubic_theta(cubic);
    silex::Element a(cubic);
    silex::Element same_a(cubic);
    silex::Element b(cubic);
    silex::Element c(cubic);
    set_cubic_element(a, 1, 1, 0);
    assert(cubic_theta.gen());
    assert(same_a.add_si(cubic_theta, 1));
    set_cubic_element(b, 3, -1, 1);
    set_cubic_element(c, -5, 2, 7);
    silex::FactoredElement cubic_value(cubic);
    assert(cubic_value.push(a, 3));
    assert(cubic_value.push(b, -2));
    assert(cubic_value.push(same_a, 5));
    assert(cubic_value.push(b, 1));
    assert(cubic_value.push(a, -4));
    assert(cubic_value.push(c, 2));
    check_factored_log_of(cubic, cubic_value, 2);

    // Values whose factors are copies made by the factored operations (in
    // place and not, from a shared list), mixed with separately pushed equal
    // factors, merge the same way.
    silex::FactoredElement cubic_power(cubic);
    assert(cubic_power.pow_si(cubic_value, 2));
    silex::FactoredElement cubic_mixed(cubic);
    assert(cubic_mixed.set(cubic_power));
    assert(cubic_mixed.multiply(cubic_mixed, cubic_value));
    assert(cubic_mixed.push(same_a, -7));
    assert(cubic_mixed.multiply(cubic_mixed, cubic_mixed));
    check_factored_log_of(cubic, cubic_mixed, 2);
    silex::FactoredElement cubic_inverse(cubic);
    assert(cubic_inverse.invert(cubic_value));
    silex::FactoredElement cubic_ratio(cubic);
    assert(cubic_ratio.divide(cubic_power, cubic_inverse));
    assert(cubic_ratio.push(a, 2));
    check_factored_log_of(cubic, cubic_ratio, 2);

    // Hundreds of entries over three distinct factors, as unit candidates
    // built from relation products have.
    silex::FactoredElement long_value(cubic);
    const silex::Element* const bases[] = {&a, &b, &same_a, &c};
    for (slong i = 0; i < 300; ++i) {
        const slong exponent = (i % 7) - 3;
        if (exponent != 0) {
            assert(long_value.push(*bases[(i * 5) % 4], exponent));
        }
    }
    check_factored_log_of(cubic, long_value, 2);
}

// Equal bases are merged before the log sum, with exponents summed exactly.
// Each case below has a result that only the merged sum gives: an entry-by-
// entry sum of irrational logs is never an exact zero, and it rounds each
// cancelled entry into the radius.
void check_merged_equal_bases() noexcept {
    silex::NumberField quadratic = quadratic_field(2);
    silex::Element unit_like(quadratic);
    silex::Element theta(quadratic);
    silex::Element same_unit(quadratic);
    silex::Element shifted(quadratic);
    set_one_plus_theta(unit_like);
    assert(theta.gen());
    assert(same_unit.add_si(theta, 1));
    assert(shifted.add_si(theta, 5));

    // Exponents that cancel completely: the value is 1 and its log is 0
    // exactly.
    silex::FactoredElement cancelled(quadratic);
    assert(cancelled.push(unit_like, 7));
    assert(cancelled.push(shifted, 3));
    assert(cancelled.push(same_unit, -4));
    assert(cancelled.push(shifted, -3));
    assert(cancelled.push(unit_like, -3));
    // Exponents that cancel on one base and add on the other: the merged
    // list is shifted^2 alone, so the log is that of a one-entry product.
    silex::FactoredElement partial(quadratic);
    assert(partial.push(unit_like, 9));
    assert(partial.push(shifted, 3));
    assert(partial.push(same_unit, -9));
    assert(partial.push(shifted, -1));
    silex::FactoredElement single(quadratic);
    assert(single.push(shifted, 2));
    for (const silex::LogEmbeddingMode mode :
         {silex::LogEmbeddingMode::plain,
          silex::LogEmbeddingMode::product}) {
        for (const slong precision : {slong(30), slong(64), slong(200)}) {
            silex::EmbeddingContext embeddings(quadratic);
            sflint::ArbVec logs(2);
            assert(cancelled.logarithmic_embedding(sflint::ArbVecRef(logs),
                                                   embeddings, mode,
                                                   precision));
            for (slong j = 0; j < 2; ++j) {
                assert(::arb_is_zero(logs.data() + j) != 0);
                assert(::arb_is_exact(logs.data() + j) != 0);
            }

            sflint::ArbVec merged(2);
            sflint::ArbVec alone(2);
            silex::EmbeddingContext other_embeddings(quadratic);
            assert(partial.logarithmic_embedding(sflint::ArbVecRef(merged),
                                                 embeddings, mode,
                                                 precision));
            assert(single.logarithmic_embedding(sflint::ArbVecRef(alone),
                                                other_embeddings, mode,
                                                precision));
            for (slong j = 0; j < 2; ++j) {
                assert(::arb_equal(merged.data() + j, alone.data() + j) !=
                       0);
            }
        }
    }
    check_factored_log_of(quadratic, cancelled, 2);
    check_factored_log_of(quadratic, partial, 2);

    // Exponents whose sum leaves slong: a^max * a^max * b has the merged
    // exponent 2 * (2^63 - 1) on a.  The result must contain the exact log,
    // (2^64 - 2) log|a| + log|b| at each place (a cubic field, with a real
    // and a complex place).
    silex::NumberField cubic = cubic_field();
    silex::Element a(cubic);
    silex::Element same_a(cubic);
    silex::Element b(cubic);
    silex::Element cubic_theta(cubic);
    set_cubic_element(a, 1, 1, 0);
    assert(cubic_theta.gen());
    assert(same_a.add_si(cubic_theta, 1));
    set_cubic_element(b, 3, -1, 1);
    const slong max = std::numeric_limits<slong>::max();
    silex::FactoredElement wide(cubic);
    assert(wide.push(a, max));
    assert(wide.push(b, 1));
    assert(wide.push(same_a, max));
    sflint::Fmpz wide_exponent;
    sflint::fmpz_set_si(wide_exponent, max);
    ::fmpz_mul_ui(wide_exponent.raw(), wide_exponent.raw(), 2);
    for (const silex::LogEmbeddingMode mode :
         {silex::LogEmbeddingMode::plain,
          silex::LogEmbeddingMode::product}) {
        silex::EmbeddingContext exact_embeddings(cubic);
        sflint::ArbVec log_a(2);
        sflint::ArbVec log_b(2);
        assert(silex::logarithmic_embedding(sflint::ArbVecRef(log_a),
                                            exact_embeddings, a, mode, 1200));
        assert(silex::logarithmic_embedding(sflint::ArbVecRef(log_b),
                                            exact_embeddings, b, mode, 1200));
        sflint::ArbVec exact(2);
        for (slong j = 0; j < 2; ++j) {
            ::arb_mul_fmpz(exact.data() + j, log_a.data() + j,
                           wide_exponent.raw(), 1200);
            ::arb_add(exact.data() + j, exact.data() + j, log_b.data() + j,
                      1200);
        }
        for (const slong precision : {slong(30), slong(64), slong(200)}) {
            silex::EmbeddingContext embeddings(cubic);
            silex::EmbeddingContext reference_embeddings(cubic);
            sflint::ArbVec logs(2);
            sflint::ArbVec expected(2);
            assert(wide.logarithmic_embedding(sflint::ArbVecRef(logs),
                                              embeddings, mode, precision));
            assert(merged_factored_log(expected, reference_embeddings, wide,
                                       mode, precision, 2));
            for (slong j = 0; j < 2; ++j) {
                assert(::arb_equal(logs.data() + j, expected.data() + j) !=
                       0);
                assert(has_radius_at_most_2exp(logs.data() + j, -precision));
                assert(::arb_contains(logs.data() + j, exact.data() + j) !=
                       0);
            }
        }
    }
}

// A context whose cache already holds the factor logs at a higher precision
// (stored by an earlier call), and a value whose first pass misses the target
// so a second precision pass runs: both must equal the merged reference,
// which doubles on its own.
void check_warm_cache_and_second_pass() noexcept {
    silex::NumberField cubic = cubic_field();
    silex::Element a(cubic);
    silex::Element b(cubic);
    set_cubic_element(a, 1, 1, 0);
    set_cubic_element(b, 3, -1, 1);
    silex::FactoredElement warm_value(cubic);
    assert(warm_value.push(a, 3));
    assert(warm_value.push(b, -2));
    assert(warm_value.push(a, 5));
    assert(warm_value.push(b, 1));
    for (const silex::LogEmbeddingMode mode :
         {silex::LogEmbeddingMode::plain,
          silex::LogEmbeddingMode::product}) {
        silex::EmbeddingContext embeddings(cubic);
        silex::EmbeddingContext reference_embeddings(cubic);
        sflint::ArbVec warmup(2);
        assert(warm_value.logarithmic_embedding(sflint::ArbVecRef(warmup),
                                                embeddings, mode, 400));
        // The reference context is warmed the same way.
        sflint::ArbVec reference_warmup(2);
        assert(merged_factored_log(reference_warmup, reference_embeddings,
                                      warm_value, mode, 400, 2));
        std::uint64_t version = 0;
        assert(silex::detail::embedding_log_cache_version(version,
                                                          embeddings));
        sflint::ArbVec logs(2);
        sflint::ArbVec expected(2);
        assert(warm_value.logarithmic_embedding(sflint::ArbVecRef(logs),
                                                embeddings, mode, 30));
        assert(merged_factored_log(expected, reference_embeddings,
                                      warm_value, mode, 30, 2));
        for (slong j = 0; j < 2; ++j) {
            assert(::arb_equal(logs.data() + j, expected.data() + j) != 0);
            assert(has_radius_at_most_2exp(logs.data() + j, -30));
        }
        // The same call again changes nothing in the cache.
        std::uint64_t after = 0;
        assert(silex::detail::embedding_log_cache_version(after,
                                                          embeddings));
        sflint::ArbVec again(2);
        assert(warm_value.logarithmic_embedding(sflint::ArbVecRef(again),
                                                embeddings, mode, 30));
        std::uint64_t repeated = 0;
        assert(silex::detail::embedding_log_cache_version(repeated,
                                                          embeddings));
        assert(repeated == after);
        for (slong j = 0; j < 2; ++j) {
            assert(::arb_equal(again.data() + j, logs.data() + j) != 0);
        }
    }

    // (1 + theta)^n has a tiny conjugate, so the factor log at the first
    // bucket is too wide and the pass must repeat at doubled precision.  The
    // easy factor comes first so its row is ready before the pass fails.
    silex::NumberField quadratic = quadratic_field(2);
    silex::Element base(quadratic);
    silex::Element power(quadratic);
    set_one_plus_theta(base);
    sflint::Fmpz n;
    sflint::fmpz_set_si(n, 300);
    assert(power.pow_fmpz(base, sflint::FmpzConstRef(n)));
    silex::Element other(quadratic);
    assert(other.gen());
    silex::FactoredElement hard(quadratic);
    assert(hard.push(other, 2));
    assert(hard.push(power, 1));
    assert(hard.push(power, -1));
    assert(hard.push(power, 2));
    for (const silex::LogEmbeddingMode mode :
         {silex::LogEmbeddingMode::plain,
          silex::LogEmbeddingMode::product}) {
        silex::EmbeddingContext embeddings(quadratic);
        silex::EmbeddingContext reference_embeddings(quadratic);
        sflint::ArbVec logs(2);
        sflint::ArbVec expected(2);
        sflint::ArbVec first_pass(2);
        assert(!merged_factored_log(first_pass, reference_embeddings, hard,
                                       mode, 30, 2, true));
        assert(hard.logarithmic_embedding(sflint::ArbVecRef(logs), embeddings,
                                          mode, 30));
        assert(merged_factored_log(expected, reference_embeddings, hard,
                                      mode, 30, 2));
        for (slong j = 0; j < 2; ++j) {
            assert(::arb_equal(logs.data() + j, expected.data() + j) != 0);
            assert(has_radius_at_most_2exp(logs.data() + j, -30));
        }
    }
}

// The factor list of `value` is exactly `expected`, in order: equal factors
// and equal exponents, nothing merged or reordered.
bool has_factor_list(const silex::FactoredElement& value,
                     const std::vector<const silex::Element*>& factors,
                     const std::vector<slong>& exponents) noexcept {
    const auto entries = value.factors();
    if (entries.size() != factors.size() ||
        factors.size() != exponents.size()) {
        return false;
    }
    for (std::size_t i = 0; i < entries.size(); ++i) {
        if (!entries[i].factor.equal(*factors[i]) ||
            entries[i].exponent != exponents[i]) {
            return false;
        }
    }
    return true;
}

// A copy made with set() may share factor storage with its source.  Every
// mutation of one value, including in-place forms whose output aliases an
// input, must leave the other value and the factors it exposes unchanged.
void check_copies_are_independent() noexcept {
    silex::NumberField field = quadratic_field(-5);
    silex::Element a(field);
    silex::Element b(field);
    silex::Element c(field);
    set_one_plus_theta(a);
    assert(b.set_si(3));
    assert(c.set_si(-7));

    silex::FactoredElement source(field);
    assert(source.push(a, 2));
    assert(source.push(b, -1));
    assert(source.push(a, 5));
    const std::vector<const silex::Element*> source_factors{&a, &b, &a};
    const std::vector<slong> source_exponents{2, -1, 5};

    silex::FactoredElement extra(field);
    assert(extra.push(c, 3));

    // Mutating the copy leaves the source unchanged.
    {
        silex::FactoredElement copy(field);
        assert(copy.set(source));
        assert(has_factor_list(copy, source_factors, source_exponents));
        const silex::Element* source_first = source.factor(0);
        assert(copy.push(c, 4));
        assert(copy.multiply(copy, extra));
        assert(copy.divide(copy, extra));
        assert(copy.pow_si(copy, 2));
        assert(has_factor_list(source, source_factors, source_exponents));
        assert(source.factor(0) == source_first);
        assert(has_factor_list(copy, {&a, &b, &a, &c, &c, &c},
                               {4, -2, 10, 8, 6, -6}));
    }

    // Mutating the source leaves the copy unchanged, also after the source
    // is cleared.
    {
        silex::FactoredElement working(field);
        silex::FactoredElement copy(field);
        assert(working.set(source));
        assert(copy.set(working));
        const silex::Element* copy_first = copy.factor(0);
        assert(working.multiply(working, working));
        assert(has_factor_list(working, {&a, &b, &a, &a, &b, &a},
                               {2, -1, 5, 2, -1, 5}));
        assert(has_factor_list(copy, source_factors, source_exponents));
        assert(copy.factor(0) == copy_first);
        assert(working.push(*working.factor(1), 7));
        assert(working.length() == 7);
        working.clear();
        assert(has_factor_list(copy, source_factors, source_exponents));
        assert(copy.factor(0) == copy_first);
    }

    // A push of one of the value's own factors copies it before the list
    // changes, whether or not the storage is shared.
    {
        silex::FactoredElement shared(field);
        silex::FactoredElement copy(field);
        assert(shared.set(source));
        assert(copy.set(shared));
        assert(shared.push(*shared.factor(2), -3));
        assert(has_factor_list(shared, {&a, &b, &a, &a}, {2, -1, 5, -3}));
        assert(has_factor_list(copy, source_factors, source_exponents));
        assert(copy.push(*copy.factor(1), 1));
        assert(has_factor_list(copy, {&a, &b, &a, &b}, {2, -1, 5, 1}));
        assert(has_factor_list(source, source_factors, source_exponents));
    }

    // Products with an empty operand and set_element do not touch the other
    // operand.
    {
        silex::FactoredElement empty(field);
        silex::FactoredElement product(field);
        assert(empty.one());
        assert(product.multiply(empty, source));
        assert(has_factor_list(product, source_factors, source_exponents));
        assert(product.multiply(product, empty));
        assert(product.push(c, 1));
        assert(product.set_element(b));
        assert(has_factor_list(product, {&b}, {1}));
        assert(has_factor_list(source, source_factors, source_exponents));
        assert(empty.length() == 0);
    }

    // A failed in-place product leaves the value unchanged.
    {
        silex::FactoredElement copy(field);
        silex::FactoredElement huge(field);
        assert(copy.set(source));
        assert(huge.push(c, std::numeric_limits<slong>::min()));
        assert(!copy.divide(copy, huge));
        assert(has_factor_list(copy, source_factors, source_exponents));
        assert(has_factor_list(source, source_factors, source_exponents));
    }

    // The same failure on an unshared list (the in-place path).
    {
        silex::FactoredElement own(field);
        silex::FactoredElement huge(field);
        assert(own.push(a, 2));
        assert(huge.push(c, std::numeric_limits<slong>::min()));
        assert(!own.divide(own, huge));
        assert(has_factor_list(own, {&a}, {2}));
    }
}

int main() {
    check_copies_are_independent();
    check_many_factor_evaluation();
    check_repeated_factor_logs();
    check_merged_equal_bases();
    check_warm_cache_and_second_pass();
    check_multi_factor_leftover_roots();
    check_definite_non_powers();
    check_reduced_root_factored_norm();
    check_reduced_root_high_degree_screen();

    silex::NumberField field = degree_one_field();
    silex::NumberField other_field = degree_one_field();

    silex::Element two(field);
    silex::Element three(field);
    silex::Element four(field);
    silex::Element eight(field);
    silex::Element zero(field);
    silex::Element other_parent(other_field);
    assert(two.set_si(2));
    assert(three.set_si(3));
    assert(four.set_si(4));
    assert(eight.set_si(8));
    assert(zero.zero());
    assert(other_parent.set_si(2));

    silex::FactoredElement product(field);
    assert(product.is_defined());
    assert(product.parent() != nullptr &&
           product.parent()->has_same_data(field));
    assert(product.length() == 0);

    slong exponent = 3;
    assert(product.push(two, exponent));
    exponent = -1;
    assert(product.push(three, exponent));
    assert(product.length() == 2);

    silex::Element evaluated(field);
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    silex::FactoredElement recursive_eval(field);
    assert(recursive_eval.push(two, 13));
    assert(recursive_eval.push(three, -11));
    assert(recursive_eval.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8192, 177147));

    silex::FactoredElement single(field);
    exponent = 2;
    assert(single.push(two, exponent));

    silex::FactoredElement combined(field);
    assert(combined.multiply(product, single));
    assert(combined.length() == 3);
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 32, 3));

    assert(combined.divide(combined, single));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    assert(combined.invert(product));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 8));

    assert(combined.multiply(combined, product));
    assert(combined.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    exponent = -2;
    assert(combined.pow_si(product, exponent));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 9, 64));

    exponent = 0;
    assert(combined.pow_si(product, exponent));
    assert(combined.length() == 0);
    assert(combined.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    silex::FactoredElement copy(field);
    assert(copy.set(product));
    assert(copy.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    silex::FactoredElement moved(std::move(copy));
    assert(moved.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    assert(!product.push(zero, exponent));
    assert(!product.push(other_parent, exponent));

    silex::FactoredElement other_product(other_field);
    assert(other_product.set_element(other_parent));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.multiply(product, other_product));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.divide(product, other_product));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.invert(other_product));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));
    assert(!product.pow_si(other_product, exponent));
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 3));

    silex::FactoredElement root_input(field);
    assert(root_input.push(two, 6));
    assert(root_input.push(three, -4));
    assert(root_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 81));

    silex::FactoredElement root(field);
    assert(root.root_si(root_input, 2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 9));
    assert(combined.pow_si(root, 2));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 81));

    assert(root.root_si(root_input, 1));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 81));

    assert(root.root_si(root, 2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8, 9));

    silex::FactoredElement empty_root(field);
    assert(empty_root.root_si(empty_root, 7));
    assert(empty_root.length() == 0);
    assert(empty_root.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    assert(root.set_element(three));
    assert(!root.root_si(root_input, 5));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!root.root_si(root_input, 0));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!root.root_si(root_input, -2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!root.root_si(other_product, 2));
    assert(root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));

    bool is_square = false;
    silex::FactoredElement square_input(field);
    silex::FactoredElement square_root(field);
    silex::FactoredElement square_power(field);
    assert(square_input.push(two, 4));
    assert(square_input.push(three, -2));
    assert(square_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 16, 9));
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_power.pow_si(square_root, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 16, 9));

    assert(square_input.is_square(is_square, square_input));
    assert(is_square);
    assert(square_power.pow_si(square_input, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 16, 9));

    assert(square_input.one());
    assert(square_input.push(four, 1));
    assert(square_input.push(three, 2));
    assert(square_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 36, 1));
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_power.pow_si(square_root, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 36, 1));

    assert(square_input.one());
    assert(square_input.push(four, -1));
    assert(square_input.push(three, 2));
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 2));
    assert(square_power.pow_si(square_root, 2));
    assert(square_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 9, 4));

    assert(square_input.one());
    assert(square_input.push(two, 1));
    assert(square_root.set_element(four));
    is_square = true;
    assert(square_input.is_square(is_square, square_root));
    assert(!is_square);
    assert(square_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    assert(square_input.one());
    assert(square_input.is_square(is_square, square_root));
    assert(is_square);
    assert(square_root.length() == 0);
    assert(!square_input.is_square(is_square, other_product));

    bool is_power = false;
    silex::FactoredElement power_input(field);
    silex::FactoredElement power_root(field);
    silex::FactoredElement power_check(field);
    assert(power_input.push(two, 6));
    assert(power_input.push(three, -3));
    assert(power_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));
    assert(power_input.is_power_si(is_power, power_root, 3,
                                   silex::FactoredRootStrategy::reduced));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));

    assert(power_input.is_power_si(is_power, power_input, 3));
    assert(is_power);
    assert(power_check.pow_si(power_input, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 27));

    assert(power_input.one());
    assert(power_input.push(four, 1));
    assert(power_input.is_power_si(is_power, power_root, 2));
    assert(is_power);
    assert(power_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 2, 1));

    assert(power_input.one());
    assert(power_input.push(two, 3));
    assert(power_input.push(three, 6));
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));

    assert(power_input.one());
    assert(power_input.push(eight, -1));
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(is_power);
    assert(power_check.pow_si(power_root, 3));
    assert(power_check.evaluate(evaluated));
    assert(element_is_rational(evaluated, 1, 8));

    assert(power_input.one());
    assert(power_input.push(two, 1));
    assert(power_root.set_element(four));
    is_power = true;
    assert(power_input.is_power_si(is_power, power_root, 3));
    assert(!is_power);
    assert(power_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    assert(power_input.one());
    assert(power_input.is_power_si(is_power, power_input, 1));
    assert(is_power);
    assert(power_input.length() == 0);

    assert(power_root.set_element(four));
    is_power = true;
    assert(!power_input.is_power_si(is_power, power_root, 0));
    assert(is_power);
    assert(power_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));
    assert(!power_input.is_power_si(is_power, other_product, 3));

    silex::CompactElement compact(field);
    assert(compact.is_defined());
    assert(compact.parent() != nullptr &&
           compact.parent()->has_same_data(field));
    assert(compact.one(3));
    assert(compact.base() == 3);
    assert(compact.length() == 0);
    assert(compact.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    silex::FactoredElement compact_input(field);
    silex::FactoredElement compact_root(field);
    silex::FactoredElement compact_power(field);
    assert(compact_input.push(two, 3));
    assert(compact_input.push(three, 6));
    assert(compact_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));
    assert(compact_input.is_power_si(
            is_power, compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(compact_root.length() == 2);
    slong compact_exp = 0;
    assert(compact_root.exponent(compact_exp, 0));
    assert(compact_exp == 1);
    assert(compact_root.exponent(compact_exp, 1));
    assert(compact_exp == 2);
    assert(compact.set_factored_element(compact_input, 3));
    assert(compact.length() == 2);
    assert(compact.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));
    assert(compact.root_base(compact_root));
    assert(compact_power.pow_si(compact_root, 3));
    assert(compact_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 8 * 729, 1));

    assert(compact_input.one());
    assert(compact_input.push(two, 4));
    assert(compact_input.push(four, 1));
    assert(!compact_root.root_si(compact_input, 3));
    assert(compact_input.is_power_si(
            is_power, compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));
    assert(compact_power.pow_si(compact_root, 3));
    assert(compact_power.evaluate(evaluated));
    assert(element_is_rational(evaluated, 64, 1));

    assert(compact_input.one());
    assert(compact_input.push(two, 5));
    assert(compact_input.push(three, -4));
    assert(compact_input.evaluate(evaluated));
    assert(element_is_rational(evaluated, 32, 81));
    assert(compact.set_factored_element(compact_input, 3));
    assert(compact.evaluate(evaluated));
    assert(element_is_rational(evaluated, 32, 81));

    assert(compact_input.one());
    assert(compact_input.push(two, 1));
    assert(compact.set_factored_element(compact_input, 3));
    assert(compact_root.set_element(four));
    assert(!compact.root_base(compact_root));
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));
    is_power = true;
    assert(compact_input.is_power_si(
            is_power, compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(!is_power);
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    assert(compact_input.one());
    assert(!compact.set_factored_element(compact_input, 1));
    assert(!compact_input.is_power_si(is_power, compact_root, -3));
    assert(!compact.root_base(other_product));
    is_power = true;
    assert(!compact_input.is_power_si(
            is_power, compact_root, 3,
            static_cast<silex::FactoredRootStrategy>(99)));
    assert(is_power);
    assert(compact_root.evaluate(evaluated));
    assert(element_is_rational(evaluated, 4, 1));

    silex::FactoredElement overflow_input(field);
    assert(overflow_input.push(two, std::numeric_limits<slong>::min()));
    assert(combined.set_element(three));
    assert(!combined.invert(overflow_input));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!combined.pow_si(overflow_input, -1));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));
    assert(!combined.divide(product, overflow_input));
    assert(combined.evaluate(evaluated));
    assert(element_is_rational(evaluated, 3, 1));

    assert(product.one());
    assert(product.length() == 0);
    assert(product.evaluate(evaluated));
    assert(evaluated.equal_si(1));

    set_rational(two, 5, 2);
    assert(product.set_element(two));
    assert(product.length() == 1);
    assert(product.evaluate(evaluated));
    assert(element_is_rational(evaluated, 5, 2));
    assert(!product.set_element(zero));

    silex::EmbeddingContext rational_embeddings(field);
    silex::Element log_two(field);
    silex::Element log_four(field);
    silex::FactoredElement guarded_logs(field);
    sflint::ArbVec guarded_log_row(1);
    const slong large_even_exponent =
            (std::numeric_limits<slong>::max() / 2048) * 2;
    assert(log_two.set_si(2));
    assert(log_four.set_si(4));
    assert(guarded_logs.push(log_two, large_even_exponent));
    assert(guarded_logs.push(log_four, -large_even_exponent / 2));
    assert(guarded_logs.logarithmic_embedding(
            sflint::ArbVecRef(guarded_log_row), rational_embeddings,
            silex::LogEmbeddingMode::product, 32));
    assert(sflint::arb_contains_si(guarded_log_row.data() + 0, 0));
    assert(has_radius_at_most_2exp(guarded_log_row.data() + 0, -32));

    silex::NumberField quadratic = quadratic_field(2);
    silex::NumberField other_quadratic = quadratic_field(3);
    silex::EmbeddingContext embeddings(quadratic);
    silex::EmbeddingContext other_embeddings(other_quadratic);
    silex::Element unit_like(quadratic);
    silex::Element expanded(quadratic);
    sflint::ArbVec elt_logs(2);
    set_one_plus_theta(unit_like);

    silex::FactoredElement quadratic_square(quadratic);
    silex::FactoredElement quadratic_root(quadratic);
    silex::FactoredElement quadratic_power(quadratic);
    silex::Element quadratic_evaluated(quadratic);
    assert(expanded.multiply(unit_like, unit_like));
    assert(quadratic_square.push(expanded, 1));
    assert(quadratic_square.is_square(is_square, quadratic_root));
    assert(is_square);
    assert(quadratic_power.pow_si(quadratic_root, 2));
    assert(quadratic_power.evaluate(quadratic_evaluated));
    assert(quadratic_evaluated.equal(expanded));

    silex::NumberField cubic = cubic_field();
    silex::Element cubic_eight(cubic);
    silex::Element cubic_evaluated(cubic);
    silex::FactoredElement cubic_power(cubic);
    silex::FactoredElement cubic_root(cubic);
    silex::FactoredElement cubic_check(cubic);
    silex::FactoredElement cubic_compact_root(cubic);
    silex::Element cubic_theta(cubic);
    silex::Element cubic_shift(cubic);
    silex::Element cubic_cube(cubic);
    assert(cubic_eight.set_si(8));
    assert(cubic_power.push(cubic_eight, 1));
    assert(cubic_power.is_power_si(is_power, cubic_root, 3));
    assert(is_power);
    assert(cubic_root.evaluate(cubic_evaluated));
    assert(cubic_evaluated.equal_si(2));
    assert(cubic_check.pow_si(cubic_root, 3));
    assert(cubic_check.evaluate(cubic_evaluated));
    assert(cubic_evaluated.equal(cubic_eight));
    assert(cubic_theta.gen());
    assert(cubic_shift.add_si(cubic_theta, 1));
    sflint::Fmpz cubic_exponent;
    sflint::fmpz_set_ui(cubic_exponent, 3);
    assert(cubic_cube.pow_fmpz(cubic_shift,
                               sflint::FmpzConstRef(cubic_exponent)));
    assert(cubic_power.one());
    assert(cubic_power.push(cubic_cube, 1));
    is_power = false;
    assert(cubic_power.is_power_si(
            is_power, cubic_compact_root, 3,
            silex::FactoredRootStrategy::compact));
    assert(is_power);
    assert(cubic_check.pow_si(cubic_compact_root, 3));
    assert(cubic_check.evaluate(cubic_evaluated));
    assert(cubic_evaluated.equal(cubic_cube));

    silex::FactoredElement compact_log_source(quadratic);
    silex::CompactElement compact_logs(quadratic);
    silex::Element compact_log_value(quadratic);
    sflint::ArbVec compact_log_row(2);
    assert(compact_log_source.push(unit_like, 5));
    assert(compact_logs.set_factored_element(compact_log_source, 3));
    assert(compact_logs.evaluate(compact_log_value));
    assert(compact_logs.logarithmic_embedding(
            sflint::ArbVecRef(compact_log_row), embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(silex::logarithmic_embedding(
            sflint::ArbVecRef(elt_logs), embeddings, compact_log_value,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_overlaps(
            sflint::ArbConstRef(compact_log_row.data() + 0),
            sflint::ArbConstRef(elt_logs.data() + 0)));
    assert(sflint::arb_overlaps(
            sflint::ArbConstRef(compact_log_row.data() + 1),
            sflint::ArbConstRef(elt_logs.data() + 1)));

    silex::FactoredElement empty_logs(quadratic);
    sflint::ArbVec fac_logs(2);
    assert(empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 0));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 0));

    exponent = 2;
    assert(empty_logs.push(unit_like, exponent));
    assert(empty_logs.evaluate(expanded));
    assert(empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(silex::logarithmic_embedding(
            sflint::ArbVecRef(elt_logs), embeddings, expanded,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_overlaps(sflint::ArbConstRef(fac_logs.data() + 0),
                                sflint::ArbConstRef(elt_logs.data() + 0)));
    assert(sflint::arb_overlaps(sflint::ArbConstRef(fac_logs.data() + 1),
                                sflint::ArbConstRef(elt_logs.data() + 1)));

    sflint::arb_set_si(fac_logs.data() + 0, 123);
    sflint::arb_set_si(fac_logs.data() + 1, 456);
    assert(!empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), other_embeddings,
            silex::LogEmbeddingMode::product, 160));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 123));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 456));
    assert(!empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            static_cast<silex::LogEmbeddingMode>(99), 160));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 123));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 456));
    assert(!empty_logs.logarithmic_embedding(
            sflint::ArbVecRef(fac_logs), embeddings,
            silex::LogEmbeddingMode::product, 0));
    assert(sflint::arb_contains_si(fac_logs.data() + 0, 123));
    assert(sflint::arb_contains_si(fac_logs.data() + 1, 456));

    return 0;
}
