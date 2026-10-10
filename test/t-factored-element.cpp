#include <silex/archimedean.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/number_field.hpp>

#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>

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
    // In this cubic field the exact cube test reports an unsupported
    // certificate (false) for this non-cube, so this case only checks that the
    // call does not crash or claim a power wrongly; the definite negative
    // cases are in check_definite_non_powers.
    is_power = true;
    const bool answered = input.is_power_si(
            is_power, root, 3, silex::FactoredRootStrategy::reduced);
    assert(!answered || !is_power);

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
    is_power = true;
    const bool square_answered = input.is_square(is_power, root);
    assert(!square_answered || !is_power);
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

// The per-entry factored log pass: the factor log of every entry, in order,
// at the bucketed precision of FactoredElement::logarithmic_embedding for the
// working precision `work`, summed entry by entry with its own exponent, then
// each coordinate rounded to the shortest form within 2^-precision.  Returns
// false when a partial sum misses the target radius (the caller doubles).
bool per_entry_factored_log_pass(sflint::ArbVec& out,
                                 silex::EmbeddingContext& embeddings,
                                 const silex::FactoredElement& value,
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
    bool first = true;
    for (const auto& entry : value.factors()) {
        if (!silex::logarithmic_embedding(sflint::ArbVecRef(row), embeddings,
                                          entry.factor, mode,
                                          factor_precision)) {
            return false;
        }
        for (slong j = 0; j < places; ++j) {
            if (first) {
                ::arb_mul_si(out.data() + j, row.data() + j, entry.exponent,
                             factor_precision);
            } else {
                ::arb_addmul_si(out.data() + j, row.data() + j,
                                entry.exponent, factor_precision);
            }
            if (!has_radius_at_most_2exp(out.data() + j, -precision)) {
                return false;
            }
        }
        first = false;
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

// The per-entry factored log with its own precision doubling: start at the
// target plus exponent and factor-count guard bits and double the working
// precision until every partial sum meets 2^-precision.  `first_pass_only`
// stops after the first pass (false when that pass misses the target).
bool per_entry_factored_log(sflint::ArbVec& out,
                            silex::EmbeddingContext& embeddings,
                            const silex::FactoredElement& value,
                            silex::LogEmbeddingMode mode,
                            slong precision,
                            slong places,
                            bool first_pass_only = false) noexcept {
    const auto bits = [](ulong x) noexcept {
        slong count = 0;
        while (x != 0) {
            ++count;
            x >>= 1;
        }
        return count;
    };
    slong exponent_bits = 0;
    for (const auto& entry : value.factors()) {
        const slong b = bits(entry.exponent < 0
                                     ? -static_cast<ulong>(entry.exponent)
                                     : static_cast<ulong>(entry.exponent));
        exponent_bits = b > exponent_bits ? b : exponent_bits;
    }
    slong work = precision + exponent_bits +
                 bits(static_cast<ulong>(value.factors().size()));
    for (int attempt = 0; attempt < 12; ++attempt) {
        if (per_entry_factored_log_pass(out, embeddings, value, mode,
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


// A factored element that repeats equal factors (as separate copies, built
// along different paths, interleaved with other factors) has the log
// embedding of the per-entry sum, bit for bit: looking each distinct factor
// up once per pass changes neither the terms nor their order.
void check_factored_log_of(const silex::NumberField& field,
                           const silex::FactoredElement& value,
                           slong places) noexcept {
    for (const silex::LogEmbeddingMode mode :
         {silex::LogEmbeddingMode::plain,
          silex::LogEmbeddingMode::product}) {
        for (const slong precision : {slong(30), slong(64), slong(200)}) {
            silex::EmbeddingContext embeddings(field);
            silex::EmbeddingContext reference_embeddings(field);
            sflint::ArbVec logs(places);
            sflint::ArbVec expected(places);
            assert(value.logarithmic_embedding(sflint::ArbVecRef(logs),
                                               embeddings, mode, precision));
            assert(per_entry_factored_log(expected, reference_embeddings,
                                          value, mode, precision, places));
            for (slong j = 0; j < places; ++j) {
                assert(::arb_equal(logs.data() + j, expected.data() + j) !=
                       0);
                assert(has_radius_at_most_2exp(logs.data() + j, -precision));
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

// A context whose cache already holds the factor logs at a higher precision
// (stored by an earlier call), and a value whose first pass misses the target
// so a second precision pass runs: both must equal the per-entry reference,
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
        // The reference context is warmed the same way, entry by entry.
        sflint::ArbVec reference_warmup(2);
        assert(per_entry_factored_log(reference_warmup, reference_embeddings,
                                      warm_value, mode, 400, 2));
        std::uint64_t version = 0;
        assert(silex::detail::embedding_log_cache_version(version,
                                                          embeddings));
        sflint::ArbVec logs(2);
        sflint::ArbVec expected(2);
        assert(warm_value.logarithmic_embedding(sflint::ArbVecRef(logs),
                                                embeddings, mode, 30));
        assert(per_entry_factored_log(expected, reference_embeddings,
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
        assert(!per_entry_factored_log(first_pass, reference_embeddings, hard,
                                       mode, 30, 2, true));
        assert(hard.logarithmic_embedding(sflint::ArbVecRef(logs), embeddings,
                                          mode, 30));
        assert(per_entry_factored_log(expected, reference_embeddings, hard,
                                      mode, 30, 2));
        for (slong j = 0; j < 2; ++j) {
            assert(::arb_equal(logs.data() + j, expected.data() + j) != 0);
            assert(has_radius_at_most_2exp(logs.data() + j, -30));
        }
    }
}

int main() {
    check_many_factor_evaluation();
    check_repeated_factor_logs();
    check_warm_cache_and_second_pass();
    check_multi_factor_leftover_roots();
    check_definite_non_powers();

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
