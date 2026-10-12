#include <silex/order_unit.hpp>

#include <silex/class_group.hpp>
#include <silex/factor_base.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/residue_field.hpp>
#include <silex/signature.hpp>
#include <silex/unit.hpp>
#include <silex/zeta.hpp>

#include "class_group/class_group_certification_internal.hpp"
#include "class_group/class_group_internal.hpp"
#include "factored_element/wide_factored_internal.hpp"
#include "order_unit/order_unit_internal.hpp"
#include "order_unit/relation_unit_internal.hpp"
#include "order_unit/class_unit_transaction_internal.hpp"
#include "order_unit/compute_internal.hpp"
#include "order_unit/precision_cap_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"

#include <flint/fmpq_poly.h>
#include <flint/ulong_extras.h>

#include <cassert>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

bool set_fmpz_si(sflint::Fmpz& out, slong value) noexcept {
    sflint::fmpz_set_si(sflint::FmpzRef(out), value);
    return true;
}

int test_unit_extraction_cache_clear() {
    silex::NumberField field = silex::test::quadratic_field(2);
    silex::detail::RelationUnitExtractionState state;
    state.cutoff_inverse.emplace(1, 1);
    state.cutoff_inverse_precision = 64;
    state.cutoff_inverse_rank = 1;
    state.reduce_mod_units_inverses.emplace_back(64, 1);
    state.reduce_mod_units_log_generators.emplace_back(field);
    assert(state.reduce_mod_units_log_generators.back().one());
    state.expand_reduce_mod_units_log_generators = true;
    state.reduce_mod_units_log_generator_expansion_attempted = true;
    state.reduce_mod_units_log_generators_expanded = true;

    state.clear_dependent_unit_cache();
    assert(!state.cutoff_inverse.has_value());
    assert(state.cutoff_inverse_precision == 0);
    assert(state.cutoff_inverse_rank == 0);
    assert(state.reduce_mod_units_inverses.empty());
    assert(state.reduce_mod_units_log_generators.empty());
    assert(state.expand_reduce_mod_units_log_generators);
    assert(!state.reduce_mod_units_log_generator_expansion_attempted);
    assert(!state.reduce_mod_units_log_generators_expanded);
    return 0;
}

void set_zmat_entry_si(sflint::FmpzMat& matrix,
                       slong row,
                       slong column,
                       slong value) noexcept {
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(matrix), row, column),
            value);
}

void set_zmat_entry_str(sflint::FmpzMat& matrix,
                        slong row,
                        slong column,
                        const char* value) noexcept {
    assert(sflint::fmpz_set_str(
            sflint::fmpz_mat_entry(sflint::FmpzMatRef(matrix), row, column),
            value));
}

void set_arb_mat_entry_si_2exp(sflint::ArbMat& matrix,
                               slong row,
                               slong column,
                               slong mantissa,
                               slong exponent) noexcept {
    arb_struct* entry = arb_mat_entry(matrix.raw(), row, column);
    ::arf_set_si(arb_midref(entry), mantissa);
    ::arf_mul_2exp_si(arb_midref(entry), arb_midref(entry), exponent);
    ::mag_zero(arb_radref(entry));
}

void add_arb_mat_entry_2exp(sflint::ArbMat& matrix,
                            slong row,
                            slong column,
                            slong exponent) noexcept {
    sflint::Arf term;
    ::arf_set_ui(term.raw(), 1);
    ::arf_mul_2exp_si(term.raw(), term.raw(), exponent);
    arb_struct* entry = arb_mat_entry(matrix.raw(), row, column);
    ::arf_add(arb_midref(entry), arb_midref(entry), term.raw(),
              ARF_PREC_EXACT, ARF_RND_NEAR);
    ::mag_zero(arb_radref(entry));
}

void set_arb_mat_entry_interval_si_2exp(
        sflint::ArbMat& matrix,
        slong row,
        slong column,
        slong lower_mantissa,
        slong lower_exponent,
        slong upper_mantissa,
        slong upper_exponent) noexcept {
    sflint::Arf lower;
    sflint::Arf upper;
    ::arf_set_si(lower.raw(), lower_mantissa);
    ::arf_mul_2exp_si(lower.raw(), lower.raw(), lower_exponent);
    ::arf_set_si(upper.raw(), upper_mantissa);
    ::arf_mul_2exp_si(upper.raw(), upper.raw(), upper_exponent);
    ::arb_set_interval_arf(arb_mat_entry(matrix.raw(), row, column),
                           lower.raw(), upper.raw(), 256);
}

void set_arb_mat_entry_huge_2exp(sflint::ArbMat& matrix,
                                 slong row,
                                 slong column,
                                 const char* exponent) noexcept {
    sflint::Fmpz power;
    assert(sflint::fmpz_set_str(sflint::FmpzRef(power), exponent));
    arb_struct* entry = arb_mat_entry(matrix.raw(), row, column);
    ::arf_one(arb_midref(entry));
    ::arf_mul_2exp_fmpz(arb_midref(entry), arb_midref(entry), power.raw());
    ::mag_zero(arb_radref(entry));
}

void assert_regulator_pivot_rows(
        const sflint::ArbMat& matrix,
        const std::vector<slong>& expected) noexcept {
    std::vector<slong> actual;
    assert(silex::detail::regulator_pivot_rows_for_testing(
            actual, sflint::ArbMatConstRef(matrix), 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(actual == expected);
}

bool mat_entry_is_si(const sflint::FmpzMat& matrix,
                     slong row,
                     slong col,
                     slong value) noexcept {
    return sflint::fmpz_equal_si(
            sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(matrix), row, col),
            value);
}

bool check_unit_proof_record(const silex::OrderUnitGroup& group,
                             slong index,
                             slong expected_ell,
                             silex::ProofState expected_status,
                             slong expected_aux_bound,
                             slong expected_local_primes,
                             bool expected_changed) noexcept {
    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    silex::ProofState status = silex::ProofState::not_checked;
    slong local_primes = -1;
    bool changed = !expected_changed;
    return group.unit_proof_record(sflint::FmpzRef(ell), status,
                                   sflint::FmpzRef(aux_bound), local_primes,
                                   changed, index) &&
           sflint::fmpz_equal_si(ell, expected_ell) &&
           status == expected_status &&
           sflint::fmpz_equal_si(aux_bound, expected_aux_bound) &&
           (expected_local_primes >= 0
                    ? local_primes == expected_local_primes
                    : local_primes > 0) &&
           changed == expected_changed;
}

bool has_unit_proof_record(const silex::OrderUnitGroup& group,
                           slong expected_ell,
                           silex::ProofState expected_status,
                           bool expected_changed) noexcept {
    for (slong i = 0; i < group.unit_proof_record_count(); ++i) {
        sflint::Fmpz ell;
        sflint::Fmpz aux_bound;
        silex::ProofState status = silex::ProofState::not_checked;
        slong local_primes = 0;
        bool changed = !expected_changed;
        if (group.unit_proof_record(sflint::FmpzRef(ell), status,
                                    sflint::FmpzRef(aux_bound), local_primes,
                                    changed, i) &&
            sflint::fmpz_equal_si(ell, expected_ell) &&
            status == expected_status && local_primes > 0 &&
            changed == expected_changed) {
            return true;
        }
    }
    return false;
}

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quadratic_field(slong radicand) noexcept {
    return silex::test::quadratic_field(radicand);
}

silex::NumberField cubic_field(slong c1, slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quartic_field(slong c3,
                                 slong c2,
                                 slong c1,
                                 slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, c3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, c2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::NumberField quintic_field(slong c3,
                                 slong c2,
                                 slong c1,
                                 slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 5, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, c3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, c2);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);

    return silex::test::field_by_polynomial(
        sflint::FmpqPolyConstRef(polynomial));
}

silex::OrderUnitGroup real_quadratic_unit_group() noexcept {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);

    silex::OrderUnitGroup group;
    assert(group.compute(order));
    return group;
}

bool order_value_has_power_one(const silex::OrderElement& generator,
                               slong exponent) noexcept {
    const silex::Order* order = generator.parent();
    if (order == nullptr || order->parent() == nullptr) {
        return false;
    }

    silex::Element value(*order->parent());
    silex::Element power(*order->parent());
    sflint::Fmpz exp;
    sflint::fmpz_set_si(sflint::FmpzRef(exp), exponent);
    return generator.get_element(value) &&
           power.pow_fmpz(value, sflint::FmpzConstRef(exp)) &&
           power.equal_si(1);
}

bool check_rank_zero_group(const silex::OrderUnitGroup& group,
                           const silex::Order& order,
                           slong expected_order) noexcept {
    if (!group.is_defined() || !group.is_set() ||
        !silex::same_order_parent(group.parent(), &order) ||
        group.free_rank() != 0 ||
        group.certification_status() != silex::CertificationMode::proven) {
        return false;
    }

    sflint::Fmpz torsion_order;
    silex::OrderElement generator(order);
    sflint::Arb regulator;
    return group.torsion_order(sflint::FmpzRef(torsion_order)) &&
           sflint::fmpz_equal_si(torsion_order, expected_order) &&
           group.torsion_generator(generator) &&
           order_value_has_power_one(generator, expected_order) &&
           group.regulator(sflint::ArbRef(regulator)) &&
           sflint::arb_is_one(regulator);
}

bool set_real_quadratic_unit(silex::Element& out) noexcept {
    const silex::NumberField* field = out.parent();
    if (field == nullptr) {
        return false;
    }

    silex::Element theta(*field);
    return theta.gen() && out.add_si(theta, 1);
}

bool check_first_free_generator(const silex::OrderUnitGroup& group,
                                const silex::Element& expected) noexcept {
    const silex::NumberField* field = expected.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           expanded.equal(expected);
}

bool check_first_free_generator_square(const silex::OrderUnitGroup& group,
                                       const silex::Element& expected_square)
        noexcept {
    const silex::NumberField* field = expected_square.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    silex::Element square(*field);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           square.multiply(expanded, expanded) &&
           square.equal(expected_square);
}

bool check_first_free_generator_is_order_unit(
        const silex::OrderUnitGroup& group,
        const silex::Order& order) noexcept {
    const silex::NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    silex::OrderElement order_element(order);
    silex::Ideal principal(order);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           order_element.set_element(expanded) &&
           principal.set_principal(order_element) &&
           principal.is_one();
}

bool first_free_generator_is_pm_unit_power_pm1(
        const silex::OrderUnitGroup& group,
        const silex::Element& unit) noexcept {
    const silex::NumberField* field = unit.parent();
    if (field == nullptr) {
        return false;
    }

    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    silex::Element inverse(*field);
    silex::Element candidate(*field);
    if (!group.free_generator(compact, 0) || !compact.evaluate(expanded) ||
        !inverse.invert(unit)) {
        return false;
    }
    const silex::Element* bases[] = {&unit, &inverse};
    for (const silex::Element* base : bases) {
        if (expanded.equal(*base)) {
            return true;
        }
        if (!candidate.negate(*base)) {
            return false;
        }
        if (expanded.equal(candidate)) {
            return true;
        }
    }
    return false;
}

bool set_quadratic_coeffs(silex::Element& out,
                          slong constant_num,
                          ulong constant_den,
                          slong linear_num,
                          ulong linear_den) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;

    sflint::fmpq_set_si(coefficient, constant_num, constant_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    sflint::fmpq_set_si(coefficient, linear_num, linear_den);
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool set_integral_quadratic_coeffs_str(silex::Element& out,
                                       const char* constant,
                                       const char* linear) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::Fmpq coefficient;
    sflint::Fmpz integer;

    if (!sflint::fmpz_set_str(sflint::FmpzRef(integer), constant)) {
        return false;
    }
    sflint::fmpq_set_fmpz(coefficient, sflint::FmpzConstRef(integer));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 0, coefficient);
    if (!sflint::fmpz_set_str(sflint::FmpzRef(integer), linear)) {
        return false;
    }
    sflint::fmpq_set_fmpz(coefficient, sflint::FmpzConstRef(integer));
    sflint::fmpq_poly_set_coeff_fmpq(polynomial, 1, coefficient);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
}

bool degree_one_residue_root(sflint::Fmpz& root,
                             const silex::PrimeIdeal& prime) noexcept {
    sflint::Fmpz p;
    sflint::FmpzPoly residue;
    sflint::Fmpz constant;
    sflint::Fmpz leading;
    sflint::Fmpz inverse_leading;
    if (prime.residue_degree() != 1 ||
        !prime.rational_prime(sflint::FmpzRef(p)) ||
        !prime.residue_polynomial(sflint::FmpzPolyRef(residue)) ||
        ::fmpz_poly_degree(residue.raw()) != 1) {
        return false;
    }

    ::fmpz_poly_get_coeff_fmpz(constant.raw(), residue.raw(), 0);
    ::fmpz_poly_get_coeff_fmpz(leading.raw(), residue.raw(), 1);
    ::fmpz_mod(leading.raw(), leading.raw(), p.raw());
    if (::fmpz_invmod(inverse_leading.raw(), leading.raw(), p.raw()) == 0) {
        return false;
    }
    ::fmpz_neg(root.raw(), constant.raw());
    ::fmpz_mod(root.raw(), root.raw(), p.raw());
    ::fmpz_mul(root.raw(), root.raw(), inverse_leading.raw());
    ::fmpz_mod(root.raw(), root.raw(), p.raw());
    return true;
}

bool check_real_quadratic_group(const silex::OrderUnitGroup& group,
                                const silex::Order& order,
                                const silex::Element& expected_unit) noexcept {
    if (!group.is_defined() || !group.is_set() ||
        !silex::same_order_parent(group.parent(), &order) ||
        group.free_rank() != 1 ||
        group.certification_status() != silex::CertificationMode::proven ||
        !check_first_free_generator(group, expected_unit)) {
        return false;
    }

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    return group.torsion_order(sflint::FmpzRef(torsion_order)) &&
           sflint::fmpz_equal_si(torsion_order, 2) &&
           group.regulator(sflint::ArbRef(regulator)) &&
           sflint::arb_is_positive(regulator);
}

struct TorsionEventCounter {
    const char* label = nullptr;
    slong count = 0;
};

void torsion_event_callback(void* user,
                            silex::DiagnosticsModule,
                            silex::ProfileEvent event,
                            const char*,
                            const char* label) noexcept {
    TorsionEventCounter* counter = static_cast<TorsionEventCounter*>(user);
    if (event == silex::ProfileEvent::event && label != nullptr &&
        std::strcmp(label, counter->label) == 0) {
        ++counter->count;
    }
}

// Forces the shared precision-doubling cap (T-111) down to one bit so each
// capped doubling loop must stop after its first undecided pass and return
// false, then checks the same call succeeds with the default cap.
struct PrecisionCapOverride {
    explicit PrecisionCapOverride(slong cap) noexcept {
        silex::detail::set_precision_doubling_cap_for_testing(cap);
    }
    ~PrecisionCapOverride() {
        silex::detail::set_precision_doubling_cap_for_testing(0);
    }
};

int test_precision_doubling_cap_default_value() {
    using silex::detail::precision_doubling_cap;
    assert(precision_doubling_cap(1) == (WORD(1) << 20));
    assert(precision_doubling_cap(1 << 14) == (WORD(1) << 20));
    assert(precision_doubling_cap(1 << 15) == (WORD(1) << 21));
    assert(precision_doubling_cap(0) == (WORD(1) << 20));
    assert(precision_doubling_cap(std::numeric_limits<slong>::max()) ==
           (WORD(1) << 20));
    assert(silex::detail::precision_doubling_allowed(
            (WORD(1) << 19), WORD(1) << 20, nullptr, "x"));
    assert(!silex::detail::precision_doubling_allowed(
            (WORD(1) << 19) + 1, WORD(1) << 20, nullptr, "x"));
    return 0;
}

int test_precision_doubling_cap_forced_failures() {
    silex::NumberField field = quadratic_field(2);
    silex::Element epsilon(field);
    assert(set_quadratic_coeffs(epsilon, 1, 1, 1, 1));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::EmbeddingContext embeddings(field);
    silex::FactoredElementSpan span(&generator, 1);

    // Every capped loop emits one profile event when the cap stops it.
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    TorsionEventCounter counter;
    counter.label = "unit_group.precision_doubling_cap";
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            torsion_event_callback, &counter);
    slong expected_events = 0;
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
#define SILEX_TEST_CAP_EVENT() (++expected_events)
#else
#define SILEX_TEST_CAP_EVENT() ((void)0)
#endif

    // Factored logarithmic embedding: in a field whose two complex pairs
    // agree to about 300 bits, the embedding is unavailable at the first
    // working precision, so the loop must double; the cap stops it.
    {
        silex::NumberField close = silex::test::close_imaginary_pairs_field();
        silex::Element gen(close);
        assert(gen.gen());
        silex::FactoredElement factored(close);
        assert(factored.set_element(gen));
        silex::EmbeddingContext close_embeddings(close);
        sflint::ArbVec logs(2);
        {
            PrecisionCapOverride cap(1);
            assert(!factored.logarithmic_embedding(
                    sflint::ArbVecRef(logs), close_embeddings,
                    silex::LogEmbeddingMode::product, 1, &diagnostics));
            SILEX_TEST_CAP_EVENT();
            assert(counter.count == expected_events);
        }
        assert(factored.logarithmic_embedding(
                sflint::ArbVecRef(logs), close_embeddings,
                silex::LogEmbeddingMode::product, 1));
    }

    // Compact independence: precision 1 is not enough to decide.
    {
        bool independent = false;
        {
            PrecisionCapOverride cap(1);
            assert(!silex::detail::compact_independent(
                    independent, embeddings, span, 1, &diagnostics));
            SILEX_TEST_CAP_EVENT();
            assert(counter.count == expected_events);
        }
        assert(silex::detail::compact_independent(
                independent, embeddings, span, 1));
        assert(independent);
    }

    // Compact regulator to a tolerance no finite precision meets quickly.
    {
        sflint::Arb regulator;
        {
            PrecisionCapOverride cap(1);
            assert(!silex::detail::compact_regulator_adaptive(
                    sflint::ArbRef(regulator), embeddings, span, 32,
                    1000000, &diagnostics));
            SILEX_TEST_CAP_EVENT();
            assert(counter.count == expected_events);
        }
        assert(silex::detail::compact_regulator_adaptive(
                sflint::ArbRef(regulator), embeddings, span, 32, 64));
    }

    // Conjugate-log cutoff inverse: in the cyclic cubic x^3 - 3x + 1 the
    // units theta and theta - 1 give a rank-two group whose 2x2 cutoff
    // matrix cannot be inverted at one bit.
    {
        silex::NumberField cubic = cubic_field(-3, 1);
        silex::Order order;
        order = silex::test::equation_order(cubic);
        silex::Element theta(cubic);
        silex::Element theta_minus_one(cubic);
        assert(theta.gen());
        assert(theta_minus_one.add_si(theta, -1));
        silex::FactoredElement generators[] = {
                silex::FactoredElement(cubic), silex::FactoredElement(cubic)};
        assert(generators[0].set_element(theta));
        assert(generators[1].set_element(theta_minus_one));
        silex::EmbeddingContext cubic_embeddings(cubic);
        silex::OrderUnitGroup group;
        assert(group.set_units(order,
                               silex::FactoredElementSpan(generators, 2),
                               cubic_embeddings, 64));
        assert(group.free_rank() == 2);
        silex::detail::RelationUnitExtractionState state;
        const sflint::ArbMat* inverse = nullptr;
        slong precision = 1;
        group.set_diagnostics(&diagnostics);
        {
            PrecisionCapOverride cap(1);
            assert(!silex::detail::conj_log_cutoff_inverse(
                    inverse, state, group, cubic_embeddings, precision));
            SILEX_TEST_CAP_EVENT();
            assert(counter.count == expected_events);
        }
        group.set_diagnostics(nullptr);
        precision = 1;
        assert(silex::detail::conj_log_cutoff_inverse(
                inverse, state, group, cubic_embeddings, precision));
        assert(inverse != nullptr);
    }
#undef SILEX_TEST_CAP_EVENT
    return 0;
}

// An exact half-integer coordinate stops the reduce-mod-units precision
// ladder.  In the cyclic cubic x^3 - 3x + 1 the units theta^2 and theta - 1
// generate a rank-two group in which theta^3 has coordinates exactly
// (3/2, 0).  No precision certifies a rounding of 3/2, and the zero
// coordinate defeats the approximate rounding, so the ladder used to double
// the precision to its cap before returning the value unreduced.  Once the
// coordinate ball is narrow (radius at most 0.1) and contains the
// half-integer, and the zero coordinate's ball is narrow too, it now returns
// the value unreduced at once.
int test_reduce_mod_units_half_integer_tie() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order order;
    order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    silex::Element theta_minus_one(cubic);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));
    silex::FactoredElement generators[] = {
            silex::FactoredElement(cubic), silex::FactoredElement(cubic)};
    assert(generators[0].push(theta, 2));
    assert(generators[1].set_element(theta_minus_one));
    silex::EmbeddingContext embeddings(cubic);
    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 64));
    assert(group.free_rank() == 2);
    const slong precision_before = embeddings.precision();

    std::vector<silex::FactoredElement> values;
    values.emplace_back(cubic);
    assert(values[0].push(theta, 3));
    assert(silex::detail::reduce_relation_units_modulo(
            values, group, embeddings, 32));

    // The value is returned unreduced: still theta^3.
    silex::Element theta_cubed(cubic);
    silex::Element value(cubic);
    silex::Element theta_squared(cubic);
    assert(theta_squared.multiply(theta, theta));
    assert(theta_cubed.multiply(theta_squared, theta));
    assert(values[0].evaluate(value));
    assert(value.equal(theta_cubed));

    // The ladder stopped well below its cap of
    // kReduceModUnitsMaxPrecision bits.
    assert(embeddings.precision() <= FLINT_MAX(precision_before, 1024));
    return 0;
}

// The reduction step value * prod_j g_j^(-c_j) as it was first written:
// copy the generator, raise it to -c_j, multiply and copy back.  The test
// below pins apply_unit_reduction to the factor list this builds, including
// its two exits: false when a scaled factor exponent overflows slong (as
// pow_si reports), and true with the values unchanged when -c_j itself does
// not fit an slong.
bool reference_apply_unit_reduction(
        std::vector<silex::FactoredElement>& values,
        const silex::OrderUnitGroup& group,
        const sflint::FmpzMat& coefficients) noexcept {
    const silex::NumberField* field = group.parent()->parent();
    const slong rank = group.free_rank();
    std::vector<silex::FactoredElement> reduced;
    for (std::size_t i = 0; i < values.size(); ++i) {
        reduced.emplace_back(*field);
        if (!reduced.back().set(values[i])) {
            return false;
        }
        for (slong j = 0; j < rank; ++j) {
            const sflint::FmpzConstRef coefficient = sflint::fmpz_mat_entry(
                    coefficients, static_cast<slong>(i), j);
            if (sflint::fmpz_is_zero(coefficient)) {
                continue;
            }
            sflint::Fmpz exponent;
            sflint::fmpz_neg(sflint::FmpzRef(exponent), coefficient);
            if (!sflint::fmpz_fits_si(sflint::FmpzConstRef(exponent))) {
                return true;
            }
            silex::FactoredElement unit(*field);
            silex::FactoredElement power(*field);
            silex::FactoredElement product(*field);
            if (!group.free_generator(unit, j) ||
                !power.pow_si(unit, sflint::fmpz_get_si(
                                            sflint::FmpzConstRef(exponent))) ||
                !product.multiply(reduced.back(), power) ||
                !reduced.back().set(product)) {
                return false;
            }
        }
    }
    values = std::move(reduced);
    return true;
}

bool same_factor_lists(const std::vector<silex::FactoredElement>& left,
                       const std::vector<silex::FactoredElement>& right) {
    if (left.size() != right.size()) {
        return false;
    }
    for (std::size_t i = 0; i < left.size(); ++i) {
        const auto a = left[i].factors();
        const auto b = right[i].factors();
        if (a.size() != b.size()) {
            return false;
        }
        for (std::size_t k = 0; k < a.size(); ++k) {
            if (a[k].exponent != b[k].exponent ||
                !a[k].factor.equal(b[k].factor)) {
                return false;
            }
        }
    }
    return true;
}

std::vector<silex::FactoredElement> copy_factored_values(
        const std::vector<silex::FactoredElement>& values) {
    std::vector<silex::FactoredElement> out;
    for (const auto& value : values) {
        out.emplace_back(*value.parent());
        assert(out.back().set(value));
    }
    return out;
}

// apply_unit_reduction appends each generator's factors, with exponent
// -c_j * e, to the reduced value in place.  It must build exactly the
// factor list of the copy/pow/multiply/set path above (no merging of equal
// factors), and keep that path's failure and early-exit behaviour.
int test_apply_unit_reduction_factor_list() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order order;
    order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    silex::Element theta_minus_one(cubic);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));
    // Multi-factor generators with a repeated factor: theta^2 (theta-1) and
    // (theta-1)^-1 theta^3 (theta-1)^2, of log coordinates (2,1) and (3,1).
    silex::FactoredElement generators[] = {
            silex::FactoredElement(cubic), silex::FactoredElement(cubic)};
    assert(generators[0].push(theta, 2));
    assert(generators[0].push(theta_minus_one, 1));
    assert(generators[1].push(theta_minus_one, -1));
    assert(generators[1].push(theta, 3));
    assert(generators[1].push(theta_minus_one, 2));
    silex::EmbeddingContext embeddings(cubic);
    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 64));
    assert(group.free_rank() == 2);
    slong largest = 0;
    bool repeated = false;
    for (slong j = 0; j < 2; ++j) {
        silex::FactoredElement generator(cubic);
        assert(group.free_generator(generator, j));
        assert(generator.length() >= 2);
        for (const auto& entry : generator.factors()) {
            const slong magnitude =
                    entry.exponent < 0 ? -entry.exponent : entry.exponent;
            largest = magnitude > largest ? magnitude : largest;
        }
        // Some free generator holds two equal factors (no merging).
        const auto fs = generator.factors();
        for (slong a = 0; a < generator.length(); ++a) {
            for (slong b = a + 1; b < generator.length(); ++b) {
                if (fs[a].factor.equal(fs[b].factor)) repeated = true;
            }
        }
    }
    assert(largest >= 2);
    assert(repeated);

    std::vector<silex::FactoredElement> values;
    values.emplace_back(cubic);
    assert(values[0].push(theta, 5));
    assert(values[0].push(theta_minus_one, -2));
    values.emplace_back(cubic);
    assert(values[1].push(theta_minus_one, 3));
    values.emplace_back(cubic);
    assert(values[2].push(theta, 1));

    // Mixed signs and a zero row and entry.
    sflint::FmpzMat coefficients(3, 2);
    set_zmat_entry_si(coefficients, 0, 0, 1);
    set_zmat_entry_si(coefficients, 0, 1, -2);
    set_zmat_entry_si(coefficients, 2, 0, -3);
    set_zmat_entry_si(coefficients, 2, 1, 4);
    std::vector<silex::FactoredElement> actual = copy_factored_values(values);
    std::vector<silex::FactoredElement> expected =
            copy_factored_values(values);
    assert(silex::detail::apply_unit_reduction(actual, group, coefficients));
    assert(reference_apply_unit_reduction(expected, group, coefficients));
    assert(same_factor_lists(actual, expected));
    assert(actual[0].length() > values[0].length());

    // A scaled exponent past slong fails, as pow_si does, with the values
    // unchanged even though row 0 was already reduced.
    sflint::FmpzMat overflow(3, 2);
    set_zmat_entry_si(overflow, 0, 0, 1);
    set_zmat_entry_si(overflow, 1, 0, -(WORD(1) << 62));
    set_zmat_entry_si(overflow, 1, 1, -(WORD(1) << 62));
    actual = copy_factored_values(values);
    expected = copy_factored_values(values);
    assert(!reference_apply_unit_reduction(expected, group, overflow));
    assert(!silex::detail::apply_unit_reduction(actual, group, overflow));
    assert(same_factor_lists(actual, values));
    assert(same_factor_lists(expected, values));

    // -c_j that does not fit an slong returns true with the values
    // unchanged.
    sflint::FmpzMat unfit(3, 2);
    set_zmat_entry_si(unfit, 0, 1, 1);
    set_zmat_entry_si(unfit, 1, 0, std::numeric_limits<slong>::min());
    actual = copy_factored_values(values);
    expected = copy_factored_values(values);
    assert(reference_apply_unit_reduction(expected, group, unfit));
    assert(silex::detail::apply_unit_reduction(actual, group, unfit));
    assert(same_factor_lists(actual, values));
    assert(same_factor_lists(expected, values));

    // A coefficient matrix of the wrong shape is rejected.
    sflint::FmpzMat wrong(2, 2);
    actual = copy_factored_values(values);
    assert(!silex::detail::apply_unit_reduction(actual, group, wrong));
    assert(same_factor_lists(actual, values));
    return 0;
}

bool arb_vecs_identical(const sflint::ArbVec& left,
                        const sflint::ArbVec& right) noexcept {
    if (left.length() != right.length()) {
        return false;
    }
    for (slong i = 0; i < left.length(); ++i) {
        if (!arb_equal(left.data() + i, right.data() + i)) {
            return false;
        }
    }
    return true;
}

// The factored log memo returns exactly what a fresh
// FactoredElement::logarithmic_embedding call computes.  It stores a result
// only under a known embedding-cache version, returns it only while that
// version holds, never confuses two candidates with the same exponents, and
// recomputes once a factor log in the cache has been upgraded.
int test_factored_log_memo_matches_direct_call() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Element theta(cubic);
    silex::Element theta_minus_one(cubic);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));

    silex::FactoredElement value(cubic);
    assert(value.push(theta, 2));
    assert(value.push(theta_minus_one, -1));
    assert(value.push(theta, 1));
    // Same length and exponents, different factors.
    silex::FactoredElement other(cubic);
    assert(other.push(theta_minus_one, 2));
    assert(other.push(theta, -1));
    assert(other.push(theta_minus_one, 1));

    silex::EmbeddingContext embeddings(cubic);
    const slong places = 3;
    const slong precision = 64;
    const auto mode = silex::LogEmbeddingMode::product;
    std::uint64_t version = 0;
    assert(!silex::detail::embedding_log_cache_version(version, embeddings));

    silex::detail::FactoredLogMemo memo;
    sflint::ArbVec first(places);
    sflint::ArbVec memoized(places);
    sflint::ArbVec direct(places);

    // No cache version before the first call, so nothing is stored.
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(first), &memo, value, embeddings, mode,
            precision));
    assert(memo.stores == 0 && memo.hits == 0);
    std::uint64_t version_before = 0;
    assert(silex::detail::embedding_log_cache_version(version_before,
                                                      embeddings));

    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(memoized), &memo, value, embeddings, mode,
            precision));
    assert(memo.stores == 1 && memo.hits == 0);
    // Two distinct factors, each stored once.
    assert(memo.factors.size() == 2);
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(memoized), &memo, value, embeddings, mode,
            precision));
    assert(memo.hits == 1);
    assert(value.logarithmic_embedding(sflint::ArbVecRef(direct), embeddings,
                                       mode, precision));
    assert(arb_vecs_identical(memoized, direct));
    assert(arb_vecs_identical(first, direct));

    // A different precision or mode is a different result.
    sflint::ArbVec plain(places);
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(plain), &memo, value, embeddings,
            silex::LogEmbeddingMode::plain, precision));
    assert(memo.hits == 1);

    // Same exponents, different factors: no hit.
    sflint::ArbVec other_memoized(places);
    sflint::ArbVec other_direct(places);
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(other_memoized), &memo, other, embeddings,
            mode, precision));
    assert(memo.hits == 1);
    assert(other.logarithmic_embedding(sflint::ArbVecRef(other_direct),
                                       embeddings, mode, precision));
    assert(arb_vecs_identical(other_memoized, other_direct));
    assert(!arb_vecs_identical(other_memoized, memoized));
    // New cache entries do not change the version.
    assert(silex::detail::embedding_log_cache_version(version, embeddings));
    assert(version == version_before);

    // Upgrading a cached factor log changes the version; the memo then
    // recomputes, matching a fresh call under the upgraded cache.
    sflint::ArbVec factor_logs(places);
    assert(silex::logarithmic_embedding(sflint::ArbVecRef(factor_logs),
                                        embeddings, theta, mode, 4096));
    assert(silex::detail::embedding_log_cache_version(version, embeddings));
    assert(version != version_before);
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(memoized), &memo, value, embeddings, mode,
            precision));
    assert(memo.hits == 1);
    assert(value.logarithmic_embedding(sflint::ArbVecRef(direct), embeddings,
                                       mode, precision));
    assert(arb_vecs_identical(memoized, direct));
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(memoized), &memo, value, embeddings, mode,
            precision));
    assert(memo.hits == 2);
    assert(arb_vecs_identical(memoized, direct));

    // Without a memo the call is the direct call.
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(memoized), nullptr, value, embeddings, mode,
            precision));
    assert(arb_vecs_identical(memoized, direct));
    return 0;
}

// A factor log that fails at the starting precision and is stored only at
// a doubled precision changes what a fresh call computes: the fresh call
// then hits the doubled-precision entry and finishes at the starting
// precision.  The memo must not return the result of the first call.
int test_factored_log_memo_after_failed_factor_log() {
    silex::NumberField field = quadratic_field(2);
    silex::Element theta(field);
    silex::Element unit(field);
    silex::Element factor(field);
    assert(theta.gen());
    assert(unit.add_si(theta, -1));
    sflint::Fmpz exponent;
    fmpz_set_si(exponent.raw(), 40);
    assert(factor.pow_fmpz(unit, sflint::FmpzConstRef(exponent)));

    silex::FactoredElement value(field);
    assert(value.push(factor, 1));

    silex::EmbeddingContext embeddings(field);
    const slong places = 2;
    const slong precision = 16;
    const auto mode = silex::LogEmbeddingMode::product;
    // Create the cache, so the memo has a version to store under.
    sflint::ArbVec theta_logs(places);
    assert(silex::logarithmic_embedding(sflint::ArbVecRef(theta_logs),
                                        embeddings, theta, mode, 4096));
    // The factor's conjugate (sqrt(2) - 1)^40 is too small to separate from
    // 0 at the starting precision.
    sflint::ArbVec factor_logs(places);
    assert(!silex::logarithmic_embedding(sflint::ArbVecRef(factor_logs),
                                         embeddings, factor, mode, precision));

    silex::detail::FactoredLogMemo memo;
    sflint::ArbVec first(places);
    sflint::ArbVec memoized(places);
    sflint::ArbVec direct(places);
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(first), &memo, value, embeddings, mode,
            precision));
    assert(silex::detail::memoized_factored_log_embedding(
            sflint::ArbVecRef(memoized), &memo, value, embeddings, mode,
            precision));
    assert(value.logarithmic_embedding(sflint::ArbVecRef(direct), embeddings,
                                       mode, precision));
    assert(arb_vecs_identical(memoized, direct));
    return 0;
}

// round_arb_approx_fmpz ports the reference approximate rounding.  Its width
// check computes floor(lower) - ceil(upper), never positive, so it never
// rejects: a wide ball rounds to div(ceil(upper) + floor(lower), 2).  The
// rounding fails only through Silex's finiteness guard on sqrt(|value|),
// which fires for a ball with 0 strictly inside it (the reference would
// return 0 there); an exact 0 or a ball with endpoint 0 still rounds.
int test_round_arb_approx_width_check() {
    const slong precision = 64;
    sflint::Fmpz out;

    // 100 +/- 48: wide (radius far above sqrt(100)), excludes 0.  Rounds
    // to div(148 + 52, 2) = 100.
    sflint::Arb wide;
    arb_set_ui(wide.raw(), 100);
    arb_add_error_2exp_si(wide.raw(), 5);
    arb_add_error_2exp_si(wide.raw(), 4);
    assert(silex::detail::round_arb_approx_fmpz(sflint::FmpzRef(out),
                                                wide.raw(), precision));
    assert(fmpz_get_si(out.raw()) == 100);

    // Narrow ball excluding 0.
    sflint::Arb narrow;
    arb_set_si(narrow.raw(), -7);
    arb_add_error_2exp_si(narrow.raw(), -20);
    assert(silex::detail::round_arb_approx_fmpz(sflint::FmpzRef(out),
                                                narrow.raw(), precision));
    assert(fmpz_get_si(out.raw()) == -7);

    // 0 +/- 1/2 has 0 strictly inside: sqrt(|value|) is indeterminate and
    // the guard rejects it.
    sflint::Arb zero_ball;
    arb_zero(zero_ball.raw());
    arb_add_error_2exp_si(zero_ball.raw(), -1);
    assert(!silex::detail::round_arb_approx_fmpz(sflint::FmpzRef(out),
                                                 zero_ball.raw(), precision));

    // 1/4 +/- 1/2 also has 0 strictly inside and is rejected although its
    // midpoint is positive.
    sflint::Arb straddle;
    arb_set_ui(straddle.raw(), 1);
    arb_mul_2exp_si(straddle.raw(), straddle.raw(), -2);
    arb_add_error_2exp_si(straddle.raw(), -1);
    assert(!silex::detail::round_arb_approx_fmpz(sflint::FmpzRef(out),
                                                 straddle.raw(), precision));

    // An exact 0 rounds to 0.
    sflint::Arb exact_zero;
    arb_zero(exact_zero.raw());
    fmpz_set_si(out.raw(), 5);
    assert(silex::detail::round_arb_approx_fmpz(sflint::FmpzRef(out),
                                                exact_zero.raw(), precision));
    assert(fmpz_is_zero(out.raw()));

    // 1/2 +/- 1/2 = [0, 1] touches 0 only at an endpoint, so sqrt(|value|)
    // stays finite and it rounds to div(1 + 0, 2) = 0.
    sflint::Arb touching;
    arb_set_ui(touching.raw(), 1);
    arb_mul_2exp_si(touching.raw(), touching.raw(), -1);
    arb_add_error_2exp_si(touching.raw(), -1);
    fmpz_set_si(out.raw(), 5);
    assert(silex::detail::round_arb_approx_fmpz(sflint::FmpzRef(out),
                                                touching.raw(), precision));
    assert(fmpz_is_zero(out.raw()));
    return 0;
}

// The half-integer stop also requires every coordinate that defeats the
// approximate rounding to be a narrow ball.  An exact tie 3/2 beside a small
// nonzero coordinate whose wide ball still contains 0 keeps the ladder going:
// more precision can exclude 0 from that ball, and then the approximate
// rounding succeeds and reduces the value.
int test_reduce_mod_units_half_integer_tie_stop_rule() {
    const slong precision = 64;
    sflint::ArbMat rows(1, 2);
    arb_struct* tie = arb_mat_entry(rows.raw(), 0, 0);
    arb_struct* other = arb_mat_entry(rows.raw(), 0, 1);
    arb_set_ui(tie, 3);
    arb_mul_2exp_si(tie, tie, -1);
    arb_add_error_2exp_si(tie, -30);

    // 1/32 +/- 1/4: a nonzero coordinate whose ball contains 0 and has
    // radius above 0.1.  The ladder must continue.
    arb_set_ui(other, 1);
    arb_mul_2exp_si(other, other, -5);
    arb_add_error_2exp_si(other, -2);
    assert(!silex::detail::coordinate_rows_stop_at_half_integer_tie(
            rows, precision));

    // 1/32 +/- 1/16: the same coordinate as a narrow ball (radius at most
    // 0.1).  The ladder stops; this near-zero case is left unreduced.
    arb_set_ui(other, 1);
    arb_mul_2exp_si(other, other, -5);
    arb_add_error_2exp_si(other, -4);
    assert(silex::detail::coordinate_rows_stop_at_half_integer_tie(
            rows, precision));

    // 0 +/- 2^-30: an exactly zero coordinate beside the tie.  The ladder
    // stops.
    arb_zero(other);
    arb_add_error_2exp_si(other, -30);
    assert(silex::detail::coordinate_rows_stop_at_half_integer_tie(
            rows, precision));

    // Without a half-integer tie the ladder continues.
    arb_set_ui(tie, 1);
    arb_add_error_2exp_si(tie, -30);
    assert(!silex::detail::coordinate_rows_stop_at_half_integer_tie(
            rows, precision));
    return 0;
}

int test_degree_one_rank_zero_group() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderUnitGroup group;
    assert(!group.is_defined());
    assert(!group.is_set());
    assert(group.free_rank() == -1);
    assert(group.certification_status() == silex::CertificationMode::unknown);

    assert(group.compute(order));
    assert(check_rank_zero_group(group, order, 2));

    silex::OrderElement generator(order);
    assert(group.torsion_generator(generator));
    assert(generator.equal_si(-1));
    return 0;
}

int test_compute_real_quadratic_proven() {
    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Order sqrt2_order;
    sqrt2_order = silex::test::equation_order(sqrt2);
    assert(sqrt2_order.is_maximal());
    silex::Element epsilon2(sqrt2);
    assert(set_quadratic_coeffs(epsilon2, 1, 1, 1, 1));

    silex::OrderUnitGroup group;
    assert(group.compute(sqrt2_order));
    assert(check_real_quadratic_group(group, sqrt2_order, epsilon2));

    silex::NumberField sqrt5 = quadratic_field(5);
    silex::Order equation;
    silex::Order maximal(sqrt5);
    equation = silex::test::equation_order(sqrt5);
    assert(!equation.is_maximal());
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    silex::Element epsilon5(sqrt5);
    assert(set_quadratic_coeffs(epsilon5, 1, 2, 1, 2));

    silex::OrderUnitGroup half_integral;
    assert(half_integral.compute(maximal));
    assert(check_real_quadratic_group(half_integral, maximal, epsilon5));

    silex::NumberField sqrt17345 = quadratic_field(17345);
    silex::Order equation17345;
    silex::Order maximal17345(sqrt17345);
    equation17345 = silex::test::equation_order(sqrt17345);
    assert(!equation17345.is_maximal());
    assert(maximal17345.maximal_order(equation17345));
    assert(maximal17345.is_maximal());
    silex::Element epsilon17345(sqrt17345);
    assert(set_integral_quadratic_coeffs_str(
            epsilon17345, "301977958012", "2292915721"));

    silex::OrderUnitGroup large_half_integral;
    assert(large_half_integral.compute(maximal17345));
    assert(check_real_quadratic_group(
            large_half_integral, maximal17345, epsilon17345));

    return 0;
}

int test_compute_with_class_group_rank_zero() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 64;
    options.max_relations = 16;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 80));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(class_group.relation_saturation_record_count() == 0);
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(!class_group.try_certify_with_units(
            units, silex::CertificationMode::grh, 80));
    assert(check_rank_zero_group(units, order, 2));

    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 1));

    silex::ClassGroupContext proven_class_group;
    silex::OrderUnitGroup proven_units;
    options.requested_certification = silex::CertificationMode::proven;
    assert(proven_units.compute_with_class_group(
            proven_class_group, order, sflint::FmpzConstRef(bound),
            options, 80));
    assert(proven_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.relation_saturation_record_count() == 0);
    assert(proven_class_group.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.unit_proof_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(check_rank_zero_group(proven_units, order, 2));

    return 0;
}

int test_compute_with_class_group_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 160));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(units.is_set());
    assert(silex::same_order_parent(units.parent(), &order));
    assert(units.free_rank() == 1);
    silex::Element epsilon(field);
    assert(set_quadratic_coeffs(epsilon, 1, 1, 1, 1));
    assert(check_real_quadratic_group(units, order, epsilon));

    silex::ClassGroupContext proven_class_group;
    silex::OrderUnitGroup proven_units;
    options.requested_certification = silex::CertificationMode::proven;
    assert(proven_units.compute_with_class_group(
            proven_class_group, order, sflint::FmpzConstRef(bound),
            options, 160));
    assert(proven_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.unit_proof_status() ==
           silex::ProofState::verified);
    assert(proven_class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(check_real_quadratic_group(proven_units, order, epsilon));

    silex::NumberField sqrt5 = quadratic_field(5);
    silex::Order equation;
    silex::Order maximal(sqrt5);
    equation = silex::test::equation_order(sqrt5);
    assert(maximal.maximal_order(equation));
    silex::Element epsilon5(sqrt5);
    assert(set_quadratic_coeffs(epsilon5, 1, 2, 1, 2));

    sflint::Fmpz sqrt5_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(sqrt5_bound), maximal));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(sqrt5_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(sqrt5_bound), 2);
    }
    silex::ClassGroupComputeOptions sqrt5_options;
    sqrt5_options.max_candidates = 256;
    sqrt5_options.max_relations = 48;
    sqrt5_options.requested_certification =
            silex::CertificationMode::proven;

    silex::ClassGroupContext sqrt5_class_group;
    silex::OrderUnitGroup sqrt5_units;
    assert(sqrt5_units.compute_with_class_group(
            sqrt5_class_group, maximal, sflint::FmpzConstRef(sqrt5_bound),
            sqrt5_options, 160));
    assert(sqrt5_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(sqrt5_class_group.relation_kernel_unit_count() >=
           sqrt5_units.free_rank());
    assert(sqrt5_class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(check_real_quadratic_group(sqrt5_units, maximal, epsilon5));

    return 0;
}

// A proven compute on a maximal order passes its relation-kernel units
// through the index-bounded saturation pass
// (set_relation_kernel_units_index_bounded_saturated) before
// prove_index_bound, and that pass includes the torsion generator when
// ell | w.  In Q(sqrt97), w = 2, so the ell = 2 pass runs with the torsion
// column.  This pins the published generator, eps = 5604 + 569 sqrt97
// (norm -1; reference GP 2.17.4: quadunit(97) = 5035 + 1138*w with
// w = (1 + sqrt97)/2), and checks that the same class-group context replayed
// through the public saturation pass has index bound at least 2 and is
// changed by it to +-eps^(+-1), so drift in that pass or the proof shows up.
int test_compute_with_class_group_proven_saturated_generator() {
    silex::NumberField field = quadratic_field(97);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::ClassGroupComputeOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 160));
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);

    silex::Element epsilon(field);
    assert(set_quadratic_coeffs(epsilon, 5604, 1, 569, 1));
    assert(check_real_quadratic_group(units, order, epsilon));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup initial(order);
    assert(initial.set_relation_kernel_units_index_bounded(
            order, class_group, embeddings, 160, 160));
    sflint::Fmpz index_bound;
    assert(initial.regulator_index_bound(sflint::FmpzRef(index_bound), 160));
    assert(sflint::fmpz_cmp_ui(sflint::FmpzConstRef(index_bound), 2) >= 0);
    // The relation-kernel candidates come from an LLL-reduced kernel basis,
    // and here the best candidate is already a fundamental unit even though
    // the regulator bound alone does not show it.
    assert(first_free_generator_is_pm_unit_power_pm1(initial, epsilon));

    sflint::Fmpz aux_start;
    sflint::Fmpz aux_max;
    assert(set_fmpz_si(aux_start, 2));
    assert(set_fmpz_si(aux_max, 31));
    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    bool stable = false;
    assert(saturated.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, class_group, embeddings, 160, 160, 1,
            sflint::FmpzConstRef(aux_start), sflint::FmpzConstRef(aux_max),
            2));
    // Saturation finds no proper superlattice of an already saturated
    // group, so the generator is kept.
    assert(!changed);
    assert(first_free_generator_is_pm_unit_power_pm1(saturated, epsilon));

    return 0;
}

int test_compute_with_class_group_real_quadratic_210() {
    silex::NumberField field = quadratic_field(210);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    silex::Element expected_unit(field);
    assert(set_quadratic_coeffs(expected_unit, 29, 1, 2, 1));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 10000;
    options.max_relations = 1000;
    options.requested_certification = silex::CertificationMode::proven;

    auto check_pair = [&](const silex::ClassGroupContext& class_group,
                          const silex::OrderUnitGroup& units) {
        assert(class_group.has_presentation());
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.regulator_proof_status() ==
               silex::ProofState::verified);

        sflint::Fmpz class_order;
        sflint::Fmpz invariant;
        assert(class_group.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 4));
        assert(class_group.invariant_count() == 2);
        assert(class_group.invariant(sflint::FmpzRef(invariant), 0));
        assert(sflint::fmpz_equal_si(invariant, 2));
        assert(class_group.invariant(sflint::FmpzRef(invariant), 1));
        assert(sflint::fmpz_equal_si(invariant, 2));
        assert(check_real_quadratic_group(units, order, expected_unit));
    };

    auto compute_pair = [&](silex::ClassGroupContext& class_group,
                            silex::OrderUnitGroup& units,
                            slong requested_bound, slong precision) {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, requested_bound));
        assert(units.compute_with_class_group(
                class_group, order, sflint::FmpzConstRef(bound), options,
                precision));
        check_pair(class_group, units);
    };

    const slong requested_bounds[] = {10, 19, 20, 50};
    for (const slong requested_bound : requested_bounds) {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        compute_pair(class_group, units, requested_bound, 128);
    }

    silex::ClassGroupContext preserved_class_group;
    silex::OrderUnitGroup preserved_units;
    compute_pair(preserved_class_group, preserved_units, 20, 192);

    const slong relation_count_before =
            preserved_class_group.relation_count();
    const auto relations_before = preserved_class_group.relations();
    const auto invariant_generators_before =
            preserved_class_group.invariant_generator_matrix();
    const auto regulator_before = preserved_units.regulator();
    const auto torsion_order_before = preserved_units.torsion_order();
    assert(relations_before.has_value());
    assert(invariant_generators_before.has_value());
    assert(regulator_before.has_value());
    assert(torsion_order_before.has_value());

    silex::ClassGroupComputeOptions constrained_options = options;
    constrained_options.max_candidates = 0;
    constrained_options.max_relations = 0;
    sflint::Fmpz constrained_bound;
    assert(set_fmpz_si(constrained_bound, 10));
    assert(!preserved_units.compute_with_class_group(
            preserved_class_group, order,
            sflint::FmpzConstRef(constrained_bound), constrained_options,
            128));

    check_pair(preserved_class_group, preserved_units);
    assert(preserved_class_group.relation_count() == relation_count_before);
    const auto relations_after = preserved_class_group.relations();
    const auto invariant_generators_after =
            preserved_class_group.invariant_generator_matrix();
    const auto regulator_after = preserved_units.regulator();
    const auto torsion_order_after = preserved_units.torsion_order();
    assert(relations_after.has_value());
    assert(invariant_generators_after.has_value());
    assert(regulator_after.has_value());
    assert(torsion_order_after.has_value());
    assert(sflint::fmpz_mat_equal(*relations_before, *relations_after));
    assert(sflint::fmpz_mat_equal(*invariant_generators_before,
                                  *invariant_generators_after));
    assert(::arb_equal(regulator_before->raw(), regulator_after->raw()) != 0);
    assert(sflint::fmpz_equal(*torsion_order_before,
                              *torsion_order_after));

    return 0;
}

int test_nonmaximal_quadratic_pair_rejection_preserves_output() {
    silex::NumberField field = quadratic_field(-47);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(!equation.is_maximal());
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(bound), maximal));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext fresh_class_group;
    silex::OrderUnitGroup fresh_units;
    assert(!fresh_units.compute_with_class_group(
            fresh_class_group, equation, sflint::FmpzConstRef(bound),
            options, 128));
    assert(!fresh_class_group.has_factor_base());
    assert(!fresh_class_group.has_presentation());
    assert(!fresh_units.is_set());

    silex::ClassGroupContext preserved_class_group;
    silex::OrderUnitGroup preserved_units;
    assert(preserved_units.compute_with_class_group(
            preserved_class_group, maximal, sflint::FmpzConstRef(bound),
            options, 128));
    assert(preserved_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(preserved_units.certification_status() ==
           silex::CertificationMode::proven);
    assert(check_rank_zero_group(preserved_units, maximal, 2));

    const slong relation_count_before =
            preserved_class_group.relation_count();
    const auto relations_before = preserved_class_group.relations();
    const auto invariant_generators_before =
            preserved_class_group.invariant_generator_matrix();
    const auto class_order_before = preserved_class_group.order();
    const auto regulator_before = preserved_units.regulator();
    const auto torsion_order_before = preserved_units.torsion_order();
    silex::OrderElement torsion_before(maximal);
    assert(relations_before.has_value());
    assert(invariant_generators_before.has_value());
    assert(class_order_before.has_value());
    assert(regulator_before.has_value());
    assert(torsion_order_before.has_value());
    assert(preserved_units.torsion_generator(torsion_before));
    assert(sflint::fmpz_equal_si(*class_order_before, 5));

    assert(!preserved_units.compute_with_class_group(
            preserved_class_group, equation, sflint::FmpzConstRef(bound),
            options, 128));
    assert(silex::same_order_parent(preserved_class_group.parent(),
                                    &maximal));
    assert(silex::same_order_parent(preserved_units.parent(), &maximal));
    assert(preserved_class_group.relation_count() == relation_count_before);
    assert(preserved_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(preserved_units.certification_status() ==
           silex::CertificationMode::proven);

    const auto relations_after = preserved_class_group.relations();
    const auto invariant_generators_after =
            preserved_class_group.invariant_generator_matrix();
    const auto class_order_after = preserved_class_group.order();
    const auto regulator_after = preserved_units.regulator();
    const auto torsion_order_after = preserved_units.torsion_order();
    silex::OrderElement torsion_after(maximal);
    assert(relations_after.has_value());
    assert(invariant_generators_after.has_value());
    assert(class_order_after.has_value());
    assert(regulator_after.has_value());
    assert(torsion_order_after.has_value());
    assert(preserved_units.torsion_generator(torsion_after));
    assert(sflint::fmpz_mat_equal(*relations_before, *relations_after));
    assert(sflint::fmpz_mat_equal(*invariant_generators_before,
                                  *invariant_generators_after));
    assert(sflint::fmpz_equal(*class_order_before, *class_order_after));
    assert(::arb_equal(regulator_before->raw(), regulator_after->raw()) != 0);
    assert(sflint::fmpz_equal(*torsion_order_before,
                              *torsion_order_after));
    assert(torsion_before.equal(torsion_after));

    return 0;
}

int test_compute_with_class_group_quintic_proven() {
    silex::NumberField field = quintic_field(0, 0, -1, -1);
    silex::Order equation;
    equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 1000;
    options.zeta_bf_max_cutoff = 20000;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 128));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    // The native reference-style unit extraction can certify this pair without
    // preserving the older incidental surplus of five kernel rows.
    assert(class_group.relation_kernel_unit_count() >= units.free_rank());
    // A degree-five hR comes only from Belabas-Friedman, which assumes GRH,
    // so `proven` is reached by generation, proven units, and relation
    // saturation at every prime dividing h_cand; no analytic check is
    // recorded as the proof.
    assert(class_group.analytic_class_regulator_certification() !=
           silex::CertificationMode::proven);
    assert(class_group.factor_base_generation_checked_status() ==
           silex::ProofState::verified);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(class_group.zeta_bf_proof_status() ==
           silex::ProofState::not_checked);

    sflint::Fmpz class_order;
    assert(class_group.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_is_one(sflint::FmpzConstRef(class_order)));
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(units.free_rank() == 2);

    return 0;
}

int test_compute_with_class_group_quintic_small_bf_budget_trivial() {
    silex::NumberField field = quintic_field(0, 0, -1, -1);
    silex::Order equation;
    equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 5000;
    options.max_relations = 1000;
    options.zeta_bf_max_cutoff = 1;
    options.requested_certification = silex::CertificationMode::proven;

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    assert(units.compute_with_class_group(class_group, order,
                                          sflint::FmpzConstRef(bound),
                                          options, 128));
    assert(class_group.has_presentation());
    assert(class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(class_group.relation_saturation_status() ==
           silex::ProofState::verified);
    assert(class_group.unit_proof_status() == silex::ProofState::verified);
    assert(class_group.regulator_proof_status() ==
           silex::ProofState::verified);
    assert(class_group.zeta_bf_proof_status() ==
           silex::ProofState::not_checked);
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(units.free_rank() == 2);

    return 0;
}

int test_compute_with_class_group_failures_preserve_output() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);

    silex::OrderUnitGroup units;
    assert(units.compute(degree_one_order));
    assert(check_rank_zero_group(units, degree_one_order, 2));

    silex::ClassGroupContext class_group(degree_one_order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(class_group.build_factor_base(sflint::FmpzConstRef(bound)));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupComputeOptions options;
    options.max_candidates = 64;
    options.max_relations = 16;

    assert(!units.compute_with_class_group(class_group, order,
                                           sflint::FmpzConstRef(bound),
                                           options, 160));
    assert(check_rank_zero_group(units, degree_one_order, 2));
    assert(silex::same_order_parent(class_group.parent(), &degree_one_order));

    assert(!units.compute_with_class_group(class_group, degree_one_order,
                                           sflint::FmpzConstRef(bound),
                                           options, 0));
    assert(check_rank_zero_group(units, degree_one_order, 2));
    assert(silex::same_order_parent(class_group.parent(), &degree_one_order));

    return 0;
}

int test_relation_kernel_rescale_log_matrix_word_rounding() {
    // Source: reference implementation polarit2.c:RgM_rescale_to_int rescales t_REAL
    // matrices using the represented, word-rounded real precision before
    // integer LLL.  Relation-kernel unit selection uses this path for compact
    // log matrices, so keep it aligned with the class-group reference boundary.
    sflint::ArbMat matrix(2, 2);
    sflint::FmpzMat rows32(2, 2);
    sflint::FmpzMat rows64(2, 2);
    sflint::FmpzMat native_rows(2, 2);
    sflint::FmpzMat expected32(2, 2);
    sflint::FmpzMat expected64(2, 2);

    set_arb_mat_entry_si_2exp(matrix, 0, 0, 1, 0);
    add_arb_mat_entry_2exp(matrix, 0, 0, -100);
    set_arb_mat_entry_si_2exp(matrix, 0, 1, 1, -1);
    set_arb_mat_entry_si_2exp(matrix, 1, 0, 1, -66);
    set_arb_mat_entry_si_2exp(matrix, 1, 1, -3, -2);

    set_zmat_entry_str(expected32, 0, 0, "73786976294838206464");
    set_zmat_entry_str(expected32, 0, 1, "36893488147419103232");
    set_zmat_entry_si(expected32, 1, 0, 1);
    set_zmat_entry_str(expected32, 1, 1, "-55340232221128654848");

    set_zmat_entry_str(expected64, 0, 0,
                       "1267650600228229401496703205377");
    set_zmat_entry_str(expected64, 0, 1,
                       "633825300114114700748351602688");
    set_zmat_entry_str(expected64, 1, 0, "17179869184");
    set_zmat_entry_str(expected64, 1, 1,
                       "-950737950171172051122527404032");

    assert(silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows32), sflint::ArbMatConstRef(matrix), 65,
            32));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(rows32),
                                  sflint::FmpzMatConstRef(expected32)));
    assert(silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), 65,
            64));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(rows64),
                                  sflint::FmpzMatConstRef(expected64)));

    assert(!silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), 65,
            16));
    assert(!silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), 0,
            64));
    assert(!silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(rows64), sflint::ArbMatConstRef(matrix), -1,
            64));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(rows64),
                                  sflint::FmpzMatConstRef(expected64)));

    sflint::ArbMat empty_matrix(0, 0);
    sflint::FmpzMat empty_rows(0, 0);
    for (const slong word_bits : {slong{32}, slong{64}}) {
        const slong maximum_safe_precision =
                std::numeric_limits<slong>::max() - (word_bits - 1);
        assert(silex::detail::
                       relation_kernel_rescale_log_matrix_for_testing(
                               sflint::FmpzMatRef(empty_rows),
                               sflint::ArbMatConstRef(empty_matrix),
                               maximum_safe_precision, word_bits));
        sflint::FmpzMat overflow_rows(2, 2);
        sflint::fmpz_mat_set(sflint::FmpzMatRef(overflow_rows),
                             sflint::FmpzMatConstRef(expected64));
        assert(!silex::detail::
                        relation_kernel_rescale_log_matrix_for_testing(
                                sflint::FmpzMatRef(overflow_rows),
                                sflint::ArbMatConstRef(matrix),
                                maximum_safe_precision + 1, word_bits));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(overflow_rows),
                sflint::FmpzMatConstRef(expected64)));
    }

    const slong native_word_bits =
            silex::detail::relation_kernel_native_word_bits();
    assert(native_word_bits == 32 || native_word_bits == 64);
    assert(silex::detail::relation_kernel_rescale_log_matrix_for_testing(
            sflint::FmpzMatRef(native_rows), sflint::ArbMatConstRef(matrix),
            65, native_word_bits));
    assert(sflint::fmpz_mat_equal(
            sflint::FmpzMatConstRef(native_rows),
            native_word_bits == 32 ? sflint::FmpzMatConstRef(expected32)
                                   : sflint::FmpzMatConstRef(expected64)));
    return 0;
}

int test_reduced_regulator_reconstruction_correction() {
    constexpr slong precision = 256;

    auto reset_sentinels = [](sflint::Arb& regulator,
                              sflint::FmpzMat& coordinates) noexcept {
        sflint::arb_set_si(regulator, 42);
        coordinates = sflint::FmpzMat(1, 2);
        set_zmat_entry_si(coordinates, 0, 0, 47);
        set_zmat_entry_si(coordinates, 0, 1, -53);
    };
    auto assert_sentinels = [](const sflint::Arb& regulator,
                               const sflint::FmpzMat& coordinates) noexcept {
        assert(::arb_equal_si(regulator.raw(), 42) != 0);
        assert(sflint::fmpz_mat_nrows(coordinates) == 1);
        assert(sflint::fmpz_mat_ncols(coordinates) == 2);
        assert(mat_entry_is_si(coordinates, 0, 0, 47));
        assert(mat_entry_is_si(coordinates, 0, 1, -53));
    };
    auto call = [](sflint::Arb& out,
                   const sflint::ArbMat& coordinates,
                   const sflint::Arb& regulator_multiple,
                   const sflint::Arb& z,
                   sflint::FmpzMat* integer_coordinates) noexcept {
        return silex::detail::reduced_regulator_from_coordinates(
                sflint::ArbRef(out), coordinates,
                sflint::ArbConstRef(regulator_multiple),
                sflint::ArbConstRef(z), integer_coordinates, precision,
                nullptr);
    };

    sflint::ArbMat one_coordinate(1, 1);
    set_arb_mat_entry_si_2exp(one_coordinate, 0, 0, 1, 0);
    sflint::Arb z;
    sflint::arb_set_ui(z, 8);
    sflint::Arb regulator_multiple;
    sflint::Arb regulator_out;
    sflint::FmpzMat integer_out(1, 2);

    sflint::arb_set_ui(regulator_multiple, 1);
    ::arb_mul_2exp_si(regulator_multiple.raw(), regulator_multiple.raw(), -4);
    reset_sentinels(regulator_out, integer_out);
    assert(!call(regulator_out, one_coordinate, regulator_multiple, z,
                 &integer_out));
    assert_sentinels(regulator_out, integer_out);

    sflint::arb_set_ui(regulator_multiple, 1);
    ::arb_mul_2exp_si(regulator_multiple.raw(), regulator_multiple.raw(), -3);
    assert(call(regulator_out, one_coordinate, regulator_multiple, z,
                &integer_out));
    assert(::arb_equal(regulator_out.raw(), regulator_multiple.raw()) != 0);
    assert(sflint::fmpz_mat_nrows(integer_out) == 1);
    assert(sflint::fmpz_mat_ncols(integer_out) == 1);
    assert(mat_entry_is_si(integer_out, 0, 0, 1));

    sflint::arb_set_ui(regulator_multiple, 1);
    ::arb_mul_2exp_si(regulator_multiple.raw(), regulator_multiple.raw(), -2);
    sflint::arb_set_si(regulator_out, 42);
    assert(call(regulator_out, one_coordinate, regulator_multiple, z,
                nullptr));
    assert(::arb_equal(regulator_out.raw(), regulator_multiple.raw()) != 0);

    sflint::Arf lower;
    sflint::Arf upper;
    ::arf_set_ui(lower.raw(), 3);
    ::arf_mul_2exp_si(lower.raw(), lower.raw(), -5);
    ::arf_set_ui(upper.raw(), 5);
    ::arf_mul_2exp_si(upper.raw(), upper.raw(), -5);
    ::arb_set_interval_arf(regulator_multiple.raw(), lower.raw(), upper.raw(),
                           precision);
    reset_sentinels(regulator_out, integer_out);
    assert(!call(regulator_out, one_coordinate, regulator_multiple, z,
                 &integer_out));
    assert_sentinels(regulator_out, integer_out);

    sflint::ArbMat huge_error(1, 1);
    set_arb_mat_entry_si_2exp(huge_error, 0, 0, 2, 0);
    sflint::Fmpz huge_error_exponent;
    sflint::fmpz_set_si(sflint::FmpzRef(huge_error_exponent),
                        std::numeric_limits<slong>::max());
    sflint::fmpz_add_ui(sflint::FmpzRef(huge_error_exponent),
                        sflint::FmpzConstRef(huge_error_exponent), 1);
    assert(!sflint::fmpz_fits_si(
            sflint::FmpzConstRef(huge_error_exponent)));
    ::arb_add_error_2exp_fmpz(
            arb_mat_entry(huge_error.raw(), 0, 0),
            huge_error_exponent.raw());
    assert(::arb_is_finite(arb_mat_entry(huge_error.raw(), 0, 0)) != 0);
    assert(::arb_contains_si(arb_mat_entry(huge_error.raw(), 0, 0), 2) != 0);
    sflint::arb_one(regulator_multiple);
    sflint::arb_one(z);
    reset_sentinels(regulator_out, integer_out);
    assert(!call(regulator_out, huge_error, regulator_multiple, z,
                 &integer_out));
    assert_sentinels(regulator_out, integer_out);

    sflint::Fmpz q;
    sflint::Fmpz k;
    sflint::Fmpz x_denominator;
    sflint::fmpz_one(sflint::FmpzRef(q));
    sflint::fmpz_mul_2exp(sflint::FmpzRef(q), sflint::FmpzConstRef(q), 40);
    sflint::fmpz_add_ui(sflint::FmpzRef(k), sflint::FmpzConstRef(q), 1);
    sflint::fmpz_mul_2exp(sflint::FmpzRef(x_denominator),
                          sflint::FmpzConstRef(q), 1);
    sflint::fmpz_add_ui(sflint::FmpzRef(x_denominator),
                        sflint::FmpzConstRef(x_denominator), 1);

    sflint::Fmpz one;
    sflint::Fmpz two;
    sflint::fmpz_one(sflint::FmpzRef(one));
    sflint::fmpz_set_si(sflint::FmpzRef(two), 2);
    sflint::Fmpq x;
    sflint::Fmpq one_over_q;
    sflint::Fmpq one_over_k;
    sflint::fmpq_set_fmpz_frac(x, sflint::FmpzConstRef(two),
                               sflint::FmpzConstRef(x_denominator));
    sflint::fmpq_set_fmpz_frac(one_over_q, sflint::FmpzConstRef(one),
                               sflint::FmpzConstRef(q));
    sflint::fmpq_set_fmpz_frac(one_over_k, sflint::FmpzConstRef(one),
                               sflint::FmpzConstRef(k));
    sflint::Fmpq crossover;
    ::fmpq_add(crossover.raw(), one_over_q.raw(), one_over_k.raw());
    ::fmpq_div_2exp(crossover.raw(), crossover.raw(), 1);
    sflint::Fmpq midpoint;
    ::fmpq_add(midpoint.raw(), x.raw(), crossover.raw());
    ::fmpq_div_2exp(midpoint.raw(), midpoint.raw(), 1);
    assert(::fmpq_cmp(x.raw(), midpoint.raw()) < 0);
    assert(::fmpq_cmp(midpoint.raw(), crossover.raw()) < 0);

    sflint::ArbMat high_denominator(1, 1);
    arb_struct* high_denominator_entry =
            arb_mat_entry(high_denominator.raw(), 0, 0);
    ::arf_set_fmpq(arb_midref(high_denominator_entry), midpoint.raw(),
                   precision, ARF_RND_NEAR);
    ::mag_set_ui_2exp_si(arb_radref(high_denominator_entry), 1, -79);
    assert(sflint::arb_contains_fmpq(high_denominator_entry, one_over_q));
    assert(sflint::arb_contains_fmpq(high_denominator_entry, one_over_k));
    sflint::arb_set_fmpz(regulator_multiple, sflint::FmpzConstRef(k));
    sflint::arb_set_ui(z, 1);
    ::arb_mul_2exp_si(z.raw(), z.raw(), -1);

    sflint::Fmpq expected_reduced;
    sflint::fmpq_set_fmpz_frac(expected_reduced, sflint::FmpzConstRef(k),
                               sflint::FmpzConstRef(q));
    assert(call(regulator_out, high_denominator, regulator_multiple, z,
                &integer_out));
    assert(sflint::arb_contains_fmpq(sflint::ArbConstRef(regulator_out),
                                     expected_reduced));
    assert(!sflint::arb_contains_si(regulator_out, 1));
    assert(mat_entry_is_si(integer_out, 0, 0, 1));

    ::arb_neg(high_denominator_entry, high_denominator_entry);
    sflint::Fmpq negative_one_over_q;
    sflint::Fmpq negative_one_over_k;
    sflint::fmpq_neg(negative_one_over_q, one_over_q);
    sflint::fmpq_neg(negative_one_over_k, one_over_k);
    assert(sflint::arb_contains_fmpq(high_denominator_entry,
                                     negative_one_over_q));
    assert(sflint::arb_contains_fmpq(high_denominator_entry,
                                     negative_one_over_k));
    assert(call(regulator_out, high_denominator, regulator_multiple, z,
                &integer_out));
    assert(sflint::arb_contains_fmpq(sflint::ArbConstRef(regulator_out),
                                     expected_reduced));
    assert(!sflint::arb_contains_si(regulator_out, 1));
    assert(mat_entry_is_si(integer_out, 0, 0, -1));
    return 0;
}

int test_hnf_regulator_multiple_certification_correction() {
    constexpr slong precision = 256;

    auto call = [](sflint::Arb& out,
                   const sflint::ArbMat& basis) noexcept {
        return silex::detail::hnf_regulator_multiple_from_basis(
                sflint::ArbRef(out), basis, 2, precision);
    };
    auto reset_sentinel = [](sflint::Arb& out) noexcept {
        sflint::arb_set_si(out, 42);
    };
    auto assert_sentinel = [](const sflint::Arb& out) noexcept {
        assert(::arb_equal_si(out.raw(), 42) != 0);
    };
    auto assert_exact_inverse_product = [](const sflint::ArbMat& basis) {
        sflint::ArbMat inverse(2, 2);
        sflint::ArbMat product(2, 2);
        assert(::arb_mat_inv(inverse.raw(), basis.raw(), precision) != 0);
        ::arb_mat_mul(product.raw(), inverse.raw(), basis.raw(), precision);
        for (slong i = 0; i < 2; ++i) {
            for (slong j = 0; j < 2; ++j) {
                assert(::arb_equal_si(
                               arb_mat_entry(product.raw(), i, j),
                               i == j ? 1 : 0) != 0);
            }
        }
    };
    auto set_symmetric_basis = [](sflint::ArbMat& basis,
                                  slong mantissa,
                                  slong exponent) noexcept {
        set_arb_mat_entry_si_2exp(basis, 0, 0, 1, 0);
        set_arb_mat_entry_si_2exp(basis, 1, 0, 1, 0);
        set_arb_mat_entry_si_2exp(basis, 0, 1, -mantissa, exponent);
        set_arb_mat_entry_si_2exp(basis, 1, 1, mantissa, exponent);
    };

    sflint::Arb out;

    sflint::ArbMat passing(2, 2);
    set_symmetric_basis(passing, 1, 0);
    assert_exact_inverse_product(passing);
    reset_sentinel(out);
    assert(call(out, passing));
    assert(::arb_equal_si(out.raw(), 1) != 0);

    // This basis has an exact inverse but kR = 1/2 while its lower-right
    // minor is 1.  Inverse-residual validation alone cannot reject it.
    sflint::ArbMat mismatch(2, 2);
    set_arb_mat_entry_si_2exp(mismatch, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(mismatch, 1, 0, 1, 0);
    set_arb_mat_entry_si_2exp(mismatch, 0, 1, 0, 0);
    set_arb_mat_entry_si_2exp(mismatch, 1, 1, 1, 0);
    assert_exact_inverse_product(mismatch);
    reset_sentinel(out);
    assert(!call(out, mismatch));
    assert_sentinel(out);

    sflint::ArbMat consistent_small(2, 2);
    set_symmetric_basis(consistent_small, 1, -4);
    reset_sentinel(out);
    assert(!call(out, consistent_small));
    assert_sentinel(out);

    sflint::ArbMat threshold(2, 2);
    set_symmetric_basis(threshold, 1, -3);
    reset_sentinel(out);
    assert(call(out, threshold));
    sflint::Arb one_eighth;
    sflint::arb_one(one_eighth);
    ::arb_mul_2exp_si(one_eighth.raw(), one_eighth.raw(), -3);
    assert(::arb_equal(out.raw(), one_eighth.raw()) != 0);

    sflint::ArbMat threshold_overlap(2, 2);
    set_arb_mat_entry_si_2exp(threshold_overlap, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(threshold_overlap, 1, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(
            threshold_overlap, 1, 1, 3, -5, 5, -5);
    ::arb_neg(arb_mat_entry(threshold_overlap.raw(), 0, 1),
              arb_mat_entry(threshold_overlap.raw(), 1, 1));
    reset_sentinel(out);
    assert(!call(out, threshold_overlap));
    assert_sentinel(out);

    // A zero-containing difference ball remains certifiable through its
    // absolute upper bound; it must not be mistaken for threshold overlap.
    sflint::ArbMat rounded(2, 2);
    set_arb_mat_entry_si_2exp(rounded, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(rounded, 1, 0, 1, 0);
    set_arb_mat_entry_si_2exp(rounded, 0, 1, -1, 0);
    ::mag_set_ui_2exp_si(
            arb_radref(arb_mat_entry(rounded.raw(), 0, 1)), 1, -80);
    set_arb_mat_entry_si_2exp(rounded, 1, 1, 1, 0);
    reset_sentinel(out);
    assert(call(out, rounded));
    assert(::arb_contains_si(out.raw(), 1) != 0);
    assert(::arb_is_exact(out.raw()) == 0);
    sflint::Arb one;
    sflint::Arb difference;
    sflint::arb_one(one);
    sflint::arb_sub(difference, one, out, precision);
    assert(::arb_contains_zero(difference.raw()) != 0);
    assert(::arb_is_zero(difference.raw()) == 0);

    return 0;
}

int test_regulator_pivot_exponent_policy() {
    // The source regulator-multiple pivot compares integer exponent bins,
    // preserves the first row on a tie, and accepts only gexpo > -32.
    sflint::ArbMat threshold(1, 1);
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 0, 0);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{-1});
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 1, -32);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{-1});
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 3, -33);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{-1});
    set_arb_mat_entry_si_2exp(threshold, 0, 0, 1, -31);
    assert_regulator_pivot_rows(threshold, std::vector<slong>{0});

    sflint::ArbMat tie(2, 1);
    set_arb_mat_entry_si_2exp(tie, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(tie, 1, 0, 7, -2);
    assert_regulator_pivot_rows(tie, std::vector<slong>{0});
    set_arb_mat_entry_si_2exp(tie, 0, 0, 7, -2);
    set_arb_mat_entry_si_2exp(tie, 1, 0, 1, 0);
    assert_regulator_pivot_rows(tie, std::vector<slong>{0});
    set_arb_mat_entry_si_2exp(tie, 1, 0, 1, 1);
    assert_regulator_pivot_rows(tie, std::vector<slong>{1});

    sflint::ArbMat same_bin_ball(2, 1);
    set_arb_mat_entry_interval_si_2exp(
            same_bin_ball, 0, 0, 5, -2, 3, -1);
    set_arb_mat_entry_si_2exp(same_bin_ball, 1, 0, 7, -2);
    assert_regulator_pivot_rows(same_bin_ball, std::vector<slong>{0});

    sflint::ArbMat huge(2, 1);
    set_arb_mat_entry_huge_2exp(
            huge, 0, 0, "92233720368547758081234567890");
    set_arb_mat_entry_huge_2exp(
            huge, 1, 0, "92233720368547758081234567891");
    assert_regulator_pivot_rows(huge, std::vector<slong>{1});
    return 0;
}

int test_regulator_pivot_arb_ambiguity_and_invalidity() {
    std::vector<slong> pivots{47};
    sflint::ArbMat exact_zero(2, 1);
    set_arb_mat_entry_si_2exp(exact_zero, 0, 0, 0, 0);
    set_arb_mat_entry_si_2exp(exact_zero, 1, 0, 1, 0);
    assert_regulator_pivot_rows(exact_zero, std::vector<slong>{1});

    sflint::ArbMat later_equal(2, 1);
    set_arb_mat_entry_si_2exp(later_equal, 0, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(later_equal, 1, 0, 3, -2, 5, -2);
    assert_regulator_pivot_rows(later_equal, std::vector<slong>{0});

    sflint::ArbMat earlier_equal(2, 1);
    set_arb_mat_entry_interval_si_2exp(earlier_equal, 0, 0, 3, -2, 5, -2);
    set_arb_mat_entry_si_2exp(earlier_equal, 1, 0, 1, 0);
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(earlier_equal), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat later_outranks(2, 1);
    set_arb_mat_entry_si_2exp(later_outranks, 0, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(
            later_outranks, 1, 0, 3, -1, 5, -1);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(later_outranks), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat contains_zero(2, 1);
    set_arb_mat_entry_si_2exp(contains_zero, 0, 0, 1, 0);
    ::arb_zero(arb_mat_entry(contains_zero.raw(), 1, 0));
    ::arb_add_error_2exp_si(arb_mat_entry(contains_zero.raw(), 1, 0), -2);
    assert_regulator_pivot_rows(contains_zero, std::vector<slong>{0});
    ::arb_add_error_2exp_si(arb_mat_entry(contains_zero.raw(), 1, 0), 2);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(contains_zero), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat earlier_contains_zero(2, 1);
    ::arb_zero(arb_mat_entry(earlier_contains_zero.raw(), 0, 0));
    ::arb_add_error_2exp_si(
            arb_mat_entry(earlier_contains_zero.raw(), 0, 0), 0);
    set_arb_mat_entry_si_2exp(earlier_contains_zero, 1, 0, 1, 0);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(earlier_contains_zero), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat threshold_crossing(1, 1);
    set_arb_mat_entry_interval_si_2exp(
            threshold_crossing, 0, 0, 3, -33, 1, -31);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(threshold_crossing), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    sflint::ArbMat below_threshold(1, 1);
    set_arb_mat_entry_interval_si_2exp(
            below_threshold, 0, 0, 3, -35, 5, -34);
    assert_regulator_pivot_rows(below_threshold, std::vector<slong>{-1});

    sflint::ArbMat nonfinite_tail(1, 2);
    set_arb_mat_entry_si_2exp(nonfinite_tail, 0, 0, 1, 0);
    ::arb_pos_inf(arb_mat_entry(nonfinite_tail.raw(), 0, 1));
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(nonfinite_tail), 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(pivots.empty());

    // The first column is certified and pivoted before the second column's
    // threshold-crossing ball makes the result inconclusive.  No partial
    // pivot trace may escape.
    sflint::ArbMat late_ambiguity(2, 2);
    set_arb_mat_entry_si_2exp(late_ambiguity, 0, 0, 1, 0);
    set_arb_mat_entry_interval_si_2exp(
            late_ambiguity, 1, 1, 3, -33, 1, -31);
    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(late_ambiguity), 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(pivots.empty());

    pivots = {47};
    assert(silex::detail::regulator_pivot_rows_for_testing(
            pivots, sflint::ArbMatConstRef(exact_zero), 0) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(pivots.empty());
    return 0;
}

int test_regulator_pivot_clean_cols_and_mapping() {
    const silex::Signature signature(2, 0);
    std::vector<slong> selected{47};
    std::vector<slong> pivots{47};

    sflint::ArbMat discarded(1, 2);
    set_arb_mat_entry_si_2exp(discarded, 0, 0, 1, -3);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, discarded, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(selected.empty());
    assert((pivots == std::vector<slong>{0}));

    sflint::ArbMat retained(1, 2);
    set_arb_mat_entry_si_2exp(retained, 0, 0, 1, -2);
    set_arb_mat_entry_si_2exp(retained, 0, 1, 1, 0);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, retained, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert((selected == std::vector<slong>{0}));

    sflint::ArbMat ambiguous(1, 2);
    set_arb_mat_entry_interval_si_2exp(ambiguous, 0, 0, 3, -4, 5, -4);
    set_arb_mat_entry_si_2exp(ambiguous, 0, 1, 0, 0);
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, ambiguous, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(selected.empty());
    assert(pivots.empty());

    ::arb_zero(arb_mat_entry(ambiguous.raw(), 0, 0));
    ::arb_add_error_2exp_si(arb_mat_entry(ambiguous.raw(), 0, 0), -2);
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, ambiguous, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::precision_inconclusive);
    assert(selected.empty());
    assert(pivots.empty());

    // An earlier possible retention is resolved by a later certified entry.
    sflint::ArbMat resolved(1, 2);
    set_arb_mat_entry_interval_si_2exp(resolved, 0, 0, 3, -4, 5, -4);
    set_arb_mat_entry_si_2exp(resolved, 0, 1, 1, 0);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, resolved, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert((selected == std::vector<slong>{0}));

    // Certified retain cannot hide a later nonfinite entry.
    set_arb_mat_entry_si_2exp(resolved, 0, 0, 1, -2);
    ::arb_pos_inf(arb_mat_entry(resolved.raw(), 0, 1));
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, resolved, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(selected.empty());
    assert(pivots.empty());

    // Invalid in a later unit row dominates earlier cleaning uncertainty.
    sflint::ArbMat late_invalid(2, 2);
    set_arb_mat_entry_interval_si_2exp(
            late_invalid, 0, 0, 3, -4, 5, -4);
    ::arb_pos_inf(arb_mat_entry(late_invalid.raw(), 1, 1));
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, late_invalid, signature, 1, 2, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(selected.empty());
    assert(pivots.empty());

    // Two retained originals separated by a discarded row keep source order.
    const silex::Signature signature3(3, 0);
    sflint::ArbMat mapped(3, 3);
    set_arb_mat_entry_si_2exp(mapped, 0, 0, 1, -2);
    set_arb_mat_entry_si_2exp(mapped, 0, 1, 1, 0);
    set_arb_mat_entry_si_2exp(mapped, 1, 0, 1, -3);
    set_arb_mat_entry_si_2exp(mapped, 2, 1, 1, -2);
    set_arb_mat_entry_si_2exp(mapped, 2, 2, 1, 0);
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, mapped, signature3, 2, 3, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert((selected == std::vector<slong>{0, 2}));

    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, mapped, silex::Signature(2, 0), 2, 3, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);
    assert(selected.empty());
    assert(pivots.empty());
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, mapped,
            silex::Signature(std::numeric_limits<slong>::max(),
                             std::numeric_limits<slong>::max()),
            2, 3, 256) ==
           silex::detail::RegulatorPivotOutcome::invalid);

    // reference's rank-zero branch bypasses filter_unit_log_columns, so a finite ambiguous log
    // is irrelevant rather than precision-inconclusive.
    sflint::ArbMat rank_zero_ambiguous(1, 1);
    ::arb_zero(arb_mat_entry(rank_zero_ambiguous.raw(), 0, 0));
    ::arb_add_error_2exp_si(
            arb_mat_entry(rank_zero_ambiguous.raw(), 0, 0), -2);
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, rank_zero_ambiguous, silex::Signature(1, 0),
            0, 1, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(selected.empty());
    assert(pivots.empty());
    ::arb_pos_inf(arb_mat_entry(rank_zero_ambiguous.raw(), 0, 0));
    selected = {47};
    pivots = {47};
    assert(silex::detail::regulator_pivot_unit_indices_for_testing(
            selected, pivots, rank_zero_ambiguous, silex::Signature(1, 0),
            0, 1, 256) ==
           silex::detail::RegulatorPivotOutcome::success);
    assert(selected.empty());
    assert(pivots.empty());
    return 0;
}

int test_regulator_pivot_multi_column_elimination() {
    sflint::ArbMat matrix(2, 2);
    set_arb_mat_entry_si_2exp(matrix, 0, 0, 1, 0);
    set_arb_mat_entry_si_2exp(matrix, 0, 1, 1, 0);
    set_arb_mat_entry_si_2exp(matrix, 1, 0, 7, -2);
    set_arb_mat_entry_si_2exp(matrix, 1, 1, 3, 0);
    sflint::ArbMat input_snapshot(2, 2);
    ::arb_mat_set(input_snapshot.raw(), matrix.raw());
    assert_regulator_pivot_rows(matrix, std::vector<slong>{0, 1});
    assert(::arb_mat_equal(input_snapshot.raw(), matrix.raw()) != 0);
    assert_regulator_pivot_rows(matrix, std::vector<slong>{0, 1});
    assert(::arb_mat_equal(input_snapshot.raw(), matrix.raw()) != 0);
    return 0;
}

void assert_regulator_pivot_route_state_preserved(
        const silex::detail::RegulatorPivotFinishTestResult& result)
        noexcept {
    assert(result.goal_reached_after == result.goal_reached_before);
    assert(result.relation_need_after == result.relation_need_before);
    assert(result.completion_old_need_after ==
           result.completion_old_need_before);
    assert(result.completion_dependent_trials_after ==
           result.completion_dependent_trials_before);
    assert(result.completion_subfactor_base_trials_after ==
           result.completion_subfactor_base_trials_before);
    assert(result.analytic_extra_relation_requests_after ==
           result.analytic_extra_relation_requests_before);
    assert(result.finish_unit_log_rotation_after ==
           result.finish_unit_log_rotation_before);
    assert(result.squash_index_after == result.squash_index_before);
    assert(result.candidates_tried_after == result.candidates_tried_before);
    assert(result.accepted_relations_after ==
           result.accepted_relations_before);
    assert(result.finish_unit_log_rotation_active_after ==
           result.finish_unit_log_rotation_active_before);
    assert(result.finish_full_rank_relation_active_after ==
           result.finish_full_rank_relation_active_before);
    assert(result.relation_control_state_unchanged);
}

int test_regulator_pivot_finish_route_matrix() {
    using Action = silex::detail::RegulatorPivotFinishAction;
    using Outcome = silex::detail::RegulatorPivotOutcome;
    using Path = silex::detail::RegulatorPivotFinishPath;
    for (const Path path : {Path::independent_unit_preprobe,
                            Path::regulator_product}) {
        const auto full =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::success, 3, 3, true, 160, 0);
        assert(full.action == Action::proceed);
        assert(full.all_control_state_unchanged);

        const auto partial =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::success, 1, 3, true, 160, 0);
        assert(partial.action == Action::request_relations);
        assert(!partial.goal_reached_after);
        assert(partial.relation_need_after == 2);
        assert(partial.completion_old_need_after == 2);
        assert(partial.completion_dependent_trials_after == 0);
        assert(partial.completion_subfactor_base_trials_after ==
               partial.completion_subfactor_base_trials_before);
        assert(partial.candidates_tried_after ==
               partial.candidates_tried_before);
        assert(partial.accepted_relations_after ==
               partial.accepted_relations_before);
        assert(partial.analytic_extra_relation_requests_after ==
               partial.analytic_extra_relation_requests_before + 1);
        assert(partial.finish_unit_log_rotation_active_after);
        assert(partial.finish_unit_log_rotation_after ==
               partial.squash_index_before);
        assert(partial.squash_index_after == partial.squash_index_before + 1);
        assert(!partial.finish_full_rank_relation_active_after);
        assert(partial.analytic_finish_precision_after ==
               partial.analytic_finish_precision_before);
        assert(partial.analytic_precision_doublings_after ==
               partial.analytic_precision_doublings_before);
        assert(partial.analytic_precision_inconclusive_after ==
               partial.analytic_precision_inconclusive_before);
        assert(partial.factor_base_restart_requests_after ==
               partial.factor_base_restart_requests_before);
        assert(partial.factor_base_restart_pending_after ==
               partial.factor_base_restart_pending_before);
        assert(partial.factor_base_restart_allow_past_half_after ==
               partial.factor_base_restart_allow_past_half_before);

        const auto precision =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::precision_inconclusive, 0, 3,
                                false, 160, 0);
        assert(precision.action == Action::retry_precision);
        assert_regulator_pivot_route_state_preserved(precision);
        assert(precision.analytic_finish_precision_after == 320);
        assert(precision.analytic_precision_doublings_after ==
               precision.analytic_precision_doublings_before + 1);
        assert(precision.analytic_precision_inconclusive_after ==
               precision.analytic_precision_inconclusive_before + 1);
        assert(precision.factor_base_restart_requests_after ==
               precision.factor_base_restart_requests_before);
        assert(!precision.factor_base_restart_pending_after);
        assert(precision.factor_base_restart_allow_past_half_after ==
               precision.factor_base_restart_allow_past_half_before);

        const auto precision_cap =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::precision_inconclusive, 0, 3,
                                false, 40960, 0);
        assert(precision_cap.action == Action::failed);
        assert(precision_cap.all_control_state_unchanged);

        const auto restart =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::precision_inconclusive, 0, 3,
                                true, 160, 7);
        assert(restart.action == Action::restart_factor_base);
        assert_regulator_pivot_route_state_preserved(restart);
        assert(restart.analytic_finish_precision_after ==
               restart.analytic_finish_precision_before);
        assert(restart.analytic_precision_doublings_after ==
               restart.analytic_precision_doublings_before);
        assert(restart.analytic_precision_inconclusive_after ==
               restart.analytic_precision_inconclusive_before);
        assert(restart.factor_base_restart_requests_after ==
               restart.factor_base_restart_requests_before + 1);
        assert(restart.factor_base_restart_pending_after);
        assert(!restart.factor_base_restart_allow_past_half_after);

        const auto invalid =
                silex::detail::relation_search::
                        regulator_pivot_finish_for_testing(
                                path, Outcome::invalid, 0, 3, true, 160, 7);
        assert(invalid.action == Action::failed);
        assert(invalid.all_control_state_unchanged);
    }
    return 0;
}

int test_relation_kernel_units_rank_one() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(bound), options));
    assert(context.relation_kernel_unit_count() >= 6);

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    assert(group.set_relation_kernel_units(equation_order, context,
                                           embeddings, 160));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    assert(group.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
    auto owned_torsion_order = group.torsion_order();
    assert(owned_torsion_order.has_value());
    assert(sflint::fmpz_equal_si(*owned_torsion_order, 2));
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_positive(regulator));
    auto owned_regulator = group.regulator();
    assert(owned_regulator.has_value());
    assert(sflint::arb_is_positive(*owned_regulator));

    return 0;
}

int test_relation_kernel_units_bounded_rank_one() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz factor_bound;
    assert(set_fmpz_si(factor_bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() >= 6);

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    sflint::Fmpz denominator_bound;
    assert(set_fmpz_si(denominator_bound, 2));
    assert(group.set_relation_kernel_units_bounded(
            equation_order, context, embeddings,
            sflint::FmpzConstRef(denominator_bound), 80, 160));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    silex::OrderUnitGroup bound_one;
    assert(set_fmpz_si(denominator_bound, 1));
    assert(bound_one.set_relation_kernel_units_bounded(
            equation_order, context, embeddings,
            sflint::FmpzConstRef(denominator_bound), 80, 160));
    assert(bound_one.is_set());
    assert(bound_one.free_rank() == 1);
    assert(check_first_free_generator_is_order_unit(bound_one,
                                                    equation_order));

    return 0;
}

int test_relation_kernel_units_index_bounded_rank_one() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz factor_bound;
    assert(set_fmpz_si(factor_bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() >= 6);

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    assert(group.set_relation_kernel_units_index_bounded(
            equation_order, context, embeddings, 80, 160));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator_is_order_unit(group, equation_order));
    // The relation-kernel candidates come from an LLL-reduced kernel basis,
    // and the best one is the fundamental unit 1 + sqrt(2).  The regulator
    // index bound of that group is 1, so the index-one check proves it.
    silex::Element epsilon(field);
    assert(set_quadratic_coeffs(epsilon, 1, 1, 1, 1));
    assert(first_free_generator_is_pm_unit_power_pm1(group, epsilon));
    assert(group.certification_status() == silex::CertificationMode::proven);

    return 0;
}

int test_relation_kernel_units_index_bounded_quartic_power_root() {
    silex::NumberField field = quartic_field(-1, -1, 0, 2);
    silex::Order equation_order;
    silex::Order order(field);
    equation_order = silex::test::equation_order(field);
    assert(order.maximal_order(equation_order));

    // This fixture's relation kernel comes from a factor base of the former
    // generation-bound size.
    sflint::Fmpz factor_bound;
    assert(silex::test::former_generation_bound(factor_bound, order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_bound), 2);
    }

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() == 6);

    silex::EmbeddingContext embeddings(field);
    sflint::Arb analytic_hR;
    assert(silex::zeta_class_regulator_product(
            sflint::ArbRef(analytic_hR), order, 128));

    silex::OrderUnitGroup direct(order);
    sflint::Fmpz direct_bound;
    assert(direct.set_relation_kernel_units(order, context, embeddings, 128));
    assert(direct.class_regulator_index_bound(
            sflint::FmpzRef(direct_bound), context,
            sflint::ArbConstRef(analytic_hR), 128));
    // The six candidates come from an LLL-reduced kernel basis; the units
    // they generate have class-regulator index bound 3 (11 for the raw
    // Smith-form kernel rows).
    assert(sflint::fmpz_equal_si(direct_bound, 3));

    silex::OrderUnitGroup index_bounded(order);
    sflint::Fmpz bounded_value;
    assert(index_bounded.set_relation_kernel_units_index_bounded(
            order, context, embeddings, 128, 128));
    assert(index_bounded.class_regulator_index_bound(
            sflint::FmpzRef(bounded_value), context,
            sflint::ArbConstRef(analytic_hR), 128));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(bounded_value),
                              sflint::FmpzConstRef(direct_bound)));

    return 0;
}

namespace {

// One grh validation saturation pass (inexact index bound 2) on `context`.
bool grh_validation_saturation_pass(
        bool& changed,
        silex::ClassGroupContext& context,
        const silex::OrderUnitGroup& units,
        const silex::ClassGroupComputeOptions& options,
        silex::detail::RelationSaturationMemo* memo) {
    sflint::Fmpz index_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(index_bound), 2);
    return silex::detail::saturate_candidate_class_relations_with_units(
            changed, context, units, options,
            sflint::FmpzConstRef(index_bound), false, memo);
}

void assert_same_saturation_state(const silex::ClassGroupContext& left,
                                  const silex::ClassGroupContext& right) {
    assert(left.relation_count() == right.relation_count());
    sflint::Fmpz left_order;
    sflint::Fmpz right_order;
    assert(left.order(sflint::FmpzRef(left_order)));
    assert(right.order(sflint::FmpzRef(right_order)));
    assert(sflint::fmpz_equal(sflint::FmpzConstRef(left_order),
                              sflint::FmpzConstRef(right_order)));
    assert(left.relation_saturation_status() ==
           right.relation_saturation_status());
    assert(left.unit_proof_status() == right.unit_proof_status());
    assert(left.regulator_proof_status() == right.regulator_proof_status());
    assert(left.relation_saturation_record_count() ==
           right.relation_saturation_record_count());
    for (slong i = 0; i < left.relation_saturation_record_count(); ++i) {
        auto left_record = left.relation_saturation_record(i);
        auto right_record = right.relation_saturation_record(i);
        assert(left_record.has_value() && right_record.has_value());
        assert(sflint::fmpz_equal(sflint::FmpzConstRef(left_record->ell),
                                  sflint::FmpzConstRef(right_record->ell)));
        assert(left_record->status == right_record->status);
    }
}

}  // namespace

// The grh validation memo (detail::RelationSaturationMemo): a bounded
// saturation pass that found nothing is skipped while h_cand, the candidate
// regulator, the auxiliary-prime bound and the append limit are unchanged,
// and the skip leaves the state of a real pass that appends nothing.
// Q(sqrt(10)) has h = 2, so with a candidate h_cand = 2 no pass at ell = 2
// can append a relation.
int test_grh_validation_saturation_memo_replays_unchanged_pass() {
    silex::NumberField field = quadratic_field(10);
    silex::Order equation_order;
    silex::Order order(field);
    equation_order = silex::test::equation_order(field);
    assert(order.maximal_order(equation_order));

    sflint::Fmpz factor_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_bound), order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_bound), 2);
    }
    silex::ClassGroupCandidateOptions candidate_options;
    candidate_options.max_candidates = 5000;
    candidate_options.max_relations = 500;

    // `real` runs every pass; `memoed` and `fresh` run through the memo.
    // `fresh` stands for the class/unit loop's context after an iteration
    // that appended relations without changing h_cand: the append reset its
    // saturation records and status.
    silex::ClassGroupContext real;
    silex::ClassGroupContext memoed;
    silex::ClassGroupContext fresh;
    assert(real.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                  candidate_options));
    assert(memoed.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                    candidate_options));
    assert(fresh.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                   candidate_options));
    sflint::Fmpz h_cand;
    assert(real.order(sflint::FmpzRef(h_cand)));
    assert(sflint::fmpz_equal_si(h_cand, 2));
    assert_same_saturation_state(real, memoed);
    assert_same_saturation_state(real, fresh);

    silex::OrderUnitGroup units(order);
    assert(units.compute(order));

    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::grh;
    options.max_relations = 500;
    silex::detail::RelationSaturationMemo memo;

    // First pass: both run it, find nothing, and the memo records it.
    bool real_changed = true;
    bool memo_changed = true;
    assert(grh_validation_saturation_pass(real_changed, real, units, options,
                                          nullptr));
    assert(grh_validation_saturation_pass(memo_changed, memoed, units,
                                          options, &memo));
    assert(!real_changed && !memo_changed);
    assert(memo.valid);
    assert(memo.skips == 0);
    assert(sflint::fmpz_equal_si(memo.class_order, 2));
    assert(real.relation_saturation_status() ==
           silex::ProofState::unavailable);
    assert(real.relation_saturation_record_count() == 1);
    auto record = real.relation_saturation_record(0);
    assert(record.has_value());
    assert(sflint::fmpz_equal_si(record->ell, 2));
    assert(record->status == silex::ProofState::unavailable);
    assert(units.certification_status() == silex::CertificationMode::proven);
    assert(real.unit_proof_status() == silex::ProofState::verified);
    assert_same_saturation_state(real, memoed);

    // Same h_cand and units on a context with no saturation state: the memo
    // skips the pass and writes the state the real pass wrote.
    assert(fresh.relation_saturation_status() ==
           silex::ProofState::not_checked);
    assert(fresh.relation_saturation_record_count() == 0);
    assert(fresh.unit_proof_status() == silex::ProofState::not_checked);
    memo_changed = true;
    assert(grh_validation_saturation_pass(memo_changed, fresh, units, options,
                                          &memo));
    assert(!memo_changed);
    assert(memo.skips == 1);
    assert_same_saturation_state(real, fresh);

    // A repeated pass on the context that recorded the memo is skipped too.
    real_changed = true;
    memo_changed = true;
    assert(grh_validation_saturation_pass(real_changed, real, units, options,
                                          nullptr));
    assert(grh_validation_saturation_pass(memo_changed, memoed, units,
                                          options, &memo));
    assert(!real_changed && !memo_changed);
    assert(memo.skips == 2);
    assert_same_saturation_state(real, memoed);

    // A different auxiliary-prime bound asks a different question: the pass
    // runs and the memo is recorded again for the new bound.
    options.relation_saturation_aux_prime_bound = 37;
    assert(grh_validation_saturation_pass(memo_changed, memoed, units,
                                          options, &memo));
    assert(!memo_changed);
    assert(memo.skips == 2);
    assert(memo.valid);
    assert(sflint::fmpz_equal_si(memo.aux_bound, 37));
    assert_same_saturation_state(real, memoed);

    // The exact imaginary-quadratic index route never uses the memo.
    sflint::Fmpz exact_bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(exact_bound), 2);
    assert(silex::detail::saturate_candidate_class_relations_with_units(
            memo_changed, memoed, units, options,
            sflint::FmpzConstRef(exact_bound), true, &memo));
    assert(memo.skips == 2);

    // A memo whose h_cand differs from the current one is not used: the
    // pass runs and the memo is recorded again for the current h_cand.
    assert(memo.valid);
    sflint::fmpz_set_ui(sflint::FmpzRef(memo.class_order), 4);
    assert(grh_validation_saturation_pass(memo_changed, fresh, units, options,
                                          &memo));
    assert(!memo_changed);
    assert(memo.skips == 2);
    assert(memo.valid);
    assert(sflint::fmpz_equal_si(memo.class_order, 2));

    // Likewise for a different factor-base generator count.
    assert(memo.generator_count == fresh.generator_count());
    memo.generator_count = fresh.generator_count() + 1;
    assert(grh_validation_saturation_pass(memo_changed, fresh, units, options,
                                          &memo));
    assert(!memo_changed);
    assert(memo.skips == 2);
    assert(memo.valid);
    assert(memo.generator_count == fresh.generator_count());

    // Units generated by eps^2, eps = 3 + sqrt(10) the fundamental unit, have
    // twice the regulator: a regulator ratio of two is not within the
    // factor-two window, so the pass runs.
    silex::EmbeddingContext embeddings(field);
    silex::Element theta(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(theta.gen());
    assert(epsilon.add_si(theta, 3));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square_generator(field);
    assert(square_generator.set_element(epsilon2));
    silex::FactoredElement square_generators[] = {
            std::move(square_generator)};
    silex::OrderUnitGroup square_units(order);
    assert(square_units.set_units(
            order, silex::FactoredElementSpan(square_generators, 1),
            embeddings, 128));
    sflint::Arb regulator;
    sflint::Arb square_regulator;
    sflint::Arb twice_regulator;
    assert(units.regulator(sflint::ArbRef(regulator)));
    assert(square_units.regulator(sflint::ArbRef(square_regulator)));
    sflint::arb_mul_2exp_si(twice_regulator, regulator, 1);
    assert(sflint::arb_overlaps(square_regulator, twice_regulator));
    assert(memo.valid);
    assert(grh_validation_saturation_pass(memo_changed, fresh, square_units,
                                          options, &memo));
    assert(memo.skips == 2);

    return 0;
}

// A bounded pass that appends a relation invalidates the memo.  In the
// degree-one field with factor base {2} and the single relation (-4), h_cand
// is 2 and the pass at ell = 2 appends the relation (2) (as in
// test_relation_saturation_bounded_append_with_units_degree_one of
// t-class-group).
int test_grh_validation_saturation_memo_invalidated_by_change() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(-4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());
    sflint::Fmpz h_cand;
    assert(context.order(sflint::FmpzRef(h_cand)));
    assert(sflint::fmpz_equal_si(h_cand, 2));

    silex::OrderUnitGroup units(order);
    assert(units.compute(order));

    silex::ClassGroupComputeOptions options;
    options.requested_certification = silex::CertificationMode::grh;
    options.max_relations = 500;

    // A memo left valid by an earlier no-op pass with another h_cand.
    silex::detail::RelationSaturationMemo memo;
    memo.valid = true;
    sflint::fmpz_set_ui(sflint::FmpzRef(memo.class_order), 4);

    bool changed = false;
    assert(grh_validation_saturation_pass(changed, context, units, options,
                                          &memo));
    assert(changed);
    assert(context.relation_count() == 2);
    assert(memo.skips == 0);
    assert(!memo.valid);

    return 0;
}

int test_relation_kernel_units_index_bounded_cubic_rank_target() {
    silex::NumberField field = cubic_field(-2, -5);
    silex::Order equation_order;
    silex::Order order(field);
    equation_order = silex::test::equation_order(field);
    assert(order.maximal_order(equation_order));

    sflint::Fmpz factor_bound;
    assert(silex::factor_base_class_group_bound(
            sflint::FmpzRef(factor_bound), order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_bound), 2);
    }

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 5000;
    options.max_relations = 500;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(order, sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() == 1);

    sflint::Fmpz aux_start;
    sflint::Fmpz aux_max;
    assert(set_fmpz_si(aux_start, 2));
    assert(set_fmpz_si(aux_max, 31));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup refined(order);
    bool changed = false;
    bool stable = false;
    assert(refined.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, context, embeddings, 128, 128, 1,
            sflint::FmpzConstRef(aux_start),
            sflint::FmpzConstRef(aux_max), 2));

    sflint::Fmpz index_bound;
    // GP: K = bnfinit(x^3-2*x-5, 1) has K.reg = 2.35861..., h = 2,
    // signature (1, 1) and d_K = -643; the refined group has regulator
    // 9.43444... (index 4).  The order is maximal and the degree prime, so
    // the Pohst-Zassenhaus (6.22) bound R >= 1.05676... (above the Friedman
    // 1989 Table 6 bound 0.28) gives floor(9.43444... / 1.05676...) = 8.
    assert(refined.regulator_index_bound(sflint::FmpzRef(index_bound), 128));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(index_bound), 8));

    return 0;
}

int test_relation_kernel_units_index_bounded_saturated() {
    silex::NumberField field = quadratic_field(2);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(field);
    assert(equation_order.is_maximal());

    sflint::Fmpz factor_bound;
    assert(set_fmpz_si(factor_bound, 2));

    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 32;

    silex::ClassGroupContext context;
    assert(context.compute_candidate(equation_order,
                                     sflint::FmpzConstRef(factor_bound),
                                     options));
    assert(context.relation_kernel_unit_count() >= 6);

    sflint::Fmpz start;
    sflint::Fmpz max;
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    bool changed = true;
    bool stable = true;
    assert(group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, equation_order, context, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 3));
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &equation_order));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator_is_order_unit(group, equation_order));

    return 0;
}

int test_relation_kernel_units_rank_zero_delegates() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(4));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group;
    assert(group.set_relation_kernel_units(order, context, embeddings, 80));
    assert(check_rank_zero_group(group, order, 2));

    silex::OrderUnitGroup bounded;
    assert(bounded.set_relation_kernel_units_bounded(
            order, context, embeddings, sflint::FmpzConstRef(bound), 40, 80));
    assert(check_rank_zero_group(bounded, order, 2));

    silex::OrderUnitGroup index_bounded;
    assert(index_bounded.set_relation_kernel_units_index_bounded(
            order, context, embeddings, 40, 80));
    assert(check_rank_zero_group(index_bounded, order, 2));

    silex::OrderUnitGroup saturated;
    bool changed = true;
    bool stable = false;
    assert(saturated.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, context, embeddings, 40, 80, 1,
            sflint::FmpzConstRef(bound), sflint::FmpzConstRef(bound), 2));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(saturated, order, 2));

    return 0;
}

int test_relation_kernel_unit_failures_preserve_output() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);

    silex::OrderUnitGroup group;
    assert(group.compute(degree_one_order));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 2));
    silex::ClassGroupContext insufficient(order);
    assert(insufficient.build_factor_base(sflint::FmpzConstRef(bound)));
    silex::EmbeddingContext embeddings(field);
    assert(!group.set_relation_kernel_units(order, insufficient, embeddings,
                                            160));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    sflint::Fmpz denominator_bound;
    assert(set_fmpz_si(denominator_bound, 0));
    assert(!group.set_relation_kernel_units_bounded(
            order, insufficient, embeddings,
            sflint::FmpzConstRef(denominator_bound), 80, 160));
    assert(set_fmpz_si(denominator_bound, 2));
    assert(!group.set_relation_kernel_units_bounded(
            order, insufficient, embeddings,
            sflint::FmpzConstRef(denominator_bound), 160, 80));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    assert(!group.set_relation_kernel_units_index_bounded(
            order, insufficient, embeddings, 80, 160));
    assert(!group.set_relation_kernel_units_index_bounded(
            order, insufficient, embeddings, 160, 80));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    silex::NumberField other_field = quadratic_field(3);
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);
    silex::ClassGroupContext wrong_parent(other_order);
    silex::EmbeddingContext other_embeddings(other_field);
    assert(!group.set_relation_kernel_units(order, wrong_parent, embeddings,
                                            160));
    assert(!group.set_relation_kernel_units(order, insufficient,
                                            other_embeddings, 160));
    assert(!group.set_relation_kernel_units(order, insufficient, embeddings,
                                            0));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    sflint::Fmpz start;
    sflint::Fmpz max;
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));
    bool changed = true;
    bool stable = false;
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(changed);
    assert(!stable);
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, other_embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 160, 80, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 0,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(set_fmpz_si(start, 1));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(set_fmpz_si(start, 5));
    assert(set_fmpz_si(max, 4));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2));
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));
    assert(!group.set_relation_kernel_units_index_bounded_saturated(
            changed, stable, order, insufficient, embeddings, 80, 160, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 0));
    assert(changed);
    assert(!stable);
    assert(check_rank_zero_group(group, degree_one_order, 2));

    return 0;
}

int test_saturate_row_real_quadratic_square() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));

    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    // ell = 2 divides the torsion order 2 of Z[sqrt2], so rows carry the
    // torsion exponent in a final column.
    sflint::FmpzMat row(1, 2);
    sflint::Fmpz ell;
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    assert(set_fmpz_si(ell, 2));

    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    assert(saturated.saturate_row(changed, group,
                                  sflint::FmpzMatConstRef(row), 0,
                                  sflint::FmpzConstRef(ell), embeddings,
                                  128));
    assert(changed);
    assert(saturated.free_rank() == 1);
    assert(saturated.certification_status() ==
           silex::CertificationMode::unknown);
    assert(check_first_free_generator_square(saturated, epsilon2));

    return 0;
}

// One field-modulus cache shared by two saturation calls, as
// prove_index_bound shares it across its passes: the first call fills it,
// the second adds no entry, and both give the result of an uncached call.
int test_saturate_local_with_kernel_shares_field_modulus_cache() {
    using OrderUnitGroupAccess = silex::detail::OrderUnitGroupAccess;

    // x^4 - x - 1 (signature (2, 1), torsion +-1) has the units x and x - 1.
    sflint::FmpqPoly polynomial;
    const slong coefficients[] = {-1, -1, 0, 0, 1};
    for (slong i = 0; i < 5; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);

    silex::Element x(field);
    silex::Element one(field);
    silex::Element x_minus_one(field);
    assert(x.gen() && one.one() && x_minus_one.subtract(x, one));

    // x^2 in a compact form with several bases, and x - 1.
    silex::FactoredElement x_squared(field);
    silex::FactoredElement other(field);
    assert(x_squared.push(x, 3) && x_squared.push(x_minus_one, 2) &&
           x_squared.push(x_minus_one, -2) && x_squared.push(x, -1));
    assert(other.set_element(x_minus_one));
    silex::FactoredElement generators[] = {std::move(x_squared),
                                           std::move(other)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 128));
    assert(group.free_rank() == 2);

    // ell = 2 divides the torsion order, so the kernel has a torsion column.
    sflint::FmpzMat kernel(1, 3);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(kernel, 0, 0), 1);
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));

    silex::OrderUnitGroup uncached(order);
    bool uncached_changed = false;
    assert(OrderUnitGroupAccess::saturate_local_with_kernel(
            uncached, uncached_changed, group,
            sflint::FmpzMatConstRef(kernel), sflint::FmpzConstRef(ell),
            embeddings, 128, nullptr));
    assert(uncached_changed);
    assert(check_first_free_generator(uncached, x));

    silex::detail::CompactFieldModulusCache cache;
    silex::OrderUnitGroup first(order);
    bool first_changed = false;
    assert(OrderUnitGroupAccess::saturate_local_with_kernel(
            first, first_changed, group, sflint::FmpzMatConstRef(kernel),
            sflint::FmpzConstRef(ell), embeddings, 128, &cache));
    const std::size_t filled = cache.entries.size();
    assert(filled > 0 && cache.field.has_same_data(field));
    assert(first_changed && check_first_free_generator(first, x));

    silex::OrderUnitGroup second(order);
    bool second_changed = false;
    assert(OrderUnitGroupAccess::saturate_local_with_kernel(
            second, second_changed, group, sflint::FmpzMatConstRef(kernel),
            sflint::FmpzConstRef(ell), embeddings, 128, &cache));
    assert(cache.entries.size() == filled);
    assert(second_changed && check_first_free_generator(second, x));

    silex::FactoredElement expected(field);
    silex::FactoredElement actual(field);
    silex::Element expected_value(field);
    silex::Element actual_value(field);
    for (const silex::OrderUnitGroup* result : {&first, &second}) {
        assert(result->free_rank() == uncached.free_rank());
        for (slong i = 0; i < uncached.free_rank(); ++i) {
            assert(uncached.free_generator(expected, i) &&
                   result->free_generator(actual, i) &&
                   expected.evaluate(expected_value) &&
                   actual.evaluate(actual_value) &&
                   expected_value.equal(actual_value));
        }
    }
    return 0;
}

int test_saturate_row_no_root_and_divisible_copy() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::FmpzMat row(1, 2);  // free exponent, torsion exponent 0
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);

    silex::OrderUnitGroup out(order);
    bool changed = true;
    assert(out.saturate_row(changed, group, sflint::FmpzMatConstRef(row), 0,
                            sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!changed);
    assert(check_first_free_generator(out, epsilon));

    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 2);
    changed = true;
    assert(group.saturate_row(changed, group, sflint::FmpzMatConstRef(row), 0,
                              sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!changed);
    assert(check_first_free_generator(group, epsilon));

    return 0;
}

void set_arb_coordinate(sflint::ArbMat& coordinates,
                        slong row,
                        slong numerator,
                        slong denominator,
                        slong error_exponent) noexcept {
    sflint::Fmpq value;
    ::fmpq_set_si(value.raw(), numerator, static_cast<ulong>(denominator));
    arb_struct* entry = arb_mat_entry(coordinates.raw(), row, 0);
    ::arb_set_fmpq(entry, value.raw(), 512);
    ::arb_add_error_2exp_si(entry, error_exponent);
}

bool candidate_row_equals(const sflint::FmpzMat& candidate,
                          const slong* expected,
                          slong length) noexcept {
    for (slong i = 0; i < length; ++i) {
        if (!sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(candidate, 0, i), expected[i])) {
            return false;
        }
    }
    return true;
}

int test_dependent_relation_rational_candidate() {
    sflint::ArbMat coordinates(2, 1);
    sflint::FmpzMat candidate(1, 3);
    bool found = false;

    // Coordinates 1/1000003 and 5/7: the common denominator is
    // 7 * 1000003 = 7000021, recovered without scanning denominators.
    set_arb_coordinate(coordinates, 0, 1, 1000003, -200);
    set_arb_coordinate(coordinates, 1, 5, 7, -200);
    assert(silex::detail::dependent_relation_rational_candidate(
            found, candidate, coordinates, 10000000));
    assert(found);
    const slong expected[] = {-7, -5000015, 7000021};
    assert(candidate_row_equals(candidate, expected, 3));

    // The same coordinates with a denominator bound below 7000021.
    found = true;
    assert(silex::detail::dependent_relation_rational_candidate(
            found, candidate, coordinates, 7000020));
    assert(!found);

    // Integral coordinates give denominator one.
    set_arb_coordinate(coordinates, 0, 3, 1, -40);
    set_arb_coordinate(coordinates, 1, -2, 1, -40);
    assert(silex::detail::dependent_relation_rational_candidate(
            found, candidate, coordinates, 2));
    assert(found);
    const slong integral[] = {-3, 2, 1};
    assert(candidate_row_equals(candidate, integral, 3));

    // A ball with radius above one does not determine a candidate.
    set_arb_coordinate(coordinates, 0, 1, 3, 1);
    found = true;
    assert(silex::detail::dependent_relation_rational_candidate(
            found, candidate, coordinates, 1000));
    assert(!found);

    // Shape errors fail.
    sflint::FmpzMat wrong(1, 2);
    assert(!silex::detail::dependent_relation_rational_candidate(
            found, wrong, coordinates, 1000));
    assert(!silex::detail::dependent_relation_rational_candidate(
            found, candidate, coordinates, 0));
    return 0;
}

int test_dependent_relation_denominator_step() {
    sflint::ArbMat coordinates(2, 1);
    bool possible = false;
    slong step = 0;

    // Narrow balls around 3/7 and 5/14: only multiples of 14 can pass the
    // scan's integer test.
    set_arb_coordinate(coordinates, 0, 3, 7, -100);
    set_arb_coordinate(coordinates, 1, 5, 14, -100);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 128));
    assert(possible);
    assert(step == 14);

    // The same balls with bound 13: no m <= 13 is a multiple of 14.
    possible = true;
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 13, 128));
    assert(!possible);

    // A ball around 1/1009 with bound 1000 holds no rational with a
    // denominator <= 1000.
    set_arb_coordinate(coordinates, 0, 1, 1009, -100);
    set_arb_coordinate(coordinates, 1, 0, 1, -100);
    possible = true;
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 1000, 128));
    assert(!possible);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 1009, 128));
    assert(possible);
    assert(step == 1009);

    // Only narrow balls constrain the step: the ball 1/2 +- 1/4 holds 1/3,
    // 2/5, 3/5, ..., so it determines nothing.
    set_arb_coordinate(coordinates, 0, 1, 3, -100);
    set_arb_coordinate(coordinates, 1, 1, 2, -2);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 50, 128));
    assert(possible);
    assert(step == 3);

    // No narrow ball: every denominator stays a candidate.
    set_arb_coordinate(coordinates, 0, 1, 3, -4);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 50, 128));
    assert(possible);
    assert(step == 1);

    // The rounding margin of the scaled balls widens the balls by about
    // |x| 2^(8 - precision): at 12 bits the 2^-100 ball around 3/7 is no
    // longer narrow for bound 20.
    set_arb_coordinate(coordinates, 0, 3, 7, -100);
    set_arb_coordinate(coordinates, 1, 0, 1, -100);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 128));
    assert(possible);
    assert(step == 7);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 12));
    assert(possible);
    assert(step == 1);

    // Negative coordinates give the same steps as positive ones.
    set_arb_coordinate(coordinates, 0, -3, 7, -100);
    set_arb_coordinate(coordinates, 1, -5, 14, -100);
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 128));
    assert(possible);
    assert(step == 14);
    possible = true;
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 13, 128));
    assert(!possible);

    // A non-finite entry constrains nothing and is skipped.
    set_arb_coordinate(coordinates, 0, -3, 7, -100);
    ::arb_pos_inf(arb_mat_entry(coordinates.raw(), 1, 0));
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 128));
    assert(possible);
    assert(step == 7);
    ::arb_indeterminate(arb_mat_entry(coordinates.raw(), 1, 0));
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 128));
    assert(possible);
    assert(step == 7);
    ::arb_indeterminate(arb_mat_entry(coordinates.raw(), 0, 0));
    assert(silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 128));
    assert(possible);
    assert(step == 1);

    // The skipped denominators rest on FLINT's rounding in
    // arb_mat_scalar_mul_si: whenever the scan's integer test (that product
    // followed by arb_get_unique_fmpz) succeeds for every entry, m must be a
    // multiple of the step.  Pin this against the linked FLINT.
    {
        const slong precisions[] = {4, 8, 12, 20, 53, 128};
        const slong errors[] = {-6, -20, -40, -100, -200};
        const slong numerators[] = {1, -3, 5, 22, -77};
        const slong denominators[] = {1, 3, 7, 14, 30};
        const slong bound = 60;
        for (slong precision : precisions) {
            for (slong error : errors) {
                for (slong a : numerators) {
                    for (slong b : denominators) {
                        sflint::ArbMat balls(2, 1);
                        set_arb_coordinate(balls, 0, a, b, error);
                        set_arb_coordinate(balls, 1, a + 1, b + 2, error);
                        assert(silex::detail::
                                       dependent_relation_denominator_step(
                                               possible, step, balls, bound,
                                               precision));
                        for (slong m = 1; m <= bound; ++m) {
                            sflint::ArbMat scaled(2, 1);
                            ::arb_mat_set(scaled.raw(), balls.raw());
                            ::arb_mat_scalar_mul_si(scaled.raw(), scaled.raw(),
                                                    m, precision);
                            sflint::Fmpz integer;
                            bool passes = true;
                            for (slong i = 0; i < 2; ++i) {
                                passes = passes &&
                                         ::arb_get_unique_fmpz(
                                                 integer.raw(),
                                                 arb_mat_entry(scaled.raw(),
                                                               i, 0));
                            }
                            if (passes) {
                                assert(possible);
                                assert(m % step == 0);
                            }
                        }
                    }
                }
            }
        }
    }

    // Shape and argument errors fail.
    sflint::ArbMat wide(2, 2);
    assert(!silex::detail::dependent_relation_denominator_step(
            possible, step, wide, 20, 128));
    assert(!silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 0, 128));
    assert(!silex::detail::dependent_relation_denominator_step(
            possible, step, coordinates, 20, 0));
    return 0;
}

struct ScopeCounter {
    const char* label = nullptr;
    slong count = 0;
};

void scope_begin_callback(void* user,
                          silex::DiagnosticsModule,
                          silex::ProfileEvent event,
                          const char*,
                          const char* label) noexcept {
    ScopeCounter* counter = static_cast<ScopeCounter*>(user);
    if (event == silex::ProfileEvent::begin_scope && label != nullptr &&
        std::strcmp(label, counter->label) == 0) {
        ++counter->count;
    }
}

// The bounded dependent-relation search (set_relation_kernel_units_bounded)
// tries only the denominators that can pass its integer test, in the order
// of the full scan from min_denominator, and recovers the same relations.
int test_dependent_relation_min_denominator_scan() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement y(field);
    assert(y.set_element(epsilon));

    const auto group_of_power = [&](silex::OrderUnitGroup& group,
                                    slong exponent) {
        silex::FactoredElement generator(field);
        assert(generator.push(epsilon, exponent));
        silex::FactoredElement generators[] = {std::move(generator)};
        assert(silex::detail::order_unit_group_set_units_internal(
                group, order, silex::FactoredElementSpan(generators, 1),
                embeddings, 128, true));
    };
    const auto search = [&](bool& recovered, sflint::FmpzMat& rel,
                            const silex::OrderUnitGroup& group,
                            const silex::FactoredElement& value,
                            slong bound) {
        silex::FactoredElement root(field);
        sflint::Fmpz torsion_exp;
        sflint::Fmpz denominator_bound;
        assert(set_fmpz_si(denominator_bound, bound));
        return silex::detail::dependent_relation_bounded_min_denominator(
                recovered, root, rel, torsion_exp, group, value, embeddings,
                sflint::FmpzConstRef(denominator_bound), 2, 128, 128, false,
                false);
    };

    // <epsilon^1009>: epsilon has coordinate 1/1009, so the scan from 2
    // reaches its first candidate at 1009 and tries no other denominator.
    {
        silex::OrderUnitGroup group(order);
        group_of_power(group, 1009);
        silex::DiagnosticsContext diagnostics;
        silex::diagnostics_context_init(diagnostics);
        ScopeCounter counter;
        counter.label = "unit_group.dependent_relation_try_denominator";
        silex::diagnostics_set_profiling(
                diagnostics, true,
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group),
                scope_begin_callback, &counter);
        group.set_diagnostics(&diagnostics);
        bool recovered = false;
        sflint::FmpzMat rel(1, 2);
        assert(search(recovered, rel, group, y, 2000));
        assert(recovered);
        const slong expected[] = {-1, 1009};
        assert(candidate_row_equals(rel, expected, 2));
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
        assert(counter.count == 1);
#else
        assert(counter.count == 0);
#endif
        group.set_diagnostics(nullptr);

        // Bound 1008 admits no candidate, so nothing is recovered.
        recovered = true;
        assert(search(recovered, rel, group, y, 1008));
        assert(!recovered);
    }

    // <epsilon^1000000007> with bound 1000000007: one candidate, where the
    // full scan tried about 10^9 denominators.
    {
        silex::OrderUnitGroup group(order);
        group_of_power(group, 1000000007);
        bool recovered = false;
        sflint::FmpzMat rel(1, 2);
        assert(search(recovered, rel, group, y, 1000000007));
        assert(recovered);
        const slong expected[] = {-1, 1000000007};
        assert(candidate_row_equals(rel, expected, 2));
    }

    // y = epsilon^3 in <epsilon>: the coordinate 3 is integral, and the
    // scan starts at min_denominator 2, so it recovers (-6, 2) as before.
    {
        silex::OrderUnitGroup group(order);
        group_of_power(group, 1);
        silex::FactoredElement cube(field);
        assert(cube.push(epsilon, 3));
        bool recovered = false;
        sflint::FmpzMat rel(1, 2);
        assert(search(recovered, rel, group, cube, 50));
        assert(recovered);
        const slong expected[] = {-6, 2};
        assert(candidate_row_equals(rel, expected, 2));
    }
    return 0;
}

// rank_zero_torsion returns the first root-of-unity power of largest order
// lying in the order, with consistent element and coordinates.
int test_rank_zero_torsion_generators() {
    const auto check = [](const silex::Order& order, slong expected_order,
                          slong root_power) {
        const silex::NumberField* field = order.parent();
        assert(field != nullptr);
        sflint::Fmpz torsion_order;
        silex::OrderElement generator(order);
        assert(silex::detail::rank_zero_torsion(
                sflint::FmpzRef(torsion_order), generator, order));
        assert(sflint::fmpz_equal_si(torsion_order, expected_order));

        sflint::Fmpz field_order;
        silex::Element root(*field);
        assert(silex::roots_of_unity(sflint::FmpzRef(field_order), root,
                                     *field));
        silex::Element expected(*field);
        assert(silex::detail::compute_power(expected, root, root_power));
        silex::Element value(*field);
        assert(generator.get_element(value));
        assert(value.equal(expected));

        silex::OrderElement reference(order);
        assert(reference.set_element(expected));
        sflint::FmpzMat coordinates(1, order.degree());
        sflint::FmpzMat reference_coordinates(1, order.degree());
        assert(generator.get_coordinates(sflint::FmpzMatRef(coordinates)));
        assert(reference.get_coordinates(
                sflint::FmpzMatRef(reference_coordinates)));
        assert(::fmpz_mat_equal(coordinates.raw(),
                                reference_coordinates.raw()) != 0);
    };

    silex::NumberField rational = degree_one_field();
    silex::Order rational_order;
    rational_order = silex::test::equation_order(rational);
    check(rational_order, 2, 1);

    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Order sqrt2_order;
    sqrt2_order = silex::test::equation_order(sqrt2);
    check(sqrt2_order, 2, 1);

    silex::NumberField gaussian = quadratic_field(-1);
    silex::Order gaussian_order;
    gaussian_order = silex::test::equation_order(gaussian);
    check(gaussian_order, 4, 1);

    // Z[sqrt(-3)] misses the primitive sixth and third roots of unity, so
    // the first power in the order is the third, -1.
    silex::NumberField eisenstein = quadratic_field(-3);
    silex::Order equation;
    equation = silex::test::equation_order(eisenstein);
    assert(!equation.is_maximal());
    check(equation, 2, 3);
    silex::Order maximal(eisenstein);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    check(maximal, 6, 1);

    silex::NumberField zeta8 = quartic_field(0, 0, 0, 1);
    silex::Order zeta8_order;
    zeta8_order = silex::test::equation_order(zeta8);
    check(zeta8_order, 8, 1);
    return 0;
}

int test_add_dependent_unit_large_denominator() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));

    // The group <epsilon^N> has index N = 1000000007 in the unit group mod
    // torsion, so adjoining epsilon needs the relation
    // epsilon^N = (epsilon^N)^1, whose denominator is N.
    constexpr slong kIndex = 1000000007;
    silex::FactoredElement generator(field);
    assert(generator.push(epsilon, kIndex));
    silex::FactoredElement generators[] = {std::move(generator)};

    // epsilon^N is a unit by construction; the trusted install skips the
    // public per-generator check, which would expand it.
    silex::OrderUnitGroup group(order);
    assert(silex::detail::order_unit_group_set_units_internal(
            group, order, silex::FactoredElementSpan(generators, 1),
            embeddings, 128, true));
    sflint::Fmpz index_bound;
    assert(group.regulator_index_bound(sflint::FmpzRef(index_bound), 128));
    assert(::fmpz_cmp_si(index_bound.raw(), kIndex) >= 0);

    silex::FactoredElement candidate(field);
    assert(candidate.set_element(epsilon));
    silex::detail::RelationUnitExtractionState state;
    bool changed = false;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128));
    assert(changed);
    assert(group.free_rank() == 1);

    silex::OrderUnitGroup expected(order);
    silex::FactoredElement fundamental(field);
    assert(fundamental.set_element(epsilon));
    silex::FactoredElement fundamental_generators[] = {
            std::move(fundamental)};
    assert(expected.set_units(
            order, silex::FactoredElementSpan(fundamental_generators, 1),
            embeddings, 128));
    sflint::Arb regulator;
    sflint::Arb expected_regulator;
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(expected.regulator(sflint::ArbRef(expected_regulator)));
    assert(::arb_overlaps(regulator.raw(), expected_regulator.raw()) != 0);

    // Adjoining epsilon again finds the denominator-one relation and leaves
    // the group unchanged.
    changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128));
    assert(!changed);
    return 0;
}

int test_add_dependent_unit_overflowing_index_bound_fallback() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));

    // <epsilon^N> has index N above 2^62 in the unit group mod torsion,
    // and its a priori regulator index bound does not fit in slong.
    constexpr slong kIndex = (WORD(1) << 62) + (WORD(1) << 61) + 135;
    silex::FactoredElement generator(field);
    assert(generator.push(epsilon, kIndex));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(silex::detail::order_unit_group_set_units_internal(
            group, order, silex::FactoredElementSpan(generators, 1),
            embeddings, 128, true));
    sflint::Fmpz index_bound;
    assert(group.regulator_index_bound(sflint::FmpzRef(index_bound), 128));
    assert(::fmpz_fits_si(index_bound.raw()) == 0);

    silex::FactoredElement candidate(field);
    assert(candidate.set_element(epsilon));
    sflint::Arb starting_regulator;
    assert(group.regulator(sflint::ArbRef(starting_regulator)));

    // Without a fallback the search does not run and the group is kept.
    silex::detail::RelationUnitExtractionState state;
    bool changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128));
    assert(!changed);
    changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128,
                                             nullptr));
    assert(!changed);

    // A fallback below the true index only misses the relation.
    silex::detail::DependentUnitIndexFallback small;
    ::fmpz_set_si(small.bound.raw(), kIndex / 2);
    small.computed = true;
    small.usable = true;
    changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128,
                                             &small));
    assert(!changed);
    sflint::Arb regulator;
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(::arb_equal(regulator.raw(), starting_regulator.raw()) != 0);

    // An unusable fallback is ignored.
    silex::detail::DependentUnitIndexFallback unusable;
    ::fmpz_set_si(unusable.bound.raw(), kIndex);
    unusable.computed = true;
    changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128,
                                             &unusable));
    assert(!changed);

    // A fallback computed from a class group context with neither verified
    // generation nor a GRH generation record is not usable, so the group is
    // kept.
    silex::ClassGroupContext bare_class_group;
    silex::detail::AnalyticClassRegulatorCache bare_cache;
    silex::detail::DependentUnitIndexFallback computed_fallback;
    computed_fallback.class_group = &bare_class_group;
    computed_fallback.analytic_cache = &bare_cache;
    computed_fallback.order = &order;
    computed_fallback.precision = 128;
    changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128,
                                             &computed_fallback));
    assert(!changed);
    assert(computed_fallback.computed);
    assert(!computed_fallback.usable);
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(::arb_equal(regulator.raw(), starting_regulator.raw()) != 0);

    // A word-size fallback at least the index replaces the overflowing
    // bound, and adjoining epsilon recovers the full unit group.
    silex::detail::DependentUnitIndexFallback fallback;
    ::fmpz_set_si(fallback.bound.raw(), kIndex);
    fallback.computed = true;
    fallback.usable = true;
    changed = false;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128,
                                             &fallback));
    assert(changed);
    assert(group.free_rank() == 1);

    silex::OrderUnitGroup expected(order);
    silex::FactoredElement fundamental(field);
    assert(fundamental.set_element(epsilon));
    silex::FactoredElement fundamental_generators[] = {
            std::move(fundamental)};
    assert(expected.set_units(
            order, silex::FactoredElementSpan(fundamental_generators, 1),
            embeddings, 128));
    sflint::Arb expected_regulator;
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(expected.regulator(sflint::ArbRef(expected_regulator)));
    assert(::arb_overlaps(regulator.raw(), expected_regulator.raw()) != 0);
    return 0;
}

int test_dependent_relation_wrong_simplest_rational_rejected() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.push(epsilon, 5));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(silex::detail::order_unit_group_set_units_internal(
            group, order, silex::FactoredElementSpan(generators, 1),
            embeddings, 128, true));

    // y = epsilon^3 has coordinate 3/5 against the generator epsilon^5.
    silex::FactoredElement y(field);
    assert(y.push(epsilon, 3));
    constexpr slong kPrecision = 128;
    sflint::ArbVec generator_log(2);
    sflint::ArbVec y_log(2);
    assert(generators[0].logarithmic_embedding(
            sflint::ArbVecRef(generator_log), embeddings,
            silex::LogEmbeddingMode::product, kPrecision));
    assert(y.logarithmic_embedding(sflint::ArbVecRef(y_log), embeddings,
                                   silex::LogEmbeddingMode::product,
                                   kPrecision));
    sflint::ArbMat inverse(1, 1);
    ::arb_inv(arb_mat_entry(inverse.raw(), 0, 0), generator_log.data(),
              kPrecision);

    // Widen the inverse so the coordinate ball is about 3/5 +- 7/100.  The
    // ball contains the true value 3/5, but its simplest rational is 2/3.
    sflint::Arb widen;
    ::arb_set_si(widen.raw(), 7);
    ::arb_div_si(widen.raw(), widen.raw(), 100, kPrecision);
    sflint::Arb magnitude;
    ::arb_abs(magnitude.raw(), y_log.data());
    ::arb_div(widen.raw(), widen.raw(), magnitude.raw(), kPrecision);
    sflint::ArbMat wide_inverse(1, 1);
    ::arb_set(arb_mat_entry(wide_inverse.raw(), 0, 0),
              arb_mat_entry(inverse.raw(), 0, 0));
    ::arb_add_error(arb_mat_entry(wide_inverse.raw(), 0, 0), widen.raw());

    sflint::ArbMat coordinates(1, 1);
    ::arb_mul(arb_mat_entry(coordinates.raw(), 0, 0), y_log.data(),
              arb_mat_entry(wide_inverse.raw(), 0, 0), kPrecision);
    sflint::Fmpq true_value;
    ::fmpq_set_si(true_value.raw(), 3, 5);
    sflint::Arb true_ball;
    ::arb_set_fmpq(true_ball.raw(), true_value.raw(), kPrecision);
    assert(::arb_contains(arb_mat_entry(coordinates.raw(), 0, 0),
                          true_ball.raw()) != 0);
    sflint::FmpzMat guess(1, 2);
    bool found = false;
    assert(silex::detail::dependent_relation_rational_candidate(
            found, guess, coordinates, 5));
    assert(found);
    const slong wrong[] = {-2, 3};
    assert(candidate_row_equals(guess, wrong, 2));

    // The exact check rejects the wrong candidate y^3 = (epsilon^5)^2.
    silex::FactoredElement root(field);
    sflint::FmpzMat relation(1, 2);
    sflint::Fmpz torsion_exp;
    bool recovered = true;
    sflint::Fmpz bound;
    ::fmpz_set_si(bound.raw(), 5);
    assert(silex::detail::dependent_relation_bounded_with_inverse(
            recovered, root, relation, torsion_exp, group, y, embeddings,
            wide_inverse, sflint::FmpzConstRef(bound), kPrecision, true,
            false));
    assert(!recovered);

    // The precise inverse gives the true relation y^5 = (epsilon^5)^3.
    assert(silex::detail::dependent_relation_bounded_with_inverse(
            recovered, root, relation, torsion_exp, group, y, embeddings,
            inverse, sflint::FmpzConstRef(bound), kPrecision, true, false));
    assert(recovered);
    const slong expected[] = {-3, 5};
    assert(candidate_row_equals(relation, expected, 2));
    return 0;
}

// On the require_y_root path the exact torsion check decides a candidate
// relation, as reference Unit/UnitGrpCtx.jl `_add_dependent_unit!` does with
// `_check_relation_mod_torsion`; no m-th root is extracted.  For
// y = -epsilon against <epsilon^3>, y^3 * (epsilon^3)^-1 = -1, so the
// relation (-1, 3) holds with torsion exponent 1 and the root stays y.
int test_dependent_relation_y_root_uses_torsion_check() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.push(epsilon, 3));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(silex::detail::order_unit_group_set_units_internal(
            group, order, silex::FactoredElementSpan(generators, 1),
            embeddings, 128, true));

    silex::Element minus_epsilon(field);
    assert(minus_epsilon.negate(epsilon));
    silex::FactoredElement y(field);
    assert(y.set_element(minus_epsilon));

    constexpr slong kPrecision = 128;
    sflint::ArbVec generator_log(2);
    assert(generators[0].logarithmic_embedding(
            sflint::ArbVecRef(generator_log), embeddings,
            silex::LogEmbeddingMode::product, kPrecision));
    sflint::ArbMat inverse(1, 1);
    ::arb_inv(arb_mat_entry(inverse.raw(), 0, 0), generator_log.data(),
              kPrecision);
    sflint::Fmpz bound;
    ::fmpz_set_si(bound.raw(), 3);

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    ScopeCounter counter;
    counter.label = "unit_group.dependent_relation_power_root";
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            scope_begin_callback, &counter);
    group.set_diagnostics(&diagnostics);

    const slong expected[] = {-1, 3};
    for (const bool require_torsion_exponent : {true, false}) {
        silex::FactoredElement root(field);
        sflint::FmpzMat relation(1, 2);
        sflint::Fmpz torsion_exp;
        ::fmpz_set_si(torsion_exp.raw(), 7);
        bool recovered = false;
        assert(silex::detail::dependent_relation_bounded_with_inverse(
                recovered, root, relation, torsion_exp, group, y, embeddings,
                inverse, sflint::FmpzConstRef(bound), kPrecision, true,
                require_torsion_exponent));
        assert(recovered);
        assert(candidate_row_equals(relation, expected, 2));
        silex::Element root_value(field);
        assert(root.evaluate(root_value));
        assert(root_value.equal(minus_epsilon));
        if (require_torsion_exponent) {
            assert(::fmpz_equal_si(torsion_exp.raw(), 1));
        } else {
            assert(::fmpz_equal_si(torsion_exp.raw(), 0) ||
                   ::fmpz_equal_si(torsion_exp.raw(), 1));
        }
    }
    assert(counter.count == 0);
    group.set_diagnostics(nullptr);
    return 0;
}

int test_saturate_row_rank_zero_and_failures() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);

    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));
    assert(check_rank_zero_group(rank_zero, degree_one_order, 2));

    // Rank zero with ell = 2 | w = 2: a row is only the torsion exponent.
    sflint::FmpzMat zero_width(1, 0);
    sflint::FmpzMat torsion_only(1, 1);
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    assert(copied.saturate_row(changed, rank_zero,
                               sflint::FmpzMatConstRef(torsion_only), 0,
                               sflint::FmpzConstRef(ell),
                               degree_one_embeddings, 80));
    assert(!changed);
    assert(check_rank_zero_group(copied, degree_one_order, 2));
    assert(!copied.saturate_row(changed, rank_zero,
                                sflint::FmpzMatConstRef(zero_width), 0,
                                sflint::FmpzConstRef(ell),
                                degree_one_embeddings, 80));
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::FmpzMat row(1, 2);
    sflint::FmpzMat free_width(1, 1);  // missing the torsion column
    sflint::FmpzMat bad_width(1, 3);
    sflint::Fmpz not_prime;
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(free_width, 0, 0), 1);
    assert(set_fmpz_si(not_prime, 4));
    // A rejected row leaves the output flag and group untouched.
    changed = true;
    assert(!copied.saturate_row(changed, group,
                                sflint::FmpzMatConstRef(free_width), 0,
                                sflint::FmpzConstRef(ell), embeddings, 128));
    assert(changed);
    assert(check_rank_zero_group(copied, degree_one_order, 2));
    assert(!copied.saturate_row(changed, group, sflint::FmpzMatConstRef(row),
                                1, sflint::FmpzConstRef(ell), embeddings,
                                128));
    assert(!copied.saturate_row(changed, group,
                                sflint::FmpzMatConstRef(bad_width), 0,
                                sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!copied.saturate_row(changed, group, sflint::FmpzMatConstRef(row),
                                0, sflint::FmpzConstRef(not_prime),
                                embeddings, 128));
    assert(!copied.saturate_row(changed, group, sflint::FmpzMatConstRef(row),
                                0, sflint::FmpzConstRef(ell), embeddings, 0));
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    return 0;
}

int test_residue_dlog_kernel_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 5));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    // 5 is inert, so the residue field is F_25.  Both eps = 1 + sqrt2 (its
    // norm to F_5 is N(eps) = -1 = 4, a square) and -1 are squares in F_25,
    // so with ell = 2 | w = 2 the kernel is all of F_2^2: the free column and
    // the torsion column.
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());
    sflint::FmpzMat kernel(0, 0);
    assert(group.residue_dlog_kernel(kernel, prime_span,
                                     sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(kernel) == 2);
    assert(sflint::fmpz_mat_nrows(kernel) == 2);
    assert(mat_entry_is_si(kernel, 0, 0, 1));
    assert(mat_entry_is_si(kernel, 0, 1, 0));
    assert(mat_entry_is_si(kernel, 1, 0, 0));
    assert(mat_entry_is_si(kernel, 1, 1, 1));

    // ell = 3 does not divide w = 2: no torsion column.  7 = 1 mod 3 splits
    // in Q(sqrt2) (3^2 = 2 mod 7).
    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(ell, 3));
    silex::PrimeIdealList split;
    assert(silex::decompose_prime(split, order, sflint::FmpzConstRef(p)));
    assert(split.size() == 2);
    assert(group.residue_dlog_kernel(
            kernel, silex::PrimeIdealSpan(split.at(0), split.size()),
            sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(kernel) == 1);

    // Failure with ell = 2 | w = 2 leaves the output untouched: the prime
    // above 2 has no ell-th-power residue character, both alone and after a
    // prime above 7 whose torsion-extended column succeeds.
    assert(set_fmpz_si(p, 2));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList above_ell;
    assert(silex::decompose_prime(above_ell, order, sflint::FmpzConstRef(p)));
    assert(above_ell.size() == 1);
    silex::PrimeIdeal mixed[2];
    assert(mixed[0].define(order) && mixed[0].set(*split.at(0)));
    assert(mixed[1].define(order) && mixed[1].set(*above_ell.at(0)));
    sflint::FmpzMat sentinel(1, 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(sentinel, 0, 0), 7);
    assert(!group.residue_dlog_kernel(
            sentinel, silex::PrimeIdealSpan(above_ell.at(0), above_ell.size()),
            sflint::FmpzConstRef(ell)));
    assert(!group.residue_dlog_kernel(sentinel,
                                      silex::PrimeIdealSpan(mixed, 2),
                                      sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(sentinel) == 1);
    assert(sflint::fmpz_mat_ncols(sentinel) == 1);
    assert(mat_entry_is_si(sentinel, 0, 0, 7));

    return 0;
}

int test_residue_dlog_proof_kernel_torsion() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::OrderUnitGroup group;
    assert(group.compute(order));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 3));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);

    sflint::FmpzMat kernel(1, 1);
    assert(group.residue_dlog_proof_kernel(
            kernel, silex::PrimeIdealSpan(primes.at(0), primes.size()),
            sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 1);

    assert(set_fmpz_si(ell, 3));
    assert(group.residue_dlog_proof_kernel(
            kernel, silex::PrimeIdealSpan(primes.at(0), primes.size()),
            sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 0);

    return 0;
}

// Independent oracle for the residue_dlog_proof_kernel regression below:
// evaluate an ambient element's alpha-polynomial at alpha = s mod q and take
// its Legendre symbol, written additively (0 for a square, 1 otherwise).
slong legendre_bit_at_alpha(const silex::Element& element,
                            slong s,
                            slong q) noexcept {
    sflint::FmpqPoly polynomial;
    assert(element.get_fmpq_poly(sflint::FmpqPolyRef(polynomial)));
    slong value = 0;
    for (slong i = fmpq_poly_degree(polynomial.raw()); i >= 0; --i) {
        sflint::Fmpq coefficient;
        fmpq_poly_get_coeff_fmpq(coefficient.raw(), polynomial.raw(), i);
        const slong numerator = static_cast<slong>(
                fmpz_fdiv_ui(fmpq_numref(coefficient.raw()),
                             static_cast<ulong>(q)));
        const slong denominator = static_cast<slong>(
                fmpz_fdiv_ui(fmpq_denref(coefficient.raw()),
                             static_cast<ulong>(q)));
        assert(denominator != 0);
        const slong inverse = static_cast<slong>(n_invmod(
                static_cast<ulong>(denominator), static_cast<ulong>(q)));
        value = (value * s + numerator * inverse) % q;
    }
    assert(value != 0);
    return n_jacobi(static_cast<ulong>(value), static_cast<ulong>(q)) == 1
            ? 0
            : 1;
}

// Returns whether the mod-2 row space of `kernel` is exactly the set of
// vectors v in F_2^2 with v . column = 0 for every oracle column.
bool kernel_matches_oracle(const sflint::FmpzMat& kernel,
                           const std::vector<std::pair<slong, slong>>& columns)
        noexcept {
    if (sflint::fmpz_mat_ncols(kernel) != 2) {
        return false;
    }
    bool in_span[4] = {false, false, false, false};
    const slong rows = sflint::fmpz_mat_nrows(kernel);
    if (rows > 2) {
        return false;
    }
    for (slong mask = 0; mask < (slong(1) << rows); ++mask) {
        slong a = 0;
        slong b = 0;
        for (slong r = 0; r < rows; ++r) {
            if ((mask >> r) & 1) {
                a += static_cast<slong>(fmpz_fdiv_ui(
                        fmpz_mat_entry(kernel.raw(), r, 0), 2));
                b += static_cast<slong>(fmpz_fdiv_ui(
                        fmpz_mat_entry(kernel.raw(), r, 1), 2));
            }
        }
        in_span[2 * (a % 2) + (b % 2)] = true;
    }
    for (slong a = 0; a < 2; ++a) {
        for (slong b = 0; b < 2; ++b) {
            bool expected = true;
            for (const auto& column : columns) {
                expected = expected &&
                           (a * column.first + b * column.second) % 2 == 0;
            }
            if (in_span[2 * a + b] != expected) {
                return false;
            }
        }
    }
    return true;
}

int test_residue_dlog_proof_kernel_nonpower_basis_order() {
    // x^2 - 12 stays on the generic backend; its maximal order [1, alpha/2]
    // is not a power basis.  epsilon = 2 + alpha/2 is the fundamental unit,
    // and the subgroup <-1, epsilon^2> is not 2-saturated.  The residue
    // images of epsilon^2 are squares at every degree-one prime, so the
    // proof kernel at the primes over 23 and 37 must contain (1, 0).
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -12);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    assert(field.backend_kind() == silex::NumberFieldBackendKind::generic);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(!order.is_equation_order());

    silex::Element epsilon_squared(field);
    {
        sflint::FmpqPoly unit;
        sflint::fmpq_poly_set_coeff_si(unit, 0, 7);
        sflint::fmpq_poly_set_coeff_si(unit, 1, 2);
        assert(epsilon_squared.set_fmpq_poly(sflint::FmpqPolyConstRef(unit)));
    }
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon_squared));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));
    assert(group.free_rank() == 1);

    silex::FactoredElement free_generator(field);
    silex::Element free_value(field);
    silex::Element minus_one(field);
    silex::Element alpha(field);
    assert(alpha.gen());
    assert(group.free_generator(free_generator, 0));
    assert(free_generator.evaluate(free_value));
    assert(minus_one.set_si(-1));

    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 2));
    std::vector<silex::PrimeIdeal> all_primes;
    std::vector<std::pair<slong, slong>> all_columns;
    for (const slong q : {23, 37}) {
        sflint::Fmpz p;
        assert(set_fmpz_si(p, q));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
        assert(primes.size() == 2);
        std::vector<silex::PrimeIdeal> over_q;
        std::vector<std::pair<slong, slong>> q_columns;
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(prime != nullptr);
            assert(prime->residue_degree() == 1);

            // Find alpha mod P from ideal membership alone: alpha - s in P.
            silex::Ideal ideal(order);
            assert(prime->get_ideal(ideal));
            slong s = -1;
            for (slong t = 0; t < q; ++t) {
                if ((t * t - 12) % q != 0) {
                    continue;
                }
                silex::Element difference(field);
                silex::OrderElement difference_order(order);
                assert(difference.add_si(alpha, -t));
                assert(difference_order.set_element(difference));
                if (ideal.contains(difference_order)) {
                    assert(s == -1);
                    s = t;
                }
            }
            assert(s >= 0);

            const std::pair<slong, slong> column{
                    legendre_bit_at_alpha(free_value, s, q),
                    legendre_bit_at_alpha(minus_one, s, q)};
            // epsilon^2 is a square in K, so its image is a square.
            assert(column.first == 0);

            // The residue-field column and the direct degree-one evaluation
            // column agree with the oracle.
            sflint::FmpzMat via_residue_field(2, 1);
            sflint::FmpzMat via_direct(2, 1);
            assert(silex::detail::residue_dlog_proof_matrix(
                    via_residue_field, group, *prime,
                    sflint::FmpzConstRef(ell)));
            assert(silex::detail::saturation_proof_prime_column_direct_degree_one(
                    via_direct, group, *prime, sflint::FmpzConstRef(ell)));
            assert(sflint::fmpz_mat_nrows(via_residue_field) == 2);
            assert(sflint::fmpz_mat_nrows(via_direct) == 2);
            assert(static_cast<slong>(fmpz_fdiv_ui(
                           fmpz_mat_entry(via_residue_field.raw(), 0, 0), 2)) ==
                   column.first);
            assert(static_cast<slong>(fmpz_fdiv_ui(
                           fmpz_mat_entry(via_residue_field.raw(), 1, 0), 2)) ==
                   column.second);
            assert(fmpz_mat_equal(via_residue_field.raw(),
                                  via_direct.raw()) != 0);

            std::vector<silex::PrimeIdeal> single;
            single.emplace_back(order);
            assert(single.back().set(*prime));
            sflint::FmpzMat kernel(0, 0);
            assert(group.residue_dlog_proof_kernel(
                    kernel, silex::PrimeIdealSpan(single.data(), 1),
                    sflint::FmpzConstRef(ell)));
            assert(kernel_matches_oracle(kernel, {column}));

            over_q.emplace_back(order);
            assert(over_q.back().set(*prime));
            all_primes.emplace_back(order);
            assert(all_primes.back().set(*prime));
            q_columns.push_back(column);
            all_columns.push_back(column);
        }

        sflint::FmpzMat kernel(0, 0);
        assert(group.residue_dlog_proof_kernel(
                kernel, silex::PrimeIdealSpan(over_q.data(), over_q.size()),
                sflint::FmpzConstRef(ell)));
        assert(kernel_matches_oracle(kernel, q_columns));
    }

    sflint::FmpzMat kernel(0, 0);
    assert(group.residue_dlog_proof_kernel(
            kernel, silex::PrimeIdealSpan(all_primes.data(), all_primes.size()),
            sflint::FmpzConstRef(ell)));
    assert(kernel_matches_oracle(kernel, all_columns));
    // (1, 0), the epsilon^2 exponent vector, lies in the kernel.
    bool contains_unit_vector = false;
    for (slong r = 0; r < sflint::fmpz_mat_nrows(kernel); ++r) {
        contains_unit_vector =
                contains_unit_vector ||
                (fmpz_fdiv_ui(fmpz_mat_entry(kernel.raw(), r, 0), 2) == 1 &&
                 fmpz_fdiv_ui(fmpz_mat_entry(kernel.raw(), r, 1), 2) == 0);
    }
    assert(contains_unit_vector);
    return 0;
}

int test_rank_one_degree_one_root_image_matches_prime_ideal() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));

    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::FactoredElement stored_generator(field);
    assert(group.free_generator(stored_generator, 0));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(ell, 2));

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p), 1));
    assert(primes.size() > 0);

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);

        bool via_prime = false;
        bool via_root = false;
        bool via_nmod_root = false;
        sflint::Fmpz root;
        assert(degree_one_residue_root(root, *prime));
        assert(silex::detail::saturation_proof_prime_known_rank_one_free_generator_nonzero(
                via_prime, group, stored_generator, *prime,
                sflint::FmpzConstRef(ell)));
        assert(silex::detail::rank_one_factored_image_nonzero_at_degree_one_root(
                via_root, stored_generator, sflint::FmpzConstRef(p),
                sflint::FmpzConstRef(root), sflint::FmpzConstRef(ell)));
        assert(silex::detail::rank_one_factored_image_nonzero_at_degree_one_root_nmod(
                via_nmod_root, stored_generator, sflint::fmpz_get_ui(
                                                        sflint::FmpzConstRef(p)),
                sflint::fmpz_get_ui(sflint::FmpzConstRef(root)),
                sflint::FmpzConstRef(ell)));
        assert(via_root == via_prime);
        assert(via_nmod_root == via_prime);
        ++checked;
    }
    assert(checked == primes.size());
    return 0;
}

int test_direct_degree_one_proof_column_matches_residue_field() {
    silex::OrderUnitGroup group = real_quadratic_unit_group();
    const silex::Order* order = group.parent();
    assert(order != nullptr);
    if (order == nullptr) {
        return 1;
    }

    sflint::Fmpz p;
    sflint::Fmpz ell;
    const bool prime_set = set_fmpz_si(p, 7) && set_fmpz_si(ell, 2);
    assert(prime_set);
    if (!prime_set) {
        return 1;
    }

    silex::PrimeIdealList primes;
    const bool decomposed = silex::decompose_prime(
            primes, *order, sflint::FmpzConstRef(p), 1);
    assert(decomposed);
    assert(primes.size() > 0);
    if (!decomposed || primes.size() == 0) {
        return 1;
    }

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        if (prime == nullptr) {
            return 1;
        }

        sflint::FmpzMat generic_column(0, 0);
        sflint::FmpzMat direct_column(0, 0);
        const bool columns_match =
                silex::detail::saturation_proof_prime_column(
                        generic_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                silex::detail::saturation_proof_prime_column_direct_degree_one(
                        direct_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(generic_column),
                        sflint::FmpzMatConstRef(direct_column));
        assert(columns_match);
        if (!columns_match) {
            return 1;
        }
        ++checked;
    }
    assert(checked == primes.size());
    return checked == primes.size() ? 0 : 1;
}

bool generic_residue_dlog_column_for_test(sflint::FmpzMat& out,
                                          const silex::OrderUnitGroup& group,
                                          const silex::PrimeIdeal& prime,
                                          sflint::FmpzConstRef ell) noexcept {
    const silex::Order* order = group.parent();
    const silex::NumberField* field =
            order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    if (!group.is_set() || field == nullptr) {
        return false;
    }

    silex::ResidueField residue_field(prime);
    if (!residue_field.is_defined()) {
        return false;
    }

    sflint::FmpzMat candidate(rank, 1);
    silex::ResidueFieldElement image(residue_field);
    silex::ResidueFieldQuotientLog quotient_log(residue_field);
    if (!quotient_log.is_defined() || !quotient_log.set_ell(ell)) {
        return false;
    }
    for (slong i = 0; i < rank; ++i) {
        silex::FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !quotient_log.apply(sflint::fmpz_mat_entry(candidate, i, 0),
                                image)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

// Per-entry oracle for residue_dlog_proof_matrix: every entry, including the
// torsion entry when ell divides the torsion order, takes its own
// ResidueFieldQuotientLog::apply.
bool per_entry_residue_dlog_proof_matrix_for_test(
        sflint::FmpzMat& out,
        const silex::OrderUnitGroup& group,
        const silex::PrimeIdeal& prime,
        sflint::FmpzConstRef ell) noexcept {
    const silex::Order* order = group.parent();
    const silex::NumberField* field =
            order == nullptr ? nullptr : order->parent();
    const slong rank = group.free_rank();
    slong proof_rank = 0;
    if (!group.is_set() || field == nullptr ||
        !silex::detail::dlog_proof_rank(proof_rank, group, ell)) {
        return false;
    }

    silex::ResidueField residue_field(prime);
    silex::ResidueFieldElement image(residue_field);
    silex::ResidueFieldQuotientLog quotient_log(residue_field);
    if (!residue_field.is_defined() || !quotient_log.is_defined() ||
        !quotient_log.set_ell(ell)) {
        return false;
    }

    sflint::FmpzMat candidate(proof_rank, 1);
    for (slong i = 0; i < rank; ++i) {
        silex::FactoredElement generator(*field);
        if (!group.free_generator(generator, i) ||
            !image.set_factored_element(generator) ||
            !quotient_log.apply(sflint::fmpz_mat_entry(candidate, i, 0),
                                image)) {
            return false;
        }
    }
    if (proof_rank > rank) {
        silex::OrderElement torsion(*order);
        if (!group.torsion_generator(torsion) ||
            !image.set_order_element(torsion) ||
            !quotient_log.apply(sflint::fmpz_mat_entry(candidate, rank, 0),
                                image)) {
            return false;
        }
    }

    out = std::move(candidate);
    return true;
}

// residue_dlog_proof_matrix shares one walk over the quotient-generator
// powers across a column; every entry must equal the per-entry quotient log.
// The cyclic cubic x^3 - 3x + 1 with units theta and theta - 1 gives two
// free entries, plus the torsion entry -1 at ell = 2.
int test_residue_dlog_proof_matrix_matches_per_entry_quotient_log() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order order;
    order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    silex::Element theta_minus_one(cubic);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));
    silex::FactoredElement generators[] = {
            silex::FactoredElement(cubic), silex::FactoredElement(cubic)};
    assert(generators[0].set_element(theta));
    assert(generators[1].push(theta_minus_one, 3));
    assert(generators[1].push(theta, -2));
    silex::EmbeddingContext embeddings(cubic);
    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 64));
    assert(group.free_rank() == 2);

    const slong ells[] = {2, 3, 5, 7, 13, 31, 101, 1009};
    slong checked = 0;
    slong torsion_columns = 0;
    for (const slong ell_si : ells) {
        sflint::Fmpz ell;
        assert(set_fmpz_si(ell, ell_si));
        slong proof_rank = 0;
        assert(silex::detail::dlog_proof_rank(
                proof_rank, group, sflint::FmpzConstRef(ell)));
        slong primes_for_ell = 0;
        for (ulong q = static_cast<ulong>(ell_si) + 1; primes_for_ell < 6;
             q += static_cast<ulong>(ell_si)) {
            if (q == 3 || n_is_prime(q) == 0) {
                continue;
            }
            sflint::Fmpz q_fmpz;
            assert(set_fmpz_si(q_fmpz, static_cast<slong>(q)));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(
                    primes, order, sflint::FmpzConstRef(q_fmpz), 1));
            for (slong i = 0; i < primes.size(); ++i) {
                const silex::PrimeIdeal* prime = primes.at(i);
                assert(prime != nullptr);
                if (prime == nullptr) {
                    return 1;
                }
                sflint::FmpzMat shared(0, 0);
                sflint::FmpzMat expected(0, 0);
                const bool matches =
                        silex::detail::residue_dlog_proof_matrix(
                                shared, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        per_entry_residue_dlog_proof_matrix_for_test(
                                expected, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        sflint::fmpz_mat_nrows(shared) == proof_rank &&
                        sflint::fmpz_mat_equal(
                                sflint::FmpzMatConstRef(shared),
                                sflint::FmpzMatConstRef(expected));
                assert(matches);
                if (!matches) {
                    return 1;
                }
                // The public column path hands its own residue field to the
                // shared-walk code; it must agree with the oracle too.
                sflint::FmpzMat via_column(0, 0);
                const bool column_matches =
                        silex::detail::saturation_proof_prime_column(
                                via_column, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        sflint::fmpz_mat_equal(
                                sflint::FmpzMatConstRef(via_column),
                                sflint::FmpzMatConstRef(expected));
                assert(column_matches);
                if (!column_matches) {
                    return 1;
                }
                ++checked;
                torsion_columns += proof_rank > group.free_rank() ? 1 : 0;
            }
            if (primes.size() > 0) {
                ++primes_for_ell;
            }
        }
    }
    // Inert primes (residue degree 3 in this cyclic cubic): the shared walk
    // declines and every entry takes ResidueFieldQuotientLog::apply.
    slong higher_degree = 0;
    const slong inert_primes[] = {2, 5, 13};
    for (const slong q : inert_primes) {
        sflint::Fmpz q_fmpz;
        assert(set_fmpz_si(q_fmpz, q));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(
                primes, order, sflint::FmpzConstRef(q_fmpz), 3));
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(prime != nullptr);
            if (prime == nullptr || prime->residue_degree() < 2) {
                continue;
            }
            for (const slong ell_si : {2, 3, 7, 31, 61}) {
                if ((q * q * q - 1) % ell_si != 0) {
                    continue;
                }
                sflint::Fmpz ell;
                assert(set_fmpz_si(ell, ell_si));
                sflint::FmpzMat shared(0, 0);
                sflint::FmpzMat expected(0, 0);
                const bool matches =
                        silex::detail::residue_dlog_proof_matrix(
                                shared, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        per_entry_residue_dlog_proof_matrix_for_test(
                                expected, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        sflint::fmpz_mat_equal(
                                sflint::FmpzMatConstRef(shared),
                                sflint::FmpzMatConstRef(expected));
                assert(matches);
                if (!matches) {
                    return 1;
                }
                ++higher_degree;
            }
        }
    }
    assert(higher_degree >= 4);
    if (higher_degree < 4) {
        return 1;
    }
    assert(checked >= 48);
    assert(torsion_columns > 0);
    return checked >= 48 && torsion_columns > 0 ? 0 : 1;
}

// From ell of 7 bits on a degree-one column takes its quotient logs by
// baby-step giant-step with one baby-step table; every entry must still
// equal the per-entry quotient log.  ell = 61 (6 bits) stays on the shared
// walk, while 67 (7 bits), 262139 and 262147 (18 and 19 bits, either side of
// the reference table cutoff) and 1048583 (21 bits) take baby-step
// giant-step.  The cyclic cubic
// x^3 - 3x + 1 splits completely at each q used, giving three degree-one
// primes per q.
int test_residue_dlog_proof_matrix_baby_giant_matches_per_entry_quotient_log() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order order;
    order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    silex::Element theta_minus_one(cubic);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));
    silex::FactoredElement generators[] = {
            silex::FactoredElement(cubic), silex::FactoredElement(cubic)};
    assert(generators[0].set_element(theta));
    assert(generators[1].push(theta_minus_one, 3));
    assert(generators[1].push(theta, -2));
    silex::EmbeddingContext embeddings(cubic);
    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 64));

    const slong ells[] = {61, 67, 262139, 262147, 1048583};
    slong checked = 0;
    for (const slong ell_si : ells) {
        sflint::Fmpz ell;
        assert(set_fmpz_si(ell, ell_si));
        slong primes_for_ell = 0;
        for (ulong q = static_cast<ulong>(ell_si) + 1; primes_for_ell < 2;
             q += static_cast<ulong>(ell_si)) {
            if (n_is_prime(q) == 0) {
                continue;
            }
            sflint::Fmpz q_fmpz;
            assert(set_fmpz_si(q_fmpz, static_cast<slong>(q)));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(
                    primes, order, sflint::FmpzConstRef(q_fmpz), 1));
            for (slong i = 0; i < primes.size(); ++i) {
                const silex::PrimeIdeal* prime = primes.at(i);
                assert(prime != nullptr);
                if (prime == nullptr) {
                    return 1;
                }
                sflint::FmpzMat shared(0, 0);
                sflint::FmpzMat expected(0, 0);
                const bool matches =
                        silex::detail::residue_dlog_proof_matrix(
                                shared, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        per_entry_residue_dlog_proof_matrix_for_test(
                                expected, group, *prime,
                                sflint::FmpzConstRef(ell)) &&
                        sflint::fmpz_mat_nrows(shared) == 2 &&
                        sflint::fmpz_mat_equal(
                                sflint::FmpzMatConstRef(shared),
                                sflint::FmpzMatConstRef(expected));
                assert(matches);
                if (!matches) {
                    return 1;
                }
                if (ell_si == 262147 && q == 1048589 && i == 0) {
                    // External spot value.  The prime has root r = 688049 of
                    // x^3 - 3x + 1 and quotient generator zeta = 16; the
                    // column's theta entry is the discrete log of
                    // r^((q-1)/ell) = r^4 to base zeta.  gp:
                    //   znlog(Mod(688049,1048589)^4, Mod(16,1048589), 262147)
                    // returns 94028 (and Mod(16,1048589)^262147 == 1).
                    silex::ResidueField rf(*prime);
                    silex::ResidueFieldQuotientLog ql(rf);
                    silex::ResidueFieldElement zeta_el(rf);
                    silex::ResidueFieldElement img(rf);
                    sflint::Fmpz zs;
                    sflint::Fmpz rs;
                    assert(ql.set_ell(sflint::FmpzConstRef(ell)) &&
                           ql.quotient_generator(zeta_el) &&
                           zeta_el.degree_one_scalar(sflint::FmpzRef(zs)) &&
                           img.set_element(theta) &&
                           img.degree_one_scalar(sflint::FmpzRef(rs)));
                    assert(sflint::fmpz_get_si(sflint::FmpzConstRef(zs)) ==
                           16);
                    assert(sflint::fmpz_get_si(sflint::FmpzConstRef(rs)) ==
                           688049);
                    const slong spot = sflint::fmpz_get_si(
                            sflint::fmpz_mat_entry(
                                    sflint::FmpzMatConstRef(shared), 0, 0));
                    assert(spot == 94028);
                    if (spot != 94028) {
                        return 1;
                    }
                }
                ++checked;
            }
            if (primes.size() > 0) {
                ++primes_for_ell;
            }
        }
    }
    assert(checked == 30);
    return checked == 30 ? 0 : 1;
}

// On a degree-one prime with word-size p the shared walk (or baby-step
// giant-step) resolves every nonzero entry, so the only entries a column
// leaves for ResidueFieldQuotientLog::apply are zero residues, which apply
// rejects; a unit group's column cannot produce one, so this drives
// residue_quotient_log_column with chosen images.  A column mixing resolved
// entries with one zero entry must fail like the per-entry quotient log,
// wherever the zero entry sits and on both routes (ell = 7 at q = 71 walks,
// ell = 67 at q = 269 and ell = 262147 at q = 1048589 take baby-step
// giant-step); the same column without it must equal the per-entry values.
int test_residue_quotient_log_column_mixes_resolved_and_apply_entries() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order order;
    order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    silex::Element theta_minus_one(cubic);
    silex::Element theta_plus_two(cubic);
    assert(theta.gen());
    assert(theta_minus_one.add_si(theta, -1));
    assert(theta_plus_two.add_si(theta, 2));
    const silex::Element* values[] = {&theta, &theta_minus_one,
                                      &theta_plus_two};
    const slong value_count = 3;

    const slong cases[][2] = {{7, 71}, {67, 269}, {262147, 1048589}};
    slong checked = 0;
    for (const auto& entry : cases) {
        sflint::Fmpz ell;
        sflint::Fmpz q;
        assert(set_fmpz_si(ell, entry[0]) && set_fmpz_si(q, entry[1]));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, order,
                                      sflint::FmpzConstRef(q), 1));
        assert(primes.size() == 3);
        for (slong i = 0; i < primes.size(); ++i) {
            const silex::PrimeIdeal* prime = primes.at(i);
            assert(prime != nullptr);
            if (prime == nullptr) {
                return 1;
            }
            silex::ResidueField residue_field(*prime);
            silex::ResidueFieldQuotientLog quotient_log(residue_field);
            assert(residue_field.is_defined() && residue_field.degree() == 1);
            assert(quotient_log.set_ell(sflint::FmpzConstRef(ell)));

            // Without a zero entry: every entry equals its own apply.
            std::vector<silex::ResidueFieldElement> images;
            for (slong k = 0; k < value_count; ++k) {
                images.emplace_back(residue_field);
                assert(images.back().set_element(*values[k]));
                assert(!images.back().is_zero());
            }
            sflint::FmpzMat column(value_count, 1);
            sflint::Fmpz expected;
            bool matches = silex::detail::residue_quotient_log_column(
                    column, images, residue_field, quotient_log,
                    sflint::FmpzConstRef(ell), nullptr);
            for (slong k = 0; matches && k < value_count; ++k) {
                matches = quotient_log.apply(
                                  sflint::FmpzRef(expected),
                                  images[static_cast<std::size_t>(k)]) &&
                          ::fmpz_equal(expected.raw(),
                                       ::fmpz_mat_entry(column.raw(), k,
                                                        0)) != 0;
            }
            assert(matches);
            if (!matches) {
                return 1;
            }

            // A zero entry at each position: apply rejects it, so the
            // column fails.
            for (slong zero_at = 0; zero_at <= value_count; ++zero_at) {
                std::vector<silex::ResidueFieldElement> mixed;
                for (slong k = 0; k <= value_count; ++k) {
                    mixed.emplace_back(residue_field);
                    if (k == zero_at) {
                        assert(mixed.back().zero());
                        assert(!quotient_log.apply(sflint::FmpzRef(expected),
                                                   mixed.back()));
                        continue;
                    }
                    const slong source = k < zero_at ? k : k - 1;
                    assert(mixed.back().set_element(*values[source]));
                }
                sflint::FmpzMat mixed_column(value_count + 1, 1);
                const bool rejected =
                        !silex::detail::residue_quotient_log_column(
                                mixed_column, mixed, residue_field,
                                quotient_log, sflint::FmpzConstRef(ell),
                                nullptr);
                assert(rejected);
                if (!rejected) {
                    return 1;
                }
            }
            ++checked;
        }
    }
    assert(checked == 9);
    return checked == 9 ? 0 : 1;
}

// word_quotient_log_resolve takes the shared walk when ell has fewer than 7
// bits (kWordQuotientLogBabyGiantBits) and baby-step giant-step from 7 bits
// on.  The route counters of ComputeTestSeams, read through TestHookAccess,
// show which one residue_quotient_log_column ran: ell = 3 and ell = 61 (6
// bits, a prime just below the cutoff 64) walk; ell = 67
// (7 bits, just above the cutoff 64) and ell = 127 (7 bits, the largest) use
// baby-step giant-step.  Each case uses the first prime q = k ell + 1 that
// splits completely in the cubic, so ell divides q - 1.
int test_residue_quotient_log_column_route_follows_the_bit_cutoff() {
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order order;
    order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    silex::Element theta_plus_two(cubic);
    assert(theta.gen());
    assert(theta_plus_two.add_si(theta, 2));
    const silex::Element* values[] = {&theta, &theta_plus_two};

    struct Case {
        slong ell;
        bool baby_giant;
    };
    const Case cases[] = {{3, false}, {61, false}, {67, true}, {127, true}};
    for (const Case& entry : cases) {
        sflint::Fmpz ell;
        assert(set_fmpz_si(ell, entry.ell));
        silex::PrimeIdealList primes;
        for (ulong k = 1; k < 10000; ++k) {
            const ulong q = k * static_cast<ulong>(entry.ell) + 1;
            if (n_is_prime(q) == 0) {
                continue;
            }
            sflint::Fmpz q_fmpz;
            assert(set_fmpz_si(q_fmpz, static_cast<slong>(q)));
            primes = silex::PrimeIdealList();
            if (silex::decompose_prime(primes, order,
                                       sflint::FmpzConstRef(q_fmpz), 1) &&
                primes.size() == 3) {
                break;
            }
            primes = silex::PrimeIdealList();
        }
        assert(primes.size() == 3);
        if (primes.size() != 3) {
            return 1;
        }
        const silex::PrimeIdeal* prime = primes.at(0);
        assert(prime != nullptr);
        if (prime == nullptr) {
            return 1;
        }
        silex::ResidueField residue_field(*prime);
        silex::ResidueFieldQuotientLog quotient_log(residue_field);
        assert(residue_field.is_defined() && residue_field.degree() == 1);
        assert(quotient_log.set_ell(sflint::FmpzConstRef(ell)));
        std::vector<silex::ResidueFieldElement> images;
        for (const silex::Element* value : values) {
            images.emplace_back(residue_field);
            assert(images.back().set_element(*value));
        }
        sflint::FmpzMat column(2, 1);

        const unsigned long walk_before =
                silex::detail::TestHookAccess::word_quotient_log_walk_calls();
        const unsigned long bsgs_before = silex::detail::TestHookAccess::
                word_quotient_log_baby_giant_calls();
        const bool ok = silex::detail::residue_quotient_log_column(
                column, images, residue_field, quotient_log,
                sflint::FmpzConstRef(ell), nullptr);
        const unsigned long walk_after =
                silex::detail::TestHookAccess::word_quotient_log_walk_calls();
        const unsigned long bsgs_after = silex::detail::TestHookAccess::
                word_quotient_log_baby_giant_calls();
        assert(ok);
        assert(walk_after - walk_before == (entry.baby_giant ? 0u : 1u));
        assert(bsgs_after - bsgs_before == (entry.baby_giant ? 1u : 0u));
        if (!ok || walk_after - walk_before != (entry.baby_giant ? 0u : 1u) ||
            bsgs_after - bsgs_before != (entry.baby_giant ? 1u : 0u)) {
            return 1;
        }
    }
    return 0;
}

int test_direct_degree_one_residue_dlog_matrix_matches_residue_field() {
    silex::OrderUnitGroup group = real_quadratic_unit_group();
    const silex::Order* order = group.parent();
    assert(order != nullptr);
    if (order == nullptr) {
        return 1;
    }

    sflint::Fmpz p;
    sflint::Fmpz ell;
    const bool prime_set = set_fmpz_si(p, 7) && set_fmpz_si(ell, 2);
    assert(prime_set);
    if (!prime_set) {
        return 1;
    }

    silex::PrimeIdealList primes;
    const bool decomposed = silex::decompose_prime(
            primes, *order, sflint::FmpzConstRef(p), 1);
    assert(decomposed);
    assert(primes.size() > 0);
    if (!decomposed || primes.size() == 0) {
        return 1;
    }

    slong checked = 0;
    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        if (prime == nullptr) {
            return 1;
        }

        sflint::FmpzMat generic_column(0, 0);
        sflint::FmpzMat direct_column(0, 0);
        sflint::FmpzMat saturation_column(0, 0);
        const bool columns_match =
                generic_residue_dlog_column_for_test(
                        generic_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                silex::detail::residue_dlog_matrix_direct_degree_one(
                        direct_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                silex::detail::saturation_prime_column(
                        saturation_column, group, *prime,
                        sflint::FmpzConstRef(ell)) &&
                sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(generic_column),
                        sflint::FmpzMatConstRef(direct_column)) &&
                sflint::fmpz_mat_equal(
                        sflint::FmpzMatConstRef(generic_column),
                        sflint::FmpzMatConstRef(saturation_column));
        assert(columns_match);
        if (!columns_match) {
            return 1;
        }
        ++checked;
    }
    assert(checked == primes.size());
    return checked == primes.size() ? 0 : 1;
}

int test_half_integral_direct_degree_one_residue_dlog_matches_residue_field() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));

    silex::OrderUnitGroup group;
    assert(group.compute(maximal));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 11));
    assert(set_fmpz_si(ell, 2));

    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(
            primes, maximal, sflint::FmpzConstRef(p), 1));
    assert(primes.size() == 2);

    for (slong i = 0; i < primes.size(); ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);

        sflint::FmpzMat generic_column(0, 0);
        sflint::FmpzMat direct_column(0, 0);
        sflint::FmpzMat saturation_column(0, 0);
        assert(generic_residue_dlog_column_for_test(
                generic_column, group, *prime,
                sflint::FmpzConstRef(ell)));
        assert(silex::detail::residue_dlog_matrix_direct_degree_one(
                direct_column, group, *prime,
                sflint::FmpzConstRef(ell)));
        assert(silex::detail::saturation_prime_column(
                saturation_column, group, *prime,
                sflint::FmpzConstRef(ell)));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(generic_column),
                sflint::FmpzMatConstRef(direct_column)));
        assert(sflint::fmpz_mat_equal(
                sflint::FmpzMatConstRef(generic_column),
                sflint::FmpzMatConstRef(saturation_column)));
    }

    return 0;
}

int test_select_saturation_proof_primes() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 31));

    silex::PrimeIdealList selected;
    sflint::FmpzMat kernel(1, 1);
    bool certified = false;
    assert(group.select_saturation_proof_primes(
            selected, certified, kernel, sflint::FmpzConstRef(ell),
            sflint::FmpzConstRef(bound)));
    assert(certified);
    assert(selected.size() > 0);
    assert(sflint::fmpz_mat_nrows(kernel) == 0);
    assert(sflint::fmpz_mat_ncols(kernel) == 2);
    for (slong i = 0; i < selected.size(); ++i) {
        assert(selected.at(i) != nullptr);
        assert(silex::same_order_parent(selected.at(i)->parent(), &order));
        assert(selected.at(i)->residue_degree() == 1);
    }

    assert(set_fmpz_si(bound, 2));
    certified = false;
    sflint::FmpzMat preserved_kernel(1, 1);
    silex::PrimeIdealList preserved;
    assert(group.select_saturation_proof_primes(
            preserved, certified, preserved_kernel,
            sflint::FmpzConstRef(ell), sflint::FmpzConstRef(bound)) == false);
    assert(!certified);
    assert(!preserved.is_defined());
    assert(sflint::fmpz_mat_nrows(preserved_kernel) == 1);
    assert(sflint::fmpz_mat_ncols(preserved_kernel) == 1);

    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));
    assert(set_fmpz_si(ell, 3));
    assert(set_fmpz_si(bound, 5));
    assert(rank_zero.select_saturation_proof_primes(
            preserved, certified, preserved_kernel,
            sflint::FmpzConstRef(ell), sflint::FmpzConstRef(bound)));
    assert(certified);
    assert(preserved.size() == 0);
    assert(sflint::fmpz_mat_nrows(preserved_kernel) == 0);
    assert(sflint::fmpz_mat_ncols(preserved_kernel) == 0);

    return 0;
}

int test_unit_proof_records() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::OrderUnitGroup group;
    assert(group.compute(order));
    assert(group.unit_proof_record_count() == 0);

    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(aux_bound, 11));
    assert(!group.unit_proof_verified(sflint::FmpzConstRef(ell)));

    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_local_saturated(
            status, changed, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 80));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.unit_proof_record_count() == 1);
    assert(proved.unit_proof_verified(sflint::FmpzConstRef(ell)));
    assert(check_unit_proof_record(proved, 0, 2,
                                   silex::ProofState::verified, 11, -1,
                                   false));
    auto proof_record = proved.unit_proof_record(0);
    assert(proof_record.has_value());
    assert(sflint::fmpz_equal_si(proof_record->ell, 2));
    assert(proof_record->status == silex::ProofState::verified);
    assert(sflint::fmpz_equal_si(proof_record->aux_prime_bound, 11));
    assert(proof_record->local_primes > 0);
    assert(!proof_record->changed);

    assert(set_fmpz_si(aux_bound, 13));
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_local_saturated(
            status, changed, proved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 80));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.unit_proof_record_count() == 1);
    assert(check_unit_proof_record(proved, 0, 2,
                                   silex::ProofState::verified, 13, -1,
                                   false));

    assert(set_fmpz_si(ell, 3));
    assert(set_fmpz_si(aux_bound, 17));
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_local_saturated(
            status, changed, proved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 80));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.unit_proof_record_count() == 2);
    assert(check_unit_proof_record(proved, 1, 3,
                                   silex::ProofState::verified, 17, 0,
                                   false));

    silex::OrderUnitGroup copied(order);
    assert(copied.set(proved));
    assert(copied.unit_proof_record_count() == 2);
    assert(copied.unit_proof_verified(sflint::FmpzConstRef(ell)));

    silex::OrderUnitGroup empty(order);
    assert(empty.compute(order));
    assert(empty.unit_proof_record_count() == 0);
    copied.swap(empty);
    assert(copied.unit_proof_record_count() == 0);
    assert(empty.unit_proof_record_count() == 2);
    assert(empty.unit_proof_verified(sflint::FmpzConstRef(ell)));

    assert(set_fmpz_si(ell, 5));
    assert(!copied.unit_proof_verified(sflint::FmpzConstRef(ell)));
    slong local_primes = 17;
    status = silex::ProofState::verified;
    changed = true;
    assert(!copied.unit_proof_record(sflint::FmpzRef(ell), status,
                                     sflint::FmpzRef(aux_bound),
                                     local_primes, changed, 99));
    assert(sflint::fmpz_equal_si(ell, 5));
    assert(sflint::fmpz_equal_si(aux_bound, 17));
    assert(status == silex::ProofState::verified);
    assert(local_primes == 17);
    assert(changed);
    assert(!copied.unit_proof_record(99).has_value());

    return 0;
}

int test_prove_local_saturated_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order,
                           silex::FactoredElementSpan(primitive_generators,
                                                      1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 31));

    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_local_saturated(
            status, changed, group, sflint::FmpzConstRef(ell), 2,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(check_first_free_generator(proved, epsilon));
    assert(proved.certification_status() == silex::CertificationMode::unknown);
    assert(check_unit_proof_record(proved, 0, 2,
                                   silex::ProofState::verified, 31, -1,
                                   false));

    silex::Element epsilon2(field);
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement square_generators[] = {std::move(square)};
    silex::OrderUnitGroup square_group(order);
    assert(square_group.set_units(order,
                                  silex::FactoredElementSpan(square_generators,
                                                             1),
                                  embeddings, 128));

    silex::OrderUnitGroup mutated(order);
    status = silex::ProofState::not_checked;
    changed = false;
    assert(mutated.prove_local_saturated(
            status, changed, square_group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::unavailable);
    assert(changed);
    assert(check_first_free_generator_square(mutated, epsilon2));
    assert(check_unit_proof_record(mutated, 0, 2,
                                   silex::ProofState::unavailable, 31, -1,
                                   true));

    assert(set_fmpz_si(bound, 2));
    silex::OrderUnitGroup unavailable(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(unavailable.prove_local_saturated(
            status, changed, square_group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::unavailable);
    assert(!changed);
    assert(check_first_free_generator(unavailable, epsilon2));
    assert(check_unit_proof_record(unavailable, 0, 2,
                                   silex::ProofState::unavailable, 2, 0,
                                   false));

    sflint::Fmpz not_prime;
    assert(set_fmpz_si(not_prime, 4));
    status = silex::ProofState::verified;
    changed = true;
    assert(!unavailable.prove_local_saturated(
            status, changed, square_group, sflint::FmpzConstRef(not_prime), 1,
            sflint::FmpzConstRef(bound), embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(changed);

    return 0;
}

int test_regulator_index_bound() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz bound;
    assert(rank_zero.regulator_index_bound(sflint::FmpzRef(bound), 128));
    assert(sflint::fmpz_equal_si(bound, 1));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup full;
    assert(full.compute(order));
    // floor(log(1 + sqrt2) / 0.48) = floor(1.836...) = 1, with 0.48 the
    // Friedman 1989 Table 6 bound for signature (2, 0).
    assert(full.regulator_index_bound(sflint::FmpzRef(bound), 160));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(bound), 1));

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(square)};
    silex::OrderUnitGroup subgroup(order);
    assert(subgroup.set_units(order, silex::FactoredElementSpan(generators, 1),
                              embeddings, 160));
    sflint::Fmpz subgroup_bound;
    assert(subgroup.regulator_index_bound(sflint::FmpzRef(subgroup_bound),
                                          160));
    // floor(2 log(1 + sqrt2) / 0.48) = floor(3.672...) = 3.
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(subgroup_bound), 3));

    // Q(sqrt7): d_K = 28 and R_K = log(8 + 3 sqrt7) = 2.76865...  For the
    // maximal order the Pohst-Zassenhaus (6.22) bound log(28/4)/2 =
    // 0.97295... exceeds Table 6's 0.48: floor(2.76865... / 0.97295...) = 2.
    // The order of conductor 3 has the same unit group, but it is not
    // maximal, so it keeps the signature bound: floor(2.76865... / 0.48) = 5.
    silex::NumberField sqrt7 = quadratic_field(7);
    silex::Order sqrt7_maximal = silex::test::equation_order(sqrt7);
    assert(sqrt7_maximal.is_maximal());
    silex::OrderUnitGroup sqrt7_units;
    assert(sqrt7_units.compute(sqrt7_maximal));
    assert(sqrt7_units.regulator_index_bound(sflint::FmpzRef(bound), 160));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(bound), 2));
    sflint::Fmpz three;
    assert(set_fmpz_si(three, 3));
    silex::Order sqrt7_conductor3 = silex::Order::quadratic_order(
            sqrt7, sflint::FmpzConstRef(three));
    assert(sqrt7_conductor3.maximality_known());
    assert(!sqrt7_conductor3.is_maximal());
    silex::Element sqrt7_unit(sqrt7);
    assert(set_quadratic_coeffs(sqrt7_unit, 8, 1, 3, 1));
    silex::FactoredElement sqrt7_generator(sqrt7);
    assert(sqrt7_generator.set_element(sqrt7_unit));
    silex::FactoredElement sqrt7_generators[] = {std::move(sqrt7_generator)};
    silex::EmbeddingContext sqrt7_embeddings(sqrt7);
    silex::OrderUnitGroup conductor3_units(sqrt7_conductor3);
    assert(conductor3_units.set_units(
            sqrt7_conductor3,
            silex::FactoredElementSpan(sqrt7_generators, 1),
            sqrt7_embeddings, 160));
    assert(conductor3_units.regulator_index_bound(sflint::FmpzRef(bound), 160));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(bound), 5));

    sflint::Fmpz sentinel;
    assert(set_fmpz_si(sentinel, 42));
    silex::OrderUnitGroup unset(order);
    assert(!unset.regulator_index_bound(sflint::FmpzRef(sentinel), 160));
    assert(sflint::fmpz_equal_si(sentinel, 42));
    assert(!rank_zero.regulator_index_bound(sflint::FmpzRef(sentinel), 0));
    assert(sflint::fmpz_equal_si(sentinel, 42));

    return 0;
}

int test_class_regulator_index_bound() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderUnitGroup units;
    assert(units.compute(order));

    auto prepare_context = [&](silex::ClassGroupContext& context) noexcept {
        sflint::Fmpz bound;
        assert(set_fmpz_si(bound, 2));
        context = silex::ClassGroupContext(order);
        assert(context.is_defined());
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));

        const silex::FactorBase* base = context.factor_base();
        assert(base != nullptr);
        silex::Relation relation(*base);
        silex::Element alpha(field);
        assert(alpha.set_si(-2));
        assert(relation.set_generator(alpha));
        assert(context.append_relation(relation));
        assert(context.publish_presentation());
        assert(context.has_presentation());
    };

    silex::ClassGroupContext class_group;
    prepare_context(class_group);

    sflint::Arb product;
    assert(units.class_regulator_product(sflint::ArbRef(product),
                                         class_group, 128));
    assert(sflint::arb_contains_si(product, 1));

    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb analytic_hR;
    sflint::arb_set_fmpq(analytic_hR, half, 128);

    sflint::Fmpz index_bound;
    assert(units.class_regulator_index_bound(
            sflint::FmpzRef(index_bound), class_group,
            sflint::ArbConstRef(analytic_hR), 128));
    assert(sflint::fmpz_equal_si(index_bound, 2));

    sflint::Fmpz sentinel;
    assert(set_fmpz_si(sentinel, 42));
    sflint::Arb zero;
    sflint::arb_zero(zero);
    assert(!units.class_regulator_index_bound(
            sflint::FmpzRef(sentinel), class_group,
            sflint::ArbConstRef(zero), 128));
    assert(sflint::fmpz_equal_si(sentinel, 42));
    assert(!units.class_regulator_index_bound(
            sflint::FmpzRef(sentinel), class_group,
            sflint::ArbConstRef(analytic_hR), 0));
    assert(sflint::fmpz_equal_si(sentinel, 42));

    // An analytic hR of 2 against the candidate product 1 puts the quotient
    // below one; the bound fails closed rather than reporting index one.
    sflint::Arb large_analytic_hR;
    sflint::arb_set_si(large_analytic_hR, 2);
    assert(!units.class_regulator_index_bound(
            sflint::FmpzRef(sentinel), class_group,
            sflint::ArbConstRef(large_analytic_hR), 128));
    assert(sflint::fmpz_equal_si(sentinel, 42));

    return 0;
}

int test_class_regulator_index_bound_interval_boundary() {
    sflint::Arb candidate_product;
    sflint::Arb analytic_product;
    sflint::Arb analytic_error;
    sflint::arb_set_si(candidate_product, 2);
    sflint::arb_set_si(analytic_product, 2);
    sflint::arb_set_si(analytic_error, 1);
    sflint::arb_add_error(analytic_product, analytic_error);
    assert(sflint::arb_contains_si(analytic_product, 1));
    assert(sflint::arb_contains_si(analytic_product, 3));

    sflint::Fmpz bound;
    assert(silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(analytic_product), 128, nullptr));
    assert(sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) >= 0);

    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb narrow_analytic_product;
    sflint::Arb narrow_error;
    sflint::arb_set_si(narrow_analytic_product, 2);
    sflint::arb_set_fmpq(narrow_error, half, 128);
    sflint::arb_add_error(narrow_analytic_product, narrow_error);
    assert(silex::detail::
                   class_regulator_index_is_one_from_candidate_product(
                           sflint::ArbConstRef(candidate_product),
                           sflint::ArbConstRef(narrow_analytic_product), 128,
                           nullptr));
    assert(silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(narrow_analytic_product), 128, nullptr));
    assert(sflint::fmpz_equal_si(bound, 2));

    return 0;
}

// The quotient h_cand R_cand / (hR)_analytic is a positive integer (the
// combined class/unit index) when the analytic enclosure is correct.  An
// upper endpoint below one leaves no positive integer in the enclosure, so
// the bound fails closed instead of reporting index one (decision
// 2026-09-26).
int test_class_regulator_index_bound_below_one_fails_closed() {
    sflint::Arb candidate_product;
    sflint::Arb analytic_product;
    sflint::Fmpz bound;

    // 1 / 2: upper endpoint below one, fail closed, output unchanged.
    sflint::arb_set_si(candidate_product, 1);
    sflint::arb_set_si(analytic_product, 2);
    assert(set_fmpz_si(bound, 42));
    assert(!silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(analytic_product), 128, nullptr));
    assert(sflint::fmpz_equal_si(bound, 42));

    // 1 / (2 +/- 1/2): the enclosure [2/5, 2/3] is still below one.
    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb error;
    sflint::arb_set_fmpq(error, half, 128);
    sflint::arb_add_error(analytic_product, error);
    assert(!silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(analytic_product), 128, nullptr));
    assert(sflint::fmpz_equal_si(bound, 42));

    // An exact quotient of one still gives one.
    sflint::arb_set_si(analytic_product, 1);
    assert(silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound),
            sflint::ArbConstRef(candidate_product),
            sflint::ArbConstRef(analytic_product), 128, nullptr));
    assert(sflint::fmpz_equal_si(bound, 1));

    return 0;
}

// Both index bounds fail closed only when the rigorous upper endpoint of the
// quotient is below one (decision 2026-09-26).  This builds a numerator whose
// quotient by an exact one is a nonzero-width enclosure [1 - 2r, 1]: midpoint
// and lower endpoint below one, upper endpoint exactly one.  It must give
// bound one, not fail.  arb_div rounds the radius up (mag_div), so the
// midpoint is set to 1 - r, where r is the radius arb_div produces from the
// numerator radius 1/4; both bounds divide exactly so at this precision.
int test_index_bounds_accept_upper_endpoint_exactly_one() {
    const slong prec = 128;
    sflint::Arb one;
    sflint::arb_set_si(one, 1);
    sflint::Fmpq quarter;
    sflint::fmpq_set_si(quarter, 1, 4);
    sflint::Arb error;
    sflint::arb_set_fmpq(error, quarter, prec);
    sflint::Arb probe;
    sflint::arb_set_si(probe, 0);
    sflint::arb_add_error(probe, error);
    sflint::Arb probe_quotient;
    sflint::arb_div(probe_quotient, probe, one, prec);

    sflint::Arf quotient_radius;
    ::arf_set_mag(quotient_radius.raw(), arb_radref(probe_quotient.raw()));
    sflint::Arf midpoint;
    ::arf_one(midpoint.raw());
    assert(::arf_sub(midpoint.raw(), midpoint.raw(), quotient_radius.raw(),
                     prec, ARF_RND_DOWN) == 0);
    sflint::Arb numerator;
    ::arb_set_arf(numerator.raw(), midpoint.raw());
    ::mag_set(arb_radref(numerator.raw()), arb_radref(probe.raw()));

    // Precondition: the quotient enclosure the bounds compute has nonzero
    // radius, contains 1/2, and has upper endpoint exactly one.
    sflint::Arb quotient;
    sflint::arb_div(quotient, numerator, one, prec);
    assert(!::arb_is_exact(quotient.raw()));
    assert(::arb_contains_si(quotient.raw(), 1));
    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    assert(::arb_contains_fmpq(quotient.raw(), half.raw()));
    sflint::Arf upper;
    sflint::arb_get_ubound_arf(upper, quotient, prec);
    assert(::arf_equal_si(upper.raw(), 1));

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 42));
    assert(silex::detail::class_regulator_index_bound_from_candidate_product(
            sflint::FmpzRef(bound), sflint::ArbConstRef(numerator),
            sflint::ArbConstRef(one), prec, nullptr));
    assert(sflint::fmpz_equal_si(bound, 1));

    assert(set_fmpz_si(bound, 42));
    assert(silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(numerator),
            sflint::ArbConstRef(one), prec));
    assert(sflint::fmpz_equal_si(bound, 1));

    return 0;
}

// The unit index bound is the floor of the upper endpoint of R_G / R_lower,
// and fails closed when that endpoint is below one (decision 2026-09-26).
int test_unit_index_bound_from_regulator_quotient() {
    sflint::Arb regulator;
    sflint::Arb lower;
    sflint::Fmpz bound;

    // 7 / 2: floor(3.5) = 3.
    sflint::arb_set_si(regulator, 7);
    sflint::arb_set_si(lower, 2);
    assert(silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 128));
    assert(sflint::fmpz_equal_si(bound, 3));

    // An exact quotient of one gives one.
    sflint::arb_set_si(regulator, 2);
    assert(silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 128));
    assert(sflint::fmpz_equal_si(bound, 1));

    // An enclosure 1 +/- 1/2 of the quotient has upper endpoint 3/2, so the
    // bound is floor(3/2) = 1.
    sflint::Fmpq half;
    sflint::fmpq_set_si(half, 1, 2);
    sflint::Arb error;
    sflint::arb_set_fmpq(error, half, 128);
    sflint::arb_set_si(regulator, 2);
    sflint::arb_add_error(regulator, error);
    sflint::arb_set_si(lower, 2);
    assert(silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 128));
    assert(sflint::fmpz_equal_si(bound, 1));

    // A subgroup regulator below the regulator lower bound (1 / 2 < 1)
    // contradicts the lower bound: fail closed, output unchanged.
    assert(set_fmpz_si(bound, 42));
    sflint::arb_set_si(regulator, 1);
    assert(!silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 128));
    assert(sflint::fmpz_equal_si(bound, 42));

    // Non-positive inputs and a nonpositive precision fail.
    sflint::arb_set_si(regulator, 0);
    assert(!silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 128));
    sflint::arb_set_si(regulator, 7);
    sflint::arb_set_si(lower, 0);
    assert(!silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 128));
    sflint::arb_set_si(lower, 2);
    assert(!silex::detail::unit_index_bound_from_regulator_quotient(
            sflint::FmpzRef(bound), sflint::ArbConstRef(regulator),
            sflint::ArbConstRef(lower), 0));
    assert(sflint::fmpz_equal_si(bound, 42));

    return 0;
}

// Gates that trust a caller-supplied analytic class-regulator product,
// saturation index bound, or exact norm are internal.  A requires-expression
// reports false for a private member, so these fail to compile if any of the
// gates becomes public again.
template <typename Context>
constexpr bool has_public_class_unit_with_units_gate =
        requires(Context& context,
                 silex::OrderUnitGroup& units,
                 const sflint::Arb& analytic_hR) {
            context.try_certify_class_unit_with_units(
                    units, sflint::ArbConstRef(analytic_hR), slong{1});
        };

template <typename Context>
constexpr bool has_public_analytic_index_bound_gate =
        requires(Context& context,
                 const silex::OrderUnitGroup& units,
                 const sflint::Arb& analytic_hR,
                 const sflint::Fmpz& aux_bound) {
            context.try_analytic_index_bound_with_units(
                    units, sflint::ArbConstRef(analytic_hR),
                    sflint::FmpzConstRef(aux_bound), slong{1});
        };

template <typename Context>
constexpr bool has_public_saturation_index_bound_gate =
        requires(Context& context,
                 const silex::OrderUnitGroup& units,
                 const sflint::Fmpz& index_bound,
                 const sflint::Fmpz& aux_bound) {
            context.try_prove_relation_saturation_index_bound_with_units(
                    units, sflint::FmpzConstRef(index_bound),
                    sflint::FmpzConstRef(aux_bound));
        };

template <typename Context>
constexpr bool has_public_generator_with_norm_gate =
        requires(Context& context,
                 bool& partial_throttle_exit,
                 const silex::Element& generator,
                 const sflint::Fmpq& norm) {
            context.try_append_generator_relation_with_norm(
                    partial_throttle_exit, generator,
                    sflint::FmpqConstRef(norm));
        };

// Control: the zeta gate, which computes hR itself, stays public.
template <typename Context>
constexpr bool has_public_class_unit_with_zeta_gate =
        requires(Context& context, silex::OrderUnitGroup& units) {
            context.try_certify_class_unit_with_zeta(units, slong{1});
        };

static_assert(
        has_public_class_unit_with_zeta_gate<silex::ClassGroupContext>);
static_assert(
        !has_public_class_unit_with_units_gate<silex::ClassGroupContext>);
static_assert(
        !has_public_analytic_index_bound_gate<silex::ClassGroupContext>);
static_assert(
        !has_public_saturation_index_bound_gate<silex::ClassGroupContext>);
static_assert(
        !has_public_generator_with_norm_gate<silex::ClassGroupContext>);

bool real_quadratic_two_candidate(silex::ClassGroupContext& out,
                                  const silex::Order& order) noexcept {
    sflint::Fmpz bound;
    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    return set_fmpz_si(bound, 2) &&
           out.compute_candidate(order, sflint::FmpzConstRef(bound),
                                 options) &&
           out.has_presentation() &&
           out.certification_status() == silex::CertificationMode::unknown;
}

bool real_quadratic_two_units(silex::OrderUnitGroup& out,
                              const silex::Order& order,
                              silex::EmbeddingContext& embeddings,
                              const silex::Element& unit) noexcept {
    const silex::NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }
    silex::FactoredElement generator(*field);
    if (!generator.set_element(unit)) {
        return false;
    }
    silex::FactoredElement generators[] = {std::move(generator)};
    return out.set_units(order, silex::FactoredElementSpan(generators, 1),
                         embeddings, 256) &&
           out.certification_status() == silex::CertificationMode::unknown;
}

// The finite-index loop of reference Clgp.jl:_class_unit_group saturates and
// then re-runs _validate_class_unit_group, which reads the class number of the
// current relation module (class_group_current_h).  On Z[sqrt(2)] with the
// factor base over the primes up to 7 and the relations sqrt(2), 3, 5,
// 11 + 6 sqrt(2) = (3 + sqrt(2))^2 and 3 - sqrt(2), the tentative class group
// is Z/2.  Saturation at 2 finds 3 + sqrt(2) and appends it as a class
// relation, which clears the stored presentation.  Both saturation helpers
// must validate against the rebuilt presentation instead of reading
// "presentation unavailable" as a result.  With factor base generation
// recorded as checked the real index bound drops to 1.  Without it the bound
// stays 2.  An outcome the small-index loop cannot resolve after a rebuild
// keeps Silex's behaviour from before the rebuild instead of failing the
// compute: a round that rebuilt the presentation sets the candidate validated
// (in the proven path the full proven re-check then decides), and a loop
// entered after saturation at 2 rebuilt it returns without an attempt.
enum class AppendedClassRelationCase {
    small_index_generation_checked,
    small_index_bound_two,
    small_index_bound_three,
    at_two_generation_checked,
    at_two_then_small_index,
};

bool appended_class_relation_case(AppendedClassRelationCase which) {
    using CertificationAccess =
            silex::detail::ClassGroupCertificationAccess;
    using silex::detail::ValidateRefineOutcome;
    const bool record_generation_check =
            which == AppendedClassRelationCase::
                             small_index_generation_checked ||
            which == AppendedClassRelationCase::at_two_generation_checked;
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());

    silex::ClassGroupContext context(order);
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 7));
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    assert(context.generator_count() == 5);
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);
    const slong coefficients[5][2] = {
            {0, 1}, {3, 0}, {5, 0}, {11, 6}, {3, -1}};
    for (const auto& pair : coefficients) {
        silex::Relation relation(*base);
        silex::Element generator(field);
        assert(set_quadratic_coeffs(generator, pair[0], 1, pair[1], 1));
        assert(relation.set_generator(generator));
        assert(context.append_relation(relation));
    }
    assert(context.publish_presentation());
    if (record_generation_check) {
        assert(CertificationAccess::record_factor_base_honesty_proof(
                context, sflint::FmpzConstRef(bound)));
    }
    sflint::Fmpz class_order;
    assert(context.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 2));

    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::OrderUnitGroup units(order);
    assert(real_quadratic_two_units(units, order, embeddings, epsilon));
    assert(units.free_rank() == 1);

    silex::detail::RelationUnitExtractionState state;
    silex::detail::ValidateRefineSummary summary;
    silex::detail::AnalyticClassRegulatorCache analytic_cache;
    const slong relations_before = context.relation_count();
    bool validated = false;
    bool attempted = false;
    bool entered_after_rebuild = false;
    if (which == AppendedClassRelationCase::at_two_generation_checked ||
        which == AppendedClassRelationCase::at_two_then_small_index) {
        bool rebuilt = false;
        assert(silex::detail::unknown_saturation_at_two(
                validated, rebuilt, state, summary, context, units,
                embeddings, analytic_cache, order, 1, 128, nullptr, false));
        assert(summary.outcome !=
               ValidateRefineOutcome::local_rank_unavailable);
        assert(state.saturation_at_two_done);
        assert(rebuilt);
        assert(context.relation_count() > relations_before);
        assert(context.has_presentation());
        assert(context.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 1));
        if (record_generation_check) {
            assert(validated);
            assert(summary.outcome == ValidateRefineOutcome::proven);
            assert(sflint::fmpz_equal_si(summary.last_index_bound, 1));
            return true;
        }
        assert(!validated);
        assert(sflint::fmpz_equal_si(summary.last_index_bound, 2));
        entered_after_rebuild = true;
    } else {
        const slong entry_bound =
                which == AppendedClassRelationCase::small_index_bound_three
                        ? 3
                        : 2;
        ::fmpz_set_si(summary.last_index_bound.raw(), entry_bound);
    }

    const slong relations_entering = context.relation_count();
    assert(silex::detail::unknown_small_index_saturation(
            validated, attempted, state, summary, context, units, embeddings,
            analytic_cache, order, 1, 128, nullptr, false,
            entered_after_rebuild));
    assert(summary.outcome != ValidateRefineOutcome::local_rank_unavailable);
    assert(context.has_presentation());
    assert(context.order(sflint::FmpzRef(class_order)));
    assert(sflint::fmpz_equal_si(class_order, 1));
    if (entered_after_rebuild) {
        // Nothing saturates at 2 or 3 any more; Silex skipped this loop
        // before the rebuild, so the caller continues.
        assert(context.relation_count() == relations_entering);
        assert(!validated);
        assert(!attempted);
        return true;
    }
    assert(attempted);
    assert(context.relation_count() > relations_before);
    assert(validated);
    if (record_generation_check) {
        assert(summary.outcome == ValidateRefineOutcome::proven);
        assert(sflint::fmpz_equal_si(summary.last_index_bound, 1));
    } else {
        // Bound 2 either did not decrease (entered at 2) or decreased from 3
        // with no saturation at 2 or 3 in the next round.
        assert(sflint::fmpz_equal_si(summary.last_index_bound, 2));
    }
    return true;
}

int test_small_index_saturation_revalidates_appended_class_relation() {
    assert(appended_class_relation_case(
            AppendedClassRelationCase::small_index_generation_checked));
    assert(appended_class_relation_case(
            AppendedClassRelationCase::small_index_bound_two));
    assert(appended_class_relation_case(
            AppendedClassRelationCase::small_index_bound_three));
    assert(appended_class_relation_case(
            AppendedClassRelationCase::at_two_generation_checked));
    assert(appended_class_relation_case(
            AppendedClassRelationCase::at_two_then_small_index));
    return 0;
}

// Class/unit validation reads the candidate hR over the analytic hR as an
// integer index, and so recognizes index one, only once factor-base
// generation is checked up to the recorded requirement.  Q(sqrt(-26))
// (h = 6) has requirement floor(sqrt(104 / 3)) = 5.  The public check
// accepts any bound up to the build bound, so a check at 2 marks the
// generation check verified without showing that the base generates the
// class group; validation must then not report index one.  The same
// candidate checked at 5 does.
int test_validation_index_one_requires_generation_to_required_bound() {
    silex::NumberField field = quadratic_field(-26);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup units(order);
    assert(units.set_units(order, silex::FactoredElementSpan(), embeddings,
                           128));
    assert(units.free_rank() == 0);

    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 5));
    silex::ClassGroupCandidateOptions options;
    options.max_candidates = 256;
    options.max_relations = 48;
    for (const slong checked : {2, 5}) {
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), options));
        sflint::Fmpz class_order;
        assert(class_group.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 6));
        sflint::Fmpz required;
        assert(class_group.factor_base_generation_bound(
                sflint::FmpzRef(required)));
        assert(sflint::fmpz_equal_si(required, 5));
        sflint::Fmpz checked_bound;
        assert(set_fmpz_si(checked_bound, checked));
        assert(class_group.check_factor_base_generation_bound(
                sflint::FmpzConstRef(checked_bound)));
        sflint::Fmpz recorded_check;
        assert(class_group.factor_base_generation_checked_bound(
                sflint::FmpzRef(recorded_check)));
        assert(sflint::fmpz_equal_si(recorded_check, checked));

        silex::detail::AnalyticClassRegulatorCache analytic_cache;
        sflint::Fmpz index_bound;
        sflint::Arb expected_regulator;
        const bool estimated = silex::detail::class_unit_validation_estimate(
                index_bound, expected_regulator, units, class_group,
                analytic_cache, order, 128, false);
        const bool index_one =
                estimated &&
                sflint::fmpz_is_one(sflint::FmpzConstRef(index_bound));
        assert(index_one == (checked >= 5));
    }
    return 0;
}

int test_class_unit_regulator_certification() {
    using CertificationAccess =
            silex::detail::ClassGroupCertificationAccess;

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));

    // Public route: the zeta gate computes hR itself and certifies the
    // fundamental unit 1 + sqrt(2), R = log(1 + sqrt(2)) = 0.88137...
    {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units(order);
        assert(real_quadratic_two_candidate(class_group, order));
        assert(real_quadratic_two_units(units, order, embeddings, epsilon));
        assert(class_group.try_certify_class_unit_with_zeta(units, 256));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::verified);
        // The quadratic L(1, chi) route is unconditional.
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::proven);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.regulator_proof_status() ==
               silex::ProofState::verified);
        // Index one against the unconditional hR proves saturation at
        // every prime; no ell-local records are kept.
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(class_group.relation_saturation_record_count() == 0);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);
        // Q(sqrt 2) has class number 1 and torsion {+1, -1}.
        sflint::Fmpz class_order;
        assert(class_group.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, 1));
        sflint::Fmpz torsion_order;
        assert(units.torsion_order(sflint::FmpzRef(torsion_order)));
        assert(sflint::fmpz_equal_si(torsion_order, 2));
        sflint::Arb regulator;
        assert(units.regulator(sflint::ArbRef(regulator)));
        sflint::Arb lower;
        sflint::Arb upper;
        sflint::arb_set_d(lower, 0.8813);
        sflint::arb_set_d(upper, 0.8814);
        assert(sflint::arb_gt(regulator, lower));
        assert(sflint::arb_lt(regulator, upper));
    }

    // Internal gate: the access-class entry point certifies when it is fed
    // the independently computed zeta value of hR, as the internal
    // validation and compute routes do.
    {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units(order);
        assert(real_quadratic_two_candidate(class_group, order));
        assert(real_quadratic_two_units(units, order, embeddings, epsilon));
        sflint::Arb zeta_hR;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(zeta_hR),
                                                   order, 256));
        // h = 1, so the zeta value is R = log(1 + sqrt(2)) =
        // 0.8813735870195430252...: it overlaps an independent Arb
        // evaluation and lies inside 0.881373587019543 +/- 1e-15.
        sflint::Arb reference;
        ::arb_sqrt_ui(reference.raw(), 2, 256);
        ::arb_add_ui(reference.raw(), reference.raw(), 1, 256);
        ::arb_log(reference.raw(), reference.raw(), 256);
        assert(::arb_overlaps(zeta_hR.raw(), reference.raw()) != 0);
        sflint::Arb decimal_window;
        assert(::arb_set_str(decimal_window.raw(),
                             "0.881373587019543 +/- 1e-15", 256) == 0);
        assert(::arb_contains(decimal_window.raw(), zeta_hR.raw()) != 0);
        // The same value reported as GRH-conditional (as when the
        // quadratic L(1, chi) evaluation fails and the zeta route falls
        // back to Belabas-Friedman) cannot publish `proven`, although this
        // order has an unconditional route: the gate honours the route that
        // produced the value, not the order type.
        assert(!CertificationAccess::try_certify_class_unit_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR), 256,
                false));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::not_checked);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(units.certification_status() ==
               silex::CertificationMode::unknown);
        assert(CertificationAccess::try_certify_class_unit_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR), 256,
                true));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::proven);

        // A later Belabas-Friedman audit is recorded but does not downgrade
        // the unconditional L(1, chi) check to `grh`.  (Precision 64: the
        // Belabas-Friedman evaluation does not reach 256 bits within the
        // cutoff on this field.)
        assert(class_group.try_certify_class_unit_with_zeta_bf(units, 20000,
                                                               64));
        assert(class_group.zeta_bf_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::proven);
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);

        // A Belabas-Friedman hR cannot prove other, unproven units, even on
        // a context that is already proven: the gate fails and changes
        // nothing.
        silex::OrderUnitGroup other(order);
        assert(real_quadratic_two_units(other, order, embeddings, epsilon));
        assert(other.certification_status() ==
               silex::CertificationMode::unknown);
        assert(!class_group.try_certify_class_unit_with_zeta_bf(other, 20000,
                                                                64));
        assert(other.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::proven);
    }

    // Negative: the index-two subgroup <3 + 2 sqrt(2)> has regulator
    // 2 log(1 + sqrt(2)) = 1.7627...  Its own hR can no longer be passed to
    // a public gate (see the static_asserts above), and every public route
    // that remains computes hR itself and rejects the subgroup.
    {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units(order);
        assert(real_quadratic_two_candidate(class_group, order));
        assert(real_quadratic_two_units(units, order, embeddings, epsilon2));

        assert(!class_group.try_certify_class_unit_with_zeta(units, 256));
        assert(!class_group.try_certify_with_units(
                units, silex::CertificationMode::proven, 256));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::not_checked);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(units.certification_status() ==
               silex::CertificationMode::unknown);

        // The internal gate, fed the zeta value rather than the subgroup's
        // own product, also rejects it and leaves both objects unchanged.
        sflint::Arb zeta_hR;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(zeta_hR),
                                                   order, 256));
        assert(!CertificationAccess::try_certify_class_unit_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR), 256,
                true));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::not_checked);
        assert(units.certification_status() ==
               silex::CertificationMode::unknown);
    }

    // Negative with proven-labelled units: mark <3 + 2 sqrt(2)> proven
    // through the internal access class, so the unit label cannot be what
    // blocks certification.  The analytic index test (ratio 2) must reject
    // it on both the public zeta gate and the internal gate.
    {
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units(order);
        assert(real_quadratic_two_candidate(class_group, order));
        assert(real_quadratic_two_units(units, order, embeddings, epsilon2));
        silex::detail::TestHookAccess::mark_certification_proven(units);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);
        assert(silex::detail::order_unit_torsion_is_computed(units));

        assert(!class_group.try_certify_class_unit_with_zeta(units, 256));
        sflint::Arb zeta_hR;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(zeta_hR),
                                                   order, 256));
        assert(!CertificationAccess::try_certify_class_unit_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR), 256,
                true));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::not_checked);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(class_group.regulator_proof_status() ==
               silex::ProofState::not_checked);
    }

    return 0;
}

// Cached torsion handed to the internal unit setter must never reach a
// proven label.  The setter rejects any cached pair that is not a cyclic
// group of exactly that order in the order, on both the trusted and the
// validating path; and proven publication recomputes torsion, so an
// under-claimed but well-formed pair (Q(i): w = 2 instead of 4) is still
// rejected by the class/unit gates.
int test_cached_torsion_never_reaches_proven() {
    using OrderUnitGroupAccess = silex::detail::OrderUnitGroupAccess;

    {
        silex::NumberField field = quadratic_field(2);
        silex::Order order;
        order = silex::test::equation_order(field);
        silex::EmbeddingContext embeddings(field);
        silex::Element epsilon(field);
        assert(set_real_quadratic_unit(epsilon));
        silex::FactoredElement generator(field);
        assert(generator.set_element(epsilon));
        silex::FactoredElement generators[] = {std::move(generator)};

        silex::Element theta(field);
        assert(theta.gen());
        struct FakeTorsion {
            slong order;
            slong generator;  // 0 means sqrt(2)
        };
        const FakeTorsion fakes[] = {
                {1, 1},   // w = 1: odd, and -1 is always torsion
                {6, -1},  // no element of order 6 in Q(sqrt 2)
                {8, 0},   // sqrt(2) is not a root of unity
                {4, -1},  // -1 has order 2, not 4
        };
        for (const FakeTorsion& fake : fakes) {
            sflint::Fmpz torsion_order;
            assert(set_fmpz_si(torsion_order, fake.order));
            silex::OrderElement torsion_generator(order);
            if (fake.generator == 0) {
                assert(torsion_generator.set_element(theta));
            } else {
                assert(torsion_generator.set_si(fake.generator));
            }
            for (const bool trusted : {true, false}) {
                silex::OrderUnitGroup units(order);
                assert(!OrderUnitGroupAccess::set_units(
                        units, order,
                        silex::FactoredElementSpan(generators, 1),
                        embeddings, 256, trusted, &torsion_order,
                        &torsion_generator));
                assert(!units.is_set());
            }
        }

        // The honest cached pair (2, -1) is accepted.
        sflint::Fmpz torsion_order;
        assert(set_fmpz_si(torsion_order, 2));
        silex::OrderElement torsion_generator(order);
        assert(torsion_generator.set_si(-1));
        silex::OrderUnitGroup units(order);
        assert(OrderUnitGroupAccess::set_units(
                units, order, silex::FactoredElementSpan(generators, 1),
                embeddings, 256, true, &torsion_order, &torsion_generator));
        assert(silex::detail::order_unit_torsion_is_computed(units));
    }

    {
        silex::NumberField field = quadratic_field(-1);
        silex::Order order;
        order = silex::test::equation_order(field);
        assert(order.is_maximal());
        silex::EmbeddingContext embeddings(field);

        sflint::Fmpz bound;
        assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                    order));
        if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
            sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
        }
        silex::ClassGroupCandidateOptions options;
        options.max_candidates = 256;
        options.max_relations = 48;

        // Under-claimed torsion: w = 2 with generator -1 passes the setter's
        // shape check but is not the torsion of Z[i] (w = 4).
        sflint::Fmpz torsion_order;
        assert(set_fmpz_si(torsion_order, 2));
        silex::OrderElement torsion_generator(order);
        assert(torsion_generator.set_si(-1));
        silex::OrderUnitGroup under(order);
        assert(OrderUnitGroupAccess::set_units(
                under, order, silex::FactoredElementSpan(), embeddings, 256,
                true, &torsion_order, &torsion_generator));
        assert(!silex::detail::order_unit_torsion_is_computed(under));

        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), options));
        assert(!class_group.try_certify_class_unit_with_zeta(under, 256));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(under.certification_status() ==
               silex::CertificationMode::unknown);

        // Control: the same gate certifies the torsion Silex computes.
        silex::OrderUnitGroup computed(order);
        assert(computed.set_units(order, silex::FactoredElementSpan(),
                                  embeddings, 256));
        assert(silex::detail::order_unit_torsion_is_computed(computed));
        sflint::Fmpz computed_order;
        assert(computed.torsion_order(sflint::FmpzRef(computed_order)));
        assert(sflint::fmpz_equal_si(computed_order, 4));
        assert(class_group.try_certify_class_unit_with_zeta(computed, 256));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(computed.certification_status() ==
               silex::CertificationMode::proven);
    }

    return 0;
}

// Belabas-Friedman (2015, Theorem 1) bounds the truncation error of its hR
// evaluation only under GRH for zeta_K and zeta_Q.  For degree >= 3 that is
// the only analytic hR, so it may never publish `proven`: the gates fail
// closed, and a check that is recorded carries the `grh` label.
int test_belabas_friedman_class_regulator_is_grh_conditional() {
    using CertificationAccess =
            silex::detail::ClassGroupCertificationAccess;

    // x^3 - x - 1: h = 1, unit rank 1, and the root theta is a fundamental
    // unit (theta (theta^2 - 1) = 1; R = log theta = 0.2812...).
    silex::NumberField field = cubic_field(-1, -1);
    silex::Order equation;
    equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());
    silex::EmbeddingContext embeddings(field);

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }

    silex::Element theta(field);
    assert(theta.gen());

    // Unproven units and a candidate class group: before this rule the
    // zeta gates published `proven` from the Belabas-Friedman hR (index
    // one).  They now fail closed and leave both objects unchanged.
    {
        silex::ClassGroupCandidateOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), options));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        silex::OrderUnitGroup units(order);
        silex::FactoredElement generator(field);
        assert(generator.set_element(theta));
        silex::FactoredElement generators[] = {std::move(generator)};
        assert(units.set_units(order,
                               silex::FactoredElementSpan(generators, 1),
                               embeddings, 256));
        assert(units.certification_status() ==
               silex::CertificationMode::unknown);

        sflint::Arb zeta_hR;
        assert(silex::zeta_class_regulator_product(sflint::ArbRef(zeta_hR),
                                                   order, 256));

        assert(!class_group.try_certify_class_unit_with_zeta(units, 256));
        assert(!class_group.try_certify_class_unit_with_zeta_bf(units, 20000,
                                                                256));
        assert(!CertificationAccess::try_certify_class_unit_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR), 256,
                false));
        assert(!class_group.try_certify_with_units(
                units, silex::CertificationMode::proven, 256, 20000));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::not_checked);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::unknown);
        assert(class_group.zeta_bf_proof_status() ==
               silex::ProofState::not_checked);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(units.certification_status() ==
               silex::CertificationMode::unknown);
    }

    // A `proven` request succeeds only through generation, proven units,
    // and relation saturation at every p | h_cand.  A Belabas-Friedman
    // audit recorded afterwards is labelled `grh`; it does not replace the
    // saturation proof.
    {
        silex::ClassGroupComputeOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        options.zeta_bf_max_cutoff = 20000;
        options.requested_certification = silex::CertificationMode::proven;
        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        assert(units.compute_with_class_group(class_group, order,
                                              sflint::FmpzConstRef(bound),
                                              options, 128));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() !=
               silex::CertificationMode::proven);

        assert(class_group.try_certify_class_unit_with_zeta_bf(units, 20000,
                                                               128));
        assert(class_group.zeta_bf_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::grh);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
    }

    return 0;
}

// x^3 - 11: h = 2 (Cl = Z/2), unit rank 1.  Oracle, GP 2.17
// `bnfinit(x^3 - 11, 1)`: `.cyc` = [2], `.reg` = 5.58720662606090776185...,
// `.tu` = [2, -1], `bnfcertify` = 1.  A `proven` request must reach
// `proven` through the saturation route, which on this field runs the ell = 2
// local proof; the Belabas-Friedman hR (GRH) supplies no proof component.
int test_saturation_proves_nontrivial_class_group_without_grh() {
    using CertificationAccess =
            silex::detail::ClassGroupCertificationAccess;

    silex::NumberField field = cubic_field(0, -11);
    silex::Order equation;
    equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    assert(order.is_maximal());
    silex::EmbeddingContext embeddings(field);

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }

    silex::ClassGroupComputeOptions proven_options;
    proven_options.max_candidates = 5000;
    proven_options.max_relations = 500;
    proven_options.zeta_bf_max_cutoff = 20000;
    proven_options.requested_certification =
            silex::CertificationMode::proven;

    silex::ClassGroupContext proven_class_group;
    silex::OrderUnitGroup proven_units;
    assert(proven_units.compute_with_class_group(
            proven_class_group, order, sflint::FmpzConstRef(bound),
            proven_options, 128));
    assert(proven_class_group.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_units.certification_status() ==
           silex::CertificationMode::proven);
    assert(proven_class_group.analytic_class_regulator_certification() !=
           silex::CertificationMode::proven);
    assert(proven_class_group.relation_saturation_status() ==
           silex::ProofState::verified);

    // Cl = Z/2, as in GP.
    assert(proven_class_group.invariant_count() == 1);
    sflint::Fmpz invariant;
    assert(proven_class_group.invariant(sflint::FmpzRef(invariant), 0));
    assert(sflint::fmpz_equal_si(invariant, 2));

    // The class group is proven by the ell = 2 local proof, and every
    // recorded saturation check is verified.
    bool saw_two = false;
    for (slong i = 0; i < proven_class_group.relation_saturation_record_count();
         ++i) {
        sflint::Fmpz ell;
        silex::ProofState status = silex::ProofState::not_checked;
        assert(proven_class_group.relation_saturation_record(
                sflint::FmpzRef(ell), status, i));
        assert(status == silex::ProofState::verified);
        if (sflint::fmpz_equal_si(ell, 2)) {
            saw_two = true;
        }
    }
    assert(saw_two);

    // R = 5.587206626060907761855..., torsion {+1, -1}, as in GP.
    sflint::Arb regulator;
    assert(proven_units.regulator(sflint::ArbRef(regulator)));
    sflint::Arb reference;
    assert(::arb_set_str(reference.raw(),
                         "5.58720662606090776185541302052 +/- 1e-28",
                         128) == 0);
    assert(::arb_overlaps(regulator.raw(), reference.raw()) != 0);
    sflint::Fmpz torsion_order;
    assert(proven_units.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));

    // Proven units, GRH-conditional hR, saturation not yet proven: every
    // analytic gate fails closed and leaves the candidate unchanged.  The
    // units are genuinely fundamental, so the analytic index-one test
    // itself passes; only the GRH rule blocks publication.
    {
        silex::ClassGroupCandidateOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), options));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        silex::OrderUnitGroup units(std::move(proven_units));
        assert(units.certification_status() ==
               silex::CertificationMode::proven);

        sflint::Arb zeta_hR;
        sflint::Arb error_bound;
        ulong cutoff = 0;
        slong work_precision = 0;
        assert(silex::zeta_class_regulator_product_bf_audit(
                sflint::ArbRef(zeta_hR), sflint::ArbRef(error_bound), cutoff,
                work_precision, order, 20000, 64));
        sflint::Fmpz aux_bound;
        assert(set_fmpz_si(aux_bound, 1000));

        const auto assert_unchanged = [&]() {
            assert(class_group.certification_status() ==
                   silex::CertificationMode::unknown);
            assert(class_group.relation_saturation_status() !=
                   silex::ProofState::verified);
            assert(class_group.analytic_class_regulator_status() ==
                   silex::ProofState::not_checked);
            assert(class_group.analytic_class_regulator_certification() ==
                   silex::CertificationMode::unknown);
            assert(class_group.zeta_bf_proof_status() ==
                   silex::ProofState::not_checked);
            assert(class_group.unit_proof_status() ==
                   silex::ProofState::not_checked);
            assert(class_group.regulator_proof_status() ==
                   silex::ProofState::not_checked);
            assert(units.certification_status() ==
                   silex::CertificationMode::proven);
        };

        assert(!class_group.try_certify_class_unit_with_zeta(units, 64));
        assert_unchanged();
        assert(!class_group.try_certify_class_unit_with_zeta_bf(units, 20000,
                                                                64));
        assert_unchanged();
        assert(!CertificationAccess::try_certify_class_unit_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR), 64,
                false));
        assert_unchanged();
        assert(!CertificationAccess::try_analytic_index_bound_with_units(
                class_group, units, sflint::ArbConstRef(zeta_hR),
                sflint::FmpzConstRef(aux_bound), 64, false));
        assert_unchanged();
        assert(!CertificationAccess::try_certify_class_unit_with_bf_audit(
                class_group, units, sflint::ArbConstRef(zeta_hR),
                sflint::ArbConstRef(error_bound), cutoff, 20000, 64,
                work_precision));
        assert_unchanged();

        // Control: the unconditional saturation route proves the same pair.
        assert(class_group.try_certify_with_units(
                units, silex::CertificationMode::proven, 256, 20000));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::unknown);

        // With saturation proven, the same GRH-conditional values pass the
        // analytic gates, which shows that the failures above came from the
        // GRH rule and not from an inconclusive enclosure.  The check is
        // recorded as `grh` and the label stays `proven`.
        assert(class_group.try_certify_class_unit_with_zeta(units, 64));
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::grh);
        assert(class_group.try_certify_class_unit_with_zeta_bf(units, 20000,
                                                               64));
        assert(class_group.zeta_bf_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::grh);
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);

        // The grh-mode acceptance record cannot claim an unconditional hR
        // for an order without an unconditional route, and it never changes
        // the certification label or the proof states.
        const silex::ProofState unit_proof_before =
                class_group.unit_proof_status();
        const silex::ProofState saturation_before =
                class_group.relation_saturation_status();
        CertificationAccess::record_grh_acceptance_analytic_check(
                class_group, true);
        assert(class_group.analytic_class_regulator_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::grh);
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.unit_proof_status() == unit_proof_before);
        assert(class_group.relation_saturation_status() ==
               saturation_before);
    }

    return 0;
}

// The exact imaginary-quadratic route proves the class group from the exact
// class number h(D) (FLINT reduced forms): with generation verified, the
// relation lattice has index h_cand / h = 1 in the full one, so it is
// saturated at every prime ell | h.  The route stores one verified
// relation-saturation record for each such ell, so a later GRH-conditional
// Belabas-Friedman audit keeps `proven` by promoting from those records.
int test_exact_imaginary_quadratic_saturation_records() {
    struct Case {
        slong radicand;
        slong class_number;
        slong prime_count;
    };
    // D = -4, -3, -20, -23, -56, -47, -87.
    const Case cases[] = {{-1, 1, 0},  {-3, 1, 0},  {-5, 2, 1},
                          {-23, 3, 1}, {-14, 4, 1}, {-47, 5, 1},
                          {-87, 6, 2}};

    for (const Case& test_case : cases) {
        silex::NumberField field = quadratic_field(test_case.radicand);
        silex::Order equation;
        equation = silex::test::equation_order(field);
        silex::Order order(field);
        assert(order.maximal_order(equation));
        assert(order.is_maximal());

        sflint::Fmpz bound;
        assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                    order));
        if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
            sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
        }

        silex::ClassGroupComputeOptions options;
        options.max_candidates = 5000;
        options.max_relations = 500;
        options.requested_certification = silex::CertificationMode::proven;

        silex::ClassGroupContext class_group;
        silex::OrderUnitGroup units;
        assert(units.compute_with_class_group(class_group, order,
                                              sflint::FmpzConstRef(bound),
                                              options, 128));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(units.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);

        sflint::Fmpz class_order;
        assert(class_group.order(sflint::FmpzRef(class_order)));
        assert(sflint::fmpz_equal_si(class_order, test_case.class_number));

        // One verified record for each distinct prime ell | h, and none
        // for h = 1.
        const auto check_records = [&]() {
            assert(class_group.relation_saturation_record_count() ==
                   test_case.prime_count);
            for (slong i = 0; i < test_case.prime_count; ++i) {
                auto record = class_group.relation_saturation_record(i);
                assert(record.has_value());
                assert(record->status == silex::ProofState::verified);
                assert(sflint::fmpz_is_prime(
                        sflint::FmpzConstRef(record->ell)));
                assert(sflint::fmpz_divisible(
                        sflint::FmpzConstRef(class_order),
                        sflint::FmpzConstRef(record->ell)));
                for (slong j = 0; j < i; ++j) {
                    auto other = class_group.relation_saturation_record(j);
                    assert(other.has_value());
                    assert(!sflint::fmpz_equal(
                            sflint::FmpzConstRef(other->ell),
                            sflint::FmpzConstRef(record->ell)));
                }
            }
        };
        check_records();

        // The BF audit records a GRH-conditional check; it succeeds by
        // promoting from the stored saturation records.
        assert(class_group.try_certify_class_unit_with_zeta_bf(units, 20000,
                                                               128));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.zeta_bf_proof_status() ==
               silex::ProofState::verified);
        assert(class_group.analytic_class_regulator_certification() ==
               silex::CertificationMode::grh);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        check_records();
    }

    return 0;
}

// Counts profile events with one label.
int test_add_dependent_unit_precision_cap() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.push(epsilon, 7));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(silex::detail::order_unit_group_set_units_internal(
            group, order, silex::FactoredElementSpan(generators, 1),
            embeddings, 128, true));

    // 2 is not a unit, so no relation with a bounded denominator exists and
    // every precision is inconclusive.  The doubling stops at the cap and
    // the group is left unchanged.
    silex::Element two(field);
    assert(two.set_si(2));
    silex::FactoredElement candidate(field);
    assert(candidate.set_element(two));
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    TorsionEventCounter counter;
    counter.label = "unit_group.dependent_relation_precision_cap";
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            torsion_event_callback, &counter);
    group.set_diagnostics(&diagnostics);
    const slong rank_before = group.free_rank();
    sflint::Arb regulator_before;
    assert(group.regulator(sflint::ArbRef(regulator_before)));

    silex::detail::RelationUnitExtractionState state;
    bool changed = true;
    assert(silex::detail::add_dependent_unit(changed, group, candidate,
                                             embeddings, state, 128));
    assert(!changed);
    assert(state.rel_add_precision == silex::detail::kRelAddStartPrecision);
    assert(group.free_rank() == rank_before);
    sflint::Arb regulator_after;
    assert(group.regulator(sflint::ArbRef(regulator_after)));
    assert(::arb_equal(regulator_after.raw(), regulator_before.raw()) != 0);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    assert(counter.count == 1);
#else
    assert(counter.count == 0);
#endif
    group.set_diagnostics(nullptr);
    return 0;
}

// The exact imaginary-quadratic index route may report index one under a
// proven request only after it has published the pair proven.  With torsion
// that Silex did not compute (Z[i] with an under-claimed w = 2), the route
// must fail closed rather than hand callers an unproven index one, and a
// failed call leaves no generation check behind.  The candidate-pair
// validation loop then recomputes the units with the computed torsion and
// publishes that pair proven; the under-claimed units are never proven.
int test_rank_zero_quadratic_index_one_requires_computed_torsion() {
    using CertificationAccess =
            silex::detail::ClassGroupCertificationAccess;
    using OrderUnitGroupAccess = silex::detail::OrderUnitGroupAccess;

    silex::NumberField field = quadratic_field(-1);
    silex::Order order;
    order = silex::test::equation_order(field);
    assert(order.is_maximal());
    silex::EmbeddingContext embeddings(field);

    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::ClassGroupCandidateOptions candidate_options;
    candidate_options.max_candidates = 256;
    candidate_options.max_relations = 48;

    sflint::Fmpz torsion_order;
    assert(set_fmpz_si(torsion_order, 2));
    silex::OrderElement torsion_generator(order);
    assert(torsion_generator.set_si(-1));
    const auto under_claimed = [&](silex::OrderUnitGroup& units) {
        assert(OrderUnitGroupAccess::set_units(
                units, order, silex::FactoredElementSpan(), embeddings, 256,
                true, &torsion_order, &torsion_generator));
        assert(!silex::detail::order_unit_torsion_is_computed(units));
    };

    // Direct: no index is returned and nothing is published.
    {
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), candidate_options));
        silex::OrderUnitGroup under(order);
        under_claimed(under);

        sflint::Fmpz index;
        assert(!CertificationAccess::rank_zero_quadratic_class_index_bound(
                sflint::FmpzRef(index), class_group, under,
                silex::CertificationMode::proven));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::not_checked);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(under.certification_status() ==
               silex::CertificationMode::unknown);
        // The failed call ran the generation check inside its transaction,
        // so that check is rolled back too.
        assert(class_group.factor_base_generation_checked_status() ==
               silex::ProofState::not_checked);
    }

    // Through the candidate-pair validation loop: the exact index route
    // fails closed and leaves no generation check, so the loop does not take
    // the integer index-one shortcut.  With the under-claimed w the analytic
    // index bound is above one, so the loop does not reach the saturation
    // fallback either (where prove_index_bound would reject the w = 2 input;
    // test_unit_proofs_require_computed_torsion covers that rejection).  It
    // recomputes the units in `scratch` with the computed torsion (w = 4),
    // proves them, and publishes that pair proven through the exact index
    // route.  The under-claimed (w = 2) units are never proven.
    {
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), candidate_options));
        silex::OrderUnitGroup under(order);
        under_claimed(under);
        silex::OrderUnitGroup scratch(order);

        silex::DiagnosticsContext diagnostics;
        silex::diagnostics_context_init(diagnostics);
        TorsionEventCounter counter;
        // The recompute runs on `scratch` and reports its proven result.
        counter.label = "unit_group.validation_recompute.result.proven";
        silex::diagnostics_set_profiling(
                diagnostics, true,
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group),
                torsion_event_callback, &counter);
        under.set_diagnostics(&diagnostics);
        scratch.set_diagnostics(&diagnostics);

        silex::ClassGroupComputeOptions options;
        options.requested_certification = silex::CertificationMode::proven;
        silex::detail::AnalyticClassRegulatorCache analytic_cache;
        silex::detail::ValidateRefineSummary summary;
        assert(silex::detail::try_validate_refine_loop(
                class_group, under, scratch, order, options, embeddings,
                analytic_cache, summary, 256, false));
        assert(summary.outcome ==
               silex::detail::ValidateRefineOutcome::proven);
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        // The loop swaps the refined units into `under`.  Whichever object
        // is proven has the computed torsion w = 4; no object with the
        // under-claimed w = 2 is proven.
        assert(under.certification_status() ==
               silex::CertificationMode::proven);
        bool saw_under_claimed = false;
        for (const silex::OrderUnitGroup* group : {&under, &scratch}) {
            if (!group->is_set()) {
                continue;
            }
            sflint::Fmpz group_torsion;
            assert(group->torsion_order(sflint::FmpzRef(group_torsion)));
            if (sflint::fmpz_equal_si(group_torsion, 2)) {
                saw_under_claimed = true;
                assert(group->certification_status() !=
                       silex::CertificationMode::proven);
            }
            if (group->certification_status() ==
                silex::CertificationMode::proven) {
                assert(sflint::fmpz_equal_si(group_torsion, 4));
                assert(silex::detail::order_unit_torsion_is_computed(
                        *group));
            }
        }
        assert(saw_under_claimed);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
        assert(counter.count >= 1);
#else
        assert(counter.count == 0);
#endif
        under.set_diagnostics(nullptr);
        scratch.set_diagnostics(nullptr);
    }

    // Control: the torsion Silex computes certifies at index one.
    {
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), candidate_options));
        silex::OrderUnitGroup computed(order);
        assert(computed.set_units(order, silex::FactoredElementSpan(),
                                  embeddings, 256));
        assert(silex::detail::order_unit_torsion_is_computed(computed));
        // Q(i) has the roots of unity {1, i, -1, -i}.
        sflint::Fmpz computed_torsion;
        assert(computed.torsion_order(sflint::FmpzRef(computed_torsion)));
        assert(sflint::fmpz_equal_si(computed_torsion, 4));

        sflint::Fmpz index;
        assert(CertificationAccess::rank_zero_quadratic_class_index_bound(
                sflint::FmpzRef(index), class_group, computed,
                silex::CertificationMode::proven));
        assert(sflint::fmpz_is_one(sflint::FmpzConstRef(index)));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.relation_saturation_status() ==
               silex::ProofState::verified);
        assert(computed.certification_status() ==
               silex::CertificationMode::proven);
    }

    return 0;
}

// The public unit proofs check the torsion of their input themselves.  On
// Z[i] with torsion under-claimed as {1, -1} (w = 2 instead of 4), the
// regulator index bound of the empty free part is one, so without the check
// prove_index_bound and the saturate methods would label the group proven,
// prove_local_saturated would record an ell = 2 result, and the grh
// publication would label the pair grh.  They must fail closed and leave the
// output unchanged; the same group with computed torsion is accepted.
int test_unit_proofs_require_computed_torsion() {
    using OrderUnitGroupAccess = silex::detail::OrderUnitGroupAccess;

    silex::NumberField field = quadratic_field(-1);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::OrderUnitGroup computed(order);
    assert(computed.set_units(order, silex::FactoredElementSpan(), embeddings,
                              256));
    assert(silex::detail::order_unit_torsion_is_computed(computed));

    sflint::Fmpz torsion_order;
    assert(set_fmpz_si(torsion_order, 2));
    silex::OrderElement torsion_generator(order);
    assert(torsion_generator.set_si(-1));
    silex::OrderUnitGroup under(order);
    assert(OrderUnitGroupAccess::set_units(
            under, order, silex::FactoredElementSpan(), embeddings, 256, true,
            &torsion_order, &torsion_generator));
    assert(!silex::detail::order_unit_torsion_is_computed(under));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 64));
    sflint::Fmpz aux_max;
    assert(set_fmpz_si(aux_max, 256));
    const auto torsion_is = [](const silex::OrderUnitGroup& units, slong w) {
        sflint::Fmpz stored;
        return units.torsion_order(sflint::FmpzRef(stored)) &&
               sflint::fmpz_equal_si(stored, w);
    };

    // Each rejection below emits its own torsion event; under profiling
    // builds the count shows which check fired.
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    TorsionEventCounter counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            torsion_event_callback, &counter);
    const auto expect_events = [&](slong expected) {
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
        assert(counter.count == expected);
#else
        (void)expected;
        assert(counter.count == 0);
#endif
    };

    // prove_index_bound.
    {
        counter.label = "unit_group.prove_index_bound.torsion_not_computed";
        counter.count = 0;
        silex::OrderUnitGroup out(order);
        out.set_diagnostics(&diagnostics);
        silex::ProofState status = silex::ProofState::not_checked;
        bool changed = true;
        assert(!out.prove_index_bound(status, changed, under, 4,
                                      sflint::FmpzConstRef(aux_bound), 2,
                                      embeddings, 256));
        assert(!out.is_set());
        assert(status == silex::ProofState::not_checked);
        assert(changed);
        expect_events(1);

        assert(out.prove_index_bound(status, changed, computed, 4,
                                     sflint::FmpzConstRef(aux_bound), 2,
                                     embeddings, 256));
        assert(status == silex::ProofState::verified);
        assert(out.certification_status() == silex::CertificationMode::proven);
        assert(torsion_is(out, 4));

        // A rejected input leaves a previously published result unchanged.
        assert(!out.prove_index_bound(status, changed, under, 4,
                                      sflint::FmpzConstRef(aux_bound), 2,
                                      embeddings, 256));
        assert(out.certification_status() == silex::CertificationMode::proven);
        assert(torsion_is(out, 4));
        expect_events(2);
    }

    // saturate_index_bounded.
    {
        counter.label = "unit_group.saturate_index_bounded.torsion_not_computed";
        counter.count = 0;
        silex::OrderUnitGroup out(order);
        out.set_diagnostics(&diagnostics);
        bool changed = true;
        bool stable = true;
        assert(!out.saturate_index_bounded(changed, stable, under, embeddings,
                                           4, sflint::FmpzConstRef(aux_bound),
                                           2, 256));
        assert(!out.is_set());
        assert(changed && stable);
        expect_events(1);

        // An unset input fails as an argument error, not a torsion event.
        silex::OrderUnitGroup unset(order);
        assert(!out.saturate_index_bounded(changed, stable, unset, embeddings,
                                           4, sflint::FmpzConstRef(aux_bound),
                                           2, 256));
        expect_events(1);

        assert(out.saturate_index_bounded(changed, stable, computed,
                                          embeddings, 4,
                                          sflint::FmpzConstRef(aux_bound), 2,
                                          256));
        assert(out.certification_status() == silex::CertificationMode::proven);
        assert(torsion_is(out, 4));

        // A rejected input leaves a previously published result unchanged.
        out.set_diagnostics(&diagnostics);
        changed = true;
        stable = true;
        assert(!out.saturate_index_bounded(changed, stable, under, embeddings,
                                           4, sflint::FmpzConstRef(aux_bound),
                                           2, 256));
        assert(out.certification_status() == silex::CertificationMode::proven);
        assert(torsion_is(out, 4));
        assert(changed && stable);
        expect_events(2);
    }

    // saturate_index_bounded_adaptive, on its own: it runs the passes
    // through the unchecked body, so its own entry check is what rejects.
    {
        counter.label =
                "unit_group.saturate_index_bounded_adaptive.torsion_not_computed";
        counter.count = 0;
        silex::OrderUnitGroup adaptive(order);
        adaptive.set_diagnostics(&diagnostics);
        bool changed = true;
        bool stable = true;
        assert(!adaptive.saturate_index_bounded_adaptive(
                changed, stable, under, embeddings, 4,
                sflint::FmpzConstRef(aux_bound),
                sflint::FmpzConstRef(aux_max), 2, 256));
        assert(!adaptive.is_set());
        assert(changed && stable);
        expect_events(1);

        assert(adaptive.saturate_index_bounded_adaptive(
                changed, stable, computed, embeddings, 4,
                sflint::FmpzConstRef(aux_bound),
                sflint::FmpzConstRef(aux_max), 2, 256));
        assert(adaptive.certification_status() ==
               silex::CertificationMode::proven);
        assert(torsion_is(adaptive, 4));

        adaptive.set_diagnostics(&diagnostics);
        changed = true;
        stable = true;
        assert(!adaptive.saturate_index_bounded_adaptive(
                changed, stable, under, embeddings, 4,
                sflint::FmpzConstRef(aux_bound),
                sflint::FmpzConstRef(aux_max), 2, 256));
        assert(adaptive.certification_status() ==
               silex::CertificationMode::proven);
        assert(torsion_is(adaptive, 4));
        assert(changed && stable);
        expect_events(2);
    }

    // prove_local_saturated at ell = 2, which divides both w = 2 and w = 4,
    // so the ell-local test uses the torsion column.  It records no proof
    // on the under-claimed torsion and leaves the output unchanged.
    {
        counter.label = "unit_group.prove_local_saturated.torsion_not_computed";
        counter.count = 0;
        sflint::Fmpz ell;
        assert(set_fmpz_si(ell, 2));
        silex::OrderUnitGroup local(order);
        local.set_diagnostics(&diagnostics);
        silex::ProofState status = silex::ProofState::not_checked;
        bool changed = true;
        assert(!local.prove_local_saturated(
                status, changed, under, sflint::FmpzConstRef(ell), 4,
                sflint::FmpzConstRef(aux_bound), embeddings, 256));
        assert(!local.is_set());
        assert(status == silex::ProofState::not_checked);
        assert(changed);
        expect_events(1);

        assert(local.prove_local_saturated(
                status, changed, computed, sflint::FmpzConstRef(ell), 4,
                sflint::FmpzConstRef(aux_bound), embeddings, 256));
        assert(status == silex::ProofState::verified);
        assert(local.unit_proof_verified(sflint::FmpzConstRef(ell)));
        assert(torsion_is(local, 4));
        const slong records = local.unit_proof_record_count();

        local.set_diagnostics(&diagnostics);
        status = silex::ProofState::not_checked;
        changed = true;
        assert(!local.prove_local_saturated(
                status, changed, under, sflint::FmpzConstRef(ell), 4,
                sflint::FmpzConstRef(aux_bound), embeddings, 256));
        assert(status == silex::ProofState::not_checked);
        assert(changed);
        assert(local.unit_proof_verified(sflint::FmpzConstRef(ell)));
        assert(local.unit_proof_record_count() == records);
        assert(torsion_is(local, 4));
        expect_events(2);
    }

    // try_certify_with_units treats a proven label on uncomputed torsion as
    // unproven: a proven request fails and publishes nothing.  The same
    // label on computed torsion certifies (h = 1).
    {
        sflint::Fmpz bound;
        assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                    order));
        if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
            sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
        }
        silex::ClassGroupCandidateOptions candidate_options;
        candidate_options.max_candidates = 256;
        candidate_options.max_relations = 48;

        silex::OrderUnitGroup labelled(order);
        assert(labelled.set(under));
        silex::detail::TestHookAccess::mark_certification_proven(labelled);
        assert(labelled.certification_status() ==
               silex::CertificationMode::proven);

        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), candidate_options));
        assert(!class_group.try_certify_with_units(
                labelled, silex::CertificationMode::proven, 256));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(class_group.unit_proof_status() ==
               silex::ProofState::not_checked);
        assert(class_group.regulator_proof_status() ==
               silex::ProofState::not_checked);

        silex::OrderUnitGroup good(order);
        assert(good.set(computed));
        silex::detail::TestHookAccess::mark_certification_proven(good);
        assert(class_group.try_certify_with_units(
                good, silex::CertificationMode::proven, 256));
        assert(class_group.certification_status() ==
               silex::CertificationMode::proven);
        assert(class_group.unit_proof_status() == silex::ProofState::verified);
    }

    // The grh publication of a completed class/unit pair checks the torsion
    // itself: GRH does not bear on w, so under-claimed torsion gets no grh
    // label and neither object changes.  Computed torsion is published.
    {
        using CertificationAccess =
                silex::detail::ClassGroupCertificationAccess;
        counter.label = "unit_group.grh_publication.torsion_not_computed";
        counter.count = 0;
        sflint::Fmpz bound;
        assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                    order));
        if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
            sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
        }
        silex::ClassGroupCandidateOptions candidate_options;
        candidate_options.max_candidates = 256;
        candidate_options.max_relations = 48;
        silex::ClassGroupContext class_group;
        assert(class_group.compute_candidate(
                order, sflint::FmpzConstRef(bound), candidate_options));

        silex::OrderUnitGroup labelled(order);
        assert(labelled.set(under));
        assert(!CertificationAccess::publish_grh_labels(class_group, labelled,
                                                        &diagnostics));
        assert(class_group.certification_status() ==
               silex::CertificationMode::unknown);
        assert(labelled.certification_status() ==
               silex::CertificationMode::unknown);
        assert(torsion_is(labelled, 2));
        expect_events(1);

        silex::OrderUnitGroup good(order);
        assert(good.set(computed));
        assert(CertificationAccess::publish_grh_labels(class_group, good,
                                                       &diagnostics));
        assert(class_group.certification_status() ==
               silex::CertificationMode::grh);
        assert(good.certification_status() == silex::CertificationMode::grh);
        expect_events(1);
    }

    return 0;
}

// Positive-rank companion of test_unit_proofs_require_computed_torsion.  In
// Q(zeta8) = Q(x), x^4 + 1 = 0 (rank 1, w = 8) the units are
// <zeta8> x <1 + sqrt2> with sqrt2 = x - x^3.  A torsion under-claim of
// {+-1} or {+-1, +-i} paired with 1 + sqrt2 generates a subgroup of index
// 4 or 2: the ell = 2 torsion column then describes the wrong torsion and
// hides zeta8.  prove_index_bound must reject both under-claims and prove
// the group with the computed torsion, w = 8.
int test_unit_proofs_require_computed_torsion_positive_rank() {
    using OrderUnitGroupAccess = silex::detail::OrderUnitGroupAccess;

    silex::NumberField field = quartic_field(0, 0, 0, 1);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element x(field);
    assert(x.gen());
    silex::Element x2(field);
    assert(x2.multiply(x, x));
    silex::Element x3(field);
    assert(x3.multiply(x2, x));
    silex::Element unit(field);
    assert(unit.subtract(x, x3));
    assert(unit.add_si(unit, 1));
    silex::FactoredElement generator(field);
    assert(generator.set_element(unit));
    silex::FactoredElement generators[] = {std::move(generator)};
    const silex::FactoredElementSpan span(generators, 1);

    const auto torsion_is = [](const silex::OrderUnitGroup& units, slong w) {
        sflint::Fmpz stored;
        return units.torsion_order(sflint::FmpzRef(stored)) &&
               sflint::fmpz_equal_si(stored, w);
    };

    silex::OrderUnitGroup computed(order);
    assert(computed.set_units(order, span, embeddings, 256));
    assert(silex::detail::order_unit_torsion_is_computed(computed));
    assert(computed.free_rank() == 1);
    assert(torsion_is(computed, 8));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 64));

    // Under-claims: w = 2 with generator -1, w = 4 with generator x^2 = i.
    for (slong w : {2, 4}) {
        sflint::Fmpz torsion_order;
        assert(set_fmpz_si(torsion_order, w));
        silex::OrderElement torsion_generator(order);
        if (w == 2) {
            assert(torsion_generator.set_si(-1));
        } else {
            assert(torsion_generator.set_element(x2));
        }
        silex::OrderUnitGroup under(order);
        assert(OrderUnitGroupAccess::set_units(
                under, order, span, embeddings, 256, true, &torsion_order,
                &torsion_generator));
        assert(!silex::detail::order_unit_torsion_is_computed(under));
        assert(under.free_rank() == 1);
        assert(torsion_is(under, w));

        silex::OrderUnitGroup out(order);
        silex::ProofState status = silex::ProofState::not_checked;
        bool changed = true;
        assert(!out.prove_index_bound(status, changed, under, 4,
                                      sflint::FmpzConstRef(aux_bound), 2,
                                      embeddings, 256));
        assert(!out.is_set());
        assert(status == silex::ProofState::not_checked);
        assert(changed);
    }

    silex::OrderUnitGroup out(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(out.prove_index_bound(status, changed, computed, 4,
                                 sflint::FmpzConstRef(aux_bound), 2,
                                 embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(out.certification_status() == silex::CertificationMode::proven);
    assert(out.free_rank() == 1);
    assert(torsion_is(out, 8));

    return 0;
}

// The cached-torsion test rejects any w with phi(w) not dividing the degree
// before it computes a power, and otherwise requires a generator of exact
// order w in the order.
int test_cached_torsion_is_valid() {
    silex::NumberField field = quadratic_field(-1);
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderElement minus_one(order);
    assert(minus_one.set_si(-1));
    silex::OrderElement i_unit(order);
    sflint::Fmpz w;
    {
        silex::OrderUnitGroup units(order);
        assert(units.compute(order));
        assert(units.torsion_generator(i_unit));
        assert(units.torsion_order(sflint::FmpzRef(w)));
        assert(sflint::fmpz_equal_si(w, 4));
    }

    const auto valid = [&](slong value, const silex::OrderElement& g) {
        sflint::Fmpz order_value;
        assert(set_fmpz_si(order_value, value));
        return silex::detail::cached_torsion_is_valid(
                order, sflint::FmpzConstRef(order_value), g);
    };
    assert(valid(4, i_unit));
    assert(valid(2, minus_one));  // a subgroup: valid, not all torsion
    assert(!valid(2, i_unit));    // wrong exact order
    assert(!valid(8, i_unit));    // phi(8) = 4 does not divide 2
    assert(!valid(6, minus_one)); // phi(6) = 2 divides 2; not order 6
    assert(!valid(3, minus_one)); // odd
    assert(!valid(0, minus_one));
    assert(!valid(-2, minus_one));
    assert(!valid(WORD(1) << 62, minus_one));  // phi = 2^61
    return 0;
}

int test_prove_index_bound() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 5));
    silex::OrderUnitGroup proved_rank_zero(degree_one_order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved_rank_zero.prove_index_bound(
            status, changed, rank_zero, 1, sflint::FmpzConstRef(aux_bound), 1,
            degree_one_embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved_rank_zero.certification_status() ==
           silex::CertificationMode::proven);
    assert(proved_rank_zero.unit_proof_record_count() == 0);

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup full_rank(order);
    assert(full_rank.set_units(order,
                               silex::FactoredElementSpan(
                                       primitive_generators, 1),
                               embeddings, 128));

    assert(set_fmpz_si(aux_bound, 31));
    silex::OrderUnitGroup proved(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_index_bound(
            status, changed, full_rank, 1, sflint::FmpzConstRef(aux_bound), 1,
            embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    assert(check_first_free_generator(proved, epsilon));
    // log(1 + sqrt2) = 0.8813... and the Friedman 1989 Table 6 bound
    // R >= 0.48 for signature (2, 0) give the index bound
    // floor(0.8813... / 0.48) = 1, so no ell-local proof is needed.
    assert(proved.unit_proof_record_count() == 0);

    silex::Element epsilon2(field);
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement square_generators[] = {std::move(square)};
    silex::OrderUnitGroup square_group(order);
    assert(square_group.set_units(order,
                                  silex::FactoredElementSpan(square_generators,
                                                             1),
                                  embeddings, 128));

    silex::OrderUnitGroup restarted(order);
    status = silex::ProofState::not_checked;
    changed = false;
    assert(restarted.prove_index_bound(
            status, changed, square_group, 1, sflint::FmpzConstRef(aux_bound),
            1, embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(changed);
    assert(restarted.certification_status() ==
           silex::CertificationMode::proven);
    assert(check_first_free_generator_square(restarted, epsilon2));
    // The square subgroup has index bound 3; the 2-local step adjoins the
    // square root and restarts, and the restarted full group has index
    // bound 1, so the only record is the 2-local step that changed it.
    assert(restarted.unit_proof_record_count() == 1);
    assert(has_unit_proof_record(restarted, 2, silex::ProofState::unavailable,
                                 true));

    silex::OrderUnitGroup restart_limited(order);
    status = silex::ProofState::not_checked;
    changed = false;
    assert(restart_limited.prove_index_bound(
            status, changed, square_group, 1, sflint::FmpzConstRef(aux_bound),
            0, embeddings, 128));
    assert(status == silex::ProofState::unavailable);
    assert(changed);
    assert(restart_limited.certification_status() ==
           silex::CertificationMode::unknown);
    assert(check_first_free_generator_square(restart_limited, epsilon2));
    assert(has_unit_proof_record(restart_limited, 2,
                                 silex::ProofState::unavailable, true));

    sflint::Fmpz bad_bound;
    assert(set_fmpz_si(bad_bound, 1));
    status = silex::ProofState::verified;
    changed = true;
    assert(!proved.prove_index_bound(
            status, changed, full_rank, 1, sflint::FmpzConstRef(bad_bound), 1,
            embeddings, 128));
    assert(status == silex::ProofState::verified);
    assert(changed);
    assert(check_first_free_generator(proved, epsilon));

    return 0;
}

// True when the first free generator of `group` is one of +/-unit^(+/-1).
// True when the regulator ball contains `expected`, a 100-digit decimal
// read at 1024 bits (its conversion error is far below the ball radius).
bool regulator_contains(const silex::OrderUnitGroup& group,
                        const char* expected) noexcept {
    sflint::Arb regulator;
    sflint::Arb reference;
    return group.regulator(sflint::ArbRef(regulator)) &&
           ::arb_set_str(reference.raw(), expected, 1024) == 0 &&
           ::arb_contains(regulator.raw(), reference.raw()) != 0;
}

bool set_single_unit_group(silex::OrderUnitGroup& out,
                           const silex::Order& order,
                           const silex::Element& unit,
                           silex::EmbeddingContext& embeddings) noexcept {
    const silex::NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }
    silex::FactoredElement generator(*field);
    if (!generator.set_element(unit)) {
        return false;
    }
    silex::FactoredElement generators[] = {std::move(generator)};
    return out.set_units(order, silex::FactoredElementSpan(generators, 1),
                         embeddings, 256);
}

// Proves `start` in the non-maximal `order` and checks that the published
// group never leaves the order.  The residue-character proof cannot see an
// l-th root that lies in the maximal order but not in `order`, so such a
// kernel vector must fail closed rather than adjoin the root.
void assert_nonmaximal_proof_fails_closed(
        const silex::Order& order,
        const silex::OrderUnitGroup& start,
        const silex::Element& expected_generator,
        const char* expected_regulator,
        slong aux_bound_value,
        silex::EmbeddingContext& embeddings) {
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, aux_bound_value));
    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = false;
    assert(proved.prove_index_bound(status, changed, start, 8,
                                    sflint::FmpzConstRef(aux_bound), 4,
                                    embeddings, 256));
    assert(status == silex::ProofState::unavailable);
    assert(proved.certification_status() ==
           silex::CertificationMode::unknown);
    assert(proved.free_rank() == 1);
    assert(check_first_free_generator_is_order_unit(proved, order));
    assert(first_free_generator_is_pm_unit_power_pm1(proved,
                                                     expected_generator));
    assert(regulator_contains(proved, expected_regulator));
    sflint::Fmpz torsion_order;
    assert(proved.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
}

// Z[sqrt5] is the non-maximal equation order of Q(sqrt5), of index 2 in
// Z[phi], phi = (1 + sqrt5)/2.  Its unit group is <-1, 2 + sqrt5> with
// 2 + sqrt5 = phi^3 and regulator 3 log(phi); reference GP 2.17.4 session:
//     quadunit(20)       \\ 2 + w, w = sqrt5
//     quadregulator(20)  \\ 1.4436354751788103424...
//     quadregulator(5)   \\ 0.4812118250596034474... (Z[phi])
// A proof must never adjoin phi, which is an l-th root in Q(sqrt5) of the
// supplied units but lies outside Z[sqrt5].
int test_prove_index_bound_nonmaximal_quadratic_sqrt5() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    assert(!equation.is_maximal());
    silex::EmbeddingContext embeddings(field);
    // GP: \p 100; quadregulator(20)
    const char* const kLogFundamental =
            "1.443635475178810342493276740273105269405553003156981558983054"
            "506520491602824665323236028287368170425";

    silex::Element order_unit(field);  // 2 + sqrt5 = phi^3
    assert(set_quadratic_coeffs(order_unit, 2, 1, 1, 1));
    silex::Element phi6(field);  // 9 + 4 sqrt5 = (2 + sqrt5)^2
    assert(phi6.multiply(order_unit, order_unit));
    silex::Element phi12(field);  // 161 + 72 sqrt5 = (2 + sqrt5)^4
    assert(phi12.multiply(phi6, phi6));

    // With aux bound 200 every prime l uses the direct proof route: l = 2
    // adjoins the square roots 9 + 4 sqrt5 and 2 + sqrt5, which lie in the
    // order, and l = 3 then meets phi and fails closed.
    const silex::Element* inputs[] = {&order_unit, &phi6, &phi12};
    for (const silex::Element* input : inputs) {
        silex::OrderUnitGroup start(equation);
        assert(set_single_unit_group(start, equation, *input, embeddings));
        assert_nonmaximal_proof_fails_closed(equation, start, order_unit,
                                             kLogFundamental, 200,
                                             embeddings);
    }
    // Aux bounds 3 and 2 send l = 3 through the in-place stable route.
    silex::OrderUnitGroup fundamental_start(equation);
    assert(set_single_unit_group(fundamental_start, equation, order_unit,
                                 embeddings));
    assert_nonmaximal_proof_fails_closed(equation, fundamental_start,
                                         order_unit, kLogFundamental, 3,
                                         embeddings);
    assert_nonmaximal_proof_fails_closed(equation, fundamental_start,
                                         order_unit, kLogFundamental, 2,
                                         embeddings);

    // The public row saturation leaves the group unchanged when the l-th
    // root (phi for l = 3) lies outside the order.
    silex::OrderUnitGroup start(equation);
    assert(set_single_unit_group(start, equation, order_unit, embeddings));
    sflint::FmpzMat row(1, 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(row, 0, 0), 1);
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 3));
    silex::OrderUnitGroup saturated(equation);
    bool changed = true;
    assert(saturated.saturate_row(changed, start,
                                  sflint::FmpzMatConstRef(row), 0,
                                  sflint::FmpzConstRef(ell), embeddings, 256));
    assert(!changed);
    assert(check_first_free_generator(saturated, order_unit));
    assert(saturated.certification_status() ==
           silex::CertificationMode::unknown);

    // Row saturation still adjoins a square root that lies in the order.
    silex::OrderUnitGroup square_start(equation);
    assert(set_single_unit_group(square_start, equation, phi6, embeddings));
    assert(set_fmpz_si(ell, 2));
    // ell = 2 divides w = 2, so the row carries a torsion exponent (0).
    sflint::FmpzMat square_row(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(square_row, 0, 0), 1);
    silex::OrderUnitGroup square_saturated(equation);
    changed = false;
    assert(square_saturated.saturate_row(
            changed, square_start, sflint::FmpzMatConstRef(square_row), 0,
            sflint::FmpzConstRef(ell), embeddings, 256));
    assert(changed);
    assert(check_first_free_generator_is_order_unit(square_saturated,
                                                    equation));
    assert(first_free_generator_is_pm_unit_power_pm1(square_saturated,
                                                     order_unit));
    return 0;
}

// Detail-level check of the dependent-relation adjoin.  With the group
// <2 + sqrt5> of Z[sqrt5] and the relation row (-1, 3), i.e.
// root^3 = 2 + sqrt5, the root phi lies outside the order and must not be
// adjoined; the in-order square root 2 + sqrt5 of <9 + 4 sqrt5> still is.
int test_adjoin_dependent_relation_nonmaximal_root() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    assert(!equation.is_maximal());
    silex::EmbeddingContext embeddings(field);

    silex::Element phi(field);
    assert(set_quadratic_coeffs(phi, 1, 2, 1, 2));
    silex::Element order_unit(field);  // 2 + sqrt5 = phi^3
    assert(set_quadratic_coeffs(order_unit, 2, 1, 1, 1));
    silex::Element order_unit2(field);
    assert(order_unit2.multiply(order_unit, order_unit));

    silex::OrderUnitGroup group(equation);
    assert(set_single_unit_group(group, equation, order_unit, embeddings));
    silex::FactoredElement root(field);
    assert(root.set_element(phi));
    sflint::FmpzMat rel(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rel, 0, 0), -1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rel, 0, 1), 3);

    silex::OrderUnitGroup out(equation);
    bool changed = true;
    assert(silex::detail::adjoin_verified_dependent_relation(
            changed, out, group, root, sflint::FmpzMatConstRef(rel), 0,
            embeddings, 256));
    assert(!changed);
    assert(check_first_free_generator(out, order_unit));

    silex::OrderUnitGroup square_group(equation);
    assert(set_single_unit_group(square_group, equation, order_unit2,
                                 embeddings));
    silex::FactoredElement square_root(field);
    assert(square_root.set_element(order_unit));
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rel, 0, 1), 2);
    silex::OrderUnitGroup square_out(equation);
    changed = false;
    assert(silex::detail::adjoin_verified_dependent_relation(
            changed, square_out, square_group, square_root,
            sflint::FmpzMatConstRef(rel), 0, embeddings, 256));
    assert(changed);
    assert(check_first_free_generator_is_order_unit(square_out, equation));
    assert(first_free_generator_is_pm_unit_power_pm1(square_out, order_unit));
    return 0;
}

// Control: in the maximal order Z[phi] the same proof adjoins phi and
// proves the group with regulator log(phi).
int test_prove_index_bound_maximal_quadratic_sqrt5_control() {
    silex::NumberField field = quadratic_field(5);
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal(field);
    assert(maximal.maximal_order(equation));
    assert(maximal.is_maximal());
    silex::EmbeddingContext embeddings(field);

    silex::Element phi(field);
    assert(set_quadratic_coeffs(phi, 1, 2, 1, 2));
    silex::Element phi3(field);
    silex::Element phi6(field);
    assert(phi3.multiply(phi, phi));
    assert(phi3.multiply(phi3, phi));
    assert(phi6.multiply(phi3, phi3));

    silex::OrderUnitGroup start(maximal);
    assert(set_single_unit_group(start, maximal, phi6, embeddings));
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 200));
    silex::OrderUnitGroup proved(maximal);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = false;
    assert(proved.prove_index_bound(status, changed, start, 8,
                                    sflint::FmpzConstRef(aux_bound), 4,
                                    embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    assert(first_free_generator_is_pm_unit_power_pm1(proved, phi));
    // GP: \p 100; quadregulator(5)
    assert(regulator_contains(
            proved,
            "0.481211825059603447497758913424368423135184334385660519661018"
            "1688401638676082217744120094291227234750"));
    return 0;
}

// x^3 - 4x^2 - 8x - 1 (discriminant 2213, totally real, Z[theta] maximal):
// GP 2.17 `K = bnfinit(x^3-4*x^2-8*x-1, 1)` gives the fundamental units
// theta and -4 theta^2 - theta + 6, class number one, and bnfcertify(K) = 1.
// The regulator is 12.6808...  On the maximal order the Pohst-Zassenhaus
// (6.22) bound R >= 2.4930... gives the regulator index bound
// floor(12.6808... / 2.4930...) = 5.  The proofs below run on the equation
// order, the same lattice but not known to be maximal, where the bound is
// the Friedman 1989 Table 6 bound R >= 0.52 for signature (3, 0) and the
// index bound is floor(12.6808... / 0.52) = 24.  With the pre-scan bound 200,
// ell = 13, 17
// and 19 each have a single degree-one prime q = 1 mod ell below the bound,
// so the unit-rank-two local proof must continue past the pre-scan
// (reference `saturate!` scan) and still verify every ell <= 24.
int test_prove_index_bound_cubic2213_past_prescan_bound() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -8);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order maximal = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::Order order = silex::test::equation_order(field);
    assert(!order.maximality_known());
    silex::EmbeddingContext embeddings(field);

    sflint::FmpqPoly theta_poly;
    sflint::fmpq_poly_set_coeff_si(theta_poly, 1, 1);
    sflint::FmpqPoly second_poly;
    sflint::fmpq_poly_set_coeff_si(second_poly, 2, -4);
    sflint::fmpq_poly_set_coeff_si(second_poly, 1, -1);
    sflint::fmpq_poly_set_coeff_si(second_poly, 0, 6);
    silex::Element theta(field);
    silex::Element second(field);
    assert(theta.set_fmpq_poly(sflint::FmpqPolyConstRef(theta_poly)));
    assert(second.set_fmpq_poly(sflint::FmpqPolyConstRef(second_poly)));
    silex::FactoredElement generators[] = {silex::FactoredElement(field),
                                           silex::FactoredElement(field)};
    assert(generators[0].set_element(theta));
    assert(generators[1].set_element(second));
    silex::OrderUnitGroup start(order);
    assert(start.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 256));
    assert(start.free_rank() == 2);

    silex::OrderUnitGroup maximal_start(maximal);
    assert(maximal_start.set_units(
            maximal, silex::FactoredElementSpan(generators, 2), embeddings,
            256));
    sflint::Fmpz maximal_index_bound;
    assert(maximal_start.regulator_index_bound(
            sflint::FmpzRef(maximal_index_bound), 256));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(maximal_index_bound), 5));

    sflint::Fmpz index_bound;
    assert(start.regulator_index_bound(sflint::FmpzRef(index_bound), 256));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(index_bound), 24));

    // The validation proof uses this pre-scan bound (kComputeProofAuxMax).
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 1000));

    // GP, f = x^3 - 4x^2 - 8x - 1:
    //     forprime(q=2, 1000, if(q%31==1 && #polrootsmod(f,q), print(q)))
    // prints only 683 (one root), so one character cannot settle rank two.
    // The public bounded local proof stays bounded by aux_bound and is
    // unavailable.
    sflint::Fmpz ell31;
    assert(set_fmpz_si(ell31, 31));
    silex::OrderUnitGroup local(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(local.prove_local_saturated(
            status, changed, start, sflint::FmpzConstRef(ell31), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 256));
    assert(status == silex::ProofState::unavailable);
    assert(!changed);
    assert(!local.unit_proof_verified(sflint::FmpzConstRef(ell31)));

    silex::OrderUnitGroup validated(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(validated.prove_index_bound(status, changed, start, 1,
                                       sflint::FmpzConstRef(aux_bound), 2,
                                       embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(validated.certification_status() ==
           silex::CertificationMode::proven);
    assert(has_unit_proof_record(validated, 23, silex::ProofState::verified,
                                 false));
    assert(!has_unit_proof_record(validated, 29, silex::ProofState::verified,
                                  false));

    // GP: forprime(ell=13, 19, forprime(q=2, 200,
    //         if(q%ell==1 && #polrootsmod(f,q), print(ell, " ", q))))
    // prints 13 157, 17 137 and 19 191, each with one root: the pre-scan
    // below 200 yields one character for each of these ell.
    sflint::Fmpz small_aux_bound;
    assert(set_fmpz_si(small_aux_bound, 200));
    silex::OrderUnitGroup proved(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(proved.prove_index_bound(status, changed, start, 1,
                                    sflint::FmpzConstRef(small_aux_bound), 2,
                                    embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    const slong kProvedElls[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    for (slong ell : kProvedElls) {
        assert(has_unit_proof_record(proved, ell, silex::ProofState::verified,
                                     false));
    }
    assert(!has_unit_proof_record(proved, 29, silex::ProofState::verified,
                                  false));
    // GP: \p 110; K.reg
    assert(regulator_contains(
            proved,
            "12.68082022713452093603368239737775544244960210761943928030154"
            "4077156925247307418251826663747356303483"));
    return 0;
}

// Q(zeta_13), f = x^12 + x^11 + ... + 1 (Z[zeta] maximal).  The cyclotomic
// units u_a = 1 + zeta + ... + zeta^(a-1), a = 2..6, are fundamental: GP 2.17
// `K = bnfinit(polcyclo(13), 1)` gives class number one, bnfcertify(K) = 1,
// and the log-embedding determinant of u_2..u_6 equals K.reg.  At ell = 17 the
// stable scan stops with a nonempty character kernel whose row roots stay
// undecided (degree 12), so the proof follows the reference `saturate!`
// (Saturate.jl:396-411): it doubles `stable` and rescans until the kernel is
// empty, instead of failing at once.
int test_prove_index_bound_cyclotomic13_undecided_root() {
    sflint::FmpqPoly polynomial;
    for (slong i = 0; i <= 12; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, 1);
    }
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);

    std::vector<silex::FactoredElement> generators;
    for (slong a = 2; a <= 6; ++a) {
        sflint::FmpqPoly unit_poly;
        for (slong i = 0; i < a; ++i) {
            sflint::fmpq_poly_set_coeff_si(unit_poly, i, 1);
        }
        silex::Element unit(field);
        assert(unit.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly)));
        generators.emplace_back(field);
        assert(generators.back().set_element(unit));
    }
    silex::OrderUnitGroup start(order);
    assert(start.set_units(order,
                           silex::FactoredElementSpan(generators.data(),
                                                      generators.size()),
                           embeddings, 256));
    assert(start.free_rank() == 5);

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 1000));
    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_index_bound(status, changed, start, 1,
                                    sflint::FmpzConstRef(aux_bound), 2,
                                    embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    assert(has_unit_proof_record(proved, 17, silex::ProofState::verified,
                                 false));
    // GP: \p 110; K.reg
    assert(regulator_contains(
            proved,
            "120.7840313626563127293175827399128565840801560246896328383867"
            "406585225412356786592960914162040943733"));
    return 0;
}

// Resets the stable proof scan limits when a test leaves its scope.
struct StableProofScanLimitsForTesting {
    StableProofScanLimitsForTesting(slong max_aux_primes,
                                    ulong max_q) noexcept {
        silex::detail::set_stable_proof_scan_limits_for_testing(
                max_aux_primes, max_q);
    }
    ~StableProofScanLimitsForTesting() {
        silex::detail::set_stable_proof_scan_limits_for_testing(0, 0);
    }
    StableProofScanLimitsForTesting(const StableProofScanLimitsForTesting&) =
            delete;
    StableProofScanLimitsForTesting& operator=(
            const StableProofScanLimitsForTesting&) = delete;
};

// Q(zeta_13) with the cyclotomic units u_2^17, u_3, ..., u_6 (see
// test_prove_index_bound_cyclotomic13_undecided_root): a subgroup of index 17
// in the full unit group, so it is not 17-saturated, and the 17th root of
// the first generator is u_2, which Silex cannot compute in degree 12.
silex::NumberField cyclotomic13_field() noexcept {
    sflint::FmpqPoly polynomial;
    for (slong i = 0; i <= 12; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, 1);
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool set_cyclotomic13_u2_power17_group(
        silex::OrderUnitGroup& out,
        const silex::Order& order,
        const silex::NumberField& field,
        silex::EmbeddingContext& embeddings) noexcept {
    std::vector<silex::FactoredElement> generators;
    for (slong a = 2; a <= 6; ++a) {
        sflint::FmpqPoly unit_poly;
        for (slong i = 0; i < a; ++i) {
            sflint::fmpq_poly_set_coeff_si(unit_poly, i, 1);
        }
        silex::Element unit(field);
        if (!unit.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly))) {
            return false;
        }
        if (a == 2) {
            // The expanded element u_2^17, so no factored exponent
            // reveals the root.
            sflint::Fmpz exponent;
            silex::Element power(field);
            if (!set_fmpz_si(exponent, 17) ||
                !power.pow_fmpz(unit, sflint::FmpzConstRef(exponent))) {
                return false;
            }
            unit = std::move(power);
        }
        generators.emplace_back(field);
        if (!generators.back().set_element(unit)) {
            return false;
        }
    }
    return out.set_units(order,
                         silex::FactoredElementSpan(generators.data(),
                                                    generators.size()),
                         embeddings, 256) &&
           out.free_rank() == 5;
}

// u_2 = 1 + zeta_13 in Q(zeta_13) (see cyclotomic13_field).
bool set_cyclotomic13_u2(silex::Element& out) noexcept {
    sflint::FmpqPoly unit_poly;
    sflint::fmpq_poly_set_coeff_si(unit_poly, 0, 1);
    sflint::fmpq_poly_set_coeff_si(unit_poly, 1, 1);
    return out.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly));
}

// Appends `column` (rows x 1) to `matrix` (rows x n).
sflint::FmpzMat append_proof_column_for_test(const sflint::FmpzMat& matrix,
                                             const sflint::FmpzMat& column) {
    const slong rows = sflint::fmpz_mat_nrows(matrix);
    const slong cols = sflint::fmpz_mat_ncols(matrix);
    sflint::FmpzMat out(rows, cols + 1);
    for (slong i = 0; i < rows; ++i) {
        for (slong j = 0; j < cols; ++j) {
            sflint::fmpz_set(sflint::fmpz_mat_entry(out, i, j),
                             sflint::fmpz_mat_entry(matrix, i, j));
        }
        sflint::fmpz_set(sflint::fmpz_mat_entry(out, i, cols),
                         sflint::fmpz_mat_entry(column, i, 0));
    }
    return out;
}

// The stable proof scan's columns at q (stable_proof_prime_columns) against
// the residue-field route the scan used before, saturation_proof_prime_column
// at each prime of decompose_prime(order, q, 1): one column per degree-one
// prime, in the same order, each the direct root evaluation at that prime.
// The direct route takes logarithms against the ell-th power of the least
// primitive root mod q, the residue-field route against the first residue
// k = 1, 2, ... with k^((q - 1) / ell) != 1, so a column is the
// residue-field column times a unit modulo ell (equal for ell = 2), and the
// proof kernel after every column, which drives the scan, is the same.
// Q(zeta_13): q = 1 mod 13 splits completely, so every q below has 12 roots,
// and the equation order takes the word-size root route, whose root order
// the scan restores to decompose_prime's.
int test_stable_proof_prime_columns_match_residue_field_route() {
    silex::NumberField field = cyclotomic13_field();
    silex::Order equation = silex::test::equation_order(field);
    silex::Order maximal = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    const silex::Order* orders[] = {&equation, &maximal};
    // (ell, q): q = 1 mod ell and q = 1 mod 13.  For ell = 13, q = 79 the
    // least primitive root is 3 and 3^6 = 18 mod 79, while the residue-field
    // route takes 2^6 = 64 mod 79; for ell = 5, q = 521 they are
    // 3^104 = 25 and 2^104 = 516.  Those columns differ by a unit factor.
    // ell = 67 (7 bits) takes the baby-step giant-step branch of the word
    // route (kWordQuotientLogBabyGiantBits); q = 10453 = 1 mod 13 * 67 has
    // least primitive root 5, and 2^156 = (5^156)^57 mod 10453.
    const slong cases[][2] = {{2, 53}, {3, 79}, {5, 521}, {13, 79},
                              {67, 10453}};
    slong scaled_columns = 0;
    for (const silex::Order* order : orders) {
        silex::EmbeddingContext embeddings(field);
        silex::OrderUnitGroup group(*order);
        assert(set_cyclotomic13_u2_power17_group(group, *order, field,
                                                 embeddings));
        for (const auto& item : cases) {
            sflint::Fmpz ell;
            sflint::Fmpz q;
            assert(set_fmpz_si(ell, item[0]));
            assert(set_fmpz_si(q, item[1]));

            // The two quotient generators: the least primitive root g
            // raised to (q - 1) / ell, and k^((q - 1) / ell) for the first
            // k = 1, 2, ... where that is not 1; `factor` is the exponent t
            // with the second equal to the first to the t.
            const ulong q_ui = static_cast<ulong>(item[1]);
            const ulong ell_ui = static_cast<ulong>(item[0]);
            const ulong exponent = (q_ui - 1) / ell_ui;
            ulong primitive = 0;
            for (ulong g = 1; g < q_ui && primitive == 0; ++g) {
                bool generates = true;
                for (ulong f = 2; f <= q_ui - 1 && generates; ++f) {
                    if ((q_ui - 1) % f == 0 && n_is_prime(f) != 0 &&
                        n_powmod(g, static_cast<slong>((q_ui - 1) / f),
                                 q_ui) == 1) {
                        generates = false;
                    }
                }
                if (generates) {
                    primitive = g;
                }
            }
            ulong first = 0;
            for (ulong k = 1; k < q_ui && first == 0; ++k) {
                if (n_powmod(k, static_cast<slong>(exponent), q_ui) != 1) {
                    first = n_powmod(k, static_cast<slong>(exponent), q_ui);
                }
            }
            const ulong zeta = n_powmod(primitive,
                                        static_cast<slong>(exponent), q_ui);
            ulong factor = 0;
            for (ulong t = 1, power = zeta; t < ell_ui; ++t) {
                if (power == first) {
                    factor = t;
                    break;
                }
                power = n_mulmod2_preinv(power, zeta, q_ui,
                                         n_preinvert_limb(q_ui));
            }
            assert(factor != 0);

            std::vector<sflint::FmpzMat> columns;
            assert(silex::detail::stable_proof_prime_columns(
                    columns, group, sflint::FmpzConstRef(q),
                    sflint::FmpzConstRef(ell)));
            silex::PrimeIdealList primes;
            assert(silex::decompose_prime(primes, *order,
                                          sflint::FmpzConstRef(q), 1));
            assert(primes.size() == 12);
            assert(static_cast<slong>(columns.size()) == primes.size());

            const slong rows = sflint::fmpz_mat_nrows(columns[0]);
            sflint::FmpzMat scan_matrix(rows, 0);
            sflint::FmpzMat field_matrix(rows, 0);
            for (slong i = 0; i < primes.size(); ++i) {
                const silex::PrimeIdeal* prime = primes.at(i);
                assert(prime != nullptr);
                sflint::FmpzMat expected(0, 0);
                assert(silex::detail::saturation_proof_prime_column(
                        expected, group, *prime, sflint::FmpzConstRef(ell)));
                const sflint::FmpzMat& column =
                        columns[static_cast<std::size_t>(i)];
                assert(sflint::fmpz_mat_nrows(expected) == rows);
                assert(sflint::fmpz_mat_nrows(column) == rows);
                assert(sflint::fmpz_mat_ncols(column) == 1);

                // The column is the direct evaluation at this prime's root.
                sflint::Fmpz root;
                assert(degree_one_residue_root(root, *prime));
                sflint::FmpzMat direct(0, 0);
                assert(silex::detail::saturation_proof_prime_column_at_degree_one_root(
                        direct, group, sflint::FmpzConstRef(q),
                        sflint::FmpzConstRef(root),
                        sflint::FmpzConstRef(ell)));
                assert(fmpz_mat_equal(direct.raw(), column.raw()) != 0);

                // column = t * expected mod ell, entry by entry, where the
                // residue-field generator is the t-th power of the direct
                // one.
                for (slong r = 0; r < rows; ++r) {
                    const ulong b = fmpz_fdiv_ui(
                            fmpz_mat_entry(expected.raw(), r, 0), ell_ui);
                    assert(fmpz_cmp_ui(fmpz_mat_entry(column.raw(), r, 0),
                                       (factor * b) % ell_ui) == 0);
                }
                if (item[0] == 2) {
                    assert(fmpz_mat_equal(column.raw(), expected.raw()) != 0);
                }
                if (fmpz_mat_equal(column.raw(), expected.raw()) == 0) {
                    ++scaled_columns;
                }

                // The proof kernel after this column is the same.
                scan_matrix = append_proof_column_for_test(scan_matrix,
                                                           column);
                field_matrix = append_proof_column_for_test(field_matrix,
                                                            expected);
                sflint::FmpzMat scan_kernel(0, 0);
                sflint::FmpzMat field_kernel(0, 0);
                assert(silex::detail::dlog_kernel_from_matrix(
                        scan_kernel, scan_matrix, sflint::FmpzConstRef(ell)));
                assert(silex::detail::dlog_kernel_from_matrix(
                        field_kernel, field_matrix,
                        sflint::FmpzConstRef(ell)));
                assert(fmpz_mat_equal(scan_kernel.raw(), field_kernel.raw()) !=
                       0);
            }
        }
    }
    // The scaled case is exercised, not only literal equality.
    assert(scaled_columns > 0);
    return 0;
}

// Checks the stable proof scan's columns at q (stable_proof_prime_columns)
// against saturation_proof_prime_column at each of the `prime_count` primes
// of decompose_prime(order, q, 1), as
// test_stable_proof_prime_columns_match_residue_field_route does, but
// without computing the scale factor beforehand.  Where the direct root
// evaluation succeeds, the scan's next column is that evaluation and equals
// t times the residue-field column modulo ell, for one unit t per q (both
// builders' quotient generators depend only on q).  Where it declines, the
// scan falls back to the residue-field route at that prime; on the root
// path q does not divide [O : Z[theta]], so that route reduces the same
// theta-polynomials (residue_variable_numerator) and declines too, and the
// scan, like the residue-field route alone, has no column there.  The proof
// kernel after every column is the same for both.  Adds the direct columns
// with t != 1 to `scaled` and the declined primes to `declined`.
void check_stable_proof_prime_columns_at_q(slong& scaled,
                                           slong& declined,
                                           const silex::OrderUnitGroup& group,
                                           const silex::Order& order,
                                           sflint::FmpzConstRef q,
                                           slong ell,
                                           slong prime_count) {
    sflint::Fmpz ell_value;
    assert(set_fmpz_si(ell_value, ell));
    const sflint::FmpzConstRef ell_ref(ell_value);
    const ulong ell_ui = static_cast<ulong>(ell);
    std::vector<sflint::FmpzMat> columns;
    assert(silex::detail::stable_proof_prime_columns(columns, group, q,
                                                     ell_ref));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, q, 1));
    assert(primes.size() == prime_count);

    slong rows = 0;
    assert(silex::detail::dlog_proof_rank(rows, group, ell_ref));
    ulong factor = 0;
    std::size_t next = 0;
    sflint::FmpzMat scan_matrix(rows, 0);
    sflint::FmpzMat field_matrix(rows, 0);
    for (slong i = 0; i < prime_count; ++i) {
        const silex::PrimeIdeal* prime = primes.at(i);
        assert(prime != nullptr);
        sflint::Fmpz root;
        assert(degree_one_residue_root(root, *prime));
        sflint::FmpzMat direct(0, 0);
        sflint::FmpzMat expected(0, 0);
        if (!silex::detail::saturation_proof_prime_column_at_degree_one_root(
                    direct, group, q, sflint::FmpzConstRef(root), ell_ref)) {
            assert(!silex::detail::saturation_proof_prime_column(
                    expected, group, *prime, ell_ref));
            ++declined;
            continue;
        }
        assert(silex::detail::saturation_proof_prime_column(expected, group,
                                                            *prime, ell_ref));
        assert(next < columns.size());
        const sflint::FmpzMat& column = columns[next++];
        assert(sflint::fmpz_mat_nrows(expected) == rows);
        assert(sflint::fmpz_mat_nrows(column) == rows);
        assert(sflint::fmpz_mat_ncols(column) == 1);
        assert(fmpz_mat_equal(direct.raw(), column.raw()) != 0);
        for (slong r = 0; r < rows; ++r) {
            const ulong b = fmpz_fdiv_ui(fmpz_mat_entry(expected.raw(), r, 0),
                                         ell_ui);
            if (factor == 0 && b != 0) {
                factor = (fmpz_fdiv_ui(fmpz_mat_entry(column.raw(), r, 0),
                                       ell_ui) *
                          n_invmod(b, ell_ui)) %
                         ell_ui;
                assert(factor != 0);
            }
            assert(fmpz_cmp_ui(fmpz_mat_entry(column.raw(), r, 0),
                               (factor * b) % ell_ui) == 0);
        }
        if (fmpz_mat_equal(column.raw(), expected.raw()) == 0) {
            ++scaled;
        }

        scan_matrix = append_proof_column_for_test(scan_matrix, column);
        field_matrix = append_proof_column_for_test(field_matrix, expected);
        sflint::FmpzMat scan_kernel(0, 0);
        sflint::FmpzMat field_kernel(0, 0);
        assert(silex::detail::dlog_kernel_from_matrix(scan_kernel, scan_matrix,
                                                      ell_ref));
        assert(silex::detail::dlog_kernel_from_matrix(field_kernel,
                                                      field_matrix, ell_ref));
        assert(fmpz_mat_equal(scan_kernel.raw(), field_kernel.raw()) != 0);
    }
    assert(next == columns.size());
}

// The stable proof scan's columns on the routes that
// test_stable_proof_prime_columns_match_residue_field_route does not reach.
// Each case has a direct column that is a nontrivial multiple of the
// residue-field column, so the direct route, not the residue-field one,
// built the columns that the check finds direct.
// - The decompose_prime fallback: in Q(zeta_13) the first generator u_2 of
//   <u_2, ..., u_6> is stored as the factors (u_2 (zeta - r))^1
//   (zeta - r)^-1, which vanish at the root r of one prime over q = 79 and
//   at no other root.  The direct evaluation declines at r only; the
//   fallback finds the prime over r, whose residue-field column declines as
//   well, so the scan has the columns of the other 11 primes, in order, as
//   the residue-field route alone has.
// - The fmpz route of saturation_proof_prime_column_at_degree_one_root, for
//   q = 18446744073709555657 > 2^64, q = 1 mod 3 * 13.
// - The maximal-order root route, for the maximal order Z[theta / 2] of
//   theta^3 + 2 theta^2 - 8 theta - 8, which is not the equation order:
//   theta = 2 alpha with alpha^3 + alpha^2 - 2 alpha - 1 = 0, so the field is
//   Q(zeta_7)^+ and alpha, alpha^2 - 2 are units (conjugates 2 cos(2 pi / 7)
//   and 2 cos(4 pi / 7)); ell = 7 and q = 43 = 1 mod 7 splits completely.
int test_stable_proof_prime_columns_fallback_and_wide_routes() {
    silex::NumberField field = cyclotomic13_field();
    silex::Order equation = silex::test::equation_order(field);
    {
        sflint::Fmpz q;
        assert(set_fmpz_si(q, 79));
        silex::PrimeIdealList primes;
        assert(silex::decompose_prime(primes, equation,
                                      sflint::FmpzConstRef(q), 1));
        assert(primes.size() == 12);
        sflint::Fmpz root;
        assert(primes.at(4) != nullptr &&
               degree_one_residue_root(root, *primes.at(4)));

        sflint::FmpqPoly linear;
        sflint::fmpq_poly_set_coeff_si(linear, 1, 1);
        sflint::Fmpq constant;
        sflint::fmpq_set_fmpz(constant, sflint::FmpzConstRef(root));
        sflint::fmpq_neg(sflint::FmpqRef(constant),
                         sflint::FmpqConstRef(constant));
        sflint::fmpq_poly_set_coeff_fmpq(linear, 0, constant);
        silex::Element vanishing(field);
        silex::Element u2(field);
        silex::Element u2_times_vanishing(field);
        assert(vanishing.set_fmpq_poly(sflint::FmpqPolyConstRef(linear)));
        assert(set_cyclotomic13_u2(u2));
        assert(u2_times_vanishing.multiply(u2, vanishing));

        std::vector<silex::FactoredElement> generators;
        for (slong a = 2; a <= 6; ++a) {
            generators.emplace_back(field);
            if (a == 2) {
                assert(generators.back().push(u2_times_vanishing, 1));
                assert(generators.back().push(vanishing, -1));
                continue;
            }
            sflint::FmpqPoly unit_poly;
            for (slong i = 0; i < a; ++i) {
                sflint::fmpq_poly_set_coeff_si(unit_poly, i, 1);
            }
            silex::Element unit(field);
            assert(unit.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly)));
            assert(generators.back().set_element(unit));
        }
        silex::EmbeddingContext embeddings(field);
        silex::OrderUnitGroup group(equation);
        assert(group.set_units(equation,
                               silex::FactoredElementSpan(generators.data(),
                                                          generators.size()),
                               embeddings, 256));
        assert(group.free_rank() == 5);

        slong scaled = 0;
        slong declined = 0;
        check_stable_proof_prime_columns_at_q(scaled, declined, group,
                                              equation,
                                              sflint::FmpzConstRef(q), 13, 12);
        assert(declined == 1);
        assert(scaled > 0);
    }
    {
        sflint::Fmpz q;
        assert(sflint::fmpz_set_str(sflint::FmpzRef(q),
                                    "18446744073709555657"));
        assert(!sflint::fmpz_abs_fits_ui(sflint::FmpzConstRef(q)));
        silex::EmbeddingContext embeddings(field);
        silex::OrderUnitGroup group(equation);
        assert(set_cyclotomic13_u2_power17_group(group, equation, field,
                                                 embeddings));
        slong scaled = 0;
        slong declined = 0;
        check_stable_proof_prime_columns_at_q(scaled, declined, group,
                                              equation,
                                              sflint::FmpzConstRef(q), 3, 12);
        assert(declined == 0);
        assert(scaled > 0);
    }
    {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 2, 2);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -8);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, -8);
        silex::NumberField cubic = silex::test::field_by_polynomial(
                sflint::FmpqPolyConstRef(polynomial));
        silex::Order cubic_equation = silex::test::equation_order(cubic);
        silex::Order maximal(cubic);
        assert(maximal.maximal_order(cubic_equation));
        assert(maximal.is_maximal() && !maximal.is_equation_order());

        // alpha = theta / 2 and alpha^2 - 2 = theta^2 / 4 - 2.
        sflint::FmpqPoly alpha_poly;
        sflint::FmpqPoly conjugate_poly;
        sflint::Fmpq coefficient;
        sflint::fmpq_set_si(coefficient, 1, 2);
        sflint::fmpq_poly_set_coeff_fmpq(alpha_poly, 1, coefficient);
        sflint::fmpq_set_si(coefficient, 1, 4);
        sflint::fmpq_poly_set_coeff_fmpq(conjugate_poly, 2, coefficient);
        sflint::fmpq_poly_set_coeff_si(conjugate_poly, 0, -2);
        silex::Element alpha(cubic);
        silex::Element conjugate(cubic);
        assert(alpha.set_fmpq_poly(sflint::FmpqPolyConstRef(alpha_poly)));
        assert(conjugate.set_fmpq_poly(
                sflint::FmpqPolyConstRef(conjugate_poly)));
        std::vector<silex::FactoredElement> generators;
        generators.emplace_back(cubic);
        assert(generators.back().set_element(alpha));
        generators.emplace_back(cubic);
        assert(generators.back().set_element(conjugate));
        silex::EmbeddingContext embeddings(cubic);
        silex::OrderUnitGroup group(maximal);
        assert(group.set_units(maximal,
                               silex::FactoredElementSpan(generators.data(),
                                                          generators.size()),
                               embeddings, 256));
        assert(group.free_rank() == 2);

        sflint::Fmpz q;
        assert(set_fmpz_si(q, 43));
        slong scaled = 0;
        slong declined = 0;
        check_stable_proof_prime_columns_at_q(scaled, declined, group,
                                              maximal,
                                              sflint::FmpzConstRef(q), 7, 3);
        assert(declined == 0);
        assert(scaled > 0);
    }
    return 0;
}

// For an odd prime ell, as set_cyclotomic13_u2_power17_group, but the first
// `powered` generators u_a (a = 2, ..., powered + 1) are replaced by
// -u_a^ell = (-u_a)^ell, stored as the factors u_a^ell and (-1)^1: a
// subgroup of index ell^powered in the full unit group.  The ell kernel row
// of such a generator is an ell-th power but not a formal one: the formal
// split leaves the rational -1, which Element::is_power roots in any
// degree, so the root -u_a is computed, unlike the root of the expanded
// u_2^17.
bool set_cyclotomic13_minus_unit_power_group(
        silex::OrderUnitGroup& out,
        const silex::Order& order,
        const silex::NumberField& field,
        silex::EmbeddingContext& embeddings,
        slong ell,
        slong powered) noexcept {
    silex::Element minus_one(field);
    sflint::Fmpz minus_one_value;
    if (!set_fmpz_si(minus_one_value, -1) ||
        !minus_one.set_fmpz(sflint::FmpzConstRef(minus_one_value))) {
        return false;
    }
    std::vector<silex::FactoredElement> generators;
    for (slong a = 2; a <= 6; ++a) {
        sflint::FmpqPoly unit_poly;
        for (slong i = 0; i < a; ++i) {
            sflint::fmpq_poly_set_coeff_si(unit_poly, i, 1);
        }
        silex::Element unit(field);
        if (!unit.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly))) {
            return false;
        }
        generators.emplace_back(field);
        const bool pushed = a <= powered + 1
                                    ? generators.back().push(unit, ell) &&
                                             generators.back().push(
                                                     minus_one, 1)
                                   : generators.back().set_element(unit);
        if (!pushed) {
            return false;
        }
    }
    return out.set_units(order,
                         silex::FactoredElementSpan(generators.data(),
                                                    generators.size()),
                         embeddings, 256) &&
           out.free_rank() == 5;
}

// Counts the stable proof's deferred and accepted kernel-row roots in the
// first prove_local_saturated_stable_in_place_ pass of a run: with
// prove_index_bound, the pass for the smallest ell >= aux_bound.
struct FirstStablePassRootEvents {
    bool active = false;
    slong passes = 0;
    slong deferred = 0;
    slong accepted = 0;
};

void first_stable_pass_root_callback(void* user,
                                     silex::DiagnosticsModule,
                                     silex::ProfileEvent event,
                                     const char*,
                                     const char* label) noexcept {
    FirstStablePassRootEvents* events =
            static_cast<FirstStablePassRootEvents*>(user);
    if (label == nullptr) {
        return;
    }
    if (std::strcmp(label,
                    "unit_group.prove_local_saturated_stable_in_place") == 0) {
        if (event == silex::ProfileEvent::begin_scope) {
            events->active = events->passes == 0;
        } else if (event == silex::ProfileEvent::end_scope) {
            events->active = false;
            ++events->passes;
        }
        return;
    }
    if (!events->active || event != silex::ProfileEvent::event) {
        return;
    }
    if (std::strcmp(label,
                    "unit_group.stable_relation_saturation.root_deferred") ==
        0) {
        ++events->deferred;
    } else if (std::strcmp(label, "unit_group.stable_relation_saturation."
                                  "accepted_root") == 0) {
        ++events->accepted;
    }
}

// Q(sqrt 118) with the fundamental unit eps = 306917 + 28254 sqrt(118).  GP:
// K = bnfinit(x^2 - 118, 1); K.fu gives eps up to sign, and eps is a cube
// at both primes above q = 13 and q = 19, the split primes q = 1 mod 3
// below 31, but not at the primes above q = 31.
bool set_quadratic118_fundamental_unit_group(
        silex::OrderUnitGroup& out,
        const silex::Order& order,
        const silex::NumberField& field,
        silex::EmbeddingContext& embeddings) noexcept {
    sflint::FmpqPoly unit_poly;
    sflint::fmpq_poly_set_coeff_si(unit_poly, 0, 306917);
    sflint::fmpq_poly_set_coeff_si(unit_poly, 1, 28254);
    silex::Element unit(field);
    silex::FactoredElement generator(field);
    return unit.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly)) &&
           generator.set_element(unit) &&
           out.set_units(order, silex::FactoredElementSpan(&generator, 1),
                         embeddings, 256) &&
           out.free_rank() == 1;
}

silex::NumberField quadratic118_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -118);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// kernel_row_root tells an error apart from an undecided root, so the proof
// scan fails closed at once on an error and rescans only an undecided row.
int test_kernel_row_root_three_way_result() {
    {
        silex::NumberField field = quadratic118_field();
        silex::Order order = silex::test::verified_maximal_order(
                silex::test::equation_order(field));
        silex::EmbeddingContext embeddings(field);
        silex::OrderUnitGroup group(order);
        assert(set_quadratic118_fundamental_unit_group(group, order, field,
                                                       embeddings));
        sflint::Fmpz ell;
        assert(set_fmpz_si(ell, 3));
        sflint::FmpzMat rows(2, 1);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 1);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 1, 0), 3);

        // eps is not a cube (it is fundamental): decided, not a power.
        silex::FactoredElement root(field);
        bool is_power = true;
        assert(silex::detail::kernel_row_root(
                       is_power, root, group, sflint::FmpzMatConstRef(rows), 0,
                       sflint::FmpzConstRef(ell)) ==
               silex::detail::KernelRowRootResult::decided);
        assert(!is_power);
        // eps^3 is a cube: decided, with an exactly verified root.
        is_power = false;
        assert(silex::detail::kernel_row_root(
                       is_power, root, group, sflint::FmpzMatConstRef(rows), 1,
                       sflint::FmpzConstRef(ell)) ==
               silex::detail::KernelRowRootResult::decided);
        assert(is_power);

        // Errors: a row out of range, a kernel width that is neither the
        // free rank nor the free rank plus one, and a composite exponent.
        assert(silex::detail::kernel_row_root(
                       is_power, root, group, sflint::FmpzMatConstRef(rows), 2,
                       sflint::FmpzConstRef(ell)) ==
               silex::detail::KernelRowRootResult::error);
        sflint::FmpzMat wide(1, 3);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(wide, 0, 0), 1);
        assert(silex::detail::kernel_row_root(
                       is_power, root, group, sflint::FmpzMatConstRef(wide), 0,
                       sflint::FmpzConstRef(ell)) ==
               silex::detail::KernelRowRootResult::error);
        sflint::Fmpz composite;
        assert(set_fmpz_si(composite, 9));
        assert(silex::detail::kernel_row_root(
                       is_power, root, group, sflint::FmpzMatConstRef(rows), 0,
                       sflint::FmpzConstRef(composite)) ==
               silex::detail::KernelRowRootResult::error);
    }

    // u_2^17 is a 17th power, but its root u_2 is not computed in degree 12:
    // undecided, not an error.
    silex::NumberField field = cyclotomic13_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group(order);
    assert(set_cyclotomic13_u2_power17_group(group, order, field, embeddings));
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 17));
    sflint::FmpzMat rows(1, 5);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 1);
    silex::FactoredElement root(field);
    bool is_power = false;
    assert(silex::detail::kernel_row_root(
                   is_power, root, group, sflint::FmpzMatConstRef(rows), 0,
                   sflint::FmpzConstRef(ell)) ==
           silex::detail::KernelRowRootResult::undecided);
    return 0;
}

// KernelRowRootMode::defer_undecidable skips the root attempt only where
// Element::is_power cannot find the root of a non-rational l-th power
// (l >= 3, degree >= 10) and the row product is not a formal l-th power.
int test_kernel_row_root_defers_undecidable_rows() {
    {
        // The premise: in degree 12, Element::is_power with exponent 3
        // leaves a non-rational cube undecided and decides a rational one.
        // If Element::is_power learns to find such roots here, update
        // element_power_test_roots_non_rational (relation_kernel.cpp).
        silex::NumberField field = cyclotomic13_field();
        sflint::FmpqPoly theta_poly;
        sflint::fmpq_poly_set_coeff_si(theta_poly, 1, 1);
        sflint::fmpq_poly_set_coeff_si(theta_poly, 0, 2);
        silex::Element base(field);
        silex::Element cube(field);
        silex::Element root(field);
        sflint::Fmpz three;
        assert(set_fmpz_si(three, 3));
        assert(base.set_fmpq_poly(sflint::FmpqPolyConstRef(theta_poly)));
        assert(cube.pow_fmpz(base, sflint::FmpzConstRef(three)));
        bool is_power = false;
        assert(!cube.is_power(is_power, root, sflint::FmpzConstRef(three)));
        silex::Element eight(field);
        sflint::Fmpz eight_value;
        assert(set_fmpz_si(eight_value, 8));
        assert(eight.set_fmpz(sflint::FmpzConstRef(eight_value)));
        assert(eight.is_power(is_power, root, sflint::FmpzConstRef(three)));
        assert(is_power);
    }

    silex::NumberField field = cyclotomic13_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group(order);
    assert(set_cyclotomic13_u2_power17_group(group, order, field, embeddings));
    sflint::Fmpz ell;
    assert(set_fmpz_si(ell, 17));
    sflint::FmpzMat rows(2, 5);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 1, 1), 17);
    silex::FactoredElement root(field);
    bool is_power = false;
    // u_2^17 (expanded): degree 12 and l = 17, so the row is deferred
    // instead of ending undecided after the full attempt.
    assert(silex::detail::kernel_row_root(
                   is_power, root, group, sflint::FmpzMatConstRef(rows), 0,
                   sflint::FmpzConstRef(ell), nullptr,
                   silex::detail::KernelRowRootMode::defer_undecidable) ==
           silex::detail::KernelRowRootResult::deferred);
    // u_3^17 is a formal 17th power: still tested, and decided with the
    // formal root u_3.
    assert(silex::detail::kernel_row_root(
                   is_power, root, group, sflint::FmpzMatConstRef(rows), 1,
                   sflint::FmpzConstRef(ell), nullptr,
                   silex::detail::KernelRowRootMode::defer_undecidable) ==
           silex::detail::KernelRowRootResult::decided);
    assert(is_power);
    silex::FactoredElement u3(field);
    assert(group.free_generator(u3, 1));
    std::vector<silex::FactoredElement> got;
    std::vector<silex::FactoredElement> want;
    got.emplace_back(field);
    want.emplace_back(field);
    assert(got.back().set(root));
    assert(want.back().set(u3));
    assert(same_factor_lists(got, want));

    // Degree 2: Element::is_power decides cubes, so the row is tested.
    silex::NumberField quadratic = quadratic118_field();
    silex::Order quadratic_order = silex::test::verified_maximal_order(
            silex::test::equation_order(quadratic));
    silex::EmbeddingContext quadratic_embeddings(quadratic);
    silex::OrderUnitGroup quadratic_group(quadratic_order);
    assert(set_quadratic118_fundamental_unit_group(
            quadratic_group, quadratic_order, quadratic,
            quadratic_embeddings));
    sflint::Fmpz ell3;
    assert(set_fmpz_si(ell3, 3));
    sflint::FmpzMat eps_row(1, 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(eps_row, 0, 0), 1);
    silex::FactoredElement quadratic_root(quadratic);
    is_power = true;
    assert(silex::detail::kernel_row_root(
                   is_power, quadratic_root, quadratic_group,
                   sflint::FmpzMatConstRef(eps_row), 0,
                   sflint::FmpzConstRef(ell3), nullptr,
                   silex::detail::KernelRowRootMode::defer_undecidable) ==
           silex::detail::KernelRowRootResult::decided);
    assert(!is_power);
    return 0;
}

// Without require_y_root (the bounded relation-kernel search), when the exact
// m-th root search cannot complete and the torsion check verifies the
// relation for y, the element returned and adjoined is y, as in reference
// Unit/UnitGrpCtx.jl `_add_dependent_unit!` (`_transform(vcat(U.units, y),
// m)`), and not whatever `root` held before.  Against <u_2^17, u_3, ...,
// u_6> in Q(zeta_13), y = u_2 satisfies y^17 = g_1, but the root of the
// expanded u_2^17 is not computed in degree 12.
int test_dependent_relation_unbounded_root_failure_adjoins_y() {
    silex::NumberField field = cyclotomic13_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup group(order);
    assert(set_cyclotomic13_u2_power17_group(group, order, field, embeddings));

    // The premise: the compact 17th-root search of g_1 = u_2^17 fails, so
    // the power-root search returns false and the torsion check decides.
    silex::FactoredElement g1(field);
    assert(group.free_generator(g1, 0));
    silex::FactoredElement g1_root(field);
    bool is_power = false;
    assert(!g1.is_power_si(is_power, g1_root, 17,
                           silex::FactoredRootStrategy::compact));

    silex::Element u2(field);
    assert(set_cyclotomic13_u2(u2));
    silex::FactoredElement y(field);
    assert(y.set_element(u2));

    // A stale value in `root`, as left by an earlier kernel generator.
    silex::FactoredElement root(field);
    assert(group.free_generator(root, 1));

    sflint::FmpzMat relation(1, 6);
    sflint::Fmpz torsion_exp;
    sflint::Fmpz bound;
    assert(set_fmpz_si(bound, 17));
    bool recovered = false;
    assert(silex::detail::dependent_relation_bounded_min_denominator(
            recovered, root, relation, torsion_exp, group, y, embeddings,
            sflint::FmpzConstRef(bound), 2, 128, 256, false, false));
    assert(recovered);
    const slong expected[] = {-1, 0, 0, 0, 0, 17};
    assert(candidate_row_equals(relation, expected, 6));
    silex::Element root_value(field);
    assert(root.evaluate(root_value));
    assert(root_value.equal(u2));

    // Adjoining it gives <U, y>, the full unit group.
    silex::OrderUnitGroup refined(order);
    bool changed = false;
    assert(silex::detail::adjoin_verified_dependent_relation(
            changed, refined, group, root, sflint::FmpzMatConstRef(relation),
            0, embeddings, 256));
    assert(changed);
    assert(refined.free_rank() == 5);
    // GP: \p 110; K.reg (see test_prove_index_bound_cyclotomic13_undecided_root)
    assert(regulator_contains(
            refined,
            "120.7840313626563127293175827399128565840801560246896328383867"
            "406585225412356786592960914162040943733"));
    return 0;
}

// A kernel row that is not an ell-th power, met when the scan has used every
// prime q of its range, fails closed as `unavailable`: a rescan would return
// the same kernel, so doubling `stable` would never end.  The test hook ends
// the q range at 30, so the ell = 3 scan sees only q = 13 and q = 19, where
// eps is a cube; without the limit the scan reaches q = 31 and proves ell = 3.
int test_prove_index_bound_non_power_range_exhausted() {
    silex::NumberField field = quadratic118_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup start(order);
    assert(set_quadratic118_fundamental_unit_group(start, order, field,
                                                   embeddings));
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 30));

    {
        silex::OrderUnitGroup proved(order);
        silex::ProofState status = silex::ProofState::not_checked;
        bool changed = true;
        assert(proved.prove_index_bound(status, changed, start, 1,
                                        sflint::FmpzConstRef(aux_bound), 2,
                                        embeddings, 256));
        assert(status == silex::ProofState::verified);
        sflint::Fmpz ell3;
        assert(set_fmpz_si(ell3, 3));
        assert(proved.unit_proof_verified(sflint::FmpzConstRef(ell3)));
    }

    const StableProofScanLimitsForTesting limits(0, 30);
    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_index_bound(status, changed, start, 1,
                                    sflint::FmpzConstRef(aux_bound), 2,
                                    embeddings, 256));
    assert(status == silex::ProofState::unavailable);
    assert(!changed);
    assert(proved.certification_status() != silex::CertificationMode::proven);
    const slong last = proved.unit_proof_record_count() - 1;
    // Two degree-one primes above each of q = 13 and q = 19.
    assert(check_unit_proof_record(proved, last, 3,
                                   silex::ProofState::unavailable, 30, 4,
                                   false));
    return 0;
}

// Counts, in prove_local_saturated, the outcomes of the stable selection
// run before the pre-scan's kernel-row root attempts, and the root passes
// (saturate_local_with_kernel scopes) that follow it.
struct PrescanStableFirstEvents {
    slong verified = 0;
    slong not_verified = 0;
    slong root_passes = 0;
};

void prescan_stable_first_callback(void* user,
                                   silex::DiagnosticsModule,
                                   silex::ProfileEvent event,
                                   const char*,
                                   const char* label) noexcept {
    PrescanStableFirstEvents* events =
            static_cast<PrescanStableFirstEvents*>(user);
    if (label == nullptr) {
        return;
    }
    if (event == silex::ProfileEvent::begin_scope &&
        std::strcmp(label, "unit_group.saturate_local_with_kernel") == 0) {
        ++events->root_passes;
    } else if (event == silex::ProfileEvent::event &&
               std::strcmp(label,
                           "unit_group.proof_prescan_stable_first_verified") ==
                       0) {
        ++events->verified;
    } else if (event == silex::ProfileEvent::event &&
               std::strcmp(label, "unit_group.proof_prescan_stable_first_"
                                  "not_verified") == 0) {
        ++events->not_verified;
    }
}

// The pre-scan of a prove_index_bound pass that leaves a nonempty kernel
// first runs the selection that opens the reference stable scan, and
// publishes its result when that kernel is empty, without testing any
// pre-scan row for a root.  Otherwise it tests the pre-scan rows and then
// runs the stable scan, as before.  In Q(sqrt 22), eps = 197 + 42 sqrt(22)
// (GP: K = bnfinit(x^2 - 22, 1); K.fu, up to sign) is a cube at both primes
// above q = 7 and at neither prime above q = 13.  With aux_bound = 10 the ell = 3 pre-scan sees only q = 7, so
// eps stays in its kernel; the stable scan continues to q = 13.  With the q
// range ended at 10 it cannot.
int test_prove_index_bound_prescan_stable_first() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -22);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup start(order);
    {
        sflint::FmpqPoly unit_poly;
        sflint::fmpq_poly_set_coeff_si(unit_poly, 0, 197);
        sflint::fmpq_poly_set_coeff_si(unit_poly, 1, 42);
        silex::Element unit(field);
        silex::FactoredElement generator(field);
        assert(unit.set_fmpq_poly(sflint::FmpqPolyConstRef(unit_poly)));
        assert(generator.set_element(unit));
        assert(start.set_units(order,
                               silex::FactoredElementSpan(&generator, 1),
                               embeddings, 256));
        assert(start.free_rank() == 1);
    }
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 10));

    for (const bool limit_range : {false, true}) {
        silex::DiagnosticsContext diagnostics;
        silex::diagnostics_context_init(diagnostics);
        PrescanStableFirstEvents events;
        silex::diagnostics_set_profiling(
                diagnostics, true,
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group),
                prescan_stable_first_callback, &events);
        const StableProofScanLimitsForTesting limits(0, limit_range ? 10 : 0);
        silex::OrderUnitGroup proved(order);
        proved.set_diagnostics(&diagnostics);
        silex::ProofState status = silex::ProofState::not_checked;
        bool changed = true;
        assert(proved.prove_index_bound(status, changed, start, 1,
                                        sflint::FmpzConstRef(aux_bound), 2,
                                        embeddings, 256));
        assert(!changed);
        // The records are those of the path without the early selection
        // (checked against it): ell = 2, then ell = 3.
        assert(proved.unit_proof_record_count() == 2);
        assert(check_unit_proof_record(proved, 0, 2,
                                       silex::ProofState::verified, 10, 3,
                                       false));
        if (!limit_range) {
            // The stable selection empties the kernel at the first prime
            // above q = 13, after both primes above q = 7, where eps is a
            // cube: ell = 3 is verified and eps is never tested for a root.
            assert(status == silex::ProofState::verified);
            assert(proved.certification_status() ==
                   silex::CertificationMode::proven);
            assert(check_unit_proof_record(proved, 1, 3,
                                           silex::ProofState::verified, 10, 3,
                                           false));
        } else {
            // The range ends before q = 13: the selection keeps eps, the
            // pre-scan tests it (not a cube), and the stable fallback fails
            // closed after the two primes above q = 7.
            assert(status == silex::ProofState::unavailable);
            assert(check_unit_proof_record(proved, 1, 3,
                                           silex::ProofState::unavailable, 10,
                                           2, false));
        }
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
        assert(events.verified == (limit_range ? 0 : 1));
        assert(events.not_verified == (limit_range ? 1 : 0));
        assert(events.root_passes == (limit_range ? 1 : 0));
#endif
    }
    return 0;
}

// A deferred row that is a true power is rooted after the rescan (decision
// for deferred root attempts).  In Q(zeta_13), the group of
// set_cyclotomic13_minus_unit_power_group with ell = 17 and one powered
// generator has index 17, and aux_bound = 17
// sends ell = 17 to the stable proof.  Degree 12 and ell = 17, so the first
// scan, at stable = 3.5, defers the non-formal row of -u_2^17 without
// testing it; the rescan at stable = 7 tests it and accepts the root -u_2,
// and the next scan has an empty kernel.  An attempt in the first scan
// would accept the same root at stable = 3.5.  The proof record then has
// the same local primes either way, since a verifying scan stops at the
// first prime that empties the kernel, whatever `stable` is; the published
// generators may differ.  The resource guard (4096 primes q per scan, far
// above what this proof needs) makes a scan that never tests the row fail
// closed quickly instead of after 2^16 primes q.
int test_prove_index_bound_roots_power_after_deferred_rescan() {
    silex::NumberField field = cyclotomic13_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup start(order);
    assert(set_cyclotomic13_minus_unit_power_group(start, order, field,
                                                   embeddings, 17, 1));
    sflint::Fmpz ell17;
    assert(set_fmpz_si(ell17, 17));

    // The row is deferred under defer_undecidable, and rooted when tested.
    sflint::FmpzMat rows(1, 5);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(rows, 0, 0), 1);
    silex::FactoredElement root(field);
    bool is_power = false;
    assert(silex::detail::kernel_row_root(
                   is_power, root, start, sflint::FmpzMatConstRef(rows), 0,
                   sflint::FmpzConstRef(ell17), nullptr,
                   silex::detail::KernelRowRootMode::defer_undecidable) ==
           silex::detail::KernelRowRootResult::deferred);
    assert(silex::detail::kernel_row_root(
                   is_power, root, start, sflint::FmpzMatConstRef(rows), 0,
                   sflint::FmpzConstRef(ell17)) ==
           silex::detail::KernelRowRootResult::decided);
    assert(is_power);
    silex::Element u2(field);
    silex::Element minus_u2(field);
    silex::Element root_value(field);
    assert(set_cyclotomic13_u2(u2));
    assert(minus_u2.negate(u2));
    assert(root.evaluate(root_value));
    assert(root_value.equal(minus_u2));

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    FirstStablePassRootEvents events;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            first_stable_pass_root_callback, &events);
    const StableProofScanLimitsForTesting limits(4096, 0);
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 17));
    silex::OrderUnitGroup proved(order);
    proved.set_diagnostics(&diagnostics);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = false;
    assert(proved.prove_index_bound(status, changed, start, 1,
                                    sflint::FmpzConstRef(aux_bound), 2,
                                    embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    // The ell = 17 record: verified, the group changed, and 13 degree-one
    // primes (12 above q = 443 and the first above q = 1327).
    bool found = false;
    for (slong i = 0; i < proved.unit_proof_record_count(); ++i) {
        found = found || check_unit_proof_record(
                                 proved, i, 17, silex::ProofState::verified,
                                 17, 13, true);
    }
    assert(found);
    // The root was adjoined: the regulator is that of the full unit group.
    // GP: \p 110; K.reg
    assert(regulator_contains(
            proved,
            "120.7840313626563127293175827399128565840801560246896328383867"
            "406585225412356786592960914162040943733"));
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The ell = 17 pass deferred the row once, then accepted its root.
    assert(events.deferred == 1);
    assert(events.accepted == 1);
#endif
    return 0;
}

// When the stable scan has used every prime q of its range, rows are
// tested at once, not deferred: a rescan would return the same kernel.  In
// Q(zeta_13), every generator of the group of
// set_cyclotomic13_minus_unit_power_group with ell = 19 and five powered
// generators is a 19th power, so the ell = 19 character kernel stays full
// (dimension 5), and the scan stops only after more than 3.5 * 5 primes
// without a change in it.  The test hook ends the q range at 2000: the
// primes ell <= 17 are still proved (ell = 17 needs the first degree-one
// prime above q = 1327), and the ell = 19 scan, sent to the stable proof by
// aux_bound = 19, sees only the 12 degree-one primes above q = 1483 (the
// first prime q = 1 mod 19 that splits in Q(zeta_13)) and ends with the
// range exhausted.  The first row, -u_2^19, is then rooted in that first
// scan; a deferral would have rescanned first.  The published result is
// the same either way, so only the profile events, counted in builds with
// profiling enabled, tell the two apart.  After the root the next row
// tested is not a 19th power, and the range is still exhausted, so ell = 19
// fails closed as `unavailable`.
int test_prove_index_bound_range_exhausted_tests_rows_at_once() {
    silex::NumberField field = cyclotomic13_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup start(order);
    assert(set_cyclotomic13_minus_unit_power_group(start, order, field,
                                                   embeddings, 19, 5));

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    FirstStablePassRootEvents events;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            first_stable_pass_root_callback, &events);
    const StableProofScanLimitsForTesting limits(0, 2000);
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 19));
    silex::OrderUnitGroup proved(order);
    proved.set_diagnostics(&diagnostics);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = false;
    assert(proved.prove_index_bound(status, changed, start, 1,
                                    sflint::FmpzConstRef(aux_bound), 2,
                                    embeddings, 256));
    assert(status == silex::ProofState::unavailable);
    assert(changed);
    assert(proved.certification_status() != silex::CertificationMode::proven);
    sflint::Fmpz ell17;
    sflint::Fmpz ell19;
    assert(set_fmpz_si(ell17, 17));
    assert(set_fmpz_si(ell19, 19));
    assert(proved.unit_proof_verified(sflint::FmpzConstRef(ell17)));
    const slong last = proved.unit_proof_record_count() - 1;
    // The 12 degree-one primes above q = 1483.
    assert(check_unit_proof_record(proved, last, 19,
                                   silex::ProofState::unavailable, 19, 12,
                                   false));
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The ell = 19 pass tested the row at once and accepted its root.
    assert(events.deferred == 0);
    assert(events.accepted == 1);
#endif
    return 0;
}

// The undecided-root rescans of a group that is not 17-saturated: the row
// for u_2^17 stays in every character kernel, so the proof fails closed as
// `unavailable` and never certifies ell = 17.  Both stops are tested through
// the test hook: the per-scan guard (here 256 primes q instead of 2^16,
// which takes minutes) and the end of the q range (here q <= 2000), where
// the undecided row fails closed at once instead of rescanning forever.
int test_prove_index_bound_undecided_root_not_saturated() {
    silex::NumberField field = cyclotomic13_field();
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);
    silex::OrderUnitGroup start(order);
    assert(set_cyclotomic13_u2_power17_group(start, order, field, embeddings));
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 1000));
    sflint::Fmpz ell17;
    assert(set_fmpz_si(ell17, 17));

    for (int run = 0; run < 2; ++run) {
        // run 0: guard; run 1: end of range.  At q <= 2000 the primes
        // q = 1 mod 17 that split in Q(zeta_13) are 443 and 1327 (q = 1
        // mod 221), 24 degree-one primes.  A scan stopped by the guard
        // records no local primes.
        const StableProofScanLimitsForTesting limits(run == 0 ? 256 : 0,
                                                     run == 0 ? 0 : 2000);
        silex::OrderUnitGroup proved(order);
        silex::ProofState status = silex::ProofState::not_checked;
        bool changed = true;
        assert(proved.prove_index_bound(status, changed, start, 1,
                                        sflint::FmpzConstRef(aux_bound), 2,
                                        embeddings, 256));
        assert(status == silex::ProofState::unavailable);
        assert(!changed);
        assert(proved.certification_status() !=
               silex::CertificationMode::proven);
        assert(!proved.unit_proof_verified(sflint::FmpzConstRef(ell17)));
        const slong last = proved.unit_proof_record_count() - 1;
        assert(check_unit_proof_record(proved, last, 17,
                                       silex::ProofState::unavailable, 1000,
                                       run == 0 ? 0 : 24, false));
    }
    return 0;
}

// The stronger regulator bound must not change the published generators.
// x^3 - 4x^2 - 8x - 1 (discriminant 2213): the units theta^2 and
// -4 theta^2 - theta + 6 span a subgroup of index 2 in the fundamental
// units (see above).  saturate_index_bounded scans every prime ell up to
// the regulator index bound: 5 on the maximal order (Pohst-Zassenhaus
// (6.22)), 24 on the equation order, which is not known to be maximal and
// keeps the signature bound.  Z[theta] is the maximal order, so both runs
// saturate the same lattice and must publish identical generators.
int test_saturate_index_bounded_generators_independent_of_bound() {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, -8);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order maximal = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::Order equation = silex::test::equation_order(field);
    assert(!equation.maximality_known());
    silex::EmbeddingContext embeddings(field);

    sflint::FmpqPoly theta_poly;
    sflint::fmpq_poly_set_coeff_si(theta_poly, 2, 1);
    sflint::FmpqPoly second_poly;
    sflint::fmpq_poly_set_coeff_si(second_poly, 2, -4);
    sflint::fmpq_poly_set_coeff_si(second_poly, 1, -1);
    sflint::fmpq_poly_set_coeff_si(second_poly, 0, 6);
    silex::Element theta_squared(field);
    silex::Element second(field);
    assert(theta_squared.set_fmpq_poly(sflint::FmpqPolyConstRef(theta_poly)));
    assert(second.set_fmpq_poly(sflint::FmpqPolyConstRef(second_poly)));
    silex::FactoredElement generators[] = {silex::FactoredElement(field),
                                           silex::FactoredElement(field)};
    assert(generators[0].set_element(theta_squared));
    assert(generators[1].set_element(second));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 1000));
    std::vector<silex::FactoredElement> published[2];
    const silex::Order* orders[2] = {&maximal, &equation};
    for (int run = 0; run < 2; ++run) {
        const silex::Order& order = *orders[run];
        silex::OrderUnitGroup start(order);
        assert(start.set_units(order, silex::FactoredElementSpan(generators, 2),
                               embeddings, 256));
        sflint::Fmpz index_bound;
        assert(start.regulator_index_bound(sflint::FmpzRef(index_bound), 256));
        assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(index_bound),
                                     run == 0 ? 10 : 48));
        silex::OrderUnitGroup out(order);
        bool changed = false;
        bool stable = false;
        assert(out.saturate_index_bounded(changed, stable, start, embeddings,
                                          4, sflint::FmpzConstRef(aux_bound),
                                          4, 256));
        assert(changed);
        assert(out.free_rank() == 2);
        for (slong j = 0; j < 2; ++j) {
            silex::FactoredElement unit(field);
            assert(out.free_generator(unit, j));
            published[run].emplace_back(field);
            assert(published[run].back().set(unit));
        }
        assert(regulator_contains(
                out,
                "12.68082022713452093603368239737775544244960210761943928030154"
                "4077156925247307418251826663747356303483"));
    }
    assert(same_factor_lists(published[0], published[1]));
    return 0;
}

// x^4 + 4x^3 + x^2 + 2x - 3 (discriminant -70640, signature (2, 1), Z[theta]
// maximal): GP 2.17 `K = bnfinit(f, 1)` gives class number one, w = 2,
// bnfcertify(K) = 1 and K.reg = 14.46098754389304356...; the units
// theta^3 + 4 theta^2 - 2 and -theta^3 - 5 theta^2 - 6 theta - 7 have norm 1
// and their log-embedding determinant equals K.reg, so they are fundamental.
// The Friedman 1989 Table 6 bound for signature (2, 1) is R >= 0.36, and
// floor(14.46098... / 0.36) = floor(40.169...) = 40.  GP:
//     forprime(q=2, 1000, if(q%31==1 && #polrootsmod(f,q), print(q)))
// prints only 311 (one root), so the 31-local proof of this rank-two group
// must continue past the pre-scan bound 1000.
int test_prove_index_bound_quartic70640() {
    silex::NumberField field = quartic_field(4, 1, 2, -3);
    silex::Order order = silex::test::verified_maximal_order(
            silex::test::equation_order(field));
    silex::EmbeddingContext embeddings(field);

    sflint::FmpqPoly first_poly;
    sflint::fmpq_poly_set_coeff_si(first_poly, 3, 1);
    sflint::fmpq_poly_set_coeff_si(first_poly, 2, 4);
    sflint::fmpq_poly_set_coeff_si(first_poly, 0, -2);
    sflint::FmpqPoly second_poly;
    sflint::fmpq_poly_set_coeff_si(second_poly, 3, -1);
    sflint::fmpq_poly_set_coeff_si(second_poly, 2, -5);
    sflint::fmpq_poly_set_coeff_si(second_poly, 1, -6);
    sflint::fmpq_poly_set_coeff_si(second_poly, 0, -7);
    silex::Element first(field);
    silex::Element second(field);
    assert(first.set_fmpq_poly(sflint::FmpqPolyConstRef(first_poly)));
    assert(second.set_fmpq_poly(sflint::FmpqPolyConstRef(second_poly)));
    silex::FactoredElement generators[] = {silex::FactoredElement(field),
                                           silex::FactoredElement(field)};
    assert(generators[0].set_element(first));
    assert(generators[1].set_element(second));
    silex::OrderUnitGroup start(order);
    assert(start.set_units(order, silex::FactoredElementSpan(generators, 2),
                           embeddings, 256));
    assert(start.free_rank() == 2);

    sflint::Fmpz index_bound;
    assert(start.regulator_index_bound(sflint::FmpzRef(index_bound), 256));
    assert(sflint::fmpz_equal_si(sflint::FmpzConstRef(index_bound), 40));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 1000));
    silex::OrderUnitGroup proved(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = true;
    assert(proved.prove_index_bound(status, changed, start, 1,
                                    sflint::FmpzConstRef(aux_bound), 2,
                                    embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(!changed);
    assert(proved.certification_status() == silex::CertificationMode::proven);
    const slong kProvedElls[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (slong ell : kProvedElls) {
        assert(has_unit_proof_record(proved, ell, silex::ProofState::verified,
                                     false));
    }
    assert(!has_unit_proof_record(proved, 41, silex::ProofState::verified,
                                  false));
    // GP: \p 110; K.reg
    assert(regulator_contains(
            proved,
            "14.46098754389304356185892633589732684899498168515774627523333"
            "6609135104418249234251372768481728608544835148228"));
    return 0;
}

// u = -(3 + 2 sqrt2) = -eps^2, eps = 1 + sqrt2, is a square in Q(sqrt2)
// only up to torsion: u * (-1) = eps^2.  Every residue character of order 2
// takes the same value on u and on -1, so the 2-saturation candidates must
// include the torsion generator, as reference
// `compute_candidates_for_saturate` and `saturate!` do (see
// docs/reference/algorithms_and_sources.rst); the candidate u * (-1) has
// the root +/-eps.  A candidate set without torsion never adjoins that root
// and fails closed.
int test_prove_torsion_twisted_saturation_root() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element twisted(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(twisted.multiply(epsilon, epsilon));
    assert(twisted.negate(twisted));

    // Witness identity r^ell = prod u^e * zeta^t for e = t = 1, ell = 2:
    // r = eps, u = twisted = -eps^2, zeta = -1 (the order's torsion
    // generator).  u * zeta = (-eps^2) * (-1) = eps^2 = r^2.
    silex::Element witness_power(field);
    silex::Element witness_product(field);
    assert(witness_power.multiply(epsilon, epsilon));
    assert(witness_product.negate(twisted));
    assert(witness_power.equal(witness_product));

    silex::OrderUnitGroup start(order);
    assert(set_single_unit_group(start, order, twisted, embeddings));

    sflint::Fmpz ell;
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(aux_bound, 31));

    // Bounded pre-scan route of the public l-local proof.
    silex::OrderUnitGroup local(order);
    silex::ProofState status = silex::ProofState::not_checked;
    bool changed = false;
    assert(local.prove_local_saturated(
            status, changed, start, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 256));
    assert(status == silex::ProofState::unavailable);
    assert(changed);
    assert(first_free_generator_is_pm_unit_power_pm1(local, epsilon));

    silex::OrderUnitGroup local_again(order);
    status = silex::ProofState::not_checked;
    changed = true;
    assert(local_again.prove_local_saturated(
            status, changed, local, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(aux_bound), embeddings, 256));
    assert(status == silex::ProofState::verified);
    assert(!changed);

    // Index-bound proof through the pre-scan (aux_bound 31) and directly
    // through the reference-stable scan (aux_bound 2 <= ell).
    const slong aux_bounds[] = {31, 2};
    for (slong bound : aux_bounds) {
        assert(set_fmpz_si(aux_bound, bound));
        silex::OrderUnitGroup proved(order);
        status = silex::ProofState::not_checked;
        changed = false;
        assert(proved.prove_index_bound(status, changed, start, 1,
                                        sflint::FmpzConstRef(aux_bound), 4,
                                        embeddings, 256));
        assert(status == silex::ProofState::verified);
        assert(changed);
        assert(proved.certification_status() ==
               silex::CertificationMode::proven);
        assert(first_free_generator_is_pm_unit_power_pm1(proved, epsilon));
        // GP: \p 100; quadregulator(8)
        assert(regulator_contains(
                proved,
                "0.8813735870195430252326093249797923090281603282616354107532"
                "956086533771842220260878337068919102560"));
        sflint::Fmpz torsion_order;
        assert(proved.torsion_order(sflint::FmpzRef(torsion_order)));
        assert(sflint::fmpz_equal_si(torsion_order, 2));
    }

    return 0;
}

// Non-proof saturation on the same u = -(3 + 2 sqrt2) = -eps^2 at ell = 2.
// At the two primes above 7 (sqrt2 = +/-3 mod 7), -1 is not a square
// (7 = 3 mod 4), so u is not a square either and the free-exponent kernel
// is empty: a free-only kernel misses the root.  The torsion-extended
// kernel is spanned by (e, t) = (1, 1), whose candidate u * (-1) = eps^2 has
// the root eps.  GP 2.17.4 cross-check:
//     kronecker(-1, 7)                               \\ -1
//     [kronecker(-9, 7), kronecker(3, 7)]            \\ [-1, -1]: u mod P
//     bnfisunit(bnfinit(x^2 - 2), -(3 + 2*x))        \\ [2, 1]~
int test_saturate_torsion_twisted_non_proof() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element twisted(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(twisted.multiply(epsilon, epsilon));
    assert(twisted.negate(twisted));

    silex::OrderUnitGroup start(order);
    assert(set_single_unit_group(start, order, twisted, embeddings));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    sflint::FmpzMat free_kernel(0, 0);
    assert(silex::detail::residue_dlog_free_kernel(
            free_kernel, start, prime_span, sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(free_kernel) == 1);
    assert(sflint::fmpz_mat_nrows(free_kernel) == 0);

    sflint::FmpzMat kernel(0, 0);
    assert(start.residue_dlog_kernel(kernel, prime_span,
                                     sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(kernel) == 2);
    assert(sflint::fmpz_mat_nrows(kernel) == 1);
    assert(mat_entry_is_si(kernel, 0, 0, 1));
    assert(mat_entry_is_si(kernel, 0, 1, 1));

    bool changed = false;
    silex::OrderUnitGroup row_out(order);
    assert(row_out.saturate_row(changed, start, sflint::FmpzMatConstRef(kernel),
                                0, sflint::FmpzConstRef(ell), embeddings,
                                256));
    assert(changed);
    assert(first_free_generator_is_pm_unit_power_pm1(row_out, epsilon));

    // The untwisted row (1, 0) is u itself, which is not a square.
    sflint::FmpzMat untwisted(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(untwisted, 0, 0), 1);
    silex::OrderUnitGroup untwisted_out(order);
    changed = true;
    assert(untwisted_out.saturate_row(
            changed, start, sflint::FmpzMatConstRef(untwisted), 0,
            sflint::FmpzConstRef(ell), embeddings, 256));
    assert(!changed);
    assert(check_first_free_generator(untwisted_out, twisted));

    silex::OrderUnitGroup local(order);
    changed = false;
    assert(local.saturate_local_once(changed, start, prime_span,
                                     sflint::FmpzConstRef(ell), embeddings,
                                     256));
    assert(changed);
    assert(first_free_generator_is_pm_unit_power_pm1(local, epsilon));

    // The bounded pass selects the two primes above 7 first (3 and 5 are
    // inert), adjoins eps, and the next pass leaves it unchanged.
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 31));
    silex::OrderUnitGroup bounded(order);
    bool stable = false;
    changed = false;
    assert(bounded.saturate_bounded(changed, stable, start,
                                    sflint::FmpzConstRef(ell), 2,
                                    sflint::FmpzConstRef(aux_bound), 4,
                                    embeddings, 256));
    assert(changed);
    assert(stable);
    assert(first_free_generator_is_pm_unit_power_pm1(bounded, epsilon));
    assert(bounded.certification_status() ==
           silex::CertificationMode::unknown);

    return 0;
}

// Z[3 sqrt2] = Z[sqrt18], the equation order of x^2 - 18, has index 3 in
// Z[sqrt2].  Its unit group is <-1, 17 + 4 sqrt18> = <-1, (1 + sqrt2)^4>,
// of index 4 in Z[sqrt2]^x; reference GP 2.17.4 session:
//     quadunit(72)       \\ 17 + 4*w, w = sqrt18
//     quadregulator(72)  \\ 3.5254943480781721009...
// The square root 3 + 2 sqrt2 = 3 + (2/3) sqrt18 of the fundamental unit is
// not in the order and must not be adjoined.
int test_prove_index_bound_nonmaximal_quadratic_sqrt18() {
    sflint::FmpzPoly polynomial;
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 2, 1);
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 0, -18);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpzPolyConstRef(polynomial));
    silex::Order equation = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    // GP: \p 100; quadregulator(72)
    const char* const kLogFundamental =
            "3.525494348078172100930437299919169236112641313046541643013182"
            "434613508736888104351334827567641024171";

    silex::Element order_unit(field);  // 17 + 4 sqrt18 = (1 + sqrt2)^4
    assert(set_quadratic_coeffs(order_unit, 17, 1, 4, 1));
    silex::Element order_unit2(field);  // 577 + 136 sqrt18
    assert(order_unit2.multiply(order_unit, order_unit));
    // -(17 + 4 sqrt18) = -eps^4 is a square only up to torsion:
    // -(17 + 4 sqrt18) * (-1) = eps^4 = (3 + 2 sqrt2)^2, and 3 + 2 sqrt2 is
    // not in the order (see the comment above).  The torsion-extended
    // kernel's twisted candidate (e, t) = (1, 1) sees this same
    // out-of-order root and must also fail closed.
    silex::Element neg_order_unit(field);
    assert(neg_order_unit.negate(order_unit));

    // l = 2 adjoins 17 + 4 sqrt18 (a square root in the order) from its
    // square, then meets 3 + 2 sqrt2 and fails closed.
    const silex::Element* inputs[] = {&order_unit, &order_unit2,
                                      &neg_order_unit};
    for (const silex::Element* input : inputs) {
        silex::OrderUnitGroup start(equation);
        assert(set_single_unit_group(start, equation, *input, embeddings));
        assert_nonmaximal_proof_fails_closed(equation, start, order_unit,
                                             kLogFundamental, 200,
                                             embeddings);
    }
    silex::OrderUnitGroup fundamental_start(equation);
    assert(set_single_unit_group(fundamental_start, equation, order_unit,
                                 embeddings));
    assert_nonmaximal_proof_fails_closed(equation, fundamental_start,
                                         order_unit, kLogFundamental, 2,
                                         embeddings);
    silex::OrderUnitGroup twisted_start(equation);
    assert(set_single_unit_group(twisted_start, equation, neg_order_unit,
                                 embeddings));
    assert_nonmaximal_proof_fails_closed(equation, twisted_start, order_unit,
                                         kLogFundamental, 2, embeddings);
    return 0;
}

// Non-proof twisted row in the nonmaximal order Z[sqrt18] (index 3 in
// Z[sqrt2]).  u = -(17 + 4 sqrt18) = -eps^4 with eps = 1 + sqrt2.  At the
// two primes above 7 (sqrt18 = +/-2 mod 7), u reduces to 3 and 5 and -1 to
// 6, none of them squares mod 7, so the torsion-extended kernel is spanned
// by (e, t) = (1, 1).  Its candidate u * (-1) = 17 + 4 sqrt18 has the square
// root 3 + 2 sqrt2 = 3 + (2/3) sqrt18, which is not in the order, so every
// non-proof route must succeed without adjoining it.  GP 2.17.4 cross-check:
//     [kronecker(3, 7), kronecker(5, 7), kronecker(-1, 7)]  \\ [-1, -1, -1]
//     kronecker(18, 7)                                      \\ 1: 7 splits
int test_saturate_torsion_twisted_non_proof_nonmaximal_sqrt18() {
    sflint::FmpzPoly polynomial;
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 2, 1);
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 0, -18);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpzPolyConstRef(polynomial));
    silex::Order order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element order_unit(field);  // 17 + 4 sqrt18 = (1 + sqrt2)^4
    assert(set_quadratic_coeffs(order_unit, 17, 1, 4, 1));
    silex::Element twisted(field);
    assert(twisted.negate(order_unit));

    silex::OrderUnitGroup start(order);
    assert(set_single_unit_group(start, order, twisted, embeddings));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 7));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 2);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    sflint::FmpzMat kernel(0, 0);
    assert(start.residue_dlog_kernel(kernel, prime_span,
                                     sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(kernel) == 2);
    assert(sflint::fmpz_mat_nrows(kernel) == 1);
    assert(mat_entry_is_si(kernel, 0, 0, 1));
    assert(mat_entry_is_si(kernel, 0, 1, 1));
    sflint::FmpzMat proof_kernel(0, 0);
    assert(start.residue_dlog_proof_kernel(proof_kernel, prime_span,
                                           sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(kernel),
                                  sflint::FmpzMatConstRef(proof_kernel)));

    bool changed = true;
    silex::OrderUnitGroup row_out(order);
    assert(row_out.saturate_row(changed, start, sflint::FmpzMatConstRef(kernel),
                                0, sflint::FmpzConstRef(ell), embeddings,
                                256));
    assert(!changed);
    assert(check_first_free_generator(row_out, twisted));
    assert(check_first_free_generator_is_order_unit(row_out, order));

    silex::OrderUnitGroup local(order);
    changed = true;
    assert(local.saturate_local_once(changed, start, prime_span,
                                     sflint::FmpzConstRef(ell), embeddings,
                                     256));
    assert(!changed);
    assert(check_first_free_generator(local, twisted));

    // The bounded pass skips 3 (decompose_prime(order, 3, 1) finds no
    // degree-one prime of this order above 3) and 5 (inert), selects the
    // two primes above 7, and adjoins nothing.
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 31));
    silex::OrderUnitGroup bounded(order);
    bool stable = false;
    changed = true;
    assert(bounded.saturate_bounded(changed, stable, start,
                                    sflint::FmpzConstRef(ell), 2,
                                    sflint::FmpzConstRef(aux_bound), 4,
                                    embeddings, 256));
    assert(!changed);
    assert(stable);
    assert(check_first_free_generator(bounded, twisted));
    assert(bounded.certification_status() ==
           silex::CertificationMode::unknown);

    return 0;
}

// True when the first free generator f of `group` satisfies
// f^exponent == unit^(+/-exponent), i.e. f = unit^(+/-1) up to a root of
// unity whose order divides `exponent`.
bool first_free_generator_matches_unit_mod_torsion(
        const silex::OrderUnitGroup& group,
        const silex::Element& unit,
        slong exponent) noexcept {
    const silex::NumberField* field = unit.parent();
    if (field == nullptr) {
        return false;
    }

    sflint::Fmpz power;
    sflint::fmpz_set_si(sflint::FmpzRef(power), exponent);
    silex::FactoredElement compact(*field);
    silex::Element expanded(*field);
    silex::Element generator_power(*field);
    silex::Element unit_power(*field);
    silex::Element inverse_power(*field);
    return group.free_generator(compact, 0) &&
           compact.evaluate(expanded) &&
           generator_power.pow_fmpz(expanded, sflint::FmpzConstRef(power)) &&
           unit_power.pow_fmpz(unit, sflint::FmpzConstRef(power)) &&
           inverse_power.invert(unit_power) &&
           (generator_power.equal(unit_power) ||
            generator_power.equal(inverse_power));
}

// Odd ell dividing the torsion order: ell = 3 in the maximal order
// Z[zeta_12] of x^4 - x^2 + 1 (w = 12, free rank 1, fundamental unit
// x - 1).  u = (1 - x)^3 * zeta_3^(-1) with zeta_3 = x^4, so u is a cube
// only up to torsion: u * x^k is a cube exactly when k = 1 mod 3.  At the
// four degree-one primes above 13, zeta_3 is not a cube (F_13^x has no
// cube of order 3), so neither is u, and the free-exponent kernel is
// empty.  The torsion-extended kernel is one-dimensional, and its row
// yields a cube root (1 - x) * zeta with zeta in mu_12.  GP 2.17.4
// cross-check (field variable y):
//     K = bnfinit(y^4 - y^2 + 1, 1); [K.tu, K.fu]  \\ [[12, y], [y - 1]]
//     u = Mod(1 - y, y^4 - y^2 + 1)^3 * Mod(y, y^4 - y^2 + 1)^8;
//     bnfisunit(K, lift(u))                         \\ [3, 2]~
//     #nfroots(K, x^3 - lift(u))                    \\ 0
//     #nfroots(K, x^3 - lift(u * Mod(y, y^4 - y^2 + 1)))  \\ 3
int test_saturate_torsion_twisted_non_proof_odd_ell() {
    sflint::FmpzPoly polynomial;
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 4, 1);
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 2, -1);
    ::fmpz_poly_set_coeff_si(polynomial.raw(), 0, 1);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpzPolyConstRef(polynomial));
    silex::Order order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element theta(field);
    silex::Element eta(field);  // 1 - x
    assert(theta.gen());
    assert(eta.negate(theta));
    assert(eta.add_si(eta, 1));
    sflint::Fmpz exponent;
    silex::Element eta_cube(field);
    silex::Element zeta3_inverse(field);  // x^8 = zeta_3^2
    silex::Element twisted(field);
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), 3);
    assert(eta_cube.pow_fmpz(eta, sflint::FmpzConstRef(exponent)));
    sflint::fmpz_set_si(sflint::FmpzRef(exponent), 8);
    assert(zeta3_inverse.pow_fmpz(theta, sflint::FmpzConstRef(exponent)));
    assert(twisted.multiply(eta_cube, zeta3_inverse));

    silex::OrderUnitGroup start(order);
    assert(set_single_unit_group(start, order, twisted, embeddings));
    sflint::Fmpz torsion_order;
    assert(start.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 12));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 13));
    assert(set_fmpz_si(ell, 3));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 4);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    sflint::FmpzMat free_kernel(0, 0);
    assert(silex::detail::residue_dlog_free_kernel(
            free_kernel, start, prime_span, sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(free_kernel) == 1);
    assert(sflint::fmpz_mat_nrows(free_kernel) == 0);

    // The torsion entry depends on which primitive 12th root is the torsion
    // generator, so only the shape and the nonzero entries are fixed.
    sflint::FmpzMat kernel(0, 0);
    assert(start.residue_dlog_kernel(kernel, prime_span,
                                     sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_ncols(kernel) == 2);
    assert(sflint::fmpz_mat_nrows(kernel) == 1);
    assert(!sflint::fmpz_is_zero(sflint::fmpz_mat_entry(kernel, 0, 0)));
    assert(!sflint::fmpz_is_zero(sflint::fmpz_mat_entry(kernel, 0, 1)));
    sflint::FmpzMat proof_kernel(0, 0);
    assert(start.residue_dlog_proof_kernel(proof_kernel, prime_span,
                                           sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_equal(sflint::FmpzMatConstRef(kernel),
                                  sflint::FmpzMatConstRef(proof_kernel)));

    bool changed = false;
    silex::OrderUnitGroup row_out(order);
    assert(row_out.saturate_row(changed, start, sflint::FmpzMatConstRef(kernel),
                                0, sflint::FmpzConstRef(ell), embeddings,
                                256));
    assert(changed);
    assert(first_free_generator_matches_unit_mod_torsion(row_out, eta, 12));

    // The untwisted row (1, 0) is u itself, which is not a cube.
    sflint::FmpzMat untwisted(1, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(untwisted, 0, 0), 1);
    silex::OrderUnitGroup untwisted_out(order);
    changed = true;
    assert(untwisted_out.saturate_row(
            changed, start, sflint::FmpzMatConstRef(untwisted), 0,
            sflint::FmpzConstRef(ell), embeddings, 256));
    assert(!changed);
    assert(check_first_free_generator(untwisted_out, twisted));

    silex::OrderUnitGroup local(order);
    changed = false;
    assert(local.saturate_local_once(changed, start, prime_span,
                                     sflint::FmpzConstRef(ell), embeddings,
                                     256));
    assert(changed);
    assert(first_free_generator_matches_unit_mod_torsion(local, eta, 12));

    // The bounded pass skips 7 (its primes have degree two) and selects the
    // four primes above 13, adjoins the cube root, and the next pass leaves
    // it unchanged.
    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 61));
    silex::OrderUnitGroup bounded(order);
    bool stable = false;
    changed = false;
    assert(bounded.saturate_bounded(changed, stable, start,
                                    sflint::FmpzConstRef(ell), 4,
                                    sflint::FmpzConstRef(aux_bound), 4,
                                    embeddings, 256));
    assert(changed);
    assert(stable);
    assert(first_free_generator_matches_unit_mod_torsion(bounded, eta, 12));
    assert(bounded.certification_status() ==
           silex::CertificationMode::unknown);

    return 0;
}

int test_saturate_index_bounded() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz aux_bound;
    assert(set_fmpz_si(aux_bound, 5));
    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    bool stable = false;
    assert(copied.saturate_index_bounded(
            changed, stable, rank_zero, degree_one_embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 2, 128));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(square)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::OrderUnitGroup saturated(order);
    changed = false;
    stable = true;
    assert(saturated.saturate_index_bounded(
            changed, stable, group, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 3, 128));
    // The square subgroup has index bound floor(2 log(1 + sqrt2) / 0.48) = 3
    // (Friedman 1989 Table 6, signature (2, 0)).  The 2-local pass adjoins
    // the square root and the 3-local pass is stable, so the result is the
    // full unit group, whose index bound 1 certifies it.
    assert(changed);
    assert(stable);
    assert(check_first_free_generator_square(saturated, epsilon2));
    assert(saturated.certification_status() ==
           silex::CertificationMode::proven);

    silex::OrderUnitGroup zero_pass(order);
    changed = true;
    stable = true;
    assert(zero_pass.saturate_index_bounded(
            changed, stable, group, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 0, 128));
    assert(!changed);
    assert(!stable);
    assert(check_first_free_generator(zero_pass, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup primitive_group(order);
    assert(primitive_group.set_units(order,
                                     silex::FactoredElementSpan(
                                             primitive_generators, 1),
                                     embeddings, 128));
    silex::Element epsilon3(field);
    assert(epsilon3.multiply(epsilon2, epsilon));
    silex::FactoredElement cube(field);
    assert(cube.set_element(epsilon3));
    silex::FactoredElement cube_generators[] = {std::move(cube)};
    silex::OrderUnitGroup no_progress(order);
    assert(no_progress.set_units(order,
                                 silex::FactoredElementSpan(cube_generators,
                                                            1),
                                 embeddings, 128));
    silex::OrderUnitGroup preserved(order);
    assert(preserved.set(primitive_group));
    // The full group has index bound floor(log(1 + sqrt2) / 0.48) = 1, so
    // no ell-local pass runs and the group is certified unchanged.
    silex::OrderUnitGroup certified(order);
    changed = true;
    stable = false;
    assert(certified.saturate_index_bounded(
            changed, stable, primitive_group, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 3, 128));
    assert(!changed);
    assert(stable);
    assert(check_first_free_generator(certified, epsilon));
    assert(certified.certification_status() ==
           silex::CertificationMode::proven);

    // The cube subgroup has index bound floor(3 log(1 + sqrt2) / 0.48) = 5.
    // With auxiliary bound 2 the 2-local pass finds no prime and nothing has
    // changed, so the call fails and preserves the destination.
    assert(set_fmpz_si(aux_bound, 2));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(aux_bound), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator(preserved, epsilon));

    sflint::Fmpz bad_bound;
    assert(set_fmpz_si(bad_bound, 1));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(bad_bound), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator(preserved, epsilon));

    return 0;
}

int test_saturate_index_bounded_adaptive() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz start;
    sflint::Fmpz max;
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 8));
    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    bool stable = false;
    assert(copied.saturate_index_bounded_adaptive(
            changed, stable, rank_zero, degree_one_embeddings, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 2, 128));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement square(field);
    assert(square.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(square)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 5));
    silex::OrderUnitGroup adaptive(order);
    changed = false;
    stable = true;
    assert(adaptive.saturate_index_bounded_adaptive(
            changed, stable, group, embeddings, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 3, 128));
    // Index bound 3 for the square subgroup (see
    // test_saturate_index_bounded): the first round adjoins the square root
    // and is stable, and the full group is certified.
    assert(changed);
    assert(stable);
    assert(check_first_free_generator_square(adaptive, epsilon2));
    assert(adaptive.certification_status() ==
           silex::CertificationMode::proven);

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup primitive_group(order);
    assert(primitive_group.set_units(order,
                                     silex::FactoredElementSpan(
                                             primitive_generators, 1),
                                     embeddings, 128));
    silex::Element epsilon3(field);
    assert(epsilon3.multiply(epsilon2, epsilon));
    silex::FactoredElement cube(field);
    assert(cube.set_element(epsilon3));
    silex::FactoredElement cube_generators[] = {std::move(cube)};
    silex::OrderUnitGroup no_progress(order);
    assert(no_progress.set_units(order,
                                 silex::FactoredElementSpan(cube_generators,
                                                            1),
                                 embeddings, 128));
    silex::OrderUnitGroup preserved(order);
    assert(preserved.set(primitive_group));
    // The cube subgroup has index bound floor(3 log(1 + sqrt2) / 0.48) = 5;
    // with auxiliary bounds capped at 2 no pass makes progress.
    assert(set_fmpz_si(start, 2));
    assert(set_fmpz_si(max, 2));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded_adaptive(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(start), sflint::FmpzConstRef(max), 3, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator(preserved, epsilon));

    sflint::Fmpz bad_start;
    assert(set_fmpz_si(bad_start, 1));
    changed = true;
    stable = false;
    assert(!preserved.saturate_index_bounded_adaptive(
            changed, stable, no_progress, embeddings, 1,
            sflint::FmpzConstRef(bad_start), sflint::FmpzConstRef(max), 3,
            128));
    assert(changed);
    assert(!stable);

    return 0;
}

int test_saturate_local_once_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 5));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, order, sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    assert(saturated.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(changed);
    assert(saturated.free_rank() == 1);
    assert(check_first_free_generator_square(saturated, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup already_saturated(order);
    assert(already_saturated.set_units(
            order, silex::FactoredElementSpan(primitive_generators, 1),
            embeddings, 128));

    silex::OrderUnitGroup copied(order);
    changed = true;
    assert(copied.saturate_local_once(
            changed, already_saturated, prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(!changed);
    assert(check_first_free_generator(copied, epsilon));

    return 0;
}

int test_saturate_local_once_rank_zero_and_failures() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz p;
    sflint::Fmpz ell;
    assert(set_fmpz_si(p, 5));
    assert(set_fmpz_si(ell, 2));
    silex::PrimeIdealList primes;
    assert(silex::decompose_prime(primes, degree_one_order,
                                  sflint::FmpzConstRef(p)));
    assert(primes.size() == 1);
    silex::PrimeIdealSpan prime_span(primes.at(0), primes.size());

    // ell = 2 divides w = 2, so the only column is the torsion generator
    // -1, a square mod 5 (-1 = 2^2): the kernel is spanned by (1).  Its
    // candidate -1 is not a square in Q, so nothing is adjoined below.
    sflint::FmpzMat kernel(5, 5);
    assert(rank_zero.residue_dlog_kernel(kernel, prime_span,
                                         sflint::FmpzConstRef(ell)));
    assert(sflint::fmpz_mat_nrows(kernel) == 1);
    assert(sflint::fmpz_mat_ncols(kernel) == 1);
    assert(mat_entry_is_si(kernel, 0, 0, 1));

    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    assert(copied.saturate_local_once(
            changed, rank_zero, prime_span, sflint::FmpzConstRef(ell),
            degree_one_embeddings, 80));
    assert(!changed);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::NumberField other_field = quadratic_field(3);
    silex::Order order;
    silex::Order other_order;
    order = silex::test::equation_order(field);
    other_order = silex::test::equation_order(other_field);
    silex::EmbeddingContext embeddings(field);
    silex::EmbeddingContext other_embeddings(other_field);
    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::PrimeIdealList other_primes;
    assert(silex::decompose_prime(other_primes, other_order,
                                  sflint::FmpzConstRef(p)));
    assert(other_primes.size() >= 1);
    silex::PrimeIdealSpan other_prime_span(other_primes.at(0), 1);

    sflint::Fmpz not_prime;
    assert(set_fmpz_si(not_prime, 4));
    assert(!copied.saturate_local_once(
            changed, group, silex::PrimeIdealSpan(),
            sflint::FmpzConstRef(ell), embeddings, 128));
    assert(!copied.saturate_local_once(
            changed, group, other_prime_span, sflint::FmpzConstRef(ell),
            embeddings, 128));
    assert(!copied.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(not_prime),
            embeddings, 128));
    assert(!copied.saturate_local_once(
            changed, group, prime_span, sflint::FmpzConstRef(ell),
            other_embeddings, 128));
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    return 0;
}

int test_select_saturation_primes() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    assert(set_real_quadratic_unit(epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 5));

    silex::PrimeIdealList selected;
    assert(group.select_saturation_primes(
            selected, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound)));
    assert(selected.size() == 1);
    assert(selected.at(0) != nullptr);
    assert(silex::same_order_parent(selected.at(0)->parent(), &order));

    assert(group.select_saturation_primes(
            selected, sflint::FmpzConstRef(ell), 0,
            sflint::FmpzConstRef(bound)));
    assert(selected.size() == 0);

    silex::PrimeIdealList preserved;
    assert(group.select_saturation_primes(
            preserved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound)));
    assert(set_fmpz_si(bound, 2));
    assert(!group.select_saturation_primes(
            preserved, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound)));
    assert(preserved.size() == 1);

    return 0;
}

// The public selector and the bounded saturation passes' selector run the
// same scans, so they pick the same primes.  Z[sqrt7] (disc 28) at ell = 3
// and Z[sqrt13] (disc 52) at ell = 3 each have a rational q = 1 mod ell
// dividing the discriminant (7 and 13, both ramified) with a degree-one
// prime above it: both selectors skip it, as the reference
// compute_candidates_for_saturate does.  Z[sqrt3] (disc 12) at ell = 2 is
// the ell = 2 skip case: q = 3 is odd, ramified, and its degree-one prime
// has an accepted character column, so removing the skip from either
// selector changes the selected primes.  Z[sqrt2] at ell = 2 and the cubic
// x^3 - 3x + 1 at ell = 2 are agreement-only cases (skipped_q = 0): Z[sqrt2]
// has no such q, and the cubic's q = 3 is not selected even without the
// skip.
void assert_saturation_selectors_agree_for(const silex::Order& order,
                                           const silex::Element* units,
                                           slong unit_count,
                                           slong ell_value,
                                           slong target_len,
                                           slong bound_value,
                                           slong skipped_q) {
    silex::EmbeddingContext embeddings(*order.parent());
    silex::OrderUnitGroup group(order);
    {
        std::vector<silex::FactoredElement> generators;
        for (slong i = 0; i < unit_count; ++i) {
            silex::FactoredElement generator(*order.parent());
            assert(generator.set_element(units[i]));
            generators.push_back(std::move(generator));
        }
        assert(group.set_units(
                order,
                silex::FactoredElementSpan(generators.data(), unit_count),
                embeddings, 256));
    }

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, ell_value));
    assert(set_fmpz_si(bound, bound_value));

    silex::PrimeIdealList public_primes;
    silex::PrimeIdealList bounded_primes;
    sflint::FmpzMat kernel(0, 0);
    assert(group.select_saturation_primes(public_primes,
                                          sflint::FmpzConstRef(ell),
                                          target_len,
                                          sflint::FmpzConstRef(bound)));
    assert(silex::detail::select_saturation_primes_with_kernel(
            bounded_primes, kernel, group, sflint::FmpzConstRef(ell),
            target_len, sflint::FmpzConstRef(bound)));
    assert(public_primes.size() == target_len);
    assert(bounded_primes.size() == target_len);
    for (slong i = 0; i < target_len; ++i) {
        const silex::PrimeIdeal* lhs = public_primes.at(i);
        const silex::PrimeIdeal* rhs = bounded_primes.at(i);
        assert(lhs != nullptr && rhs != nullptr);
        assert(lhs->equal(*rhs));
        sflint::Fmpz q;
        assert(lhs->rational_prime(sflint::FmpzRef(q)));
        assert(::fmpz_cmp_si(q.raw(), skipped_q) != 0);
    }
}

void assert_saturation_selectors_agree(slong radicand,
                                       const char* unit_constant,
                                       const char* unit_linear,
                                       slong ell_value,
                                       slong target_len,
                                       slong bound_value,
                                       slong skipped_q) {
    silex::NumberField field = quadratic_field(radicand);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::Element unit(field);
    assert(set_integral_quadratic_coeffs_str(unit, unit_constant,
                                             unit_linear));
    assert_saturation_selectors_agree_for(order, &unit, 1, ell_value,
                                          target_len, bound_value, skipped_q);
}

int test_saturation_selectors_agree() {
    assert_saturation_selectors_agree(7, "8", "3", 3, 2, 100, 7);
    assert_saturation_selectors_agree(13, "649", "180", 3, 2, 100, 13);
    // ell = 2 skip case: 2 + sqrt3 is the fundamental unit, and 3 | 12.
    assert_saturation_selectors_agree(3, "2", "1", 2, 2, 200, 3);
    // Agreement only: no skipped prime (skipped_q = 0).
    assert_saturation_selectors_agree(2, "1", "1", 2, 3, 100, 0);

    // Agreement only (skipped_q = 0).  The root theta of x^3 - 3x + 1 and
    // its conjugate theta^2 - 2 are independent units (rank 2).
    silex::NumberField cubic = cubic_field(-3, 1);
    silex::Order cubic_order;
    cubic_order = silex::test::equation_order(cubic);
    silex::Element theta(cubic);
    assert(theta.gen());
    silex::Element conjugate(cubic);
    assert(conjugate.multiply(theta, theta));
    assert(conjugate.add_si(conjugate, -2));
    silex::Element cubic_units[] = {std::move(theta), std::move(conjugate)};
    assert_saturation_selectors_agree_for(cubic_order, cubic_units, 2, 2, 2,
                                          200, 0);
    return 0;
}

// Forces the source_shape_fallback pass of both saturation selectors.  In
// Z[sqrt2] at ell = 3 the first q = 1 mod 3 is 7 (split, degree one).  With
// bound 5 the first pass finds no candidate, so each selector falls back to
// every prime ideal above every p <= 5: 2 is ramified (norm 2) and 3 inert
// (norm 9), neither with 3 | N - 1, while 5 is inert with norm 25, so the
// fallback selects the degree-two prime above 5.  With bound 7 the first
// pass selects the prime above 7 and the fallback does not run.
int test_saturation_selectors_source_shape_fallback() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element unit(field);
    assert(set_integral_quadratic_coeffs_str(unit, "1", "1"));
    silex::FactoredElement generator(field);
    assert(generator.set_element(unit));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 256));

    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    TorsionEventCounter counter;
    counter.label = "unit_group.saturation_selector.source_shape_fallback";
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::unit_group),
            torsion_event_callback, &counter);
    group.set_diagnostics(&diagnostics);
    const auto expect_fallback_events = [&](slong expected) {
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
        assert(counter.count == expected);
#else
        (void)expected;
        assert(counter.count == 0);
#endif
    };
    const auto selected_prime_is = [](const silex::PrimeIdealList& primes,
                                      slong q_value, slong degree) {
        sflint::Fmpz q;
        return primes.size() == 1 && primes.at(0) != nullptr &&
               primes.at(0)->rational_prime(sflint::FmpzRef(q)) &&
               sflint::fmpz_equal_si(q, q_value) &&
               primes.at(0)->residue_degree() == degree;
    };

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 3));
    for (slong bound_value : {5, 7}) {
        const bool fallback = bound_value == 5;
        const slong q_value = fallback ? 5 : 7;
        const slong degree = fallback ? 2 : 1;
        assert(set_fmpz_si(bound, bound_value));

        counter.count = 0;
        silex::PrimeIdealList public_primes;
        assert(group.select_saturation_primes(public_primes,
                                              sflint::FmpzConstRef(ell), 1,
                                              sflint::FmpzConstRef(bound)));
        assert(selected_prime_is(public_primes, q_value, degree));
        expect_fallback_events(fallback ? 1 : 0);

        counter.count = 0;
        silex::PrimeIdealList kernel_primes;
        sflint::FmpzMat kernel(0, 0);
        assert(silex::detail::select_saturation_primes_with_kernel(
                kernel_primes, kernel, group, sflint::FmpzConstRef(ell), 1,
                sflint::FmpzConstRef(bound)));
        assert(selected_prime_is(kernel_primes, q_value, degree));
        assert(kernel_primes.at(0)->equal(*public_primes.at(0)));
        expect_fallback_events(fallback ? 1 : 0);
    }

    group.set_diagnostics(nullptr);
    return 0;
}

int test_saturate_bounded_real_quadratic() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon2));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 5));

    silex::OrderUnitGroup saturated(order);
    bool changed = false;
    bool stable = false;
    assert(saturated.saturate_bounded(
            changed, stable, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, embeddings, 128));
    assert(changed);
    assert(stable);
    assert(saturated.free_rank() == 1);
    assert(check_first_free_generator_square(saturated, epsilon2));

    silex::FactoredElement primitive(field);
    assert(primitive.set_element(epsilon));
    silex::FactoredElement primitive_generators[] = {std::move(primitive)};
    silex::OrderUnitGroup already(order);
    assert(already.set_units(order,
                             silex::FactoredElementSpan(primitive_generators,
                                                        1),
                             embeddings, 128));

    silex::OrderUnitGroup unchanged(order);
    changed = true;
    stable = false;
    assert(unchanged.saturate_bounded(
            changed, stable, already, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, embeddings, 128));
    assert(!changed);
    assert(stable);
    assert(check_first_free_generator(unchanged, epsilon));

    return 0;
}

int test_saturate_bounded_rank_zero_zero_pass_partial_and_failures() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);
    silex::EmbeddingContext degree_one_embeddings(degree_one);
    silex::OrderUnitGroup rank_zero;
    assert(rank_zero.compute(degree_one_order));

    sflint::Fmpz ell;
    sflint::Fmpz bound;
    assert(set_fmpz_si(ell, 2));
    assert(set_fmpz_si(bound, 5));

    silex::OrderUnitGroup copied(degree_one_order);
    bool changed = true;
    bool stable = false;
    assert(copied.saturate_bounded(
            changed, stable, rank_zero, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, degree_one_embeddings, 80));
    assert(!changed);
    assert(stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    changed = true;
    stable = true;
    assert(copied.saturate_bounded(
            changed, stable, rank_zero, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 0, degree_one_embeddings, 80));
    assert(!changed);
    assert(!stable);
    assert(check_rank_zero_group(copied, degree_one_order, 2));

    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::Element epsilon(field);
    silex::Element epsilon2(field);
    silex::Element epsilon4(field);
    assert(set_real_quadratic_unit(epsilon));
    assert(epsilon2.multiply(epsilon, epsilon));
    assert(epsilon4.multiply(epsilon2, epsilon2));
    silex::FactoredElement generator(field);
    assert(generator.set_element(epsilon4));
    silex::FactoredElement generators[] = {std::move(generator)};
    silex::OrderUnitGroup group(order);
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));

    silex::OrderUnitGroup partial(order);
    changed = false;
    stable = true;
    assert(partial.saturate_bounded(
            changed, stable, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 1, embeddings, 128));
    assert(changed);
    assert(!stable);
    assert(check_first_free_generator_square(partial, epsilon4));

    silex::OrderUnitGroup preserved(degree_one_order);
    assert(preserved.set(rank_zero));
    assert(set_fmpz_si(bound, 2));
    changed = false;
    stable = false;
    assert(!preserved.saturate_bounded(
            changed, stable, group, sflint::FmpzConstRef(ell), 1,
            sflint::FmpzConstRef(bound), 3, embeddings, 128));
    assert(check_rank_zero_group(preserved, degree_one_order, 2));

    return 0;
}

int test_supplied_real_quadratic_units() {
    silex::NumberField field = quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::Element unit(field);
    assert(set_real_quadratic_unit(unit));
    silex::FactoredElement generator(field);
    assert(generator.set_element(unit));
    silex::FactoredElement generators[] = {std::move(generator)};

    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(generators, 1),
                           embeddings, 128));
    assert(group.is_defined());
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &order));
    assert(group.free_rank() == 1);
    assert(group.certification_status() == silex::CertificationMode::unknown);

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    assert(group.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_positive(regulator));
    assert(check_first_free_generator(group, unit));

    silex::OrderUnitGroup copy;
    assert(copy.set(group));
    assert(copy.free_rank() == 1);
    assert(copy.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_positive(regulator));
    assert(check_first_free_generator(copy, unit));

    silex::OrderUnitGroup moved(std::move(copy));
    assert(moved.free_rank() == 1);
    assert(check_first_free_generator(moved, unit));
    assert(!copy.is_defined());

    return 0;
}

int test_supplied_rank_zero_units() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);

    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(), embeddings,
                           80));
    assert(group.is_defined());
    assert(group.is_set());
    assert(silex::same_order_parent(group.parent(), &order));
    assert(group.free_rank() == 0);
    assert(group.certification_status() == silex::CertificationMode::unknown);

    sflint::Fmpz torsion_order;
    sflint::Arb regulator;
    assert(group.torsion_order(sflint::FmpzRef(torsion_order)));
    assert(sflint::fmpz_equal_si(torsion_order, 2));
    assert(group.regulator(sflint::ArbRef(regulator)));
    assert(sflint::arb_is_one(regulator));

    return 0;
}

int test_supplied_unit_failures_preserve_output() {
    silex::NumberField field = quadratic_field(2);
    silex::NumberField other_field = quadratic_field(3);
    silex::Order order;
    order = silex::test::equation_order(field);
    silex::EmbeddingContext embeddings(field);
    silex::EmbeddingContext other_embeddings(other_field);

    silex::Element unit(field);
    assert(set_real_quadratic_unit(unit));
    silex::FactoredElement compact_unit(field);
    assert(compact_unit.set_element(unit));
    silex::FactoredElement units[] = {std::move(compact_unit)};

    silex::OrderUnitGroup group;
    assert(group.set_units(order, silex::FactoredElementSpan(units, 1),
                           embeddings, 128));
    sflint::Arb regulator;
    assert(group.regulator(sflint::ArbRef(regulator)));

    silex::Element two(field);
    assert(two.set_si(2));
    silex::FactoredElement nonunit(field);
    assert(nonunit.set_element(two));
    silex::FactoredElement nonunits[] = {std::move(nonunit)};
    assert(!group.set_units(order, silex::FactoredElementSpan(nonunits, 1),
                            embeddings, 128));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator(group, unit));

    silex::Element minus_one(field);
    assert(minus_one.set_si(-1));
    silex::FactoredElement torsion(field);
    assert(torsion.set_element(minus_one));
    silex::FactoredElement torsion_generators[] = {std::move(torsion)};
    assert(!group.set_units(
            order, silex::FactoredElementSpan(torsion_generators, 1),
            embeddings, 128));
    assert(group.free_rank() == 1);
    assert(check_first_free_generator(group, unit));

    assert(!group.set_units(order, silex::FactoredElementSpan(), embeddings,
                            128));
    assert(!group.set_units(order, silex::FactoredElementSpan(units, 1),
                            embeddings, 0));
    assert(!group.set_units(order, silex::FactoredElementSpan(units, 1),
                            other_embeddings, 128));

    silex::Element foreign_unit(other_field);
    assert(set_real_quadratic_unit(foreign_unit));
    silex::FactoredElement foreign(other_field);
    assert(foreign.set_element(foreign_unit));
    silex::FactoredElement foreign_generators[] = {std::move(foreign)};
    assert(!group.set_units(
            order, silex::FactoredElementSpan(foreign_generators, 1),
            embeddings, 128));
    assert(group.is_set());
    assert(group.free_rank() == 1);
    assert(check_first_free_generator(group, unit));

    return 0;
}

int test_imaginary_quadratic_torsion() {
    silex::NumberField qi = quadratic_field(-1);
    silex::Order equation_order;
    equation_order = silex::test::equation_order(qi);

    silex::OrderUnitGroup group;
    assert(group.compute(equation_order));
    assert(check_rank_zero_group(group, equation_order, 4));

    return 0;
}

int test_nonmaximal_torsion_filtering() {
    silex::NumberField qsqrt_minus3 = quadratic_field(-3);
    silex::Order equation_order;
    silex::Order maximal_order(qsqrt_minus3);
    equation_order = silex::test::equation_order(qsqrt_minus3);
    assert(maximal_order.maximal_order(equation_order));

    silex::OrderUnitGroup nonmaximal;
    assert(nonmaximal.compute(equation_order));
    assert(check_rank_zero_group(nonmaximal, equation_order, 2));

    silex::OrderUnitGroup maximal;
    assert(maximal.compute(maximal_order));
    assert(check_rank_zero_group(maximal, maximal_order, 6));

    return 0;
}

int test_positive_rank_failure_preserves_output() {
    silex::NumberField degree_one = degree_one_field();
    silex::Order degree_one_order;
    degree_one_order = silex::test::equation_order(degree_one);

    silex::OrderUnitGroup group;
    assert(group.compute(degree_one_order));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    silex::NumberField real_quadratic = quadratic_field(5);
    silex::Order real_order;
    real_order = silex::test::equation_order(real_quadratic);
    assert(!real_order.is_maximal());
    assert(!group.compute(real_order));
    assert(check_rank_zero_group(group, degree_one_order, 2));

    return 0;
}

// A monic integral x^2+b*x+c whose equation order is not maximal: compute()
// publishes the continued-fraction unit for the maximal order only, and the
// equation order (which does not contain that unit) keeps today's failure.
void assert_nonmaximal_monic_quadratic_units(slong linear,
                                             slong constant,
                                             slong constant_num,
                                             ulong constant_den,
                                             slong linear_num,
                                             ulong linear_den,
                                             const char* regulator) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, linear);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, constant);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order equation_order = silex::test::equation_order(field);
    silex::Order maximal_order(field);
    assert(maximal_order.maximal_order(equation_order));
    assert(!equation_order.is_maximal());
    assert(maximal_order.is_maximal());

    silex::OrderUnitGroup nonmaximal;
    assert(!nonmaximal.compute(equation_order));
    assert(!nonmaximal.is_set());

    silex::OrderUnitGroup maximal;
    assert(maximal.compute(maximal_order));
    assert(maximal.free_rank() == 1);
    assert(maximal.certification_status() == silex::CertificationMode::proven);

    silex::FactoredElement generator(field);
    silex::Element value(field);
    silex::Element expected(field);
    sflint::FmpqPoly expected_polynomial;
    sflint::Fmpq coefficient;
    sflint::fmpq_set_si(coefficient, constant_num, constant_den);
    sflint::fmpq_poly_set_coeff_fmpq(expected_polynomial, 0, coefficient);
    sflint::fmpq_set_si(coefficient, linear_num, linear_den);
    sflint::fmpq_poly_set_coeff_fmpq(expected_polynomial, 1, coefficient);
    assert(expected.set_fmpq_poly(
            sflint::FmpqPolyConstRef(expected_polynomial)));
    assert(maximal.free_generator(generator, 0));
    assert(generator.evaluate(value));
    assert(value.equal(expected));

    sflint::Arb published;
    sflint::Arb reference;
    assert(maximal.regulator(sflint::ArbRef(published)));
    arb_set_str(reference.raw(), regulator, 128);
    assert(arb_overlaps(published.raw(), reference.raw()));
}

int test_nonmaximal_monic_quadratic_units() {
    // x^2 - 12: index 2 in Z[(sqrt 3)]; unit 2 + sqrt(3) = 2 + theta/2.
    assert_nonmaximal_monic_quadratic_units(
            0, -12, 2, 1, 1, 2, "1.316957896924816708625046347 +/- 1e-25");
    // x^2 + 3x - 9: D = 45, index 3; unit (1 + sqrt 5)/2 = 1 + theta/3.
    assert_nonmaximal_monic_quadratic_units(
            3, -9, 1, 1, 1, 3, "0.4812118250596034474977589134 +/- 1e-25");
    // x^2 + 2x - 44: D = 180, index 6; unit (1 + sqrt 5)/2 = (theta + 4)/6.
    assert_nonmaximal_monic_quadratic_units(
            2, -44, 2, 3, 1, 6, "0.4812118250596034474977589134 +/- 1e-25");
    return 0;
}

int test_group_keeps_parent_order_alive() {
    silex::OrderUnitGroup group = real_quadratic_unit_group();
    assert(group.is_defined());
    assert(group.is_set());
    assert(group.parent() != nullptr);
    assert(group.parent()->parent() != nullptr);
    assert(group.free_rank() == 1);

    silex::OrderElement torsion(*group.parent());
    assert(group.torsion_generator(torsion));

    silex::FactoredElement generator(*group.parent()->parent());
    silex::Element expanded(*group.parent()->parent());
    silex::OrderElement order_element(*group.parent());
    silex::Ideal principal(*group.parent());
    assert(group.free_generator(generator, 0));
    assert(generator.evaluate(expanded));
    assert(order_element.set_element(expanded));
    assert(principal.set_principal(order_element));
    assert(principal.is_one());

    return 0;
}

int test_set_move_and_access_failures() {
    silex::NumberField field = degree_one_field();
    silex::Order order;
    order = silex::test::equation_order(field);

    silex::OrderUnitGroup group;
    assert(group.compute(order));

    silex::OrderUnitGroup copy;
    assert(copy.set(group));
    assert(check_rank_zero_group(copy, order, 2));

    silex::OrderUnitGroup moved(std::move(copy));
    assert(check_rank_zero_group(moved, order, 2));
    assert(!copy.is_defined());

    sflint::Fmpz order_out;
    silex::NumberField other_field = quadratic_field(-1);
    silex::Order other_order;
    other_order = silex::test::equation_order(other_field);
    silex::OrderElement wrong_generator(other_order);
    assert(!moved.torsion_generator(wrong_generator));
    assert(moved.torsion_order(sflint::FmpzRef(order_out)));

    silex::FactoredElement free_generator(field);
    assert(!moved.free_generator(free_generator, 0));

    moved.clear();
    assert(!moved.is_defined());
    assert(moved.define(order));
    assert(moved.is_defined());
    assert(!moved.is_set());

    silex::NumberField sqrt2 = quadratic_field(2);
    silex::Order sqrt2_order;
    sqrt2_order = silex::test::equation_order(sqrt2);
    silex::Element epsilon(sqrt2);
    assert(set_quadratic_coeffs(epsilon, 1, 1, 1, 1));

    silex::OrderUnitGroup positive_rank;
    assert(positive_rank.compute(sqrt2_order));
    assert(check_real_quadratic_group(positive_rank, sqrt2_order, epsilon));

    silex::OrderUnitGroup positive_rank_copy;
    assert(positive_rank_copy.set(positive_rank));
    assert(check_real_quadratic_group(positive_rank_copy, sqrt2_order,
                                      epsilon));

    silex::OrderUnitGroup positive_rank_move;
    positive_rank_move = std::move(positive_rank_copy);
    assert(check_real_quadratic_group(positive_rank_move, sqrt2_order,
                                      epsilon));
    assert(!positive_rank_copy.is_defined());

    return 0;
}

// Wide-exponent unit-basis transform (T-157).  In the cyclic cubic
// x^3 + x^2 - 2x - 1 let v1 = theta and v2 = theta^2 - 2 (conjugate units).
// With w = v2, w2 = v2^2 (a distinct base) and K = 2k, k = 2^40, the stored
// units u1 = v1 w^(-2k) and u2 = v1^-1 w2^k w under T = [[K+1, K], [1, 1]]
// give u1^(K+1) u2^K = v1 w^(-4k^2) w2^(2k^2) = v1 and u1 u2 = v2, with
// merged exponents of size 2^82 that overflow slong.
struct WideTransformCase {
    silex::NumberField field;
    silex::Order order;
    silex::Element v1;
    silex::Element v2;
    std::vector<silex::FactoredElement> units;
    sflint::FmpzMat transform{2, 2};

    explicit WideTransformCase(bool with_w2)
        : field(make_wide_case_field()),
          order(field),
          v1(field),
          v2(field) {
        silex::Order equation = silex::test::equation_order(field);
        assert(order.maximal_order(equation));
        silex::Element theta(field);
        silex::Element two(field);
        silex::Element theta_squared(field);
        assert(theta.gen());
        assert(v1.set(theta));
        assert(theta_squared.multiply(theta, theta));
        assert(v2.add_si(theta_squared, -2));
        silex::Element w2(field);
        assert(w2.multiply(v2, v2));
        const slong k = slong(1) << 40;
        silex::FactoredElement u1(field);
        silex::FactoredElement u2(field);
        assert(u1.push(v1, 1));
        assert(u1.push(v2, -2 * k));
        assert(u2.push(v1, -1));
        if (with_w2) {
            assert(u2.push(w2, k));
            assert(u2.push(v2, 1));
        } else {
            assert(u2.push(v2, 2 * k + 1));
        }
        units.push_back(std::move(u1));
        units.push_back(std::move(u2));
        const slong big = 2 * k;
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(transform, 0, 0), big + 1);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(transform, 0, 1), big);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(transform, 1, 0), 1);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(transform, 1, 1), 1);
    }

    static silex::NumberField make_wide_case_field() {
        sflint::FmpqPoly polynomial;
        sflint::fmpq_poly_zero(polynomial);
        sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
        sflint::fmpq_poly_set_coeff_si(polynomial, 1, -2);
        sflint::fmpq_poly_set_coeff_si(polynomial, 0, -1);
        return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    }
};

bool factored_equals(const silex::FactoredElement& value,
                     const silex::Element& expected) {
    silex::Element evaluated(*value.parent());
    return value.evaluate(evaluated) && evaluated.equal(expected);
}

int test_transform_factored_units_wide_exponents() {
    WideTransformCase c(true);
    silex::EmbeddingContext embeddings(c.field);

    // The word-size transform overflows on these exponents and leaves the
    // units unchanged.
    std::vector<silex::FactoredElement> word_units;
    for (const auto& unit : c.units) {
        word_units.emplace_back(c.field);
        assert(word_units.back().set(unit));
    }
    bool overflowed = false;
    assert(!silex::detail::transform_factored_units(word_units, c.transform,
                                                    overflowed));
    assert(overflowed);

    silex::detail::WideFactoredUnits wide;
    bool has_wide = false;
    assert(silex::detail::transform_factored_units_wide(
            c.units, wide, has_wide, c.transform, c.order, embeddings,
            nullptr));
    assert(has_wide);
    assert(silex::detail::wide_factored_units_exceed_slong(wide));
    slong max_bits = 0;
    for (slong i = 0; i < sflint::fmpz_mat_nrows(wide.exponents); ++i) {
        for (slong j = 0; j < sflint::fmpz_mat_ncols(wide.exponents); ++j) {
            max_bits = FLINT_MAX(
                    max_bits,
                    slong(::fmpz_bits(
                            ::fmpz_mat_entry(wide.exponents.raw(), i, j))));
        }
    }
    assert(max_bits == 83);  // 4k^2 = 2^82
    assert(c.units.size() == 2);
    assert(factored_equals(c.units[0], c.v1));
    assert(factored_equals(c.units[1], c.v2));
    return 0;
}

// With w only (no w2) the merged exponents fit in slong: the products
// themselves are the units and nothing wide is stored.
int test_transform_factored_units_wide_merge_fits() {
    WideTransformCase c(false);
    silex::EmbeddingContext embeddings(c.field);
    silex::detail::WideFactoredUnits wide;
    bool has_wide = true;
    assert(silex::detail::transform_factored_units_wide(
            c.units, wide, has_wide, c.transform, c.order, embeddings,
            nullptr));
    assert(!has_wide);
    assert(c.units.size() == 2);
    assert(factored_equals(c.units[0], c.v1));
    assert(factored_equals(c.units[1], c.v2));
    return 0;
}

int test_group_wide_free_generators() {
    WideTransformCase c(true);
    silex::EmbeddingContext embeddings(c.field);
    silex::detail::WideFactoredUnits wide;
    bool has_wide = false;
    assert(silex::detail::transform_factored_units_wide(
            c.units, wide, has_wide, c.transform, c.order, embeddings,
            nullptr));
    assert(has_wide);

    silex::OrderUnitGroup group;
    assert(silex::detail::OrderUnitGroupAccess::set_units_with_wide_form(
            group, c.order,
            silex::FactoredElementSpan(c.units.data(), 2), wide, embeddings,
            256));
    assert(group.is_set());
    assert(group.free_rank() == 2);
    const silex::detail::WideFactoredUnits* stored =
            silex::detail::OrderUnitGroupAccess::wide_free_generators(group);
    assert(stored != nullptr);
    assert(silex::detail::wide_factored_units_exceed_slong(*stored));

    // A copy keeps the wide form.
    silex::OrderUnitGroup copy;
    assert(copy.set(group));
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(copy) !=
           nullptr);

    for (slong i = 0; i < 2; ++i) {
        silex::FactoredElement generator(c.field);
        assert(group.free_generator(generator, i));
        assert(generator.length() == 1);
        slong exponent = 0;
        assert(generator.exponent(exponent, 0));
        assert(exponent == 1);
        assert(factored_equals(generator, i == 0 ? c.v1 : c.v2));
    }
    return 0;
}

bool same_factored(const silex::FactoredElement& left,
                   const silex::FactoredElement& right) {
    if (left.length() != right.length()) {
        return false;
    }
    for (slong i = 0; i < left.length(); ++i) {
        const auto& a = left.factors()[static_cast<std::size_t>(i)];
        const auto& b = right.factors()[static_cast<std::size_t>(i)];
        if (a.exponent != b.exponent || !a.factor.equal(b.factor)) {
            return false;
        }
    }
    return true;
}

// Only an exponent overflow sends the transform to the wide path: invalid
// inputs fail with `overflowed` false.
int test_transform_factored_units_overflow_flag() {
    WideTransformCase c(true);
    bool overflowed = true;

    // Transform of the wrong size.
    sflint::FmpzMat wrong_size(1, 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(wrong_size, 0, 0), 1);
    assert(!silex::detail::transform_factored_units(c.units, wrong_size,
                                                    overflowed));
    assert(!overflowed);

    // A unit from another field: pow_si fails without an overflow.
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);
    silex::NumberField other = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Element root(other);
    assert(root.gen());
    std::vector<silex::FactoredElement> mixed;
    mixed.emplace_back(c.field);
    assert(mixed.back().push(c.v1, 1));
    mixed.emplace_back(other);
    assert(mixed.back().push(root, 1));
    sflint::FmpzMat swap(2, 2);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(swap, 0, 1), 1);
    sflint::fmpz_set_si(sflint::fmpz_mat_entry(swap, 1, 0), 1);
    overflowed = true;
    assert(!silex::detail::transform_factored_units(mixed, swap, overflowed));
    assert(!overflowed);
    return 0;
}

// A failing wide transform leaves units, wide and has_wide unchanged.  The
// equation order's maximality is not known, so the rebuild step rejects it.
int test_transform_factored_units_wide_failure() {
    WideTransformCase c(true);
    silex::Order equation = silex::test::equation_order(c.field);
    assert(!equation.maximality_known());
    silex::EmbeddingContext embeddings(c.field);
    std::vector<silex::FactoredElement> before;
    for (const auto& unit : c.units) {
        before.emplace_back(c.field);
        assert(before.back().set(unit));
    }
    silex::detail::WideFactoredUnits wide;
    bool has_wide = false;
    assert(!silex::detail::transform_factored_units_wide(
            c.units, wide, has_wide, c.transform, equation, embeddings,
            nullptr));
    assert(!has_wide);
    assert(wide.bases.empty());
    assert(sflint::fmpz_mat_nrows(wide.exponents) == 0);
    assert(c.units.size() == before.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
        assert(same_factored(c.units[i], before[i]));
    }
    return 0;
}

// set_units_with_wide_form fails on mismatched sizes and leaves the group
// unchanged; move and clear() carry or drop the wide form with the group.
int test_group_wide_form_failure_move_clear() {
    WideTransformCase c(true);
    silex::EmbeddingContext embeddings(c.field);
    silex::detail::WideFactoredUnits wide;
    bool has_wide = false;
    assert(silex::detail::transform_factored_units_wide(
            c.units, wide, has_wide, c.transform, c.order, embeddings,
            nullptr));
    assert(has_wide);
    silex::detail::WideFactoredUnits spare;
    assert(silex::detail::wide_factored_units_copy(spare, wide));

    silex::OrderUnitGroup group;
    assert(silex::detail::OrderUnitGroupAccess::set_units_with_wide_form(
            group, c.order,
            silex::FactoredElementSpan(c.units.data(), 2), wide, embeddings,
            256));
    const silex::detail::WideFactoredUnits* stored =
            silex::detail::OrderUnitGroupAccess::wide_free_generators(group);
    assert(stored != nullptr);

    // One rebuilt unit against a two-row wide form.
    assert(!silex::detail::OrderUnitGroupAccess::set_units_with_wide_form(
            group, c.order,
            silex::FactoredElementSpan(c.units.data(), 1), spare, embeddings,
            256));
    assert(group.is_set());
    assert(group.free_rank() == 2);
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(group) ==
           stored);
    // Two rebuilt units against a one-row wide form.
    silex::detail::WideFactoredUnits one_row;
    assert(silex::detail::wide_factored_units_copy(one_row, spare));
    sflint::FmpzMat short_exponents(
            1, sflint::fmpz_mat_ncols(spare.exponents));
    one_row.exponents.swap(short_exponents);
    assert(!silex::detail::OrderUnitGroupAccess::set_units_with_wide_form(
            group, c.order,
            silex::FactoredElementSpan(c.units.data(), 2), one_row,
            embeddings, 256));
    assert(group.is_set());
    assert(group.free_rank() == 2);
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(group) ==
           stored);
    for (slong i = 0; i < 2; ++i) {
        silex::FactoredElement generator(c.field);
        assert(group.free_generator(generator, i));
        assert(factored_equals(generator, i == 0 ? c.v1 : c.v2));
    }

    silex::OrderUnitGroup moved(std::move(group));
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(moved) ==
           stored);
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(group) ==
           nullptr);
    silex::OrderUnitGroup assigned;
    assigned = std::move(moved);
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(
                   assigned) == stored);
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(moved) ==
           nullptr);

    assigned.clear();
    assert(!assigned.is_set());
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(
                   assigned) == nullptr);
    return 0;
}

// reduce_stored_relation_units on stored units whose unit-basis transform
// overflows slong: u1 = v1 v2^M w2^(-M/2) = v1 with M = 2^62 + 2 and
// u2 = v1^3 v2^2.  Reduction combines u2 with a multiple of u1 whose merged
// exponents exceed slong, so the group keeps the wide form and its free
// generators are the rebuilt units.
int test_reduce_stored_units_keeps_wide_form() {
    WideTransformCase c(true);
    silex::EmbeddingContext embeddings(c.field);
    silex::Element w2(c.field);
    assert(w2.multiply(c.v2, c.v2));
    const slong m = (slong(1) << 62) + 2;
    std::vector<silex::FactoredElement> units;
    units.emplace_back(c.field);
    assert(units.back().push(c.v1, 1));
    assert(units.back().push(c.v2, m));
    assert(units.back().push(w2, -(m / 2)));
    units.emplace_back(c.field);
    assert(units.back().push(c.v1, 3));
    assert(units.back().push(c.v2, 2));

    silex::OrderUnitGroup group(c.order);
    assert(silex::detail::order_unit_group_set_units_internal(
            group, c.order, silex::FactoredElementSpan(units.data(), 2),
            embeddings, 256, true));
    assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(group) ==
           nullptr);

    bool changed = false;
    assert(silex::detail::reduce_stored_relation_units(changed, group,
                                                       embeddings, 256));
    assert(changed);
    assert(group.free_rank() == 2);
    const silex::detail::WideFactoredUnits* stored =
            silex::detail::OrderUnitGroupAccess::wide_free_generators(group);
    assert(stored != nullptr);
    assert(silex::detail::wide_factored_units_exceed_slong(*stored));
    for (slong i = 0; i < 2; ++i) {
        silex::FactoredElement generator(c.field);
        assert(group.free_generator(generator, i));
        assert(generator.length() == 1);
    }
    return 0;
}

// Q(zeta_7), f = x^6 + x^5 + ... + 1 (Z[zeta] maximal): signature (0, 3),
// unit rank 2, torsion of order 14.  The cyclotomic units u2 = 1 + zeta and
// u3 = 1 + zeta + zeta^2 are independent.
struct Cyclotomic7Units {
    silex::NumberField field;
    silex::Order order;
    silex::Element u2;
    silex::Element u3;

    Cyclotomic7Units()
        : field(make_field()),
          order(silex::test::equation_order(field)),
          u2(field),
          u3(field) {
        assert(set_partial_sum(u2, 2));
        assert(set_partial_sum(u3, 3));
    }

    static silex::NumberField make_field() {
        sflint::FmpqPoly polynomial;
        for (slong i = 0; i <= 6; ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial, i, 1);
        }
        return silex::test::field_by_polynomial(
                sflint::FmpqPolyConstRef(polynomial));
    }

    // 1 + zeta + ... + zeta^(a-1).
    static bool set_partial_sum(silex::Element& out, slong a) {
        sflint::FmpqPoly polynomial;
        for (slong i = 0; i < a; ++i) {
            sflint::fmpq_poly_set_coeff_si(polynomial, i, 1);
        }
        return out.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial));
    }

    // {u2^e2, u2^e21 u3}.
    std::vector<silex::FactoredElement> basis(slong e2, slong e21) const {
        std::vector<silex::FactoredElement> out;
        out.emplace_back(field);
        assert(out.back().push(u2, e2));
        out.emplace_back(field);
        if (e21 != 0) {
            assert(out.back().push(u2, e21));
        }
        assert(out.back().push(u3, 1));
        return out;
    }
};

bool same_regulator(const silex::OrderUnitGroup& left,
                    const silex::OrderUnitGroup& right) {
    const auto a = left.regulator();
    const auto b = right.regulator();
    return a.has_value() && b.has_value() &&
           ::arb_equal(a->raw(), b->raw()) != 0;
}

// The internal re-basis steps keep the torsion stored in the group instead
// of recomputing it.  The cached-torsion hook accepts the under-claimed pair
// (2, -1) in Q(zeta_7), whose torsion has order 14, so a recomputation would
// show as order 14.
int test_rebasis_keeps_group_torsion() {
    Cyclotomic7Units c;
    silex::EmbeddingContext embeddings(c.field);
    sflint::Fmpz two;
    sflint::fmpz_set_si(sflint::FmpzRef(two), 2);
    silex::OrderElement minus_one(c.order);
    assert(minus_one.set_si(-1));

    const auto torsion_order_is = [](const silex::OrderUnitGroup& group,
                                     slong expected) {
        const auto order = group.torsion_order();
        return order.has_value() &&
               sflint::fmpz_equal_si(sflint::FmpzConstRef(*order), expected);
    };

    // Reduction of a skewed basis.
    {
        const std::vector<silex::FactoredElement> skewed = c.basis(1, 5);
        const silex::FactoredElementSpan span(skewed.data(), skewed.size());
        silex::OrderUnitGroup computed(c.order);
        assert(silex::detail::order_unit_group_set_units_internal(
                computed, c.order, span, embeddings, 256, true));
        assert(torsion_order_is(computed, 14));

        silex::OrderUnitGroup group(c.order);
        assert(silex::detail::order_unit_group_set_units_internal(
                group, c.order, span, embeddings, 256, true, &two,
                &minus_one));
        assert(torsion_order_is(group, 2));
        bool changed = false;
        assert(silex::detail::reduce_stored_relation_units(changed, group,
                                                           embeddings, 256));
        assert(changed);
        assert(torsion_order_is(group, 2));
        // The propagated under-claimed torsion is rejected by the guard
        // that proven publication relies on; computed torsion passes.
        assert(silex::detail::order_unit_torsion_is_computed(computed));
        assert(!silex::detail::order_unit_torsion_is_computed(group));
    }

    // Adjoining u2 = sqrt(u2^2) to <u2^2, u3>: the relation row
    // (-1, 0, 2) says root^2 = u2^2.
    {
        const std::vector<silex::FactoredElement> square = c.basis(2, 0);
        const silex::FactoredElementSpan span(square.data(), square.size());
        silex::FactoredElement root(c.field);
        assert(root.set_element(c.u2));
        sflint::FmpzMat rel(1, 3);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(rel, 0, 0), -1);
        sflint::fmpz_set_si(sflint::fmpz_mat_entry(rel, 0, 2), 2);

        silex::OrderUnitGroup computed(c.order);
        assert(silex::detail::order_unit_group_set_units_internal(
                computed, c.order, span, embeddings, 256, true));
        silex::OrderUnitGroup computed_out(c.order);
        bool changed = false;
        assert(silex::detail::adjoin_verified_dependent_relation(
                changed, computed_out, computed, root,
                sflint::FmpzMatConstRef(rel), 0, embeddings, 256));
        assert(changed);
        assert(torsion_order_is(computed_out, 14));

        silex::OrderUnitGroup group(c.order);
        assert(silex::detail::order_unit_group_set_units_internal(
                group, c.order, span, embeddings, 256, true, &two,
                &minus_one));
        silex::OrderUnitGroup out(c.order);
        changed = false;
        assert(silex::detail::adjoin_verified_dependent_relation(
                changed, out, group, root, sflint::FmpzMatConstRef(rel), 0,
                embeddings, 256));
        assert(changed);
        assert(torsion_order_is(out, 2));
        assert(same_regulator(out, computed_out));
        assert(silex::detail::order_unit_torsion_is_computed(computed_out));
        assert(!silex::detail::order_unit_torsion_is_computed(out));
    }
    return 0;
}

// Stored units past the reduction's factor-count cutoff (20000 factors) are
// rebuilt as single field elements equal to the units the unpadded basis
// reduces to.  The padding (u2 u2^-1)^10000 is a formal product equal to 1
// that the factored form does not merge, so u2 is stored with 20001
// factors.  The reduced basis {u2, u3} takes the identity-transform path,
// the skewed basis {u2, u2^5 u3} the transform path.  The rebuild needs an
// order known to be maximal; on Z[zeta_7] not marked maximal it fails, and
// the reduction keeps the factored units (identity transform) or the merged
// transformed products (skewed basis).
int test_reduce_stored_units_rebuilds_large_factored_units() {
    Cyclotomic7Units c;
    silex::EmbeddingContext embeddings(c.field);
    silex::Order maximal(c.field);
    assert(maximal.maximal_order(c.order));
    constexpr slong kPaddingPairs = 10000;
    for (const bool known_maximal : {true, false}) {
        const silex::Order& order = known_maximal ? maximal : c.order;
        for (const slong e21 : {slong(0), slong(5)}) {
            const std::vector<silex::FactoredElement> plain = c.basis(1, e21);
            std::vector<silex::FactoredElement> padded = c.basis(1, e21);
            for (slong k = 0; k < kPaddingPairs; ++k) {
                assert(padded[0].push(c.u2, 1));
                assert(padded[0].push(c.u2, -1));
            }
            assert(padded[0].length() == 2 * kPaddingPairs + 1);

            silex::OrderUnitGroup plain_group(order);
            silex::OrderUnitGroup padded_group(order);
            assert(silex::detail::order_unit_group_set_units_internal(
                    plain_group, order,
                    silex::FactoredElementSpan(plain.data(), plain.size()),
                    embeddings, 256, true));
            assert(silex::detail::order_unit_group_set_units_internal(
                    padded_group, order,
                    silex::FactoredElementSpan(padded.data(), padded.size()),
                    embeddings, 256, true));
            bool plain_changed = false;
            bool padded_changed = false;
            assert(silex::detail::reduce_stored_relation_units(
                    plain_changed, plain_group, embeddings, 256));
            assert(silex::detail::reduce_stored_relation_units(
                    padded_changed, padded_group, embeddings, 256));
            assert(plain_changed == (e21 != 0));
            assert(padded_changed == (known_maximal || e21 != 0));
            assert(padded_group.free_rank() == 2);
            assert(silex::detail::OrderUnitGroupAccess::wide_free_generators(
                           padded_group) == nullptr);
            slong padded_length = 0;
            for (slong i = 0; i < 2; ++i) {
                silex::FactoredElement expected(c.field);
                silex::FactoredElement stored(c.field);
                assert(plain_group.free_generator(expected, i));
                assert(padded_group.free_generator(stored, i));
                padded_length += stored.length();
                silex::Element expected_value(c.field);
                silex::Element stored_value(c.field);
                assert(expected.evaluate(expected_value));
                assert(stored.evaluate(stored_value));
                assert(stored_value.equal(expected_value));
            }
            if (known_maximal) {
                assert(padded_length == 2);
            } else if (e21 == 0) {
                assert(padded_length > 2 * kPaddingPairs);
            }
            const auto plain_regulator = plain_group.regulator();
            const auto padded_regulator = padded_group.regulator();
            assert(plain_regulator.has_value() &&
                   padded_regulator.has_value());
            assert(::arb_overlaps(plain_regulator->raw(),
                                  padded_regulator->raw()) != 0);
        }
    }
    return 0;
}


// Runs the paired class/unit transaction for the monic irreducible
// polynomial with coefficients `coefficients` (constant term first) on its
// maximal order, with the factor-base bound of the CLI.
struct TransactionRun {
    bool ok = false;
    silex::detail::ClassUnitTransactionReport audit;
};

TransactionRun run_transaction(const std::vector<slong>& coefficients,
                               silex::CertificationMode mode,
                               slong max_candidates, slong max_relations) {
    sflint::FmpzPoly polynomial;
    slong degree = 0;
    for (const slong coefficient : coefficients) {
        sflint::fmpz_poly_set_coeff_si(polynomial, degree++, coefficient);
    }
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpzPolyConstRef(polynomial));
    silex::Order equation = silex::test::equation_order(field);
    silex::Order order(field);
    assert(order.maximal_order(equation));
    silex::ClassGroupComputeOptions options;
    options.max_candidates = max_candidates;
    options.max_relations = max_relations;
    options.requested_certification = mode;
    sflint::Fmpz bound;
    assert(silex::factor_base_class_group_bound(sflint::FmpzRef(bound),
                                                order));
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 2);
    }
    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    TransactionRun run;
    run.ok = silex::detail::compute_class_unit_transaction(
            units, class_group, order, sflint::FmpzConstRef(bound), options,
            128, run.audit);
    return run;
}

// compute_with_class_group records the last validation outcome of a failed
// attempt only when validation ran (outcome other than not_run).  No public
// input leaves compute_ok true with the summary at not_run, so the seam
// clears the summary after the failed attempt.  Q(sqrt 146) in proven mode
// with max_relations 9 fails after validation ends in no_progress.
int test_compute_failure_records_validation_outcome_only_when_it_ran() {
    const std::vector<slong> sqrt146 = {-146, 0, 1};
    const TransactionRun natural = run_transaction(
            sqrt146, silex::CertificationMode::proven, 5000, 9);
    assert(!natural.ok);
    assert(std::strcmp(natural.audit.failure_reason,
                       "class_unit_computation_failed") == 0);
    assert(natural.audit.last_validation_outcome != nullptr);
    assert(std::strcmp(natural.audit.last_validation_outcome,
                       "no_progress") == 0);

    silex::detail::TestHookAccess::set_forget_validation_summary(true);
    const TransactionRun forgotten = run_transaction(
            sqrt146, silex::CertificationMode::proven, 5000, 9);
    silex::detail::TestHookAccess::set_forget_validation_summary(false);
    assert(!forgotten.ok);
    assert(std::strcmp(forgotten.audit.failure_reason,
                       "class_unit_computation_failed") == 0);
    assert(forgotten.audit.last_validation_outcome == nullptr);
    assert(!forgotten.audit.final_result_published);

    // The seam is off again: the run repeats the natural result.
    const TransactionRun again = run_transaction(
            sqrt146, silex::CertificationMode::proven, 5000, 9);
    assert(!again.ok && again.audit.last_validation_outcome != nullptr);
    return 0;
}

// The grh continuation fails closed when extend_tentative_relations_ reports
// success but relation_count() did not grow.  No public cap reaches this, so
// the seam arms one continuation that skips the extension.  x^3 - 4x + 1 in
// grh mode needs one continuation to succeed.  The exit leaves no cap
// recorded: the latest relation-search slice was not stopped by a cap, so
// the failure is the generic class_unit_computation_failed.
int test_grh_continuation_without_relation_progress_fails_closed() {
    const std::vector<slong> cubic = {1, -4, 0, 1};
    const TransactionRun natural = run_transaction(
            cubic, silex::CertificationMode::grh, 5000, 500);
    assert(natural.ok && natural.audit.final_result_published);

    silex::detail::TestHookAccess::set_continuation_appends_nothing(true);
    const TransactionRun stalled = run_transaction(
            cubic, silex::CertificationMode::grh, 5000, 500);
    silex::detail::TestHookAccess::set_continuation_appends_nothing(false);
    assert(!stalled.ok);
    assert(!stalled.audit.final_result_published);
    assert(std::strcmp(stalled.audit.failure_reason,
                       "class_unit_computation_failed") == 0);
    assert(stalled.audit.relation_search_budget_cap ==
           silex::detail::RelationSearchBudgetCap::none);

    // The arming was consumed: a repeat succeeds like the natural run.
    const TransactionRun repeat = run_transaction(
            cubic, silex::CertificationMode::grh, 5000, 500);
    assert(repeat.ok);
    return 0;
}

}  // namespace

int main() {
    test_compute_failure_records_validation_outcome_only_when_it_ran();
    test_grh_continuation_without_relation_progress_fails_closed();
    test_grh_validation_saturation_memo_replays_unchanged_pass();
    test_grh_validation_saturation_memo_invalidated_by_change();
    test_unit_extraction_cache_clear();
    test_precision_doubling_cap_default_value();
    test_precision_doubling_cap_forced_failures();
    test_reduce_mod_units_half_integer_tie();
    test_apply_unit_reduction_factor_list();
    test_reduce_mod_units_half_integer_tie_stop_rule();
    test_factored_log_memo_matches_direct_call();
    test_factored_log_memo_after_failed_factor_log();
    test_round_arb_approx_width_check();
    test_degree_one_rank_zero_group();
    test_compute_real_quadratic_proven();
    test_compute_with_class_group_rank_zero();
    test_compute_with_class_group_real_quadratic();
    test_compute_with_class_group_proven_saturated_generator();
    test_compute_with_class_group_real_quadratic_210();
    test_nonmaximal_quadratic_pair_rejection_preserves_output();
    test_compute_with_class_group_quintic_proven();
    test_compute_with_class_group_quintic_small_bf_budget_trivial();
    test_compute_with_class_group_failures_preserve_output();
    test_relation_kernel_rescale_log_matrix_word_rounding();
    test_reduced_regulator_reconstruction_correction();
    test_hnf_regulator_multiple_certification_correction();
    test_regulator_pivot_exponent_policy();
    test_regulator_pivot_arb_ambiguity_and_invalidity();
    test_regulator_pivot_clean_cols_and_mapping();
    test_regulator_pivot_multi_column_elimination();
    test_regulator_pivot_finish_route_matrix();
    test_relation_kernel_units_rank_one();
    test_relation_kernel_units_bounded_rank_one();
    test_relation_kernel_units_index_bounded_rank_one();
    test_relation_kernel_units_index_bounded_quartic_power_root();
    test_relation_kernel_units_index_bounded_cubic_rank_target();
    test_relation_kernel_units_index_bounded_saturated();
    test_relation_kernel_units_rank_zero_delegates();
    test_small_index_saturation_revalidates_appended_class_relation();
    test_validation_index_one_requires_generation_to_required_bound();
    test_relation_kernel_unit_failures_preserve_output();
    test_saturate_row_real_quadratic_square();
    test_saturate_row_no_root_and_divisible_copy();
    test_saturate_local_with_kernel_shares_field_modulus_cache();
    test_dependent_relation_rational_candidate();
    test_dependent_relation_denominator_step();
    test_dependent_relation_min_denominator_scan();
    test_rank_zero_torsion_generators();
    test_add_dependent_unit_large_denominator();
    test_add_dependent_unit_overflowing_index_bound_fallback();
    test_dependent_relation_wrong_simplest_rational_rejected();
    test_dependent_relation_y_root_uses_torsion_check();
    test_add_dependent_unit_precision_cap();
    test_saturate_row_rank_zero_and_failures();
    test_residue_dlog_kernel_real_quadratic();
    test_residue_dlog_proof_kernel_torsion();
    test_residue_dlog_proof_kernel_nonpower_basis_order();
    test_rank_one_degree_one_root_image_matches_prime_ideal();
    test_direct_degree_one_proof_column_matches_residue_field();
    test_direct_degree_one_residue_dlog_matrix_matches_residue_field();
    test_residue_dlog_proof_matrix_matches_per_entry_quotient_log();
    test_residue_dlog_proof_matrix_baby_giant_matches_per_entry_quotient_log();
    test_residue_quotient_log_column_mixes_resolved_and_apply_entries();
    test_residue_quotient_log_column_route_follows_the_bit_cutoff();
    test_half_integral_direct_degree_one_residue_dlog_matches_residue_field();
    test_select_saturation_proof_primes();
    test_unit_proof_records();
    test_prove_local_saturated_real_quadratic();
    test_regulator_index_bound();
    test_class_regulator_index_bound();
    test_class_regulator_index_bound_interval_boundary();
    test_class_regulator_index_bound_below_one_fails_closed();
    test_index_bounds_accept_upper_endpoint_exactly_one();
    test_unit_index_bound_from_regulator_quotient();
    test_class_unit_regulator_certification();
    test_belabas_friedman_class_regulator_is_grh_conditional();
    test_saturation_proves_nontrivial_class_group_without_grh();
    test_exact_imaginary_quadratic_saturation_records();
    test_cached_torsion_never_reaches_proven();
    test_rank_zero_quadratic_index_one_requires_computed_torsion();
    test_unit_proofs_require_computed_torsion();
    test_unit_proofs_require_computed_torsion_positive_rank();
    test_cached_torsion_is_valid();
    test_prove_index_bound();
    test_prove_index_bound_nonmaximal_quadratic_sqrt5();
    test_adjoin_dependent_relation_nonmaximal_root();
    test_prove_index_bound_maximal_quadratic_sqrt5_control();
    test_prove_index_bound_cubic2213_past_prescan_bound();
    test_saturate_index_bounded_generators_independent_of_bound();
    test_prove_index_bound_quartic70640();
    test_prove_index_bound_cyclotomic13_undecided_root();
    test_stable_proof_prime_columns_match_residue_field_route();
    test_stable_proof_prime_columns_fallback_and_wide_routes();
    test_kernel_row_root_three_way_result();
    test_kernel_row_root_defers_undecidable_rows();
    test_dependent_relation_unbounded_root_failure_adjoins_y();
    test_prove_index_bound_non_power_range_exhausted();
    test_prove_index_bound_prescan_stable_first();
    test_prove_index_bound_roots_power_after_deferred_rescan();
    test_prove_index_bound_range_exhausted_tests_rows_at_once();
    test_prove_index_bound_undecided_root_not_saturated();
    test_prove_torsion_twisted_saturation_root();
    test_saturate_torsion_twisted_non_proof();
    test_prove_index_bound_nonmaximal_quadratic_sqrt18();
    test_saturate_torsion_twisted_non_proof_nonmaximal_sqrt18();
    test_saturate_torsion_twisted_non_proof_odd_ell();
    test_saturate_index_bounded();
    test_saturate_index_bounded_adaptive();
    test_saturate_local_once_real_quadratic();
    test_saturate_local_once_rank_zero_and_failures();
    test_select_saturation_primes();
    test_saturation_selectors_agree();
    test_saturation_selectors_source_shape_fallback();
    test_saturate_bounded_real_quadratic();
    test_saturate_bounded_rank_zero_zero_pass_partial_and_failures();
    test_supplied_real_quadratic_units();
    test_supplied_rank_zero_units();
    test_supplied_unit_failures_preserve_output();
    test_imaginary_quadratic_torsion();
    test_nonmaximal_torsion_filtering();
    test_positive_rank_failure_preserves_output();
    test_nonmaximal_monic_quadratic_units();
    test_group_keeps_parent_order_alive();
    test_set_move_and_access_failures();
    test_transform_factored_units_wide_exponents();
    test_transform_factored_units_wide_merge_fits();
    test_transform_factored_units_overflow_flag();
    test_transform_factored_units_wide_failure();
    test_group_wide_form_failure_move_clear();
    test_reduce_stored_units_keeps_wide_form();
    test_rebasis_keeps_group_torsion();
    test_reduce_stored_units_rebuilds_large_factored_units();
    test_group_wide_free_generators();
    return 0;
}
