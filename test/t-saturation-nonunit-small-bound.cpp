// Regression test for saturation_process_candidates_once on a non-unit
// candidate whose residual conjugate-log bound is small.  The residual
// (unit-only) presentation reduction assumes an empty ideal decomposition,
// so it must be taken only when the formal leftover is a unit; before the
// fix a non-unit candidate with log bound at most 1 took it, failed the
// compact evaluation invariant, and the call returned false instead of
// "power" or "not a power".  run_unit_candidate_nonunit_leftover covers a
// unit candidate whose leftover is not a unit.
//
// Fixture: Q(sqrt 2), factor base bound 2 (the prime P over 2, 2 = P^2),
// one relation with a given generator, prime n = 2, candidate = that
// relation.  2 = (sqrt 2)^2 is a square; -2 is not (it is negative at both
// real places).  Both have row (2) and log bound ceil(log 2) = 1.

#include <silex/class_group.hpp>
#include <silex/relation.hpp>

#include "class_group/class_group_storage_internal.hpp"
#include "class_group/relation_saturation_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {
namespace sflint = silex::flint;

bool run(slong generator, bool expect_power) {
    silex::NumberField field = silex::test::quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 2);
    silex::ClassGroupContext context(order);
    assert(context.is_defined());
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);

    silex::Relation relation(*base);
    silex::Element alpha(field);
    assert(alpha.set_si(generator));
    assert(relation.set_generator(alpha));
    assert(context.append_relation(relation));
    assert(context.publish_presentation());

    sflint::FmpzMat candidates(1, 1);
    sflint::fmpz_set_si(
            sflint::FmpzRef(sflint::fmpz_mat_entry(candidates, 0, 0)), 1);

    silex::detail::SaturationCandidateProcessingResult result;
    const bool ok = silex::detail::saturation_process_candidates_once(
            result, context, sflint::FmpzMatConstRef(candidates), 2, true, nullptr, nullptr,
            nullptr, 128);
    std::printf("generator=%ld process_candidates_once ok=%d saw_power=%d "
                "saw_non_power=%d tested=%ld compact=%ld structural=%ld\n",
                static_cast<long>(generator), ok ? 1 : 0, result.saw_power ? 1 : 0,
                result.saw_non_power ? 1 : 0,
                static_cast<long>(result.tested_power_candidates),
                static_cast<long>(result.compact_power_candidates),
                static_cast<long>(result.structural_power_candidates));
    std::fflush(stdout);
    if (!ok || result.saw_power != expect_power ||
        result.saw_non_power != !expect_power) {
        return false;
    }
    // The compact branch ran, not the structural shortcut.
    if (result.compact_power_candidates != 1 ||
        result.structural_power_candidates != 0) {
        return false;
    }
    if (!expect_power) {
        return result.appended_class_relations == 0;
    }
    // The appended relation's generator is an n-th root of the candidate.
    if (result.appended_class_relations != 1) {
        return false;
    }
    silex::Element root(field);
    silex::Element square(field);
    if (context.relation_count() < 2 ||
        !context.relation_generator(root, context.relation_count() - 1) ||
        !square.multiply(root, root)) {
        return false;
    }
    return square.equal(alpha);
}

// A unit candidate whose formal leftover is not a unit.  Relations
// gamma = eps^20 sqrt 2 (row (1)) and delta = 2 / eps (row (2)), with
// eps = 1 + sqrt 2, and candidate coefficients 2 and -1 outside [0, n):
// gamma^2 delta^-1 = eps^41 has a zero row and is not a square (eps is
// fundamental of norm -1).  The formal split leaves delta, which is not a
// unit, with log bound 2 <= isqrt(37).  Before the residual branch was
// gated on a unit leftover, the unit-only reduction of delta failed the
// compact evaluation invariant and the call returned false.
struct ResidualNotUnitCounter {
    long events = 0;
};

void residual_not_unit_callback(void* user,
                                silex::DiagnosticsModule,
                                silex::ProfileEvent event,
                                const char*,
                                const char* label) noexcept {
    if (event == silex::ProfileEvent::event && label != nullptr &&
        std::strcmp(label, "class_group.saturation_residual_not_unit") == 0) {
        ++static_cast<ResidualNotUnitCounter*>(user)->events;
    }
}

bool run_unit_candidate_nonunit_leftover() {
    silex::NumberField field = silex::test::quadratic_field(2);
    silex::Order order;
    order = silex::test::equation_order(field);

    sflint::Fmpz bound;
    sflint::fmpz_set_si(sflint::FmpzRef(bound), 2);
    silex::ClassGroupContext context(order);
    assert(context.is_defined());
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    ResidualNotUnitCounter counter;
    silex::diagnostics_set_profiling(
            diagnostics, true,
            silex::diagnostics_module_bit(silex::DiagnosticsModule::class_group),
            residual_not_unit_callback, &counter);
    context.set_diagnostics(&diagnostics);
    assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
    const silex::FactorBase* base = context.factor_base();
    assert(base != nullptr);

    silex::Element sqrt2(field);
    silex::Element eps(field);
    silex::Element eps_power(field);
    silex::Element gamma(field);
    silex::Element delta(field);
    sflint::Fmpz twenty;
    sflint::fmpz_set_si(sflint::FmpzRef(twenty), 20);
    assert(sqrt2.gen() && eps.add_si(sqrt2, 1) &&
           eps_power.pow_fmpz(eps, sflint::FmpzConstRef(twenty)) &&
           gamma.multiply(eps_power, sqrt2));
    // 2 / eps = 2 (sqrt 2 - 1) = 2 sqrt 2 - 2.
    assert(delta.add(sqrt2, sqrt2) && delta.add_si(delta, -2));

    // delta's row is twice gamma's, so the dependent relation is appended
    // with the relation search's keep policy.
    for (const silex::Element* generator : {&gamma, &delta}) {
        silex::Relation relation(*base);
        assert(relation.set_generator(*generator));
        bool skipped = false;
        assert(silex::detail::ClassGroupRelationSearchAccess::
                       append_search_relation(skipped, context, relation));
        assert(!skipped);
    }
    assert(context.publish_presentation());

    // Coefficients follow the relation order of saturation_relation_data.
    std::vector<silex::FactoredElement> relations;
    sflint::FmpzMat rows(context.relation_count(), context.generator_count());
    assert(silex::detail::saturation_relation_data(
            relations, sflint::FmpzMatRef(rows), context));
    assert(context.relation_count() == 2 && context.generator_count() == 1);
    sflint::FmpzMat candidates(2, 1);
    for (slong j = 0; j < 2; ++j) {
        const slong row = sflint::fmpz_get_si(sflint::FmpzConstRef(
                sflint::fmpz_mat_entry(sflint::FmpzMatConstRef(rows), j, 0)));
        assert(row == 1 || row == 2);
        sflint::fmpz_set_si(
                sflint::FmpzRef(sflint::fmpz_mat_entry(candidates, j, 0)),
                row == 1 ? 2 : -1);
    }

    silex::detail::SaturationCandidateProcessingResult result;
    const bool ok = silex::detail::saturation_process_candidates_once(
            result, context, sflint::FmpzMatConstRef(candidates), 2, true,
            nullptr, nullptr, nullptr, 128);
    std::printf("unit candidate, non-unit leftover: ok=%d saw_power=%d "
                "saw_non_power=%d compact=%ld\n",
                ok ? 1 : 0, result.saw_power ? 1 : 0,
                result.saw_non_power ? 1 : 0,
                static_cast<long>(result.compact_power_candidates));
    std::fflush(stdout);
#if defined(SILEX_ENABLE_PROFILING) && SILEX_ENABLE_PROFILING
    // The residual branch was considered and rejected for the non-unit
    // leftover.
    assert(counter.events == 1);
#endif
    return ok && !result.saw_power && result.saw_non_power &&
           result.compact_power_candidates == 1 &&
           result.structural_power_candidates == 0 &&
           result.appended_class_relations == 0;
}
}  // namespace

int main() {
    // Control: 4 has log bound ceil(log 4) = 2 > sqrt(2) and was handled
    // before the fix.
    const bool control = run(4, true);
    assert(control);
    // 2 has log bound 1: a non-unit square.
    const bool square = run(2, true);
    assert(square);
    // -2 has log bound 1: a non-unit non-square.
    const bool non_square = run(-2, false);
    assert(non_square);
    const bool unit_candidate = run_unit_candidate_nonunit_leftover();
    assert(unit_candidate);
    return 0;
}
