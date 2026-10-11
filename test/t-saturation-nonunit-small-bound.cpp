// Regression test for saturation_process_candidates_once on a non-unit
// candidate whose residual conjugate-log bound is small.  The residual
// (unit-only) presentation reduction assumes an empty ideal decomposition,
// so it must be taken only for unit candidates; before the fix a non-unit
// candidate with log bound at most 1 took it, failed the compact
// evaluation invariant, and the call returned false instead of "power" or
// "not a power".
//
// Fixture: Q(sqrt 2), factor base bound 2 (the prime P over 2, 2 = P^2),
// one relation with a given generator, prime n = 2, candidate = that
// relation.  2 = (sqrt 2)^2 is a square; -2 is not (it is negative at both
// real places).  Both have row (2) and log bound ceil(log 2) = 1.

#include <silex/class_group.hpp>
#include <silex/relation.hpp>

#include "class_group/relation_saturation_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstdio>

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
    return 0;
}
