// Reproducer for a suspected defect in saturation_process_candidates_once:
// a non-unit candidate whose conjugate-log bound is at most 1 makes
// compact_use_residual_power_test choose the unit-only presentation
// reduction (empty ideal decomposition), whose invariant check then fails
// for a non-unit, so the call returns false instead of "square" or "not a
// square".
//
// Fixture: Q(sqrt 2), factor base bound 2 (the prime P over 2, 2 = P^2),
// one relation with generator 2 (row (2)), prime n = 2, candidate = that
// relation.  The element 2 = (sqrt 2)^2 is a square, so a correct run
// returns true with saw_power set.  The log bound of 2 is
// ceil(log 2) = 1 and sqrt(1) = 1, so the residual branch is taken.
//
// The test is registered DISABLED while the defect stands; re-enable it by
// removing the DISABLED property in CMakeLists.txt (search for
// silex-saturation-nonunit-small-bound) once the library is fixed.

#include <silex/class_group.hpp>
#include <silex/relation.hpp>

#include "class_group/relation_saturation_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstdio>

namespace {
namespace sflint = silex::flint;

bool run(slong generator) {
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
            result, context, sflint::FmpzMatConstRef(candidates), 2, false, nullptr, nullptr,
            nullptr, 128);
    std::printf("generator=%ld process_candidates_once ok=%d saw_power=%d "
                "saw_non_power=%d tested=%ld compact=%ld structural=%ld\n",
                static_cast<long>(generator), ok ? 1 : 0, result.saw_power ? 1 : 0,
                result.saw_non_power ? 1 : 0,
                static_cast<long>(result.tested_power_candidates),
                static_cast<long>(result.compact_power_candidates),
                static_cast<long>(result.structural_power_candidates));
    std::fflush(stdout);
    // Correct behavior: the candidate is a square and is reported as one.
    return ok && result.saw_power;
}
}  // namespace

int main() {
    // Control: 4 has log bound ceil(log 4) = 2 > sqrt(2), so the original
    // (non-residual) branch is taken.
    const bool control = run(4);
    assert(control);
    // Defect: 2 has log bound 1, so the residual branch is taken.
    const bool defect = run(2);
    assert(defect);
    return 0;
}
