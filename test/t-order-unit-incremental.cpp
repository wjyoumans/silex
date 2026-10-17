#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>

#include "order_unit/compute_internal.hpp"
#include "order_unit/order_unit_internal.hpp"
#include "order_unit/precision_cap_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct PrecisionCapOverride {
    explicit PrecisionCapOverride(slong cap) noexcept {
        silex::detail::set_precision_doubling_cap_for_testing(cap);
    }
    ~PrecisionCapOverride() {
        silex::detail::set_precision_doubling_cap_for_testing(0);
    }
};

// The real cyclotomic quartic x^4 - x^3 - 4x^2 + 4x + 1 (Q(zeta_15)^+,
// unit rank 3), in which theta is a unit (norm f(0) = 1).
silex::NumberField real_quartic_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 4, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, -1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 2, -4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 4);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, 1);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

bool has_state_at(const silex::detail::IncrementalUnitContext& context,
                  slong precision) noexcept {
    for (const auto& state : context.log_states) {
        if (state.precision == precision) {
            return true;
        }
    }
    return false;
}

// The dependent candidate theta^2 after theta.  Its Gram determinant is
// zero, so the incremental test decides it only once the ball lies below
// compact_independent's lower bound (about 2^-42 for degree 4, r = 3).  From
// 8 bits the ball radius is far above that bound, and the test goes straight
// to a precision near where the radius reaches it instead of doubling
// through 16 and 32 bits.  The decision is compact_independent's, the
// candidate is rejected, and nothing falls back.
int test_dependent_candidate_skips_to_deciding_precision() {
    silex::NumberField field = real_quartic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::FactoredElement first(field);
    silex::FactoredElement square(field);
    assert(first.push(theta, 1));
    assert(square.push(theta, 2));
    silex::EmbeddingContext embeddings(field);

    constexpr slong precision = 8;
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 3));
    bool changed = false;
    assert(context.add_independent_unit(changed, first, embeddings,
                                        precision, nullptr));
    assert(changed);
    assert(context.log_states.size() == 1);

    std::vector<silex::FactoredElement> trial;
    trial.emplace_back(field);
    assert(trial.back().set(first));
    trial.emplace_back(field);
    assert(trial.back().set(square));
    bool expected = true;
    assert(silex::detail::compact_independent(
            expected, embeddings,
            silex::FactoredElementSpan(trial.data(), trial.size()),
            precision));
    assert(!expected);

    changed = true;
    assert(context.add_independent_unit(changed, square, embeddings,
                                        precision, nullptr));
    assert(!changed);
    assert(context.units.size() == 1);
    assert(context.incremental_decisions == 2);
    assert(context.fallback_decisions == 0);
    assert(context.log_states.size() == 2);
    assert(has_state_at(context, precision));
    assert(!has_state_at(context, 2 * precision));
    assert(context.log_states.back().precision > 4 * precision);
    return 0;
}

// With the cap forced to 32 bits, from 8 bits the test skips 16, cannot
// decide theta^2 at 32 (the radius is still above the lower bound), and then
// tries the skipped 16 bits before giving up, so it tries every precision
// the plain doubling loop tries.  The candidate falls back to
// compact_independent, which also stops at the cap, and is not kept, as
// with plain doubling.
int test_skipped_precisions_tried_before_fallback() {
    silex::NumberField field = real_quartic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::FactoredElement first(field);
    silex::FactoredElement square(field);
    assert(first.push(theta, 1));
    assert(square.push(theta, 2));
    silex::EmbeddingContext embeddings(field);

    constexpr slong precision = 8;
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 3));
    PrecisionCapOverride cap(32);
    bool changed = false;
    assert(context.add_independent_unit(changed, first, embeddings,
                                        precision, nullptr));
    assert(changed);

    changed = true;
    assert(context.add_independent_unit(changed, square, embeddings,
                                        precision, nullptr));
    assert(!changed);
    assert(context.units.size() == 1);
    assert(context.incremental_decisions == 1);
    assert(context.fallback_decisions == 1);
    assert(context.log_states.size() == 3);
    assert(has_state_at(context, precision));
    assert(has_state_at(context, 4 * precision));
    assert(context.log_states.back().precision == 2 * precision);
    return 0;
}

}  // namespace

int main() {
    assert(test_dependent_candidate_skips_to_deciding_precision() == 0);
    assert(test_skipped_precisions_tried_before_fallback() == 0);
    return 0;
}
