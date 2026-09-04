#include "order_unit_internal.hpp"

namespace silex::detail {

bool unit_regulator_matches_reconstruction(
        flint::ArbConstRef actual_regulator,
        flint::ArbConstRef expected_regulator,
        slong precision) noexcept {
    if (precision <= 0 ||
        !flint::arb_is_finite(actual_regulator) ||
        !flint::arb_is_positive(actual_regulator) ||
        !flint::arb_is_finite(expected_regulator) ||
        !flint::arb_is_positive(expected_regulator)) {
        return false;
    }

    flint::Arb difference;
    flint::Arb abs_difference;
    flint::Arb one;
    flint::arb_sub(difference, actual_regulator.raw(),
                   expected_regulator.raw(), precision);
    flint::arb_abs(abs_difference, difference);
    flint::arb_set_ui(one, 1);
    return flint::arb_lt(abs_difference, one);
}

bool unit_regulator_matches_reconstruction(
        const OrderUnitGroup& units,
        flint::ArbConstRef expected_regulator,
        slong precision) noexcept {
    flint::Arb actual_regulator;
    return units.regulator(flint::ArbRef(actual_regulator)) &&
           unit_regulator_matches_reconstruction(
                   flint::ArbConstRef(actual_regulator), expected_regulator,
                   precision);
}

}  // namespace silex::detail
