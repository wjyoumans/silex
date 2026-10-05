#pragma once

#include <silex/flint/fmpz.hpp>

#include <flint/flint.h>

namespace silex::detail {

// A required prime that the factor-base honesty scan could not witness: the
// prime, the stage of its last witness search against the last stage
// allowed, and that stage's caps.  The lattice route uses radius,
// max_twists and random_tries; the T2 route uses random_tries,
// factor_attempts and element_steps.  The relation-search honesty
// checkpoint moves it into a class/unit transaction report when the
// failure ends the candidate attempt.  It is move-only, as its prime is.
struct FactorBaseHonestyFailure {
    bool recorded = false;
    flint::Fmpz p;
    slong residue_degree = 0;
    bool direct_witness_search = false;
    slong stage = 0;
    slong max_stage = 0;
    slong radius = 0;
    slong max_twists = 0;
    slong random_tries = 0;
    slong factor_attempts = 0;
    slong element_steps = 0;
};

}  // namespace silex::detail
