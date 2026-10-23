#include <flint/arb_mat.h>
#include <flint/fmpq_mat.h>

#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arb_vec.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/order_unit.hpp>

#include "order_unit/compute_internal.hpp"
#include "order_unit/order_unit_internal.hpp"
#include "order_unit/precision_cap_internal.hpp"
#include "test_hook_access.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct PrecisionCapOverride {
    explicit PrecisionCapOverride(slong cap) noexcept {
        silex::detail::TestHookAccess::set_precision_doubling_cap(cap);
    }
    ~PrecisionCapOverride() {
        silex::detail::TestHookAccess::set_precision_doubling_cap(0);
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

struct IncrementalTargetSlackOverride {
    explicit IncrementalTargetSlackOverride(slong bits) noexcept {
        silex::detail::TestHookAccess::set_incremental_target_slack(bits);
    }
    ~IncrementalTargetSlackOverride() {
        silex::detail::TestHookAccess::set_incremental_target_slack(
                silex::detail::IncrementalTargetTestSeams::kOff);
    }
};

using RelationMode = silex::detail::IncrementalRelationTestSeams::Mode;

// Forces the relation step of the incremental test (see
// incremental_relation_dependent in compute.cpp) for the scope.  The tests
// of the lower-bound ladder below force it to fail_reconstruction, so their
// dependent candidates are decided by the lower bound B, as before the
// relation step.
struct IncrementalRelationOverride {
    explicit IncrementalRelationOverride(RelationMode mode) noexcept {
        silex::detail::TestHookAccess::set_incremental_relation_mode(mode);
    }
    ~IncrementalRelationOverride() {
        silex::detail::TestHookAccess::set_incremental_relation_mode(
                RelationMode::off);
    }
};

// The real cyclotomic octic Q(zeta_17)^+, the field of 2 cos(2 pi / 17),
// with minimal polynomial x^8 + x^7 - 7x^6 - 6x^5 + 15x^4 + 10x^3 - 10x^2
// - 4x + 1 (unit rank 7), in which theta is a unit (norm f(0) = 1).
silex::NumberField real_octic_field() noexcept {
    constexpr slong coefficients[9] = {1, -4, -10, 10, 15, -6, -7, 1, 1};
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (slong i = 0; i < 9; ++i) {
        sflint::fmpq_poly_set_coeff_si(polynomial, i, coefficients[i]);
    }
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

std::vector<slong> state_precisions(
        const silex::detail::IncrementalUnitContext& context) {
    std::vector<slong> precisions;
    for (const auto& state : context.log_states) {
        precisions.push_back(state.precision);
    }
    return precisions;
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

// Integer rows of an invertible 5 x 6 matrix, divided by 3 so that the
// row balls are inexact at a finite precision.
constexpr slong kBorderRows = 5;
constexpr slong kBorderColumns = 6;
constexpr slong kBorderNumerators[kBorderRows][kBorderColumns] = {
        {3, -1, 4, 1, -5, 9},
        {2, 6, -5, 3, 5, -8},
        {-9, 7, 9, 3, 2, 3},
        {8, -4, 6, 2, -6, 4},
        {3, 3, -8, 3, 2, 7}};

// Each step extends the previous step's bordered enclosures by one row, as
// add_independent_unit does after each accepted unit, starting from the
// empty Gram matrix (inverse 0 x 0, determinant 1).  After every step each
// entry of the bordered inverse must contain the exact rational entry of
// the inverse of the Gram matrix of the rows so far, and the determinant
// ball the exact determinant; they must also overlap arb_mat_inv and
// arb_mat_det of the ball Gram matrix, the full recompute, since both
// enclose the same values.
int test_bordered_inverse_encloses_exact_inverse() {
    constexpr slong precision = 40;
    std::vector<sflint::ArbVec> rows;
    sflint::ArbMat inverse(0, 0);
    sflint::Arb determinant;
    arb_one(determinant.raw());
    sflint::Fmpq entry;
    for (slong k = 0; k < kBorderRows; ++k) {
        sflint::ArbVec row(kBorderColumns);
        for (slong j = 0; j < kBorderColumns; ++j) {
            fmpq_set_si(entry.raw(), kBorderNumerators[k][j], 3);
            arb_set_fmpq(row.data() + j, entry.raw(), precision);
        }
        // projection = L v, inverse_projection = G^{-1} L v,
        // complement = v.v - projection . inverse_projection, computed as
        // incremental_independence_at does.
        sflint::ArbVec projection(k);
        sflint::ArbVec inverse_projection(k);
        sflint::Arb complement;
        for (slong i = 0; i < k; ++i) {
            arb_dot(projection.data() + i, nullptr, 0,
                    rows[static_cast<std::size_t>(i)].data(), 1, row.data(),
                    1, kBorderColumns, precision);
        }
        for (slong i = 0; i < k; ++i) {
            arb_dot(inverse_projection.data() + i, nullptr, 0,
                    arb_mat_entry(inverse.raw(), i, 0), 1, projection.data(),
                    1, k, precision);
        }
        sflint::Arb norm;
        arb_dot(norm.raw(), nullptr, 0, row.data(), 1, row.data(), 1,
                kBorderColumns, precision);
        arb_dot(complement.raw(), norm.raw(), 1, projection.data(), 1,
                inverse_projection.data(), 1, k, precision);
        rows.push_back(std::move(row));

        sflint::ArbMat next(0, 0);
        sflint::Arb next_determinant;
        assert(silex::detail::incremental_gram_bordered_inverse(
                next, next_determinant, inverse, determinant,
                inverse_projection, complement, precision));
        inverse.swap(next);
        determinant.swap(next_determinant);

        const slong size = k + 1;
        assert(arb_mat_nrows(inverse.raw()) == size);
        assert(arb_mat_ncols(inverse.raw()) == size);
        sflint::FmpqMat gram(size, size);
        sflint::FmpqMat exact_inverse(size, size);
        sflint::Fmpq exact_determinant;
        sflint::ArbMat ball_gram(size, size);
        for (slong i = 0; i < size; ++i) {
            for (slong j = 0; j < size; ++j) {
                slong dot = 0;
                for (slong c = 0; c < kBorderColumns; ++c) {
                    dot += kBorderNumerators[i][c] * kBorderNumerators[j][c];
                }
                fmpq_set_si(fmpq_mat_entry(gram.raw(), i, j), dot, 9);
                arb_dot(arb_mat_entry(ball_gram.raw(), i, j), nullptr, 0,
                        rows[static_cast<std::size_t>(i)].data(), 1,
                        rows[static_cast<std::size_t>(j)].data(), 1,
                        kBorderColumns, precision);
            }
        }
        assert(fmpq_mat_inv(exact_inverse.raw(), gram.raw()) != 0);
        fmpq_mat_det(exact_determinant.raw(), gram.raw());
        sflint::ArbMat full_inverse(size, size);
        sflint::Arb full_determinant;
        assert(arb_mat_inv(full_inverse.raw(), ball_gram.raw(), precision) !=
               0);
        arb_mat_det(full_determinant.raw(), ball_gram.raw(), precision);
        for (slong i = 0; i < size; ++i) {
            for (slong j = 0; j < size; ++j) {
                assert(arb_contains_fmpq(
                               arb_mat_entry(inverse.raw(), i, j),
                               fmpq_mat_entry(exact_inverse.raw(), i, j)) !=
                       0);
            }
        }
        assert(arb_mat_overlaps(inverse.raw(), full_inverse.raw()) != 0);
        assert(arb_contains_fmpq(determinant.raw(),
                                 exact_determinant.raw()) != 0);
        assert(arb_overlaps(determinant.raw(), full_determinant.raw()) != 0);
        assert(sflint::arb_is_positive(determinant));
    }
    return 0;
}

// The update declines when it cannot prove the bordered matrix invertible
// (a complement ball containing zero), on a determinant ball that is not
// positive, and on mismatched shapes, leaving the caller to recompute.
int test_bordered_inverse_declines_undecided_pivot() {
    constexpr slong precision = 64;
    sflint::ArbMat inverse(1, 1);
    arb_set_si(arb_mat_entry(inverse.raw(), 0, 0), 1);
    sflint::Arb determinant;
    arb_one(determinant.raw());
    sflint::ArbVec inverse_projection(1);
    arb_set_si(inverse_projection.data(), 1);
    sflint::Arb complement;
    sflint::ArbMat out(0, 0);
    sflint::Arb out_determinant;

    // complement = 0 +/- 1/2
    arb_zero(complement.raw());
    mag_set_d(arb_radref(complement.raw()), 0.5);
    assert(!silex::detail::incremental_gram_bordered_inverse(
            out, out_determinant, inverse, determinant, inverse_projection,
            complement, precision));

    arb_set_si(complement.raw(), 2);
    assert(silex::detail::incremental_gram_bordered_inverse(
            out, out_determinant, inverse, determinant, inverse_projection,
            complement, precision));

    sflint::Arb not_positive;
    arb_zero(not_positive.raw());
    mag_set_d(arb_radref(not_positive.raw()), 0.5);
    assert(!silex::detail::incremental_gram_bordered_inverse(
            out, out_determinant, inverse, not_positive, inverse_projection,
            complement, precision));

    sflint::ArbVec short_projection(0);
    assert(!silex::detail::incremental_gram_bordered_inverse(
            out, out_determinant, inverse, determinant, short_projection,
            complement, precision));
    return 0;
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
    IncrementalRelationOverride relation(RelationMode::fail_reconstruction);

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
    IncrementalRelationOverride relation(RelationMode::fail_reconstruction);

    constexpr slong precision = 8;
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 3));
    PrecisionCapOverride cap(32);
    bool changed = false;
    assert(context.add_independent_unit(changed, first, embeddings,
                                        precision, nullptr));
    assert(changed);
    // The accept extended the Gram enclosures at the decision precision.
    assert(context.log_states.size() == 1);
    assert(context.log_states.front().gram_bordered);
    assert(context.log_states.front().gram_units == 1);

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
    // Before the fallback the test ran again on full recomputes, which
    // replaced the bordered enclosures.
    for (const auto& state : context.log_states) {
        assert(!state.gram_bordered);
    }
    return 0;
}

// Accepting theta, theta - 1 and theta + 2 (all units: norms f(0), f(1),
// f(-2) are 1) with theta^2 in between.  After each accept short of full
// rank, the log state at the decision precision holds bordered Gram
// enclosures for every kept unit, and they overlap the full recompute
// (arb_mat_inv and arb_mat_det of the Gram matrix of its rows), which
// encloses the same values.  The decisions are compact_independent's.
int test_accept_borders_gram_enclosures() {
    silex::NumberField field = real_quartic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::Element minus_one(field);
    assert(minus_one.add_si(theta, -1));
    silex::Element plus_two(field);
    assert(plus_two.add_si(theta, 2));
    std::vector<silex::FactoredElement> candidates;
    std::vector<bool> independent;
    const auto add = [&](const silex::Element& base, slong exponent,
                         bool expected) {
        candidates.emplace_back(field);
        assert(candidates.back().push(base, exponent));
        independent.push_back(expected);
    };
    add(theta, 1, true);
    add(minus_one, 1, true);
    add(theta, 2, false);
    add(plus_two, 1, true);
    silex::EmbeddingContext embeddings(field);

    constexpr slong precision = 64;
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 3));
    for (std::size_t c = 0; c < candidates.size(); ++c) {
        bool changed = !independent[c];
        assert(context.add_independent_unit(changed, candidates[c],
                                            embeddings, precision, nullptr));
        assert(changed == independent[c]);
        if (!changed || context.full_rank) {
            continue;
        }
        const slong kept = static_cast<slong>(context.units.size());
        slong bordered = 0;
        for (const auto& state : context.log_states) {
            if (!state.gram_bordered) {
                continue;
            }
            ++bordered;
            assert(state.gram_units == kept);
            assert(static_cast<slong>(state.rows.size()) == kept);
            sflint::ArbMat gram(kept, kept);
            for (slong i = 0; i < kept; ++i) {
                for (slong j = 0; j < kept; ++j) {
                    const sflint::ArbVec& left =
                            state.rows[static_cast<std::size_t>(i)];
                    arb_dot(arb_mat_entry(gram.raw(), i, j), nullptr, 0,
                            left.data(), 1,
                            state.rows[static_cast<std::size_t>(j)].data(), 1,
                            left.length(), state.precision);
                }
            }
            sflint::ArbMat full_inverse(kept, kept);
            sflint::Arb full_determinant;
            assert(arb_mat_inv(full_inverse.raw(), gram.raw(),
                               state.precision) != 0);
            arb_mat_det(full_determinant.raw(), gram.raw(), state.precision);
            assert(arb_mat_overlaps(state.gram_inverse.raw(),
                                    full_inverse.raw()) != 0);
            assert(arb_overlaps(state.gram_determinant.raw(),
                                full_determinant.raw()) != 0);
        }
        assert(bordered == 1);
    }
    assert(context.full_rank);
    assert(context.units.size() == 3);
    assert(context.fallback_decisions == 0);
    return 0;
}

// The independent candidate theta - 1 (also a unit: norm f(1) = 1) after
// theta, started at 1 bit.  The Gram determinant ball is undecided at 1 bit
// and the test decides at 4 bits.  The accepted unit's log row must land in
// the log state at that decision precision, and only there.  (This does not
// separate the new schedule from plain doubling: the determinant of an
// independent pair is far above the lower bound, so both decide at 4 bits.)
int test_independent_candidate_accepted_after_skip() {
    silex::NumberField field = real_quartic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::Element shifted(field);
    assert(shifted.add_si(theta, -1));
    silex::FactoredElement first(field);
    silex::FactoredElement second(field);
    assert(first.push(theta, 1));
    assert(second.push(shifted, 1));
    silex::EmbeddingContext embeddings(field);

    constexpr slong precision = 1;
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 3));
    bool changed = false;
    assert(context.add_independent_unit(changed, first, embeddings,
                                        precision, nullptr));
    assert(changed);

    changed = false;
    assert(context.add_independent_unit(changed, second, embeddings,
                                        precision, nullptr));
    assert(changed);
    assert(context.units.size() == 2);
    assert(context.incremental_decisions == 2);
    assert(context.fallback_decisions == 0);
    for (const auto& state : context.log_states) {
        assert(state.rows.size() ==
               (state.precision == 4 * precision ? 2u : 1u));
    }
    assert(context.log_states.back().precision == 4 * precision);
    assert(context.log_states.back().rows.size() == 2);
    return 0;
}

// Adds theta (accepted at 7 bits) and then `candidate` to a context over
// the real octic at 7 bits, where the doubling precisions are 7, 14, ...,
// 112, 224, ...  Returns the context; `changed` reports whether the
// candidate was kept, which must agree with compact_independent.
silex::detail::IncrementalUnitContext octic_context_after(
        bool& changed,
        const silex::NumberField& field,
        silex::EmbeddingContext& embeddings,
        const silex::FactoredElement& first,
        const silex::FactoredElement& candidate,
        slong precision) noexcept {
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 7));
    changed = false;
    assert(context.add_independent_unit(changed, first, embeddings,
                                        precision, nullptr));
    assert(changed);
    assert(state_precisions(context) == std::vector<slong>{precision});

    std::vector<silex::FactoredElement> trial;
    trial.emplace_back(field);
    assert(trial.back().set(first));
    trial.emplace_back(field);
    assert(trial.back().set(candidate));
    bool expected = true;
    assert(silex::detail::compact_independent(
            expected, embeddings,
            silex::FactoredElementSpan(trial.data(), trial.size()),
            precision));

    changed = !expected;
    assert(context.add_independent_unit(changed, candidate, embeddings,
                                        precision, nullptr));
    assert(changed == expected);
    return context;
}

// The dependent candidate theta^2 after theta in the real octic, whose
// lower bound B is about 2^-125.  From 7 bits the ball is undecided and its
// target precision lies above 112, so plain doubling would decide at 224.
// The test instead tries the target plus 32 bits rounded up to a whole
// number of limbs, 192 bits, which lies below 224 and decides there.  The
// next dependent candidate, theta^3, starts at that precision and decides
// at once, so no doubling precision is ever tried above 7.
int test_dependent_candidate_decides_at_target_step() {
    silex::NumberField field = real_octic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::FactoredElement first(field);
    silex::FactoredElement square(field);
    silex::FactoredElement cube(field);
    assert(first.push(theta, 1));
    assert(square.push(theta, 2));
    assert(cube.push(theta, 3));
    silex::EmbeddingContext embeddings(field);
    IncrementalRelationOverride relation(RelationMode::fail_reconstruction);

    constexpr slong precision = 7;
    bool changed = true;
    silex::detail::IncrementalUnitContext context = octic_context_after(
            changed, field, embeddings, first, square, precision);
    assert(!changed);
    assert(context.units.size() == 1);
    assert(context.incremental_decisions == 2);
    assert(context.fallback_decisions == 0);
    assert(state_precisions(context) ==
           (std::vector<slong>{precision, 192}));
    assert(192 % FLINT_BITS == 0);
    assert(16 * precision < 192 && 192 < 32 * precision);
    assert(context.reject_precision == 192);

    changed = true;
    assert(context.add_independent_unit(changed, cube, embeddings,
                                        precision, nullptr));
    assert(!changed);
    assert(context.incremental_decisions == 3);
    assert(context.fallback_decisions == 0);
    assert(state_precisions(context) ==
           (std::vector<slong>{precision, 192}));
    assert(context.reject_precision == 192);
    return 0;
}

// The same candidate with the slack forced to -64 bits, so the target step
// lands short of the target and its ball stays undecided.  The test then
// takes the doubling precision it would have tried without the step, 224
// bits, and decides there.  The decision is compact_independent's, and
// nothing falls back to it.
int test_undecided_target_step_falls_back_to_doubling() {
    silex::NumberField field = real_octic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::FactoredElement first(field);
    silex::FactoredElement square(field);
    assert(first.push(theta, 1));
    assert(square.push(theta, 2));
    silex::EmbeddingContext embeddings(field);
    IncrementalRelationOverride relation(RelationMode::fail_reconstruction);

    constexpr slong precision = 7;
    IncrementalTargetSlackOverride slack(-64);
    bool changed = true;
    silex::detail::IncrementalUnitContext context = octic_context_after(
            changed, field, embeddings, first, square, precision);
    assert(!changed);
    assert(context.units.size() == 1);
    assert(context.incremental_decisions == 2);
    assert(context.fallback_decisions == 0);
    // 128 bits is ceil64(target - 64); 224 is the next doubling precision.
    assert(state_precisions(context) ==
           (std::vector<slong>{precision, 128, 32 * precision}));
    assert(context.reject_precision == 32 * precision);
    return 0;
}

// A candidate that starts at a reject precision between doublings and stays
// undecided up to the cap.  The dependent candidate theta^3 after theta, with
// the cap forced to 112 bits, starts at a reject precision of 96 bits, which
// lies between the doubling precisions 56 and 112 of the ladder 7, 14, ...,
// 112.  Its ball needs about 192 bits, so it is undecided at 96 and 112.  The
// test then tries the doubling precisions below the start last, so a log
// state exists at each of 14, 28 and 56 bits, and compact_independent, which
// also stops at the cap, leaves the candidate undecided.
int test_precisions_below_a_between_doubling_start_are_tried() {
    silex::NumberField field = real_octic_field();
    silex::Element theta(field);
    assert(theta.gen());
    silex::FactoredElement first(field);
    silex::FactoredElement cube(field);
    assert(first.push(theta, 1));
    assert(cube.push(theta, 3));
    silex::EmbeddingContext embeddings(field);
    IncrementalRelationOverride relation(RelationMode::fail_reconstruction);

    constexpr slong precision = 7;
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, precision));
    bool changed = false;
    assert(context.add_independent_unit(changed, first, embeddings,
                                        precision, nullptr));
    assert(changed);
    context.reject_precision = 96;

    PrecisionCapOverride cap(112);
    changed = true;
    assert(context.add_independent_unit(changed, cube, embeddings,
                                        precision, nullptr));
    assert(context.units.size() == 1);
    const std::vector<slong> states = state_precisions(context);
    for (const slong below : {14, 28, 56, 96, 112}) {
        assert(std::find(states.begin(), states.end(), below) !=
               states.end());
    }
    return 0;
}

// The sextic x^6 - 2, of signature (2, 2) and unit rank 3, in which
// theta - 1 is a unit (norm f(1) = -1).  It has complex places, whose
// product-mode log entries the torsion check halves.
silex::NumberField sextic_with_complex_places() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 6, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -2);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

// Adds u^2 and then u^3 to a context of rank `rank` at `precision`, with
// `diagnostics` (nullable).  u^3 is dependent (z = (u^3)^2 (u^2)^-3 = 1,
// so the reconstructed relation has m = 2), as compact_independent
// confirms.  Returns the context.
silex::detail::IncrementalUnitContext square_then_cube(
        const silex::NumberField& field,
        silex::EmbeddingContext& embeddings,
        const silex::Element& u,
        slong rank,
        slong precision,
        const silex::DiagnosticsContext* diagnostics) noexcept {
    silex::FactoredElement square(field);
    silex::FactoredElement cube(field);
    assert(square.push(u, 2));
    assert(cube.push(u, 3));

    std::vector<silex::FactoredElement> trial;
    trial.emplace_back(field);
    assert(trial.back().set(square));
    trial.emplace_back(field);
    assert(trial.back().set(cube));
    bool expected = true;
    assert(silex::detail::compact_independent(
            expected, embeddings,
            silex::FactoredElementSpan(trial.data(), trial.size()),
            precision));
    assert(!expected);

    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, rank));
    bool changed = false;
    assert(context.add_independent_unit(changed, square, embeddings,
                                        precision, diagnostics));
    assert(changed);
    changed = true;
    assert(context.add_independent_unit(changed, cube, embeddings, precision,
                                        diagnostics));
    assert(!changed);
    assert(context.units.size() == 1);
    assert(context.incremental_decisions == 2);
    assert(context.fallback_decisions == 0);
    return context;
}

// The relation step decides a dependent candidate at the starting
// precision, where the Gram determinant ball is still undecided: in the
// real octic (B about 2^-125) at 64 bits, and in the sextic x^6 - 2 (B
// about 2^-46, complex places) at 32 bits.  The relation u^3 against u^2
// has m = 2, and its unit passes the torsion check, so no other precision
// is tried.
int test_dependent_candidate_decided_by_relation() {
    struct Case {
        silex::NumberField field;
        slong unit_shift;
        slong rank;
        slong precision;
    };
    Case cases[] = {{real_octic_field(), 0, 7, 64},
                    {sextic_with_complex_places(), -1, 3, 32}};
    for (Case& c : cases) {
        silex::Element theta(c.field);
        silex::Element u(c.field);
        assert(theta.gen());
        assert(u.add_si(theta, c.unit_shift));
        silex::EmbeddingContext embeddings(c.field);
        const silex::detail::IncrementalUnitContext context =
                square_then_cube(c.field, embeddings, u, c.rank, c.precision,
                                 nullptr);
        assert(context.relation_decisions == 1);
        assert(state_precisions(context) ==
               std::vector<slong>{c.precision});
        assert(context.reject_precision == c.precision);
    }
    return 0;
}

// A forced failure of the relation step, either no reconstructed relation
// or a perturbed one (m + 1, so the relation unit is u^3, not a root of
// unity, and the torsion check cannot pass), decides nothing: the dependent
// candidate goes on to the lower-bound path and is decided there, at a
// higher precision, as before the relation step.  The decision is the same
// and nothing falls back to compact_independent.
int test_forced_relation_failure_takes_lower_bound_path() {
    silex::NumberField field = real_octic_field();
    silex::Element theta(field);
    assert(theta.gen());
    constexpr slong precision = 64;
    for (const RelationMode mode : {RelationMode::fail_reconstruction,
                                    RelationMode::perturb_relation}) {
        IncrementalRelationOverride relation(mode);
        silex::EmbeddingContext embeddings(field);
        const silex::detail::IncrementalUnitContext context =
                square_then_cube(field, embeddings, theta, 7, precision,
                                 nullptr);
        assert(context.relation_decisions == 0);
        assert(context.log_states.size() > 1);
        assert(context.reject_precision > precision);
    }
    return 0;
}

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
void count_debug_failure(void* user,
                         silex::DiagnosticsModule,
                         silex::DebugLevel,
                         const char*,
                         const char*,
                         int,
                         const char*,
                         const char*) noexcept {
    ++*static_cast<slong*>(user);
}
#endif

// With expensive debug checks on, a candidate decided dependent by the
// relation step is checked against compact_independent, which must decide
// it and call it dependent (add_independent_unit's debug check); no check
// may fail.  Without debug checks this only repeats the decision.
int test_relation_decision_debug_check() {
    silex::NumberField field = sextic_with_complex_places();
    silex::Element theta(field);
    silex::Element u(field);
    assert(theta.gen());
    assert(u.add_si(theta, -1));
    silex::DiagnosticsContext diagnostics;
    silex::diagnostics_context_init(diagnostics);
    slong failures = 0;
#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    silex::diagnostics_set_debug_checks(
            diagnostics, silex::DebugLevel::expensive,
            silex::diagnostics_all_modules, count_debug_failure, &failures);
#endif
    silex::EmbeddingContext embeddings(field);
    const silex::detail::IncrementalUnitContext context =
            square_then_cube(field, embeddings, u, 3, 32, &diagnostics);
    assert(context.relation_decisions == 1);
    assert(failures == 0);
    return 0;
}

// A dependent candidate whose relation unit is a root of unity other than
// 1: -u against the kept u gives m = 1 and z = (-u) u^-1 = -1.  The relation
// step must accept this through the torsion check on z != 1.
int test_relation_unit_nontrivial_root_of_unity() {
    silex::NumberField field = sextic_with_complex_places();
    silex::Element theta(field);
    silex::Element u(field);
    silex::Element minus_one(field);
    assert(theta.gen());
    assert(u.add_si(theta, -1));
    assert(minus_one.set_si(-1));
    silex::FactoredElement kept(field);
    silex::FactoredElement negated(field);
    assert(kept.push(u, 1));
    assert(negated.push(minus_one, 1));
    assert(negated.push(u, 1));

    constexpr slong precision = 32;
    silex::EmbeddingContext embeddings(field);
    silex::detail::IncrementalUnitContext context;
    assert(context.define(field, 3));
    bool changed = false;
    assert(context.add_independent_unit(changed, kept, embeddings, precision,
                                        nullptr));
    assert(changed);
    changed = true;
    assert(context.add_independent_unit(changed, negated, embeddings,
                                        precision, nullptr));
    assert(!changed);
    assert(context.units.size() == 1);
    assert(context.relation_decisions == 1);
    assert(context.fallback_decisions == 0);
    return 0;
}

}  // namespace

int main() {
    assert(test_bordered_inverse_encloses_exact_inverse() == 0);
    assert(test_bordered_inverse_declines_undecided_pivot() == 0);
    assert(test_accept_borders_gram_enclosures() == 0);
    assert(test_dependent_candidate_skips_to_deciding_precision() == 0);
    assert(test_skipped_precisions_tried_before_fallback() == 0);
    assert(test_independent_candidate_accepted_after_skip() == 0);
    assert(test_dependent_candidate_decides_at_target_step() == 0);
    assert(test_undecided_target_step_falls_back_to_doubling() == 0);
    assert(test_dependent_candidate_decided_by_relation() == 0);
    assert(test_forced_relation_failure_takes_lower_bound_path() == 0);
    assert(test_relation_decision_debug_check() == 0);
    assert(test_relation_unit_nontrivial_root_of_unity() == 0);
    assert(test_precisions_below_a_between_doubling_start_are_tried() == 0);
    return 0;
}
