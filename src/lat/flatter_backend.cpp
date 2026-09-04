#if defined(SILEX_WITH_FLATTER) && SILEX_WITH_FLATTER
#include <gmp.h>
#endif

#include "flatter_backend_internal.hpp"

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/flint.h>

#include <limits>

#if defined(SILEX_WITH_FLATTER) && SILEX_WITH_FLATTER
#include <flatter/computation_context.h>
#include <flatter/data/lattice.h>
#include <flatter/flatter.h>
#include <flatter/problems.h>

#endif

namespace silex::lat::detail {

namespace {

slong matrix_rank(flint::FmpzMatConstRef input) noexcept {
    return fmpz_mat_rank(input.raw());
}

bool valid_output_dimensions(flint::FmpzMatRef reduced,
                             flint::FmpzMatRef transform,
                             flint::FmpzMatConstRef input,
                             slong rank,
                             double root_hermite_factor) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    return rows >= 0 && cols >= 0 && rank >= 0 &&
           fmpz_mat_nrows(reduced.raw()) == rows &&
           fmpz_mat_ncols(reduced.raw()) == rank &&
           fmpz_mat_nrows(transform.raw()) == cols &&
           fmpz_mat_ncols(transform.raw()) == rank &&
           root_hermite_factor > 1.0;
}

bool dimensions_fit_uint(flint::FmpzMatConstRef input) noexcept {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    return rows <= std::numeric_limits<unsigned int>::max() &&
           cols <= std::numeric_limits<unsigned int>::max();
}

#if defined(SILEX_WITH_FLATTER) && SILEX_WITH_FLATTER
void set_mpz_from_fmpz(mpz_t out, const fmpz* value) noexcept {
    fmpz_get_mpz(out, value);
}

void set_fmpz_from_mpz(fmpz* out, const mpz_t value) noexcept {
    fmpz_set_mpz(out, value);
}
#endif

}  // namespace

FlatterLllResult flatter_column_lll_transform(
        flint::FmpzMatRef reduced,
        flint::FmpzMatRef transform,
        flint::FmpzMatConstRef input,
        double root_hermite_factor,
        unsigned int max_threads) noexcept {
    const slong rank = matrix_rank(input);
    if (!valid_output_dimensions(
                reduced, transform, input, rank, root_hermite_factor)) {
        return {FlatterBackendStatus::invalid_dimensions, rank};
    }
    if (!dimensions_fit_uint(input)) {
        return {FlatterBackendStatus::dimension_overflow, rank};
    }

#if defined(SILEX_WITH_FLATTER) && SILEX_WITH_FLATTER
    const auto rows = static_cast<unsigned int>(fmpz_mat_nrows(input.raw()));
    const auto cols = static_cast<unsigned int>(fmpz_mat_ncols(input.raw()));
    const auto output_rank = static_cast<unsigned int>(rank);

    flatter::Lattice lattice(cols, rows);
    auto lattice_basis = lattice.basis().data<mpz_t>();
    for (unsigned int row = 0; row < rows; ++row) {
        for (unsigned int col = 0; col < cols; ++col) {
            set_mpz_from_fmpz(
                    lattice_basis(row, col),
                    fmpz_mat_entry(input.raw(), row, col));
        }
    }

    flatter::Matrix backend_transform(
            flatter::ElementType::MPZ, cols, lattice.rank());
    flatter::initialize();
    flatter::ComputationContext context(max_threads == 0 ? 1 : max_threads);
    flatter::LatticeReductionParams params(
            lattice, backend_transform, root_hermite_factor);
    if (rank == static_cast<slong>(cols) && cols <= 32 &&
        lattice.basis().prec() <= 128) {
        params.phase = 3;
    }
    flatter::LatticeReduction reduction(params, context);
    reduction.solve();
    lattice.update_rank();
    flatter::finalize();

    if (lattice.rank() != output_rank) {
        return {FlatterBackendStatus::reduction_failed,
                static_cast<slong>(lattice.rank())};
    }

    const auto reduced_basis = lattice.basis().data<mpz_t>();
    const auto backend_transform_data = backend_transform.data<mpz_t>();
    for (unsigned int row = 0; row < rows; ++row) {
        for (unsigned int col = 0; col < output_rank; ++col) {
            set_fmpz_from_mpz(
                    fmpz_mat_entry(reduced.raw(), row, col),
                    reduced_basis(row, col));
        }
    }
    for (unsigned int row = 0; row < cols; ++row) {
        for (unsigned int col = 0; col < output_rank; ++col) {
            set_fmpz_from_mpz(
                    fmpz_mat_entry(transform.raw(), row, col),
                    backend_transform_data(row, col));
        }
    }
    flint::FmpzMat expected(rows, output_rank);
    flint::fmpz_mat_mul(flint::FmpzMatRef(expected), input,
                        flint::FmpzMatConstRef(transform.raw()));
    if (!flint::fmpz_mat_equal(flint::FmpzMatConstRef(expected),
                               flint::FmpzMatConstRef(reduced.raw()))) {
        return {FlatterBackendStatus::transform_unavailable, rank};
    }
    return {FlatterBackendStatus::success, rank};
#else
    (void)reduced;
    (void)transform;
    (void)input;
    (void)root_hermite_factor;
    (void)max_threads;
    return {FlatterBackendStatus::unavailable, rank};
#endif
}

}  // namespace silex::lat::detail
