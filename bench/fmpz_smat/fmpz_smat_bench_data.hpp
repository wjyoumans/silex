#ifndef SILEX_BENCH_FMPZ_SMAT_DATA_HPP
#define SILEX_BENCH_FMPZ_SMAT_DATA_HPP

#include <silex/fmpz_smat.hpp>

namespace silex::bench::fmpz_smat {

inline void set_entry_si(fmpz_mat_t matrix,
        slong row,
        slong col,
        slong value) noexcept {
    fmpz_set_si(fmpz_mat_entry(matrix, row, col), value);
}

inline void set_entry_ui_2exp(fmpz_mat_t matrix,
        slong row,
        slong col,
        ulong value,
        ulong exponent) noexcept {
    fmpz_set_ui(fmpz_mat_entry(matrix, row, col), value);
    fmpz_mul_2exp(fmpz_mat_entry(matrix, row, col),
            fmpz_mat_entry(matrix, row, col),
            exponent);
}

inline void set_sparse_kernel_matrix(silex::fmpz_smat::SparseMat& matrix) noexcept {
    fmpz_mat_t dense;
    fmpz_mat_init(dense, 3, 4);

    set_entry_ui_2exp(dense, 0, 0, 1, 80);
    set_entry_si(dense, 0, 3, -7);
    set_entry_si(dense, 2, 1, 5);
    set_entry_si(dense, 2, 2, -2);

    matrix.set_fmpz_mat(dense);
    fmpz_mat_clear(dense);
}

inline void set_sparse_kernel_right(fmpz_mat_t right) noexcept {
    set_entry_si(right, 0, 0, 3);
    set_entry_si(right, 0, 1, -1);
    set_entry_si(right, 1, 0, 4);
    set_entry_si(right, 1, 1, 2);
    set_entry_si(right, 2, 0, -6);
    set_entry_si(right, 2, 1, 5);
    set_entry_si(right, 3, 0, 8);
    set_entry_si(right, 3, 1, -3);
}

inline void set_rank_matrix(fmpz_mat_t dense) noexcept {
    set_entry_si(dense, 0, 0, 2);
    set_entry_si(dense, 1, 1, 2);
}

inline void set_mod_rank_context_rows(fmpz_mat_t dense) noexcept {
    set_entry_si(dense, 0, 3, 1);
    set_entry_si(dense, 1, 1, -2);
    set_entry_si(dense, 1, 4, 5);
    set_entry_si(dense, 2, 0, 9);
    set_entry_si(dense, 3, 0, 6);
    set_entry_si(dense, 3, 2, 3);
    set_entry_si(dense, 4, 1, 1);
    set_entry_si(dense, 4, 3, 4);
}

inline void set_hnf_context_rows(fmpz_mat_t dense) noexcept {
    set_entry_si(dense, 0, 0, 2);
    set_entry_si(dense, 1, 1, 5);
    set_entry_si(dense, 2, 0, 4);
    set_entry_si(dense, 2, 1, 10);
}

}  // namespace silex::bench::fmpz_smat

#endif  // SILEX_BENCH_FMPZ_SMAT_DATA_HPP
