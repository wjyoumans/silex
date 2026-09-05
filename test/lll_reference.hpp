#pragma once

#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz_mat.hpp>

namespace silex::test {

inline void exact_binary_parameter(flint::Fmpq& out, double value) {
    flint::Arf exact;
    arf_set_d(exact.raw(), value);
    arf_get_fmpq(out.raw(), exact.raw());
}

// Test-only rational Gram-Schmidt oracle; see the source-reference document.
// This does not call either of the production reducedness checkers.
inline bool rational_lll_reduced(flint::FmpzMatConstRef input,
                                 double delta = 0.99, double eta = 0.51) {
    const slong rows = fmpz_mat_nrows(input.raw());
    const slong cols = fmpz_mat_ncols(input.raw());
    flint::FmpqMat orthogonal(rows, cols), mu(rows, rows), norms(rows, 1);
    flint::Fmpq d, e, dot, term, coefficient;
    exact_binary_parameter(d, delta);
    exact_binary_parameter(e, eta);
    for (slong i = 0; i < rows; ++i) {
        for (slong k = 0; k < cols; ++k) {
            fmpq_set_fmpz(fmpq_mat_entry(orthogonal.raw(), i, k),
                    fmpz_mat_entry(input.raw(), i, k));
        }
        for (slong j = 0; j < i; ++j) {
            fmpq_zero(dot.raw());
            for (slong k = 0; k < cols; ++k) {
                fmpq_set_fmpz(coefficient.raw(), fmpz_mat_entry(input.raw(), i, k));
                fmpq_mul(term.raw(), coefficient.raw(),
                        fmpq_mat_entry(orthogonal.raw(), j, k));
                fmpq_add(dot.raw(), dot.raw(), term.raw());
            }
            fmpq_div(fmpq_mat_entry(mu.raw(), i, j), dot.raw(),
                    fmpq_mat_entry(norms.raw(), j, 0));
            fmpq_abs(term.raw(), fmpq_mat_entry(mu.raw(), i, j));
            if (fmpq_cmp(term.raw(), e.raw()) > 0) return false;
            for (slong k = 0; k < cols; ++k) {
                fmpq_mul(term.raw(), fmpq_mat_entry(mu.raw(), i, j),
                        fmpq_mat_entry(orthogonal.raw(), j, k));
                fmpq_sub(fmpq_mat_entry(orthogonal.raw(), i, k),
                        fmpq_mat_entry(orthogonal.raw(), i, k), term.raw());
            }
        }
        auto* norm = fmpq_mat_entry(norms.raw(), i, 0);
        for (slong k = 0; k < cols; ++k) {
            fmpq_mul(term.raw(), fmpq_mat_entry(orthogonal.raw(), i, k),
                    fmpq_mat_entry(orthogonal.raw(), i, k));
            fmpq_add(norm, norm, term.raw());
        }
        if (fmpq_sgn(norm) <= 0) return false;
        if (i != 0) {
            fmpq_mul(term.raw(), fmpq_mat_entry(mu.raw(), i, i - 1),
                    fmpq_mat_entry(mu.raw(), i, i - 1));
            fmpq_sub(term.raw(), d.raw(), term.raw());
            fmpq_mul(term.raw(), term.raw(), fmpq_mat_entry(norms.raw(), i - 1, 0));
            if (fmpq_cmp(term.raw(), norm) > 0) return false;
        }
    }
    return true;
}

}  // namespace silex::test
