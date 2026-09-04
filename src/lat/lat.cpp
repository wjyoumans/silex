#include <silex/lat.hpp>

#include <silex/diagnostics.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_lll.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/flint/fmpz_mod_ctx.hpp>
#include <silex/flint/fmpz_mod_mat.hpp>

#include <flint/arb_mat.h>
#include <flint/arf.h>
#include <flint/flint.h>

#include <cmath>
#include <vector>

namespace silex::lat {
namespace {

constexpr slong lat_enum_double_max_dim = 32;
constexpr slong lat_enum_double_max_coord = WORD(10000);
constexpr double lat_enum_double_min_diag = 1e-150;
constexpr double lat_enum_double_max_entry = 1e150;

struct EnumStateArb {
    arb_mat_struct* cholesky = nullptr;
    const arb_struct* bound_sq = nullptr;
    fmpz_mat_struct* coeffs = nullptr;
    slong max_coord = -1;
    slong prec = 0;
    bool precision_failure = false;
    bool aborted_by_callback = false;
    ShortVectorCallback callback = nullptr;
    void* user = nullptr;
};

struct EnumStateDouble {
    double* cholesky = nullptr;
    double* partial = nullptr;
    fmpz_mat_struct* coeffs = nullptr;
    slong dimension = 0;
    double bound = 0.0;
    slong max_coord = -1;
    bool aborted_by_callback = false;
    ShortVectorCallback callback = nullptr;
    void* user = nullptr;
};

slong nonzero_rows(const fmpz_mat_t matrix) noexcept {
    slong rows = fmpz_mat_nrows(matrix);
    while (rows > 0 && fmpz_mat_is_zero_row(matrix, rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (fmpz_mat_is_zero_row(matrix, i) != 0) {
            return -1;
        }
    }
    return rows;
}

bool hnf_basis(fmpz_mat_t out, const fmpz_mat_t basis) noexcept {
    flint::FmpzMat full(fmpz_mat_nrows(basis), fmpz_mat_ncols(basis));
    fmpz_mat_hnf(full.raw(), basis);

    const slong rows = nonzero_rows(full.raw());
    if (rows < 0) {
        return false;
    }

    flint::FmpzMatConstWindow window(full, 0, 0, rows, fmpz_mat_ncols(basis));
    fmpz_mat_clear(out);
    fmpz_mat_init_set(out, window.raw());
    return true;
}

bool integral_divexact_matrix(fmpz_mat_t out,
        const fmpz_mat_t matrix,
        const fmpz_t denominator) noexcept {
    if (fmpz_is_zero(denominator) != 0) {
        return false;
    }

    for (slong i = 0; i < fmpz_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < fmpz_mat_ncols(matrix); ++j) {
            if (fmpz_divisible(fmpz_mat_entry(matrix, i, j), denominator) == 0) {
                return false;
            }
        }
    }

    for (slong i = 0; i < fmpz_mat_nrows(matrix); ++i) {
        for (slong j = 0; j < fmpz_mat_ncols(matrix); ++j) {
            fmpz_divexact(fmpz_mat_entry(out, i, j),
                    fmpz_mat_entry(matrix, i, j),
                    denominator);
        }
    }
    return true;
}

bool same_ambient(const Lat& left, const Lat& right) noexcept {
    return left.ambient_dim() == right.ambient_dim();
}

#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
bool lat_check_failed(const DiagnosticsContext* diagnostics,
        DebugLevel level,
        const char* function,
        const char* file,
        int line,
        const char* label,
        const char* expression) noexcept {
    debug_check_failed(diagnostics, DiagnosticsModule::lattice, level, function,
            file, line, label, expression);
    return false;
}
#endif

bool enum_indeterminate_arb(EnumStateArb& state, const arb_t value) noexcept {
    if (flint::arb_is_negative(value)) {
        return false;
    }
    if (flint::arb_is_nonnegative(value)) {
        return true;
    }

    state.precision_failure = true;
    return true;
}

bool enum_recurse_arb(EnumStateArb& state,
        slong coordinate,
        const arb_t accumulated) noexcept {
    if (state.precision_failure || state.aborted_by_callback) {
        return false;
    }

    const slong dimension = arb_mat_nrows(state.cholesky);

    flint::Arb remaining;
    flint::arb_sub(remaining, state.bound_sq, accumulated, state.prec);
    if (!enum_indeterminate_arb(state, remaining.raw())) {
        return true;
    }

    flint::Arb delta;
    flint::Arb diag;
    flint::Arb abs_diag;
    flint::Arb tmp;

    flint::arb_set(diag, arb_mat_entry(state.cholesky, coordinate, coordinate));
    flint::arb_abs(abs_diag, diag);

    if (!flint::arb_is_positive(abs_diag)) {
        state.precision_failure = true;
        return false;
    }

    flint::arb_zero(delta);
    for (slong i = coordinate + 1; i < dimension; ++i) {
        flint::arb_mul_fmpz(tmp, arb_mat_entry(state.cholesky, i, coordinate),
                fmpz_mat_entry(state.coeffs, 0, i), state.prec);
        flint::arb_add(delta, delta, tmp, state.prec);
    }

    flint::Arb center;
    flint::Arb radius;
    flint::Arb half_width;
    flint::Arb lower_arb;
    flint::Arb upper_arb;
    flint::Arb coord_arb;
    flint::Arb y_coord;
    flint::Arb y_coord_sq;
    flint::Arb new_accumulated;

    flint::arb_neg(center, delta);
    flint::arb_div(center, center, diag, state.prec);
    flint::arb_sqrt(radius, remaining, state.prec);
    flint::arb_div(half_width, radius, abs_diag, state.prec);
    flint::arb_sub(lower_arb, center, half_width, state.prec);
    flint::arb_add(upper_arb, center, half_width, state.prec);

    flint::Fmpz lower;
    flint::Fmpz upper;
    flint::Fmpz current;

    flint::Arf lower_bound;
    flint::Arf upper_bound;
    flint::arb_get_lbound_arf(lower_bound, lower_arb, state.prec);
    flint::arb_get_ubound_arf(upper_bound, upper_arb, state.prec);

    if (flint::arf_is_nan(lower_bound) || flint::arf_is_inf(lower_bound) ||
        flint::arf_is_nan(upper_bound) || flint::arf_is_inf(upper_bound)) {
        state.precision_failure = true;
        return false;
    }

    flint::arf_get_fmpz(lower, lower_bound, ARF_RND_CEIL);
    flint::arf_get_fmpz(upper, upper_bound, ARF_RND_FLOOR);

    if (state.max_coord >= 0) {
        flint::Fmpz cap;
        fmpz_set_si(cap.raw(), -state.max_coord);
        if (fmpz_cmp(lower.raw(), cap.raw()) < 0) {
            fmpz_set(lower.raw(), cap.raw());
        }
        fmpz_set_si(cap.raw(), state.max_coord);
        if (fmpz_cmp(upper.raw(), cap.raw()) > 0) {
            fmpz_set(upper.raw(), cap.raw());
        }
    }

    if (fmpz_cmp(lower.raw(), upper.raw()) > 0) {
        return true;
    }

    fmpz_set(current.raw(), lower.raw());
    while (fmpz_cmp(current.raw(), upper.raw()) <= 0) {
        if (state.precision_failure || state.aborted_by_callback) {
            return !state.precision_failure;
        }

        flint::arb_set_fmpz(coord_arb, current);
        flint::arb_mul(y_coord, diag, coord_arb, state.prec);
        flint::arb_add(y_coord, y_coord, delta, state.prec);
        flint::arb_sqr(y_coord_sq, y_coord, state.prec);
        flint::arb_add(new_accumulated, accumulated, y_coord_sq, state.prec);

        flint::arb_sub(tmp, state.bound_sq, new_accumulated, state.prec);
        if (!flint::arb_is_negative(tmp)) {
            fmpz_set(fmpz_mat_entry(state.coeffs, 0, coordinate), current.raw());

            if (coordinate == 0) {
                bool nonzero = false;
                for (slong i = 0; i < dimension; ++i) {
                    if (fmpz_is_zero(fmpz_mat_entry(state.coeffs, 0, i)) == 0) {
                        nonzero = true;
                        break;
                    }
                }
                if (nonzero && state.callback(state.coeffs, state.user) == 0) {
                    state.aborted_by_callback = true;
                    return true;
                }
            } else {
                if (!enum_recurse_arb(state, coordinate - 1, new_accumulated.raw())) {
                    return false;
                }
                if (state.aborted_by_callback) {
                    return true;
                }
            }
        }

        fmpz_add_ui(current.raw(), current.raw(), UWORD(1));
    }

    return true;
}

bool enum_recurse_double(EnumStateDouble& state,
        slong coordinate,
        double accumulated) noexcept {
    if (state.aborted_by_callback) {
        return true;
    }

    const double remaining = state.bound - accumulated;
    if (remaining < 0.0) {
        return true;
    }

    const double diag = state.cholesky[coordinate * state.dimension + coordinate];
    if (!(diag > 0.0)) {
        return false;
    }

    double delta = 0.0;
    for (slong i = coordinate + 1; i < state.dimension; ++i) {
        delta += state.cholesky[i * state.dimension + coordinate] *
                 state.partial[i];
    }

    const double center = -delta / diag;
    const double half_width = std::sqrt(remaining) / diag;
    const double lower_d = center - half_width - 1e-9;
    const double upper_d = center + half_width + 1e-9;

    if (lower_d < static_cast<double>(WORD_MIN + 1) ||
        upper_d > static_cast<double>(WORD_MAX - 1)) {
        return false;
    }

    slong lower = static_cast<slong>(std::ceil(lower_d));
    slong upper = static_cast<slong>(std::floor(upper_d));
    if (state.max_coord >= 0) {
        if (lower < -state.max_coord) {
            lower = -state.max_coord;
        }
        if (upper > state.max_coord) {
            upper = state.max_coord;
        }
    }

    if (lower > upper) {
        return true;
    }

    slong center_coord = static_cast<slong>(std::floor(center + 0.5));
    if (center_coord < lower) {
        center_coord = lower;
    }
    if (center_coord > upper) {
        center_coord = upper;
    }

    const slong max_offset =
            FLINT_MAX(center_coord - lower, upper - center_coord);
    for (slong offset = 0; offset <= max_offset; ++offset) {
        const slong sides = offset == 0 ? 1 : 2;
        for (slong side = 0; side < sides; ++side) {
            const slong value =
                    side == 0 ? center_coord - offset : center_coord + offset;
            if (value < lower || value > upper) {
                continue;
            }
            if (state.aborted_by_callback) {
                return true;
            }

            const double y = diag * static_cast<double>(value) + delta;
            const double new_accumulated = accumulated + y * y;
            if (new_accumulated > state.bound) {
                continue;
            }

            state.partial[coordinate] = static_cast<double>(value);
            fmpz_set_si(fmpz_mat_entry(state.coeffs, 0, coordinate), value);

            if (coordinate == 0) {
                bool nonzero = false;
                for (slong i = 0; i < state.dimension; ++i) {
                    if (fmpz_is_zero(fmpz_mat_entry(state.coeffs, 0, i)) == 0) {
                        nonzero = true;
                        break;
                    }
                }
                if (nonzero && state.callback(state.coeffs, state.user) == 0) {
                    state.aborted_by_callback = true;
                    return true;
                }
            } else {
                if (!enum_recurse_double(state, coordinate - 1,
                            new_accumulated)) {
                    return false;
                }
                if (state.aborted_by_callback) {
                    return true;
                }
            }
        }
    }

    return true;
}

bool enum_downconvert_cholesky(double* out,
        const arb_mat_t cholesky,
        slong dimension) noexcept {
    flint::Arf mid;

    bool ok = true;
    for (slong i = 0; ok && i < dimension; ++i) {
        for (slong j = 0; j <= i; ++j) {
            const arb_struct* entry = arb_mat_entry(cholesky, i, j);
            flint::arf_set(mid, arb_midref(entry));
            const double value = flint::arf_get_d(mid, ARF_RND_NEAR);

            if (!std::isfinite(value)) {
                ok = false;
                break;
            }

            if (i == j) {
                if (value < lat_enum_double_min_diag ||
                    value > lat_enum_double_max_entry) {
                    ok = false;
                    break;
                }
            } else if (std::fabs(value) > lat_enum_double_max_entry) {
                ok = false;
                break;
            }

            out[i * dimension + j] = value;
        }
    }

    return ok;
}

bool enum_dispatch_double(const arb_mat_t cholesky,
        const arb_t bound_sq,
        slong dimension,
        slong max_coord,
        slong prec,
        ShortVectorCallback callback,
        void* user) noexcept {
    if (dimension > lat_enum_double_max_dim ||
        (max_coord >= 0 && max_coord > lat_enum_double_max_coord)) {
        return false;
    }

    std::vector<double> cholesky_d(
            static_cast<std::size_t>(dimension * dimension), 0.0);

    if (!enum_downconvert_cholesky(cholesky_d.data(), cholesky, dimension)) {
        return false;
    }

    flint::Arf upper_bound;
    flint::arb_get_ubound_arf(upper_bound, bound_sq, prec);
    if (flint::arf_is_nan(upper_bound) || flint::arf_is_inf(upper_bound)) {
        return false;
    }

    double bound = flint::arf_get_d(upper_bound, ARF_RND_UP);
    if (!std::isfinite(bound) || bound < 0.0) {
        return false;
    }

    bound *= 1.0 + std::ldexp(1.0, -40);

    std::vector<double> partial(static_cast<std::size_t>(dimension), 0.0);

    flint::FmpzMat coeffs(1, dimension);

    EnumStateDouble state;
    state.cholesky = cholesky_d.data();
    state.partial = partial.data();
    state.coeffs = coeffs.raw();
    state.dimension = dimension;
    state.bound = bound;
    state.max_coord = max_coord;
    state.callback = callback;
    state.user = user;

    return enum_recurse_double(state, dimension - 1, 0.0);
}

bool append_p_division_generators(fmpz_mat_t out,
        const fmpz_mat_t hnf,
        const fmpz_mod_mat_t nullspace,
        slong nullity,
        const fmpz_t p,
        const fmpz_mod_ctx_t ctx) noexcept {
    const slong r = fmpz_mat_nrows(hnf);
    const slong n = fmpz_mat_ncols(hnf);

    flint::FmpzMat tmp(r + nullity, n);
    flint::Fmpz coeff;
    flint::Fmpz lift;

    // Raw entries keep this p-saturation lifting kernel aligned with FLINT nullspace output.
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(tmp.raw(), i, j), fmpz_mat_entry(hnf, i, j));
        }
    }

    bool ok = true;
    for (slong k = 0; ok && k < nullity; ++k) {
        for (slong j = 0; ok && j < n; ++j) {
            fmpz_zero(lift.raw());
            for (slong i = 0; i < r; ++i) {
                fmpz_mod_mat_get_entry(coeff.raw(), nullspace, i, k, ctx);
                fmpz_addmul(lift.raw(), coeff.raw(), fmpz_mat_entry(hnf, i, j));
            }
            if (fmpz_divisible(lift.raw(), p) == 0) {
                ok = false;
            } else {
                fmpz_divexact(fmpz_mat_entry(tmp.raw(), r + k, j), lift.raw(), p);
            }
        }
    }

    if (ok) {
        fmpz_mat_swap(out, tmp.raw());
    }

    return ok;
}

}  // namespace

Lat::Lat(slong ambient_dim)
    : basis_(0, ambient_dim < 0 ? 0 : ambient_dim), ambient_dim_(ambient_dim) {
    if (ambient_dim < 0) {
        flint_throw(FLINT_ERROR, "Lat::Lat: negative ambient dimension\n");
    }
}

Lat::~Lat() noexcept {
    ambient_dim_ = 0;
    is_hnf_ = true;
}

void Lat::swap(Lat& other) noexcept {
    basis_.swap(other.basis_);
    std::swap(ambient_dim_, other.ambient_dim_);
    std::swap(is_hnf_, other.is_hnf_);
}

bool Lat::set(const Lat& other) noexcept {
    if (this == &other) {
        return true;
    }
    ambient_dim_ = other.ambient_dim_;
    set_basis_direct(other.basis_.raw(), other.is_hnf_);
    return true;
}

bool Lat::set_basis(flint::FmpzMatConstRef matrix) noexcept {
    if (matrix.raw() == nullptr || flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }
    set_basis_direct(matrix.raw(), fmpz_mat_is_in_hnf(matrix.raw()) != 0);
    return true;
}

bool Lat::get_basis(flint::FmpzMatRef matrix) const noexcept {
    if (matrix.raw() == nullptr ||
        flint::fmpz_mat_nrows(matrix) != flint::fmpz_mat_nrows(basis_) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }
    flint::fmpz_mat_set(matrix, basis_);
    return true;
}

flint::FmpzMat Lat::basis() const noexcept {
    flint::FmpzMat out(flint::fmpz_mat_nrows(basis_), ambient_dim_);
    get_basis(flint::FmpzMatRef(out));
    return out;
}

bool Lat::hnf(Lat& out) const noexcept {
    Lat tmp(ambient_dim_);

    if (is_hnf_) {
        const slong rows = nonzero_rows(basis_.raw());
        if (rows < 0) {
            return false;
        }
        flint::FmpzMat trimmed(rows, ambient_dim_);
        flint::FmpzMatConstWindow window(basis_, 0, 0, rows, ambient_dim_);
        flint::fmpz_mat_set(trimmed, window.const_ref());
        tmp.set_basis_direct(trimmed.raw(), true);
    } else {
        flint::FmpzMat hnf_matrix(flint::fmpz_mat_nrows(basis_), ambient_dim_);
        if (!hnf_basis(hnf_matrix.raw(), basis_.raw())) {
            return false;
        }
        tmp.set_basis_direct(hnf_matrix.raw(), true);
    }

    out.swap(tmp);
    return true;
}

bool Lat::hnf_transform(Lat& out, flint::FmpzMatRef transform) const noexcept {
    if (transform.raw() == nullptr ||
        flint::fmpz_mat_nrows(transform) != flint::fmpz_mat_nrows(basis_) ||
        flint::fmpz_mat_ncols(transform) != flint::fmpz_mat_nrows(basis_)) {
        return false;
    }

    Lat tmp(ambient_dim_);
    flint::FmpzMat hnf_matrix(flint::fmpz_mat_nrows(basis_), ambient_dim_);
    fmpz_mat_hnf_transform(hnf_matrix.raw(), transform.raw(), basis_.raw());
    tmp.set_basis_direct(hnf_matrix.raw(),
            fmpz_mat_is_in_hnf(hnf_matrix.raw()) != 0 &&
                    nonzero_rows(hnf_matrix.raw()) == flint::fmpz_mat_nrows(hnf_matrix));
    out.swap(tmp);
    return true;
}

bool Lat::contains(const fmpz* vector) const noexcept {
    if (vector == nullptr) {
        return false;
    }

    Lat hnf_lattice(ambient_dim_);
    Lat hnf_candidate(ambient_dim_);
    if (!hnf(hnf_lattice)) {
        return false;
    }

    flint::FmpzMat candidate(hnf_lattice.nrows() + 1, ambient_dim_);
    if (hnf_lattice.nrows() > 0) {
        flint::FmpzMatWindow top(
                candidate, 0, 0, hnf_lattice.nrows(), ambient_dim_);
        flint::fmpz_mat_set(top.ref(), hnf_lattice.basis_);
    }
    for (slong j = 0; j < ambient_dim_; ++j) {
        flint::fmpz_set(
                flint::fmpz_mat_entry(candidate, hnf_lattice.nrows(), j), vector + j);
    }
    hnf_candidate.set_basis(candidate);
    if (!hnf_candidate.hnf(hnf_candidate)) {
        return false;
    }
    return fmpz_mat_equal(hnf_lattice.basis_.raw(), hnf_candidate.basis_.raw()) != 0;
}

bool Lat::contains_row(flint::FmpzMatConstRef matrix, slong row) const noexcept {
    if (matrix.raw() == nullptr || row < 0 || row >= flint::fmpz_mat_nrows(matrix) ||
        flint::fmpz_mat_ncols(matrix) != ambient_dim_) {
        return false;
    }
    return contains(fmpz_mat_row(matrix.raw(), row));
}

bool Lat::sum(Lat& out, const Lat& other) const noexcept {
    if (!same_ambient(*this, other)) {
        return false;
    }

    Lat tmp(ambient_dim_);
    flint::FmpzMat combined(
            flint::fmpz_mat_nrows(basis_) + flint::fmpz_mat_nrows(other.basis_),
            ambient_dim_);
    if (flint::fmpz_mat_nrows(basis_) == 0) {
        flint::FmpzMatWindow bottom(combined, 0, 0,
                flint::fmpz_mat_nrows(other.basis_), ambient_dim_);
        flint::fmpz_mat_set(bottom.ref(), other.basis_);
    } else if (flint::fmpz_mat_nrows(other.basis_) == 0) {
        flint::FmpzMatWindow top(combined, 0, 0, flint::fmpz_mat_nrows(basis_),
                ambient_dim_);
        flint::fmpz_mat_set(top.ref(), basis_);
    } else {
        fmpz_mat_concat_vertical(combined.raw(), basis_.raw(), other.basis_.raw());
    }

    tmp.set_basis(combined);
    if (!tmp.hnf(tmp)) {
        return false;
    }
    out.swap(tmp);
    return true;
}

bool Lat::intersection(Lat& out, const Lat& other) const noexcept {
    if (!same_ambient(*this, other)) {
        return false;
    }

    Lat left_hnf(ambient_dim_);
    Lat right_hnf(ambient_dim_);
    Lat tmp(ambient_dim_);
    if (!hnf(left_hnf) || !other.hnf(right_hnf)) {
        return false;
    }

    const slong r = left_hnf.nrows();
    const slong s = right_hnf.nrows();
    const slong n = ambient_dim_;
    if (r == 0 || s == 0) {
        out.swap(tmp);
        return true;
    }

    flint::FmpzMat columns(n, r + s);
    flint::FmpzMat snf(n, r + s);
    flint::FmpzMat right_transform(r + s, r + s);

    // Matrix-entry loops here build the SNF input and recover the intersection basis.
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(columns.raw(), j, i),
                    fmpz_mat_entry(left_hnf.basis_.raw(), i, j));
        }
    }
    for (slong i = 0; i < s; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_neg(fmpz_mat_entry(columns.raw(), j, r + i),
                    fmpz_mat_entry(right_hnf.basis_.raw(), i, j));
        }
    }

    fmpz_mat_snf_transform(snf.raw(), nullptr, right_transform.raw(), columns.raw());
    slong rank = 0;
    while (rank < n && rank < r + s &&
            !flint::fmpz_is_zero(flint::fmpz_mat_entry(snf, rank, rank))) {
        ++rank;
    }

    const slong nullity = r + s - rank;
    flint::FmpzMat intersection_basis(nullity, n);
    for (slong k = 0; k < nullity; ++k) {
        for (slong i = 0; i < r; ++i) {
            for (slong j = 0; j < n; ++j) {
                fmpz_addmul(fmpz_mat_entry(intersection_basis.raw(), k, j),
                        fmpz_mat_entry(right_transform.raw(), i, rank + k),
                        fmpz_mat_entry(left_hnf.basis_.raw(), i, j));
            }
        }
    }

    tmp.set_basis(intersection_basis);

    if (!tmp.hnf(tmp)) {
        return false;
    }
    out.swap(tmp);
    return true;
}

bool Lat::index(flint::FmpzRef index, const Lat& sublattice) const noexcept {
    if (index.raw() == nullptr || !same_ambient(*this, sublattice)) {
        return false;
    }

    Lat ambient_hnf(ambient_dim_);
    Lat sub_hnf(ambient_dim_);
    if (!hnf(ambient_hnf) || !sublattice.hnf(sub_hnf)) {
        return false;
    }

    const slong r = ambient_hnf.nrows();
    const slong n = ambient_dim_;
    if (r != sub_hnf.nrows()) {
        return false;
    }
    if (r == 0) {
        flint::fmpz_one(index);
        return true;
    }
    for (slong i = 0; i < r; ++i) {
        if (!ambient_hnf.contains(fmpz_mat_row(sub_hnf.basis_.raw(), i))) {
            return false;
        }
    }

    flint::FmpzMat left(n, r);
    flint::FmpzMat right(n, r);
    flint::FmpzMat solution(r, r);
    flint::FmpzMat integral_solution(r, r);
    flint::Fmpz denominator;
    flint::Fmpz det;

    // Transposed solve inputs are assembled entrywise for the FLINT linear solve.
    for (slong i = 0; i < r; ++i) {
        for (slong j = 0; j < n; ++j) {
            fmpz_set(fmpz_mat_entry(left.raw(), j, i),
                    fmpz_mat_entry(ambient_hnf.basis_.raw(), i, j));
            fmpz_set(fmpz_mat_entry(right.raw(), j, i),
                    fmpz_mat_entry(sub_hnf.basis_.raw(), i, j));
        }
    }

    bool ok = false;
    if (fmpz_mat_can_solve(solution.raw(), denominator.raw(), left.raw(),
                right.raw()) != 0 &&
        integral_divexact_matrix(integral_solution.raw(), solution.raw(),
                denominator.raw())) {
        fmpz_mat_det(det.raw(), integral_solution.raw());
        if (!flint::fmpz_is_zero(det)) {
            flint::fmpz_abs(index, det);
            ok = true;
        }
    }

    return ok;
}

bool Lat::saturate(Lat& out, flint::FmpzConstRef p) const noexcept {
    if (p.raw() == nullptr || !flint::fmpz_is_prime(p)) {
        return false;
    }

    Lat tmp(ambient_dim_);
    if (!hnf(tmp)) {
        return false;
    }

    silex::flint::FmpzModCtx ctx(p.raw());
    flint::FmpzMat lifted(0, ambient_dim_);

    bool ok = true;
    for (;;) {
        const slong r = tmp.nrows();
        const slong n = ambient_dim_;
        if (r == 0) {
            break;
        }

        silex::flint::FmpzModMat hmod(r, n, ctx);
        silex::flint::FmpzModMat transpose(n, r, ctx);
        silex::flint::FmpzModMat nullspace(r, r, ctx);

        fmpz_mod_mat_set_fmpz_mat(hmod.raw(), tmp.basis_.raw(), ctx.raw());
        fmpz_mod_mat_transpose(transpose.raw(), hmod.raw(), ctx.raw());
        const slong nullity =
                fmpz_mod_mat_nullspace(nullspace.raw(), transpose.raw(), ctx.raw());
        if (nullity == 0) {
            break;
        }

        lifted = flint::FmpzMat(0, n);
        if (!append_p_division_generators(lifted.raw(), tmp.basis_.raw(),
                    nullspace.raw(), nullity, p.raw(), ctx.raw()) ||
                !tmp.set_basis(lifted) || !tmp.hnf(tmp)) {
            ok = false;
            break;
        }
    }

    if (!ok) {
        return false;
    }

    out.swap(tmp);
    return true;
}

bool Lat::lll_reduce(Lat& out) const noexcept {
    Lat hnf_lattice(ambient_dim_);
    Lat tmp(ambient_dim_);
    if (!hnf(hnf_lattice)) {
        return false;
    }

    const slong r = hnf_lattice.nrows();
    if (r == 0) {
        out.swap(tmp);
        return true;
    }

    flint::FmpzMat reduced(flint::fmpz_mat_nrows(hnf_lattice.basis_),
            flint::fmpz_mat_ncols(hnf_lattice.basis_));
    flint::FmpzMat transform(r, r);
    fmpz_mat_set(reduced.raw(), hnf_lattice.basis_.raw());

    silex::flint::FmpzLll config;
    fmpz_lll(reduced.raw(), transform.raw(), config.raw());
    tmp.set_basis_direct(reduced.raw(), false);

    out.swap(tmp);
    return true;
}

bool Lat::enum_short_vectors_arb(flint::ArbConstRef bound_sq,
        slong max_coord,
        slong prec,
        ShortVectorCallback callback,
        void* user) const noexcept {
    if (bound_sq.raw() == nullptr || callback == nullptr || prec < 2) {
        return false;
    }
    if (flint::arb_is_negative(bound_sq.raw())) {
        return true;
    }

    const slong row_count = flint::fmpz_mat_nrows(basis_);
    const slong col_count = flint::fmpz_mat_ncols(basis_);
    if (row_count == 0) {
        return true;
    }

    flint::ArbMat gram(row_count, row_count);
    flint::ArbMat cholesky(row_count, row_count);
    flint::Arb entry;
    flint::Arb accumulated;

    for (slong i = 0; i < row_count; ++i) {
        for (slong j = i; j < row_count; ++j) {
            flint::arb_zero(accumulated);
            for (slong k = 0; k < col_count; ++k) {
                flint::arb_set_fmpz(entry, fmpz_mat_entry(basis_.raw(), i, k));
                flint::arb_mul_fmpz(entry, entry,
                        fmpz_mat_entry(basis_.raw(), j, k), prec);
                flint::arb_add(accumulated, accumulated, entry, prec);
            }
            flint::arb_set(arb_mat_entry(gram.raw(), i, j), accumulated);
            if (i != j) {
                flint::arb_set(arb_mat_entry(gram.raw(), j, i), accumulated);
            }
        }
    }

    if (flint::arb_mat_cho(cholesky, gram, prec) == 0) {
        return false;
    }

    if (enum_dispatch_double(cholesky.raw(), bound_sq.raw(), row_count, max_coord,
                prec, callback, user)) {
        return true;
    }

    flint::FmpzMat coeffs(1, row_count);

    EnumStateArb state;
    state.cholesky = cholesky.raw();
    state.bound_sq = bound_sq.raw();
    state.coeffs = coeffs.raw();
    state.max_coord = max_coord;
    state.prec = prec;
    state.callback = callback;
    state.user = user;

    flint::arb_zero(accumulated);
    const bool ok = enum_recurse_arb(state, row_count - 1, accumulated.raw());
    return ok && (!state.precision_failure || state.aborted_by_callback);
}

bool Lat::check(const DiagnosticsContext* diagnostics) const noexcept {
#if defined(SILEX_ENABLE_DEBUG_CHECKS) && SILEX_ENABLE_DEBUG_CHECKS
    if (!debug_check_enabled(diagnostics, DiagnosticsModule::lattice,
                DebugLevel::cheap)) {
        return true;
    }

    if (ambient_dim_ < 0) {
        return lat_check_failed(diagnostics, DebugLevel::cheap, __func__,
                __FILE__, __LINE__, "ambient dimension",
                "ambient_dim_ >= 0");
    }

    if (flint::fmpz_mat_ncols(basis_) != ambient_dim_) {
        return lat_check_failed(diagnostics, DebugLevel::cheap, __func__,
                __FILE__, __LINE__, "basis ambient dimension",
                "fmpz_mat_ncols(basis_) == ambient_dim_");
    }

    if (!debug_check_enabled(diagnostics, DiagnosticsModule::lattice,
                DebugLevel::normal) ||
        !is_hnf_) {
        return true;
    }

    flint::FmpzMat hnf_basis(flint::fmpz_mat_nrows(basis_),
            flint::fmpz_mat_ncols(basis_));
    fmpz_mat_hnf(hnf_basis.raw(), basis_.raw());
    const bool ok = fmpz_mat_equal(hnf_basis.raw(), basis_.raw()) != 0;
    if (!ok) {
        return lat_check_failed(diagnostics, DebugLevel::normal, __func__,
                __FILE__, __LINE__, "hnf marker",
                "is_hnf_ implies basis_ is in HNF");
    }
#else
    (void)diagnostics;
#endif
    return true;
}

void Lat::set_basis_direct(const fmpz_mat_t matrix, bool is_hnf) noexcept {
    flint::FmpzMat tmp(fmpz_mat_nrows(matrix), fmpz_mat_ncols(matrix));
    fmpz_mat_set(tmp.raw(), matrix);
    basis_.swap(tmp);
    is_hnf_ = is_hnf;
}

}  // namespace silex::lat
