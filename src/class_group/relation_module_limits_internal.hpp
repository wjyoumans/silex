#pragma once

#include <flint/flint.h>

namespace silex::detail {

// Largest factor base, in generators, for which a class-group context
// builds its relation modules: 2^14.  Each reset allocates three dense
// ncols x ncols mod-p matrices (the row module's and the relation module's
// fmpz_smat::ModRankContext and the relation module's checkpoint, FLINT
// nmod_mat), 2 GiB each at this limit, and FLINT aborts the process when
// such an allocation fails (flint_calloc).  A larger base fails closed
// before allocating, with reason factor_base_size_limit_exceeded in a
// class/unit transaction, so it fails the same way on every machine.  The
// largest base of the test suite and the benchmark corpus has fewer than
// 400 generators; the proven base of an imaginary quadratic field, up to
// sqrt(|D| / 3), passes the limit between |D| = 4e11 and 5e11 (about 2.4e4
// generators at |D| = 1e12).  Past it, a `proven` imaginary quadratic
// request keeps its compact relation base and is certified by the class
// group's structure certificate instead (relation_factor_base_plan.cpp).
inline constexpr slong kRelationModuleMaxColumns = WORD(1) << 14;

}  // namespace silex::detail
