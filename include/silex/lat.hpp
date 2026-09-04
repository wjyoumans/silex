#pragma once

#include <flint/arb.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <silex/flint/arb.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>

#include <utility>

namespace silex {
struct DiagnosticsContext;
}

namespace silex::lat {

// Low-level callback shape for enumeration kernels. The coefficient matrix is
// borrowed for the duration of the callback.
using ShortVectorCallback = int (*)(const fmpz_mat_t coeffs, void* user);

class Lat {
public:
    explicit Lat(slong ambient_dim);
    ~Lat() noexcept;

    Lat(const Lat&) = delete;
    Lat& operator=(const Lat&) = delete;

    Lat(Lat&& other) noexcept : Lat(0) { swap(other); }
    Lat& operator=(Lat&& other) noexcept {
        if (this != &other) {
            Lat tmp(0);
            tmp.swap(*this);
            swap(other);
        }
        return *this;
    }

    void swap(Lat& other) noexcept;
    bool set(const Lat& other) noexcept;
    bool set_basis(flint::FmpzMatConstRef matrix) noexcept;
    bool get_basis(flint::FmpzMatRef matrix) const noexcept;
    flint::FmpzMat basis() const noexcept;

    slong ambient_dim() const noexcept { return ambient_dim_; }
    slong nrows() const noexcept { return fmpz_mat_nrows(basis_.raw()); }
    bool is_hnf() const noexcept { return is_hnf_; }
    flint::FmpzMatConstRef basis_ref() const noexcept { return basis_; }
    // Low-level FLINT interop for kernels and bridge code; prefer basis_ref()
    // or get_basis() when direct raw FLINT mutation is not required.
    const fmpz_mat_struct* raw_basis() const noexcept { return basis_.raw(); }
    fmpz_mat_struct* raw_basis() noexcept { return basis_.raw(); }

    bool hnf(Lat& out) const noexcept;
    bool hnf_transform(Lat& out, flint::FmpzMatRef transform) const noexcept;
    bool contains(const fmpz* vector) const noexcept;
    bool contains_row(flint::FmpzMatConstRef matrix, slong row) const noexcept;
    bool sum(Lat& out, const Lat& other) const noexcept;
    bool intersection(Lat& out, const Lat& other) const noexcept;
    bool index(flint::FmpzRef index, const Lat& sublattice) const noexcept;
    bool saturate(Lat& out, flint::FmpzConstRef p) const noexcept;
    bool lll_reduce(Lat& out) const noexcept;
    bool enum_short_vectors_arb(flint::ArbConstRef bound_sq,
            slong max_coord,
            slong prec,
            ShortVectorCallback callback,
            void* user) const noexcept;
    bool check(const DiagnosticsContext* diagnostics) const noexcept;

private:
    void set_basis_direct(const fmpz_mat_t matrix, bool is_hnf) noexcept;

    flint::FmpzMat basis_;
    slong ambient_dim_ = 0;
    bool is_hnf_ = true;
};

inline void swap(Lat& left, Lat& right) noexcept { left.swap(right); }

}  // namespace silex::lat
