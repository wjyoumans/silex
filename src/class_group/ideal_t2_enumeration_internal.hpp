#pragma once

#include <silex/class_group.hpp>
#include <silex/flint/arb_mat.hpp>

#include <vector>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;

}  // namespace silex::detail

namespace silex::detail::relation_search {

inline constexpr slong kMaxFactorAttempts = 500;
inline constexpr slong kMaxElementSteps =
        4 * kMaxFactorAttempts * kMaxFactorAttempts;

struct FiniteIdealT2EnumerationData {
    flint::FmpzMat basis{0, 0};
    std::vector<double> quadratic_form_data;
    double initial_bound_value = 0.0;
};

class FiniteQuadraticFormEnumerationContext {
public:
    bool reset(const std::vector<double>& quadratic_form_data, slong dimension) noexcept;
    bool start(double bound, slong max_element_steps,
               bool skip_first_scalar) noexcept;
    bool next() noexcept;
    bool current_row(flint::FmpzMatRef out) const noexcept;
    slong element_steps() const noexcept;
    bool element_step_limit_reached() const noexcept;

private:
    void clear() noexcept;
    double q(slong row, slong column) const noexcept;
    double v(slong index) const noexcept;
    bool step(slong index) noexcept;

    slong dimension_ = 0;
    std::vector<double> quadratic_form_data_;
    std::vector<slong> coefficients_;
    std::vector<slong> increments_;
    std::vector<double> lengths_;
    std::vector<double> tails_;
    slong level_ = 0;
    slong max_element_steps_ = 0;
    slong skip_first_ = 0;
    slong element_steps_ = 0;
    double bound_ = 0.0;
    bool started_ = false;
    bool exhausted_ = true;
    bool need_outer_step_ = false;
    bool element_step_limit_reached_ = false;
};

bool fmpz_mat_single_row_is_primitive(flint::FmpzMatConstRef row) noexcept;

bool build_finite_ideal_t2_enumeration_data_with_retry(
        FiniteIdealT2EnumerationData& out, const Ideal& ideal,
        const DiagnosticsContext* diagnostics,
        detail::OrderMinkowskiEmbeddingCache* embedding_cache) noexcept;

double arb_midpoint_double(flint::ArbConstRef value) noexcept;

}  // namespace silex::detail::relation_search
