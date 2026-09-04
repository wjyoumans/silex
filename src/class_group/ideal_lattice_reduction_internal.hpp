#pragma once

#include <silex/diagnostics.hpp>
#include <silex/element.hpp>
#include <silex/factored_element.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>

namespace silex::detail {

struct OrderMinkowskiEmbeddingCache;

bool ideal_lattice_short_element(
        Element& out,
        const Ideal& ideal,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool weighted_ideal_lattice_short_element(
        Element& out,
        const Ideal& ideal,
        flint::FmpzMatConstRef weights,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr,
        OrderMinkowskiEmbeddingCache* cache = nullptr) noexcept;

bool weighted_ideal_lattice_short_element(
        Element& out,
        const FractionalIdeal& ideal,
        flint::FmpzMatConstRef weights,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr,
        OrderMinkowskiEmbeddingCache* cache = nullptr) noexcept;

bool reduce_ideal_lattice(
        Ideal& reduced,
        Element& multiplier,
        const Ideal& ideal,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool reduce_ideal_product(
        Ideal& reduced,
        Element& multiplier,
        const Ideal& left,
        const Ideal& right,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

bool reduce_ideal_signed_power(
        Ideal& reduced,
        FactoredElement& multiplier,
        const Ideal& ideal,
        slong exponent,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;

}  // namespace silex::detail
