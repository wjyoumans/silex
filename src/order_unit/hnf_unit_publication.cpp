#include "order_unit_internal.hpp"

namespace silex::detail {

bool publish_validated_hnf_units(
        OrderUnitGroup& out,
        const Order& order,
        std::vector<FactoredElement>& units,
        EmbeddingContext& embeddings,
        flint::ArbConstRef expected_regulator,
        slong rank,
        slong precision,
        const flint::Fmpz* cached_torsion_order,
        const OrderElement* cached_torsion_generator) noexcept {
    const DiagnosticsContext* const diagnostics = out.diagnostics();
    if (static_cast<slong>(units.size()) != rank) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit publication failed: rank mismatch");
        return false;
    }

    OrderUnitGroup candidate(order);
    candidate.set_diagnostics(diagnostics);
    bool units_set = false;
    {
        SILEX_PROFILE_SCOPE(candidate.diagnostics(),
                            DiagnosticsModule::unit_group,
                            "unit_group.set_units");
        units_set = order_unit_group_set_units_internal(
                candidate, order,
                FactoredElementSpan(units.data(), units.size()), embeddings,
                precision, false, cached_torsion_order,
                cached_torsion_generator);
    }
    if (!units_set) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit publication failed: exact unit set");
        return false;
    }
    if (!unit_regulator_matches_reconstruction(
                candidate, expected_regulator, precision)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::unit_group,
                  LogLevel::detail,
                  "HNF unit publication failed: regulator mismatch");
        return false;
    }

    out.swap(candidate);
    return true;
}

}  // namespace silex::detail
