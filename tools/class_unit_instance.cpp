#include <silex/class_group.hpp>
#include <silex/diagnostics.hpp>
#include <silex/embedding.hpp>
#include <silex/factor_base.hpp>
#include <silex/element.hpp>
#include <silex/flint/fmpq.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>
#include <silex/order_unit.hpp>
#include <silex/prime_ideal.hpp>
#include <silex/signature.hpp>
#include <silex/sunit.hpp>

#include "order_unit/class_unit_transaction_internal.hpp"

#include <flint/arf.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {
namespace sflint = silex::flint;

struct ProfileAggregate {
    std::uint64_t count = 0;
    std::uint64_t events = 0;
    std::chrono::nanoseconds inclusive{0};
    std::chrono::nanoseconds exclusive{0};
};

struct ActiveProfileScope {
    std::string key;
    std::chrono::steady_clock::time_point start;
    std::chrono::nanoseconds child_time{0};
};

struct SPrimeSelector {
    slong rational_prime = 0;
    slong canonical_index = -1;
    slong selection_index = -1;
    std::vector<std::string> beta_power_basis;
};

struct ProfileCollector {
    std::unordered_map<std::string, ProfileAggregate> aggregates;
    std::vector<ActiveProfileScope> stack;

    static std::string key(silex::DiagnosticsModule module,
                           const char* label) {
        std::string out = silex::diagnostics_module_name(module);
        out += ":";
        out += label == nullptr ? "<null>" : label;
        return out;
    }

    static void callback(void* user,
                         silex::DiagnosticsModule module,
                         silex::ProfileEvent event,
                         const char*,
                         const char* label) noexcept {
        auto* collector = static_cast<ProfileCollector*>(user);
        if (collector == nullptr) {
            return;
        }

        const std::string profile_key = key(module, label);
        if (event == silex::ProfileEvent::begin_scope) {
            collector->stack.push_back(ActiveProfileScope{
                    profile_key, std::chrono::steady_clock::now(), {}});
            return;
        }
        if (event == silex::ProfileEvent::end_scope) {
            if (collector->stack.empty()) {
                return;
            }
            ActiveProfileScope scope = std::move(collector->stack.back());
            collector->stack.pop_back();
            const auto elapsed =
                    std::chrono::steady_clock::now() - scope.start;
            auto& aggregate = collector->aggregates[scope.key];
            ++aggregate.count;
            aggregate.inclusive += elapsed;
            aggregate.exclusive += elapsed - scope.child_time;
            if (!collector->stack.empty()) {
                collector->stack.back().child_time += elapsed;
            }
            return;
        }

        auto& aggregate = collector->aggregates[profile_key];
        ++aggregate.events;
    }

    void configure(silex::DiagnosticsContext& diagnostics) noexcept {
        const auto modules =
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::field) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::element) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::order) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::ideal) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::fractional_ideal) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::prime_ideal) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::lattice) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::residue) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::relation) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::class_group) |
                silex::diagnostics_module_bit(
                        silex::DiagnosticsModule::unit_group);
        silex::diagnostics_set_profiling(
                diagnostics, true, modules, callback, this);
    }

    std::vector<std::pair<std::string, ProfileAggregate>> sorted_rows() const {
        std::vector<std::pair<std::string, ProfileAggregate>> rows;
        rows.reserve(aggregates.size());
        for (const auto& row : aggregates) {
            rows.push_back(row);
        }
        std::sort(rows.begin(), rows.end(), [](const auto& left,
                                               const auto& right) {
            return left.second.inclusive > right.second.inclusive;
        });
        return rows;
    }
};

struct Options {
    std::vector<slong> coefficients;
    std::vector<slong> warmup_coefficients;
    std::string mode = "proven";
    slong precision = 128;
    slong max_candidates = -1;
    slong max_relations = -1;
    ulong zeta_bf_max_cutoff = 20000;
    std::optional<slong> factor_base_bound_override;
    std::optional<slong> expect_class_order;
    std::optional<slong> expect_unit_rank;
    std::optional<bool> expect_success;
    std::vector<SPrimeSelector> s_prime_selectors;
    bool compute_sunit = false;
    bool logging = false;
    bool trace = false;
    bool verbose = false;
    bool profiling = false;
    bool marked_protocol = false;
};

constexpr const char* ready_marker = "__SILEX_BENCH_SILEX_READY__";
constexpr const char* target_done_marker =
        "__SILEX_BENCH_SILEX_TARGET_DONE__";
constexpr std::size_t max_protocol_phase_bytes = 1U << 20U;
constexpr std::size_t target_nonce_bytes = 32;

void write_json_string(std::ostream& out, std::string_view value) {
    out << '"';
    for (const char ch : value) {
        switch (ch) {
        case '"':
            out << "\\\"";
            break;
        case '\\':
            out << "\\\\";
            break;
        case '\b':
            out << "\\b";
            break;
        case '\f':
            out << "\\f";
            break;
        case '\n':
            out << "\\n";
            break;
        case '\r':
            out << "\\r";
            break;
        case '\t':
            out << "\\t";
            break;
        default:
            if (static_cast<unsigned char>(ch) < 0x20) {
                out << "\\u00";
                const char* hex = "0123456789abcdef";
                out << hex[(static_cast<unsigned char>(ch) >> 4) & 0xf]
                    << hex[static_cast<unsigned char>(ch) & 0xf];
            } else {
                out << ch;
            }
            break;
        }
    }
    out << '"';
}

const char* json_bool(bool value) noexcept {
    return value ? "true" : "false";
}

const char* certification_name(silex::CertificationMode mode) noexcept {
    switch (mode) {
    case silex::CertificationMode::unknown:
        return "unknown";
    case silex::CertificationMode::heuristic:
        return "heuristic";
    case silex::CertificationMode::grh:
        return "grh";
    case silex::CertificationMode::proven:
        return "proven";
    }
    return "unknown";
}

const char* proof_state_name(silex::ProofState state) noexcept {
    switch (state) {
    case silex::ProofState::not_checked:
        return "not_checked";
    case silex::ProofState::unavailable:
        return "unavailable";
    case silex::ProofState::verified:
        return "verified";
    }
    return "not_checked";
}

const char* membership_outcome_name(
        silex::SUnitMembershipOutcome outcome) noexcept {
    switch (outcome) {
    case silex::SUnitMembershipOutcome::unknown:
        return "unknown";
    case silex::SUnitMembershipOutcome::not_sunit:
        return "not_sunit";
    case silex::SUnitMembershipOutcome::verified:
        return "verified";
    }
    return "unknown";
}

bool compute_class_unit_from_options(
        silex::OrderUnitGroup& units,
        silex::ClassGroupContext& class_group,
        const silex::Order& order,
        sflint::FmpzConstRef factor_base_bound,
        const silex::ClassGroupComputeOptions& compute_options,
        slong precision,
        const Options&,
        silex::detail::ClassUnitTransactionReport& audit) noexcept {
    return silex::detail::compute_class_unit_transaction(
            units, class_group, order, factor_base_bound, compute_options,
            precision, audit);
}

std::string fmpz_string(sflint::FmpzConstRef value) {
    char* raw = fmpz_get_str(nullptr, 10, value.raw());
    if (raw == nullptr) {
        return {};
    }
    std::string out(raw);
    flint_free(raw);
    return out;
}

std::string fmpq_string(sflint::FmpqConstRef value) {
    char* raw = fmpq_get_str(nullptr, 10, value.raw());
    if (raw == nullptr) {
        return {};
    }
    std::string out(raw);
    flint_free(raw);
    return out;
}

bool parse_canonical_rational_text(const std::string& value,
                                   sflint::Fmpq* out = nullptr) noexcept {
    constexpr std::size_t max_rational_text_bytes = 8192;
    if (value.empty() || value.size() > max_rational_text_bytes) {
        return false;
    }
    sflint::Fmpq parsed;
    if (::fmpq_set_str(parsed.raw(), value.c_str(), 10) != 0 ||
        fmpq_string(sflint::FmpqConstRef(parsed)) != value) {
        return false;
    }
    if (out != nullptr) {
        sflint::fmpq_set(*out, parsed);
    }
    return true;
}

bool parse_slong_value(const char* value, slong& out) noexcept {
    if (value == nullptr || *value == '\0') {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    const long long parsed = std::strtoll(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0') {
        return false;
    }
    if (parsed < static_cast<long long>(std::numeric_limits<slong>::min()) ||
        parsed > static_cast<long long>(std::numeric_limits<slong>::max())) {
        return false;
    }
    out = static_cast<slong>(parsed);
    return true;
}

bool parse_ulong_value(const char* value, ulong& out) noexcept {
    if (value == nullptr || *value == '\0' || value[0] == '-') {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    const unsigned long long parsed = std::strtoull(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0') {
        return false;
    }
    if (parsed >
        static_cast<unsigned long long>(std::numeric_limits<ulong>::max())) {
        return false;
    }
    out = static_cast<ulong>(parsed);
    return true;
}

bool parse_bool_value(const char* value, bool& out) noexcept {
    if (value == nullptr) {
        return false;
    }
    if (std::strcmp(value, "1") == 0 ||
        std::strcmp(value, "true") == 0 ||
        std::strcmp(value, "yes") == 0) {
        out = true;
        return true;
    }
    if (std::strcmp(value, "0") == 0 ||
        std::strcmp(value, "false") == 0 ||
        std::strcmp(value, "no") == 0) {
        out = false;
        return true;
    }
    return false;
}

bool parse_s_prime_witness(const std::string& value,
                           SPrimeSelector& out) noexcept {
    const std::size_t first_colon = value.find(':');
    const std::size_t second_colon =
            first_colon == std::string::npos
                    ? std::string::npos
                    : value.find(':', first_colon + 1);
    const std::size_t third_colon =
            second_colon == std::string::npos
                    ? std::string::npos
                    : value.find(':', second_colon + 1);
    if (first_colon == std::string::npos ||
        second_colon == std::string::npos ||
        third_colon == std::string::npos ||
        value.find(':', third_colon + 1) != std::string::npos) {
        return false;
    }
    slong rational_prime = 0;
    slong canonical_index = -1;
    if (!parse_slong_value(
                value.substr(0, first_colon).c_str(), rational_prime) ||
        !parse_slong_value(
                value.substr(first_colon + 1,
                             second_colon - first_colon - 1)
                        .c_str(),
                canonical_index) ||
        rational_prime < 2) {
        return false;
    }
    slong selection_index = -1;
    if (canonical_index < 0 ||
        !parse_slong_value(
                value.substr(second_colon + 1,
                             third_colon - second_colon - 1)
                        .c_str(),
                selection_index) ||
        selection_index < -1) {
        return false;
    }
    out = {};
    out.rational_prime = rational_prime;
    out.canonical_index = canonical_index;
    out.selection_index = selection_index;
    const std::string coefficients = value.substr(third_colon + 1);
    std::size_t start = 0;
    while (start <= coefficients.size()) {
        const std::size_t comma = coefficients.find(',', start);
        const std::size_t end =
                comma == std::string::npos ? coefficients.size() : comma;
        const std::string coefficient =
                coefficients.substr(start, end - start);
        if (!parse_canonical_rational_text(coefficient)) {
            return false;
        }
        out.beta_power_basis.push_back(coefficient);
        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return !out.beta_power_basis.empty();
}

bool parse_coefficients(const std::string& value,
                        std::vector<slong>& out) noexcept {
    out.clear();
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t comma = value.find(',', start);
        const std::size_t end =
                comma == std::string::npos ? value.size() : comma;
        const std::string token = value.substr(start, end - start);
        slong coeff = 0;
        if (!parse_slong_value(token.c_str(), coeff)) {
            return false;
        }
        out.push_back(coeff);
        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return out.size() >= 2;
}

void print_usage(std::ostream& out) {
    out << "Usage: silex-class-unit-instance --coeffs c0,c1,...,1 "
           "[options]\n"
        << "\n"
        << "The coefficient list is low-to-high and must currently be "
           "monic.\n"
        << "\n"
        << "Options:\n"
        << "  --mode proven|grh\n"
        << "  --warmup-coeffs c0,c1,...,1\n"
        << "  --precision N\n"
        << "  --max-candidates N\n"
        << "  --max-relations N\n"
        << "  --bf-cutoff N\n"
        << "  --factor-base-bound N\n"
        << "  --expect-class-order N\n"
        << "  --expect-unit-rank N\n"
        << "  --expect-success 0|1\n"
        << "  --compute-sunit\n"
        << "  --s-prime-witness P:INDEX:SELECTION_INDEX:BETA0,...,BETAn "
           "(repeatable; complete manifest witness order)\n"
        << "  --marked-protocol\n"
        << "  --log --trace --verbose --profile\n";
}

bool take_value(int& index,
                int argc,
                char** argv,
                const std::string& arg,
                const char* name,
                std::string& value) {
    const std::string prefix = std::string(name) + "=";
    if (arg == name) {
        if (index + 1 >= argc) {
            return false;
        }
        value = argv[++index];
        return true;
    }
    if (arg.rfind(prefix, 0) == 0) {
        value = arg.substr(prefix.size());
        return true;
    }
    return false;
}

bool parse_options(int argc,
                   char** argv,
                   Options& options,
                   std::string& error) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        std::string value;
        if (arg == "--help" || arg == "-h") {
            print_usage(std::cout);
            std::exit(0);
        }
        if (take_value(i, argc, argv, arg, "--coeffs", value)) {
            if (!parse_coefficients(value, options.coefficients)) {
                error = "invalid --coeffs value";
                return false;
            }
            continue;
        }
        if (take_value(i, argc, argv, arg, "--warmup-coeffs", value)) {
            if (!parse_coefficients(value, options.warmup_coefficients)) {
                error = "invalid --warmup-coeffs value";
                return false;
            }
            continue;
        }
        if (take_value(i, argc, argv, arg, "--mode", value)) {
            options.mode = value;
            continue;
        }
        if (take_value(i, argc, argv, arg, "--precision", value)) {
            if (!parse_slong_value(value.c_str(), options.precision) ||
                options.precision <= 0) {
                error = "invalid --precision value";
                return false;
            }
            continue;
        }
        if (take_value(i, argc, argv, arg, "--max-candidates", value)) {
            if (!parse_slong_value(value.c_str(), options.max_candidates) ||
                options.max_candidates < 0) {
                error = "invalid --max-candidates value";
                return false;
            }
            continue;
        }
        if (take_value(i, argc, argv, arg, "--max-relations", value)) {
            if (!parse_slong_value(value.c_str(), options.max_relations) ||
                options.max_relations < 0) {
                error = "invalid --max-relations value";
                return false;
            }
            continue;
        }
        if (take_value(i, argc, argv, arg, "--bf-cutoff", value)) {
            if (!parse_ulong_value(value.c_str(),
                                   options.zeta_bf_max_cutoff)) {
                error = "invalid --bf-cutoff value";
                return false;
            }
            continue;
        }
        if (take_value(i, argc, argv, arg, "--factor-base-bound", value)) {
            slong bound = 0;
            if (!parse_slong_value(value.c_str(), bound) || bound < 2) {
                error = "invalid --factor-base-bound value";
                return false;
            }
            options.factor_base_bound_override = bound;
            continue;
        }
        if (take_value(i, argc, argv, arg, "--expect-class-order", value)) {
            slong expected = 0;
            if (!parse_slong_value(value.c_str(), expected)) {
                error = "invalid --expect-class-order value";
                return false;
            }
            options.expect_class_order = expected;
            continue;
        }
        if (take_value(i, argc, argv, arg, "--expect-unit-rank", value)) {
            slong expected = 0;
            if (!parse_slong_value(value.c_str(), expected)) {
                error = "invalid --expect-unit-rank value";
                return false;
            }
            options.expect_unit_rank = expected;
            continue;
        }
        if (take_value(i, argc, argv, arg, "--expect-success", value)) {
            bool expected = false;
            if (!parse_bool_value(value.c_str(), expected)) {
                error = "invalid --expect-success value";
                return false;
            }
            options.expect_success = expected;
            continue;
        }
        if (take_value(i, argc, argv, arg, "--s-prime-witness", value)) {
            SPrimeSelector selector;
            if (!parse_s_prime_witness(value, selector)) {
                error = "invalid --s-prime-witness value";
                return false;
            }
            options.s_prime_selectors.push_back(selector);
            options.compute_sunit = true;
            continue;
        }
        if (arg == "--compute-sunit") {
            options.compute_sunit = true;
            continue;
        }
        if (arg == "--log") {
            options.logging = true;
            continue;
        }
        if (arg == "--trace") {
            options.trace = true;
            continue;
        }
        if (arg == "--verbose") {
            options.verbose = true;
            continue;
        }
        if (arg == "--profile") {
            options.profiling = true;
            continue;
        }
        if (arg == "--marked-protocol") {
            options.marked_protocol = true;
            continue;
        }
        error = "unknown option: " + arg;
        return false;
    }

    if (options.coefficients.empty()) {
        error = "missing --coeffs";
        return false;
    }
    if (options.coefficients.back() != 1) {
        error = "only monic input is supported by this diagnostic runner";
        return false;
    }
    if (!options.warmup_coefficients.empty() &&
        options.warmup_coefficients.back() != 1) {
        error = "only monic warmup input is supported by this diagnostic runner";
        return false;
    }
    if (!options.warmup_coefficients.empty() &&
        options.warmup_coefficients.size() != options.coefficients.size()) {
        error = "warmup and target degrees differ";
        return false;
    }
    if (!options.warmup_coefficients.empty() &&
        options.warmup_coefficients == options.coefficients) {
        error = "warmup and target polynomials must differ";
        return false;
    }
    if (options.mode != "proven" && options.mode != "grh") {
        error = "invalid --mode value";
        return false;
    }
    return true;
}

void diagnostic_log_callback(void*,
                             silex::DiagnosticsModule module,
                             silex::LogLevel level,
                             const char* function,
                             const char* message,
                             const char* detail) noexcept {
    std::cerr << "[log "
              << silex::diagnostics_module_name(module) << ":"
              << silex::log_level_name(level) << "] ";
    if (function != nullptr) {
        std::cerr << function << ": ";
    }
    std::cerr << (message == nullptr ? "" : message);
    if (detail != nullptr) {
        std::cerr << ": " << detail;
    }
    std::cerr << "\n";
}

void diagnostic_verbose_callback(void*,
                                 silex::DiagnosticsModule module,
                                 silex::VerboseLevel level,
                                 const char* function,
                                 const char* message,
                                 const char* detail) noexcept {
    std::cerr << "[verbose "
              << silex::diagnostics_module_name(module) << ":"
              << silex::verbose_level_name(level) << "] ";
    if (function != nullptr) {
        std::cerr << function << ": ";
    }
    std::cerr << (message == nullptr ? "" : message);
    if (detail != nullptr) {
        std::cerr << ": " << detail;
    }
    std::cerr << "\n";
}

silex::CertificationMode requested_certification(const Options& options) {
    if (options.mode == "proven") {
        return silex::CertificationMode::proven;
    }
    if (options.mode == "grh") {
        return silex::CertificationMode::grh;
    }
    return silex::CertificationMode::unknown;
}

bool set_polynomial(sflint::FmpqPoly& out,
                    const std::vector<slong>& coefficients) noexcept {
    if (coefficients.size() < 2) {
        return false;
    }
    sflint::fmpq_poly_zero(out);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        if (coefficients[i] != 0) {
            sflint::fmpq_poly_set_coeff_si(
                    out, static_cast<slong>(i), coefficients[i]);
        }
    }
    return true;
}

bool configure_compute_options(silex::ClassGroupComputeOptions& out,
                               sflint::Fmpz& factor_base_bound,
                               const silex::Order& order,
                               const Options& options) noexcept {
    out = silex::ClassGroupComputeOptions{};
    if (!silex::factor_base_class_group_bound(
                sflint::FmpzRef(factor_base_bound), order)) {
        return false;
    }
    if (sflint::fmpz_cmp_ui(sflint::FmpzConstRef(factor_base_bound), 2) < 0) {
        sflint::fmpz_set_ui(sflint::FmpzRef(factor_base_bound), 2);
    }
    if (options.factor_base_bound_override.has_value()) {
        sflint::fmpz_set_si(sflint::FmpzRef(factor_base_bound),
                            *options.factor_base_bound_override);
    }

    if (order.degree() >= 3) {
        out.max_candidates = 5000;
        out.max_relations = order.degree() >= 5 ? 1000 : 500;
    }

    if (options.max_candidates >= 0) {
        out.max_candidates = options.max_candidates;
    }
    if (options.max_relations >= 0) {
        out.max_relations = options.max_relations;
    }
    out.zeta_bf_max_cutoff = options.zeta_bf_max_cutoff;
    out.requested_certification = requested_certification(options);
    return true;
}

struct CanonicalPrimeCandidate {
    explicit CanonicalPrimeCandidate(const silex::Order& order) noexcept
        : prime(order), ideal(order) {
    }

    silex::PrimeIdeal prime;
    silex::Ideal ideal;
    sflint::FmpzMat hnf{0, 0};
};

struct SelectedPrimeDescriptor {
    sflint::Fmpz rational_prime;
    slong canonical_index = -1;
    slong ramification_index = 0;
    slong residue_degree = 0;
    std::vector<std::string> beta_power_basis;
    sflint::FmpzMat hnf{0, 0};
};

bool construct_witness_ideal(silex::Ideal& out,
                             const silex::Order& order,
                             sflint::FmpzConstRef rational_prime,
                             const std::vector<std::string>& beta) noexcept {
    const silex::NumberField* field = order.parent();
    if (field == nullptr ||
        beta.size() != static_cast<std::size_t>(order.degree())) {
        return false;
    }
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    for (std::size_t index = 0; index < beta.size(); ++index) {
        sflint::Fmpq coefficient;
        if (!parse_canonical_rational_text(beta[index], &coefficient)) {
            return false;
        }
        sflint::fmpq_poly_set_coeff_fmpq(
                polynomial, static_cast<slong>(index),
                sflint::FmpqConstRef(coefficient));
    }
    silex::Element element(*field);
    silex::OrderElement order_element(order);
    return element.is_defined() && order_element.is_defined() &&
           element.set_fmpq_poly(sflint::FmpqPolyConstRef(polynomial)) &&
           order_element.set_element(element) &&
           silex::detail::set_known_two_generator_ideal(
                   out, rational_prime, order_element);
}

void append_prime_descriptor(
        std::vector<SelectedPrimeDescriptor>& descriptors,
        const CanonicalPrimeCandidate& candidate,
        sflint::FmpzConstRef rational_prime,
        const SPrimeSelector& witness) noexcept {
    SelectedPrimeDescriptor descriptor;
    sflint::fmpz_set(sflint::FmpzRef(descriptor.rational_prime),
                     rational_prime);
    descriptor.canonical_index = witness.canonical_index;
    descriptor.ramification_index = candidate.prime.ramification_index();
    descriptor.residue_degree = candidate.prime.residue_degree();
    descriptor.beta_power_basis = witness.beta_power_basis;
    descriptor.hnf = sflint::FmpzMat(
            sflint::fmpz_mat_nrows(candidate.hnf),
            sflint::fmpz_mat_ncols(candidate.hnf));
    sflint::fmpz_mat_set(sflint::FmpzMatRef(descriptor.hnf),
                         sflint::FmpzMatConstRef(candidate.hnf));
    descriptors.push_back(std::move(descriptor));
}

bool append_selected_prime(
        std::vector<silex::PrimeIdeal>& primes,
        std::vector<SelectedPrimeDescriptor>& descriptors,
        const silex::Order& order,
        const CanonicalPrimeCandidate& candidate,
        sflint::FmpzConstRef rational_prime,
        const SPrimeSelector& witness) noexcept {
    for (const silex::PrimeIdeal& existing : primes) {
        if (existing.equal(candidate.prime)) {
            return false;
        }
    }

    primes.emplace_back(order);
    if (!primes.back().is_defined() ||
        !primes.back().set(candidate.prime)) {
        return false;
    }
    append_prime_descriptor(descriptors, candidate, rational_prime, witness);
    return true;
}

bool build_selected_primes(
        std::vector<silex::PrimeIdeal>& primes,
        std::vector<SelectedPrimeDescriptor>& descriptors,
        std::vector<SelectedPrimeDescriptor>& canonical_decompositions,
        const silex::Order& order,
        const std::vector<SPrimeSelector>& selectors,
        const silex::DiagnosticsContext* diagnostics) noexcept {
    primes.clear();
    descriptors.clear();
    canonical_decompositions.clear();
    std::vector<slong> selection_indices;
    std::vector<slong> processed_rational_primes;
    for (const SPrimeSelector& selector : selectors) {
        if (std::find(processed_rational_primes.begin(),
                      processed_rational_primes.end(),
                      selector.rational_prime) !=
            processed_rational_primes.end()) {
            continue;
        }
        processed_rational_primes.push_back(selector.rational_prime);
        sflint::Fmpz rational_prime;
        sflint::fmpz_set_si(sflint::FmpzRef(rational_prime),
                            selector.rational_prime);
        silex::PrimeIdealList decomposition;
        if (!silex::decompose_prime(
                    decomposition, order,
                    sflint::FmpzConstRef(rational_prime), order.degree(),
                    diagnostics) ||
            decomposition.size() <= 0) {
            return false;
        }

        std::vector<CanonicalPrimeCandidate> candidates;
        candidates.reserve(static_cast<std::size_t>(decomposition.size()));
        for (slong i = 0; i < decomposition.size(); ++i) {
            const silex::PrimeIdeal* source = decomposition.at(i);
            candidates.emplace_back(order);
            CanonicalPrimeCandidate& candidate = candidates.back();
            candidate.hnf = sflint::FmpzMat(order.degree(), order.degree());
            if (source == nullptr || !candidate.prime.is_defined() ||
                !candidate.prime.set(*source) || !candidate.ideal.is_defined() ||
                !candidate.prime.get_ideal(candidate.ideal) ||
                !candidate.ideal.get_hnf(
                        sflint::FmpzMatRef(candidate.hnf))) {
                return false;
            }
        }

        std::vector<const SPrimeSelector*> witnesses;
        for (const SPrimeSelector& grouped_witness : selectors) {
            if (grouped_witness.rational_prime == selector.rational_prime) {
                witnesses.push_back(&grouped_witness);
            }
        }
        std::sort(witnesses.begin(), witnesses.end(),
                  [](const SPrimeSelector* left,
                     const SPrimeSelector* right) {
                      return left->canonical_index < right->canonical_index;
                  });
        if (witnesses.size() != candidates.size()) {
            return false;
        }
        std::vector<bool> matched(candidates.size(), false);
        for (std::size_t witness_index = 0;
             witness_index < witnesses.size(); ++witness_index) {
            const SPrimeSelector& witness = *witnesses[witness_index];
            if (witness.canonical_index != static_cast<slong>(witness_index)) {
                return false;
            }
            silex::Ideal witness_ideal(order);
            if (!witness_ideal.is_defined() ||
                !construct_witness_ideal(
                        witness_ideal, order,
                        sflint::FmpzConstRef(rational_prime),
                        witness.beta_power_basis)) {
                return false;
            }
            std::size_t match_index = candidates.size();
            for (std::size_t candidate_index = 0;
                 candidate_index < candidates.size(); ++candidate_index) {
                if (candidates[candidate_index].ideal.equal(witness_ideal)) {
                    if (match_index != candidates.size()) {
                        return false;
                    }
                    match_index = candidate_index;
                }
            }
            if (match_index == candidates.size() || matched[match_index]) {
                return false;
            }
            matched[match_index] = true;
            const CanonicalPrimeCandidate& candidate = candidates[match_index];
            append_prime_descriptor(canonical_decompositions, candidate,
                                    sflint::FmpzConstRef(rational_prime),
                                    witness);
            if (witness.selection_index >= 0 &&
                !append_selected_prime(
                        primes, descriptors, order, candidate,
                        sflint::FmpzConstRef(rational_prime), witness)) {
                return false;
            }
            if (witness.selection_index >= 0) {
                selection_indices.push_back(witness.selection_index);
            }
        }
        if (std::find(matched.begin(), matched.end(), false) != matched.end()) {
            return false;
        }
    }
    for (std::size_t expected = 0; expected < selection_indices.size();
         ++expected) {
        const auto match = std::find(
                selection_indices.begin() + static_cast<std::ptrdiff_t>(expected),
                selection_indices.end(), static_cast<slong>(expected));
        if (match == selection_indices.end()) {
            return false;
        }
        const std::size_t actual = static_cast<std::size_t>(
                std::distance(selection_indices.begin(), match));
        if (actual != expected) {
            std::swap(selection_indices[expected], selection_indices[actual]);
            std::swap(primes[expected], primes[actual]);
            std::swap(descriptors[expected], descriptors[actual]);
        }
    }
    return true;
}

bool coordinates_equal(const silex::SUnitCoordinates& left,
                       const silex::SUnitCoordinates& right) noexcept {
    return left.defined && right.defined &&
           sflint::fmpz_equal(
                   sflint::FmpzConstRef(left.torsion_exponent),
                   sflint::FmpzConstRef(right.torsion_exponent)) &&
           sflint::fmpz_mat_equal(
                   sflint::FmpzMatConstRef(left.ordinary_free_exponents),
                   sflint::FmpzMatConstRef(right.ordinary_free_exponents)) &&
           sflint::fmpz_mat_equal(
                   sflint::FmpzMatConstRef(left.nonunit_exponents),
                   sflint::FmpzMatConstRef(right.nonunit_exponents));
}

struct SUnitMembershipAudit {
    bool mixed_round_trip_verified = false;
    bool outside_support_rejected = false;
    slong verified_round_trip_count = 0;
    silex::SUnitMembershipResult mixed_result;
    silex::SUnitMembershipResult outside_result;
};

bool audit_sunit_membership(
        SUnitMembershipAudit& audit,
        const silex::SUnitGroup& group,
        const silex::NumberField& field,
        const std::vector<SelectedPrimeDescriptor>& descriptors,
        slong precision) noexcept {
    audit = {};
    if (!group.is_defined() || precision <= 0) {
        return false;
    }

    silex::SUnitCoordinates input;
    input.ordinary_free_exponents =
            sflint::FmpzMat(1, group.ordinary_free_rank());
    input.nonunit_exponents = sflint::FmpzMat(1, group.nonunit_rank());
    sflint::fmpz_one(sflint::FmpzRef(input.torsion_exponent));
    for (slong i = 0; i < group.ordinary_free_rank(); ++i) {
        sflint::fmpz_set_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(input.ordinary_free_exponents), 0,
                        i),
                i + 1);
    }
    for (slong i = 0; i < group.nonunit_rank(); ++i) {
        sflint::fmpz_set_si(
                sflint::fmpz_mat_entry(
                        sflint::FmpzMatRef(input.nonunit_exponents), 0, i),
                i + 1);
    }
    input.defined = true;

    silex::Element value(field);
    silex::EmbeddingContext embeddings(field);
    silex::SUnitCoordinates recovered;
    const slong max_precision = precision <= WORD_MAX / 4
            ? std::max<slong>(256, 4 * precision)
            : precision;
    if (!value.is_defined() || !embeddings.is_defined() ||
        !group.image(value, input) ||
        !group.preimage(audit.mixed_result, recovered, value, embeddings, 16,
                        max_precision) ||
        audit.mixed_result.outcome !=
                silex::SUnitMembershipOutcome::verified ||
        !coordinates_equal(input, recovered)) {
        return false;
    }
    audit.mixed_round_trip_verified = true;
    audit.verified_round_trip_count = 1;

    slong outside_prime = 0;
    ulong candidate = 2;
    while (candidate <= static_cast<ulong>(WORD_MAX)) {
        bool selected = false;
        for (const SelectedPrimeDescriptor& descriptor : descriptors) {
            if (sflint::fmpz_equal_si(
                        sflint::FmpzConstRef(descriptor.rational_prime),
                        static_cast<slong>(candidate))) {
                selected = true;
                break;
            }
        }
        if (!selected) {
            outside_prime = static_cast<slong>(candidate);
            break;
        }
        const ulong next_candidate = n_nextprime(candidate, 1);
        if (next_candidate <= candidate) {
            return false;
        }
        candidate = next_candidate;
    }
    if (outside_prime == 0) {
        return false;
    }

    silex::Element outside(field);
    silex::SUnitCoordinates preserved;
    preserved.defined = true;
    if (!outside.set_si(outside_prime) ||
        !group.preimage(audit.outside_result, preserved, outside, embeddings,
                        16, max_precision) ||
        audit.outside_result.outcome !=
                silex::SUnitMembershipOutcome::not_sunit) {
        return false;
    }
    audit.outside_support_rejected = true;
    return true;
}

bool run_warmup(const Options& input_options) noexcept {
    if (input_options.warmup_coefficients.empty()) {
        return true;
    }

    Options options = input_options;
    options.coefficients = input_options.warmup_coefficients;
    options.warmup_coefficients.clear();
    options.expect_class_order.reset();
    options.expect_unit_rank.reset();
    options.expect_success.reset();
    options.s_prime_selectors.clear();
    options.compute_sunit = false;
    options.logging = false;
    options.trace = false;
    options.verbose = false;
    options.profiling = false;

    sflint::FmpqPoly polynomial;
    if (!set_polynomial(polynomial, options.coefficients)) {
        return false;
    }
    silex::NumberField field = silex::NumberField::by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    silex::Order equation_order = field.is_defined()
            ? silex::Order::equation_order(field)
            : silex::Order{};
    silex::Order maximal_order(field);
    if (!equation_order.is_defined() ||
        !maximal_order.maximal_order(equation_order)) {
        return false;
    }

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions compute_options;
    if (!configure_compute_options(compute_options, factor_base_bound,
                                   maximal_order, options)) {
        return false;
    }

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport audit;
    audit.reset();
    return compute_class_unit_from_options(
            units, class_group, maximal_order,
            sflint::FmpzConstRef(factor_base_bound), compute_options,
            options.precision, options, audit);
}

int print_error_json(const std::string& error, int exit_code) {
    std::cout << "{\n  \"success\": false,\n  \"error\": ";
    write_json_string(std::cout, error);
    std::cout << "\n}\n";
    return exit_code;
}

bool protocol_nonce_is_valid(std::string_view nonce) {
    return nonce.size() == target_nonce_bytes &&
           std::all_of(nonce.begin(), nonce.end(), [](unsigned char ch) {
               return (ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f');
           });
}

bool read_protocol_phase(
        std::string* phase = nullptr,
        std::size_t max_bytes = max_protocol_phase_bytes) {
    std::string value;
    value.reserve(std::min<std::size_t>(max_bytes, 4096));
    while (true) {
        const int next = std::cin.get();
        if (next == std::char_traits<char>::eof()) {
            if (value.empty()) {
                return false;
            }
            break;
        }
        if (next == '\n') {
            break;
        }
        if (value.size() >= max_bytes) {
            return false;
        }
        value.push_back(static_cast<char>(next));
    }
    if (phase != nullptr) {
        *phase = std::move(value);
    }
    return true;
}

void emit_protocol_marker(const char* marker,
                          std::string_view nonce = {}) {
    std::cout << marker;
    if (!nonce.empty()) {
        std::cout << ':' << nonce;
    }
    std::cout << '\n';
    std::cout.flush();
}

void print_profile_json(const ProfileCollector& collector) {
    const auto rows = collector.sorted_rows();
    std::cout << "[";
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        const double inclusive_ms =
                static_cast<double>(row.second.inclusive.count()) / 1.0e6;
        const double exclusive_ms =
                static_cast<double>(row.second.exclusive.count()) / 1.0e6;
        if (i != 0) {
            std::cout << ",";
        }
        std::cout << "\n    {\"label\": ";
        write_json_string(std::cout, row.first);
        std::cout << ", \"count\": " << row.second.count
                  << ", \"events\": " << row.second.events
                  << ", \"inclusive_ms\": " << inclusive_ms
                  << ", \"exclusive_ms\": " << exclusive_ms << "}";
    }
    if (!rows.empty()) {
        std::cout << "\n  ";
    }
    std::cout << "]";
}

void print_phase_timings_json(
        const silex::detail::ClassUnitTransactionReport& audit) {
    std::cout << "{";
    bool first = true;
    for (std::size_t i = 1;
         i < silex::detail::kClassUnitStageCount;
         ++i) {
        const auto stage =
                static_cast<silex::detail::ClassUnitStage>(i);
        if (stage == silex::detail::ClassUnitStage::count) {
            continue;
        }
        if (!first) {
            std::cout << ",";
        }
        first = false;
        std::cout << "\n    ";
        write_json_string(
                std::cout, silex::detail::class_unit_stage_name(stage));
        std::cout << ": ";
        if (audit.timing_recorded[i]) {
            std::cout << audit.timing_ms[i];
        } else {
            std::cout << "null";
        }
    }
    if (!first) {
        std::cout << "\n  ";
    }
    std::cout << "}";
}

void print_optional_class_order(const silex::ClassGroupContext& class_group) {
    sflint::Fmpz class_order;
    if (class_group.order(sflint::FmpzRef(class_order))) {
        write_json_string(std::cout,
                          fmpz_string(sflint::FmpzConstRef(class_order)));
    } else {
        std::cout << "null";
    }
}

void print_class_invariants(const silex::ClassGroupContext& class_group) {
    std::cout << "[";
    for (slong i = 0; i < class_group.invariant_count(); ++i) {
        sflint::Fmpz invariant;
        if (!class_group.invariant(sflint::FmpzRef(invariant), i)) {
            continue;
        }
        if (i != 0) {
            std::cout << ", ";
        }
        write_json_string(
                std::cout, fmpz_string(sflint::FmpzConstRef(invariant)));
    }
    std::cout << "]";
}

void print_s_class_invariants(const silex::SClassGroup& class_group) {
    std::cout << "[";
    for (slong i = 0; i < class_group.invariant_count(); ++i) {
        sflint::Fmpz invariant;
        if (!class_group.invariant(sflint::FmpzRef(invariant), i)) {
            continue;
        }
        if (i != 0) {
            std::cout << ", ";
        }
        write_json_string(
                std::cout, fmpz_string(sflint::FmpzConstRef(invariant)));
    }
    std::cout << "]";
}

void print_fmpz_matrix(sflint::FmpzMatConstRef matrix) {
    std::cout << "[";
    for (slong i = 0; i < sflint::fmpz_mat_nrows(matrix); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << "[";
        for (slong j = 0; j < sflint::fmpz_mat_ncols(matrix); ++j) {
            if (j != 0) {
                std::cout << ", ";
            }
            write_json_string(
                    std::cout,
                    fmpz_string(sflint::fmpz_mat_entry(matrix, i, j)));
        }
        std::cout << "]";
    }
    std::cout << "]";
}

void print_fmpq_matrix(sflint::FmpqMatConstRef matrix) {
    std::cout << "[";
    for (slong i = 0; i < ::fmpq_mat_nrows(matrix.raw()); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        std::cout << "[";
        for (slong j = 0; j < ::fmpq_mat_ncols(matrix.raw()); ++j) {
            if (j != 0) {
                std::cout << ", ";
            }
            write_json_string(
                    std::cout,
                    fmpq_string(sflint::fmpq_mat_entry(matrix, i, j)));
        }
        std::cout << "]";
    }
    std::cout << "]";
}

void print_selected_prime_descriptors(
        const std::vector<SelectedPrimeDescriptor>& descriptors) {
    std::cout << "[";
    for (std::size_t i = 0; i < descriptors.size(); ++i) {
        if (i != 0) {
            std::cout << ",";
        }
        const SelectedPrimeDescriptor& descriptor = descriptors[i];
        std::cout << "\n      {\"p\": ";
        write_json_string(
                std::cout,
                fmpz_string(sflint::FmpzConstRef(
                        descriptor.rational_prime)));
        std::cout << ", \"canonical_index\": "
                  << static_cast<long long>(descriptor.canonical_index)
                  << ", \"e\": "
                  << static_cast<long long>(descriptor.ramification_index)
                  << ", \"f\": "
                  << static_cast<long long>(descriptor.residue_degree)
                  << ", \"beta_power_basis\": [";
        for (std::size_t coefficient = 0;
             coefficient < descriptor.beta_power_basis.size(); ++coefficient) {
            if (coefficient != 0) {
                std::cout << ", ";
            }
            write_json_string(std::cout,
                              descriptor.beta_power_basis[coefficient]);
        }
        std::cout << "], \"hnf\": ";
        print_fmpz_matrix(sflint::FmpzMatConstRef(descriptor.hnf));
        std::cout << "}";
    }
    if (!descriptors.empty()) {
        std::cout << "\n    ";
    }
    std::cout << "]";
}

bool expectation_passed(const Options& options,
                        bool compute_success,
                        const silex::ClassGroupContext& class_group,
                        const silex::OrderUnitGroup& units) {
    bool passed = true;
    if (options.expect_success.has_value()) {
        passed = passed && compute_success == *options.expect_success;
    }
    if (options.expect_class_order.has_value()) {
        sflint::Fmpz class_order;
        passed = passed &&
                 class_group.order(sflint::FmpzRef(class_order)) &&
                 sflint::fmpz_equal_si(
                         sflint::FmpzConstRef(class_order),
                         *options.expect_class_order);
    }
    if (options.expect_unit_rank.has_value()) {
        passed = passed && units.is_set() &&
                 units.free_rank() == *options.expect_unit_rank;
    }
    return passed;
}

}  // namespace

int main(int argc, char** argv) {
    Options input_options;
    std::string parse_error;
    std::string target_nonce;
    if (!parse_options(argc, argv, input_options, parse_error)) {
        return print_error_json(parse_error, 2);
    }

    if (input_options.marked_protocol && !read_protocol_phase()) {
        return print_error_json(
                "marked protocol missing warmup phase input", 5);
    }

    if (!run_warmup(input_options)) {
        return print_error_json("warmup_computation_failed", 3);
    }

    if (input_options.marked_protocol) {
        emit_protocol_marker(ready_marker);
        if (!read_protocol_phase(&target_nonce, target_nonce_bytes) ||
            !protocol_nonce_is_valid(target_nonce)) {
            return print_error_json(
                    "marked protocol target nonce is invalid", 5);
        }
    }

    sflint::FmpqPoly polynomial;
    if (!set_polynomial(polynomial, input_options.coefficients)) {
        return print_error_json("failed to construct polynomial", 2);
    }

    silex::DiagnosticsContext diagnostics;
    ProfileCollector profile_collector;
    silex::diagnostics_context_init(diagnostics);
    const bool use_diagnostics = input_options.logging ||
                                 input_options.trace ||
                                 input_options.verbose ||
                                 input_options.profiling;
    if (input_options.logging || input_options.trace) {
        silex::diagnostics_set_logging(
                diagnostics,
                input_options.trace ? silex::LogLevel::trace
                                    : silex::LogLevel::detail,
                silex::diagnostics_all_modules, diagnostic_log_callback,
                nullptr);
    }
    if (input_options.verbose) {
        silex::diagnostics_set_verbose(
                diagnostics, silex::VerboseLevel::detail,
                silex::diagnostics_all_modules, diagnostic_verbose_callback,
                nullptr);
    }
    if (input_options.profiling) {
        profile_collector.configure(diagnostics);
    }

    const std::clock_t process_cpu_start = std::clock();
    const auto start = std::chrono::steady_clock::now();
    const auto field_setup_start = start;

    silex::NumberField field =
            silex::NumberField::by_polynomial(
                    sflint::FmpqPolyConstRef(polynomial));
    silex::Order equation_order =
            field.is_defined() ? silex::Order::equation_order(field)
                               : silex::Order{};
    const auto field_setup_end = std::chrono::steady_clock::now();
    const auto maximal_order_start = field_setup_end;
    silex::Order maximal_order(field);
    const bool maximal_defined =
            equation_order.is_defined() &&
            maximal_order.maximal_order(equation_order);
    const auto maximal_order_end = std::chrono::steady_clock::now();

    sflint::Fmpz factor_base_bound;
    silex::ClassGroupComputeOptions compute_options;
    const bool options_defined =
            maximal_defined &&
            configure_compute_options(compute_options, factor_base_bound,
                                      maximal_order, input_options);
    if (use_diagnostics) {
        compute_options.diagnostics = &diagnostics;
    }

    silex::ClassGroupContext class_group;
    silex::OrderUnitGroup units;
    silex::detail::ClassUnitTransactionReport transaction_report;
    transaction_report.reset();
    bool compute_success = false;
    const auto class_unit_start = std::chrono::steady_clock::now();
    if (options_defined) {
        compute_success = compute_class_unit_from_options(
                units, class_group, maximal_order,
                sflint::FmpzConstRef(factor_base_bound), compute_options,
                input_options.precision, input_options, transaction_report);
    } else {
        transaction_report.failure_stage =
                silex::detail::ClassUnitStage::factor_base_bound;
        transaction_report.failure_reason = "input_or_options_unavailable";
    }
    const auto class_unit_end = std::chrono::steady_clock::now();

    std::vector<silex::PrimeIdeal> selected_primes;
    std::vector<SelectedPrimeDescriptor> selected_prime_descriptors;
    std::vector<SelectedPrimeDescriptor> canonical_prime_decompositions;
    silex::SClassGroup s_class_group;
    silex::SUnitGroup s_unit_group;
    silex::SUnitComputeResult sunit_result;
    SUnitMembershipAudit membership_audit;
    bool prime_selection_success = false;
    bool sunit_success = false;
    bool membership_success = false;
    const auto prime_selection_start = class_unit_end;
    if (input_options.compute_sunit && compute_success) {
        prime_selection_success = build_selected_primes(
                selected_primes, selected_prime_descriptors,
                canonical_prime_decompositions, maximal_order,
                input_options.s_prime_selectors,
                use_diagnostics ? &diagnostics : nullptr);
    }
    const auto prime_selection_end = std::chrono::steady_clock::now();
    const auto sunit_start = prime_selection_end;
    if (input_options.compute_sunit && compute_success &&
        prime_selection_success) {
        silex::SUnitComputeOptions sunit_options;
        sunit_options.regulator_precision = input_options.precision;
        sunit_options.diagnostics =
                use_diagnostics ? &diagnostics : nullptr;
        sunit_success = silex::compute_sunit_groups(
                sunit_result, s_class_group, s_unit_group, class_group, units,
                silex::PrimeIdealSpan(selected_primes.data(),
                                      selected_primes.size()),
                sunit_options);
    }
    const auto sunit_end = std::chrono::steady_clock::now();
    const auto membership_start = sunit_end;
    if (input_options.compute_sunit && sunit_success) {
        membership_success = audit_sunit_membership(
                membership_audit, s_unit_group, field,
                selected_prime_descriptors, input_options.precision);
    }
    const auto membership_end = std::chrono::steady_clock::now();

    const auto end = membership_end;
    const std::clock_t process_cpu_end = std::clock();
    const bool overall_success =
            compute_success &&
            (!input_options.compute_sunit ||
             (prime_selection_success && sunit_success &&
              membership_success));
    const double elapsed_ms =
            std::chrono::duration<double, std::milli>(end - start).count();
    const bool have_process_cpu_time =
            process_cpu_start != static_cast<std::clock_t>(-1) &&
            process_cpu_end != static_cast<std::clock_t>(-1) &&
            process_cpu_end >= process_cpu_start;
    const double process_cpu_ms = have_process_cpu_time
            ? 1000.0 * static_cast<double>(
                               process_cpu_end - process_cpu_start) /
                      static_cast<double>(CLOCKS_PER_SEC)
            : 0.0;
    const double field_setup_ms =
            std::chrono::duration<double, std::milli>(
                    field_setup_end - field_setup_start)
                    .count();
    const double maximal_order_ms =
            std::chrono::duration<double, std::milli>(
                    maximal_order_end - maximal_order_start)
                    .count();
    const double class_unit_ms =
            std::chrono::duration<double, std::milli>(
                    class_unit_end - class_unit_start)
                    .count();
    const double prime_selection_ms =
            std::chrono::duration<double, std::milli>(
                    prime_selection_end - prime_selection_start)
                    .count();
    const double sunit_ms =
            std::chrono::duration<double, std::milli>(sunit_end - sunit_start)
                    .count();
    const double membership_ms =
            std::chrono::duration<double, std::milli>(
                    membership_end - membership_start)
                    .count();
    sflint::Fmpz maximal_order_discriminant;
    const bool have_maximal_order_discriminant =
            maximal_defined && maximal_order.discriminant(
                                       sflint::FmpzRef(
                                               maximal_order_discriminant));
    const slong field_degree = static_cast<slong>(
            input_options.coefficients.size() - 1);
    sflint::FmpqMat maximal_order_basis(field_degree, field_degree);
    const bool have_maximal_order_basis =
            maximal_defined && maximal_order.get_basis(
                                       sflint::FmpqMatRef(
                                               maximal_order_basis));
    sflint::Fmpz equation_order_index;
    const bool have_equation_order_index =
            maximal_defined && equation_order.index_in(
                                       sflint::FmpzRef(equation_order_index),
                                       maximal_order);
    silex::Signature field_signature;
    const bool have_signature = field.is_defined() &&
                                silex::signature(field_signature, field);
    const bool expectations_passed =
            expectation_passed(input_options, overall_success, class_group,
                               units);

    if (input_options.marked_protocol) {
        emit_protocol_marker(target_done_marker, target_nonce);
        if (!read_protocol_phase()) {
            return print_error_json(
                    "marked protocol missing final phase input", 5);
        }
    }

    std::cout << "{\n";
    std::cout << "  \"success\": " << json_bool(overall_success) << ",\n";
    std::cout << "  \"expectations_passed\": "
              << json_bool(expectations_passed) << ",\n";
    std::cout << "  \"engine_thread_count\": "
              << flint_get_num_threads() << ",\n";
    std::cout << "  \"elapsed_ms\": " << elapsed_ms << ",\n";
    std::cout << "  \"measurement_timing\": {\n";
    std::cout << "    \"algorithm_clock\": \"std_clock_process_cpu\",\n";
    std::cout << "    \"algorithm_scope\": \"post_warmup_target\",\n";
    std::cout << "    \"component_clock\": \"steady_clock\",\n";
    std::cout << "    \"warmup_excluded\": true,\n";
    std::cout << "    \"target_cpu_ms\": ";
    if (have_process_cpu_time) {
        std::cout << process_cpu_ms;
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "    \"target_wall_ms\": " << elapsed_ms << "\n";
    std::cout << "  },\n";
    std::cout << "  \"component_timing_ms\": {\n";
    std::cout << "    \"field_setup\": " << field_setup_ms << ",\n";
    std::cout << "    \"maximal_order\": " << maximal_order_ms << ",\n";
    std::cout << "    \"class_unit\": " << class_unit_ms << ",\n";
    if (input_options.compute_sunit) {
        std::cout << "    \"s_prime_selection\": "
                  << prime_selection_ms << ",\n";
        std::cout << "    \"sunit\": " << sunit_ms << ",\n";
        std::cout << "    \"sunit_membership\": " << membership_ms
                  << ",\n";
    }
    std::cout << "    \"total\": " << elapsed_ms << "\n";
    std::cout << "  },\n";
    std::cout << "  \"warmup\": {\n";
    std::cout << "    \"used\": "
              << json_bool(!input_options.warmup_coefficients.empty())
              << ",\n";
    std::cout << "    \"degree\": ";
    if (input_options.warmup_coefficients.empty()) {
        std::cout << "null\n";
    } else {
        std::cout << static_cast<long long>(
                             input_options.warmup_coefficients.size() - 1)
                  << "\n";
    }
    std::cout << "  },\n";
    std::cout << "  \"mode\": ";
    write_json_string(std::cout, input_options.mode);
    std::cout << ",\n";
    std::cout << "  \"engine\": ";
    write_json_string(std::cout, "silex");
    std::cout << ",\n";
    std::cout << "  \"failure_stage\": ";
    if (transaction_report.failure_stage ==
        silex::detail::ClassUnitStage::none) {
        std::cout << "null";
    } else {
        write_json_string(
                std::cout,
                silex::detail::class_unit_stage_name(
                        transaction_report.failure_stage));
    }
    std::cout << ",\n";
    std::cout << "  \"failure_reason\": ";
    if (transaction_report.failure_reason == nullptr) {
        std::cout << "null";
    } else {
        write_json_string(std::cout, transaction_report.failure_reason);
    }
    std::cout << ",\n";
    std::cout << "  \"certification_status\": ";
    write_json_string(
            std::cout,
            certification_name(transaction_report.class_group_certification));
    std::cout << ",\n";
    std::cout << "  \"class_group_proof_status\": ";
    write_json_string(
            std::cout,
            certification_name(transaction_report.class_group_certification));
    std::cout << ",\n";
    std::cout << "  \"unit_group_proof_status\": ";
    write_json_string(
            std::cout,
            certification_name(transaction_report.unit_group_certification));
    std::cout << ",\n";
    std::cout << "  \"regulator_proof_status\": ";
    write_json_string(
            std::cout,
            proof_state_name(transaction_report.regulator_proof_status));
    std::cout << ",\n";
    std::cout << "  \"final_result_published\": "
              << json_bool(transaction_report.final_result_published) << ",\n";
    std::cout << "  \"phase_timing_ms\": ";
    print_phase_timings_json(transaction_report);
    std::cout << ",\n";
    std::cout << "  \"degree\": "
              << static_cast<long long>(input_options.coefficients.size() - 1)
              << ",\n";
    std::cout << "  \"signature\": ";
    if (have_signature) {
        std::cout << "[" << static_cast<long long>(field_signature.r1())
                  << ", " << static_cast<long long>(field_signature.r2())
                  << "]";
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "  \"maximal_order_discriminant\": ";
    if (have_maximal_order_discriminant) {
        write_json_string(
                std::cout,
                fmpz_string(sflint::FmpzConstRef(
                        maximal_order_discriminant)));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "  \"maximal_order_basis_power\": ";
    if (have_maximal_order_basis) {
        print_fmpq_matrix(
                sflint::FmpqMatConstRef(maximal_order_basis));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "  \"equation_order_index\": ";
    if (have_equation_order_index) {
        write_json_string(
                std::cout,
                fmpz_string(sflint::FmpzConstRef(equation_order_index)));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "  \"field_defined\": "
              << json_bool(field.is_defined()) << ",\n";
    std::cout << "  \"equation_order_defined\": "
              << json_bool(equation_order.is_defined()) << ",\n";
    std::cout << "  \"maximal_order_defined\": "
              << json_bool(maximal_defined) << ",\n";
    std::cout << "  \"precision\": "
              << static_cast<long long>(input_options.precision) << ",\n";
    std::cout << "  \"factor_base_bound\": ";
    if (options_defined) {
        write_json_string(std::cout,
                          fmpz_string(sflint::FmpzConstRef(
                                  factor_base_bound)));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";

    std::cout << "  \"options\": {\n";
    std::cout << "    \"max_candidates\": "
              << static_cast<long long>(compute_options.max_candidates)
              << ",\n";
    std::cout << "    \"max_relations\": "
              << static_cast<long long>(compute_options.max_relations)
              << ",\n";
    std::cout << "    \"zeta_bf_max_cutoff\": "
              << static_cast<unsigned long long>(
                         compute_options.zeta_bf_max_cutoff) << ",\n";
    std::cout << "    \"requested_certification\": ";
    write_json_string(
            std::cout,
            certification_name(compute_options.requested_certification));
    std::cout << "\n";
    std::cout << "  },\n";

    std::cout << "  \"class_group\": {\n";
    std::cout << "    \"has_factor_base\": "
              << json_bool(class_group.has_factor_base()) << ",\n";
    std::cout << "    \"generator_count\": "
              << static_cast<long long>(class_group.generator_count())
              << ",\n";
    std::cout << "    \"has_presentation\": "
              << json_bool(class_group.has_presentation()) << ",\n";
    std::cout << "    \"certification\": ";
    write_json_string(std::cout,
                      certification_name(
                              class_group.certification_status()));
    std::cout << ",\n";
    std::cout << "    \"factor_base_generation_status\": ";
    write_json_string(std::cout,
                      proof_state_name(
                              class_group.factor_base_generation_status()));
    std::cout << ",\n";
    std::cout << "    \"relation_saturation_status\": ";
    write_json_string(std::cout,
                      proof_state_name(
                              class_group.relation_saturation_status()));
    std::cout << ",\n";
    std::cout << "    \"analytic_class_regulator_status\": ";
    write_json_string(
            std::cout,
            proof_state_name(
                    class_group.analytic_class_regulator_status()));
    std::cout << ",\n";
    std::cout << "    \"zeta_bf_proof_status\": ";
    write_json_string(std::cout,
                      proof_state_name(class_group.zeta_bf_proof_status()));
    std::cout << ",\n";
    std::cout << "    \"unit_proof_status\": ";
    write_json_string(std::cout,
                      proof_state_name(class_group.unit_proof_status()));
    std::cout << ",\n";
    std::cout << "    \"regulator_proof_status\": ";
    write_json_string(std::cout,
                      proof_state_name(class_group.regulator_proof_status()));
    std::cout << ",\n";
    std::cout << "    \"relation_count\": "
              << static_cast<long long>(class_group.relation_count())
              << ",\n";
    std::cout << "    \"relation_rank\": "
              << static_cast<long long>(class_group.relation_rank()) << ",\n";
    std::cout << "    \"skipped_dependent_relation_count\": "
              << static_cast<long long>(
                         class_group.skipped_dependent_relation_count())
              << ",\n";
    std::cout << "    \"relation_kernel_unit_count\": "
              << static_cast<long long>(
                         class_group.relation_kernel_unit_count())
              << ",\n";
    std::cout << "    \"order\": ";
    print_optional_class_order(class_group);
    std::cout << ",\n";
    std::cout << "    \"invariants\": ";
    print_class_invariants(class_group);
    std::cout << ",\n";
    std::cout << "    \"relation_sources\": {\n";
    std::cout << "      \"search\": "
              << static_cast<long long>(
                         class_group.relation_source_count(
                                 silex::ClassGroupRelationSource::Search))
              << ",\n";
    std::cout << "      \"random_product\": "
              << static_cast<long long>(
                         class_group.relation_source_count(
                                 silex::ClassGroupRelationSource
                                         ::RandomProduct))
              << ",\n";
    std::cout << "      \"supplied\": "
              << static_cast<long long>(
                         class_group.relation_source_count(
                                 silex::ClassGroupRelationSource::Supplied))
              << ",\n";
    std::cout << "      \"saturation\": "
              << static_cast<long long>(
                         class_group.relation_source_count(
                                 silex::ClassGroupRelationSource::Saturation))
              << ",\n";
    std::cout << "      \"large_prime_match\": "
              << static_cast<long long>(
                         class_group.relation_source_count(
                                 silex::ClassGroupRelationSource
                                         ::LargePrimeMatch))
              << ",\n";
    std::cout << "      \"class_generator\": "
              << static_cast<long long>(
                         class_group.relation_source_count(
                                 silex::ClassGroupRelationSource
                                         ::ClassGenerator))
              << "\n";
    std::cout << "    }\n";
    std::cout << "  },\n";

    std::cout << "  \"unit_group\": {\n";
    std::cout << "    \"is_set\": " << json_bool(units.is_set()) << ",\n";
    std::cout << "    \"free_rank\": "
              << static_cast<long long>(units.free_rank()) << ",\n";
    std::cout << "    \"certification\": ";
    write_json_string(std::cout,
                      certification_name(units.certification_status()));
    std::cout << ",\n";
    std::cout << "    \"unit_proof_record_count\": "
              << static_cast<long long>(units.unit_proof_record_count())
              << "\n";
    std::cout << "  },\n";

    std::cout << "  \"sunit\": {\n";
    std::cout << "    \"requested\": "
              << json_bool(input_options.compute_sunit) << ",\n";
    std::cout << "    \"success\": " << json_bool(sunit_success) << ",\n";
    std::cout << "    \"failure_stage\": ";
    if (!input_options.compute_sunit) {
        std::cout << "null";
    } else if (!compute_success) {
        write_json_string(std::cout, "source_class_unit");
    } else if (!prime_selection_success) {
        write_json_string(std::cout, "s_prime_selection");
    } else if (!sunit_success) {
        write_json_string(std::cout,
                          silex::sunit_compute_stage_name(
                                  sunit_result.stage));
    } else if (!membership_success) {
        write_json_string(std::cout, "membership");
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "    \"failure_selected_index\": ";
    if (!sunit_success && sunit_result.selected_index >= 0) {
        std::cout << static_cast<long long>(sunit_result.selected_index);
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "    \"final_result_published\": "
              << json_bool(sunit_success) << ",\n";
    std::cout << "    \"selected_primes\": ";
    print_selected_prime_descriptors(selected_prime_descriptors);
    std::cout << ",\n";
    std::cout << "    \"canonical_prime_decompositions\": ";
    print_selected_prime_descriptors(canonical_prime_decompositions);
    std::cout << ",\n";
    std::cout << "    \"s_class_group\": {\n";
    std::cout << "      \"order\": ";
    sflint::Fmpz s_class_order;
    if (s_class_group.order(sflint::FmpzRef(s_class_order))) {
        write_json_string(
                std::cout,
                fmpz_string(sflint::FmpzConstRef(s_class_order)));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "      \"invariants\": ";
    print_s_class_invariants(s_class_group);
    std::cout << ",\n";
    std::cout << "      \"certification_status\": ";
    write_json_string(
            std::cout,
            certification_name(s_class_group.certification_status()));
    std::cout << ",\n";
    std::cout << "      \"proof_status\": ";
    write_json_string(std::cout,
                      proof_state_name(s_class_group.proof_status()));
    std::cout << "\n    },\n";
    std::cout << "    \"s_unit_group\": {\n";
    sflint::Fmpz torsion_order;
    std::cout << "      \"torsion_order\": ";
    if (s_unit_group.torsion_order(sflint::FmpzRef(torsion_order))) {
        write_json_string(
                std::cout,
                fmpz_string(sflint::FmpzConstRef(torsion_order)));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "      \"ordinary_free_rank\": "
              << static_cast<long long>(s_unit_group.ordinary_free_rank())
              << ",\n";
    std::cout << "      \"nonunit_rank\": "
              << static_cast<long long>(s_unit_group.nonunit_rank())
              << ",\n";
    std::cout << "      \"free_rank\": "
              << static_cast<long long>(s_unit_group.free_rank()) << ",\n";
    sflint::FmpzMat valuation_matrix(
            sunit_success ? s_unit_group.nonunit_rank() : 0,
            sunit_success ? s_unit_group.selected_prime_count() : 0);
    const bool have_valuation_matrix =
            sunit_success && s_unit_group.nonunit_valuation_matrix(
                                     sflint::FmpzMatRef(valuation_matrix));
    std::cout << "      \"valuation_matrix\": ";
    if (have_valuation_matrix) {
        print_fmpz_matrix(sflint::FmpzMatConstRef(valuation_matrix));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    sflint::Fmpz valuation_index;
    const bool have_valuation_index =
            have_valuation_matrix &&
            sflint::fmpz_mat_nrows(valuation_matrix) ==
                    sflint::fmpz_mat_ncols(valuation_matrix);
    if (have_valuation_index) {
        ::fmpz_mat_det(valuation_index.raw(), valuation_matrix.raw());
        ::fmpz_abs(valuation_index.raw(), valuation_index.raw());
    }
    std::cout << "      \"valuation_lattice_index\": ";
    if (have_valuation_index) {
        write_json_string(
                std::cout,
                fmpz_string(sflint::FmpzConstRef(valuation_index)));
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    sflint::Arb s_regulator;
    const bool have_s_regulator =
            s_unit_group.regulator(sflint::ArbRef(s_regulator));
    std::cout << "      \"regulator_midpoint\": ";
    if (have_s_regulator) {
        std::cout << std::setprecision(17)
                  << ::arf_get_d(arb_midref(s_regulator.raw()),
                                 ARF_RND_NEAR);
    } else {
        std::cout << "null";
    }
    std::cout << ",\n";
    std::cout << "      \"certification_status\": ";
    write_json_string(
            std::cout,
            certification_name(s_unit_group.certification_status()));
    std::cout << ",\n";
    std::cout << "      \"proof_status\": ";
    write_json_string(std::cout,
                      proof_state_name(s_unit_group.proof_status()));
    std::cout << ",\n";
    std::cout << "      \"regulator_proof_status\": ";
    write_json_string(
            std::cout,
            proof_state_name(s_unit_group.regulator_proof_status()));
    std::cout << "\n    },\n";
    std::cout << "    \"membership\": {\n";
    std::cout << "      \"status\": ";
    write_json_string(
            std::cout,
            membership_success ? "verified" : "unknown");
    std::cout << ",\n";
    std::cout << "      \"mixed_round_trip_verified\": "
              << json_bool(
                         membership_audit.mixed_round_trip_verified)
              << ",\n";
    std::cout << "      \"verified_round_trip_count\": "
              << static_cast<long long>(
                         membership_audit.verified_round_trip_count)
              << ",\n";
    std::cout << "      \"mixed_outcome\": ";
    write_json_string(
            std::cout,
            membership_outcome_name(
                    membership_audit.mixed_result.outcome));
    std::cout << ",\n";
    std::cout << "      \"outside_support_rejected\": "
              << json_bool(
                         membership_audit.outside_support_rejected)
              << ",\n";
    std::cout << "      \"outside_outcome\": ";
    write_json_string(
            std::cout,
            membership_outcome_name(
                    membership_audit.outside_result.outcome));
    std::cout << "\n    },\n";
    std::cout << "    \"timing_ms\": {\n";
    std::cout << "      \"prime_selection\": " << prime_selection_ms
              << ",\n";
    std::cout << "      \"construction\": " << sunit_ms << ",\n";
    std::cout << "      \"membership\": " << membership_ms << ",\n";
    std::cout << "      \"total\": "
              << prime_selection_ms + sunit_ms + membership_ms << "\n";
    std::cout << "    }\n";
    std::cout << "  },\n";

    std::cout << "  \"expectations\": {\n";
    std::cout << "    \"passed\": " << json_bool(expectations_passed);
    if (input_options.expect_success.has_value()) {
        std::cout << ",\n    \"success\": "
                  << json_bool(*input_options.expect_success);
    }
    if (input_options.expect_class_order.has_value()) {
        std::cout << ",\n    \"class_order\": "
                  << static_cast<long long>(
                             *input_options.expect_class_order);
    }
    if (input_options.expect_unit_rank.has_value()) {
        std::cout << ",\n    \"unit_rank\": "
                  << static_cast<long long>(
                             *input_options.expect_unit_rank);
    }
    std::cout << "\n  },\n";

    std::cout << "  \"profile\": ";
    print_profile_json(profile_collector);
    std::cout << "\n";
    std::cout << "}\n";

    return expectations_passed &&
                    (!input_options.compute_sunit || overall_success)
            ? 0
            : 1;
}
