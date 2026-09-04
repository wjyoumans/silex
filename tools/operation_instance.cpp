#include <silex/element.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_poly.hpp>
#include <silex/ideal.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>
#include <silex/order_element.hpp>

#include <flint/flint.h>
#include <flint/fmpz_poly.h>

#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace sflint = silex::flint;

enum class Operation {
    unknown,
    maximal_order,
    ideal_multiply,
    element_square_root,
};

struct Options {
    std::vector<std::string> coefficients;
    std::vector<std::string> warmup_coefficients;
    Operation operation = Operation::unknown;
    bool marked_protocol = false;
};

constexpr const char* ready_marker = "__SILEX_BENCH_SILEX_READY__";
constexpr const char* target_done_marker =
        "__SILEX_BENCH_SILEX_TARGET_DONE__";
constexpr std::size_t max_protocol_phase_bytes = 1U << 20U;
constexpr std::size_t target_nonce_bytes = 32;

struct OperationResult {
    bool success = false;
    std::string error;
    std::optional<double> target_cpu_ms;
    std::optional<double> target_wall_ms;
    std::optional<std::string> maximal_order_discriminant;
    std::optional<std::string> ideal_norm;
    std::optional<bool> root_found;
    std::optional<bool> root_verified;
};

struct TargetTimer {
    std::clock_t cpu_start = std::clock();
    std::chrono::steady_clock::time_point wall_start =
            std::chrono::steady_clock::now();

    void finish(OperationResult& result) const noexcept {
        const auto wall_end = std::chrono::steady_clock::now();
        const std::clock_t cpu_end = std::clock();
        result.target_wall_ms =
                std::chrono::duration<double, std::milli>(
                        wall_end - wall_start)
                        .count();
        if (cpu_start != static_cast<std::clock_t>(-1) &&
            cpu_end != static_cast<std::clock_t>(-1) &&
            cpu_end >= cpu_start) {
            result.target_cpu_ms =
                    1000.0 * static_cast<double>(cpu_end - cpu_start) /
                    static_cast<double>(CLOCKS_PER_SEC);
        }
    }
};

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
                const char* hex = "0123456789abcdef";
                out << "\\u00"
                    << hex[(static_cast<unsigned char>(ch) >> 4) & 0xf]
                    << hex[static_cast<unsigned char>(ch) & 0xf];
            } else {
                out << ch;
            }
            break;
        }
    }
    out << '"';
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

const char* operation_name(Operation operation) noexcept {
    switch (operation) {
    case Operation::maximal_order:
        return "maximal_order";
    case Operation::ideal_multiply:
        return "ideal_multiply";
    case Operation::element_square_root:
        return "element_square_root";
    case Operation::unknown:
        break;
    }
    return "unknown";
}

const char* timing_scope(Operation operation) noexcept {
    switch (operation) {
    case Operation::maximal_order:
        return "field_construction+equation_order+maximal_order";
    case Operation::ideal_multiply:
        return "ideal_multiply_only";
    case Operation::element_square_root:
        return "element_is_square_with_root_only";
    case Operation::unknown:
        break;
    }
    return "unknown";
}

bool parse_operation(std::string_view value, Operation& operation) noexcept {
    if (value == "maximal_order") {
        operation = Operation::maximal_order;
        return true;
    }
    if (value == "ideal_multiply") {
        operation = Operation::ideal_multiply;
        return true;
    }
    if (value == "element_square_root") {
        operation = Operation::element_square_root;
        return true;
    }
    return false;
}

bool parse_coefficients(const std::string& value,
                        std::vector<std::string>& out) {
    out.clear();
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t comma = value.find(',', start);
        const std::size_t end =
                comma == std::string::npos ? value.size() : comma;
        const std::string token = value.substr(start, end - start);
        sflint::Fmpz coefficient;
        if (token.empty() ||
            !sflint::fmpz_set_str(sflint::FmpzRef(coefficient),
                                  token.c_str())) {
            return false;
        }
        std::string normalized =
                fmpz_string(sflint::FmpzConstRef(coefficient));
        if (normalized.empty()) {
            return false;
        }
        out.push_back(std::move(normalized));
        if (comma == std::string::npos) {
            break;
        }
        start = comma + 1;
    }
    return out.size() >= 2 && out.back() == "1";
}

bool take_value(int& index,
                int argc,
                char** argv,
                std::string_view argument,
                std::string_view option,
                std::string& value,
                std::string& error) {
    if (argument == option) {
        if (index + 1 >= argc) {
            error = "missing value for ";
            error += option;
            return false;
        }
        value = argv[++index];
        return true;
    }
    std::string prefix(option);
    prefix += "=";
    if (argument.starts_with(prefix)) {
        value = std::string(argument.substr(prefix.size()));
        if (value.empty()) {
            error = "missing value for ";
            error += option;
            return false;
        }
        return true;
    }
    return false;
}

bool parse_options(int argc,
                   char** argv,
                   Options& options,
                   std::string& error) {
    bool have_coefficients = false;
    bool have_operation = false;
    for (int i = 1; i < argc; ++i) {
        const std::string_view argument(argv[i]);
        std::string value;
        if (argument == "--coeffs" || argument.starts_with("--coeffs=")) {
            if (!take_value(i, argc, argv, argument, "--coeffs", value,
                            error) ||
                !parse_coefficients(value, options.coefficients)) {
                if (error.empty()) {
                    error = "invalid --coeffs value; expected low-to-high "
                            "monic integer coefficients";
                }
                return false;
            }
            have_coefficients = true;
            continue;
        }
        if (argument == "--warmup-coeffs" ||
            argument.starts_with("--warmup-coeffs=")) {
            if (!take_value(i, argc, argv, argument, "--warmup-coeffs",
                            value, error) ||
                !parse_coefficients(value, options.warmup_coefficients)) {
                if (error.empty()) {
                    error = "invalid --warmup-coeffs value; expected "
                            "low-to-high monic integer coefficients";
                }
                return false;
            }
            continue;
        }
        if (argument == "--operation" ||
            argument.starts_with("--operation=")) {
            if (!take_value(i, argc, argv, argument, "--operation", value,
                            error) ||
                !parse_operation(value, options.operation)) {
                if (error.empty()) {
                    error = "invalid --operation value";
                }
                return false;
            }
            have_operation = true;
            continue;
        }
        if (argument == "--marked-protocol") {
            options.marked_protocol = true;
            continue;
        }
        error = "unknown argument: ";
        error += argument;
        return false;
    }

    if (!have_coefficients) {
        error = "missing --coeffs";
        return false;
    }
    if (!have_operation) {
        error = "missing --operation";
        return false;
    }
    if (!options.warmup_coefficients.empty()) {
        if (options.warmup_coefficients.size() !=
            options.coefficients.size()) {
            error = "--warmup-coeffs must have the same degree as --coeffs";
            return false;
        }
        if (options.warmup_coefficients == options.coefficients) {
            error = "--warmup-coeffs must define a distinct field";
            return false;
        }
    }
    return true;
}

bool set_polynomial(sflint::FmpzPoly& polynomial,
                    const std::vector<std::string>& coefficients) {
    fmpz_poly_zero(polynomial.raw());
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        sflint::Fmpz coefficient;
        if (!sflint::fmpz_set_str(sflint::FmpzRef(coefficient),
                                  coefficients[i].c_str())) {
            return false;
        }
        // This executable is direct low-level bridge code: the owning
        // polynomial and coefficient remain Silex RAII values while FLINT
        // receives borrowed handles for exact coefficient ingestion.
        fmpz_poly_set_coeff_fmpz(
                polynomial.raw(), static_cast<slong>(i), coefficient.raw());
    }
    return fmpz_poly_degree(polynomial.raw()) ==
           static_cast<slong>(coefficients.size() - 1);
}

bool setup_field(const std::vector<std::string>& coefficients,
                 silex::NumberField& field,
                 std::string& error) {
    sflint::FmpzPoly polynomial;
    if (!set_polynomial(polynomial, coefficients)) {
        error = "failed to construct defining polynomial";
        return false;
    }
    silex::NumberField candidate = silex::NumberField::by_polynomial(
            sflint::FmpzPolyConstRef(polynomial));
    if (!candidate.is_defined()) {
        error = "failed to construct number field";
        return false;
    }
    field = std::move(candidate);
    return true;
}

bool setup_maximal_order(const std::vector<std::string>& coefficients,
                         silex::NumberField& field,
                         silex::Order& maximal_order,
                         std::string& error) {
    if (!setup_field(coefficients, field, error)) {
        return false;
    }
    silex::Order equation_order = silex::Order::equation_order(field);
    if (!equation_order.is_defined()) {
        error = "failed to construct equation order";
        return false;
    }
    silex::Order candidate(field);
    if (!candidate.is_defined() ||
        !candidate.maximal_order(equation_order)) {
        error = "failed to construct maximal order";
        return false;
    }
    maximal_order = std::move(candidate);
    return true;
}

OperationResult run_maximal_order(
        const std::vector<std::string>& coefficients) {
    OperationResult result;
    silex::NumberField field;
    silex::Order maximal_order;

    const TargetTimer timer;
    const bool constructed =
            setup_maximal_order(coefficients, field, maximal_order,
                                result.error);
    timer.finish(result);
    if (!constructed) {
        return result;
    }

    sflint::Fmpz discriminant;
    if (!maximal_order.discriminant(sflint::FmpzRef(discriminant))) {
        result.error = "failed to read maximal-order discriminant";
        return result;
    }
    result.maximal_order_discriminant =
            fmpz_string(sflint::FmpzConstRef(discriminant));
    if (result.maximal_order_discriminant->empty()) {
        result.error = "failed to serialize maximal-order discriminant";
        result.maximal_order_discriminant.reset();
        return result;
    }
    result.success = true;
    return result;
}

OperationResult run_ideal_multiply(
        const std::vector<std::string>& coefficients) {
    OperationResult result;
    silex::NumberField field;
    silex::Order maximal_order;
    if (!setup_maximal_order(coefficients, field, maximal_order,
                             result.error)) {
        return result;
    }

    silex::OrderElement two(maximal_order);
    silex::OrderElement three(maximal_order);
    silex::Ideal left(maximal_order);
    silex::Ideal right(maximal_order);
    silex::Ideal product(maximal_order);
    if (!two.is_defined() || !three.is_defined() || !left.is_defined() ||
        !right.is_defined() || !product.is_defined() || !two.set_si(2) ||
        !three.set_si(3) || !left.set_principal(two) ||
        !right.set_principal(three)) {
        result.error = "failed to construct principal ideal inputs";
        return result;
    }

    const TargetTimer timer;
    const bool multiplied = product.multiply(left, right);
    timer.finish(result);
    if (!multiplied) {
        result.error = "ideal multiplication failed";
        return result;
    }

    sflint::Fmpz norm;
    if (!product.norm(sflint::FmpzRef(norm))) {
        result.error = "failed to compute product ideal norm";
        return result;
    }
    result.ideal_norm = fmpz_string(sflint::FmpzConstRef(norm));
    if (result.ideal_norm->empty()) {
        result.error = "failed to serialize product ideal norm";
        result.ideal_norm.reset();
        return result;
    }
    result.success = true;
    return result;
}

OperationResult run_element_square_root(
        const std::vector<std::string>& coefficients) {
    OperationResult result;
    silex::NumberField field;
    if (!setup_field(coefficients, field, result.error)) {
        return result;
    }

    silex::Element generator(field);
    silex::Element one_plus_generator(field);
    silex::Element square(field);
    silex::Element root(field);
    silex::Element check(field);
    if (!generator.is_defined() || !one_plus_generator.is_defined() ||
        !square.is_defined() || !root.is_defined() || !check.is_defined() ||
        !generator.gen() || !one_plus_generator.add_si(generator, 1) ||
        !square.multiply(one_plus_generator, one_plus_generator)) {
        result.error = "failed to construct square input";
        return result;
    }

    bool root_found = false;
    const TargetTimer timer;
    const bool completed = square.is_square(root_found, root);
    timer.finish(result);
    result.root_found = root_found;
    result.root_verified =
            completed && root_found && check.multiply(root, root) &&
            check.equal(square);
    if (!completed) {
        result.error = "square-root operation could not decide the input";
        return result;
    }
    if (!root_found) {
        result.error = "known-square input was reported nonsquare";
        return result;
    }
    if (!*result.root_verified) {
        result.error = "returned square root failed exact verification";
        return result;
    }
    result.success = true;
    return result;
}

OperationResult run_operation(
        Operation operation,
        const std::vector<std::string>& coefficients) {
    switch (operation) {
    case Operation::maximal_order:
        return run_maximal_order(coefficients);
    case Operation::ideal_multiply:
        return run_ideal_multiply(coefficients);
    case Operation::element_square_root:
        return run_element_square_root(coefficients);
    case Operation::unknown:
        break;
    }
    OperationResult result;
    result.error = "unknown operation";
    return result;
}

void write_optional_string(
        std::ostream& out,
        const std::optional<std::string>& value) {
    if (value.has_value()) {
        write_json_string(out, *value);
    } else {
        out << "null";
    }
}

void write_optional_bool(std::ostream& out,
                         const std::optional<bool>& value) {
    if (value.has_value()) {
        out << (*value ? "true" : "false");
    } else {
        out << "null";
    }
}

void write_optional_double(std::ostream& out,
                           const std::optional<double>& value) {
    if (value.has_value()) {
        out << *value;
    } else {
        out << "null";
    }
}

void write_result(const Options& options, const OperationResult& result) {
    std::cout << std::setprecision(17);
    std::cout << "{\n";
    std::cout << "  \"engine\": \"silex\",\n";
    std::cout << "  \"operation\": ";
    write_json_string(std::cout, operation_name(options.operation));
    std::cout << ",\n";
    std::cout << "  \"success\": "
              << (result.success ? "true" : "false") << ",\n";
    std::cout << "  \"engine_thread_count\": "
              << flint_get_num_threads() << ",\n";
    std::cout << "  \"error\": ";
    if (result.error.empty()) {
        std::cout << "null";
    } else {
        write_json_string(std::cout, result.error);
    }
    std::cout << ",\n";
    std::cout << "  \"target_cpu_ms\": ";
    write_optional_double(std::cout, result.target_cpu_ms);
    std::cout << ",\n";
    std::cout << "  \"target_wall_ms\": ";
    write_optional_double(std::cout, result.target_wall_ms);
    std::cout << ",\n";
    std::cout << "  \"maximal_order_discriminant\": ";
    write_optional_string(std::cout, result.maximal_order_discriminant);
    std::cout << ",\n";
    std::cout << "  \"ideal_norm\": ";
    write_optional_string(std::cout, result.ideal_norm);
    std::cout << ",\n";
    std::cout << "  \"root_found\": ";
    write_optional_bool(std::cout, result.root_found);
    std::cout << ",\n";
    std::cout << "  \"root_verified\": ";
    write_optional_bool(std::cout, result.root_verified);
    std::cout << ",\n";
    std::cout << "  \"source\": \"silex_public_api\",\n";
    std::cout << "  \"timing_scope\": ";
    write_json_string(std::cout, timing_scope(options.operation));
    std::cout << ",\n";
    std::cout << "  \"timing_clock\": {\n";
    std::cout << "    \"cpu\": \"std_clock_process_cpu\",\n";
    std::cout << "    \"wall\": \"steady_clock\"\n";
    std::cout << "  },\n";
    std::cout << "  \"coefficients_low_to_high\": [";
    for (std::size_t i = 0; i < options.coefficients.size(); ++i) {
        if (i != 0) {
            std::cout << ", ";
        }
        write_json_string(std::cout, options.coefficients[i]);
    }
    std::cout << "],\n";
    std::cout << "  \"warmup\": {\n";
    std::cout << "    \"used\": "
              << (options.warmup_coefficients.empty() ? "false" : "true")
              << ",\n";
    std::cout << "    \"degree\": ";
    if (options.warmup_coefficients.empty()) {
        std::cout << "null\n";
    } else {
        std::cout << static_cast<unsigned long long>(
                             options.warmup_coefficients.size() - 1)
                  << "\n";
    }
    std::cout << "  }\n";
    std::cout << "}\n";
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

}  // namespace

int main(int argc, char** argv) {
    Options options;
    std::string error;
    std::string target_nonce;
    if (!parse_options(argc, argv, options, error)) {
        OperationResult result;
        result.error = std::move(error);
        write_result(options, result);
        return 2;
    }

    if (options.marked_protocol && !read_protocol_phase()) {
        OperationResult result;
        result.error = "marked protocol missing warmup phase input";
        write_result(options, result);
        return 5;
    }

    if (!options.warmup_coefficients.empty()) {
        OperationResult warmup =
                run_operation(options.operation, options.warmup_coefficients);
        if (!warmup.success) {
            OperationResult result;
            result.error = "warmup failed: ";
            result.error += warmup.error;
            write_result(options, result);
            return 3;
        }
    }

    if (options.marked_protocol) {
        emit_protocol_marker(ready_marker);
        if (!read_protocol_phase(&target_nonce, target_nonce_bytes) ||
            !protocol_nonce_is_valid(target_nonce)) {
            OperationResult result;
            result.error = "marked protocol target nonce is invalid";
            write_result(options, result);
            return 5;
        }
    }

    OperationResult result =
            run_operation(options.operation, options.coefficients);
    if (options.marked_protocol) {
        emit_protocol_marker(target_done_marker, target_nonce);
        if (!read_protocol_phase()) {
            OperationResult protocol_result;
            protocol_result.error = "marked protocol missing final phase input";
            write_result(options, protocol_result);
            return 5;
        }
    }
    write_result(options, result);
    return result.success ? 0 : 4;
}
