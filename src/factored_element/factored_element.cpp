#include <silex/factored_element.hpp>

#include <flint/arb.h>

#include <limits>
#include <utility>

namespace silex {
namespace {

bool copy_factor(Element& out, const Element& input) noexcept {
    if (input.parent() == nullptr) {
        return false;
    }
    out = Element(*input.parent());
    return out.is_defined() && out.set(input);
}

bool negate_slong(slong& out, slong value) noexcept {
    if (value == std::numeric_limits<slong>::min()) {
        return false;
    }
    out = -value;
    return true;
}

bool multiply_slong(slong& out, slong left, slong right) noexcept {
    const slong min = std::numeric_limits<slong>::min();
    const slong max = std::numeric_limits<slong>::max();

    if (left == 0 || right == 0) {
        out = 0;
        return true;
    }
    if ((left == -1 && right == min) ||
        (right == -1 && left == min)) {
        return false;
    }

    if (left > 0) {
        if (right > 0) {
            if (left > max / right) {
                return false;
            }
        } else if (right < min / left) {
            return false;
        }
    } else {
        if (right > 0) {
            if (left < min / right) {
                return false;
            }
        } else if (left < max / right) {
            return false;
        }
    }

    out = left * right;
    return true;
}

flint_bitcnt_t abs_slong_bits(slong value) noexcept {
    if (value == 0) {
        return 0;
    }

    ulong magnitude = 0;
    if (value < 0) {
        magnitude = static_cast<ulong>(-(value + 1));
        ++magnitude;
    } else {
        magnitude = static_cast<ulong>(value);
    }

    flint_bitcnt_t bits = 0;
    while (magnitude != 0) {
        ++bits;
        magnitude >>= 1;
    }
    return bits;
}

flint_bitcnt_t size_bits(std::size_t value) noexcept {
    flint_bitcnt_t bits = 0;
    while (value != 0) {
        ++bits;
        value >>= 1;
    }
    return bits;
}

bool factored_log_work_precision(slong& out,
                                 slong target_precision,
                                 flint_bitcnt_t max_exponent_bits,
                                 flint_bitcnt_t factor_count_bits) noexcept {
    const slong max = std::numeric_limits<slong>::max();
    if (max_exponent_bits > static_cast<flint_bitcnt_t>(max) ||
        factor_count_bits > static_cast<flint_bitcnt_t>(max)) {
        return false;
    }

    const slong exponent_guard = static_cast<slong>(max_exponent_bits);
    const slong count_guard = static_cast<slong>(factor_count_bits);
    if (target_precision > max - exponent_guard ||
        target_precision + exponent_guard > max - count_guard) {
        return false;
    }

    out = target_precision + exponent_guard + count_guard;
    return true;
}

bool factored_log_bucket_precision(slong& out,
                                   slong requested_precision) noexcept {
    const slong max = std::numeric_limits<slong>::max();
    slong bucket = 1;
    while (bucket <= requested_precision) {
        if (bucket > max / 2) {
            return false;
        }
        bucket *= 2;
    }
    out = bucket;
    return true;
}

bool arb_has_radius_at_most_2exp(const arb_struct* value,
                                slong exponent) noexcept {
    return arb_is_finite(value) != 0 &&
           mag_cmp_2exp_si(arb_radref(value), exponent) <= 0;
}

void expand_arb_to_radius(arb_struct* value, slong max_radius_2exp) noexcept {
    if (arb_rel_accuracy_bits(value) < 0 || arb_is_exact(value) != 0) {
        return;
    }

    slong precision = arb_bits(value) / 2;
    if (precision < 2) {
        return;
    }

    flint::Arb rounded;
    arb_set_round(rounded.raw(), value, precision);
    while (arb_has_radius_at_most_2exp(rounded.raw(), max_radius_2exp) &&
           precision > 4) {
        arb_set(value, rounded.raw());
        precision /= 2;
        arb_set_round(rounded.raw(), rounded.raw(), precision);
    }
}

void profile_factored_log_shape(const DiagnosticsContext* diagnostics,
                                std::size_t factor_count,
                                flint_bitcnt_t max_exponent_bits) noexcept {
    if (factor_count <= 1) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.factor_count.le1");
    } else if (factor_count <= 4) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.factor_count.le4");
    } else {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.factor_count.gt4");
    }

    if (max_exponent_bits <= 4) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.max_exp_bits.le4");
    } else if (max_exponent_bits <= 8) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.max_exp_bits.le8");
    } else if (max_exponent_bits <= 16) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.max_exp_bits.le16");
    } else if (max_exponent_bits <= 32) {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.max_exp_bits.le32");
    } else {
        SILEX_PROFILE_EVENT(
                diagnostics, DiagnosticsModule::element,
                "element.factored.logarithmic_embedding.max_exp_bits.gt32");
    }
}

slong floor_half_slong(slong value) noexcept {
    slong quotient = value / 2;
    if (value < 0 && value % 2 != 0) {
        --quotient;
    }
    return quotient;
}

bool is_odd_slong(slong value) noexcept {
    return value % 2 != 0;
}

void floor_divrem_positive_slong(
        slong& quotient,
        slong& remainder,
        slong value,
        slong divisor) noexcept {
    quotient = value / divisor;
    remainder = value % divisor;
    if (remainder < 0) {
        --quotient;
        remainder += divisor;
    }
}

bool split_formal_power(FactoredElement& quotient,
                        FactoredElement& remainder,
                        const FactoredElement& input,
                        slong exponent) noexcept {
    if (!quotient.is_defined() || !remainder.is_defined() ||
        input.parent() == nullptr ||
        quotient.parent() == nullptr ||
        remainder.parent() == nullptr ||
        !quotient.parent()->has_same_data(*input.parent()) ||
        !remainder.parent()->has_same_data(*input.parent()) ||
        exponent <= 0 || !quotient.one() || !remainder.one()) {
        return false;
    }

    for (const auto& entry : input.factors()) {
        slong q = 0;
        slong r = 0;
        floor_divrem_positive_slong(q, r, entry.exponent, exponent);
        if (q != 0 && !quotient.push(entry.factor, q)) {
            return false;
        }
        if (r != 0 && !remainder.push(entry.factor, r)) {
            return false;
        }
    }

    quotient.normalize();
    remainder.normalize();
    return true;
}

bool evaluate_entries(Element& out,
                            const NumberField& parent,
                            FactorSpan<Element> entries) noexcept {
    if (!out.has_parent(parent)) {
        return false;
    }

    if (entries.empty()) {
        Element one(parent);
        return one.one() && out.set(one);
    }

    if (entries.size() == 1) {
        flint::Fmpz exponent;
        fmpz_set_si(exponent.raw(), entries.front().exponent);
        Element power(parent);
        return power.pow_fmpz(entries.front().factor,
                              flint::FmpzConstRef(exponent)) &&
               out.set(power);
    }

    Element accumulator(parent);
    Element power(parent);
    Element product(parent);
    if (!accumulator.one()) {
        return false;
    }

    // reference `evaluate(FacElem)` recursively halves large exponents and
    // multiplies the odd residual factors before squaring the recursive part.
    FactoredElement halved(parent);
    for (const auto& entry : entries) {
        if (entry.exponent > -10 && entry.exponent < 10) {
            if (entry.exponent == 0) {
                continue;
            }
            flint::Fmpz exponent;
            fmpz_set_si(exponent.raw(), entry.exponent);
            if (!power.pow_fmpz(entry.factor,
                                flint::FmpzConstRef(exponent)) ||
                !product.multiply(accumulator, power) ||
                !accumulator.set(product)) {
                return false;
            }
            continue;
        }

        const slong residual = is_odd_slong(entry.exponent) ? 1 : 0;
        const slong divided = (entry.exponent - residual) / 2;
        if (divided != 0 && !halved.push(entry.factor, divided)) {
            return false;
        }
        if (residual != 0 &&
            (!product.multiply(accumulator, entry.factor) ||
             !accumulator.set(product))) {
            return false;
        }
    }

    if (halved.length() == 0) {
        return out.set(accumulator);
    }

    Element half_value(parent);
    Element half_square(parent);
    Element result(parent);
    if (!evaluate_entries(half_value, parent, halved.factors()) ||
        !half_square.multiply(half_value, half_value) ||
        !result.multiply(half_square, accumulator)) {
        return false;
    }
    return out.set(result);
}

}  // namespace

FactoredElement::FactoredElement(const NumberField& parent) noexcept {
    define(parent);
}

FactoredElement::~FactoredElement() noexcept = default;

FactoredElement::FactoredElement(FactoredElement&& other) noexcept {
    swap(other);
}

FactoredElement& FactoredElement::operator=(
        FactoredElement&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FactoredElement::swap(FactoredElement& other) noexcept {
    parent_.swap(other.parent_);
    std::swap(value_, other.value_);
}

void FactoredElement::clear() noexcept {
    value_.clear();
    parent_.clear();
}

bool FactoredElement::define(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    clear();
    parent_ = parent;
    return true;
}

bool FactoredElement::push_copy(
        const Element& factor,
        slong exponent) noexcept {
    if (!is_defined() || !factor.has_parent(parent_) || factor.equal_si(0)) {
        return false;
    }
    if (exponent == 0) {
        return true;
    }

    Element copy;
    if (!copy_factor(copy, factor)) {
        return false;
    }
    value_.push(std::move(copy), exponent);
    return true;
}

bool FactoredElement::set(const FactoredElement& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }
    if (!is_defined()) {
        *this = FactoredElement(other.parent_);
        if (!is_defined()) {
            return false;
        }
    }
    if (!parent_.has_same_data(other.parent_)) {
        return false;
    }

    FactoredElement copy(parent_);
    copy.value_.reserve(other.value_.factors().size());
    for (const auto& entry : other.value_.factors()) {
        if (!copy.push_copy(entry.factor, entry.exponent)) {
            return false;
        }
    }

    swap(copy);
    return true;
}

bool FactoredElement::is_defined() const noexcept {
    return parent_.is_defined();
}

const NumberField* FactoredElement::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong FactoredElement::length() const noexcept {
    return is_defined() ? static_cast<slong>(value_.factors().size()) : 0;
}

FactorSpan<Element> FactoredElement::factors() const noexcept {
    return is_defined() ? value_.factors() : FactorSpan<Element>();
}

bool FactoredElement::one() noexcept {
    if (!is_defined()) {
        return false;
    }
    value_.clear();
    return true;
}

bool FactoredElement::set_element(const Element& element) noexcept {
    if (!is_defined() || !element.has_parent(parent_) || element.equal_si(0)) {
        return false;
    }

    FactoredElement next(parent_);
    if (!next.push_copy(element, 1)) {
        return false;
    }
    swap(next);
    return true;
}

bool FactoredElement::push(
        const Element& factor,
        slong exponent) noexcept {
    return push_copy(factor, exponent);
}

void FactoredElement::normalize() noexcept {
    if (is_defined()) {
        value_.normalize();
    }
}

bool FactoredElement::multiply(
        const FactoredElement& left,
        const FactoredElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_)) {
        return false;
    }

    FactoredElement result(parent_);
    result.value_.reserve(left.value_.factors().size() +
                          right.value_.factors().size());
    for (const auto& entry : left.value_.factors()) {
        if (!result.push_copy(entry.factor, entry.exponent)) {
            return false;
        }
    }
    for (const auto& entry : right.value_.factors()) {
        if (!result.push_copy(entry.factor, entry.exponent)) {
            return false;
        }
    }
    result.normalize();
    swap(result);
    return true;
}

bool FactoredElement::divide(
        const FactoredElement& left,
        const FactoredElement& right) noexcept {
    if (!is_defined() || !left.is_defined() || !right.is_defined() ||
        !parent_.has_same_data(left.parent_) ||
        !left.parent_.has_same_data(right.parent_)) {
        return false;
    }

    FactoredElement result(parent_);
    result.value_.reserve(left.value_.factors().size() +
                          right.value_.factors().size());
    for (const auto& entry : left.value_.factors()) {
        if (!result.push_copy(entry.factor, entry.exponent)) {
            return false;
        }
    }
    for (const auto& entry : right.value_.factors()) {
        slong negated = 0;
        if (!negate_slong(negated, entry.exponent) ||
            !result.push_copy(entry.factor, negated)) {
            return false;
        }
    }
    result.normalize();
    swap(result);
    return true;
}

bool FactoredElement::invert(const FactoredElement& input) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.has_same_data(input.parent_)) {
        return false;
    }

    FactoredElement result(parent_);
    result.value_.reserve(input.value_.factors().size());
    for (const auto& entry : input.value_.factors()) {
        slong negated = 0;
        if (!negate_slong(negated, entry.exponent) ||
            !result.push_copy(entry.factor, negated)) {
            return false;
        }
    }
    result.normalize();
    swap(result);
    return true;
}

bool FactoredElement::pow_si(
        const FactoredElement& input,
        slong exponent) noexcept {
    if (!is_defined() || !input.is_defined() ||
        !parent_.has_same_data(input.parent_)) {
        return false;
    }

    FactoredElement result(parent_);
    result.value_.reserve(input.value_.factors().size());
    for (const auto& entry : input.value_.factors()) {
        slong scaled = 0;
        if (!multiply_slong(scaled, entry.exponent, exponent)) {
            return false;
        }
        if (scaled != 0 && !result.push_copy(entry.factor, scaled)) {
            return false;
        }
    }
    result.normalize();
    swap(result);
    return true;
}

bool FactoredElement::root_si(
        const FactoredElement& input,
        slong exponent,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "factored.root_si");
    if (!is_defined() || !input.is_defined() || exponent <= 0 ||
        !parent_.has_same_data(input.parent_)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::element, LogLevel::detail,
                  "factored root rejected");
        return false;
    }

    FactoredElement result(parent_);
    result.value_.reserve(input.value_.factors().size());
    for (const auto& entry : input.value_.factors()) {
        if (entry.exponent % exponent != 0) {
            return false;
        }

        const slong divided = entry.exponent / exponent;
        if (divided != 0 && !result.push_copy(entry.factor, divided)) {
            return false;
        }
    }
    result.normalize();
    swap(result);
    return true;
}

bool FactoredElement::is_square(
        bool& is_square,
        FactoredElement& root,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "factored.is_square");
    if (!is_defined() || !root.is_defined() ||
        !parent_.has_same_data(root.parent_)) {
        return false;
    }

    FactoredElement structural(parent_);
    FactoredElement leftover(parent_);
    structural.value_.reserve(value_.factors().size());
    leftover.value_.reserve(value_.factors().size());
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.square_split");
        for (const auto& entry : value_.factors()) {
            const slong quotient = floor_half_slong(entry.exponent);
            if (quotient != 0 &&
                !structural.push_copy(entry.factor, quotient)) {
                return false;
            }
            if (is_odd_slong(entry.exponent) &&
                !leftover.push_copy(entry.factor, 1)) {
                return false;
            }
        }
    }

    if (leftover.length() == 0) {
        structural.normalize();
        is_square = true;
        root.swap(structural);
        return true;
    }

    Element leftover_value(parent_);
    Element leftover_root(parent_);
    bool factor_is_square = false;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.square_leftover_evaluate");
        if (!leftover.evaluate(leftover_value)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.square_leftover_is_square");
        if (!leftover_value.is_square(factor_is_square, leftover_root,
                                      diagnostics)) {
            return false;
        }
    }

    if (!factor_is_square) {
        is_square = false;
        return true;
    }

    FactoredElement factor_root(parent_);
    FactoredElement candidate(parent_);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.square_combine_root");
        if (!factor_root.set_element(leftover_root) ||
            !candidate.multiply(structural, factor_root)) {
            return false;
        }
    }

    is_square = true;
    root.swap(candidate);
    return true;
}

bool FactoredElement::is_power_si(
        bool& is_power,
        FactoredElement& root,
        slong exponent,
        const DiagnosticsContext* diagnostics) const noexcept {
    return is_power_si(
            is_power, root, exponent, FactoredRootStrategy::automatic,
            diagnostics);
}

bool FactoredElement::is_power_si(
        bool& is_power,
        FactoredElement& root,
        slong exponent,
        FactoredRootStrategy strategy,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "factored.is_power_si");
    if (!is_defined() || !root.is_defined() || exponent <= 0 ||
        !parent_.has_same_data(root.parent_)) {
        SILEX_LOG(diagnostics, DiagnosticsModule::element, LogLevel::detail,
                  "factored power check rejected");
        return false;
    }

    if (strategy == FactoredRootStrategy::automatic ||
        strategy == FactoredRootStrategy::reduced) {
        strategy = FactoredRootStrategy::reduced;
    } else if (strategy != FactoredRootStrategy::compact) {
        return false;
    }

    if (exponent == 1) {
        FactoredElement candidate(parent_);
        if (!candidate.set(*this)) {
            return false;
        }
        is_power = true;
        root.swap(candidate);
        return true;
    }

    if (strategy == FactoredRootStrategy::compact) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "factored.compact");
        FactoredElement structural(parent_);
        FactoredElement leftover(parent_);
        // reference's factored-element is_power first removes the formal
        // exponent quotient and only constructs a compact presentation for
        // the residual factors.
        {
            SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                                "factored.compact_split");
            if (!split_formal_power(structural, leftover, *this, exponent)) {
                return false;
            }
        }
        if (leftover.length() == 0) {
            SILEX_PROFILE_EVENT(
                    diagnostics, DiagnosticsModule::element,
                    "factored.compact_structural_root");
            is_power = true;
            root.swap(structural);
            return true;
        }

        CompactElement compact(parent_);
        if (!compact.set_factored_element(leftover, exponent, diagnostics)) {
            return false;
        }
        FactoredElement leftover_root(parent_);
        const CompactElement::RootStatus status =
                compact.root_base_status(leftover_root, diagnostics);
        if (status == CompactElement::RootStatus::unsupported) {
            return false;
        }
        if (status == CompactElement::RootStatus::not_power) {
            is_power = false;
            return true;
        }

        FactoredElement candidate(parent_);
        if (!candidate.multiply(structural, leftover_root)) {
            return false;
        }
        is_power = true;
        root.swap(candidate);
        return true;
    }

    if (exponent == 2) {
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::element,
                            "factored.square_root");
        return is_square(is_power, root, diagnostics);
    }

    FactoredElement structural(parent_);
    FactoredElement leftover(parent_);
    structural.value_.reserve(value_.factors().size());
    leftover.value_.reserve(value_.factors().size());
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.reduced_split");
        for (const auto& entry : value_.factors()) {
            slong quotient = 0;
            slong remainder = 0;
            floor_divrem_positive_slong(
                    quotient, remainder, entry.exponent, exponent);
            if (quotient != 0 &&
                !structural.push_copy(entry.factor, quotient)) {
                return false;
            }
            if (remainder != 0 &&
                !leftover.push_copy(entry.factor, remainder)) {
                return false;
            }
        }
    }

    if (leftover.length() == 0) {
        structural.normalize();
        is_power = true;
        root.swap(structural);
        return true;
    }

    Element leftover_value(parent_);
    Element leftover_root(parent_);
    flint::Fmpz exponent_fmpz;
    fmpz_set_si(exponent_fmpz.raw(), exponent);
    bool factor_is_power = false;
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.reduced_leftover_evaluate");
        if (!leftover.evaluate(leftover_value)) {
            return false;
        }
    }
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.reduced_leftover_is_power");
        if (!leftover_value.is_power(factor_is_power, leftover_root,
                                     flint::FmpzConstRef(exponent_fmpz),
                                     diagnostics)) {
            return false;
        }
    }

    if (!factor_is_power) {
        is_power = false;
        return true;
    }

    FactoredElement factor_root(parent_);
    FactoredElement candidate(parent_);
    {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                            "factored.reduced_combine_root");
        if (!factor_root.set_element(leftover_root) ||
            !candidate.multiply(structural, factor_root)) {
            return false;
        }
    }

    is_power = true;
    root.swap(candidate);
    return true;
}

const Element* FactoredElement::factor(slong index) const noexcept {
    const auto entries = factors();
    if (index < 0 || index >= static_cast<slong>(entries.size())) {
        return nullptr;
    }
    return &entries[static_cast<std::size_t>(index)].factor;
}

bool FactoredElement::exponent(slong& out, slong index) const noexcept {
    const auto entries = factors();
    if (index < 0 || index >= static_cast<slong>(entries.size())) {
        return false;
    }
    out = entries[static_cast<std::size_t>(index)].exponent;
    return true;
}

bool FactoredElement::evaluate(Element& out) const noexcept {
    if (!is_defined() || !out.has_parent(parent_)) {
        return false;
    }

    return evaluate_entries(out, parent_, value_.factors());
}

bool FactoredElement::logarithmic_embedding(
        flint::ArbVecRef out,
        EmbeddingContext& embeddings,
        LogEmbeddingMode mode,
        slong precision,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "element.factored.logarithmic_embedding");
    if (!is_defined() || embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(parent_) ||
        precision <= 0 ||
        (mode != LogEmbeddingMode::plain &&
         mode != LogEmbeddingMode::product)) {
        return false;
    }

    Signature sig;
    if (embeddings.is_set()) {
        sig = embeddings.signature();
    } else if (!signature(sig, parent_)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    if (out.length() != places) {
        return false;
    }

    const auto entries = value_.factors();
    flint_bitcnt_t max_exponent_bits = 0;
    for (const auto& entry : entries) {
        const flint_bitcnt_t bits = abs_slong_bits(entry.exponent);
        if (bits > max_exponent_bits) {
            max_exponent_bits = bits;
        }
    }
    profile_factored_log_shape(diagnostics, entries.size(),
                               max_exponent_bits);

    flint::ArbVec sum(places);
    flint::ArbVec row(places);
    if (entries.empty()) {
        _arb_vec_zero(sum.data(), places);
        return out.set_from(flint::ArbVecConstRef(sum));
    }

    // Source trace: reference NfOrd/FacElem.jl:116-120 and :149-205 bucket factor
    // logs, add exponent and factor-count guard bits, and double working
    // precision until every partial sum meets the requested absolute error.
    slong work_precision = 0;
    if (!factored_log_work_precision(
                work_precision, precision, max_exponent_bits,
                size_bits(entries.size()))) {
        return false;
    }

    const slong target_radius_2exp = -precision;
    for (;;) {
        slong factor_precision = 0;
        if (!factored_log_bucket_precision(factor_precision,
                                           work_precision)) {
            return false;
        }
        _arb_vec_zero(sum.data(), places);
        bool first = true;
        bool precision_too_low = false;
        for (const auto& entry : entries) {
            {
                SILEX_PROFILE_SCOPE(
                        diagnostics, DiagnosticsModule::element,
                        "element.factored.logarithmic_embedding_factor");
                if (!silex::logarithmic_embedding(
                            flint::ArbVecRef(row), embeddings, entry.factor,
                            mode, factor_precision, diagnostics)) {
                    precision_too_low = true;
                    break;
                }
            }

            for (slong j = 0; j < places; ++j) {
                if (first) {
                    arb_mul_si(sum.data() + j, row.data() + j,
                               entry.exponent, factor_precision);
                } else {
                    arb_addmul_si(sum.data() + j, row.data() + j,
                                  entry.exponent, factor_precision);
                }
                if (!arb_has_radius_at_most_2exp(
                            sum.data() + j, target_radius_2exp)) {
                    precision_too_low = true;
                    break;
                }
            }
            if (precision_too_low) {
                break;
            }
            first = false;
        }

        if (!precision_too_low) {
            for (slong j = 0; j < places; ++j) {
                expand_arb_to_radius(sum.data() + j,
                                     target_radius_2exp);
            }
            return out.set_from(flint::ArbVecConstRef(sum));
        }

        if (work_precision > std::numeric_limits<slong>::max() / 2) {
            return false;
        }
        work_precision *= 2;
    }
}

CompactElement::CompactElement(const NumberField& parent) noexcept {
    define(parent);
}

CompactElement::~CompactElement() noexcept = default;

CompactElement::CompactElement(CompactElement&& other) noexcept {
    swap(other);
}

CompactElement& CompactElement::operator=(
        CompactElement&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void CompactElement::swap(CompactElement& other) noexcept {
    parent_.swap(other.parent_);
    std::swap(base_, other.base_);
    coefficients_.swap(other.coefficients_);
}

void CompactElement::clear() noexcept {
    coefficients_.clear();
    base_ = 0;
    parent_.clear();
}

bool CompactElement::define(const NumberField& parent) noexcept {
    if (!parent.is_defined()) {
        return false;
    }

    clear();
    parent_ = parent;
    return true;
}

bool CompactElement::is_defined() const noexcept {
    return parent_.is_defined();
}

const NumberField* CompactElement::parent() const noexcept {
    return is_defined() ? &parent_ : nullptr;
}

slong CompactElement::base() const noexcept {
    return is_defined() ? base_ : 0;
}

slong CompactElement::length() const noexcept {
    return is_defined() ? static_cast<slong>(coefficients_.size()) : 0;
}

bool CompactElement::ensure_length(slong length) noexcept {
    if (!is_defined() || length < 0) {
        return false;
    }

    while (static_cast<slong>(coefficients_.size()) < length) {
        FactoredElement coefficient(parent_);
        if (!coefficient.is_defined()) {
            return false;
        }
        coefficients_.push_back(std::move(coefficient));
    }
    return true;
}

void CompactElement::normalize() noexcept {
    if (!is_defined()) {
        return;
    }
    while (!coefficients_.empty() && coefficients_.back().length() == 0) {
        coefficients_.pop_back();
    }
}

bool CompactElement::one(slong base) noexcept {
    if (!is_defined() || base <= 1) {
        return false;
    }
    base_ = base;
    coefficients_.clear();
    return true;
}

bool CompactElement::set_factored_element(
        const FactoredElement& input,
        slong base,
        const DiagnosticsContext* diagnostics) noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "compact.set_factored_element");
    if (!is_defined() || input.parent() == nullptr ||
        !input.parent()->has_same_data(parent_) || base <= 1) {
        SILEX_LOG(diagnostics, DiagnosticsModule::element, LogLevel::detail,
                  "compact representation rejected");
        return false;
    }

    CompactElement result(parent_);
    if (!result.one(base)) {
        return false;
    }

    for (const auto& entry : input.factors()) {
        slong exponent = entry.exponent;
        slong level = 0;
        while (exponent != 0) {
            const slong quotient = exponent / base;
            const slong remainder = exponent % base;
            if (remainder != 0) {
                if (!result.ensure_length(level + 1) ||
                    !result.coefficients_[
                            static_cast<std::size_t>(level)].push(
                            entry.factor, remainder)) {
                    return false;
                }
            }
            exponent = quotient;
            ++level;
        }
    }

    result.normalize();
    swap(result);
    return true;
}

bool CompactElement::evaluate(Element& out) const noexcept {
    if (!is_defined() || !out.has_parent(parent_) || base_ <= 1) {
        return false;
    }

    Element accumulator(parent_);
    Element coefficient_value(parent_);
    Element powered(parent_);
    Element product(parent_);
    flint::Fmpz power;
    if (!accumulator.one()) {
        return false;
    }
    fmpz_one(power.raw());

    for (const auto& coefficient : coefficients_) {
        if (!coefficient.evaluate(coefficient_value) ||
            !powered.pow_fmpz(coefficient_value,
                              flint::FmpzConstRef(power)) ||
            !product.multiply(accumulator, powered) ||
            !accumulator.set(product)) {
            return false;
        }
        fmpz_mul_si(power.raw(), power.raw(), base_);
    }

    return out.set(accumulator);
}

bool CompactElement::logarithmic_embedding(
        flint::ArbVecRef out,
        EmbeddingContext& embeddings,
        LogEmbeddingMode mode,
        slong precision) const noexcept {
    if (!is_defined() || embeddings.parent() == nullptr ||
        !embeddings.parent()->has_same_data(parent_) ||
        base_ <= 1 || precision <= 0 ||
        (mode != LogEmbeddingMode::plain &&
         mode != LogEmbeddingMode::product)) {
        return false;
    }

    Signature sig;
    if (embeddings.is_set()) {
        sig = embeddings.signature();
    } else if (!signature(sig, parent_)) {
        return false;
    }
    const slong places = sig.r1() + sig.r2();
    if (out.length() != places) {
        return false;
    }

    flint::ArbVec sum(places);
    flint::ArbVec row(places);
    flint::Fmpz power;
    _arb_vec_zero(sum.data(), places);
    fmpz_one(power.raw());

    for (const auto& coefficient : coefficients_) {
        if (!coefficient.logarithmic_embedding(
                    flint::ArbVecRef(row), embeddings, mode, precision)) {
            return false;
        }
        for (slong j = 0; j < places; ++j) {
            arb_addmul_fmpz(sum.data() + j, row.data() + j,
                            power.raw(), precision);
        }
        fmpz_mul_si(power.raw(), power.raw(), base_);
    }

    return out.set_from(flint::ArbVecConstRef(sum));
}

bool CompactElement::root_base(
        FactoredElement& root,
        const DiagnosticsContext* diagnostics) const noexcept {
    return root_base_status(root, diagnostics) == RootStatus::power;
}

CompactElement::RootStatus CompactElement::root_base_status(
        FactoredElement& root,
        const DiagnosticsContext* diagnostics) const noexcept {
    SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::element,
                        "compact.root_base");
    if (!is_defined() || root.parent() == nullptr ||
        !root.parent()->has_same_data(parent_) || base_ <= 1) {
        SILEX_LOG(diagnostics, DiagnosticsModule::element, LogLevel::detail,
                  "compact root rejected");
        return RootStatus::unsupported;
    }

    FactoredElement result(parent_);
    if (coefficients_.empty()) {
        if (!result.one()) {
            return RootStatus::unsupported;
        }
        root.swap(result);
        return RootStatus::power;
    }

    bool lead_is_power = false;
    if (!coefficients_[0].is_power_si(lead_is_power, result, base_,
                                      diagnostics)) {
        return RootStatus::unsupported;
    }
    if (!lead_is_power) {
        return RootStatus::not_power;
    }

    FactoredElement term(parent_);
    FactoredElement product(parent_);
    slong power = 1;
    for (std::size_t i = 1; i < coefficients_.size(); ++i) {
        if (!term.pow_si(coefficients_[i], power) ||
            !product.multiply(result, term)) {
            return RootStatus::unsupported;
        }
        result.swap(product);
        if (i + 1 < coefficients_.size() &&
            !multiply_slong(power, power, base_)) {
            return RootStatus::unsupported;
        }
    }

    root.swap(result);
    return RootStatus::power;
}

}  // namespace silex
