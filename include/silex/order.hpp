#pragma once

#include <flint/flint.h>

#include <silex/element.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpz.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>

#include <memory>
#include <optional>
#include <variant>

namespace silex {

struct DiagnosticsContext;
class Ideal;
class Order;

namespace detail {
struct OrderData;
bool order_minkowski_embedding_rows(
        flint::ArbMatRef out,
        const Order& order,
        slong precision,
        const DiagnosticsContext* diagnostics = nullptr) noexcept;
}

struct EquationOrderData {};

struct QuadraticOrderData {
    QuadraticOrderData() noexcept = default;

    explicit QuadraticOrderData(flint::FmpzConstRef conductor_value) noexcept {
        flint::fmpz_set(flint::FmpzRef(conductor), conductor_value);
    }

    QuadraticOrderData(const QuadraticOrderData& other) noexcept {
        flint::fmpz_set(flint::FmpzRef(conductor),
                        flint::FmpzConstRef(other.conductor));
    }

    QuadraticOrderData& operator=(const QuadraticOrderData& other) noexcept {
        if (this != &other) {
            flint::fmpz_set(flint::FmpzRef(conductor),
                            flint::FmpzConstRef(other.conductor));
        }
        return *this;
    }

    QuadraticOrderData(QuadraticOrderData&& other) noexcept = default;
    QuadraticOrderData& operator=(QuadraticOrderData&& other) noexcept = default;

    flint::Fmpz conductor;
};

struct CyclotomicOrderData {
    slong n = 0;
};

using OrderSpecialization = std::variant<std::monostate,
                                         EquationOrderData,
                                         QuadraticOrderData,
                                         CyclotomicOrderData>;

class Order {
public:
    Order() noexcept = default;
    explicit Order(const NumberField& parent) noexcept;
    ~Order() noexcept;

    Order(const Order& other) noexcept;
    Order& operator=(const Order& other) noexcept;

    Order(Order&& other) noexcept;
    Order& operator=(Order&& other) noexcept;

    void swap(Order& other) noexcept;
    void clear() noexcept;
    // Preferred public construction. Invalid input returns an undefined order.
    static Order equation_order(const NumberField& parent) noexcept;
    static Order from_basis(const NumberField& parent,
                            flint::FmpqMatConstRef basis) noexcept;
    static Order quadratic_order(const NumberField& parent,
                                 flint::FmpzConstRef conductor) noexcept;
    // Compatibility/scratch-object construction helpers. Prefer factories in
    // ordinary user code unless mutation/failure preservation is required.
    bool define(const NumberField& parent) noexcept;
    bool define_equation_order(const NumberField& parent) noexcept;
    bool set(const Order& other) noexcept;

    bool is_defined() const noexcept;
    bool has_same_data(const Order& other) const noexcept;
    const NumberField* parent() const noexcept;
    slong degree() const noexcept;
    bool has_basis() const noexcept;
    bool is_equation_order() const noexcept;
    bool maximality_known() const noexcept;
    bool is_maximal() const noexcept;
    void set_maximality(bool is_maximal) noexcept;

    bool set_basis(flint::FmpqMatConstRef basis) noexcept;
    bool get_basis(flint::FmpqMatRef out) const noexcept;
    bool coordinates(flint::FmpqMatRef out, const Element& element) const noexcept;
    std::optional<flint::FmpqMat> basis() const noexcept;
    std::optional<flint::FmpqMat> coordinates(
            const Element& element) const noexcept;
    bool contains(const Element& element) const noexcept;
    bool trace_matrix(flint::FmpzMatRef out) const noexcept;
    bool discriminant(flint::FmpzRef out) const noexcept;
    bool multiplication_table(flint::FmpzMatRef out) noexcept;
    bool quadratic_conductor(flint::FmpzRef out) const noexcept;
    bool index_in(flint::FmpzRef out, const Order& overorder) const noexcept;
    bool p_radical(Ideal& out, flint::FmpzConstRef prime) const noexcept;
    bool pmaximal_overorder(const Order& input, flint::FmpzConstRef prime) noexcept;
    bool maximal_order(const Order& input) noexcept;

private:
    bool set_quadratic_metadata(flint::FmpzConstRef conductor) noexcept;
    bool set_quadratic_order(flint::FmpzConstRef conductor) noexcept;
    void clear_maximality() noexcept;

    friend class Ideal;
    friend bool order_index(flint::FmpzRef out,
                            const Order& suborder,
                            const Order& overorder) noexcept;
    friend bool detail::order_minkowski_embedding_rows(
            flint::ArbMatRef out,
            const Order& order,
            slong precision,
            const DiagnosticsContext* diagnostics) noexcept;

    std::shared_ptr<detail::OrderData> data_;
};

bool order_index(flint::FmpzRef out,
                 const Order& suborder,
                 const Order& overorder) noexcept;

inline bool same_order_parent(const Order* left,
                              const Order* right) noexcept {
    return left != nullptr && right != nullptr && left->has_same_data(*right);
}

inline void swap(Order& left, Order& right) noexcept {
    left.swap(right);
}

}  // namespace silex
