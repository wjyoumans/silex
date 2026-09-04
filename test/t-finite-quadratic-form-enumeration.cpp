#include "class_group/ideal_t2_enumeration_internal.hpp"
#include "class_group/relation_candidate_internal.hpp"

#include <silex/flint/fmpz_mat.hpp>

#include <cstddef>
#include <vector>

namespace {
namespace sflint = silex::flint;
namespace relation_search = silex::detail::relation_search;

bool row_equals(const sflint::FmpzMat& matrix,
                slong row,
                const slong* expected,
                slong columns) noexcept {
    if (expected == nullptr || row < 0 || columns < 0 ||
        sflint::fmpz_mat_nrows(matrix) <= row ||
        sflint::fmpz_mat_ncols(matrix) != columns) {
        return false;
    }
    for (slong column = 0; column < columns; ++column) {
        if (!sflint::fmpz_equal_si(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(matrix), row, column),
                    expected[column])) {
            return false;
        }
    }
    return true;
}

bool set_matrix(sflint::FmpzMat& matrix,
                const slong* entries,
                slong rows,
                slong columns) noexcept {
    if (entries == nullptr || rows < 0 || columns < 0 ||
        sflint::fmpz_mat_nrows(matrix) != rows ||
        sflint::fmpz_mat_ncols(matrix) != columns) {
        return false;
    }
    for (slong row = 0; row < rows; ++row) {
        for (slong column = 0; column < columns; ++column) {
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatRef(matrix), row, column),
                    entries[row * columns + column]);
        }
    }
    return true;
}

int test_reference_recurrence_and_row_orientation() {
    constexpr slong dimension = 4;
    const std::vector<double> quadratic_form{
            19.152046648820685,
            -0.33512222311789463,
            -0.20225622165085556,
            0.49149612536977508,
            0.0,
            17.935588583123879,
            -0.4801719928624274,
            0.3343357376274302,
            0.0,
            0.0,
            15.26713342600937,
            -0.079050565175163263,
            0.0,
            0.0,
            0.0,
            15.595379702497857,
    };
    constexpr double bound = 340.44831774250531;
    constexpr slong raw_rows[][dimension] = {
            {0, 0, 0, 0},
            {1, 0, 0, 0},
            {2, 0, 0, 0},
            {3, 0, 0, 0},
            {4, 0, 0, 0},
            {0, 1, 0, 0},
            {1, 1, 0, 0},
            {-1, 1, 0, 0},
    };
    constexpr slong primitive_rows[][dimension] = {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {1, 1, 0, 0},
            {-1, 1, 0, 0},
    };

    // A deliberately nonsymmetric row basis makes left-row multiplication
    // distinguishable from either column-oriented convention.
    constexpr slong row_basis_entries[dimension * dimension] = {
            2, 1, 0, 0,
            0, 3, 1, 0,
            0, 0, 5, 1,
            1, 0, 0, 7,
    };
    constexpr slong expected_coordinates[][dimension] = {
            {2, 1, 0, 0},
            {0, 3, 1, 0},
            {2, 4, 1, 0},
            {-2, 2, 1, 0},
    };

    relation_search::FiniteQuadraticFormEnumerationContext enumeration;
    sflint::FmpzMat coefficient_row(1, dimension);
    sflint::FmpzMat coordinate_row(1, dimension);
    sflint::FmpzMat row_basis(dimension, dimension);
    if (!set_matrix(row_basis, row_basis_entries, dimension, dimension) ||
        !enumeration.reset(quadratic_form, dimension) ||
        !enumeration.start(bound, 100000, false)) {
        return 1;
    }

    slong primitive_count = 0;
    constexpr slong raw_count =
            static_cast<slong>(sizeof(raw_rows) / sizeof(raw_rows[0]));
    constexpr slong expected_primitive_count = static_cast<slong>(
            sizeof(primitive_rows) / sizeof(primitive_rows[0]));
    for (slong row = 0; row < raw_count; ++row) {
        if (!enumeration.next() ||
            !enumeration.current_row(sflint::FmpzMatRef(coefficient_row)) ||
            !row_equals(coefficient_row, 0, raw_rows[row], dimension)) {
            return 1;
        }
        if (!relation_search::fmpz_mat_single_row_is_primitive(
                    sflint::FmpzMatConstRef(coefficient_row))) {
            continue;
        }
        if (primitive_count >= expected_primitive_count ||
            !row_equals(coefficient_row, 0,
                        primitive_rows[primitive_count], dimension)) {
            return 1;
        }
        relation_search::coordinates_from_lattice_combination(
                coordinate_row, sflint::FmpzMatConstRef(coefficient_row),
                sflint::FmpzMatConstRef(row_basis));
        if (!row_equals(coordinate_row, 0,
                        expected_coordinates[primitive_count], dimension)) {
            return 1;
        }
        ++primitive_count;
    }

    return primitive_count != expected_primitive_count ||
                   enumeration.element_steps() != 12 ||
                   enumeration.element_step_limit_reached()
            ? 1
            : 0;
}

int test_exhaustion_restart_and_element_step_cap() {
    const std::vector<double> quadratic_form{1.0};
    constexpr slong expected_rows[] = {0, 1, 2};
    relation_search::FiniteQuadraticFormEnumerationContext enumeration;
    sflint::FmpzMat row(1, 1);
    if (!enumeration.reset(quadratic_form, 1) ||
        !enumeration.start(4.0, 100, false)) {
        return 1;
    }
    for (slong expected : expected_rows) {
        if (!enumeration.next() ||
            !enumeration.current_row(sflint::FmpzMatRef(row)) ||
            !row_equals(row, 0, &expected, 1)) {
            return 1;
        }
    }
    if (enumeration.next() || enumeration.next() ||
        enumeration.element_steps() != 4 ||
        enumeration.element_step_limit_reached()) {
        return 1;
    }

    constexpr slong zero = 0;
    if (!enumeration.start(4.0, 100, false) || !enumeration.next() ||
        !enumeration.current_row(sflint::FmpzMatRef(row)) ||
        !row_equals(row, 0, &zero, 1) || enumeration.element_steps() != 1 ||
        enumeration.element_step_limit_reached()) {
        return 1;
    }

    if (!enumeration.start(4.0, 1, false) || !enumeration.next() ||
        !enumeration.current_row(sflint::FmpzMatRef(row)) ||
        !row_equals(row, 0, &zero, 1) || enumeration.next() ||
        enumeration.next() || enumeration.element_steps() != 2 ||
        !enumeration.element_step_limit_reached()) {
        return 1;
    }
    return 0;
}

int test_skip_first_scalar_line() {
    constexpr slong dimension = 2;
    const std::vector<double> quadratic_form{
            1.0, 0.0,
            0.0, 1.0,
    };
    constexpr slong expected_rows[][dimension] = {
            {0, 1},
            {1, 1},
            {-1, 1},
            {0, 2},
    };

    relation_search::FiniteQuadraticFormEnumerationContext enumeration;
    sflint::FmpzMat row(1, dimension);
    if (!enumeration.reset(quadratic_form, dimension) ||
        !enumeration.start(4.0, 100, true)) {
        return 1;
    }
    for (const auto& expected : expected_rows) {
        if (!enumeration.next() ||
            !enumeration.current_row(sflint::FmpzMatRef(row)) ||
            !row_equals(row, 0, expected, dimension)) {
            return 1;
        }
    }
    return enumeration.element_steps() != 7 ||
                   enumeration.element_step_limit_reached()
            ? 1
            : 0;
}

}  // namespace

int main() {
    return test_reference_recurrence_and_row_orientation() != 0 ||
                   test_exhaustion_restart_and_element_step_cap() != 0 ||
                   test_skip_first_scalar_line() != 0
            ? 1
            : 0;
}
