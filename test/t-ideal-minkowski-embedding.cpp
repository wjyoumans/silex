#include "class_group/ideal_minkowski_embedding_internal.hpp"
#include "embedding/embedding_internal.hpp"
#include "order/order_internal.hpp"
#include "test_support.hpp"

#include <silex/embedding.hpp>
#include <silex/flint/acb.hpp>
#include <silex/flint/arb_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include <flint/acb.h>

#include <cassert>

namespace {
namespace sflint = silex::flint;

constexpr slong kDegree = 3;
constexpr slong kInitialPrecision = 128;
constexpr slong kHigherPrecision = 256;

struct FieldSetup {
    silex::NumberField field;
    silex::Order maximal_order;
};

FieldSetup pure_cubic(slong radicand) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, 1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, -radicand);

    FieldSetup setup;
    setup.field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    const silex::Order equation_order =
            silex::test::equation_order(setup.field);
    setup.maximal_order = silex::Order(setup.field);
    assert(setup.maximal_order.maximal_order(equation_order));
    assert(setup.maximal_order.is_maximal());
    assert(setup.maximal_order.degree() == kDegree);
    return setup;
}

void set_nontrivial_basis(sflint::FmpzMat& basis) noexcept {
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(basis));
    const slong entries[kDegree][kDegree] = {
            {2, 1, 0},
            {0, -1, 0},
            {1, 0, 1},
    };
    for (slong row = 0; row < kDegree; ++row) {
        for (slong column = 0; column < kDegree; ++column) {
            sflint::fmpz_set_si(
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatRef(basis), row, column),
                    entries[row][column]);
        }
    }
}

void set_arb_sentinels(sflint::ArbMat& matrix, slong first) noexcept {
    const slong rows = sflint::arb_mat_nrows_value(matrix);
    const slong columns = sflint::arb_mat_ncols_value(matrix);
    for (slong row = 0; row < rows; ++row) {
        for (slong column = 0; column < columns; ++column) {
            sflint::arb_set_si(
                    sflint::arb_mat_entry_ref(matrix, row, column),
                    first + row * columns + column);
        }
    }
}

void copy_arb_matrix(sflint::ArbMat& out,
                     const sflint::ArbMat& input) noexcept {
    assert(sflint::arb_mat_nrows_value(out) ==
           sflint::arb_mat_nrows_value(input));
    assert(sflint::arb_mat_ncols_value(out) ==
           sflint::arb_mat_ncols_value(input));
    sflint::arb_mat_set(sflint::ArbMatRef(out),
                        sflint::ArbMatConstRef(input));
}

bool arb_matrices_equal(const sflint::ArbMat& left,
                        const sflint::ArbMat& right) noexcept {
    return sflint::arb_mat_nrows_value(left) ==
                           sflint::arb_mat_nrows_value(right) &&
           sflint::arb_mat_ncols_value(left) ==
                           sflint::arb_mat_ncols_value(right) &&
           ::arb_mat_equal(left.raw(), right.raw()) != 0;
}

bool arb_matrices_overlap(const sflint::ArbMat& left,
                          const sflint::ArbMat& right) noexcept {
    return sflint::arb_mat_nrows_value(left) ==
                           sflint::arb_mat_nrows_value(right) &&
           sflint::arb_mat_ncols_value(left) ==
                           sflint::arb_mat_ncols_value(right) &&
           ::arb_mat_overlaps(left.raw(), right.raw()) != 0;
}

void integer_matrix_to_arb(sflint::ArbMat& out,
                           const sflint::FmpzMat& input) noexcept {
    assert(sflint::arb_mat_nrows_value(out) ==
           sflint::fmpz_mat_nrows(input));
    assert(sflint::arb_mat_ncols_value(out) ==
           sflint::fmpz_mat_ncols(input));
    for (slong row = 0; row < sflint::fmpz_mat_nrows(input); ++row) {
        for (slong column = 0; column < sflint::fmpz_mat_ncols(input);
             ++column) {
            ::arb_set_fmpz(
                    sflint::arb_mat_entry_ref(out, row, column).raw(),
                    sflint::fmpz_mat_entry(
                            sflint::FmpzMatConstRef(input), row, column)
                            .raw());
        }
    }
}

// The arb_addmul_fmpz loop that multiply_integer_arb_matrices used before it
// summed each entry with arb_dot_fmpz; kept here as a reference enclosure.
void addmul_loop_product(sflint::ArbMat& out,
                         const sflint::FmpzMat& left,
                         const sflint::ArbMat& right,
                         slong precision) noexcept {
    for (slong row = 0; row < sflint::arb_mat_nrows_value(out); ++row) {
        for (slong column = 0; column < sflint::arb_mat_ncols_value(out);
             ++column) {
            arb_struct* entry = arb_mat_entry(out.raw(), row, column);
            ::arb_zero(entry);
            for (slong k = 0; k < sflint::fmpz_mat_ncols(left); ++k) {
                ::arb_addmul_fmpz(entry, arb_mat_entry(right.raw(), k, column),
                                  ::fmpz_mat_entry(left.raw(), row, k),
                                  precision);
            }
        }
    }
}

// right[k][c] = pi * (k + 1) / (c + 2) + sqrt(k + c + 2) at `precision`, or,
// when `dyadic`, the exact value (2k + 1) * 2^(c - 40) + k.
void set_test_right_matrix(sflint::ArbMat& right, bool dyadic,
                           slong precision) noexcept {
    sflint::Arb pi;
    sflint::Arb term;
    ::arb_const_pi(pi.raw(), precision);
    for (slong k = 0; k < sflint::arb_mat_nrows_value(right); ++k) {
        for (slong c = 0; c < sflint::arb_mat_ncols_value(right); ++c) {
            arb_struct* entry = arb_mat_entry(right.raw(), k, c);
            if (dyadic) {
                ::arb_set_si(entry, 2 * k + 1);
                ::arb_mul_2exp_si(entry, entry, c - 40);
                ::arb_add_si(entry, entry, k, ARF_PREC_EXACT);
                continue;
            }
            ::arb_mul_si(entry, pi.raw(), k + 1, precision);
            ::arb_div_si(entry, entry, c + 2, precision);
            ::arb_sqrt_ui(term.raw(), static_cast<ulong>(k + c + 2),
                          precision);
            ::arb_add(entry, entry, term.raw(), precision);
        }
    }
}

// Entries of mixed sign, with some beyond one limb.
void set_test_left_matrix(sflint::FmpzMat& left) noexcept {
    for (slong row = 0; row < sflint::fmpz_mat_nrows(left); ++row) {
        for (slong k = 0; k < sflint::fmpz_mat_ncols(left); ++k) {
            fmpz* entry = ::fmpz_mat_entry(left.raw(), row, k);
            ::fmpz_set_si(entry, (row + 1) * (k % 3 == 0 ? -1 : 1) +
                                         3 * k - 2 * row);
            if ((row + k) % 4 == 1) {
                ::fmpz_mul_2exp(entry, entry, 70);
                ::fmpz_add_ui(entry, entry, 12345);
            }
        }
    }
}

void assert_product_matches_flint(
        const sflint::ArbMat& product,
        const sflint::FmpzMat& left,
        const sflint::ArbMat& right,
        slong precision) {
    sflint::ArbMat left_arb(kDegree, kDegree);
    sflint::ArbMat expected(kDegree, kDegree);
    integer_matrix_to_arb(left_arb, left);
    ::arb_mat_mul(expected.raw(), left_arb.raw(), right.raw(), precision);
    assert(arb_matrices_overlap(product, expected));
}

int test_owner_lifetime_reuse_and_matrix_product() {
    silex::detail::OrderMinkowskiEmbeddingCache cache;
    sflint::ArbMat first_output(kDegree, kDegree);
    {
        FieldSetup setup = pure_cubic(2);
        sflint::FmpzMat identity(kDegree, kDegree);
        sflint::fmpz_mat_one(sflint::FmpzMatRef(identity));
        assert(silex::detail::build_ideal_minkowski_embedding_rows(
                first_output, sflint::FmpzMatConstRef(identity),
                setup.maximal_order, kInitialPrecision, &cache));
        assert(cache.has_order);
        assert(cache.has_order_embedding_rows);
        assert(cache.order.has_same_data(setup.maximal_order));
        assert(cache.embedding_precision == kInitialPrecision);
        assert(sflint::arb_mat_nrows_value(cache.order_embedding_rows) ==
               kDegree);
        assert(sflint::arb_mat_ncols_value(cache.order_embedding_rows) ==
               kDegree);
        assert_product_matches_flint(first_output, identity,
                                     cache.order_embedding_rows,
                                     kInitialPrecision);
    }

    // The cache owns the order data and remains usable after the source
    // field/order setup leaves scope.
    assert(cache.order.is_defined());
    silex::Order first_owner = cache.order;
    sflint::ArbMat cached_rows_before(kDegree, kDegree);
    copy_arb_matrix(cached_rows_before, cache.order_embedding_rows);
    const arb_struct* cached_first_entry =
            arb_mat_entry(cache.order_embedding_rows.raw(), 0, 0);

    sflint::FmpzMat basis(kDegree, kDegree);
    set_nontrivial_basis(basis);
    sflint::ArbMat product(kDegree, kDegree);
    assert(silex::detail::build_ideal_minkowski_embedding_rows(
            product, sflint::FmpzMatConstRef(basis), first_owner,
            kInitialPrecision, &cache));
    assert(cache.order.has_same_data(first_owner));
    assert(cache.embedding_precision == kInitialPrecision);
    assert(cached_first_entry ==
           arb_mat_entry(cache.order_embedding_rows.raw(), 0, 0));
    assert(arb_matrices_equal(cache.order_embedding_rows,
                              cached_rows_before));
    assert_product_matches_flint(product, basis, cached_rows_before,
                                 kInitialPrecision);

    sflint::ArbMat direct_product(kDegree, kDegree);
    assert(silex::detail::multiply_integer_arb_matrices(
            direct_product, sflint::FmpzMatConstRef(basis),
            cached_rows_before, kInitialPrecision));
    assert(arb_matrices_equal(direct_product, product));
    return 0;
}

int test_precision_and_order_misses() {
    FieldSetup first_setup = pure_cubic(2);
    sflint::FmpzMat identity(kDegree, kDegree);
    sflint::fmpz_mat_one(sflint::FmpzMatRef(identity));
    silex::detail::OrderMinkowskiEmbeddingCache cache;
    sflint::ArbMat initial(kDegree, kDegree);
    assert(silex::detail::build_ideal_minkowski_embedding_rows(
            initial, sflint::FmpzMatConstRef(identity),
            first_setup.maximal_order, kInitialPrecision, &cache));
    silex::Order first_owner = cache.order;

    sflint::ArbMat higher_precision(kDegree, kDegree);
    assert(silex::detail::build_ideal_minkowski_embedding_rows(
            higher_precision, sflint::FmpzMatConstRef(identity), first_owner,
            kHigherPrecision, &cache));
    assert(cache.has_order);
    assert(cache.has_order_embedding_rows);
    assert(cache.order.has_same_data(first_owner));
    assert(cache.embedding_precision == kHigherPrecision);
    assert_product_matches_flint(higher_precision, identity,
                                 cache.order_embedding_rows,
                                 kHigherPrecision);

    sflint::ArbMat lower_precision_again(kDegree, kDegree);
    assert(silex::detail::build_ideal_minkowski_embedding_rows(
            lower_precision_again, sflint::FmpzMatConstRef(identity),
            first_owner, kInitialPrecision, &cache));
    assert(cache.has_order);
    assert(cache.has_order_embedding_rows);
    assert(cache.order.has_same_data(first_owner));
    assert(cache.embedding_precision == kInitialPrecision);
    assert_product_matches_flint(lower_precision_again, identity,
                                 cache.order_embedding_rows,
                                 kInitialPrecision);

    cache.order_embedding_rows = sflint::ArbMat(2, 2);
    assert(cache.has_order);
    assert(cache.has_order_embedding_rows);
    sflint::ArbMat repaired_shape(kDegree, kDegree);
    assert(silex::detail::build_ideal_minkowski_embedding_rows(
            repaired_shape, sflint::FmpzMatConstRef(identity), first_owner,
            kInitialPrecision, &cache));
    assert(cache.has_order);
    assert(cache.has_order_embedding_rows);
    assert(cache.order.has_same_data(first_owner));
    assert(cache.embedding_precision == kInitialPrecision);
    assert(sflint::arb_mat_nrows_value(cache.order_embedding_rows) ==
           kDegree);
    assert(sflint::arb_mat_ncols_value(cache.order_embedding_rows) ==
           kDegree);
    assert_product_matches_flint(repaired_shape, identity,
                                 cache.order_embedding_rows,
                                 kInitialPrecision);

    {
        FieldSetup other_setup = pure_cubic(3);
        sflint::ArbMat other_output(kDegree, kDegree);
        assert(!other_setup.maximal_order.has_same_data(first_owner));
        assert(silex::detail::build_ideal_minkowski_embedding_rows(
                other_output, sflint::FmpzMatConstRef(identity),
                other_setup.maximal_order, kHigherPrecision, &cache));
        assert(cache.order.has_same_data(other_setup.maximal_order));
        assert(!cache.order.has_same_data(first_owner));
        assert(cache.embedding_precision == kHigherPrecision);
        assert_product_matches_flint(other_output, identity,
                                     cache.order_embedding_rows,
                                     kHigherPrecision);
    }
    assert(cache.order.is_defined());
    assert(!cache.order.has_same_data(first_owner));
    return 0;
}

int test_invalid_inputs_preserve_outputs_and_cache() {
    FieldSetup setup = pure_cubic(2);
    sflint::FmpzMat identity(kDegree, kDegree);
    sflint::fmpz_mat_one(sflint::FmpzMatRef(identity));
    silex::detail::OrderMinkowskiEmbeddingCache cache;
    sflint::ArbMat valid_output(kDegree, kDegree);
    assert(silex::detail::build_ideal_minkowski_embedding_rows(
            valid_output, sflint::FmpzMatConstRef(identity),
            setup.maximal_order, kInitialPrecision, &cache));

    const silex::Order cache_owner = cache.order;
    const slong cache_precision = cache.embedding_precision;
    sflint::ArbMat cache_rows_before(kDegree, kDegree);
    copy_arb_matrix(cache_rows_before, cache.order_embedding_rows);

    sflint::ArbMat output(kDegree, kDegree);
    set_arb_sentinels(output, 11);
    sflint::ArbMat output_before(kDegree, kDegree);
    copy_arb_matrix(output_before, output);
    assert(!silex::detail::build_ideal_minkowski_embedding_rows(
            output, sflint::FmpzMatConstRef(identity), setup.maximal_order, 0,
            &cache));
    assert(arb_matrices_equal(output, output_before));
    assert(cache.has_order);
    assert(cache.has_order_embedding_rows);
    assert(cache.order.has_same_data(cache_owner));
    assert(cache.embedding_precision == cache_precision);
    assert(arb_matrices_equal(cache.order_embedding_rows, cache_rows_before));

    set_arb_sentinels(output, 31);
    copy_arb_matrix(output_before, output);
    assert(!silex::detail::multiply_integer_arb_matrices(
            output, sflint::FmpzMatConstRef(identity), cache_rows_before, 0));
    assert(arb_matrices_equal(output, output_before));

    silex::detail::OrderMinkowskiEmbeddingCache empty_cache;
    set_arb_sentinels(output, 51);
    copy_arb_matrix(output_before, output);
    assert(!silex::detail::build_ideal_minkowski_embedding_rows(
            output, sflint::FmpzMatConstRef(identity), setup.maximal_order, 0,
            &empty_cache));
    assert(arb_matrices_equal(output, output_before));
    assert(!empty_cache.has_order);
    assert(!empty_cache.has_order_embedding_rows);
    assert(empty_cache.embedding_precision == 0);
    return 0;
}

// The product's balls contain the exact products: right is computed at a low
// precision (or exactly), and the reference product of the same exact values
// at 2048 bits is a ball, around the exact product, inside each entry.  The
// old arb_addmul_fmpz loop's enclosure contains the same reference.
int test_product_contains_exact_values() {
    constexpr slong kLowPrecision = 64;
    constexpr slong kReferencePrecision = 2048;
    const slong shapes[][3] = {{4, 5, 3}, {3, 1, 2}, {2, 0, 3}, {6, 6, 6}};
    for (const auto& shape : shapes) {
        for (const bool dyadic : {false, true}) {
            const slong rows = shape[0];
            const slong middle = shape[1];
            const slong cols = shape[2];
            sflint::FmpzMat left(rows, middle);
            set_test_left_matrix(left);
            sflint::ArbMat right(middle, cols);
            sflint::ArbMat right_reference(middle, cols);
            set_test_right_matrix(right, dyadic, kLowPrecision);
            set_test_right_matrix(right_reference, dyadic,
                                  kReferencePrecision);

            sflint::ArbMat product(rows, cols);
            set_arb_sentinels(product, 7);
            assert(silex::detail::multiply_integer_arb_matrices(
                    product, sflint::FmpzMatConstRef(left), right,
                    kLowPrecision));
            sflint::ArbMat old_product(rows, cols);
            addmul_loop_product(old_product, left, right, kLowPrecision);
            sflint::ArbMat reference(rows, cols);
            addmul_loop_product(reference, left, right_reference,
                                kReferencePrecision);

            for (slong row = 0; row < rows; ++row) {
                for (slong column = 0; column < cols; ++column) {
                    const arb_struct* exact =
                            arb_mat_entry(reference.raw(), row, column);
                    if (dyadic) {
                        assert(::arb_is_exact(exact) != 0);
                    }
                    assert(::arb_contains(
                                   arb_mat_entry(product.raw(), row, column),
                                   exact) != 0);
                    assert(::arb_contains(arb_mat_entry(old_product.raw(),
                                                        row, column),
                                          exact) != 0);
                    if (middle == 0) {
                        assert(::arb_is_zero(arb_mat_entry(
                                       product.raw(), row, column)) != 0);
                    }
                }
            }
        }
    }
    return 0;
}

}  // namespace

// The order's own embedding context, read through OrderAccess, gains roots
// when the order's Minkowski rows are built, and a separate context for the
// same field can be seeded from it with the same certified roots.
int test_order_embedding_context_seeds_other_context() {
    FieldSetup setup = pure_cubic(2);
    const silex::EmbeddingContext* order_embeddings =
            silex::detail::OrderAccess::minkowski_embedding_context(
                    setup.maximal_order);
    assert(order_embeddings != nullptr);
    assert(!order_embeddings->is_set());

    const slong precision = 264;
    sflint::ArbMat rows(kDegree, kDegree);
    assert(silex::detail::order_minkowski_embedding_rows(
            sflint::ArbMatRef(rows), setup.maximal_order, precision));
    assert(silex::detail::OrderAccess::minkowski_embedding_context(
                   setup.maximal_order) == order_embeddings);
    assert(order_embeddings->is_set());
    assert(order_embeddings->precision() >= precision);

    silex::EmbeddingContext seeded(*setup.maximal_order.parent());
    assert(silex::detail::EmbeddingContextAccess::seed_roots(
            seeded, *order_embeddings));
    assert(seeded.precision() == order_embeddings->precision());
    assert(seeded.refine(128));
    for (slong i = 0; i < kDegree; ++i) {
        sflint::Acb seeded_root;
        sflint::Acb order_root;
        assert(seeded.get_root(sflint::AcbRef(seeded_root), i));
        assert(order_embeddings->get_root(sflint::AcbRef(order_root), i));
        assert(::acb_equal(seeded_root.raw(), order_root.raw()) != 0);
        // The root of x^3 - 2 lies in the ball.
        sflint::Acb cube;
        ::acb_pow_ui(cube.raw(), seeded_root.raw(), 3, precision);
        ::acb_sub_ui(cube.raw(), cube.raw(), 2, precision);
        assert(::acb_contains_zero(cube.raw()) != 0);
    }

    silex::Order undefined;
    assert(silex::detail::OrderAccess::minkowski_embedding_context(
                   undefined) == nullptr);
    return 0;
}

int main() {
    assert(test_owner_lifetime_reuse_and_matrix_product() == 0);
    assert(test_precision_and_order_misses() == 0);
    assert(test_invalid_inputs_preserve_outputs_and_cache() == 0);
    assert(test_product_contains_exact_values() == 0);
    assert(test_order_embedding_context_seeds_other_context() == 0);
    return 0;
}
