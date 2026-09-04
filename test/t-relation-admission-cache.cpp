#include <silex/class_group.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/relation.hpp>

#include "class_group/relation_admission_cache_internal.hpp"
#include "test_support.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

namespace {
namespace sflint = silex::flint;

constexpr ulong kAdmissionModulus = UWORD(27449);

silex::NumberField degree_one_field() noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, 1);
    return silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
}

struct DegreeOneFixture {
    silex::NumberField field = degree_one_field();
    silex::Order order = silex::test::equation_order(field);
    silex::ClassGroupContext context{order};

    DegreeOneFixture() noexcept {
        sflint::Fmpz bound;
        sflint::fmpz_set_ui(sflint::FmpzRef(bound), 3);
        assert(context.build_factor_base(sflint::FmpzConstRef(bound)));
        assert(context.generator_count() == 2);
    }
};

silex::detail::RelationAdmissionCache empty_cache(
        slong length,
        slong surplus = 0) {
    assert(length >= 0);
    silex::detail::RelationAdmissionCache cache;
    cache.target_relation_count = length + surplus;
    cache.relation_surplus = surplus;
    cache.missing = length;
    cache.modular_basis.assign(
            static_cast<std::size_t>(length * length), UWORD(0));
    return cache;
}

bool set_prime_power_relation(silex::Relation& relation,
                              const silex::FactorBase& base,
                              const silex::NumberField& field,
                              ulong prime,
                              ulong exponent,
                              slong column) noexcept {
    if (column < 0 || column >= base.length()) {
        return false;
    }

    sflint::Fmpz prime_value;
    sflint::Fmpz generator_value;
    sflint::fmpz_set_ui(sflint::FmpzRef(prime_value), prime);
    sflint::fmpz_pow_ui(sflint::FmpzRef(generator_value),
                        sflint::FmpzConstRef(prime_value), exponent);

    silex::Element generator(field);
    sflint::FmpzMat row(1, base.length());
    sflint::fmpz_mat_zero(sflint::FmpzMatRef(row));
    sflint::fmpz_set_ui(
            sflint::fmpz_mat_entry(row, 0, column), exponent);
    return generator.set_fmpz(sflint::FmpzConstRef(generator_value)) &&
           silex::detail::set_relation_from_known_row(
                   relation, base, generator,
                   sflint::FmpzMatConstRef(row));
}

bool prepare_pending_integral_relation(silex::Relation& relation,
                                       slong value,
                                       sflint::FmpzMat& coordinates) noexcept {
    if (sflint::fmpz_mat_nrows(sflint::FmpzMatConstRef(coordinates)) != 1 ||
        sflint::fmpz_mat_ncols(sflint::FmpzMatConstRef(coordinates)) != 1) {
        return false;
    }
    sflint::fmpz_set_si(
            sflint::fmpz_mat_entry(coordinates, 0, 0), value);
    sflint::Fmpq norm;
    sflint::fmpq_set_si(norm, value, 1);
    bool handled = false;
    bool smooth = false;
    return silex::detail::factor_relation_row_from_integral_coordinates_and_norm(
                   relation, handled, smooth,
                   sflint::FmpzMatConstRef(coordinates),
                   sflint::FmpqConstRef(norm)) &&
           handled && smooth && !relation.is_set();
}

bool row_is(const silex::ClassGroupContext& context,
            slong row,
            slong first,
            slong second) noexcept {
    sflint::FmpzMat rows(context.relation_count(),
                        context.generator_count());
    return context.relations(sflint::FmpzMatRef(rows)) &&
           sflint::fmpz_equal_si(
                   sflint::fmpz_mat_entry(rows, row, 0), first) &&
           sflint::fmpz_equal_si(
                   sflint::fmpz_mat_entry(rows, row, 1), second);
}

int test_modular_rank_duplicate_and_surplus_admission() {
    DegreeOneFixture fixture;
    const silex::FactorBase* base = fixture.context.factor_base();
    assert(base != nullptr);
    auto cache = empty_cache(base->length());
    silex::Relation relation(*base);
    bool retained = true;

    // The exact modulus multiple is a nonzero relation, but contributes no
    // modular pivot and is skipped while no surplus is available.
    assert(set_prime_power_relation(
            relation, *base, fixture.field, 2, kAdmissionModulus, 0));
    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, false, true));
    assert(!retained);
    assert(cache.missing == 2);
    assert(cache.skipped_dependent_relation_count == 1);
    assert(fixture.context.relation_count() == 0);

    assert(set_prime_power_relation(
            relation, *base, fixture.field, 2, 1, 0));
    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, false, true));
    assert(retained);
    assert(cache.missing == 1);
    assert(cache.retained_relation_count == 1);
    assert(cache.modular_basis.size() == 4);
    assert(cache.modular_basis[0] == UWORD(1));
    assert(fixture.context.relation_count() == 1);
    assert(fixture.context.relation_rank() == 1);

    cache.relation_row_hashes.clear();
    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, false, true));
    assert(!retained);
    assert(cache.duplicate_relation_count == 1);
    assert(cache.relation_row_hashes.size() == 1);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 1);

    assert(set_prime_power_relation(
            relation, *base, fixture.field, 2, 2, 0));
    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, false, true));
    assert(!retained);
    assert(cache.skipped_dependent_relation_count == 2);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 1);

    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, true, true));
    assert(retained);
    assert(cache.relation_surplus == 0);
    assert(cache.retained_relation_count == 2);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 2);
    assert(fixture.context.relation_rank() == 1);

    cache.relation_surplus = 1;
    assert(set_prime_power_relation(
            relation, *base, fixture.field, 2, 3, 0));
    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, false, true));
    assert(retained);
    assert(cache.relation_surplus == 0);
    assert(cache.retained_relation_count == 3);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 3);
    assert(fixture.context.relation_rank() == 1);

    assert(set_prime_power_relation(
            relation, *base, fixture.field, 3, 1, 1));
    assert(silex::detail::try_admit_relation(
            fixture.context, retained, cache, relation, false, true));
    assert(retained);
    assert(cache.missing == 0);
    assert(cache.retained_relation_count == 4);
    assert(cache.modular_basis[0] == UWORD(1));
    assert(cache.modular_basis[3] == UWORD(1));
    assert(fixture.context.relation_count() == 4);
    assert(fixture.context.relation_rank() == 2);
    return 0;
}

int test_complete_block_initialization_and_pivot_restoration() {
    DegreeOneFixture fixture;
    silex::detail::RelationAdmissionCache cache;
    cache.target_relation_count = 99;
    cache.modular_basis.assign(1, UWORD(99));

    assert(silex::detail::initialize_relation_admission_cache(
            fixture.context, cache, fixture.order, 2));
    assert(cache.target_relation_count == 4);
    assert(cache.relation_surplus == 2);
    assert(cache.relation_count_before_init == 0);
    assert(cache.retained_relation_count == 2);
    assert(cache.trivial_relation_count == 2);
    assert(cache.duplicate_relation_count == 0);
    assert(cache.skipped_dependent_relation_count == 0);
    assert(cache.missing == 0);
    assert(cache.modular_basis.size() == 4);
    assert(cache.relation_row_hashes.size() == 2);
    assert(cache.modular_basis[0] == UWORD(1));
    assert(cache.modular_basis[1] == UWORD(0));
    assert(cache.modular_basis[2] == UWORD(0));
    assert(cache.modular_basis[3] == UWORD(1));
    assert(fixture.context.relation_count() == 2);
    assert(fixture.context.relation_rank() == 2);
    assert(row_is(fixture.context, 0, 1, 0));
    assert(row_is(fixture.context, 1, 0, 1));

    const std::vector<ulong> complete_basis = cache.modular_basis;
    const std::vector<slong> selected{0, 0, 1};
    assert(silex::detail::relation_search::begin_selected_pivot_recollection(
            cache, selected, 2));
    assert(cache.missing == 2);
    assert(cache.modular_basis[0] == UWORD(0));
    assert(cache.modular_basis[3] == UWORD(0));
    silex::detail::relation_search::end_selected_pivot_recollection(
            cache, selected, 2);
    assert(cache.missing == 0);
    assert(cache.modular_basis == complete_basis);
    return 0;
}

int test_deferred_integral_admission_is_atomic_until_retained() {
    DegreeOneFixture fixture;
    const silex::FactorBase* base = fixture.context.factor_base();
    assert(base != nullptr);
    auto cache = empty_cache(base->length());

    silex::Relation relation(*base);
    sflint::FmpzMat coordinates(1, 1);
    assert(prepare_pending_integral_relation(relation, 4, coordinates));
    const std::vector<ulong> initial_basis = cache.modular_basis;
    const slong initial_missing = cache.missing;
    bool retained = true;

    sflint::FmpzMat invalid_coordinates(1, 2);
    assert(!silex::detail::try_admit_deferred_integral_relation(
            fixture.context, retained, cache, relation,
            sflint::FmpzMatConstRef(invalid_coordinates), false, true));
    assert(!retained);
    assert(!relation.is_set());
    assert(fixture.context.relation_count() == 0);
    assert(cache.retained_relation_count == 0);
    assert(cache.missing == initial_missing);
    assert(cache.modular_basis == initial_basis);

    assert(silex::detail::try_admit_deferred_integral_relation(
            fixture.context, retained, cache, relation,
            sflint::FmpzMatConstRef(coordinates), false, true));
    assert(retained);
    assert(relation.is_set());
    assert(cache.retained_relation_count == 1);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 1);
    assert(row_is(fixture.context, 0, 2, 0));
    silex::Element generator(fixture.field);
    assert(fixture.context.relation_generator(generator, 0));
    assert(generator.equal_si(4));

    silex::Relation duplicate(*base);
    sflint::FmpzMat duplicate_coordinates(1, 1);
    assert(prepare_pending_integral_relation(
            duplicate, 4, duplicate_coordinates));
    const std::vector<ulong> basis_before_duplicate = cache.modular_basis;
    assert(silex::detail::try_admit_deferred_integral_relation(
            fixture.context, retained, cache, duplicate,
            sflint::FmpzMatConstRef(duplicate_coordinates), false, true));
    assert(!retained);
    assert(!duplicate.is_set());
    assert(cache.duplicate_relation_count == 1);
    assert(cache.modular_basis == basis_before_duplicate);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 1);

    silex::Relation dependent(*base);
    sflint::FmpzMat dependent_coordinates(1, 1);
    assert(prepare_pending_integral_relation(
            dependent, 8, dependent_coordinates));
    const std::vector<ulong> basis_before_skip = cache.modular_basis;
    assert(silex::detail::try_admit_deferred_integral_relation(
            fixture.context, retained, cache, dependent,
            sflint::FmpzMatConstRef(dependent_coordinates), false, true));
    assert(!retained);
    assert(!dependent.is_set());
    assert(cache.skipped_dependent_relation_count == 1);
    assert(cache.modular_basis == basis_before_skip);
    assert(cache.missing == 1);
    assert(fixture.context.relation_count() == 1);
    return 0;
}

}  // namespace

int main() {
    assert(test_modular_rank_duplicate_and_surplus_admission() == 0);
    assert(test_complete_block_initialization_and_pivot_restoration() == 0);
    assert(test_deferred_integral_admission_is_atomic_until_retained() == 0);
    return 0;
}
