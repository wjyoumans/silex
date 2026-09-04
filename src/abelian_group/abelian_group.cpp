#include <silex/abelian_group.hpp>

#include <utility>

namespace silex {
namespace {

flint::FmpzConstRef entry_const(const flint::FmpzMat& matrix,
                                slong row,
                                slong col) noexcept {
    return flint::fmpz_mat_entry(flint::FmpzMatConstRef(matrix), row, col);
}

bool snf_full_column_rank(flint::FmpzMatConstRef snf, slong columns) noexcept {
    if (flint::fmpz_mat_nrows(snf) < columns ||
        flint::fmpz_mat_ncols(snf) < columns) {
        return false;
    }

    for (slong i = 0; i < columns; ++i) {
        if (flint::fmpz_is_zero(flint::fmpz_mat_entry(snf, i, i))) {
            return false;
        }
    }
    return true;
}

slong trimmed_nonzero_rows(const fmpz_mat_t matrix) noexcept {
    slong rows = fmpz_mat_nrows(matrix);
    while (rows > 0 && fmpz_mat_is_zero_row(matrix, rows - 1) != 0) {
        --rows;
    }
    for (slong i = 0; i < rows; ++i) {
        if (fmpz_mat_is_zero_row(matrix, i) != 0) {
            return -1;
        }
    }
    return rows;
}

bool inverse_unimodular(flint::FmpzMat& inverse,
                        flint::FmpzMatConstRef matrix) noexcept {
    flint::Fmpz den;
    if (!flint::fmpz_mat_inv(flint::FmpzMatRef(inverse),
                             flint::FmpzRef(den), matrix)) {
        return false;
    }

    if (flint::fmpz_equal_si(den, -1)) {
        flint::fmpz_mat_neg(flint::FmpzMatRef(inverse),
                            flint::FmpzMatConstRef(inverse));
    } else if (!flint::fmpz_is_one(den)) {
        return false;
    }

    return true;
}

bool set_invariants(flint::FmpzVec& invariants,
                    flint::FmpzMatConstRef snf,
                    slong columns) noexcept {
    slong invariant_count = 0;
    for (slong i = 0; i < columns; ++i) {
        if (!flint::fmpz_is_pm1(flint::fmpz_mat_entry(snf, i, i))) {
            ++invariant_count;
        }
    }

    flint::FmpzVec next(invariant_count);
    slong index = 0;
    for (slong i = 0; i < columns; ++i) {
        if (flint::fmpz_is_pm1(flint::fmpz_mat_entry(snf, i, i))) {
            continue;
        }

        flint::fmpz_abs(flint::FmpzRef(next.data() + index),
                        flint::fmpz_mat_entry(snf, i, i));
        ++index;
    }

    invariants.swap(next);
    return true;
}

}  // namespace

FiniteAbelianGroup::~FiniteAbelianGroup() noexcept = default;

FiniteAbelianGroup::FiniteAbelianGroup(FiniteAbelianGroup&& other) noexcept {
    swap(other);
}

FiniteAbelianGroup& FiniteAbelianGroup::operator=(
        FiniteAbelianGroup&& other) noexcept {
    if (this != &other) {
        clear();
        swap(other);
    }
    return *this;
}

void FiniteAbelianGroup::swap(FiniteAbelianGroup& other) noexcept {
    relations_.swap(other.relations_);
    snf_.swap(other.snf_);
    left_transform_.swap(other.left_transform_);
    right_transform_.swap(other.right_transform_);
    right_transform_inv_.swap(other.right_transform_inv_);
    invariants_.swap(other.invariants_);
    std::swap(generator_count_, other.generator_count_);
    std::swap(defined_, other.defined_);
    std::swap(has_left_transform_, other.has_left_transform_);
}

void FiniteAbelianGroup::clear() noexcept {
    relations_ = flint::FmpzMat(0, 0);
    snf_ = flint::FmpzMat(0, 0);
    left_transform_ = flint::FmpzMat(0, 0);
    right_transform_ = flint::FmpzMat(0, 0);
    right_transform_inv_ = flint::FmpzMat(0, 0);
    invariants_ = flint::FmpzVec(0);
    generator_count_ = 0;
    defined_ = false;
    has_left_transform_ = false;
}

bool FiniteAbelianGroup::set(const FiniteAbelianGroup& other) noexcept {
    if (this == &other) {
        return true;
    }
    if (!other.is_defined()) {
        clear();
        return true;
    }

    FiniteAbelianGroup copy;
    copy.relations_ = flint::FmpzMat(other.relation_count(),
                                     other.generator_count());
    copy.snf_ = flint::FmpzMat(flint::fmpz_mat_nrows(other.snf_),
                               other.generator_count());
    copy.left_transform_ = flint::FmpzMat(
            flint::fmpz_mat_nrows(other.left_transform_),
            flint::fmpz_mat_ncols(other.left_transform_));
    copy.right_transform_ = flint::FmpzMat(other.generator_count(),
                                           other.generator_count());
    copy.right_transform_inv_ = flint::FmpzMat(other.generator_count(),
                                               other.generator_count());
    copy.invariants_ = flint::FmpzVec(other.invariant_count());

    flint::fmpz_mat_set(flint::FmpzMatRef(copy.relations_),
                        flint::FmpzMatConstRef(other.relations_));
    flint::fmpz_mat_set(flint::FmpzMatRef(copy.snf_),
                        flint::FmpzMatConstRef(other.snf_));
    flint::fmpz_mat_set(flint::FmpzMatRef(copy.left_transform_),
                        flint::FmpzMatConstRef(other.left_transform_));
    flint::fmpz_mat_set(flint::FmpzMatRef(copy.right_transform_),
                        flint::FmpzMatConstRef(other.right_transform_));
    flint::fmpz_mat_set(
            flint::FmpzMatRef(copy.right_transform_inv_),
            flint::FmpzMatConstRef(other.right_transform_inv_));
    for (slong i = 0; i < other.invariant_count(); ++i) {
        flint::fmpz_set(flint::FmpzRef(copy.invariants_.data() + i),
                        flint::FmpzConstRef(other.invariants_.data() + i));
    }
    copy.generator_count_ = other.generator_count_;
    copy.defined_ = true;
    copy.has_left_transform_ = other.has_left_transform_;

    swap(copy);
    return true;
}

bool FiniteAbelianGroup::is_defined() const noexcept {
    return defined_ && flint::fmpz_mat_ncols(relations_) == generator_count_ &&
           flint::fmpz_mat_nrows(snf_) >= generator_count_ &&
           flint::fmpz_mat_ncols(snf_) == generator_count_ &&
           (!has_left_transform_ ||
            (flint::fmpz_mat_nrows(left_transform_) == relation_count() &&
             flint::fmpz_mat_ncols(left_transform_) == relation_count())) &&
           flint::fmpz_mat_nrows(right_transform_) == generator_count_ &&
           flint::fmpz_mat_ncols(right_transform_) == generator_count_ &&
           flint::fmpz_mat_nrows(right_transform_inv_) == generator_count_ &&
           flint::fmpz_mat_ncols(right_transform_inv_) == generator_count_;
}

bool FiniteAbelianGroup::set_relation_matrix(
        flint::FmpzMatConstRef relations) noexcept {
    const slong rows = flint::fmpz_mat_nrows(relations);
    const slong columns = flint::fmpz_mat_ncols(relations);
    if (rows < columns) {
        return false;
    }

    FiniteAbelianGroup candidate;
    candidate.relations_ = flint::FmpzMat(rows, columns);
    candidate.snf_ = flint::FmpzMat(columns, columns);
    candidate.left_transform_ = flint::FmpzMat(0, 0);
    candidate.right_transform_ = flint::FmpzMat(columns, columns);
    candidate.right_transform_inv_ = flint::FmpzMat(columns, columns);
    candidate.generator_count_ = columns;

    flint::fmpz_mat_set(flint::FmpzMatRef(candidate.relations_), relations);

    // reference implementations use HNF/index checkpoints during relation collection and
    // defer full relation-map construction.  Publish invariants from an HNF
    // basis of the relation lattice; compute the full left transform lazily
    // only when witness rows are requested.
    flint::FmpzMat hnf_full(rows, columns);
    ::fmpz_mat_hnf(hnf_full.raw(), relations.raw());
    const slong hnf_rows = trimmed_nonzero_rows(hnf_full.raw());
    if (hnf_rows != columns) {
        return false;
    }

    flint::FmpzMat hnf_basis(columns, columns);
    {
        flint::FmpzMatConstWindow hnf_window(
                hnf_full, 0, 0, columns, columns);
        flint::fmpz_mat_set(flint::FmpzMatRef(hnf_basis),
                            hnf_window.const_ref());
    }

    return set_relation_matrix_with_hnf_basis(relations,
                                              flint::FmpzMatConstRef(hnf_basis));
}

bool FiniteAbelianGroup::set_relation_matrix_with_hnf_basis(
        flint::FmpzMatConstRef relations,
        flint::FmpzMatConstRef hnf_basis) noexcept {
    const slong rows = flint::fmpz_mat_nrows(relations);
    const slong columns = flint::fmpz_mat_ncols(relations);
    if (rows < columns ||
        flint::fmpz_mat_nrows(hnf_basis) != columns ||
        flint::fmpz_mat_ncols(hnf_basis) != columns) {
        return false;
    }

    FiniteAbelianGroup candidate;
    candidate.relations_ = flint::FmpzMat(rows, columns);
    candidate.snf_ = flint::FmpzMat(columns, columns);
    candidate.left_transform_ = flint::FmpzMat(0, 0);
    candidate.right_transform_ = flint::FmpzMat(columns, columns);
    candidate.right_transform_inv_ = flint::FmpzMat(columns, columns);
    candidate.generator_count_ = columns;

    flint::fmpz_mat_set(flint::FmpzMatRef(candidate.relations_), relations);

    flint::FmpzMat left_reduced(columns, columns);
    flint::fmpz_mat_snf_transform(
            flint::FmpzMatRef(candidate.snf_),
            flint::FmpzMatRef(left_reduced),
            flint::FmpzMatRef(candidate.right_transform_),
            hnf_basis);

    if (!snf_full_column_rank(flint::FmpzMatConstRef(candidate.snf_), columns) ||
        !inverse_unimodular(
                candidate.right_transform_inv_,
                flint::FmpzMatConstRef(candidate.right_transform_)) ||
        !set_invariants(candidate.invariants_,
                        flint::FmpzMatConstRef(candidate.snf_), columns)) {
        return false;
    }
    candidate.defined_ = true;
    candidate.has_left_transform_ = false;

    swap(candidate);
    return true;
}

bool FiniteAbelianGroup::ensure_left_transform() const noexcept {
    if (!is_defined()) {
        return false;
    }
    if (has_left_transform_) {
        return true;
    }

    // Expensive fallback: this reconstructs the original transformed SNF
    // relation map for APIs that need witnesses in terms of stored rows.
    flint::FmpzMat full_snf(relation_count(), generator_count());
    flint::FmpzMat full_left(relation_count(), relation_count());
    flint::FmpzMat full_right(generator_count(), generator_count());
    flint::FmpzMat full_right_inv(generator_count(), generator_count());
    flint::FmpzVec full_invariants(0);

    flint::fmpz_mat_snf_transform(
            flint::FmpzMatRef(full_snf),
            flint::FmpzMatRef(full_left),
            flint::FmpzMatRef(full_right),
            flint::FmpzMatConstRef(relations_));

    if (!snf_full_column_rank(flint::FmpzMatConstRef(full_snf),
                              generator_count()) ||
        !inverse_unimodular(full_right_inv,
                            flint::FmpzMatConstRef(full_right)) ||
        !set_invariants(full_invariants, flint::FmpzMatConstRef(full_snf),
                        generator_count())) {
        return false;
    }

    snf_.swap(full_snf);
    left_transform_.swap(full_left);
    right_transform_.swap(full_right);
    right_transform_inv_.swap(full_right_inv);
    invariants_.swap(full_invariants);
    has_left_transform_ = true;
    return true;
}

slong FiniteAbelianGroup::relation_count() const noexcept {
    return defined_ ? flint::fmpz_mat_nrows(relations_) : 0;
}

slong FiniteAbelianGroup::generator_count() const noexcept {
    return defined_ ? generator_count_ : 0;
}

slong FiniteAbelianGroup::invariant_count() const noexcept {
    return defined_ ? invariants_.length() : 0;
}

slong FiniteAbelianGroup::relation_kernel_count() const noexcept {
    if (!is_defined() || relation_count() <= generator_count()) {
        return 0;
    }
    return relation_count() - generator_count();
}

bool FiniteAbelianGroup::relations(flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || flint::fmpz_mat_nrows(out) != relation_count() ||
        flint::fmpz_mat_ncols(out) != generator_count()) {
        return false;
    }
    flint::fmpz_mat_set(out, flint::FmpzMatConstRef(relations_));
    return true;
}

std::optional<flint::FmpzMat> FiniteAbelianGroup::relations() const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpzMat out(relation_count(), generator_count());
    if (!relations(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool FiniteAbelianGroup::invariant(flint::FmpzRef out,
                                   slong index) const noexcept {
    if (!is_defined() || index < 0 || index >= invariant_count()) {
        return false;
    }
    flint::fmpz_set(out, flint::FmpzConstRef(invariants_.data() + index));
    return true;
}

std::optional<flint::Fmpz> FiniteAbelianGroup::invariant(
        slong index) const noexcept {
    flint::Fmpz out;
    if (!invariant(flint::FmpzRef(out), index)) {
        return std::nullopt;
    }
    return out;
}

bool FiniteAbelianGroup::invariants(flint::FmpzVecRef out) const noexcept {
    if (!is_defined() || out.length() != invariant_count()) {
        return false;
    }

    for (slong i = 0; i < invariant_count(); ++i) {
        flint::fmpz_set(flint::FmpzRef(out.data() + i),
                        flint::FmpzConstRef(invariants_.data() + i));
    }
    return true;
}

bool FiniteAbelianGroup::order(flint::FmpzRef out) const noexcept {
    if (!is_defined()) {
        return false;
    }

    flint::Fmpz product;
    flint::fmpz_one(flint::FmpzRef(product));
    for (slong i = 0; i < invariant_count(); ++i) {
        flint::fmpz_mul(flint::FmpzRef(product),
                        flint::FmpzConstRef(product),
                        flint::FmpzConstRef(invariants_.data() + i));
    }
    flint::fmpz_set(out, flint::FmpzConstRef(product));
    return true;
}

std::optional<flint::Fmpz> FiniteAbelianGroup::order() const noexcept {
    flint::Fmpz out;
    if (!order(flint::FmpzRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool FiniteAbelianGroup::invariant_generator_matrix(
        flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || flint::fmpz_mat_nrows(out) != invariant_count() ||
        flint::fmpz_mat_ncols(out) != generator_count()) {
        return false;
    }

    slong k = 0;
    for (slong i = 0; i < generator_count(); ++i) {
        if (flint::fmpz_is_pm1(entry_const(snf_, i, i))) {
            continue;
        }

        for (slong j = 0; j < generator_count(); ++j) {
            flint::fmpz_set(
                    flint::fmpz_mat_entry(out, k, j),
                    entry_const(right_transform_inv_, i, j));
        }
        ++k;
    }
    return true;
}

std::optional<flint::FmpzMat>
FiniteAbelianGroup::invariant_generator_matrix() const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpzMat out(invariant_count(), generator_count());
    if (!invariant_generator_matrix(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool FiniteAbelianGroup::invariant_generator_relation_matrix(
        flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !ensure_left_transform() ||
        flint::fmpz_mat_nrows(out) != invariant_count() ||
        flint::fmpz_mat_ncols(out) != relation_count()) {
        return false;
    }

    slong k = 0;
    for (slong i = 0; i < generator_count(); ++i) {
        if (flint::fmpz_is_pm1(entry_const(snf_, i, i))) {
            continue;
        }

        for (slong j = 0; j < relation_count(); ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(out, k, j),
                            entry_const(left_transform_, i, j));
            if (flint::fmpz_sgn(entry_const(snf_, i, i)) < 0) {
                flint::fmpz_neg(flint::fmpz_mat_entry(out, k, j),
                                flint::FmpzConstRef(
                                        flint::fmpz_mat_entry(out, k, j).raw()));
            }
        }
        ++k;
    }
    return true;
}

bool FiniteAbelianGroup::relation_kernel_matrix(
        flint::FmpzMatRef out) const noexcept {
    if (!is_defined() || !ensure_left_transform() ||
        flint::fmpz_mat_nrows(out) != relation_kernel_count() ||
        flint::fmpz_mat_ncols(out) != relation_count()) {
        return false;
    }

    for (slong i = 0; i < relation_kernel_count(); ++i) {
        const slong source_row = generator_count() + i;
        for (slong j = 0; j < relation_count(); ++j) {
            flint::fmpz_set(flint::fmpz_mat_entry(out, i, j),
                            entry_const(left_transform_, source_row, j));
        }
    }
    return true;
}

std::optional<flint::FmpzMat>
FiniteAbelianGroup::relation_kernel_matrix() const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpzMat out(relation_kernel_count(), relation_count());
    if (!relation_kernel_matrix(flint::FmpzMatRef(out))) {
        return std::nullopt;
    }
    return out;
}

bool FiniteAbelianGroup::relation_kernel_row(
        flint::FmpzMatRef out,
        slong index) const noexcept {
    if (!is_defined() || !ensure_left_transform() ||
        index < 0 || index >= relation_kernel_count() ||
        flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != relation_count()) {
        return false;
    }

    const slong source_row = generator_count() + index;
    for (slong j = 0; j < relation_count(); ++j) {
        flint::fmpz_set(flint::fmpz_mat_entry(out, 0, j),
                        entry_const(left_transform_, source_row, j));
    }
    return true;
}

bool FiniteAbelianGroup::reduce(flint::FmpzMatRef row) const noexcept {
    if (!is_defined() || flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != generator_count()) {
        return false;
    }

    flint::FmpzMat smith_row(1, generator_count());
    flint::FmpzMat reduced(1, generator_count());
    flint::Fmpz modulus;
    flint::fmpz_mat_mul(flint::FmpzMatRef(smith_row),
                        flint::FmpzMatConstRef(row.raw()),
                        flint::FmpzMatConstRef(right_transform_));

    for (slong i = 0; i < generator_count(); ++i) {
        if (flint::fmpz_is_pm1(entry_const(snf_, i, i))) {
            continue;
        }

        flint::fmpz_abs(flint::FmpzRef(modulus),
                        entry_const(snf_, i, i));
        flint::fmpz_fdiv_r(flint::fmpz_mat_entry(smith_row, 0, i),
                           flint::FmpzConstRef(
                                   flint::fmpz_mat_entry(smith_row, 0, i).raw()),
                           flint::FmpzConstRef(modulus));
    }

    flint::fmpz_mat_mul(flint::FmpzMatRef(reduced),
                        flint::FmpzMatConstRef(smith_row),
                        flint::FmpzMatConstRef(right_transform_inv_));
    flint::fmpz_mat_set(row, flint::FmpzMatConstRef(reduced));
    return true;
}

bool FiniteAbelianGroup::invariant_coordinates(
        flint::FmpzMatRef out,
        flint::FmpzMatConstRef row) const noexcept {
    if (!is_defined() || flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) != generator_count() ||
        flint::fmpz_mat_nrows(out) != 1 ||
        flint::fmpz_mat_ncols(out) != invariant_count()) {
        return false;
    }

    flint::FmpzMat smith_row(1, generator_count());
    flint::FmpzMat candidate(1, invariant_count());
    flint::Fmpz modulus;
    flint::fmpz_mat_mul(flint::FmpzMatRef(smith_row), row,
                        flint::FmpzMatConstRef(right_transform_));

    slong k = 0;
    for (slong i = 0; i < generator_count(); ++i) {
        if (flint::fmpz_is_pm1(entry_const(snf_, i, i))) {
            continue;
        }

        flint::fmpz_abs(flint::FmpzRef(modulus),
                        entry_const(snf_, i, i));
        flint::fmpz_fdiv_r(flint::fmpz_mat_entry(candidate, 0, k),
                           flint::FmpzConstRef(
                                   flint::fmpz_mat_entry(smith_row, 0, i).raw()),
                           flint::FmpzConstRef(modulus));
        ++k;
    }

    flint::fmpz_mat_set(out, flint::FmpzMatConstRef(candidate));
    return true;
}

std::optional<flint::FmpzMat> FiniteAbelianGroup::invariant_coordinates(
        flint::FmpzMatConstRef row) const noexcept {
    if (!is_defined()) {
        return std::nullopt;
    }
    flint::FmpzMat out(1, invariant_count());
    if (!invariant_coordinates(flint::FmpzMatRef(out), row)) {
        return std::nullopt;
    }
    return out;
}

}  // namespace silex
