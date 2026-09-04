#include <silex/class_group.hpp>

#include <silex/prime_ideal.hpp>
#include <silex/relation.hpp>

#include "relation_admission_cache_internal.hpp"
#include "../factor_base/factor_base_internal.hpp"

#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

namespace silex {
namespace detail::relation_search {
namespace {

constexpr ulong kRelationAdmissionModulus = UWORD(27449);

bool set_complete_prime_block_relation(
        Relation& relation,
        const FactorBase& base,
        const Order& order,
        flint::FmpzConstRef rational_prime,
        const DiagnosticsContext* diagnostics) noexcept {
    const NumberField* field = order.parent();
    if (field == nullptr) {
        return false;
    }

    Element generator(*field);
    Relation candidate(base);
    if (!generator.is_defined() || !candidate.is_defined() ||
        !generator.set_fmpz(rational_prime) ||
        !candidate.set_generator(generator, diagnostics)) {
        return false;
    }

    relation.swap(candidate);
    return true;
}

bool relation_mod_row(std::vector<ulong>& out,
                      slong& first_nonzero,
                      flint::FmpzMatConstRef row) noexcept;

bool relation_row_already_known(const RelationMatrix& matrix,
                                flint::FmpzMatConstRef row,
                                slong first_nonzero) noexcept {
    if (!matrix.is_defined() ||
        flint::fmpz_mat_nrows(row) != 1 ||
        matrix.ncols() != flint::fmpz_mat_ncols(row) ||
        first_nonzero < 0 || first_nonzero > matrix.ncols()) {
        return false;
    }

    for (slong i = 0; i < matrix.length(); ++i) {
        slong stored_first_nonzero = matrix.ncols();
        if (first_nonzero < matrix.ncols() &&
            (!matrix.row_first_nonzero(stored_first_nonzero, i) ||
             stored_first_nonzero != first_nonzero)) {
            continue;
        }
        if (matrix.row_equal(row, i)) {
            return true;
        }
    }
    return false;
}

ulong mix_relation_row_hash(ulong hash, ulong value) noexcept {
    hash ^= value;
    hash *= UWORD(1099511628211);
    return hash;
}

ulong relation_mod_row_hash(const std::vector<ulong>& row) noexcept {
    ulong hash = UWORD(1469598103934665603);
    hash = mix_relation_row_hash(
            hash, static_cast<ulong>(row.size()));
    for (std::size_t col = 0; col < row.size(); ++col) {
        hash = mix_relation_row_hash(hash,
                                          static_cast<ulong>(col + 1));
        hash = mix_relation_row_hash(hash, row[col]);
    }
    return hash;
}

bool sync_relation_row_hashes(
        detail::RelationAdmissionCache& cache,
        const RelationMatrix& matrix) noexcept {
    if (!matrix.is_defined()) {
        return false;
    }

    cache.relation_row_hashes.clear();
    cache.relation_row_hashes.reserve(
            static_cast<std::size_t>(matrix.length()));
    flint::FmpzMat row(1, matrix.ncols());
    std::vector<ulong> mod_row;
    for (slong i = 0; i < matrix.length(); ++i) {
        slong first_nonzero = matrix.ncols();
        if (!matrix.row(flint::FmpzMatRef(row), i) ||
            !relation_mod_row(mod_row, first_nonzero,
                              flint::FmpzMatConstRef(row))) {
            cache.relation_row_hashes.clear();
            return false;
        }
        cache.relation_row_hashes.push_back(
                relation_mod_row_hash(mod_row));
    }
    return true;
}

bool relation_row_already_known_with_cache(
        bool& known,
        detail::RelationAdmissionCache& cache,
        const RelationMatrix& matrix,
        flint::FmpzMatConstRef row,
        slong first_nonzero,
        ulong row_hash,
        const DiagnosticsContext* diagnostics) noexcept {
    known = false;
    if (!matrix.is_defined()) {
        return false;
    }

    if (cache.relation_row_hashes.size() !=
        static_cast<std::size_t>(matrix.length())) {
        SILEX_PROFILE_SCOPE(diagnostics, DiagnosticsModule::class_group,
                            "class_group.relation_admission.duplicate_hash_sync");
        if (!sync_relation_row_hashes(cache, matrix)) {
            return false;
        }
    }

    if (cache.relation_row_hashes.size() !=
        static_cast<std::size_t>(matrix.length())) {
        known = relation_row_already_known(matrix, row, first_nonzero);
        return true;
    }

    for (slong i = 0; i < matrix.length(); ++i) {
        if (cache.relation_row_hashes[static_cast<std::size_t>(i)] !=
            row_hash) {
            continue;
        }
        SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                            "class_group.relation_admission.duplicate_hash_hit");
        slong stored_first_nonzero = matrix.ncols();
        if (first_nonzero < matrix.ncols() &&
            (!matrix.row_first_nonzero(stored_first_nonzero, i) ||
             stored_first_nonzero != first_nonzero)) {
            continue;
        }
        if (matrix.row_equal(row, i)) {
            known = true;
            return true;
        }
    }

    SILEX_PROFILE_EVENT(diagnostics, DiagnosticsModule::class_group,
                        "class_group.relation_admission.duplicate_hash_skip");
    return true;
}

bool append_trace_detail(char* buffer,
                         std::size_t size,
                         std::size_t& used,
                         const char* format,
                         ...) noexcept {
    if (used >= size || size == 0) {
        return false;
    }

    va_list args;
    va_start(args, format);
    const int written =
            std::vsnprintf(buffer + used, size - used, format, args);
    va_end(args);
    if (written < 0) {
        return false;
    }

    const std::size_t count = static_cast<std::size_t>(written);
    if (count >= size - used) {
        used = size - 1;
        buffer[used] = '\0';
        return false;
    }
    used += count;
    return true;
}

void append_relation_row_trace_detail(char* detail,
                                      std::size_t size,
                                      std::size_t& used,
                                      flint::FmpzMatConstRef row) noexcept {
    if (flint::fmpz_mat_nrows(row) != 1 ||
        flint::fmpz_mat_ncols(row) < 0) {
        append_trace_detail(detail, size, used, " row=<unavailable>");
        return;
    }

    const slong length = flint::fmpz_mat_ncols(row);
    append_trace_detail(detail, size, used, " row=[");
    const slong shown = length < 12 ? length : 12;
    for (slong i = 0; i < shown; ++i) {
        flint::FmpzConstRef entry =
                flint::fmpz_mat_entry(row, 0, i);
        if (i != 0) {
            append_trace_detail(detail, size, used, ",");
        }
        if (flint::fmpz_fits_si(entry)) {
            append_trace_detail(detail, size, used, "%ld",
                                static_cast<long>(flint::fmpz_get_si(entry)));
        } else {
            append_trace_detail(detail, size, used, "?");
        }
    }
    if (shown < length) {
        append_trace_detail(detail, size, used, ",...");
    }
    append_trace_detail(detail, size, used, "]");
}

void append_relation_row_trace_detail(char* detail,
                                      std::size_t size,
                                      std::size_t& used,
                                      const Relation& relation) noexcept {
    if (!relation.is_set() || relation.length() < 0) {
        append_trace_detail(detail, size, used, " row=<unset>");
        return;
    }

    flint::FmpzMat row(1, relation.length());
    if (!relation.exponents(flint::FmpzMatRef(row))) {
        append_trace_detail(detail, size, used, " row=<unavailable>");
        return;
    }

    append_relation_row_trace_detail(detail, size, used,
                                     flint::FmpzMatConstRef(row));
}

void log_relation_admission_trace(const DiagnosticsContext* diagnostics,
                                 const char* outcome,
                                 slong first_nonzero,
                                 bool modular_progress,
                                 const detail::RelationAdmissionCache& cache,
                                 slong relation_count,
                                 bool in_random_relation,
                                 flint::FmpzMatConstRef row) noexcept {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (!log_enabled(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::trace)) {
        return;
    }
    char detail[512] = {};
    std::size_t used = 0;
    append_trace_detail(detail, sizeof(detail), used,
                        "outcome=%s first_nonzero=%ld modular_progress=%d "
                        "random=%d missing=%ld relsup=%ld relations=%ld",
                        outcome, static_cast<long>(first_nonzero),
                        modular_progress ? 1 : 0,
                        in_random_relation ? 1 : 0,
                        static_cast<long>(cache.missing),
                        static_cast<long>(cache.relation_surplus),
                        static_cast<long>(relation_count));
    append_relation_row_trace_detail(detail, sizeof(detail), used, row);
    log_emit(diagnostics, DiagnosticsModule::class_group, LogLevel::trace,
             __func__, "relation admission", detail);
#else
    (void)diagnostics;
    (void)outcome;
    (void)first_nonzero;
    (void)modular_progress;
    (void)cache;
    (void)relation_count;
    (void)in_random_relation;
    (void)row;
#endif
}

void log_relation_admission_trace(const DiagnosticsContext* diagnostics,
                                 const char* outcome,
                                 slong first_nonzero,
                                 bool modular_progress,
                                 const detail::RelationAdmissionCache& cache,
                                 slong relation_count,
                                 bool in_random_relation,
                                 const Relation& relation) noexcept {
#if defined(SILEX_ENABLE_LOGGING) && SILEX_ENABLE_LOGGING
    if (!log_enabled(diagnostics, DiagnosticsModule::class_group,
                     LogLevel::trace)) {
        return;
    }
    char detail[512] = {};
    std::size_t used = 0;
    append_trace_detail(detail, sizeof(detail), used,
                        "outcome=%s first_nonzero=%ld modular_progress=%d "
                        "random=%d missing=%ld relsup=%ld relations=%ld",
                        outcome, static_cast<long>(first_nonzero),
                        modular_progress ? 1 : 0,
                        in_random_relation ? 1 : 0,
                        static_cast<long>(cache.missing),
                        static_cast<long>(cache.relation_surplus),
                        static_cast<long>(relation_count));
    append_relation_row_trace_detail(detail, sizeof(detail), used, relation);
    log_emit(diagnostics, DiagnosticsModule::class_group, LogLevel::trace,
             __func__, "relation admission", detail);
#else
    (void)diagnostics;
    (void)outcome;
    (void)first_nonzero;
    (void)modular_progress;
    (void)cache;
    (void)relation_count;
    (void)in_random_relation;
    (void)relation;
#endif
}

ulong& relation_basis_entry(detail::RelationAdmissionCache& cache,
                              slong row,
                              slong column,
                              slong n) noexcept {
    return cache.modular_basis[
            static_cast<std::size_t>(column * n + row)];
}

bool relation_basis_size(std::size_t& size, slong n) noexcept {
    if (n < 0) {
        return false;
    }
    const std::size_t width = static_cast<std::size_t>(n);
    if (width != 0 &&
        width > std::numeric_limits<std::size_t>::max() / width) {
        return false;
    }
    size = width * width;
    return true;
}

bool relation_basis_is_sized(
        const detail::RelationAdmissionCache& cache,
        slong n) noexcept {
    std::size_t size = 0;
    return relation_basis_size(size, n) &&
           cache.modular_basis.size() == size;
}

bool modular_basis_row_has_progress(bool& progress,
                                         const std::vector<ulong>& basis,
                                         slong missing,
                                         const std::vector<ulong>& row)
        noexcept {
    progress = false;
    if (row.size() >
        static_cast<std::size_t>(std::numeric_limits<slong>::max())) {
        return false;
    }
    const slong n = static_cast<slong>(row.size());
    std::size_t basis_size = 0;
    if (!relation_basis_size(basis_size, n) ||
        basis.size() != basis_size || missing < 0 || missing > n) {
        return false;
    }
    if (missing == 0) {
        progress = true;
        return true;
    }

    auto entry = [&](slong row_index, slong column_index) -> ulong {
        return basis[static_cast<std::size_t>(
                column_index * n + row_index)];
    };

    std::vector<ulong> reduced = row;
    slong pivot = n - 1;
    while (pivot >= 0 && reduced[static_cast<std::size_t>(pivot)] == 0) {
        --pivot;
    }

    while (pivot >= 0) {
        const ulong diag = entry(pivot, pivot);
        if (diag != 0) {
            const ulong scale = reduced[static_cast<std::size_t>(pivot)];
            for (slong i = 0; i < pivot; ++i) {
                const ulong basis_value = entry(i, pivot);
                if (basis_value != 0) {
                    reduced[static_cast<std::size_t>(i)] =
                            (reduced[static_cast<std::size_t>(i)] +
                             scale * (kRelationAdmissionModulus - basis_value)) %
                            kRelationAdmissionModulus;
                }
            }
            reduced[static_cast<std::size_t>(pivot)] = 0;
            do {
                --pivot;
            } while (pivot >= 0 &&
                     reduced[static_cast<std::size_t>(pivot)] == 0);
            continue;
        }

        progress = true;
        return true;
    }

    return true;
}

bool reset_modular_basis(
        detail::RelationAdmissionCache& cache,
        slong n) noexcept {
    std::size_t size = 0;
    if (!relation_basis_size(size, n)) {
        return false;
    }
    cache.modular_basis.assign(size, 0);
    cache.missing = n;
    return true;
}

ulong modular_inverse(ulong value) noexcept {
    slong t = 0;
    slong next_t = 1;
    slong r = static_cast<slong>(kRelationAdmissionModulus);
    slong next_r = static_cast<slong>(value % kRelationAdmissionModulus);

    while (next_r != 0) {
        const slong quotient = r / next_r;
        const slong old_t = t;
        t = next_t;
        next_t = old_t - quotient * next_t;
        const slong old_r = r;
        r = next_r;
        next_r = old_r - quotient * next_r;
    }
    if (r != 1) {
        return 0;
    }
    if (t < 0) {
        t += static_cast<slong>(kRelationAdmissionModulus);
    }
    return static_cast<ulong>(t);
}

bool relation_mod_row(std::vector<ulong>& out,
                      slong& first_nonzero,
                      flint::FmpzMatConstRef row) noexcept {
    const slong length = flint::fmpz_mat_ncols(row);
    first_nonzero = length;
    if (flint::fmpz_mat_nrows(row) != 1 || length < 0) {
        return false;
    }

    out.assign(static_cast<std::size_t>(length), 0);
    for (slong i = 0; i < length; ++i) {
        flint::FmpzConstRef entry =
                flint::fmpz_mat_entry(flint::FmpzMatConstRef(row), 0, i);
        if (!flint::fmpz_is_zero(entry) && first_nonzero == length) {
            first_nonzero = i;
        }
        out[static_cast<std::size_t>(i)] =
                flint::fmpz_fdiv_ui(entry, kRelationAdmissionModulus);
    }
    return true;
}

bool insert_modular_basis_row(bool& progress,
                               std::vector<ulong>& basis,
                               slong& missing,
                               const std::vector<ulong>& row) noexcept {
    progress = false;
    if (row.size() >
        static_cast<std::size_t>(std::numeric_limits<slong>::max())) {
        return false;
    }
    const slong n = static_cast<slong>(row.size());
    std::size_t basis_size = 0;
    if (!relation_basis_size(basis_size, n) ||
        basis.size() != basis_size || missing < 0 || missing > n) {
        return false;
    }
    if (missing == 0) {
        progress = true;
        return true;
    }

    auto entry = [&](slong row_index, slong column_index) -> ulong& {
        return basis[static_cast<std::size_t>(
                column_index * n + row_index)];
    };

    std::vector<ulong> reduced = row;
    slong pivot = n - 1;
    while (pivot >= 0 && reduced[static_cast<std::size_t>(pivot)] == 0) {
        --pivot;
    }

    while (pivot >= 0) {
        ulong& diag = entry(pivot, pivot);
        if (diag != 0) {
            const ulong scale = reduced[static_cast<std::size_t>(pivot)];
            for (slong i = 0; i < pivot; ++i) {
                const ulong basis_value = entry(i, pivot);
                if (basis_value != 0) {
                    reduced[static_cast<std::size_t>(i)] =
                            (reduced[static_cast<std::size_t>(i)] +
                             scale * (kRelationAdmissionModulus - basis_value)) %
                            kRelationAdmissionModulus;
                }
            }
            reduced[static_cast<std::size_t>(pivot)] = 0;
            do {
                --pivot;
            } while (pivot >= 0 &&
                     reduced[static_cast<std::size_t>(pivot)] == 0);
            continue;
        }

        const ulong inverse = modular_inverse(
                reduced[static_cast<std::size_t>(pivot)]);
        if (inverse == 0) {
            return false;
        }

        for (slong i = pivot; i-- > 0;) {
            const ulong value = reduced[static_cast<std::size_t>(i)];
            if (value == 0 || entry(i, i) == 0) {
                continue;
            }
            const ulong scale = kRelationAdmissionModulus - value;
            for (slong j = 0; j < i; ++j) {
                const ulong basis_value = entry(j, i);
                if (basis_value != 0) {
                    reduced[static_cast<std::size_t>(j)] =
                            (reduced[static_cast<std::size_t>(j)] +
                             scale * basis_value) %
                            kRelationAdmissionModulus;
                }
            }
            reduced[static_cast<std::size_t>(i)] = 0;
        }

        for (slong i = 0; i < pivot; ++i) {
            const ulong value = reduced[static_cast<std::size_t>(i)];
            entry(i, pivot) = value == 0
                    ? UWORD(0)
                    : (value * inverse) % kRelationAdmissionModulus;
        }
        entry(pivot, pivot) = 1;

        for (slong i = pivot + 1; i < n; ++i) {
            const ulong value = entry(pivot, i);
            if (value == 0) {
                continue;
            }
            const ulong scale = kRelationAdmissionModulus - value;
            for (slong j = 0; j < pivot; ++j) {
                const ulong inserted = entry(j, pivot);
                if (inserted != 0) {
                    entry(j, i) =
                            (entry(j, i) + scale * inserted) %
                            kRelationAdmissionModulus;
                }
            }
            entry(pivot, i) = 0;
        }

        --missing;
        progress = true;
        return true;
    }

    return true;
}

}  // namespace

bool begin_selected_pivot_recollection(
        detail::RelationAdmissionCache& cache,
        const std::vector<slong>& indices,
        slong n) noexcept {
    if (!relation_basis_is_sized(cache, n) || cache.missing != 0) {
        return false;
    }

    slong missing = 0;
    for (slong index : indices) {
        if (index < 0 || index >= n) {
            return false;
        }
        ulong& diag = relation_basis_entry(cache, index, index, n);
        if (diag != 0) {
            diag = 0;
            ++missing;
        }
    }
    cache.missing = missing;
    return true;
}

void end_selected_pivot_recollection(detail::RelationAdmissionCache& cache,
                                     const std::vector<slong>& indices,
                                     slong n) noexcept {
    if (!relation_basis_is_sized(cache, n)) {
        return;
    }
    for (slong index : indices) {
        if (index >= 0 && index < n) {
            relation_basis_entry(cache, index, index, n) = 1;
        }
    }
    cache.missing = 0;
}

}  // namespace detail::relation_search

using namespace detail::relation_search;

namespace detail {

bool try_admit_relation(
        ClassGroupContext& context,
        bool& retained,
        RelationAdmissionCache& cache,
        const Relation& relation,
        bool in_random_relation,
        bool factor_base_verified) noexcept {
    retained = false;
    if (!relation.is_set() || !context.has_factor_base()) {
        return false;
    }
    const FactorBase* relation_base = relation.factor_base();
    const FactorBase* context_base = context.factor_base();
    if (relation_base == nullptr || context_base == nullptr ||
        (factor_base_verified
                 ? relation.length() != context_base->length()
                 : !relation_base->equal(*context_base))) {
        return false;
    }

    const flint::FmpzMatConstRef relation_row = relation.exponents_ref();
    std::vector<ulong> mod_row;
    slong first_nonzero = relation.length();
    ulong relation_row_hash = 0;
    {
        SILEX_PROFILE_SCOPE(context.diagnostics(),
                            DiagnosticsModule::class_group,
                            "class_group.relation_admission.mod_row");
        if (!relation_mod_row(mod_row, first_nonzero,
                              relation_row)) {
            return false;
        }
        relation_row_hash = relation_mod_row_hash(mod_row);
    }
    {
        SILEX_PROFILE_SCOPE(context.diagnostics(),
                            DiagnosticsModule::class_group,
                            "class_group.relation_admission.duplicate_check");
        bool already_known = false;
        if (!relation_row_already_known_with_cache(
                    already_known, cache,
                    ClassGroupRelationSearchAccess::relations(context), relation_row,
                    first_nonzero, relation_row_hash,
                    context.diagnostics())) {
            return false;
        }
        if (already_known) {
            ++cache.duplicate_relation_count;
            return true;
        }
    }
    if (!relation_basis_is_sized(cache, relation.length())) {
        return false;
    }

    bool modular_progress = false;
    {
        SILEX_PROFILE_SCOPE(context.diagnostics(),
                            DiagnosticsModule::class_group,
                            "class_group.relation_admission.modular_classify");
        if (!modular_basis_row_has_progress(
                    modular_progress, cache.modular_basis, cache.missing,
                    mod_row)) {
            return false;
        }
    }
    std::vector<ulong> modular_basis_candidate;
    slong missing_candidate = cache.missing;
    bool modular_basis_candidate_valid = false;
    if (modular_progress && cache.missing > 0) {
        {
            SILEX_PROFILE_SCOPE(context.diagnostics(),
                                DiagnosticsModule::class_group,
                                "class_group.relation_admission.modular_basis_copy");
            modular_basis_candidate = cache.modular_basis;
        }
        {
            SILEX_PROFILE_SCOPE(context.diagnostics(),
                                DiagnosticsModule::class_group,
                                "class_group.relation_admission.modular_insert");
            if (!insert_modular_basis_row(modular_progress,
                                           modular_basis_candidate,
                                           missing_candidate, mod_row)) {
                return false;
            }
        }
        modular_basis_candidate_valid = modular_progress;
    }

    const bool skip_dependent_before_full_rank =
            cache.missing > 0 && !modular_progress &&
        cache.relation_surplus <= 0 &&
        !in_random_relation;
    if (skip_dependent_before_full_rank) {
        ++cache.skipped_dependent_relation_count;
        ClassGroupRelationSearchAccess::record_skipped_dependent(context);
        SILEX_PROFILE_EVENT(
                context.diagnostics(), DiagnosticsModule::class_group,
                "class_group.relation_admission.modular_skip");
        log_relation_admission_trace(context.diagnostics(), "skipped",
                                    first_nonzero, modular_progress, cache,
                                    context.relation_count(),
                                    in_random_relation, relation);
        return true;
    }

    const slong before = context.relation_count();
    bool skipped_dependent = false;
    if (!ClassGroupRelationSearchAccess::append_search_relation(
                skipped_dependent, context, relation)) {
        return false;
    }
    if (skipped_dependent) {
        ++cache.skipped_dependent_relation_count;
        log_relation_admission_trace(context.diagnostics(), "skipped",
                                    first_nonzero, modular_progress, cache,
                                    context.relation_count(),
                                    in_random_relation, relation);
        return true;
    }
    if (context.relation_count() <= before) {
        return true;
    }

    retained = true;
    if (cache.relation_row_hashes.size() ==
        static_cast<std::size_t>(before)) {
        cache.relation_row_hashes.push_back(relation_row_hash);
    } else if (!sync_relation_row_hashes(
                       cache, ClassGroupRelationSearchAccess::relations(context))) {
        return false;
    }
    if (modular_progress) {
        if (modular_basis_candidate_valid) {
            cache.modular_basis = std::move(modular_basis_candidate);
            cache.missing = missing_candidate;
        }
    }
    ++cache.retained_relation_count;
    if (!modular_progress && cache.missing > 0 &&
        cache.relation_surplus > 0 && first_nonzero < relation.length()) {
        --cache.relation_surplus;
    }
    log_relation_admission_trace(context.diagnostics(), "retained",
                                first_nonzero, modular_progress, cache,
                                context.relation_count(),
                                in_random_relation, relation);
    return true;
}

bool try_admit_deferred_integral_relation(
        ClassGroupContext& context,
        bool& retained,
        RelationAdmissionCache& cache,
        Relation& relation,
        flint::FmpzMatConstRef integral_coordinates,
        bool in_random_relation,
        bool factor_base_verified) noexcept {
    retained = false;
    if (!relation.is_defined() || !context.has_factor_base()) {
        return false;
    }
    const FactorBase* relation_base = relation.factor_base();
    const FactorBase* context_base = context.factor_base();
    if (relation_base == nullptr || context_base == nullptr ||
        (factor_base_verified
                 ? relation.length() != context_base->length()
                 : !relation_base->equal(*context_base))) {
        return false;
    }

    const flint::FmpzMatConstRef relation_row =
            detail::pending_relation_exponents_ref(relation);
    if (flint::fmpz_mat_nrows(relation_row) != 1 ||
        flint::fmpz_mat_ncols(relation_row) != relation.length()) {
        return false;
    }

    std::vector<ulong> mod_row;
    slong first_nonzero = relation.length();
    ulong relation_row_hash = 0;
    {
        SILEX_PROFILE_SCOPE(context.diagnostics(),
                            DiagnosticsModule::class_group,
                            "class_group.relation_admission.mod_row");
        if (!relation_mod_row(mod_row, first_nonzero, relation_row)) {
            return false;
        }
        relation_row_hash = relation_mod_row_hash(mod_row);
    }
    {
        SILEX_PROFILE_SCOPE(context.diagnostics(),
                            DiagnosticsModule::class_group,
                            "class_group.relation_admission.duplicate_check");
        bool already_known = false;
        if (!relation_row_already_known_with_cache(
                    already_known, cache,
                    ClassGroupRelationSearchAccess::relations(context), relation_row,
                    first_nonzero, relation_row_hash,
                    context.diagnostics())) {
            return false;
        }
        if (already_known) {
            ++cache.duplicate_relation_count;
            return true;
        }
    }
    if (!relation_basis_is_sized(cache, relation.length())) {
        return false;
    }

    bool modular_progress = false;
    {
        SILEX_PROFILE_SCOPE(context.diagnostics(),
                            DiagnosticsModule::class_group,
                            "class_group.relation_admission.modular_classify");
        if (!modular_basis_row_has_progress(
                    modular_progress, cache.modular_basis, cache.missing,
                    mod_row)) {
            return false;
        }
    }
    std::vector<ulong> modular_basis_candidate;
    slong missing_candidate = cache.missing;
    bool modular_basis_candidate_valid = false;
    if (modular_progress && cache.missing > 0) {
        {
            SILEX_PROFILE_SCOPE(context.diagnostics(),
                                DiagnosticsModule::class_group,
                                "class_group.relation_admission.modular_basis_copy");
            modular_basis_candidate = cache.modular_basis;
        }
        {
            SILEX_PROFILE_SCOPE(context.diagnostics(),
                                DiagnosticsModule::class_group,
                                "class_group.relation_admission.modular_insert");
            if (!insert_modular_basis_row(modular_progress,
                                           modular_basis_candidate,
                                           missing_candidate, mod_row)) {
                return false;
            }
        }
        modular_basis_candidate_valid = modular_progress;
    }

    const bool skip_dependent_before_full_rank =
            cache.missing > 0 && !modular_progress &&
            cache.relation_surplus <= 0 && !in_random_relation;
    if (skip_dependent_before_full_rank) {
        ++cache.skipped_dependent_relation_count;
        ClassGroupRelationSearchAccess::record_skipped_dependent(context);
        SILEX_PROFILE_EVENT(
                context.diagnostics(), DiagnosticsModule::class_group,
                "class_group.relation_admission.modular_skip");
        log_relation_admission_trace(context.diagnostics(), "skipped",
                                    first_nonzero, modular_progress, cache,
                                    context.relation_count(),
                                    in_random_relation, relation_row);
        return true;
    }

    if (!detail::commit_relation_generator_from_integral_coordinates(
                relation, integral_coordinates, context.diagnostics())) {
        return false;
    }

    const slong before = context.relation_count();
    bool skipped_dependent = false;
    if (!ClassGroupRelationSearchAccess::append_search_relation(
                skipped_dependent, context, relation)) {
        return false;
    }
    if (skipped_dependent) {
        ++cache.skipped_dependent_relation_count;
        log_relation_admission_trace(context.diagnostics(), "skipped",
                                    first_nonzero, modular_progress, cache,
                                    context.relation_count(),
                                    in_random_relation, relation);
        return true;
    }
    if (context.relation_count() <= before) {
        return true;
    }

    retained = true;
    if (cache.relation_row_hashes.size() ==
        static_cast<std::size_t>(before)) {
        cache.relation_row_hashes.push_back(relation_row_hash);
    } else if (!sync_relation_row_hashes(
                       cache, ClassGroupRelationSearchAccess::relations(context))) {
        return false;
    }
    if (modular_progress) {
        if (modular_basis_candidate_valid) {
            cache.modular_basis = std::move(modular_basis_candidate);
            cache.missing = missing_candidate;
        }
    }
    ++cache.retained_relation_count;
    if (!modular_progress && cache.missing > 0 &&
        cache.relation_surplus > 0 && first_nonzero < relation.length()) {
        --cache.relation_surplus;
    }
    log_relation_admission_trace(context.diagnostics(), "retained",
                                first_nonzero, modular_progress, cache,
                                context.relation_count(),
                                in_random_relation, relation);
    return true;
}

}  // namespace detail


bool detail::initialize_relation_admission_cache(
        ClassGroupContext& context,
        detail::RelationAdmissionCache& cache,
        const Order& order,
        slong add_need) noexcept {
    SILEX_PROFILE_SCOPE(context.diagnostics(), DiagnosticsModule::class_group,
                        "class_group.relation_admission.initialize");
    const FactorBase* base = context.factor_base();
    if (base == nullptr || !base->is_defined() || add_need < 0 ||
        add_need > WORD_MAX - base->length()) {
        return false;
    }

    cache = detail::RelationAdmissionCache{};
    cache.target_relation_count =
            context.relation_count() + base->length() + add_need;
    cache.relation_surplus = add_need;
    cache.relation_count_before_init = context.relation_count();
    if (!reset_modular_basis(cache, base->length())) {
        return false;
    }

    flint::Fmpz rational_prime;
    for (slong i = 0; i < base->rational_prime_block_count(); ++i) {
        slong start = 0;
        slong length = 0;
        bool complete = false;
        if (!base->rational_prime_block(flint::FmpzRef(rational_prime),
                                        start, length, i) ||
            !detail::FactorBaseBlockAccess::rational_prime_block_is_complete(
                    complete, *base, i)) {
            return false;
        }
        if (!complete) {
            continue;
        }

        // A rational-prime relation is valid for initialization only when the
        // factor base contains the complete block above that prime. Exact
        // relation validation remains authoritative.
        Relation relation(*base);
        bool retained = false;
        if (!set_complete_prime_block_relation(
                    relation, *base, order,
                    flint::FmpzConstRef(rational_prime),
                    context.diagnostics()) ||
            !detail::try_admit_relation(
                    context, retained, cache, relation, false)) {
            return false;
        }
        if (retained) {
            ++cache.trivial_relation_count;
        }
    }

    return true;
}

}  // namespace silex
