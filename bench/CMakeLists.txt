find_package(benchmark CONFIG QUIET)

if(NOT benchmark_FOUND)
    message(FATAL_ERROR
        "SILEX_BUILD_BENCHMARKS=ON requires Google Benchmark. Install its "
        "CMake package or set CMAKE_PREFIX_PATH/benchmark_DIR to that package.")
endif()

add_executable(b-silex-status
    b-status.cpp)

target_compile_features(b-silex-status PRIVATE cxx_std_20)
target_link_libraries(b-silex-status
    PRIVATE
        Silex::silex
        benchmark::benchmark)
target_include_directories(b-silex-status PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(b-silex-status PRIVATE
    SILEX_BENCHMARK_EXPECTED_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    SILEX_BENCHMARK_EXPECTED_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    SILEX_BENCHMARK_EXPECTED_VERSION_PATCH=${PROJECT_VERSION_PATCH}
    SILEX_BENCHMARK_EXPECTED_VERSION_STRING="${PROJECT_VERSION}")
silex_apply_target_warnings(b-silex-status)

set(SILEX_NATIVE_BENCHMARK_TARGETS b-silex-status)

function(silex_add_benchmark target source)
    add_executable(${target} ${source})
    target_compile_features(${target} PRIVATE cxx_std_20)
    target_link_libraries(${target}
        PRIVATE
            Silex::silex
            benchmark::benchmark)
    target_include_directories(${target} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})
    silex_apply_target_warnings(${target})
    list(APPEND SILEX_NATIVE_BENCHMARK_TARGETS ${target})
    set(SILEX_NATIVE_BENCHMARK_TARGETS
        "${SILEX_NATIVE_BENCHMARK_TARGETS}" PARENT_SCOPE)
endfunction()

function(silex_benchmark_release_rows target)
    if(NOT TARGET ${target})
        message(FATAL_ERROR
            "Cannot classify rows for unknown benchmark target ${target}")
    endif()
    if(ARGC LESS 2)
        message(FATAL_ERROR
            "Benchmark target ${target} needs at least one release row")
    endif()
    set_property(TARGET ${target} PROPERTY
        SILEX_RELEASE_BENCHMARK_ROWS "${ARGN}")
endfunction()

function(silex_benchmark_excluded_rows target reason)
    if(NOT TARGET ${target})
        message(FATAL_ERROR
            "Cannot classify rows for unknown benchmark target ${target}")
    endif()
    if(ARGC LESS 3 OR reason STREQUAL "")
        message(FATAL_ERROR
            "Excluded benchmark rows need a nonempty reason")
    endif()
    foreach(row IN LISTS ARGN)
        set_property(TARGET ${target} APPEND PROPERTY
            SILEX_EXCLUDED_BENCHMARK_ROWS "${row}")
        set_property(TARGET ${target} APPEND PROPERTY
            SILEX_EXCLUDED_BENCHMARK_REASONS "${reason}")
    endforeach()
endfunction()

silex_add_benchmark(b-silex-lat_hnf lat/b-lat_hnf.cpp)
silex_add_benchmark(b-silex-lat_lll lat/b-lat_lll.cpp)
target_include_directories(b-silex-lat_lll PRIVATE
        ${PROJECT_SOURCE_DIR}/src)
silex_add_benchmark(b-silex-lat_intersection lat/b-lat_intersection.cpp)
silex_add_benchmark(b-silex-fmpz_smat_kernels
        fmpz_smat/b-fmpz_smat_kernels.cpp)
silex_add_benchmark(b-silex-fmpz_smat_rank
        fmpz_smat/b-fmpz_smat_rank.cpp)
silex_add_benchmark(b-silex-fmpz_smat_hnf_ctx
        fmpz_smat/b-fmpz_smat_hnf_ctx.cpp)
silex_add_benchmark(b-silex-nf_ord_maximal_order
        nf_ord/b-order_maximal_order.cpp)
silex_add_benchmark(b-silex-nf_ord_pmaximal_overorder
        nf_ord/b-order_pmaximal_overorder.cpp)
silex_add_benchmark(b-silex-nf_clgp
        nf_clgp/b-class_group_pipeline.cpp)
if(SILEX_BUILD_TESTS)
    add_test(NAME silex-class-unit-google-benchmark-success
        COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/test/test_class_unit_google_benchmark.py
                --exe $<TARGET_FILE:b-silex-nf_clgp>)
    set_tests_properties(silex-class-unit-google-benchmark-success
        PROPERTIES
            LABELS release
            RESOURCE_LOCK silex_native_benchmark_semantic
            TIMEOUT 300)
endif()
silex_add_benchmark(b-silex-nf_clgp_factor_base_honesty
        nf_clgp/b-factor_base_honesty_witness.cpp)
target_include_directories(b-silex-nf_clgp_factor_base_honesty
        PRIVATE ${PROJECT_SOURCE_DIR}/src)
silex_add_benchmark(b-silex-nf_fac_elt
        nf_fac_elt/b-roots_compact.cpp)
silex_add_benchmark(b-silex-nf_fac_elt_compact_reconstruction
        nf_fac_elt/b-bounded_compact_reconstruction.cpp)
target_include_directories(b-silex-nf_fac_elt_compact_reconstruction
        PRIVATE ${PROJECT_SOURCE_DIR}/src)
silex_add_benchmark(b-silex-nf_idl
        nf_idl/b-factor_over_base.cpp)
silex_add_benchmark(b-silex-nf_prime_idl
        nf_prime_idl/b-local.cpp)
silex_add_benchmark(b-silex-zeta_bf_linear_factor_count
        zeta/b-bf_linear_factor_count.cpp)
silex_add_benchmark(b-silex-zeta_bf_prime_scratch
        zeta/b-bf_prime_scratch.cpp)

# Release semantic coverage is intentionally explicit.  The role audit below
# compares this classification with each executable's --benchmark_list_tests
# output, so a newly registered row cannot silently miss the release gate.
silex_benchmark_release_rows(b-silex-status
    BM_status_ok
    BM_version_string)
silex_benchmark_release_rows(b-silex-lat_hnf
    BM_lat_hnf)
silex_benchmark_release_rows(b-silex-lat_lll
    BM_lat_lll
    BM_lat_lll_hnf_degree14
    BM_native_ideal_hnf_degree14_flint)
silex_benchmark_excluded_rows(b-silex-lat_lll
    "requires the optional fplll backend"
    BM_native_ideal_hnf_degree14_fplll
    BM_native_ideal_hnf_degree14_fplll_bkz/block_size:4
    BM_native_ideal_hnf_degree14_fplll_bkz/block_size:8
    BM_native_ideal_hnf_degree14_fplll_bkz/block_size:10
    BM_native_ideal_hnf_degree14_fplll_bkz/block_size:12
    BM_native_ideal_hnf_degree14_fplll_bkz/block_size:14)
silex_benchmark_excluded_rows(b-silex-lat_lll
    "requires the optional flatter backend"
    BM_native_ideal_hnf_degree14_flatter_rhf_1_02)
silex_benchmark_release_rows(b-silex-lat_intersection
    BM_lat_intersection)
silex_benchmark_release_rows(b-silex-fmpz_smat_kernels
    BM_fmpz_smat_mul_fmpz_mat
    BM_fmpz_smat_transpose)
silex_benchmark_release_rows(b-silex-fmpz_smat_rank
    BM_fmpz_smat_rank_mod_prime_ui
    BM_fmpz_smat_mod_rank_add_rows)
silex_benchmark_release_rows(b-silex-fmpz_smat_hnf_ctx
    BM_fmpz_smat_hnf_ctx_add_rows)
silex_benchmark_release_rows(b-silex-nf_ord_maximal_order
    BM_order_maximal_order_quadratic_fast
    BM_order_maximal_order_generic_quadratic
    BM_order_maximal_order_stable_cubic)
silex_benchmark_release_rows(b-silex-nf_ord_pmaximal_overorder
    BM_order_pmaximal_overorder_quadratic_fast
    BM_order_pmaximal_overorder_generic_quadratic)
silex_benchmark_release_rows(b-silex-nf_clgp
    BM_class_group_trivial_cubic_candidate
    BM_class_group_nontrivial_cubic_candidate
    BM_class_group_nontrivial_cubic_candidate_fb2x
    BM_class_group_nontrivial_cubic_candidate_fb4x
    BM_class_group_mixed_cubic_candidate
    BM_class_group_cubic_3_1_candidate
    BM_class_group_cubic_minus3_minus1_candidate
    BM_class_group_small_quartic_candidate
    BM_class_group_small_quartic_candidate_fb2x
    BM_class_group_small_quartic_candidate_fb4x
    BM_class_group_quartic_1_0_minus2_1_candidate
    BM_class_group_quartic_minus1_minus1_0_2_candidate
    BM_class_group_quartic_minus1_minus1_0_2_candidate_fb2x
    BM_class_group_quartic_minus1_minus1_0_2_candidate_fb4x
    BM_class_group_hard_cubic_0_minus4_minus7_candidate
    BM_class_group_hard_quartic_minus1_minus1_0_2_candidate
    BM_class_group_hard_quartic_minus1_minus1_0_2_fb2x_candidate
    BM_class_group_quintic_minus1_minus1_candidate
    BM_class_group_sextic_minus1_minus1_candidate
    BM_class_unit_0_1_0_degree_one_proven
    BM_class_unit_0_1_0_real_quadratic_proven
    BM_class_unit_0_1_0_imag_quadratic_proven
    BM_class_unit_0_1_0_cubic_trivial_proven
    BM_class_unit_0_1_0_cubic_nontrivial_proven
    BM_class_unit_0_1_0_quartic_cyclotomic_proven
    BM_class_unit_0_1_0_quartic_noncyclotomic_proven
    BM_class_unit_0_1_0_quartic_disc70640_proven
    BM_class_unit_0_1_0_quartic_disc223479_proven
    BM_class_unit_0_1_0_quartic_disc35019_proven
    BM_class_unit_0_1_0_quartic_disc1412343_proven
    BM_class_unit_0_1_0_quartic_x4_minus_x_minus_1_proven
    BM_class_unit_0_1_0_quintic_proven
    BM_class_unit_0_1_0_quintic_disc11119_proven
    BM_class_unit_0_1_0_quintic_disc401370255_proven
    BM_class_unit_0_1_0_quintic_disc57895_proven
    BM_class_unit_0_1_0_sextic_proven
    BM_class_unit_random_matrix/3/4/0/1/0/iterations:1)
silex_benchmark_excluded_rows(b-silex-nf_clgp
    "bounded diagnostic rows intentionally include unsuccessful outcomes"
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap/100/100
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap/250/250
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap/500/500
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap/750/750
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap/1000/1000
    BM_diagnostic_class_unit_0_1_0_quintic_disc57895_proven_cap/1500/1500)
silex_benchmark_excluded_rows(b-silex-nf_clgp
    "random diagnostic sweep rows characterize unsupported inputs"
    BM_class_unit_random_sweep/2/4/0/1/0
    BM_class_unit_random_sweep/2/8/0/1/0
    BM_class_unit_random_sweep/3/2/0/1/0)
silex_benchmark_excluded_rows(b-silex-nf_clgp
    "random matrix diagnostics are not tagged-release promises"
    BM_class_unit_random_matrix/2/2/0/0/0/iterations:1
    BM_class_unit_random_matrix/2/4/0/0/0/iterations:1
    BM_class_unit_random_matrix/2/8/0/0/0/iterations:1
    BM_class_unit_random_matrix/3/2/0/0/0/iterations:1
    BM_class_unit_random_matrix/3/4/0/0/0/iterations:1
    BM_class_unit_random_matrix/3/8/0/0/0/iterations:1
    BM_class_unit_random_matrix/4/2/0/0/0/iterations:1
    BM_class_unit_random_matrix/4/4/0/0/0/iterations:1
    BM_class_unit_random_matrix/5/2/0/0/0/iterations:1
    BM_class_unit_random_matrix/6/2/2/0/0/iterations:1
    BM_class_unit_random_matrix/2/4/0/1/0/iterations:1
    BM_class_unit_random_matrix/2/8/0/1/0/iterations:1
    BM_class_unit_random_matrix/3/2/0/1/0/iterations:1
    BM_class_unit_random_matrix/4/2/0/1/0/iterations:1
    BM_class_unit_random_matrix/5/2/0/1/0/iterations:1)
silex_benchmark_excluded_rows(b-silex-nf_clgp
    "candidate-cap boundary diagnostics permit expected failures"
    BM_class_unit_random_boundary_clgp/250
    BM_class_unit_random_boundary_clgp/1000
    BM_class_unit_random_boundary_clgp/5000)
silex_benchmark_release_rows(b-silex-nf_clgp_factor_base_honesty
    BM_honesty_scalar_full_factor_quartic
    BM_honesty_scalar_required_prime_quartic
    BM_honesty_scalar_full_factor_quintic
    BM_honesty_scalar_required_prime_quintic
    BM_honesty_current_materializing_quartic
    BM_honesty_order_element_direct_quartic
    BM_honesty_current_materializing_quintic
    BM_honesty_order_element_direct_quintic)
silex_benchmark_release_rows(b-silex-nf_fac_elt
    BM_nf_elt_is_square_quadratic
    BM_nf_elt_is_power_quadratic_square
    BM_nf_fac_elt_evaluate_cubic
    BM_nf_fac_elt_is_power_reduced_cubic
    BM_nf_fac_elt_is_power_compact_cubic
    BM_nf_compact_elt_set_fac_elt_cubic
    BM_nf_compact_elt_root_base_cubic
    BM_nf_fac_elt_log_embed_cubic)
silex_benchmark_release_rows(b-silex-nf_fac_elt_compact_reconstruction
    BM_bound_basis_leading_one_quadratic/manual_time
    BM_bound_basis_reordered_quadratic/manual_time
    BM_bound_basis_index_three_quartic/manual_time
    BM_direct_quartic_disc1412343/iterations:1/manual_time
    BM_crt_quartic_disc1412343/iterations:1/manual_time
    BM_bound_quartic_disc1412343/iterations:1/manual_time)
silex_benchmark_excluded_rows(
    b-silex-nf_fac_elt_compact_reconstruction
    "unbounded strict class/unit fixture generation is not a routine release-semantic gate"
    BM_direct_random_d15_h2_s43/iterations:1/manual_time
    BM_crt_random_d15_h2_s43/iterations:1/manual_time
    BM_bound_random_d15_h2_s43/iterations:1/manual_time
    BM_direct_random_d9_h4_s8/iterations:1/manual_time
    BM_crt_random_d9_h4_s8/iterations:1/manual_time
    BM_bound_random_d9_h4_s8/iterations:1/manual_time)
silex_benchmark_release_rows(b-silex-nf_idl
    BM_ideal_factor_over_base_quadratic)
silex_benchmark_release_rows(b-silex-nf_prime_idl
    BM_prime_ideal_reduce_nf_elem_quadratic
    BM_prime_ideal_valuation_quadratic
    BM_prime_ideal_valuation_order_element_quadratic_split)
silex_benchmark_release_rows(b-silex-zeta_bf_linear_factor_count
    BM_bf_linear_factor_count/cubic_blocker_current
    BM_bf_linear_factor_count/cubic_blocker_array_preinverse
    BM_bf_linear_factor_count/cubic_blocker_array_binexp
    BM_bf_linear_factor_count/cubic_blocker_roots
    BM_bf_linear_factor_count/cubic_blocker_distinct_degree
    BM_bf_linear_factor_count/quartic_blocker_current
    BM_bf_linear_factor_count/quartic_blocker_array_preinverse
    BM_bf_linear_factor_count/quartic_blocker_array_binexp
    BM_bf_linear_factor_count/quartic_blocker_roots
    BM_bf_linear_factor_count/quartic_blocker_distinct_degree
    BM_bf_linear_factor_count/cubic_control_current
    BM_bf_linear_factor_count/cubic_control_array_preinverse
    BM_bf_linear_factor_count/cubic_control_array_binexp
    BM_bf_linear_factor_count/cubic_control_roots
    BM_bf_linear_factor_count/cubic_control_distinct_degree
    BM_bf_linear_factor_count/quartic_control_current
    BM_bf_linear_factor_count/quartic_control_array_preinverse
    BM_bf_linear_factor_count/quartic_control_array_binexp
    BM_bf_linear_factor_count/quartic_control_roots
    BM_bf_linear_factor_count/quartic_control_distinct_degree
    BM_bf_linear_factor_count/quintic_control_current
    BM_bf_linear_factor_count/quintic_control_array_preinverse
    BM_bf_linear_factor_count/quintic_control_array_binexp
    BM_bf_linear_factor_count/quintic_control_roots
    BM_bf_linear_factor_count/quintic_control_distinct_degree
    BM_bf_linear_factor_count/sextic_control_current
    BM_bf_linear_factor_count/sextic_control_array_preinverse
    BM_bf_linear_factor_count/sextic_control_array_binexp
    BM_bf_linear_factor_count/sextic_control_roots
    BM_bf_linear_factor_count/sextic_control_distinct_degree)
silex_benchmark_release_rows(b-silex-zeta_bf_prime_scratch
    BM_bf_prime_scratch/cubic_blocker_per_prime
    BM_bf_prime_scratch/cubic_blocker_per_term
    BM_bf_prime_scratch/quartic_blocker_per_prime
    BM_bf_prime_scratch/quartic_blocker_per_term
    BM_bf_prime_scratch/cubic_control_per_prime
    BM_bf_prime_scratch/cubic_control_per_term
    BM_bf_prime_scratch/quartic_control_per_prime
    BM_bf_prime_scratch/quartic_control_per_term
    BM_bf_prime_scratch/quintic_control_per_prime
    BM_bf_prime_scratch/quintic_control_per_term
    BM_bf_prime_scratch/sextic_control_per_prime
    BM_bf_prime_scratch/sextic_control_per_term)

get_filename_component(SILEX_BENCHMARK_PACKAGE_PREFIX
    "${benchmark_DIR}/../../.." ABSOLUTE)
find_file(SILEX_GOOGLE_BENCHMARK_COMPARE_PY
    NAMES compare.py
    HINTS ${SILEX_BENCHMARK_PACKAGE_PREFIX}
    PATH_SUFFIXES
        share/googlebenchmark/tools
        share/benchmark/tools)
if(NOT SILEX_GOOGLE_BENCHMARK_COMPARE_PY)
    set(SILEX_GOOGLE_BENCHMARK_COMPARE_PY "")
endif()
mark_as_advanced(SILEX_GOOGLE_BENCHMARK_COMPARE_PY)

function(silex_benchmark_json_escape output value)
    string(REPLACE "\\" "\\\\" value "${value}")
    string(REPLACE "\"" "\\\"" value "${value}")
    foreach(code RANGE 1 31)
        string(ASCII ${code} control_character)
        string(HEX "${control_character}" hexadecimal_code)
        string(REPLACE "${control_character}"
            "\\u00${hexadecimal_code}" value "${value}")
    endforeach()
    set(${output} "${value}" PARENT_SCOPE)
endfunction()

function(silex_benchmark_require_direct_json_safe description value)
    foreach(code RANGE 1 31)
        string(ASCII ${code} control_character)
        string(FIND "${value}" "${control_character}" control_index)
        if(NOT control_index EQUAL -1)
            string(HEX "${control_character}" hexadecimal_code)
            message(FATAL_ERROR
                "${description} contains JSON control character 0x${hexadecimal_code}. "
                "Benchmark executable paths cannot contain JSON control "
                "characters, quotes, or backslashes.")
        endif()
    endforeach()
    foreach(forbidden_character IN ITEMS "\"" "\\")
        string(FIND "${value}" "${forbidden_character}" forbidden_index)
        if(NOT forbidden_index EQUAL -1)
            message(FATAL_ERROR
                "${description} contains a quote or backslash. Benchmark "
                "executable paths cannot contain JSON control characters, "
                "quotes, or backslashes.")
        endif()
    endforeach()
endfunction()

function(silex_benchmark_json_string_array output)
    set(content "")
    foreach(value IN LISTS ARGN)
        silex_benchmark_json_escape(escaped_value "${value}")
        string(APPEND content "\"${escaped_value}\",")
    endforeach()
    string(REGEX REPLACE ",$" "" content "${content}")
    set(${output} "[${content}]" PARENT_SCOPE)
endfunction()

function(silex_benchmark_file_digests output)
    set(digests "")
    foreach(candidate IN LISTS ARGN)
        if(IS_ABSOLUTE "${candidate}" AND
           EXISTS "${candidate}" AND
           NOT IS_DIRECTORY "${candidate}")
            file(SHA256 "${candidate}" digest)
            list(APPEND digests "${candidate}=${digest}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES digests)
    string(JOIN ";" digests ${digests})
    set(${output} "${digests}" PARENT_SCOPE)
endfunction()

set(SILEX_BENCHMARK_FPLLL_LIBRARIES ${FPLLL_LINK_LIBRARIES})
set(SILEX_BENCHMARK_FPLLL_INCLUDE_DIRS ${FPLLL_INCLUDE_DIRS})
if(SILEX_FPLLL_LIBRARY)
    list(APPEND SILEX_BENCHMARK_FPLLL_LIBRARIES
        "${SILEX_FPLLL_LIBRARY}")
endif()
if(SILEX_FPLLL_INCLUDE_DIR)
    list(APPEND SILEX_BENCHMARK_FPLLL_INCLUDE_DIRS
        "${SILEX_FPLLL_INCLUDE_DIR}")
endif()
list(REMOVE_DUPLICATES SILEX_BENCHMARK_FPLLL_LIBRARIES)
list(REMOVE_DUPLICATES SILEX_BENCHMARK_FPLLL_INCLUDE_DIRS)
set(SILEX_BENCHMARK_FLINT_LIBRARIES ${FLINT_LINK_LIBRARIES})
set(SILEX_BENCHMARK_FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIRS})
set(SILEX_BENCHMARK_FLINT_LIBRARY_DIRS ${FLINT_LIBRARY_DIRS})
set(SILEX_BENCHMARK_FLATTER_LIBRARIES "${SILEX_FLATTER_LIBRARY}")
set(SILEX_BENCHMARK_FLATTER_INCLUDE_DIRS "${SILEX_FLATTER_INCLUDE_DIR}")
silex_benchmark_file_digests(SILEX_BENCHMARK_FLINT_LIBRARY_DIGESTS
    ${FLINT_LINK_LIBRARIES})
silex_benchmark_file_digests(SILEX_BENCHMARK_FPLLL_LIBRARY_DIGESTS
    ${FPLLL_LINK_LIBRARIES} ${SILEX_FPLLL_LIBRARY})
silex_benchmark_file_digests(SILEX_BENCHMARK_FLATTER_LIBRARY_DIGESTS
    ${SILEX_FLATTER_LIBRARY})

set(SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARIES "")
set(SILEX_BENCHMARK_GOOGLE_BENCHMARK_LOCATION_PROPERTIES
        IMPORTED_LOCATION_RELEASE
        IMPORTED_LOCATION_RELWITHDEBINFO
        IMPORTED_LOCATION_NOCONFIG
        IMPORTED_LOCATION
        IMPORTED_IMPLIB_RELEASE
        IMPORTED_IMPLIB_RELWITHDEBINFO
        IMPORTED_IMPLIB_NOCONFIG
        IMPORTED_IMPLIB)
get_target_property(
    SILEX_BENCHMARK_GOOGLE_BENCHMARK_IMPORTED_CONFIGURATIONS
    benchmark::benchmark IMPORTED_CONFIGURATIONS)
if(SILEX_BENCHMARK_GOOGLE_BENCHMARK_IMPORTED_CONFIGURATIONS
        AND NOT SILEX_BENCHMARK_GOOGLE_BENCHMARK_IMPORTED_CONFIGURATIONS
                MATCHES "-NOTFOUND$")
    foreach(imported_configuration IN LISTS
            SILEX_BENCHMARK_GOOGLE_BENCHMARK_IMPORTED_CONFIGURATIONS)
        string(TOUPPER "${imported_configuration}"
            imported_configuration_upper)
        list(APPEND SILEX_BENCHMARK_GOOGLE_BENCHMARK_LOCATION_PROPERTIES
            "IMPORTED_LOCATION_${imported_configuration_upper}"
            "IMPORTED_IMPLIB_${imported_configuration_upper}")
    endforeach()
endif()
list(REMOVE_DUPLICATES
    SILEX_BENCHMARK_GOOGLE_BENCHMARK_LOCATION_PROPERTIES)
foreach(location_property IN LISTS
        SILEX_BENCHMARK_GOOGLE_BENCHMARK_LOCATION_PROPERTIES)
    get_target_property(location benchmark::benchmark ${location_property})
    if(location AND NOT location MATCHES "-NOTFOUND$")
        list(APPEND SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARIES "${location}")
    endif()
endforeach()
list(REMOVE_DUPLICATES SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARIES)
silex_benchmark_file_digests(
    SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARY_DIGESTS
    ${SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARIES})

foreach(variable IN ITEMS
        PROJECT_SOURCE_DIR
        PROJECT_BINARY_DIR
        CMAKE_GENERATOR
        CMAKE_GENERATOR_PLATFORM
        CMAKE_GENERATOR_TOOLSET
        CMAKE_CXX_COMPILER
        CMAKE_CXX_COMPILER_ID
        CMAKE_CXX_COMPILER_VERSION
        CMAKE_CXX_COMPILER_LAUNCHER
        CMAKE_CXX_COMPILER_TARGET
        CMAKE_CXX_FLAGS
        CMAKE_CXX_FLAGS_DEBUG
        CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_RELWITHDEBINFO
        CMAKE_CXX_FLAGS_MINSIZEREL
        CMAKE_EXE_LINKER_FLAGS
        CMAKE_EXE_LINKER_FLAGS_DEBUG
        CMAKE_EXE_LINKER_FLAGS_RELEASE
        CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
        CMAKE_EXE_LINKER_FLAGS_MINSIZEREL
        CMAKE_INTERPROCEDURAL_OPTIMIZATION
        CMAKE_INTERPROCEDURAL_OPTIMIZATION_DEBUG
        CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE
        CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO
        CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL
        CMAKE_POSITION_INDEPENDENT_CODE
        CMAKE_UNITY_BUILD
        CMAKE_UNITY_BUILD_BATCH_SIZE
        BUILD_SHARED_LIBS
        CMAKE_SYSROOT
        CMAKE_TOOLCHAIN_FILE
        SILEX_BENCHMARK_FLINT_LIBRARIES
        SILEX_BENCHMARK_FLINT_INCLUDE_DIRS
        SILEX_BENCHMARK_FLINT_LIBRARY_DIRS
        SILEX_BENCHMARK_FLINT_LIBRARY_DIGESTS
        SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARIES
        SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARY_DIGESTS
        SILEX_BENCHMARK_FPLLL_LIBRARIES
        SILEX_BENCHMARK_FPLLL_INCLUDE_DIRS
        SILEX_BENCHMARK_FPLLL_LIBRARY_DIGESTS
        SILEX_BENCHMARK_FLATTER_LIBRARIES
        SILEX_BENCHMARK_FLATTER_INCLUDE_DIRS
        SILEX_BENCHMARK_FLATTER_LIBRARY_DIGESTS
        SILEX_FPLLL_ROOT
        SILEX_FLATTER_ROOT
        SILEX_GOOGLE_BENCHMARK_COMPARE_PY
        SILEX_WITH_FPLLL
        SILEX_WITH_FLATTER
        SILEX_ENABLE_LOGGING
        SILEX_ENABLE_DEBUG_CHECKS
        SILEX_ENABLE_PROFILING
        SILEX_ENABLE_SANITIZERS
        SILEX_ENABLE_FRAME_POINTERS
        SILEX_ENABLE_EXCEPTIONS
        SILEX_ENABLE_RTTI
        FLINT_VERSION
        FPLLL_VERSION
        benchmark_VERSION)
    silex_benchmark_json_escape(escaped_value "${${variable}}")
    set("SILEX_BENCHMARK_JSON_${variable}" "${escaped_value}")
endforeach()

set(SILEX_NATIVE_BENCHMARK_MANIFEST_TARGETS "")
set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CALLS "")
# file(GENERATE resolves TARGET_FILE after configure-time JSON escaping is
# available.  Keep every target-path contributor safe for direct JSON output,
# and include the same inputs when selecting the generated CMake bracket
# delimiter below.
set(SILEX_NATIVE_BENCHMARK_PATH_INPUTS
    "${PROJECT_BINARY_DIR}"
    "${CMAKE_BUILD_TYPE}"
    "${CMAKE_CONFIGURATION_TYPES}"
    "${CMAKE_EXECUTABLE_SUFFIX}")
set(SILEX_NATIVE_BENCHMARK_CONFIGURATIONS
    ${CMAKE_CONFIGURATION_TYPES}
    ${CMAKE_BUILD_TYPE}
    Debug
    Release
    RelWithDebInfo
    MinSizeRel)
list(REMOVE_DUPLICATES SILEX_NATIVE_BENCHMARK_CONFIGURATIONS)
foreach(SILEX_NATIVE_BENCHMARK_TARGET IN LISTS
        SILEX_NATIVE_BENCHMARK_TARGETS)
    set(path_properties
        RUNTIME_OUTPUT_DIRECTORY
        RUNTIME_OUTPUT_NAME
        OUTPUT_NAME
        PREFIX
        SUFFIX)
    foreach(configuration IN LISTS SILEX_NATIVE_BENCHMARK_CONFIGURATIONS)
        string(TOUPPER "${configuration}" configuration_upper)
        list(APPEND path_properties
            "RUNTIME_OUTPUT_DIRECTORY_${configuration_upper}"
            "RUNTIME_OUTPUT_NAME_${configuration_upper}"
            "OUTPUT_NAME_${configuration_upper}"
            "${configuration_upper}_POSTFIX")
    endforeach()
    list(REMOVE_DUPLICATES path_properties)
    foreach(path_property IN LISTS path_properties)
        get_target_property(path_property_value
            ${SILEX_NATIVE_BENCHMARK_TARGET} ${path_property})
        if(path_property_value AND
           NOT path_property_value MATCHES "-NOTFOUND$")
            list(APPEND SILEX_NATIVE_BENCHMARK_PATH_INPUTS
                "${path_property_value}")
            silex_benchmark_require_direct_json_safe(
                "${SILEX_NATIVE_BENCHMARK_TARGET} ${path_property}"
                "${path_property_value}")
        endif()
    endforeach()
endforeach()
foreach(path_input IN LISTS SILEX_NATIVE_BENCHMARK_PATH_INPUTS)
    silex_benchmark_require_direct_json_safe(
        "Benchmark executable path input" "${path_input}")
endforeach()

set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_EQUALS "==")
while(TRUE)
    set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE
        "]${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_EQUALS}]")
    set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_DELIMITER_FOUND FALSE)
    foreach(path_input IN LISTS SILEX_NATIVE_BENCHMARK_PATH_INPUTS)
        string(FIND "${path_input}"
            "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE}"
            SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE_INDEX)
        if(NOT SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE_INDEX EQUAL -1)
            set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_DELIMITER_FOUND TRUE)
            break()
        endif()
    endforeach()
    if(NOT SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_DELIMITER_FOUND)
        break()
    endif()
    string(APPEND SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_EQUALS "=")
endwhile()
set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_OPEN
    "[${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_EQUALS}[")
foreach(SILEX_NATIVE_BENCHMARK_TARGET IN LISTS
        SILEX_NATIVE_BENCHMARK_TARGETS)
    get_target_property(release_rows ${SILEX_NATIVE_BENCHMARK_TARGET}
        SILEX_RELEASE_BENCHMARK_ROWS)
    if(NOT release_rows OR release_rows MATCHES "-NOTFOUND$")
        message(FATAL_ERROR
            "Benchmark target ${SILEX_NATIVE_BENCHMARK_TARGET} has no "
            "classified release rows")
    endif()
    get_target_property(excluded_rows ${SILEX_NATIVE_BENCHMARK_TARGET}
        SILEX_EXCLUDED_BENCHMARK_ROWS)
    get_target_property(excluded_reasons ${SILEX_NATIVE_BENCHMARK_TARGET}
        SILEX_EXCLUDED_BENCHMARK_REASONS)
    if(NOT excluded_rows OR excluded_rows MATCHES "-NOTFOUND$")
        set(excluded_rows "")
        set(excluded_reasons "")
    endif()

    set(classified_rows ${release_rows} ${excluded_rows})
    set(unique_classified_rows ${classified_rows})
    list(REMOVE_DUPLICATES unique_classified_rows)
    list(LENGTH classified_rows classified_row_count)
    list(LENGTH unique_classified_rows unique_classified_row_count)
    if(NOT classified_row_count EQUAL unique_classified_row_count)
        message(FATAL_ERROR
            "Benchmark target ${SILEX_NATIVE_BENCHMARK_TARGET} has a row "
            "classified more than once")
    endif()
    foreach(row IN LISTS classified_rows)
        if(NOT row MATCHES "^[A-Za-z0-9_:/-]+$")
            message(FATAL_ERROR
                "Benchmark row ${row} contains a regex metacharacter; add "
                "escaping before using it in a release semantic filter")
        endif()
    endforeach()

    list(LENGTH excluded_rows excluded_row_count)
    list(LENGTH excluded_reasons excluded_reason_count)
    if(NOT excluded_row_count EQUAL excluded_reason_count)
        message(FATAL_ERROR
            "Benchmark target ${SILEX_NATIVE_BENCHMARK_TARGET} has "
            "incomplete exclusion reasons")
    endif()

    silex_benchmark_json_string_array(release_rows_json ${release_rows})
    set(excluded_rows_json "")
    if(excluded_row_count GREATER 0)
        math(EXPR excluded_last "${excluded_row_count} - 1")
        foreach(index RANGE 0 ${excluded_last})
            list(GET excluded_rows ${index} excluded_row)
            list(GET excluded_reasons ${index} excluded_reason)
            silex_benchmark_json_escape(
                excluded_row_json "${excluded_row}")
            silex_benchmark_json_escape(
                excluded_reason_json "${excluded_reason}")
            string(APPEND excluded_rows_json
                "{\"name\":\"${excluded_row_json}\","
                "\"reason\":\"${excluded_reason_json}\"},")
        endforeach()
        string(REGEX REPLACE ",$" "" excluded_rows_json
            "${excluded_rows_json}")
    endif()
    string(APPEND SILEX_NATIVE_BENCHMARK_MANIFEST_TARGETS
        "    {\"target\":\"${SILEX_NATIVE_BENCHMARK_TARGET}\","
        "\"path\":\"$<TARGET_FILE:${SILEX_NATIVE_BENCHMARK_TARGET}>\","
        "\"release_semantic_rows\":${release_rows_json},"
        "\"excluded_semantic_rows\":[${excluded_rows_json}]},\n")
    string(APPEND SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CALLS
        "silex_audit_benchmark_roles("
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_OPEN}"
        "${SILEX_NATIVE_BENCHMARK_TARGET}"
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE} "
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_OPEN}"
        "$<TARGET_FILE:${SILEX_NATIVE_BENCHMARK_TARGET}>"
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE} "
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_OPEN}${release_rows}"
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE} "
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_OPEN}${excluded_rows}"
        "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CLOSE})\n")
endforeach()
string(REGEX REPLACE ",\n$" "\n"
    SILEX_NATIVE_BENCHMARK_MANIFEST_TARGETS
    "${SILEX_NATIVE_BENCHMARK_MANIFEST_TARGETS}")
set(SILEX_NATIVE_BENCHMARK_MANIFEST
    "${CMAKE_CURRENT_BINARY_DIR}/silex-native-benchmarks-$<CONFIG>.json")
file(GENERATE
    OUTPUT "${SILEX_NATIVE_BENCHMARK_MANIFEST}"
    CONTENT "{
  \"schema_version\": 1,
  \"kind\": \"silex-native-benchmark-executables\",
  \"source_dir\": \"${SILEX_BENCHMARK_JSON_PROJECT_SOURCE_DIR}\",
  \"build_dir\": \"${SILEX_BENCHMARK_JSON_PROJECT_BINARY_DIR}\",
  \"build\": {
    \"configuration\": \"$<CONFIG>\",
    \"cmake_generator\": \"${SILEX_BENCHMARK_JSON_CMAKE_GENERATOR}\",
    \"cmake_generator_platform\": \"${SILEX_BENCHMARK_JSON_CMAKE_GENERATOR_PLATFORM}\",
    \"cmake_generator_toolset\": \"${SILEX_BENCHMARK_JSON_CMAKE_GENERATOR_TOOLSET}\",
    \"cxx_compiler\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_COMPILER}\",
    \"cxx_compiler_id\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_COMPILER_ID}\",
    \"cxx_compiler_version\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_COMPILER_VERSION}\",
    \"cxx_compiler_launcher\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_COMPILER_LAUNCHER}\",
    \"cxx_compiler_target\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_COMPILER_TARGET}\",
    \"cxx_flags\": {
      \"all\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_FLAGS}\",
      \"debug\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_FLAGS_DEBUG}\",
      \"release\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_FLAGS_RELEASE}\",
      \"relwithdebinfo\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_FLAGS_RELWITHDEBINFO}\",
      \"minsizerel\": \"${SILEX_BENCHMARK_JSON_CMAKE_CXX_FLAGS_MINSIZEREL}\"
    },
    \"executable_linker_flags\": {
      \"all\": \"${SILEX_BENCHMARK_JSON_CMAKE_EXE_LINKER_FLAGS}\",
      \"debug\": \"${SILEX_BENCHMARK_JSON_CMAKE_EXE_LINKER_FLAGS_DEBUG}\",
      \"release\": \"${SILEX_BENCHMARK_JSON_CMAKE_EXE_LINKER_FLAGS_RELEASE}\",
      \"relwithdebinfo\": \"${SILEX_BENCHMARK_JSON_CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}\",
      \"minsizerel\": \"${SILEX_BENCHMARK_JSON_CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}\"
    },
    \"interprocedural_optimization\": {
      \"all\": \"${SILEX_BENCHMARK_JSON_CMAKE_INTERPROCEDURAL_OPTIMIZATION}\",
      \"debug\": \"${SILEX_BENCHMARK_JSON_CMAKE_INTERPROCEDURAL_OPTIMIZATION_DEBUG}\",
      \"release\": \"${SILEX_BENCHMARK_JSON_CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE}\",
      \"relwithdebinfo\": \"${SILEX_BENCHMARK_JSON_CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO}\",
      \"minsizerel\": \"${SILEX_BENCHMARK_JSON_CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL}\"
    },
    \"sysroot\": \"${SILEX_BENCHMARK_JSON_CMAKE_SYSROOT}\",
    \"toolchain_file\": \"${SILEX_BENCHMARK_JSON_CMAKE_TOOLCHAIN_FILE}\",
    \"flint_version\": \"${SILEX_BENCHMARK_JSON_FLINT_VERSION}\",
    \"google_benchmark_compare\": \"${SILEX_BENCHMARK_JSON_SILEX_GOOGLE_BENCHMARK_COMPARE_PY}\",
    \"dependencies\": {
      \"flint\": {
        \"version\": \"${SILEX_BENCHMARK_JSON_FLINT_VERSION}\",
        \"include_dirs\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLINT_INCLUDE_DIRS}\",
        \"library_dirs\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLINT_LIBRARY_DIRS}\",
        \"libraries\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLINT_LIBRARIES}\",
        \"library_digests\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLINT_LIBRARY_DIGESTS}\"
      },
      \"google_benchmark\": {
        \"version\": \"${SILEX_BENCHMARK_JSON_benchmark_VERSION}\",
        \"libraries\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARIES}\",
        \"library_digests\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_GOOGLE_BENCHMARK_LIBRARY_DIGESTS}\"
      }
    },
    \"optional_dependencies\": {
      \"fplll\": {
        \"enabled\": \"${SILEX_BENCHMARK_JSON_SILEX_WITH_FPLLL}\",
        \"version\": \"${SILEX_BENCHMARK_JSON_FPLLL_VERSION}\",
        \"root\": \"${SILEX_BENCHMARK_JSON_SILEX_FPLLL_ROOT}\",
        \"include_dirs\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FPLLL_INCLUDE_DIRS}\",
        \"libraries\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FPLLL_LIBRARIES}\",
        \"library_digests\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FPLLL_LIBRARY_DIGESTS}\"
      },
      \"flatter\": {
        \"enabled\": \"${SILEX_BENCHMARK_JSON_SILEX_WITH_FLATTER}\",
        \"root\": \"${SILEX_BENCHMARK_JSON_SILEX_FLATTER_ROOT}\",
        \"include_dirs\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLATTER_INCLUDE_DIRS}\",
        \"libraries\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLATTER_LIBRARIES}\",
        \"library_digests\": \"${SILEX_BENCHMARK_JSON_SILEX_BENCHMARK_FLATTER_LIBRARY_DIGESTS}\"
      }
    },
    \"options\": {
      \"logging\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_LOGGING}\",
      \"debug_checks\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_DEBUG_CHECKS}\",
      \"profiling\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_PROFILING}\",
      \"sanitizers\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_SANITIZERS}\",
      \"frame_pointers\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_FRAME_POINTERS}\",
      \"exceptions\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_EXCEPTIONS}\",
      \"rtti\": \"${SILEX_BENCHMARK_JSON_SILEX_ENABLE_RTTI}\",
      \"shared_libraries\": \"${SILEX_BENCHMARK_JSON_BUILD_SHARED_LIBS}\",
      \"position_independent_code\": \"${SILEX_BENCHMARK_JSON_CMAKE_POSITION_INDEPENDENT_CODE}\",
      \"unity_build\": \"${SILEX_BENCHMARK_JSON_CMAKE_UNITY_BUILD}\",
      \"unity_build_batch_size\": \"${SILEX_BENCHMARK_JSON_CMAKE_UNITY_BUILD_BATCH_SIZE}\",
      \"fplll\": \"${SILEX_BENCHMARK_JSON_SILEX_WITH_FPLLL}\",
      \"flatter\": \"${SILEX_BENCHMARK_JSON_SILEX_WITH_FLATTER}\"
    }
  },
  \"targets\": [
${SILEX_NATIVE_BENCHMARK_MANIFEST_TARGETS}  ]
}
")

set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT
    "${CMAKE_CURRENT_BINARY_DIR}/silex-native-benchmark-role-audit-$<CONFIG>.cmake")
set(SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CONTENT [=[
function(silex_audit_benchmark_roles target executable release_rows excluded_rows)
    execute_process(
        COMMAND "${executable}" --benchmark_list_tests
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR
            "${target}: --benchmark_list_tests failed (${result}): ${error}")
    endif()
    string(REPLACE "\r\n" "\n" output "${output}")
    string(REGEX REPLACE "\n+$" "" output "${output}")
    if(output STREQUAL "")
        set(actual_rows "")
    else()
        string(REPLACE "\n" ";" actual_rows "${output}")
    endif()

    set(classified_rows ${release_rows} ${excluded_rows})
    set(unique_actual_rows ${actual_rows})
    set(unique_classified_rows ${classified_rows})
    list(REMOVE_DUPLICATES unique_actual_rows)
    list(REMOVE_DUPLICATES unique_classified_rows)
    list(LENGTH actual_rows actual_count)
    list(LENGTH unique_actual_rows unique_actual_count)
    list(LENGTH classified_rows classified_count)
    list(LENGTH unique_classified_rows unique_classified_count)
    if(NOT actual_count EQUAL unique_actual_count)
        message(FATAL_ERROR
            "${target}: --benchmark_list_tests contains duplicate rows")
    endif()
    if(NOT classified_count EQUAL unique_classified_count)
        message(FATAL_ERROR
            "${target}: release/excluded row roles overlap or contain duplicates")
    endif()
    list(SORT actual_rows)
    list(SORT classified_rows)
    if(NOT actual_rows STREQUAL classified_rows)
        message(FATAL_ERROR
            "${target}: registered rows do not equal the explicit release/excluded roles\n"
            "registered: ${actual_rows}\nclassified: ${classified_rows}")
    endif()
endfunction()

]=])
string(APPEND SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CONTENT
    "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CALLS}")
file(GENERATE
    OUTPUT "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT}"
    CONTENT "${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT_CONTENT}")

function(silex_add_release_benchmark_semantic_gate target)
    get_target_property(release_rows ${target}
        SILEX_RELEASE_BENCHMARK_ROWS)
    string(JOIN "|" release_filter_body ${release_rows})
    set(gate_timeout 120)
    set(test_timeout 150)
    if(target STREQUAL "b-silex-nf_clgp")
        set(gate_timeout 600)
        set(test_timeout 660)
    endif()
    set(gate_arguments
        gate
        --manifest ${SILEX_NATIVE_BENCHMARK_MANIFEST}
        --target ${target}
        --filter "^(${release_filter_body})$"
        --require-counter success=1
        --require-counter failure_reason=0
        --timeout ${gate_timeout}
        --output
        ${CMAKE_CURRENT_BINARY_DIR}/semantic-gate/${target}.json
        --overwrite)
    foreach(row IN LISTS release_rows)
        list(APPEND gate_arguments --expect-row "${row}")
    endforeach()
    set(test_name "silex-${target}-release-semantic-gate")
    add_test(NAME ${test_name}
        COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tools/bench/native-benchmark.py
                ${gate_arguments})
    set_tests_properties(${test_name}
        PROPERTIES
            LABELS "release;benchmark-semantic"
            RESOURCE_LOCK silex_native_benchmark_semantic
            TIMEOUT ${test_timeout})
endfunction()

if(SILEX_BUILD_TESTS)
    add_test(NAME silex-native-benchmark-manifest
        COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/test/test_native_benchmark.py
                --tool
                ${PROJECT_SOURCE_DIR}/tools/bench/native-benchmark.py
                --manifest ${SILEX_NATIVE_BENCHMARK_MANIFEST})
    set_tests_properties(silex-native-benchmark-manifest
        PROPERTIES
            LABELS release
            TIMEOUT 60)

    add_test(NAME silex-native-benchmark-role-audit
        COMMAND ${CMAKE_COMMAND} -P
                ${SILEX_NATIVE_BENCHMARK_ROLE_AUDIT})
    set_tests_properties(silex-native-benchmark-role-audit
        PROPERTIES
            LABELS "release;benchmark-semantic"
            TIMEOUT 60)

    foreach(target IN LISTS SILEX_NATIVE_BENCHMARK_TARGETS)
        silex_add_release_benchmark_semantic_gate(${target})
    endforeach()

    add_test(NAME silex-compact-bound-contract-benchmark-semantic-gate
        COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tools/bench/native-benchmark.py
                gate
                --manifest ${SILEX_NATIVE_BENCHMARK_MANIFEST}
                --target b-silex-nf_fac_elt_compact_reconstruction
                --filter
                "^(BM_bound_basis_(leading_one_quadratic|reordered_quadratic|index_three_quartic)/manual_time|BM_bound_quartic_disc1412343/iterations:1/manual_time)$"
                --expect-row BM_bound_basis_leading_one_quadratic/manual_time
                --expect-row BM_bound_basis_reordered_quadratic/manual_time
                --expect-row BM_bound_basis_index_three_quartic/manual_time
                --expect-row
                BM_bound_quartic_disc1412343/iterations:1/manual_time
                --require-counter success=1
                --require-counter failure_reason=0
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:degree=2
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:unit_rank=1
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:equation_index=2
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:basis_leading_one=1
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:normalization_required=0
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:coordinate_bound_exact_match=1
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:bit_bound_exact_match=1
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:bound_attempts_exact_match=1
                --require-row-counter
                BM_bound_basis_leading_one_quadratic/manual_time:bound_precision_exact_match=1
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:degree=2
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:unit_rank=1
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:equation_index=2
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:basis_leading_one=0
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:normalization_required=1
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:coordinate_bound_exact_match=1
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:bit_bound_exact_match=1
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:bound_attempts_exact_match=1
                --require-row-counter
                BM_bound_basis_reordered_quadratic/manual_time:bound_precision_exact_match=1
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:degree=4
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:unit_rank=2
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:equation_index=3
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:basis_leading_one=0
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:normalization_required=1
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:coordinate_bound_exact_match=1
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:bit_bound_exact_match=1
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:bound_attempts_exact_match=1
                --require-row-counter
                BM_bound_basis_index_three_quartic/manual_time:bound_precision_exact_match=1
                --require-row-counter
                BM_bound_quartic_disc1412343/iterations:1/manual_time:coordinate_bound_exact_match=1
                --require-row-counter
                BM_bound_quartic_disc1412343/iterations:1/manual_time:bit_bound_exact_match=1
                --require-row-counter
                BM_bound_quartic_disc1412343/iterations:1/manual_time:bound_attempts_exact_match=1
                --require-row-counter
                BM_bound_quartic_disc1412343/iterations:1/manual_time:bound_precision_exact_match=1
                --timeout 45
                --output
                ${CMAKE_CURRENT_BINARY_DIR}/semantic-gate/compact-bound-contract.json
                --overwrite)
    set_tests_properties(silex-compact-bound-contract-benchmark-semantic-gate
        PROPERTIES
            LABELS "release;benchmark-semantic"
            RESOURCE_LOCK silex_native_benchmark_semantic
            TIMEOUT 60)
endif()

add_custom_target(benchmarks
    DEPENDS ${SILEX_NATIVE_BENCHMARK_TARGETS})

function(silex_audit_benchmark_directory directory)
    get_property(directory_targets DIRECTORY "${directory}"
        PROPERTY BUILDSYSTEM_TARGETS)
    foreach(directory_target IN LISTS directory_targets)
        get_target_property(directory_target_type ${directory_target} TYPE)
        if(directory_target_type STREQUAL "EXECUTABLE")
            list(FIND SILEX_NATIVE_BENCHMARK_TARGETS ${directory_target}
                managed_target_index)
            if(managed_target_index EQUAL -1)
                message(FATAL_ERROR
                    "Benchmark executable ${directory_target} is absent from "
                    "the managed native inventory and has no release-row audit")
            endif()
        endif()
    endforeach()
    get_property(subdirectories DIRECTORY "${directory}"
        PROPERTY SUBDIRECTORIES)
    foreach(subdirectory IN LISTS subdirectories)
        silex_audit_benchmark_directory("${subdirectory}")
    endforeach()
endfunction()

silex_audit_benchmark_directory("${CMAKE_CURRENT_SOURCE_DIR}")
