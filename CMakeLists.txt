cmake_minimum_required(VERSION 3.20)
project(Silex
    VERSION 0.1.0
    DESCRIPTION "Native C++20 library for computational algebraic number theory"
    HOMEPAGE_URL "https://github.com/wjyoumans/silex"
    LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(SILEX_BUILD_TESTS "Build Silex tests" ON)
option(SILEX_BUILD_RELEASE_TESTS
    "Build Silex package and source-archive release tests" OFF)
option(SILEX_BUILD_EXAMPLES "Build Silex examples" ON)
option(SILEX_BUILD_BENCHMARKS "Build Silex Google Benchmark targets" OFF)
option(SILEX_BUILD_BENCHMARK_ADAPTERS
    "Build Silex source-tree benchmark adapter executables" OFF)
option(SILEX_BUILD_DOCS "Build Silex Sphinx docs" OFF)
option(SILEX_ENABLE_LOGGING "Enable logging facade" ON)
option(SILEX_ENABLE_DEBUG_CHECKS "Enable debug-check facade" OFF)
option(SILEX_ENABLE_PROFILING "Enable profiling facade" OFF)
option(SILEX_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer for supported compilers" OFF)
option(SILEX_ENABLE_FRAME_POINTERS "Keep frame pointers for external profiling tools" OFF)
option(SILEX_ENABLE_EXCEPTIONS "Enable C++ exceptions in the Silex core target" OFF)
option(SILEX_ENABLE_RTTI "Enable C++ RTTI in the Silex core target" OFF)
option(SILEX_ENABLE_WERROR "Treat compiler warnings as errors for Silex-owned targets" OFF)
option(SILEX_WITH_FPLLL "Use optional fplll backend" OFF)
option(SILEX_WITH_FLATTER "Use optional flatter backend" OFF)
set(SILEX_DOCS_RELEASE "${PROJECT_VERSION}" CACHE STRING
    "Release string rendered in Silex documentation")
set(SILEX_DOCS_CHANNEL "dev" CACHE STRING
    "Documentation channel (dev or stable)")
set_property(CACHE SILEX_DOCS_CHANNEL PROPERTY STRINGS dev stable)
foreach(SILEX_UNIMPLEMENTED_BACKEND_OPTION IN ITEMS
        SILEX_WITH_SPDLOG SILEX_WITH_TRACY)
    if(DEFINED ${SILEX_UNIMPLEMENTED_BACKEND_OPTION} AND
       ${SILEX_UNIMPLEMENTED_BACKEND_OPTION})
        message(FATAL_ERROR
            "${SILEX_UNIMPLEMENTED_BACKEND_OPTION} is not implemented; "
            "use the Silex diagnostics facade until a concrete backend exists")
    endif()
endforeach()
set(SILEX_FPLLL_ROOT "$ENV{SILEX_FPLLL_ROOT}" CACHE PATH
    "Optional fplll install prefix used when SILEX_WITH_FPLLL=ON")
set(SILEX_FPLLL_INCLUDE_DIR "" CACHE PATH
    "Optional explicit fplll include directory used when SILEX_WITH_FPLLL=ON")
set(SILEX_FPLLL_LIBRARY "" CACHE FILEPATH
    "Optional explicit fplll library used when SILEX_WITH_FPLLL=ON")
set(SILEX_FLATTER_ROOT "$ENV{SILEX_FLATTER_ROOT}" CACHE PATH
    "Optional flatter install prefix used when SILEX_WITH_FLATTER=ON")
set(SILEX_FLATTER_INCLUDE_DIR "" CACHE PATH
    "Optional explicit flatter include directory used when SILEX_WITH_FLATTER=ON")
set(SILEX_FLATTER_LIBRARY "" CACHE FILEPATH
    "Optional explicit flatter library used when SILEX_WITH_FLATTER=ON")
set(SILEX_REQUIRED_FLINT_VERSION "3.0.0" CACHE STRING
    "Minimum FLINT version required by Silex")

find_package(PkgConfig REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint>=${SILEX_REQUIRED_FLINT_VERSION})
set(SILEX_HAVE_FPLLL 0)
set(SILEX_HAVE_FLATTER 0)
set(SILEX_FPLLL_TARGET "")
set(SILEX_FLATTER_TARGET "")
if(SILEX_WITH_FPLLL)
    if(SILEX_FPLLL_ROOT)
        list(PREPEND CMAKE_PREFIX_PATH "${SILEX_FPLLL_ROOT}")
    endif()
    pkg_check_modules(FPLLL QUIET IMPORTED_TARGET fplll)
    if(TARGET PkgConfig::FPLLL)
        add_library(SilexExternal::fplll INTERFACE IMPORTED)
        set_target_properties(SilexExternal::fplll PROPERTIES
            INTERFACE_LINK_LIBRARIES PkgConfig::FPLLL)
        set(SILEX_FPLLL_TARGET SilexExternal::fplll)
    else()
        if(NOT SILEX_FPLLL_INCLUDE_DIR)
            unset(SILEX_FPLLL_INCLUDE_DIR CACHE)
        endif()
        if(NOT SILEX_FPLLL_LIBRARY)
            unset(SILEX_FPLLL_LIBRARY CACHE)
        endif()
        find_path(SILEX_FPLLL_INCLUDE_DIR
            NAMES fplll/wrapper.h
            HINTS ${SILEX_FPLLL_ROOT} ENV SILEX_FPLLL_ROOT
            PATH_SUFFIXES include)
        find_library(SILEX_FPLLL_LIBRARY
            NAMES fplll
            HINTS ${SILEX_FPLLL_ROOT} ENV SILEX_FPLLL_ROOT
            PATH_SUFFIXES lib lib64)
        if(NOT SILEX_FPLLL_INCLUDE_DIR OR NOT SILEX_FPLLL_LIBRARY)
            message(FATAL_ERROR
                "SILEX_WITH_FPLLL=ON requires fplll. Set SILEX_FPLLL_ROOT, "
                "SILEX_FPLLL_INCLUDE_DIR, SILEX_FPLLL_LIBRARY, "
                "CMAKE_PREFIX_PATH, or PKG_CONFIG_PATH.")
        endif()
        add_library(SilexExternal::fplll UNKNOWN IMPORTED)
        set_target_properties(SilexExternal::fplll PROPERTIES
            IMPORTED_LOCATION "${SILEX_FPLLL_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${SILEX_FPLLL_INCLUDE_DIR}"
            INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${SILEX_FPLLL_INCLUDE_DIR}")
        set(SILEX_FPLLL_TARGET SilexExternal::fplll)
    endif()
    set(SILEX_HAVE_FPLLL 1)
endif()
if(SILEX_WITH_FLATTER)
    if(SILEX_FLATTER_ROOT)
        list(PREPEND CMAKE_PREFIX_PATH "${SILEX_FLATTER_ROOT}")
    endif()
    if(NOT SILEX_FLATTER_INCLUDE_DIR)
        unset(SILEX_FLATTER_INCLUDE_DIR CACHE)
    endif()
    if(NOT SILEX_FLATTER_LIBRARY)
        unset(SILEX_FLATTER_LIBRARY CACHE)
    endif()
    find_path(SILEX_FLATTER_INCLUDE_DIR
        NAMES flatter/problems.h
        HINTS
            ${SILEX_FLATTER_ROOT}
            ENV SILEX_FLATTER_ROOT
            "$ENV{HOME}/.local"
        PATH_SUFFIXES include)
    find_library(SILEX_FLATTER_LIBRARY
        NAMES flatter
        HINTS
            ${SILEX_FLATTER_ROOT}
            ENV SILEX_FLATTER_ROOT
            "$ENV{HOME}/.local"
        PATH_SUFFIXES lib lib64 build/lib)
    if(NOT SILEX_FLATTER_INCLUDE_DIR OR
       NOT EXISTS "${SILEX_FLATTER_INCLUDE_DIR}/flatter/problems.h" OR
       NOT EXISTS "${SILEX_FLATTER_INCLUDE_DIR}/flatter/data/lattice.h" OR
       NOT SILEX_FLATTER_LIBRARY)
        message(FATAL_ERROR
            "SILEX_WITH_FLATTER=ON requires flatter headers and library. Set "
            "SILEX_FLATTER_ROOT, SILEX_FLATTER_INCLUDE_DIR, "
            "SILEX_FLATTER_LIBRARY, or CMAKE_PREFIX_PATH.")
    endif()
    add_library(SilexExternal::flatter UNKNOWN IMPORTED)
    set_target_properties(SilexExternal::flatter PROPERTIES
        IMPORTED_LOCATION "${SILEX_FLATTER_LIBRARY}")
    set(SILEX_FLATTER_TARGET SilexExternal::flatter)
    set(SILEX_HAVE_FLATTER 1)
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/silex/build_config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/silex/build_config.hpp
    @ONLY)

if(SILEX_ENABLE_SANITIZERS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address,undefined)
    else()
        message(WARNING "SILEX_ENABLE_SANITIZERS=ON is only configured for Clang and GNU compilers")
    endif()
endif()

if(SILEX_ENABLE_FRAME_POINTERS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        add_compile_options(-fno-omit-frame-pointer)
    else()
        message(WARNING "SILEX_ENABLE_FRAME_POINTERS=ON is only configured for Clang and GNU compilers")
    endif()
endif()

function(silex_apply_target_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(SILEX_ENABLE_WERROR)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(SILEX_ENABLE_WERROR)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

add_library(silex
    src/abelian_group/abelian_group.cpp
    src/diagnostics/diagnostics.cpp
    src/fmpz_smat/fmpz_smat.cpp
    src/lat/lat.cpp
    src/lat/fplll_backend.cpp
    src/lat/flatter_backend.cpp
    src/class_group/class_relation_hnf_witness.cpp
    src/class_group/class_relation_module_context.cpp
    src/class_group/class_group.cpp
    src/class_group/factor_base_honesty.cpp
    src/class_group/factor_base_proof_targets.cpp
    src/class_group/ideal_lattice_lll.cpp
    src/class_group/ideal_lattice_reduction.cpp
    src/class_group/ideal_minkowski_embedding.cpp
    src/class_group/ideal_relation_witness.cpp
    src/class_group/lll_relation_search.cpp
    src/class_group/ideal_t2_enumeration.cpp
    src/class_group/relation_production.cpp
    src/class_group/maximal_imaginary_quadratic_relation_backend.cpp
    src/class_group/native_relation_search.cpp
    src/class_group/relation_admission_cache.cpp
    src/class_group/relation_completion_scheduler.cpp
    src/class_group/relation_factor_base_rebuild.cpp
    src/class_group/relation_factor_base_plan.cpp
    src/class_group/analytic_finish.cpp
    src/class_group/hnf_finish_workspace.cpp
    src/class_group/random_product_relation_search.cpp
    src/class_group/relation_search_driver.cpp
    src/class_group/small_norm_relation_search.cpp
    src/class_group/relation_candidate.cpp
    src/class_group/relation_search.cpp
    src/factored_element/compact_reconstruction.cpp
    src/factored_element/factored_element.cpp
    src/factor_base/factor_base.cpp
    src/archimedean/archimedean.cpp
    src/aut/aut.cpp
    src/embedding/embedding.cpp
    src/element/element.cpp
    src/hom/hom.cpp
    src/number_field/number_field.cpp
    src/signature/signature.cpp
    src/ideal/ideal.cpp
    src/ideal_factorization/ideal_factorization.cpp
    src/order/order.cpp
    src/order_element/order_element.cpp
    src/order_unit/class_unit_transaction.cpp
    src/order_unit/compact_reconstruction_bound.cpp
    src/order_unit/compute.cpp
    src/order_unit/order_unit.cpp
    src/order_unit/hnf_relation_coefficient_reduction.cpp
    src/order_unit/hnf_relation_coefficient_logs.cpp
    src/order_unit/hnf_unit_materialization.cpp
    src/order_unit/hnf_unit_log_sum_validation.cpp
    src/order_unit/hnf_regulator_acceptance.cpp
    src/order_unit/hnf_regulator_pivots.cpp
    src/order_unit/hnf_unit_publication.cpp
    src/order_unit/hnf_unit_log_matrix.cpp
    src/order_unit/hnf_unit_lattice_transform.cpp
    src/order_unit/hnf_reduced_regulator_reconstruction.cpp
    src/order_unit/hnf_class_regulator_product.cpp
    src/order_unit/relation_kernel.cpp
    src/order_unit/saturation.cpp
    src/order_unit/validation.cpp
    src/prime_ideal/prime_ideal.cpp
    src/relation/relation.cpp
    src/residue_field/degree_one_helpers.cpp
    src/residue_field/residue_field.cpp
    src/residue_ring/residue_ring.cpp
    src/status/status.cpp
    src/sunit/sunit.cpp
    src/sunit/sunit_class_context.cpp
    src/unit/roots_of_unity.cpp
    src/unit/unit.cpp
    src/version/version.cpp
    src/zeta/zeta.cpp
)
add_library(Silex::silex ALIAS silex)
set_target_properties(silex PROPERTIES EXPORT_NAME silex)

if(SILEX_WITH_FLATTER)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        set_source_files_properties(src/lat/flatter_backend.cpp PROPERTIES
            COMPILE_OPTIONS "-idirafter;${SILEX_FLATTER_INCLUDE_DIR}")
    else()
        set_source_files_properties(src/lat/flatter_backend.cpp PROPERTIES
            INCLUDE_DIRECTORIES "${SILEX_FLATTER_INCLUDE_DIR}")
    endif()
endif()

target_include_directories(silex
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(silex PUBLIC cxx_std_20)

# A static sanitizer-instrumented library requires the sanitizer runtime when
# linked by an installed consumer as well as by source-tree targets.
if(SILEX_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_link_options(silex INTERFACE -fsanitize=address,undefined)
endif()

target_link_libraries(silex
    PUBLIC
        PkgConfig::FLINT)
if(SILEX_FPLLL_TARGET)
    target_link_libraries(silex PRIVATE ${SILEX_FPLLL_TARGET})
endif()
if(SILEX_FLATTER_TARGET)
    target_link_libraries(silex PRIVATE ${SILEX_FLATTER_TARGET})
endif()

target_compile_definitions(silex
    PUBLIC
        $<$<BOOL:${SILEX_ENABLE_LOGGING}>:SILEX_ENABLE_LOGGING=1>
        $<$<BOOL:${SILEX_ENABLE_DEBUG_CHECKS}>:SILEX_ENABLE_DEBUG_CHECKS=1>
        $<$<BOOL:${SILEX_ENABLE_PROFILING}>:SILEX_ENABLE_PROFILING=1>
    PRIVATE
        $<$<BOOL:${SILEX_WITH_FPLLL}>:SILEX_WITH_FPLLL=1>
        $<$<BOOL:${SILEX_WITH_FLATTER}>:SILEX_WITH_FLATTER=1>
)

silex_apply_target_warnings(silex)

if(MSVC)
    if(NOT SILEX_ENABLE_EXCEPTIONS)
        target_compile_options(silex PRIVATE /EHs-c-)
    endif()
    if(NOT SILEX_ENABLE_RTTI)
        target_compile_options(silex PRIVATE /GR-)
    endif()
else()
    if(NOT SILEX_ENABLE_EXCEPTIONS)
        target_compile_options(silex PRIVATE -fno-exceptions)
    endif()
    if(NOT SILEX_ENABLE_RTTI)
        target_compile_options(silex PRIVATE -fno-rtti)
    endif()
endif()

if(SILEX_BUILD_BENCHMARK_ADAPTERS OR SILEX_BUILD_TESTS)
    add_executable(silex-class-unit-instance
        tools/class_unit_instance.cpp)
    target_link_libraries(silex-class-unit-instance PRIVATE Silex::silex)
    target_include_directories(silex-class-unit-instance PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(silex-class-unit-instance)

    add_executable(silex-operation-instance
        tools/operation_instance.cpp)
    target_link_libraries(silex-operation-instance PRIVATE Silex::silex)
    silex_apply_target_warnings(silex-operation-instance)
endif()

include(CTest)
if(SILEX_BUILD_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    set(SILEX_TEST_TIMEOUT_SECONDS 60)
    set(SILEX_TEST_NAMES)
    set(SILEX_RELEASE_TEST_NAMES)

    if(SILEX_BUILD_RELEASE_TESTS)
        add_test(NAME silex-install-consumer-default
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_installed_package.py
                    --cmake ${CMAKE_COMMAND}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                    --config=$<CONFIG>
                    --backend default
                    --layout default)
        add_test(NAME silex-install-consumer-nondefault-layout
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_installed_package.py
                    --cmake ${CMAKE_COMMAND}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                    --config=$<CONFIG>
                    --backend default
                    --layout nondefault)
        add_test(NAME silex-add-subdirectory-consumer
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_source_tree_consumer.py
                    --cmake ${CMAKE_COMMAND}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                    --config=$<CONFIG>)
        list(APPEND SILEX_RELEASE_TEST_NAMES
            silex-install-consumer-default
            silex-install-consumer-nondefault-layout
            silex-add-subdirectory-consumer)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
            add_test(NAME silex-source-archive
                COMMAND ${Python3_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/test/test_source_archive.py
                        --cmake ${CMAKE_COMMAND}
                        --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                        --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                        --project-version ${PROJECT_VERSION}
                        --config=$<CONFIG>)
            list(APPEND SILEX_RELEASE_TEST_NAMES silex-source-archive)
        else()
            message(STATUS
                "Silex source-archive release test not registered: "
                "the source root has no Git index")
        endif()
        list(APPEND SILEX_TEST_NAMES ${SILEX_RELEASE_TEST_NAMES})
    endif()

    foreach(SILEX_UNSUPPORTED_BACKEND IN ITEMS spdlog tracy)
        add_test(NAME silex-unsupported-${SILEX_UNSUPPORTED_BACKEND}-option
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_unsupported_cmake_options.py
                    --cmake ${CMAKE_COMMAND}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                    --unsupported-backend ${SILEX_UNSUPPORTED_BACKEND})
        list(APPEND SILEX_TEST_NAMES
            silex-unsupported-${SILEX_UNSUPPORTED_BACKEND}-option)
    endforeach()

    if(SILEX_WITH_FPLLL)
        add_test(NAME silex-install-consumer-fplll
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_installed_package.py
                    --cmake ${CMAKE_COMMAND}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                    --config=$<CONFIG>
                    --backend fplll
                    --fplll-root=${SILEX_FPLLL_ROOT}
                    --flatter-root=${SILEX_FLATTER_ROOT})
        list(APPEND SILEX_TEST_NAMES silex-install-consumer-fplll)
    endif()
    if(SILEX_WITH_FLATTER)
        add_test(NAME silex-install-consumer-flatter
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/test/test_installed_package.py
                    --cmake ${CMAKE_COMMAND}
                    --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
                    --build-dir ${CMAKE_CURRENT_BINARY_DIR}
                    --config=$<CONFIG>
                    --backend flatter
                    --fplll-root=${SILEX_FPLLL_ROOT}
                    --flatter-root=${SILEX_FLATTER_ROOT})
        list(APPEND SILEX_TEST_NAMES silex-install-consumer-flatter)
    endif()

    # Several source-boundary tests intentionally use assert for fixture setup.
    # Keep those checks active in optimized test builds as well.
    if(MSVC)
        add_compile_options(
            $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>:/UNDEBUG>)
    else()
        add_compile_options(
            $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>:-UNDEBUG>)
    endif()

    add_test(NAME silex-profile-check
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/check-cpp-profile.py
                --root ${CMAKE_CURRENT_SOURCE_DIR})
    list(APPEND SILEX_TEST_NAMES silex-profile-check)

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
        add_test(NAME silex-source-name-boundary
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/check-source-name-boundary.py
                    --root ${CMAKE_CURRENT_SOURCE_DIR})
        list(APPEND SILEX_TEST_NAMES silex-source-name-boundary)
    else()
        message(STATUS
            "Silex source-name boundary test not registered: "
            "the source root has no Git index")
    endif()

    add_test(NAME silex-native-benchmark-tool
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_native_benchmark.py
                --tool
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/bench/native-benchmark.py)
    list(APPEND SILEX_TEST_NAMES silex-native-benchmark-tool)

    add_test(NAME silex-class-unit-manifest-contract
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_class_unit_benchmark_manifest.py
                --manifest
                ${CMAKE_CURRENT_SOURCE_DIR}/test/data/class_unit_fields.json)
    list(APPEND SILEX_TEST_NAMES silex-class-unit-manifest-contract)

    add_test(NAME silex-class-unit-instance-protocol
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_class_unit_instance.py
                --exe $<TARGET_FILE:silex-class-unit-instance>
                --manifest
                ${CMAKE_CURRENT_SOURCE_DIR}/test/data/class_unit_fields.json)
    list(APPEND SILEX_TEST_NAMES silex-class-unit-instance-protocol)

    add_test(NAME silex-operation-instance
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_operation_instance.py
                --exe $<TARGET_FILE:silex-operation-instance>)
    list(APPEND SILEX_TEST_NAMES silex-operation-instance)

    add_test(NAME silex-sunit-instance-protocol
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_sunit_instance.py
                --exe $<TARGET_FILE:silex-class-unit-instance>
                --manifest
                ${CMAKE_CURRENT_SOURCE_DIR}/test/data/sunit_fields.json)
    list(APPEND SILEX_TEST_NAMES silex-sunit-instance-protocol)

    add_executable(t-silex-status test/t-status.cpp)
    target_link_libraries(t-silex-status PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-status)
    add_test(NAME silex-status COMMAND t-silex-status)
    list(APPEND SILEX_TEST_NAMES silex-status)

    add_executable(t-silex-flint test/t-flint.cpp)
    target_link_libraries(t-silex-flint PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-flint)
    add_test(NAME silex-flint COMMAND t-silex-flint)
    list(APPEND SILEX_TEST_NAMES silex-flint)

    add_executable(t-silex-flint-raii test/t-flint-raii.cpp)
    target_link_libraries(t-silex-flint-raii PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-flint-raii)
    add_test(NAME silex-flint-raii COMMAND t-silex-flint-raii)
    list(APPEND SILEX_TEST_NAMES silex-flint-raii)

    add_executable(t-silex-fmpz-smat test/t-fmpz-smat.cpp)
    target_link_libraries(t-silex-fmpz-smat PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-fmpz-smat)
    add_test(NAME silex-fmpz-smat COMMAND t-silex-fmpz-smat)
    list(APPEND SILEX_TEST_NAMES silex-fmpz-smat)

    add_executable(t-silex-lat test/t-lat.cpp)
    target_link_libraries(t-silex-lat PRIVATE Silex::silex)
    target_include_directories(t-silex-lat PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-lat)
    add_test(NAME silex-lat COMMAND t-silex-lat)
    list(APPEND SILEX_TEST_NAMES silex-lat)

    add_executable(t-silex-abelian-group test/t-abelian-group.cpp)
    target_link_libraries(t-silex-abelian-group PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-abelian-group)
    add_test(NAME silex-abelian-group COMMAND t-silex-abelian-group)
    list(APPEND SILEX_TEST_NAMES silex-abelian-group)

    add_executable(t-silex-class-group test/t-class-group.cpp)
    target_link_libraries(t-silex-class-group PRIVATE Silex::silex)
    target_include_directories(t-silex-class-group PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-class-group)
    add_test(NAME silex-class-group COMMAND t-silex-class-group)
    list(APPEND SILEX_TEST_NAMES silex-class-group)

    add_executable(t-silex-ideal-minkowski-embedding
        test/t-ideal-minkowski-embedding.cpp)
    target_link_libraries(t-silex-ideal-minkowski-embedding
        PRIVATE Silex::silex)
    target_include_directories(t-silex-ideal-minkowski-embedding PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-ideal-minkowski-embedding)
    add_test(NAME silex-ideal-minkowski-embedding
        COMMAND t-silex-ideal-minkowski-embedding)
    list(APPEND SILEX_TEST_NAMES
        silex-ideal-minkowski-embedding)

    add_executable(t-silex-finite-quadratic-form-enumeration
        test/t-finite-quadratic-form-enumeration.cpp)
    target_link_libraries(t-silex-finite-quadratic-form-enumeration
        PRIVATE Silex::silex)
    target_include_directories(
        t-silex-finite-quadratic-form-enumeration PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(
        t-silex-finite-quadratic-form-enumeration)
    add_test(NAME silex-finite-quadratic-form-enumeration
        COMMAND t-silex-finite-quadratic-form-enumeration)
    list(APPEND SILEX_TEST_NAMES
        silex-finite-quadratic-form-enumeration)

    add_executable(t-silex-ideal-lattice-reduction
        test/t-ideal-lattice-reduction.cpp)
    target_link_libraries(t-silex-ideal-lattice-reduction
        PRIVATE Silex::silex)
    target_include_directories(t-silex-ideal-lattice-reduction PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-ideal-lattice-reduction)
    add_test(NAME silex-ideal-lattice-reduction
        COMMAND t-silex-ideal-lattice-reduction)
    list(APPEND SILEX_TEST_NAMES silex-ideal-lattice-reduction)

    add_executable(t-silex-factor-base-honesty
        test/t-factor-base-honesty.cpp)
    target_link_libraries(t-silex-factor-base-honesty PRIVATE Silex::silex)
    target_include_directories(t-silex-factor-base-honesty PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-factor-base-honesty)
    add_test(NAME silex-factor-base-honesty
        COMMAND t-silex-factor-base-honesty)
    list(APPEND SILEX_TEST_NAMES silex-factor-base-honesty)

    add_executable(t-silex-relation-admission-cache
        test/t-relation-admission-cache.cpp)
    target_link_libraries(t-silex-relation-admission-cache
        PRIVATE Silex::silex)
    target_include_directories(t-silex-relation-admission-cache PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-relation-admission-cache)
    add_test(NAME silex-relation-admission-cache
        COMMAND t-silex-relation-admission-cache)
    list(APPEND SILEX_TEST_NAMES silex-relation-admission-cache)

    add_executable(t-silex-relation-factor-base-rebuild
        test/t-relation-factor-base-rebuild.cpp)
    target_link_libraries(t-silex-relation-factor-base-rebuild
        PRIVATE Silex::silex)
    target_include_directories(t-silex-relation-factor-base-rebuild
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-relation-factor-base-rebuild)
    add_test(NAME silex-relation-factor-base-rebuild
        COMMAND t-silex-relation-factor-base-rebuild)
    list(APPEND SILEX_TEST_NAMES
        silex-relation-factor-base-rebuild)

    add_executable(t-silex-relation-factor-base-plan
        test/t-relation-factor-base-plan.cpp)
    target_link_libraries(t-silex-relation-factor-base-plan
        PRIVATE Silex::silex)
    target_include_directories(t-silex-relation-factor-base-plan PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-relation-factor-base-plan)
    add_test(NAME silex-relation-factor-base-plan
        COMMAND t-silex-relation-factor-base-plan)
    list(APPEND SILEX_TEST_NAMES silex-relation-factor-base-plan)

    add_executable(t-silex-relation-completion-scheduler
        test/t-relation-completion-scheduler.cpp)
    target_link_libraries(t-silex-relation-completion-scheduler
        PRIVATE Silex::silex)
    target_include_directories(t-silex-relation-completion-scheduler
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-relation-completion-scheduler)
    add_test(NAME silex-relation-completion-scheduler
        COMMAND t-silex-relation-completion-scheduler)
    list(APPEND SILEX_TEST_NAMES
        silex-relation-completion-scheduler)

    add_executable(t-silex-class-group-ideal-relation
        test/t-class-group-ideal-relation.cpp)
    target_link_libraries(t-silex-class-group-ideal-relation
        PRIVATE Silex::silex)
    target_include_directories(t-silex-class-group-ideal-relation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-class-group-ideal-relation)
    add_test(NAME silex-class-group-ideal-relation
        COMMAND t-silex-class-group-ideal-relation)
    list(APPEND SILEX_TEST_NAMES
        silex-class-group-ideal-relation)

    add_executable(t-silex-class-group-hnf-witness
        test/t-class-group-hnf-witness.cpp)
    target_link_libraries(t-silex-class-group-hnf-witness
        PRIVATE Silex::silex)
    target_include_directories(t-silex-class-group-hnf-witness PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-class-group-hnf-witness)
    add_test(NAME silex-class-group-hnf-witness
        COMMAND t-silex-class-group-hnf-witness)
    list(APPEND SILEX_TEST_NAMES silex-class-group-hnf-witness)

    add_executable(t-silex-sunit-class-context
        test/t-sunit-class-context.cpp)
    target_link_libraries(t-silex-sunit-class-context PRIVATE Silex::silex)
    target_include_directories(t-silex-sunit-class-context PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-sunit-class-context)
    add_test(NAME silex-sunit-class-context
        COMMAND t-silex-sunit-class-context)
    list(APPEND SILEX_TEST_NAMES silex-sunit-class-context)

    add_executable(t-silex-sunit test/t-sunit.cpp)
    target_link_libraries(t-silex-sunit PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-sunit)
    add_test(NAME silex-sunit COMMAND t-silex-sunit)
    list(APPEND SILEX_TEST_NAMES silex-sunit)

    add_executable(t-silex-class-relation-module
        test/t-class-relation-module.cpp)
    target_link_libraries(t-silex-class-relation-module PRIVATE Silex::silex)
    target_include_directories(t-silex-class-relation-module PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-class-relation-module)
    add_test(NAME silex-class-relation-module
        COMMAND t-silex-class-relation-module)
    list(APPEND SILEX_TEST_NAMES silex-class-relation-module)

    add_executable(t-silex-class-unit-matrix
        test/t-class-unit-matrix.cpp)
    target_link_libraries(t-silex-class-unit-matrix PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-class-unit-matrix)
    add_test(NAME silex-class-unit-matrix
        COMMAND t-silex-class-unit-matrix)
    list(APPEND SILEX_TEST_NAMES silex-class-unit-matrix)

    add_executable(t-silex-number-field test/t-number-field.cpp)
    target_link_libraries(t-silex-number-field PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-number-field)
    add_test(NAME silex-number-field COMMAND t-silex-number-field)
    list(APPEND SILEX_TEST_NAMES silex-number-field)

    add_executable(t-silex-element test/t-element.cpp)
    target_link_libraries(t-silex-element PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-element)
    add_test(NAME silex-element COMMAND t-silex-element)
    list(APPEND SILEX_TEST_NAMES silex-element)

    add_executable(t-silex-factored-element test/t-factored-element.cpp)
    target_link_libraries(t-silex-factored-element PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-factored-element)
    add_test(NAME silex-factored-element COMMAND t-silex-factored-element)
    list(APPEND SILEX_TEST_NAMES silex-factored-element)

    add_executable(t-silex-compact-reconstruction
        test/t-compact-reconstruction.cpp)
    target_link_libraries(t-silex-compact-reconstruction
        PRIVATE Silex::silex)
    target_include_directories(t-silex-compact-reconstruction PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-compact-reconstruction)
    add_test(NAME silex-compact-reconstruction
        COMMAND t-silex-compact-reconstruction)
    list(APPEND SILEX_TEST_NAMES silex-compact-reconstruction)

    add_executable(t-silex-hom test/t-hom.cpp)
    target_link_libraries(t-silex-hom PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-hom)
    add_test(NAME silex-hom COMMAND t-silex-hom)
    list(APPEND SILEX_TEST_NAMES silex-hom)

    add_executable(t-silex-aut test/t-aut.cpp)
    target_link_libraries(t-silex-aut PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-aut)
    add_test(NAME silex-aut COMMAND t-silex-aut)
    list(APPEND SILEX_TEST_NAMES silex-aut)

    add_executable(t-silex-signature test/t-signature.cpp)
    target_link_libraries(t-silex-signature PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-signature)
    add_test(NAME silex-signature COMMAND t-silex-signature)
    list(APPEND SILEX_TEST_NAMES silex-signature)

    add_executable(t-silex-embedding test/t-embedding.cpp)
    target_link_libraries(t-silex-embedding PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-embedding)
    add_test(NAME silex-embedding COMMAND t-silex-embedding)
    list(APPEND SILEX_TEST_NAMES silex-embedding)

    add_executable(t-silex-archimedean test/t-archimedean.cpp)
    target_link_libraries(t-silex-archimedean PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-archimedean)
    add_test(NAME silex-archimedean COMMAND t-silex-archimedean)
    list(APPEND SILEX_TEST_NAMES silex-archimedean)

    add_executable(t-silex-unit test/t-unit.cpp)
    target_link_libraries(t-silex-unit PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-unit)
    add_test(NAME silex-unit COMMAND t-silex-unit)
    list(APPEND SILEX_TEST_NAMES silex-unit)

    add_executable(t-silex-order test/t-order.cpp)
    target_link_libraries(t-silex-order PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-order)
    add_test(NAME silex-order COMMAND t-silex-order)
    list(APPEND SILEX_TEST_NAMES silex-order)

    add_executable(t-silex-order-element test/t-order-element.cpp)
    target_link_libraries(t-silex-order-element PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-order-element)
    add_test(NAME silex-order-element COMMAND t-silex-order-element)
    list(APPEND SILEX_TEST_NAMES silex-order-element)

    add_executable(t-silex-order-unit test/t-order-unit.cpp)
    target_link_libraries(t-silex-order-unit PRIVATE Silex::silex)
    target_include_directories(t-silex-order-unit PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-order-unit)
    add_test(NAME silex-order-unit COMMAND t-silex-order-unit)
    list(APPEND SILEX_TEST_NAMES silex-order-unit)

    add_executable(t-silex-order-unit-coordinates
        test/t-order-unit-coordinates.cpp)
    target_link_libraries(t-silex-order-unit-coordinates PRIVATE Silex::silex)
    target_include_directories(t-silex-order-unit-coordinates PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-order-unit-coordinates)
    add_test(NAME silex-order-unit-coordinates
        COMMAND t-silex-order-unit-coordinates)
    list(APPEND SILEX_TEST_NAMES silex-order-unit-coordinates)

    add_executable(t-silex-zeta test/t-zeta.cpp)
    target_link_libraries(t-silex-zeta PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-zeta)
    add_test(NAME silex-zeta COMMAND t-silex-zeta)
    list(APPEND SILEX_TEST_NAMES silex-zeta)

    add_executable(t-silex-zeta-validation
        test/t-zeta-validation.cpp)
    target_link_libraries(t-silex-zeta-validation PRIVATE Silex::silex)
    target_include_directories(t-silex-zeta-validation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-zeta-validation)
    add_test(NAME silex-zeta-validation
        COMMAND t-silex-zeta-validation)
    list(APPEND SILEX_TEST_NAMES silex-zeta-validation)

    add_executable(t-silex-ideal test/t-ideal.cpp)
    target_link_libraries(t-silex-ideal PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-ideal)
    add_test(NAME silex-ideal COMMAND t-silex-ideal)
    list(APPEND SILEX_TEST_NAMES silex-ideal)

    add_executable(t-silex-prime-ideal test/t-prime-ideal.cpp)
    target_link_libraries(t-silex-prime-ideal PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-prime-ideal)
    add_test(NAME silex-prime-ideal COMMAND t-silex-prime-ideal)
    list(APPEND SILEX_TEST_NAMES silex-prime-ideal)

    add_executable(t-silex-factor-base test/t-factor-base.cpp)
    target_link_libraries(t-silex-factor-base PRIVATE Silex::silex)
    target_include_directories(t-silex-factor-base
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-factor-base)
    add_test(NAME silex-factor-base COMMAND t-silex-factor-base)
    list(APPEND SILEX_TEST_NAMES silex-factor-base)

    add_executable(t-silex-ideal-factorization test/t-ideal-factorization.cpp)
    target_link_libraries(t-silex-ideal-factorization PRIVATE Silex::silex)
    target_include_directories(t-silex-ideal-factorization
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    silex_apply_target_warnings(t-silex-ideal-factorization)
    add_test(NAME silex-ideal-factorization COMMAND t-silex-ideal-factorization)
    list(APPEND SILEX_TEST_NAMES silex-ideal-factorization)

    add_executable(t-silex-relation test/t-relation.cpp)
    target_link_libraries(t-silex-relation PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-relation)
    add_test(NAME silex-relation COMMAND t-silex-relation)
    list(APPEND SILEX_TEST_NAMES silex-relation)

    add_executable(t-silex-residue-field test/t-residue-field.cpp)
    target_link_libraries(t-silex-residue-field PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-residue-field)
    add_test(NAME silex-residue-field COMMAND t-silex-residue-field)
    list(APPEND SILEX_TEST_NAMES silex-residue-field)

    add_executable(t-silex-residue-ring test/t-residue-ring.cpp)
    target_link_libraries(t-silex-residue-ring PRIVATE Silex::silex)
    silex_apply_target_warnings(t-silex-residue-ring)
    add_test(NAME silex-residue-ring COMMAND t-silex-residue-ring)
    list(APPEND SILEX_TEST_NAMES silex-residue-ring)

    add_executable(t-silex-diagnostics test/t-diagnostics.cpp)
    target_link_libraries(t-silex-diagnostics PRIVATE Silex::silex)
    target_compile_definitions(t-silex-diagnostics
        PRIVATE
            SILEX_ENABLE_LOGGING=1)
    silex_apply_target_warnings(t-silex-diagnostics)
    add_test(NAME silex-diagnostics COMMAND t-silex-diagnostics)
    list(APPEND SILEX_TEST_NAMES silex-diagnostics)

    add_executable(t-silex-diagnostics-enabled test/t-diagnostics-enabled.cpp)
    target_link_libraries(t-silex-diagnostics-enabled PRIVATE Silex::silex)
    target_compile_definitions(t-silex-diagnostics-enabled
        PRIVATE
            SILEX_ENABLE_LOGGING=1
            SILEX_ENABLE_DEBUG_CHECKS=1
            SILEX_ENABLE_PROFILING=1)
    silex_apply_target_warnings(t-silex-diagnostics-enabled)
    add_test(NAME silex-diagnostics-enabled COMMAND t-silex-diagnostics-enabled)
    list(APPEND SILEX_TEST_NAMES silex-diagnostics-enabled)

    set_tests_properties(${SILEX_TEST_NAMES}
        PROPERTIES TIMEOUT ${SILEX_TEST_TIMEOUT_SECONDS})
    if(SILEX_RELEASE_TEST_NAMES)
        set_tests_properties(${SILEX_RELEASE_TEST_NAMES}
            PROPERTIES TIMEOUT 300)
    endif()
endif()

if(SILEX_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(SILEX_BUILD_DOCS)
    find_program(SPHINX_BUILD NAMES sphinx-build REQUIRED)
    if(NOT SILEX_DOCS_RELEASE)
        message(FATAL_ERROR
            "SILEX_DOCS_RELEASE must be nonempty when SILEX_BUILD_DOCS=ON")
    endif()
    if(NOT SILEX_DOCS_CHANNEL MATCHES "^(dev|stable)$")
        message(FATAL_ERROR
            "SILEX_DOCS_CHANNEL must be either 'dev' or 'stable'")
    endif()

    add_custom_target(docs-html
        COMMAND ${CMAKE_COMMAND} -E env
                "SILEX_DOCS_RELEASE=${SILEX_DOCS_RELEASE}"
                "SILEX_DOCS_CHANNEL=${SILEX_DOCS_CHANNEL}"
                ${SPHINX_BUILD} -n -W --keep-going
                -d ${CMAKE_CURRENT_BINARY_DIR}/docs/doctrees/html -b html
                ${CMAKE_CURRENT_SOURCE_DIR}/docs
                ${CMAKE_CURRENT_BINARY_DIR}/docs/html
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Building Silex HTML documentation"
        VERBATIM)
    add_custom_target(docs-linkcheck
        COMMAND ${CMAKE_COMMAND} -E env
                "SILEX_DOCS_RELEASE=${SILEX_DOCS_RELEASE}"
                "SILEX_DOCS_CHANNEL=${SILEX_DOCS_CHANNEL}"
                ${SPHINX_BUILD} -n -W --keep-going
                -d ${CMAKE_CURRENT_BINARY_DIR}/docs/doctrees/linkcheck
                -b linkcheck
                ${CMAKE_CURRENT_SOURCE_DIR}/docs
                ${CMAKE_CURRENT_BINARY_DIR}/docs/linkcheck
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Checking links in Silex documentation"
        VERBATIM)
endif()

if(SILEX_BUILD_BENCHMARKS)
    add_subdirectory(bench)
endif()

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SilexConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY ExactVersion)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SilexConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SilexConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Silex)

install(TARGETS silex
    EXPORT SilexTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/silex
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/silex/build_config.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/silex)

install(EXPORT SilexTargets
    NAMESPACE Silex::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Silex)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SilexConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SilexConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Silex)

install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
    ${CMAKE_CURRENT_SOURCE_DIR}/NOTICE.md
    ${CMAKE_CURRENT_SOURCE_DIR}/CITATION.cff
    ${CMAKE_CURRENT_SOURCE_DIR}/THIRD_PARTY_NOTICES.md
    DESTINATION ${CMAKE_INSTALL_DOCDIR})
