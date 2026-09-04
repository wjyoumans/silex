cmake_minimum_required(VERSION 3.20)
project(SilexAddSubdirectoryConsumer LANGUAGES CXX)

if(NOT SILEX_SOURCE_DIR)
    message(FATAL_ERROR "SILEX_SOURCE_DIR must name the Silex source tree")
endif()

set(SILEX_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SILEX_BUILD_RELEASE_TESTS OFF CACHE BOOL "" FORCE)
set(SILEX_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(SILEX_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(SILEX_BUILD_BENCHMARK_ADAPTERS OFF CACHE BOOL "" FORCE)
set(SILEX_BUILD_DOCS OFF CACHE BOOL "" FORCE)
add_subdirectory("${SILEX_SOURCE_DIR}" silex-build EXCLUDE_FROM_ALL)

add_executable(silex-add-subdirectory-consumer main.cpp)
target_link_libraries(silex-add-subdirectory-consumer PRIVATE Silex::silex)
