cmake_minimum_required(VERSION 3.20)
project(SilexInstallConsumer LANGUAGES CXX)

find_package(Silex 0.1.1 EXACT CONFIG REQUIRED)

set(_silex_incompatible_source "${CMAKE_CURRENT_BINARY_DIR}/incompatible-source")
set(_silex_incompatible_build "${CMAKE_CURRENT_BINARY_DIR}/incompatible-build")
file(REMOVE_RECURSE "${_silex_incompatible_source}" "${_silex_incompatible_build}")
file(MAKE_DIRECTORY "${_silex_incompatible_source}")
file(WRITE "${_silex_incompatible_source}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.20)
project(SilexIncompatibleConsumer LANGUAGES CXX)
find_package(Silex 0.0.9 CONFIG REQUIRED)
]=])
execute_process(
    COMMAND ${CMAKE_COMMAND}
        -S "${_silex_incompatible_source}"
        -B "${_silex_incompatible_build}"
        "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
    RESULT_VARIABLE _silex_incompatible_result
    OUTPUT_QUIET
    ERROR_QUIET)
if(_silex_incompatible_result EQUAL 0)
    message(FATAL_ERROR
        "Silex 0.1.1 incorrectly accepted an incompatible 0.0.9 request")
endif()

add_executable(silex-install-consumer main.cpp)
target_link_libraries(silex-install-consumer PRIVATE Silex::silex)
