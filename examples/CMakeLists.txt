add_executable(example-factored-compact-roots factored_compact_roots.cpp)
target_link_libraries(example-factored-compact-roots PRIVATE Silex::silex)
silex_apply_target_warnings(example-factored-compact-roots)

add_executable(example-field-order-ideal-basics field_order_ideal_basics.cpp)
target_link_libraries(example-field-order-ideal-basics PRIVATE Silex::silex)
silex_apply_target_warnings(example-field-order-ideal-basics)

add_executable(example-element-arithmetic-basics element_arithmetic_basics.cpp)
target_link_libraries(example-element-arithmetic-basics PRIVATE Silex::silex)
silex_apply_target_warnings(example-element-arithmetic-basics)

add_executable(example-class-unit-computation class_unit_computation.cpp)
target_link_libraries(example-class-unit-computation PRIVATE Silex::silex)
silex_apply_target_warnings(example-class-unit-computation)

add_executable(example-diagnostics-hot-paths diagnostics_hot_paths.cpp)
target_link_libraries(example-diagnostics-hot-paths PRIVATE Silex::silex)
silex_apply_target_warnings(example-diagnostics-hot-paths)

add_executable(example-log-unit-lattice log_unit_lattice.cpp)
target_link_libraries(example-log-unit-lattice PRIVATE Silex::silex)
silex_apply_target_warnings(example-log-unit-lattice)

add_executable(example-local-relation-basics local_relation_basics.cpp)
target_link_libraries(example-local-relation-basics PRIVATE Silex::silex)
silex_apply_target_warnings(example-local-relation-basics)

add_executable(example-maximal-order-zeta maximal_order_zeta.cpp)
target_link_libraries(example-maximal-order-zeta PRIVATE Silex::silex)
silex_apply_target_warnings(example-maximal-order-zeta)

add_executable(example-rank-zero-order-units rank_zero_order_units.cpp)
target_link_libraries(example-rank-zero-order-units PRIVATE Silex::silex)
silex_apply_target_warnings(example-rank-zero-order-units)

add_executable(example-real-quadratic-order-units real_quadratic_order_units.cpp)
target_link_libraries(example-real-quadratic-order-units PRIVATE Silex::silex)
silex_apply_target_warnings(example-real-quadratic-order-units)

add_executable(example-relation-kernel-order-units relation_kernel_order_units.cpp)
target_link_libraries(example-relation-kernel-order-units PRIVATE Silex::silex)
silex_apply_target_warnings(example-relation-kernel-order-units)

add_executable(example-supplied-order-units supplied_order_units.cpp)
target_link_libraries(example-supplied-order-units PRIVATE Silex::silex)
silex_apply_target_warnings(example-supplied-order-units)

if(SILEX_BUILD_TESTS)
    set(SILEX_EXAMPLE_TEST_NAMES)

    add_test(NAME silex-example-factored-compact-roots
        COMMAND example-factored-compact-roots)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-factored-compact-roots)
    add_test(NAME silex-example-field-order-ideal-basics
        COMMAND example-field-order-ideal-basics)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-field-order-ideal-basics)
    add_test(NAME silex-example-element-arithmetic-basics
        COMMAND example-element-arithmetic-basics)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-element-arithmetic-basics)
    add_test(NAME silex-example-class-unit-computation
        COMMAND example-class-unit-computation)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-class-unit-computation)
    add_test(NAME silex-example-diagnostics-hot-paths
        COMMAND example-diagnostics-hot-paths)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-diagnostics-hot-paths)
    add_test(NAME silex-example-log-unit-lattice
        COMMAND example-log-unit-lattice)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-log-unit-lattice)
    add_test(NAME silex-example-local-relation-basics
        COMMAND example-local-relation-basics)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-local-relation-basics)
    add_test(NAME silex-example-maximal-order-zeta
        COMMAND example-maximal-order-zeta)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-maximal-order-zeta)
    add_test(NAME silex-example-rank-zero-order-units
        COMMAND example-rank-zero-order-units)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-rank-zero-order-units)
    add_test(NAME silex-example-real-quadratic-order-units
        COMMAND example-real-quadratic-order-units)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-real-quadratic-order-units)
    add_test(NAME silex-example-relation-kernel-order-units
        COMMAND example-relation-kernel-order-units)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-relation-kernel-order-units)
    add_test(NAME silex-example-supplied-order-units
        COMMAND example-supplied-order-units)
    list(APPEND SILEX_EXAMPLE_TEST_NAMES
        silex-example-supplied-order-units)

    set_tests_properties(${SILEX_EXAMPLE_TEST_NAMES}
        PROPERTIES TIMEOUT ${SILEX_TEST_TIMEOUT_SECONDS})
endif()
