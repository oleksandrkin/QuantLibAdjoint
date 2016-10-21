file (GLOB TEST_SUITE_ADJOINT_FILES "*.hpp" "*.cpp")

list(REMOVE_ITEM TEST_SUITE_ADJOINT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/adjointarraytest.cpp)
list(REMOVE_ITEM TEST_SUITE_ADJOINT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/adjointcomplextest.cpp)
list(REMOVE_ITEM TEST_SUITE_ADJOINT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/adjointcurvefittingtest.cpp)
list(REMOVE_ITEM TEST_SUITE_ADJOINT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/adjointgreektest.cpp)
list(REMOVE_ITEM TEST_SUITE_ADJOINT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/adjointgreekstest.cpp)
list(REMOVE_ITEM TEST_SUITE_ADJOINT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/adjointspecialfunctionstest.cpp)

set (TEST_ADJOINT quantlib-test-suite-adjoint)
add_executable (${TEST_ADJOINT} ${TEST_SUITE_ADJOINT_FILES})
target_link_libraries (${TEST_ADJOINT} QuantLib BoostLib)
enable_testing ()
add_test (${TEST_ADJOINT} ${TEST_ADJOINT})
