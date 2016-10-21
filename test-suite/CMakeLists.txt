file (GLOB TEST_SUITE_FILES "*.hpp" "*.cpp")

# these do not appear in vcproj
list(REMOVE_ITEM TEST_SUITE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/quantlibbenchmark.cpp)

set (TEST quantlib-test-suite)
add_executable (${TEST} ${TEST_SUITE_FILES})
target_link_libraries (${TEST} QuantLib BoostLib)
enable_testing ()
add_test (${TEST} ${TEST})
