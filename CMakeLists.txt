cmake_minimum_required(VERSION 3.0.2)

option(no_warnings      	"disable all warnings"  OFF)
option(no_pragmas       	"disable pragmas"       OFF)
option(debug_symbols    	"force debug symbols"   OFF)
option(use_complex      	"turn complex math on"  OFF)
option(x64              	"force x64 mode"        ON)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++1y" COMPILER_SUPPORTS_CXX1Y)
if(COMPILER_SUPPORTS_CXX14)
    add_compile_options(-std=c++14)
elseif(COMPILER_SUPPORTS_CXX1Y)
    add_compile_options(-std=c++1y)
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

if(x64)
    add_compile_options(-m64)
else()
    add_compile_options(-m32)
endif()

if(no_warnings)
    add_compile_options(-w)
endif()

if(no_pragmas)
    add_definitions(-DNO_PRAGMA_MSG)
endif()

if(debug_symbols)
    add_compile_options(-g)
endif()


add_compile_options(-Wno-placement-ne -Wno-deprecated-declarations -fpermissive)

add_definitions(-DBOOST_TEST_STATIC_LINK)
add_definitions(-DNDEBUG)
add_definitions(-DBOOST_FILESYSTEM_NO_LIB -DBOOST_SYSTEM_NO_LIB -DQL_ADJOINT -DCL_TAPE -DCL_TAPE_CPPAD -DCL_TAPE_CAN_GET_VALUE -DCL_TAPE_COMPLEX_ENABLED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(dependencies/cpp tapescript/cpp/ tapescript/cpp/cl/tape/impl/)

# system build-in boost
# find_package(Boost REQUIRED COMPONENTS serialization filesystem system unit_test_framework)
# include_directories( ${Boost_INCLUDE_DIR} )

add_subdirectory(dependencies/cpp/boost)
add_subdirectory(ql)
add_subdirectory(Examples)
add_subdirectory(test-suite)
add_subdirectory(test-suite-adjoint)
