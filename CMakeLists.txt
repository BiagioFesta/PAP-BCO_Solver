cmake_minimum_required (VERSION 2.6)

# project name
PROJECT( pap-bco_solver )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 1 )

# Find Boost Graph
find_package(Boost COMPONENTS graph serialization)
link_directories(${Boost_LIBRARY_DIRS})
include_directories(BEFORE ${Boost_INCLUDE_DIRS})

# Exec
add_executable(${PROJECT_NAME} "./src/pap-bco_solver.cpp")

# Linker Libraries
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU"   OR   "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
        set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS "-std=c++11 ")
    elseif(COMPILER_SUPPORTS_CXX0X)
        set_property(TARGET ${PROJECT_NAME} APPEND_STRING PROPERTY COMPILE_FLAGS "-std=c++0x")
    else()
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -D_DEBUG")
endif()

# Option for testing
option(BUILD_TEST "Build testing projects" OFF)
if (BUILD_TEST)
   add_subdirectory(test)
endif()
