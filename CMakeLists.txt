cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR GCP)
project(${PROJECT_NAME_STR} C CXX)

# C++11 required
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else ()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

#normal includes: warnings will be reported for these
include_directories(
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/ext/ranker
        ${source_dir}/googletest/include
        ${source_dir}/googlemock/include
)

file(GLOB_RECURSE SRC_FILES
    ${PROJECT_SOURCE_DIR}/src/*.h
    ${PROJECT_SOURCE_DIR}/ext/ranker/*.h
)

add_executable(ExampleGenotypeConfidencePercentiler ${SRC_FILES} ${PROJECT_SOURCE_DIR}/examples/Percentiler_only.cpp)
add_executable(FullExample ${SRC_FILES} ${PROJECT_SOURCE_DIR}/examples/Full.cpp)

enable_testing()
add_subdirectory(tests)
