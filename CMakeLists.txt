cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CountIntegerPointsSimplePolytopes VERSION 1.0)

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(Count ${SOURCES})

include_directories(include)

target_link_libraries(Count PUBLIC gmp gmpxx cdd flint)
