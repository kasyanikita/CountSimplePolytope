cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(CountIntegerPointsSimplePolytopes VERSION 1.0)

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(Count ${SOURCES})

include_directories(include)

target_link_libraries(Count PUBLIC gmp gmpxx cdd flint)
# target_link_libraries(Main PUBLIC gmp gmpxx flint mpfr)
# target_link_libraries(UnitExperiment PUBLIC gmp gmpxx flint mpfr)
# target_link_libraries(TimeExperimentSimpex PUBLIC gmp gmpxx flint mpfr)
# target_link_libraries(TimeExperimentNParallel PUBLIC gmp gmpxx flint mpfr)
# target_link_libraries(StressTest PUBLIC gmp gmpxx flint mpfr)
# target_link_libraries(Triangulation PUBLIC gmp gmpxx flint mpfr cdd qhullcpp qhullstatic_r)
# target_link_libraries(Triangulation_experiment_3d PUBLIC gmp gmpxx flint mpfr cdd qhullcpp qhullstatic_r)
# cdd qhullcpp qhullstatic_r 
