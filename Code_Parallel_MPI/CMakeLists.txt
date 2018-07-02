cmake_minimum_required(VERSION 3.5)
project(Original_vedit)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp bath_setup.cpp propagation.cpp transition_matrix.cpp density.cpp opt_parser.cpp)
add_executable(Original_vedit ${SOURCE_FILES})
find_package(GSL REQUIRED)
target_link_libraries(Original_vedit GSL::gsl GSL::gslcblas)

include_directories(${GSL_INCLUDE_DIR})
