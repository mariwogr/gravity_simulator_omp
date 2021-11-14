cmake_minimum_required(VERSION 3.20)
project(gravity_simulator)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)
add_executable(psoa psoa.cpp)
target_link_libraries (psoa PUBLIC OpenMP::OpenMP_CXX)
add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG -DCMAKE_BUILD_TYPE=Release -fopenmp -openmp)

add_executable(gravity_simulator psoa.cpp psoa.hpp)
