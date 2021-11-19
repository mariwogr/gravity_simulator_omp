cmake_minimum_required(VERSION 3.20)
project(gravity_simulator)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG -fopenmp)
find_package(OpenMP REQUIRED)

add_executable(paos paos.cpp)
target_link_libraries (paos PUBLIC OpenMP::OpenMP_CXX)
