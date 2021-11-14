cmake_minimum_required(VERSION 3.20)
project(gravity_simulator)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG -DCMAKE_BUILD_TYPE=Release)

add_executable(gravity_simulator sim_soa.cpp sim_soa.hpp)
