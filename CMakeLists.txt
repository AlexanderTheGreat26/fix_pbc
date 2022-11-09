cmake_minimum_required(VERSION 3.23)
project(fix_pbc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

add_executable(fix_pbc main.cpp)
