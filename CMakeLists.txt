cmake_minimum_required(VERSION 3.10)
project(con2prim)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Conserved -> primitive, with the fork (funcs.hpp) and/or KHARMA (head.hpp) inverter
add_executable(con2prim src/main.cpp)
target_include_directories(con2prim PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Primitive -> conserved -> primitive round trip with both inverters
add_executable(prim2con src/prim2con.cpp)
target_include_directories(prim2con PRIVATE ${CMAKE_SOURCE_DIR}/include)
