cmake_minimum_required(VERSION 3.16)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project(wave
    LANGUAGES CXX)

add_executable(wave
    src/main.cpp
    src/csv_printer.cpp
)

target_compile_features(wave PRIVATE cxx_std_14) # probably C++11 would be enough

target_compile_options(wave PRIVATE -Wall -Wextra -Wpedantic)