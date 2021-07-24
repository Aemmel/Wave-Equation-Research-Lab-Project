cmake_minimum_required(VERSION 3.9)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project(wave LANGUAGES CXX)

find_package(OpenMP REQUIRED)

add_executable(wave
    src/main.cpp
    src/numerics.cpp
    src/csv_printer.cpp
)

# add_executable(wave
#     src/preliminaries/main.cpp
#     src/csv_printer.cpp
# )

target_link_libraries(wave PRIVATE OpenMP::OpenMP_CXX)

target_compile_features(wave PRIVATE cxx_std_14) # probably C++11 would be enough

target_compile_options(wave PRIVATE -Wall -Wextra -Wpedantic)