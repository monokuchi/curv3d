
cmake_minimum_required(VERSION 3.13)
project(curv3d)

add_compile_options(-g -O0)

set(SRC_FILES src/bezier.cpp)

#############################################################

add_executable(curv3d ${SRC_FILES})

target_include_directories(curv3d PUBLIC
                           "${CMAKE_SOURCE_DIR}/src"
                          )
