
cmake_minimum_required(VERSION 3.13)
project(curv3d)

add_compile_options(-g -O0)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# set(SRC_FILES src/bezier.cpp src/util.cpp)
set(SRC_FILES src/rasterize.cpp src/util.cpp src/bezier.cpp)

#############################################################

add_executable(curv3d ${SRC_FILES})

target_link_libraries (curv3d Eigen3::Eigen)

target_include_directories(curv3d PUBLIC
                           "${CMAKE_SOURCE_DIR}/src"
                          )
