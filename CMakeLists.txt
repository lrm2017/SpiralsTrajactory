cmake_minimum_required(VERSION 3.16)
project(spiral_constraint)

set(CMAKE_CXX_STANDARD 17)

find_package(casadi REQUIRED)   # 引入包，需要安装到系统环境才有

#set(CMAKE_CXX_COMPILER  "clang++" ) # 指定编译器
#如果要生成debug版本
SET(CMAKE_BUILD_TYPE "Debug")
#如果要生成release版本
#SET(CMAKE_BUILD_TYPE "Release")
#添加的支持debug和release支持
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -Os -Wall -g")

find_package(Eigen3 REQUIRED)

find_package(Matplot++ REQUIRED)

include(CheckIPOSupported)
check_ipo_supported(RESULT _IsIPOSupported)
if(_IsIPOSupported)
    message(STATUS "Turn on INTERPROCEDURAL_OPTIMIZATION")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()


option(BUILD_SHARED_LIBS   "Build as shared library." OFF)
option(BUILD_EXAMPLES      "Build example programs."   OFF)
option(BUILD_TESTS         "Build test programs."      OFF)
option(BUILD_DOCS          "Build documentation."      OFF)

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

#add_subdirectory(datavisual)    # qt 画图， 非必要

add_library(cpp src/directCollocationSolver.cpp)
target_link_libraries(cpp casadi)

add_executable(spiral_constraint main.cpp )
target_link_libraries(spiral_constraint cpp Eigen3::Eigen
        Matplot++::matplot)    # PlotLib 是画图库，非必要

add_executable(circeMP circleMPs.cpp)
target_link_libraries(
        circeMP cpp Eigen3::Eigen Matplot++::matplot
)

add_executable(casatest casatest.cpp)
target_link_libraries(
        casatest cpp Eigen3::Eigen Matplot++::matplot
)

