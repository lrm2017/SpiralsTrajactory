cmake_minimum_required(VERSION 3.19)
project(spiralsTrajactory)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(lib)

add_executable(spiralsTrajactory main.cpp)
target_link_libraries(spiralsTrajactory PolyTrajLib)