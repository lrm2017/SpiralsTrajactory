#cmake_minimum_required(VERSION 3.0)

# Find necessary libs.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(Qt5 REQUIRED Core Gui OpenGL Charts)

# Build ploy library.
add_library(PlotLib

        include/plot2d.hpp
        src/plot2d.cpp

        include/plot3d.hpp
        src/plot3d.cpp
        )

find_package(Qt5 REQUIRED COMPONENTS Charts Gui Core OpenGL)

# Grab headers from the include directory.
target_include_directories(PlotLib PUBLIC include)

# Grab headers from this src directory. Don't share them with clients.
target_include_directories(PlotLib PRIVATE src)

# Pull in Eigen for linear algebra.
# Pull in CppOptLib for optimization.
#target_link_libraries(PlotLib )

#LINK_DIRECTORIES(/lib/x86_64-linux-gnu)
# Set options for the compiler.
target_compile_options(PlotLib PRIVATE -Wall)

target_link_libraries(PlotLib
        Qt5::Gui Qt5::Core Qt5::OpenGL Qt5::Charts
        Eigen3::Eigen
        )
