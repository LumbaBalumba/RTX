cmake_minimum_required(VERSION 3.25)
project(RTX)

set(CMAKE_CXX_STANDARD 23)

find_package(SFML 2.5
        COMPONENTS
        system window graphics REQUIRED)
find_package(OpenMP)
add_executable(RTX src/main.cpp src/Vec3.h src/Vec3.cpp src/Color.h src/Color.cpp src/Figure.h src/Figure.cpp src/Sphere.h src/Sphere.cpp src/Box.h src/Box.cpp src/Tetra.h src/Tetra.cpp src/Ray.h src/Ray.cpp src/Image.h src/Image.cpp)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
target_link_libraries(RTX sfml-graphics sfml-window sfml-system)
