cmake_minimum_required(VERSION 3.25)
project(RTX)

set(CMAKE_CXX_STANDARD 23)

find_package(SFML 2.5
        COMPONENTS
        system window graphics REQUIRED)
find_package(OpenMP)
add_executable(RTX main.cpp)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
target_link_libraries(RTX sfml-graphics sfml-window sfml-system)