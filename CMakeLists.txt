cmake_minimum_required(VERSION 3.25)
project(RTX)

set(CMAKE_CXX_STANDARD 23)

find_package(SFML 2.5
        COMPONENTS
        system window graphics REQUIRED)
add_executable(RTX main.cpp)
target_link_libraries(RTX sfml-graphics sfml-window sfml-system)
