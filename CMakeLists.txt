cmake_minimum_required(VERSION 3.27)
project(Deblur)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(Deblur src/main.cpp src/image.cpp)
