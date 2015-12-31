cmake_minimum_required(VERSION 3.3)
project(satori)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_definitions(-DDEBUG)

set(SOURCE_FILES main.cpp)
add_executable(solution ${SOURCE_FILES})