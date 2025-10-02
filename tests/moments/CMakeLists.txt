cmake_minimum_required(VERSION 3.20.1)
project(test-moments)

add_executable(test-moments test_moments.cpp)

target_link_libraries(test-moments PRIVATE HighFive)
target_include_directories(test-moments PRIVATE ${CMAKE_SOURCE_DIR}/src)

