cmake_minimum_required(VERSION 3.20.1)
project(test-ampere)

add_executable(test-ampere test_ampere.cpp)

target_link_libraries(test-ampere PRIVATE HighFive)
target_include_directories(test-ampere PRIVATE ${CMAKE_SOURCE_DIR}/src)

