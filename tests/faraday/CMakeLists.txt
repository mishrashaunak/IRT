cmake_minimum_required(VERSION 3.20.1)
project(test-faraday)

add_executable(test-faraday test_faraday.cpp)

target_link_libraries(test-faraday PRIVATE HighFive)
target_include_directories(test-faraday PRIVATE ${CMAKE_SOURCE_DIR}/src)

