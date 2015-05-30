cmake_minimum_required(VERSION 2.6)
project(freq_seq)

set(CMAKE_CXX_FLAGS  "-std=c++11")

include_directories("${CMAKE_SOURCE_DIR}/third-party/bamtools/src")
include_directories("${CMAKE_SOURCE_DIR}/third-party/")
include_directories("./")

add_executable(freq_seq freq_seq.cc)
target_link_libraries(freq_seq bamtools bamtools-utils)


add_executable(get_raf get_raf.cc)
target_link_libraries(get_raf bamtools bamtools-utils)

