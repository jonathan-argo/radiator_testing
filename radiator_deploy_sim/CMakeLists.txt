cmake_minimum_required(VERSION 3.10)
project(radiator_deploy_sim)

find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/include/)

file(GLOB SOURCES "src/*.cpp")

add_executable(radiator_deploy_sim ${SOURCES})

target_link_libraries(radiator_deploy_sim Eigen3::Eigen)