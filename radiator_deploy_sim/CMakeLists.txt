cmake_minimum_required(VERSION 3.10)
project(radiator_deploy_sim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Eigen ===
find_package(Eigen3 REQUIRED)

# === NLopt (manual config since C++ header not installed via apt) ===
set(NLOPT_INCLUDE_DIR "/usr/local/include")
set(NLOPT_LIBRARY "/usr/local/lib/libnlopt.so")

# === Source files ===
file(GLOB SOURCES "src/*.cpp")

add_executable(radiator_deploy_sim ${SOURCES})

# === Includes and Linking ===
target_include_directories(radiator_deploy_sim
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${NLOPT_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(radiator_deploy_sim
    PRIVATE
    Eigen3::Eigen
    ${NLOPT_LIBRARY}
)

