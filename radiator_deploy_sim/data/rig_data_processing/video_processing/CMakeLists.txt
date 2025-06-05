cmake_minimum_required(VERSION 3.10)
project(process_video)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===== Adjust to wherever your OpenCV build is located =====
set(OpenCV_DIR /home/jonathanvollrath/libraries/opencv/build)

# ===== Find OpenCV package =====
find_package(OpenCV REQUIRED)

# ===== Source files =====
file(GLOB SOURCES "src/*.cpp")

# ===== Executable =====
add_executable(process_video ${SOURCES})

# ===== Link OpenCV =====
target_link_libraries(process_video ${OpenCV_LIBS})

target_compile_definitions(process_video PRIVATE GENERATE_MARKERS)