cmake_minimum_required(VERSION 3.16)
project(CS171-hw5 LANGUAGES C CXX)

add_compile_options(-g -msse -msse2 -msse3 -mssse3 -msse4 -msse4.1 -msse4.2 -O3)
set(CMAKE_BUILD_TYPE Release)

option(USE_FILESYSTEM "use std::filesystem" off)

# glfw
option(GLFW_BUILD_TESTS off)
option(GLFW_BUILD_EXAMPLES off)
option(GLFW_BUILD_DOCS off)
add_subdirectory(libs/glfw)
include_directories(libs/glfw/include)
list(APPEND LIBRARIES glfw ${glfw_LIBRARIES})

# stb
add_subdirectory(libs/stb)
include_directories(libs/stb)

# glad
set(GLAD_INCLUDE libs/glad/include)
include_directories(${GLAD_INCLUDE})

# glm
add_subdirectory(libs/glm)
list(APPEND LIBRARIES glm::glm)

# add filesystem
if (USE_FILESYSTEM)
    list(APPEND DEFINES -DUSE_FILESYSTEM)
    if (NOT MSVC)
        list(APPEND LIBRARIES stdc++fs)
    endif()
endif()

# set source files and add executable
file(GLOB SRC_FILES src/*.cpp)

# openmp
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    list(APPEND LIBRARIES OpenMP::OpenMP_CXX)
endif()

# CUDA
find_package(CUDA)
if (CUDA_FOUND)
    message("CUDA found")
    enable_language(CUDA)
    list(APPEND LIBRARIES ${CUDA_LIBRARIES} cuda)
    set(libcudacxx_DIR /opt/cuda/targets/x86_64-linux/lib/cmake/libcudacxx)
    set(CUB_DIR /opt/cuda/targets/x86_64-linux/lib/cmake/cub/)
    set(Thrust_DIR /opt/cuda/targets/x86_64-linux/lib/cmake/thrust)
    find_package(Thrust REQUIRED CONFIG)
    thrust_create_target(Thrust)
    list(APPEND LIBRARIES Thrust)
else()
    message("CUDA not found")
endif()

# PCISPH
add_subdirectory(pcisph)
include_directories(pcisph)

add_executable(main_cpu ${SRC_FILES} ${INTERNAL_SRC_FILES} src/glad.c)
target_include_directories(main_cpu PRIVATE include)
target_compile_features(main_cpu PRIVATE cxx_std_20)
target_compile_definitions(main_cpu PRIVATE -DUSE_CPU ${DEFINES})
target_link_libraries(main_cpu PRIVATE PCISPH_CPU ${LIBRARIES})

add_executable(main_surface ${SRC_FILES} ${INTERNAL_SRC_FILES} src/glad.c)
target_include_directories(main_surface PRIVATE include)
target_compile_features(main_surface PRIVATE cxx_std_20)
target_compile_definitions(main_surface PRIVATE -DUSE_CPU -DSURFACE_EXTRACTION ${DEFINES})
target_link_libraries(main_surface PRIVATE PCISPH_CPU ${LIBRARIES})

if (CUDA_FOUND)
    add_executable(main_cuda ${SRC_FILES} ${INTERNAL_SRC_FILES} src/glad.c)
    target_include_directories(main_cuda PRIVATE include)
    target_compile_features(main_cuda PRIVATE cxx_std_20)
    target_compile_definitions(main_cuda PRIVATE ${DEFINES})
    target_link_libraries(main_cuda PRIVATE PCISPH_CUDA ${LIBRARIES})
endif()
