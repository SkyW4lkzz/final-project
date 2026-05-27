cmake_minimum_required(VERSION 3.16)

project(lbm_final_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(
    lbm_solver
    src/main.cpp
    src/lbm.cpp
    src/poiseuille.cpp
    src/cylinder.cpp
)

target_include_directories(lbm_solver PRIVATE include)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lbm_solver PRIVATE OpenMP::OpenMP_CXX)
elseif(APPLE)
    foreach(LIBOMP_PREFIX /opt/homebrew/opt/libomp /usr/local/opt/libomp)
        if(EXISTS "${LIBOMP_PREFIX}")
            target_compile_options(lbm_solver PRIVATE -Xpreprocessor -fopenmp)
            target_include_directories(lbm_solver PRIVATE "${LIBOMP_PREFIX}/include")
            target_link_libraries(lbm_solver PRIVATE "${LIBOMP_PREFIX}/lib/libomp.dylib")
            message(STATUS "Using Homebrew libomp: ${LIBOMP_PREFIX}")
            break()
        endif()
    endforeach()
endif()

if(MSVC)
    target_compile_options(lbm_solver PRIVATE /W4)
else()
    target_compile_options(lbm_solver PRIVATE -Wall -Wextra -Wpedantic)
endif()
