cmake_minimum_required(VERSION 3.10)
project(GMRES_Solver Fortran)

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fopenmp -march=native -funroll-loops -g")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fopenmp -march=native -funroll-loops -g -fno-omit-frame-pointer")

# Configure optimization flags for GNU
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Active optimization Flags: ${CMAKE_Fortran_FLAGS}")
endif()

# Search OpenMP
find_package(OpenMP REQUIRED)
# BLAS and LAPACK not needed
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)

#Add GMRES Library
add_library(gmres_lib src/gmres_mgsr.f90 src/gmres_hh.f90 src/interfaces.f90)
#Add Poisson 2D Libray
add_library(poisson_lib src/problems/poisson.f90)
#ASS Hilbert's Library
add_library(hilbert_lib src/problems/hilbert.f90)
#Add Preconds Library
add_library(precond_lib src/preconds/chebyshev.f90)
#Add Utils Library
add_library(utils_lib src/utils/utils.f90)
#Dense Matrix Poisson Test
add_executable(test_dp tests/test_poisson.f90)
target_link_libraries(test_dp gmres_lib poisson_lib)
set_target_properties(test_dp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test_dp PRIVATE ${CMAKE_BINARY_DIR})

#Dense Matrix Hilbert Test
add_executable(test_hilbert tests/test_hilbert.f90)
target_link_libraries(test_hilbert gmres_lib hilbert_lib)
set_target_properties(test_hilbert PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test_hilbert PRIVATE ${CMAKE_BINARY_DIR})

#Matrix Free Poisson Test
add_executable(test_mfp tests/test_poisson_mf.f90)
target_link_libraries(test_mfp gmres_lib poisson_lib precond_lib)
set_target_properties(test_mfp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test_mfp PRIVATE ${CMAKE_BINARY_DIR})

#Matrix Free Poisson Test
add_executable(test1 tests/test1.f90)
target_link_libraries(test1 gmres_lib poisson_lib precond_lib utils_lib)
set_target_properties(test1 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test1 PRIVATE ${CMAKE_BINARY_DIR})

# Link OpenMP Libraries
target_link_libraries(test_mfp OpenMP::OpenMP_Fortran)
#target_link_libraries(test_hilbert OpenMP::OpenMP_Fortran)
#target_link_libraries(test_poisson OpenMP::OpenMP_Fortran)
