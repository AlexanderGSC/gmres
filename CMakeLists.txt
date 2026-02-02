cmake_minimum_required(VERSION 3.10)
project(GMRES_Solver Fortran)

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fopenmp -march=native -funroll-loops -g")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fopenmp -march=native -funroll-loops -g -fno-omit-frame-pointer")
# Search OpenMP
find_package(OpenMP REQUIRED)
# Configure optimization flags for GNU
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Flags de optimización activos: ${CMAKE_Fortran_FLAGS}")
endif()

# BLAS and LAPACK not needed
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)

#Add GMRES Library
add_library(gmres_lib src/gmres_mgsr.f90 src/gmres_hh.f90 src/matrix_utils.f90 src/precond.f90)
#Add GMRES Matrix Free Library
add_library(gmresmf_lib src/gmres_mf.f90 src/matrix_utils.f90)

#Poisson
add_executable(test_poisson tests/test_poisson.f90)
target_link_libraries(test_poisson gmres_lib)
set_target_properties(test_poisson PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test_poisson PRIVATE ${CMAKE_BINARY_DIR})

#Hilbert
add_executable(test_hilbert tests/test_hilbert.f90)
target_link_libraries(test_hilbert gmres_lib)
set_target_properties(test_hilbert PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test_hilbert PRIVATE ${CMAKE_BINARY_DIR})

#Matrix Free
add_executable(test_mf tests/test_poisson_mf.f90)
target_link_libraries(test_mf gmresmf_lib)
set_target_properties(test_hilbert PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(test_hilbert PRIVATE ${CMAKE_BINARY_DIR})

# Link OpenMP Libraries
target_link_libraries(test_mf OpenMP::OpenMP_Fortran)
target_link_libraries(test_hilbert OpenMP::OpenMP_Fortran)
target_link_libraries(test_poisson OpenMP::OpenMP_Fortran)
