cmake_minimum_required(VERSION 3.10)
project(GMRES_Solver Fortran)

# Buscar OpenMP
find_package(OpenMP REQUIRED)

# Buscar BLAS y LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Indicar dónde están los archivos de módulos y fuentes
set(SOURCES
    src/givens_mod.f90 
    src/arnoldi_mod.f90
    src/householder_mod.f90
    src/gmres_mod.f90 
    src/main.f90
)

# Configurar flags de optimización para GNU
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -funroll-loops")
    message(STATUS "Flags de optimización activos: ${CMAKE_Fortran_FLAGS}")
endif()

# Si usas C++, haz lo mismo para CXX
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")
endif()

# Crear el ejecutable
add_executable(gmres_exe ${SOURCES})

# Forzar que los archivos .mod se guarden en la carpeta build
set_target_properties(gmres_exe PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})
target_include_directories(gmres_exe PRIVATE ${CMAKE_BINARY_DIR})

# Enlazar las librerías
target_link_libraries(gmres_exe OpenMP::OpenMP_Fortran ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
