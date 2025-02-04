# CMake project file for x2dhf
cmake_minimum_required (VERSION 3.5)
project (x2dhf Fortran)

# Optional support for pthreads
option(USE_PTHREADS "Include support for pthreads" OFF)
# Support for OpenMP
option(USE_OPENMP "Include support for OpenMP" ON)

# Parse options
if(USE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()
if(USE_PTHREADS)
  enable_language(C)
endif()

# Compiler specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fmax-array-constructor=90000")
endif()

# Recurse into src
add_subdirectory(src)
