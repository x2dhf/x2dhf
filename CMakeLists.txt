# CMake project file for x2dhf
cmake_minimum_required (VERSION 3.5)
project (x2dhf Fortran)

# Include CMake macros
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Optional support for pthreads
option(USE_PTHREAD "Enable pthreads pool option" OFF)
option(USE_TPOOL   "Enable pthreads pool option" OFF)
# Support for OpenMP
option(USE_OPENMP "Include support for OpenMP" ON)
# Support for Libxc
option(USE_LIBXC "Include support for Libxc" OFF)
# Support for BLAS
option(USE_BLAS "Use external blas library" OFF)

# Detect OpenMP
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Detect pthreads
if(USE_PTHREADS OR USE_TPOOL)
  enable_language(C)

  # Test if pthreads is available on the platform
  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
  try_compile(ENABLE_PTHREADS
    SOURCES ${CMAKE_SOURCE_DIR}/src/c/mcsor_single_colour_pthread.c)
  if(NOT ENABLE_PTHREADS)
    message(FATAL_ERROR "Pthreads support not available on your platform")
  endif()
endif()

# Detect Libxc
if(USE_LIBXC)
  find_package(Libxc REQUIRED)
  message("Libxc found ${LIBXC_FOUND}: ${LIBXC_LIBRARY}")
endif()

# Detect BLAS
if(USE_BLAS)
  find_package(BLAS REQUIRED)
else()
  set(BLAS_FOUND 0)
endif()

# Compiler specific flags
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fmax-array-constructor=90000")
endif()

# Recurse into src
add_subdirectory(src)
