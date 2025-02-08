# CMake project file for x2dhf
cmake_minimum_required (VERSION 3.5)
project (x2dhf Fortran)
# Include CMake macros
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### Input options
# Optional support for pthreads
option(USE_PTHREAD "Enable pthreads pool option" OFF)
option(USE_TPOOL   "Enable pthreads pool option" OFF)
# Support for OpenMP
option(USE_OPENMP "Include support for OpenMP" ON)
# Support for Libxc
option(USE_LIBXC "Include support for Libxc (use environment variable LIBXC_DIR to guide detection)" OFF)
# Support for BLAS
option(USE_BLAS "Use external blas library" ON)

### Parse the options
# Check for conflicting options
if((${USE_OPENMP} AND ${USE_TPOOL}) OR (${USE_OPENMP} AND ${USE_PTHREAD}))
  message(FATAL_ERROR "OpenMP and pthreads support is mutually exclusive")
endif()
if(${USE_TPOOL} AND ${USE_PTHREAD})
  message(FATAL_ERROR "tpool and pthreads support is mutually exclusive")
endif()

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
  # Restore try compile target
  set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
endif()

# Detect Libxc
if(USE_LIBXC)
  find_package(LIBXC REQUIRED)
  if(NOT ${LIBXC_FOUND})
    message(FATAL_ERROR "Did not find Libxc")
  else()
    message("Libxc include directory ${LIBXC_INCLUDE_DIRS}")
  endif()
endif()

# Detect BLAS
if(USE_BLAS)
  find_package(BLAS REQUIRED)
else()
  set(BLAS_FOUND 0)
endif()

# Set the default build type if one was not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Recurse into src
add_subdirectory(src)
