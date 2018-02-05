# CMake project file for x2dhf
cmake_minimum_required (VERSION 2.6)
project (x2dhf)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the x2dhf directory."
      FORCE)

# FFLAGS depend on the compiler
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -funroll-all-loops")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall")
elseif (Fortran_COMPILER_NAME MATCHES "Intel")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
else()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif()

# Recurse into src
add_subdirectory(src)
