# CMake project file for x2dhf
cmake_minimum_required (VERSION 2.6)
project (x2dhf Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# If flags haven't been defined, set defaults
if (NOT CMAKE_Fortran_FLAGS_DEBUG)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  else()
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("No debug flags are known, we just try -O0 -g...")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif()
endif()

if (NOT CMAKE_Fortran_FLAGS_RELEASE)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -funroll-all-loops")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  else()
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("No optimized flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif()
endif()

# Recurse into src
add_subdirectory(src)
