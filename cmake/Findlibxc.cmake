include(FindPackageHandleStandardArgs)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)

if (DEFINED ENV{LIBXC_DIR})
  set(LIBXC_DIR "$ENV{LIBXC_DIR}")
endif()

find_library(LIBXC_LIBRARY
  NAMES xc
  PATH_SUFFIXES lib
  HINTS ${LIBXC_DIR}
  )

find_path(LIBXC_INCLUDE_DIRS
  NAMES xc.h
  PATH_SUFFIXES include
  HINTS ${LIBXC_DIR}
  )

find_package_handle_standard_args(LIBXC
  REQUIRED_VARS LIBXC_LIBRARY LIBXC_INCLUDE_DIRS
  )

mark_as_advanced(LIBXC_LIBRARY LIBXC_INCLUDE_DIRS)
if (NOT LIBXC_FOUND)
  set(LIBXC_DIR "${LIBXC_DIR}" CACHE STRING "Directory containing the libxc library.")
endif()
