cmake_minimum_required(VERSION $ENV{CMVERSION})
PROJECT( x2dhf Fortran )

SET(CMAKE_VERBOSE_MAKEFILE $ENV{VERBOSE})
# by default CMAKE_INSTALL_PREFIX=/usr/local/
#set(CMAKE_INSTALL_PREFIX $ENV{HOME}/git/x2dhf/)
set(CMAKE_INSTALL_PREFIX $ENV{INSTALLPREFIX})

set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90 F90)

# ld on macOS does not support --start|end-group options
#link_libraries("$ENV{LDOPT1}$ENV{LDOPT2}$ENV{LDOPT3}$ENV{LDOPT4}-Wl,--start-group \
#$ENV{LIB1} \
#$ENV{LIB2} \
#$ENV{LIB3} \
#-Wl,--end-group")
link_libraries("$ENV{LDOPT1}$ENV{LDOPT2}$ENV{LDOPT3}$ENV{LDOPT4}$ENV{LIB1}$ENV{LIB2}$ENV{LIB3}")

include (${CMAKE_CURRENT_SOURCE_DIR}/$ENV{SETTXT})

ADD_EXECUTABLE(x2dhf ${x2dhf_SRCS})

# Install
install (TARGETS x2dhf DESTINATION build )
