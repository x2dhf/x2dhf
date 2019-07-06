#!/bin/bash

# Installation directory
export target=$(pwd) # this installs x2dhf under bin/ in the present directory
export X2DHF=$(echo $target | sed -e 's/\//\\\//g' )
export INCLUDE4LXC=$target/include
export INCLUDES="$INCLUDE4LXC "

# Fortran compiler to use
export FC=gfortran
# Compiler flags to use
export FFLAGS="-O3 -Wall -ffree-line-length-none -I $INCLUDE4LXC "
#export FFLAGS=" -O -g  -fbounds-check -I $INCLUDE4LXC "

# Number of parallel processes in build
export npar=10

# No changes should be necessary hereafter.
if [[ ! -d objdir ]]; then
    mkdir objdir
fi
cd objdir
cmake .. \
      -DCMAKE_INSTALL_PREFIX=${target} \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-DNDEBUG" 

sed -i -e "s/\-o x2dhf *$/\-o x2dhf \-L $X2DHF\/lib \-lxcf90 \-lxc \n/" src/CMakeFiles/x2dhf.dir/link.txt
#tail src/CMakeFiles/x2dhf.dir/link.txt
make -j ${npar} install VERBOSE=1
cd ..


