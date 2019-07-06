#!/bin/bash

export X2DHF=$HOME/github/x2dhf
export INCLUDE4LXC=$X2DHF/include

# Fortran compiler to use
export FC=gfortran
# Compiler flags to use
export FFLAGS="-O2 -Wall -I $INCLUDE4LXC -ffree-line-length-none"

# Installation directory
export target=$(pwd) # this installs x2dhf under bin/ in the present directory
#export target=${HOME}  # this installs x2dhf in ${HOME}/bin

# Number of parallel processes in build
export npar=4

# No changes should be necessary hereafter.
if [[ ! -d objdir ]]; then
    mkdir objdir
fi
cd objdir
cmake .. \
      -DCMAKE_INSTALL_PREFIX=${target} \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-DNDEBUG"
make -j ${npar} install VERBOSE=1
cd ..
