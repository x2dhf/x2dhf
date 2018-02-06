#!/bin/bash

# Fortran compiler to use
export FC=gfortran
# Compiler flags to use
export FFLAGS="-O2 -Wall"
# Installation directory (x2dhf will be installed under bin/)
export target=$(pwd)

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

