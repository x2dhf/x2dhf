#!/bin/bash

install_prefix=$(pwd)

if [[ ! -d libxc ]]
then
    git clone https://gitlab.com/libxc/libxc.git 2>/dev/null
fi

cd libxc
git pull 2>/dev/null  

export npar=6

autoconf >/dev/null
./configure --prefix=$install_prefix/libxc

# ENABLE_FORTRAN ON
sed -ie 's/ENABLE_FORTRAN "Build Fortran 90 interface" OFF/ENABLE_FORTRAN "Build Fortran 90 interface" ON/' CMakeLists.txt

cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=$install_prefix
cd objdir && make -j $npar
#make test
make install






