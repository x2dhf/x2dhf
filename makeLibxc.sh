#!/bin/bash

install_prefix=$(pwd)

if [[ ! -d libxc ]]
then
    git clone https://gitlab.com/libxc/libxc.git 2>/dev/null
fi


cd libxc
git pull 2>/dev/null  

./configure --prefix=$install_prefix/libxc

cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=$install_prefix
cd objdir && make
#make test
make install






