#!/bin/bash

build_dir=build
rm -rf $build_dir
mkdir $build_dir
cd $build_dir

# built from source
# export PETSC_DIR=/Users/yuechengyu/Desktop/Stanford/Cardiac/petsc-balay/petsc
# export PETSC_ARCH=arch-darwin-c-opt

# export DYLD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH
# export LIBRARY_PATH=$PETSC_DIR/lib:$LIBRARY_PATH

cmake \
-DSV_USE_PETSC=ON \
-DCMAKE_BUILD_TYPE=Debug \
.. 

make -j4 CC=gcc-14 FC=gfortran-13 CXX=g++-14
