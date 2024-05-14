#!/usr/bin/env bash

# Set names/directories
target_dir=`pwd`"/hypre"
hypre_tarname="hypre-2.31.0.tar.gz"
hypre_dirname="hypre-2.31.0"
build_dir="build"

# Do compilation etc. in build directory
mkdir -p ${build_dir}
cd ${build_dir}

# Extract
if [ ! -d ${hypre_dirname} ]; then
    tar -xzf ${hypre_tarname}
fi

# Configure
cd ${hypre_dirname}/src

./configure \
    --with-openmp\
    --without-MPI\
    --with-print-errors\
    --prefix=${target_dir}

make -j
make install
