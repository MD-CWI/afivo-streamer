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

# A bit of a hack to add the -fPIC flag without compiling a shared library
# This flag is necessary when using the library with f2py
sed -i 's/$(FC_COMPILE_FLAGS)/-fPIC $(FC_COMPILE_FLAGS)/g' config/Makefile.config
sed -i 's/$(C_COMPILE_FLAGS)/-fPIC $(C_COMPILE_FLAGS)/g' config/Makefile.config
sed -i 's/$(CXX_COMPILE_FLAGS)/-fPIC $(CXX_COMPILE_FLAGS)/g' config/Makefile.config

make -j
make install
