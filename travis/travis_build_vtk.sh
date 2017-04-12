#!/bin/bash

set -e
MAKE="make --jobs=$NUM_THREADS --keep-going"

cmake_arg_str=" -DBUILD_TESTING=1 -DBUILD_SHARED_LIBS=1  -DBUILD_EXAMPLES=0 "
vtk_repo_str=""
if [ "$VTK_VERSION" == "6.0" ]; then
    vtk_repo_str=" --branch v6.0.0"
elif [ "$VTK_VERSION" == "6.2" ]; then
    vtk_repo_str=" --branch v6.2.0"
elif [ "$VTK_VERSION" == "7.0" ]; then
    vtk_repo_str=" --branch v7.0.0"
fi
if [ -d $VTK_SOURCE_DIR ]; then
    echo $VTK_SOURCE_DIR exists
    if [ ! -f $VTK_SOURCE_DIR/CMakeLists.txt ]; then
        echo $VTK_SOURCE_DIR does not contain CMakeList.txt
        rm -rf $VTK_SOURCE_DIR
    fi
fi
if [ ! -d "$VTK_SOURCE_DIR" ]; then
    echo "Using git to clone source"
    git clone $vtk_repo_str $VTK_SOURCE_DIR
fi
mkdir -p $VTK_DIR
cd $VTK_DIR
cmake $cmake_arg_str $VTK_SOURCE_DIR
$MAKE
