#!/bin/bash

set -e
MAKE="make --jobs=$NUM_THREADS --keep-going"
MAKE_TEST=ctest

mkdir -p $BUILD_DIR
cd $BUILD_DIR
CMAKE_BUILD_ARGS="$CMAKE_BUILD_ARGS -DVTK_DIR:PATH=$VTK_DIR -DBUILD_TESTING:BOOL=ON"
echo CMAKE_BUILD_ARGS: $CMAKE_BUILD_ARGS
cmake $CMAKE_BUILD_ARGS ../
$MAKE
$MAKE_TEST

