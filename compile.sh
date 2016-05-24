#!/bin/sh
# clean build dirs
rm -rf bin
rm -rf build
rm -rf lib
mkdir build
# compile with default settings
cd build
cmake ..
make
