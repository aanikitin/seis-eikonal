#!/bin/sh
# clean build dirs
rm -rf bin
rm -rf build
rm -rf lib
rm -rf Testing
rm -rf test/EIKONAL_EX1/bin
rm -rf CMakeLists.txt.user
cd wrappers/fortran
sh clean.sh
cd ../matlab
rm -rf *.mex
