#!/bin/sh
set -e
root_dir=`pwd`
# clean build dirs
rm -rf bin
rm -rf build
rm -rf lib
rm -rf CMakeLists.txt.user
# clean matlab dir
cd ${root_dir}/wrappers/matlab
rm -rf *.mex
rm -rf *.mexw64
rm -rf *.asv
