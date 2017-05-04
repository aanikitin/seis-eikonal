#!/bin/sh
set -e
# clean build dirs
rm -rf bin
rm -rf build
rm -rf lib
mkdir build
cd build
# compile with LTO enabled (requires up-to-date gcc, gcc-ar and gcc-ranlib)
# tested with gcc 5.4.0
# switch to OFF if unsuccessful
# build double precision version (default)
cmake -DGCC_LTO:BOOL=ON ..
# or to build single precision version, use:
# cmake -DGCC_LTO:BOOL=ON -DOPENST_FLOAT_PRECISION_SINGLE= ..
make
ctest ..
