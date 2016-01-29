#!/bin/sh
# clean build dirs
rm -rf build
mkdir build
cd build
# optimization flags to use during compilation
C_FLAGS_OPT=-O3
# compile all implementations of FSM3D
for f in `find ../src/fsm -type f -name "*.c" -printf "%f\n"`
do
printf "\n=== BUILDING: $f ===\n"
cmake -DIMP_FSM3D="$f" -DC_FLAGS_OPT=$C_FLAGS_OPT ..
make
done
printf "\n=== FINISHED ===\n"
