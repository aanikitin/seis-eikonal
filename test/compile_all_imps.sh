#!/bin/sh
cd ..

# clean build dirs
rm -rf build
rm -rf lib
mkdir build
cd build

# compile all implementations of FSM3D
for f in `find ../src/openst/eikonal/fsm/fsm3d_imp -type f -name "*.c" -printf "%f\n"`
do
printf "\n\n=== BUILDING: $f ===\n"
cmake -DFSM3D_IMP="$f" -DOPENST_PATH_LIB_SUFFIX="$f" ..
make
done
printf "\n=== FINISHED ===\n"

# compile all implementations of LSM3D
for f in `find ../src/openst/eikonal/lsm/lsm3d_imp -type f -name "*.c" -printf "%f\n"`
do
printf "\n\n=== BUILDING: $f ===\n"
cmake -DLSM3D_IMP="$f" -DOPENST_PATH_LIB_SUFFIX="$f" ..
make
done
printf "\n=== FINISHED ===\n"
