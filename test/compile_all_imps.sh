#!/bin/sh
cd ..

# clean build dirs
rm -rf build
rm -rf lib
mkdir build
cd build

# compile all implementations of FSM3D
for f in `find ../src/openst/eikonal/fsm/fsm3d_computepartial -type f -name "*.c" -printf "%f\n"`
do
for f2 in `find ../src/openst/eikonal/fsm/fsm3d_blockserial -type f -name "*.c" -printf "%f\n"`
do
printf "\n\n=== BUILDING: $f/${f2} ===\n"
cmake -DGCC_LTO:BOOL=ON -DFSM3D_IMP="$f" -DFSM3D_BLOCKSERIAL_IMP="${f2}" -DOPENST_PATH_LIB_SUFFIX="$f/${f2}" ..
make
done
done
printf "\n=== FINISHED ===\n"

# compile all implementations of LSM3D
for f in `find ../src/openst/eikonal/lsm/lsm3d_computepartial -type f -name "*.c" -printf "%f\n"`
do
for f2 in `find ../src/openst/eikonal/lsm/lsm3d_blockserial -type f -name "*.c" -printf "%f\n"`
do
printf "\n\n=== BUILDING: $f/${f2} ===\n"
cmake -DGCC_LTO:BOOL=ON -DLSM3D_IMP="$f" -DLSM3D_BLOCKSERIAL_IMP="${f2}" -DOPENST_PATH_LIB_SUFFIX="$f/${f2}" ..
make
done
done
printf "\n=== FINISHED ===\n"
