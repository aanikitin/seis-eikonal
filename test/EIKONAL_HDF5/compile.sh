#!/bin/sh

# compile all implementations of FSM3D
for f in `find ../../src/openst/eikonal/fsm/fsm3d_imp -type f -name "*.c" -printf "%f\n"`
do
gcc -O3 -Wall -Wextra -fopenmp -o EIKONAL_HDF5_${f} main.c fio_hdf5.c -L ../../lib/${f} -I /usr/include/hdf5/serial/ -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -I ../../include -lhdf5_hl -lhdf5 -l:libopenst.a -lm  
done
printf "\n=== FINISHED ===\n"

# compile all implementations of LSM3D
for f in `find ../../src/openst/eikonal/lsm/lsm3d_imp -type f -name "*.c" -printf "%f\n"`
do
gcc -O3 -Wall -Wextra -DLSM3D_IMP -fopenmp -o EIKONAL_HDF5_${f} main.c fio_hdf5.c -I /usr/include/hdf5/serial/ -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -L ../../lib/${f} -I ../../include -lhdf5_hl -lhdf5 -l:libopenst.a -lm 
done
