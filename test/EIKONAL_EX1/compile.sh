#!/bin/sh
echo "CC=${CC}"
mkdir bin
# compile all implementations of FSM3D
for f in `find ../../src/openst/eikonal/fsm/fsm3d_computepartial -type f -name "*.c" -printf "%f\n"`
do
for f2 in `find ../../src/openst/eikonal/fsm/fsm3d_blockserial -type f -name "*.c" -printf "%f\n"`
do
BIN_NAME="EIKONAL_EX1_${f%%.*}_${f2%%.*}"
printf "\n=== BUILDING: ${BIN_NAME} ===\n"
${CC} -DTEST_FSM -O3 -Wall -Wextra -fopenmp -o bin/${BIN_NAME} src/main.c -L ../../lib/${f}/${f2} -I ../../include -Wl,-Bstatic -lopenst -Wl,-Bdynamic -lm
done
done

# compile all implementations of LSM3D
for f in `find ../../src/openst/eikonal/lsm/lsm3d_computepartial -type f -name "*.c" -printf "%f\n"`
do
for f2 in `find ../../src/openst/eikonal/lsm/lsm3d_blockserial -type f -name "*.c" -printf "%f\n"`
do
BIN_NAME="EIKONAL_EX1_${f%%.*}_${f2%%.*}"
printf "\n=== BUILDING: ${BIN_NAME} ===\n"
${CC} -O3 -Wall -Wextra -fopenmp -o bin/${BIN_NAME} src/main.c -L ../../lib/${f}/${f2} -I ../../include -Wl,-Bstatic -lopenst -Wl,-Bdynamic -lm
done
done

printf "\n=== FINISHED ===\n"
