#!/bin/sh
set -e
RES_DIR="RES_EIKONAL_EX1_`date +%Y-%m-%d_%H-%M-%S_%N`"

mkdir ${RES_DIR}
hostname > ${RES_DIR}/cpuinfo.txt
echo >> ${RES_DIR}/cpuinfo.txt
cat /proc/cpuinfo >> ${RES_DIR}/cpuinfo.txt

BIN_DIR=bin

DIM1=400
DIM2=400
DIM3=400

NTH1=8
NTH2=1
NTH3=8

BI1=1
BI2=1
BI3=20

BJ1=1
BJ2=1
BJ3=20

BK1=100
BK2=100
BK3=400

NTEST1=1
NTEST2=1
NTEST3=1

. ./omp_env.sh

./test_serial.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_fsm3d_fsm_serial_v1 ${DIM1} ${DIM2} ${DIM3} ${NTEST1} ${NTEST2} ${NTEST3}
./test_serial.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_lsm3d_lsm_serial_v1 ${DIM1} ${DIM2} ${DIM3} ${NTEST1} ${NTEST2} ${NTEST3}

./test_block.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_fsm3d_bfsm_openmp_v1 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${BI1} ${BI2} ${BI3} ${BJ1} ${BJ2} ${BJ3} ${BK1} ${BK2} ${BK3} ${NTEST1} ${NTEST2} ${NTEST3}
./test_block.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_lsm3d_blsm_openmp_v1 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${BI1} ${BI2} ${BI3} ${BJ1} ${BJ2} ${BJ3} ${BK1} ${BK2} ${BK3} ${NTEST1} ${NTEST2} ${NTEST3}

./test_block.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_fsm3d_bfsm_openmp_v2 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${BI1} ${BI2} ${BI3} ${BJ1} ${BJ2} ${BJ3} ${BK1} ${BK2} ${BK3} ${NTEST1} ${NTEST2} ${NTEST3}
./test_block.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_lsm3d_blsm_openmp_v2 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${BI1} ${BI2} ${BI3} ${BJ1} ${BJ2} ${BJ3} ${BK1} ${BK2} ${BK3} ${NTEST1} ${NTEST2} ${NTEST3}

./test_block.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_fsm3d_bfsm_openmp_v3 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${BI1} ${BI2} ${BI3} ${BJ1} ${BJ2} ${BJ3} ${BK1} ${BK2} ${BK3} ${NTEST1} ${NTEST2} ${NTEST3}
./test_block.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_lsm3d_blsm_openmp_v3 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${BI1} ${BI2} ${BI3} ${BJ1} ${BJ2} ${BJ3} ${BK1} ${BK2} ${BK3} ${NTEST1} ${NTEST2} ${NTEST3}

./test_serial.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_fsm3d_dfsm_serial_v1 ${DIM1} ${DIM2} ${DIM3} ${NTEST1} ${NTEST2} ${NTEST3}
./test_serial.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_lsm3d_dlsm_serial_v1 ${DIM1} ${DIM2} ${DIM3} ${NTEST1} ${NTEST2} ${NTEST3}

./test_openmp.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_fsm3d_dfsm_openmp_v1 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${NTEST1} ${NTEST2} ${NTEST3}
./test_openmp.sh ${RES_DIR} ${BIN_DIR}/EIKONAL_EX1_lsm3d_dlsm_openmp_v1 ${DIM1} ${DIM2} ${DIM3} ${NTH1} ${NTH2} ${NTH3} ${NTEST1} ${NTEST2} ${NTEST3}
