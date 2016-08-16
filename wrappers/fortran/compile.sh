#!/bin/sh
export FFLAGS="-O3 -std=f2003 -freal-4-real-8 -fdiagnostics-color=always -Winline -Wall -Wextra -fmessage-length=0 -pedantic -fopenmp"
gfortran ${FFLAGS} -c OpenST.f90
gfortran -o BRT3D_EX1.exe ${FFLAGS} BRT3D_EX1.f90 OpenST.o ../../lib/libopenst.a -lm
gfortran -o EIKONAL_EX1.exe ${FFLAGS} EIKONAL_EX1.f90 OpenST.o ../../lib/libopenst.a -lm
