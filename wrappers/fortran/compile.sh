#!/bin/sh
gfortran -fdiagnostics-color=always -O3 -Winline -Wall -Wextra -fmessage-length=0 -pedantic -fopenmp -c module.f95
gfortran -fdiagnostics-color=always -O3 -Winline -Wall -Wextra -fmessage-length=0 -pedantic -fopenmp fort.f95 module.o ../../lib/libopenst.a -lm
