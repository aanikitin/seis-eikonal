#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "fsm.h"

#define BSIZE_DEF 8u
size_t BSIZE_I;
size_t BSIZE_J;
size_t BSIZE_K;

void ex1_init(double *F, size_t NI, size_t NJ, size_t NK){
    size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                F[M_MEM_ADR_3D(i,j,k,NI,NJ,NK)] = 1.0;
            }
        }
    }
}

void ex1_check(double *U, double H, size_t NI, size_t NJ, size_t NK,
              size_t SRCI, size_t SRCJ, size_t SRCK,
               double *L1, double *L2, double *Linf){
    size_t i, j, k, NN;
    double di, dj, dk, dist;
    double diff, l1, l2, linf;
    l1 = 0.0;
    l2 = 0.0;
    linf = DBL_MIN;
    NN = NI * NJ * NK;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                di = ((double)SRCI - (double)i) * H;
                dj = ((double)SRCJ - (double)j) * H;
                dk = ((double)SRCK - (double)k) * H;
                dist = sqrt(di * di + dj * dj + dk * dk);
                diff = U[M_MEM_ADR_3D(i,j,k,NI,NJ,NK)] - dist;
                l1 += fabs(diff);
                l2 += (diff * diff);
                linf = fmax(linf,fabs(diff));
            }
        }
    }
    *L1 = l1/NN;
    *L2 = l2/NN;
    *Linf = linf;
}

void build_info_print(){
#ifdef BUILD_INFO_CMAKE_C_COMPILER
    printf("BUILD_INFO_CMAKE_C_COMPILER: %s\n", BUILD_INFO_CMAKE_C_COMPILER);
#else
    printf("BUILD_INFO_CMAKE_C_COMPILER: [UNDEFINED]\n");
#endif
#ifdef BUILD_INFO_CMAKE_C_COMPILER_ID
    printf("BUILD_INFO_CMAKE_C_COMPILER_ID: %s\n", BUILD_INFO_CMAKE_C_COMPILER_ID);
#else
    printf("BUILD_INFO_CMAKE_C_COMPILER_ID: [UNDEFINED]\n");
#endif
#ifdef BUILD_INFO_CMAKE_C_COMPILER_VERSION
    printf("BUILD_INFO_CMAKE_C_COMPILER_VERSION: %s\n", BUILD_INFO_CMAKE_C_COMPILER_VERSION);
#else
    printf("BUILD_INFO_CMAKE_C_COMPILER_VERSION: [UNDEFINED]\n");
#endif
#ifdef BUILD_INFO_CMAKE_C_FLAGS
    printf("BUILD_INFO_CMAKE_C_FLAGS: %s\n", BUILD_INFO_CMAKE_C_FLAGS);
#else
    printf("BUILD_INFO_CMAKE_C_FLAGS: [UNDEFINED]\n");
#endif
}

void app_info(int usage){
    printf("SEIS-EIKONAL v0.1-POC (Proof of Concept) \nAuthor: Alexandr Nikitin, IPGG SB RAS\nTest: constant velocity, source in the center of domain\n");
#ifndef BIN_NAME
    char *BIN_NAME = "[UNDEFINED]";
#endif
    printf("BIN_NAME: %s\n",BIN_NAME);
#ifndef IMP_FSM3D
    char *IMP_FSM3D = "[UNDEFINED]";
#endif
    printf("IMP_FSM3D: %s\n",IMP_FSM3D);
    build_info_print();
    if(usage){
        printf("\nUsage: %s NI NJ NK BSIZE_I BSIZE_J BSIZE_K\n\t[NI x NJ x NK] - grid size\n\t[BSIZE_I x BSIZE_J x BSIZE_K] - task size\n",BIN_NAME);
    }
}

int main(int argc, char *argv[]){
    double *U, *F, H;
    double t1,t2;
    double L1, L2, Linf;
    int it, converged;
    int max_iter;
    size_t NI, NJ, NK;
    size_t SRCI, SRCJ, SRCK;
    int OMP_MAX_THREADS;

    if(argc >= 4){
        NI = atoi(argv[1]);
        NJ = atoi(argv[2]);
        NK = atoi(argv[3]);
        app_info(0);
    } else {
        app_info(1);
        return EXIT_FAILURE;
    }

    if(argc > 4 && argc < 8){
        BSIZE_I = atoi(argv[4]);
        BSIZE_J = atoi(argv[5]);
        BSIZE_K = atoi(argv[6]);
    } else {
        BSIZE_I = 1;
        BSIZE_J = 1;
        BSIZE_K = NK;
    }

    U = malloc(NI * NJ * NK * sizeof(double));
    assert(U);
    F = malloc(NI * NJ * NK * sizeof(double));
    assert(F);

    SRCI = NI/2;
    SRCJ = NJ/2;
    SRCK = NK/2;
    max_iter = 10;
    H = 2.0/NI;

    OMP_MAX_THREADS = omp_get_max_threads();
    printf("OMP_MAX_THREADS: %i\n",OMP_MAX_THREADS);
    printf("SRC: [%zu,%zu,%zu]\n",SRCI,SRCJ,SRCK);

    ex1_init(F, NI, NJ, NK);
    FSM3DInit(U, NI, NJ, NK, SRCI,SRCJ,SRCK);

    t1 = omp_get_wtime();
    it = FSM3D(U,F,H,NI,NJ,NK,SRCI,SRCJ,SRCK,max_iter,&converged);
    t2 = omp_get_wtime();
    ex1_check(U, H, NI, NJ, NK, SRCI,SRCJ,SRCK, &L1, &L2, &Linf);

    printf("====================================================\n");
    printf("NI,NJ,NK,OMP_MAX_THREADS,max_iter,it,converged,sec,L1(MAE),L2(MSE),L_inf,BSIZE_I,BSIZE_J,BSIZE_K\n");
    printf("%zu,%zu,%zu,%i,%i,%i,%i,%f,%e,%e,%e,%zu,%zu,%zu\n",
           NI,NJ,NK,OMP_MAX_THREADS,max_iter,it,converged,
           (t2 - t1),L1,L2,Linf,BSIZE_I,BSIZE_J,BSIZE_K);

    return EXIT_SUCCESS;
}
