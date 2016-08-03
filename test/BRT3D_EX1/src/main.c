#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"

#define TEST_ID "BRT3D_EX1"

#define DEFAULT_DOMAIN_SIZE 1.0
#define DEFAULT_SRC 0.5
#define DEFAULT_RCV 0.1
#define DEFAULT_GRID_SIZE 100u


void ex_init(double *U, double *V, size_t NI, size_t NJ, size_t NK,
             double HI, double HJ, double HK,
             double SRCI, double SRCJ, double SRCK){
    size_t i, j, k;
    double di, dj, dk, dist;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                V[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = 1.0;
                di = SRCI - (double)i * HI;
                dj = SRCJ - (double)j * HJ;
                dk = SRCK - (double)k * HK;
                dist = sqrt(di * di + dj * dj + dk * dk);
                U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = dist;
            }
        }
    }
}


void app_info(char *BIN_NAME,int usage){
    printf("TEST_ID: %s\n",TEST_ID);
    if(usage){
        printf("Usage: %s [NI NJ NK]\n" \
               "\t[NI NJ NK] - grid size\n"
               "Running using default values...\n\n",BIN_NAME);
    }

    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);
}


int main(int argc, char *argv[]){
    double *U, *V, HI, HJ, HK;
    double t1,t2;
    size_t NI, NJ, NK;
    double SRCI, SRCJ, SRCK, RCVI, RCVJ, RCVK;
    double BRT3D_TSTEP;
    int usage_flag, errcode;
    double *RAY;
    size_t RAY_NI, RAY_NJ;
    size_t i;
    size_t MAX_SEG;

    if(argc > 1){
        usage_flag = 0;
        NI = atoi(argv[1]);
        NJ = atoi(argv[2]);
        NK = atoi(argv[3]);
    } else {
        usage_flag = 1;
        NI = DEFAULT_GRID_SIZE;
        NJ = DEFAULT_GRID_SIZE;
        NK = DEFAULT_GRID_SIZE;
    }

    app_info(argv[0],usage_flag);

    U = malloc(NI * NJ * NK * sizeof(double));
    assert(U);
    V = malloc(NI * NJ * NK * sizeof(double));
    assert(V);

    RCVI = DEFAULT_RCV;
    RCVJ = DEFAULT_RCV;
    RCVK = DEFAULT_RCV;
    SRCI = DEFAULT_SRC;
    SRCJ = DEFAULT_SRC;
    SRCK = DEFAULT_SRC;
    HI = DEFAULT_DOMAIN_SIZE / (double)(NI - 1);
    HJ = DEFAULT_DOMAIN_SIZE / (double)(NJ - 1);
    HK = DEFAULT_DOMAIN_SIZE / (double)(NK - 1);

    printf("HI = %e; HJ = %e; HK = %e\n",HI,HJ,HK);

    ex_init(U, V, NI, NJ, NK, HI, HJ, HK, SRCI, SRCJ, SRCK);

    BRT3D_TSTEP = OpenST_BRT3D_SuggestTSTEP(1.0, HI, HJ, HK);
    MAX_SEG = (NI * NJ * NK);
    printf("TSTEP:\t\t%e\nMAX_SEG:\t%zu\n",BRT3D_TSTEP,MAX_SEG);

    t1 = omp_get_wtime();
    errcode = OpenST_BRT3D_Trace(U, V, NI, NJ, NK, HI, HJ, HK, BRT3D_TSTEP,
                                 RCVI, RCVJ, RCVK, SRCI, SRCJ, SRCK,
                                 MAX_SEG,
                                 &RAY, &RAY_NI, &RAY_NJ);
    t2 = omp_get_wtime();

    if(errcode != OPENST_ERR_SUCCESS){
        fprintf(stderr,"errcode: %i\n",errcode);
        goto EXIT;
    }

    printf("BRT3D time:\t%e sec\n", t2 - t1);
    printf("RAY:\t\t%zu line segments\n",RAY_NI);
    printf("RCV:\t\t{%e; %e; %e}\n", RCVI, RCVJ, RCVK);
    printf("RAY[start]:\t{%e; %e; %e}\n",
           RAY[OPENST_MEMADR_2D(0,0,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(0,1,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(0,2,RAY_NI,RAY_NJ)]);
    printf("RAY[end]:\t{%e; %e; %e}\n",
           RAY[OPENST_MEMADR_2D(RAY_NI - 1,0,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(RAY_NI - 1,1,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(RAY_NI - 1,2,RAY_NI,RAY_NJ)]);
    printf("SRC:\t\t{%e; %e; %e}\n", SRCI, SRCJ, SRCK);

    for(i = 0; i < RAY_NI; ++i){
        printf("%zu %e %e %e\n",i,
               RAY[OPENST_MEMADR_2D(i,0,RAY_NI,RAY_NJ)],
                RAY[OPENST_MEMADR_2D(i,1,RAY_NI,RAY_NJ)],
                RAY[OPENST_MEMADR_2D(i,2,RAY_NI,RAY_NJ)]);
    }

EXIT:
    return errcode;
}
