#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"

#define TEST_ID "BRT3D_EX1"

#define DEFAULT_DOMAIN_SIZE OPENST_FLOAT_1_0
#define DEFAULT_SRC OPENST_FLOAT_0_5
#define DEFAULT_RCV OPENST_FLOAT_1_0
#define DEFAULT_GRID_SIZE 50u
#define DEFAULT_V OPENST_FLOAT_1_0


void ex_init(OPENST_FLOAT *U, OPENST_FLOAT *V, size_t NI, size_t NJ, size_t NK,
             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
             OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK){
    size_t i, j, k;
    OPENST_FLOAT di, dj, dk, dist;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                V[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = DEFAULT_V;
                di = SRCI - (OPENST_FLOAT)i * HI;
                dj = SRCJ - (OPENST_FLOAT)j * HJ;
                dk = SRCK - (OPENST_FLOAT)k * HK;
                dist = OPENST_FLOAT_SQRT(di * di + dj * dj + dk * dk);
                U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = dist / DEFAULT_V;
            }
        }
    }
}


void app_info(char *BIN_NAME,int usage){
    printf("TEST_ID: %s\n",TEST_ID);
    if(usage){
        printf("Usage: %s [NI NJ NK] [RCVI RCVJ RCVK]\n" \
               "\t[NI NJ NK] - grid size\n" \
               "\t[RCVI RCVJ RCVK] - receiver coordinates\n" \
               "Running using default values...\n\n",BIN_NAME);
    }

    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);
}


int main(int argc, char *argv[]){
    OPENST_FLOAT *U, *V, HI, HJ, HK;
    double t1,t2;
    size_t NI, NJ, NK;
    OPENST_FLOAT SRCI, SRCJ, SRCK, RCVI, RCVJ, RCVK;
    OPENST_FLOAT BRT3D_TSTEP;
    int usage_flag, errcode;
    OPENST_FLOAT *RAY;
    size_t RAY_NI, RAY_NJ;
    size_t i;
    size_t MAX_SEG;

    if(argc > 1){
        usage_flag = 0;
        /* This is only a test program, ignore safe handling of inputs */
        NI = (size_t) atoi(argv[1]);
        NJ = (size_t) atoi(argv[2]);
        NK = (size_t) atoi(argv[3]);
    } else {
        usage_flag = 1;
        NI = DEFAULT_GRID_SIZE;
        NJ = DEFAULT_GRID_SIZE;
        NK = DEFAULT_GRID_SIZE;
    }

    app_info(argv[0],usage_flag);

    U = (OPENST_FLOAT *) malloc(NI * NJ * NK * sizeof(OPENST_FLOAT));
    assert(U);
    V = (OPENST_FLOAT *) malloc(NI * NJ * NK * sizeof(OPENST_FLOAT));
    assert(V);

    if(argc > 4){
        RCVI = atof(argv[4]);
        RCVJ = atof(argv[5]);
        RCVK = atof(argv[6]);
    } else {
        RCVI = DEFAULT_RCV;
        RCVJ = DEFAULT_SRC;
        RCVK = DEFAULT_SRC;
    }

    SRCI = DEFAULT_SRC;
    SRCJ = DEFAULT_SRC;
    SRCK = DEFAULT_SRC;
    HI = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NI - 1);
    HJ = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NJ - 1);
    HK = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NK - 1);

    printf("V = %e\n",DEFAULT_V);
    printf("HI = %e; HJ = %e; HK = %e\n",HI,HJ,HK);

    ex_init(U, V, NI, NJ, NK, HI, HJ, HK, SRCI, SRCJ, SRCK);

    BRT3D_TSTEP = OpenST_BRT3D_SuggestTSTEP(DEFAULT_V, HI, HJ, HK);
    MAX_SEG = (NI * NJ * NK);

#ifdef _MSC_VER
    printf("TSTEP:\t\t%e\nMAX_SEG:\t%Iu\n",
#else
    printf("TSTEP:\t\t%e\nMAX_SEG:\t%zu\n",
#endif
		BRT3D_TSTEP,MAX_SEG);

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

#ifdef _MSC_VER
    printf("RAY:\t\t%Iu line segments\n",
#else
    printf("RAY:\t\t%zu line segments\n",
#endif
    RAY_NI);

    printf("RCV:\t\t{%e; %e; %e}\n", RCVI, RCVJ, RCVK);
    printf("RAY[start]:\t{%e; %e; %e}\n",
           RAY[OPENST_MEMADR_2D(0,0,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(0,1,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(0,2,RAY_NI,RAY_NJ)]);
    printf("RAY[end]:\t\t{%e; %e; %e}\n",
           RAY[OPENST_MEMADR_2D(RAY_NI - 1,0,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(RAY_NI - 1,1,RAY_NI,RAY_NJ)],
            RAY[OPENST_MEMADR_2D(RAY_NI - 1,2,RAY_NI,RAY_NJ)]);
    printf("SRC:\t\t{%e; %e; %e}\n", SRCI, SRCJ, SRCK);

    for(i = 0; i < RAY_NI; ++i){
#ifdef _MSC_VER
        printf("%Iu %e %e %e\n",
#else
        printf("%zu %e %e %e\n",
#endif
               i,
               RAY[OPENST_MEMADR_2D(i,0,RAY_NI,RAY_NJ)],
               RAY[OPENST_MEMADR_2D(i,1,RAY_NI,RAY_NJ)],
               RAY[OPENST_MEMADR_2D(i,2,RAY_NI,RAY_NJ)]);
    }

EXIT:
    return errcode;
}
