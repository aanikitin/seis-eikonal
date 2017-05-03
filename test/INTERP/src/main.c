#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"

#define TEST_ID "INTERP"

#define DEFAULT_DOMAIN_SIZE 1.0
#define DEFAULT_GRID_SIZE 2u


void ex_init(double *V, double HI, double HJ, double HK, size_t NI, size_t NJ, size_t NK){
    size_t i, j, k;
    double PI,PJ,PK;
    printf("Input array:\nCOORD,VALUE\n");
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                PI = (double)i * HI;
                PJ = (double)j * HJ;
                PK = (double)k * HK;
                V[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = (double)OPENST_MEMADR_3D(i,j,k,NI,NJ,NK);
                printf("[%e,%e,%e],%e\n",PI,PJ,PK,V[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)]);
            }
        }
    }
}


void app_info(char *BIN_NAME,int usage){
    printf("TEST_ID: %s\n",TEST_ID);
    if(usage){
        printf("Usage: %s [NI NJ NK] [BSIZE_I BSIZE_J BSIZE_K] [EPS_MULT] [MAX_ITER]\n" \
               "\t[NI NJ NK] - grid size\n" \
               "\t[BSIZE_I x BSIZE_J x BSIZE_K] - task size\n" \
               "\t[EPS_MULT] - EPS convergence parameter multiplier\n" \
               "\t[MAX_ITER] - maximum number of iterations\n"
               "Running using default values...\n\n",BIN_NAME);
    }

    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);
}


int main(int argc, char *argv[]){

    OPENST_ERR errcode;
    double *V;
    size_t NI, NJ, NK, i, j, k;
    double HI, HJ, HK, PI, PJ, PK, value, t1, t2;

    NI = DEFAULT_GRID_SIZE;
    NJ = DEFAULT_GRID_SIZE;
    NK = DEFAULT_GRID_SIZE;

    app_info(argv[0],0);

    V = (double *)malloc(NI * NJ * NK * sizeof(double));
    assert(V);

    HI = DEFAULT_DOMAIN_SIZE / (double)(NI - 1);
    HJ = DEFAULT_DOMAIN_SIZE / (double)(NJ - 1);
    HK = DEFAULT_DOMAIN_SIZE / (double)(NK - 1);

    ex_init(V,HI,HJ,HK,NI,NJ,NK);

    printf("\nPerforming default (linear) interpolation...\n");
    printf("====================================================\n");
    printf("COORD,INTERP_VALUE,TIME\n");
    for(i = 0; i < (2 * NI - 1); ++i){
        for(j = 0; j < (2 * NJ - 1); ++j){
            for(k = 0; k < (2 * NK - 1); ++k){
                PI = (double)i * HI/2.0;
                PJ = (double)j * HJ/2.0;
                PK = (double)k * HK/2.0;
                t1 = omp_get_wtime();
                errcode = OpenST_INTERP_3D(V,NI,NJ,NK,HI,HJ,HK,PI,PJ,PK,&value);
                t2 = omp_get_wtime();
                if(!errcode){
                    printf("[%e,%e,%e],%e,%e\n",PI,PJ,PK,value,t2 - t1);
                } else {
                    printf("[%e,%e,%e],[error %i],%e\n",PI,PJ,PK,errcode,t2 - t1);
                }
            }
        }
    }

    return EXIT_SUCCESS;
}
