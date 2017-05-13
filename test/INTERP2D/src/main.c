#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"

#define TEST_ID "INTERP2D"

#define DEFAULT_DOMAIN_SIZE OPENST_FLOAT_1_0
#define DEFAULT_GRID_SIZE 2u


void ex_init(OPENST_FLOAT *V, OPENST_FLOAT HI, OPENST_FLOAT HJ,
             size_t NI, size_t NJ){
    size_t i, j;
    OPENST_FLOAT PI,PJ;
    printf("Input array:\nCOORD,VALUE\n");
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            PI = (OPENST_FLOAT)i * HI;
            PJ = (OPENST_FLOAT)j * HJ;
            V[OPENST_MEMADR_2D(i,j,NI,NJ)] = (OPENST_FLOAT)OPENST_MEMADR_2D(i,j,NI,NJ);
            printf("[%e,%e],%e\n",PI,PJ,V[OPENST_MEMADR_2D(i,j,NI,NJ)]);
        }
    }
}


void app_info(char *BIN_NAME,int usage){
    printf("TEST_ID: %s\n",TEST_ID);

    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);
}


int main(int argc, char *argv[]){

    OPENST_ERR errcode;
    OPENST_FLOAT *V;
    size_t NI, NJ, i, j;
    OPENST_FLOAT HI, HJ, PI, PJ, value;
    double t1, t2;

    NI = DEFAULT_GRID_SIZE;
    NJ = DEFAULT_GRID_SIZE;

    app_info(argv[0],0);

    V = (OPENST_FLOAT *)malloc(NI * NJ * sizeof(OPENST_FLOAT));
    assert(V);

    HI = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NI - 1);
    HJ = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NJ - 1);
    
    ex_init(V,HI,HJ,NI,NJ);

    printf("\nPerforming default (linear) interpolation...\n");
    printf("====================================================\n");
    printf("COORD,INTERP_VALUE,TIME\n");
    for(i = 0; i < (2 * NI - 1); ++i){
        for(j = 0; j < (2 * NJ - 1); ++j){
            PI = (OPENST_FLOAT)i * HI/OPENST_FLOAT_2_0;
            PJ = (OPENST_FLOAT)j * HJ/OPENST_FLOAT_2_0;
            t1 = omp_get_wtime();
            errcode = OpenST_INTERP_2D(V,NI,NJ,HI,HJ,PI,PJ,&value);
            t2 = omp_get_wtime();
            if(!errcode){
                printf("[%e,%e],%e,%e\n",PI,PJ,value,t2 - t1);
            } else {
                printf("[%e,%e],[error %i],%e\n",PI,PJ,errcode,t2 - t1);
            }
        }
    }

    return EXIT_SUCCESS;
}
