#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"
#include "fio_hdf5.h"

#define TEST_ID "EIKONAL_HDF5"

#define DEFAULT_SRC 0.5
#define DEFAULT_MAX_ITER 1000
#define DEFAULT_EPS_MULT 0.01
#define DEFAULT_V_DATASET_NAME "V"


void app_info(char *BIN_NAME,int usage){
    printf("TEST_ID: %s\n",TEST_ID);
    if(usage){
        printf("Usage: %s FILE_IN FILE_OUT SRCI SRCJ SRCK [DATASET] [EPS_MULT] [MAX_ITER] [BSIZE_I BSIZE_J BSIZE_K]\n" \
               "\tFILE_IN - HDF5 input file with datasets:\n" \
               "\t\t\"V\" or \"[DATASET]\" - array representing 3D regular velocity grid stored in row major format\n" \
               "\t\t\"HI\", \"HJ\", \"HK\" - grid steps in each dimension\n" \
               "\tFILE_OUT - HDF5 output file for eikonal's solution\n" \
               "\tSRCI SRCJ SRCK - source coordinates, with (0,0,0) corresponding to the first (i,j,k) element of velocity array\n" \
               "\t[DATASET] - velocity grid dataset name, defaults to \"V\"\n" \
               "\t[EPS_MULT] - EPS convergence parameter multiplier\n" \
               "\t[MAX_ITER] - maximum number of iterations\n" \
               "\t[BSIZE_I x BSIZE_J x BSIZE_K] - task size\n\n",BIN_NAME);
    }

    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);
}


int main(int argc, char *argv[]){

    OPENST_ERR errcode;
    char *FILE_IN, *FILE_OUT, *V_DATASET_NAME;
    double *U, *V, vmin, vmean,vmax, HI, HJ, HK;
    int max_iter;
    int it, converged;
    double t1,t2;
    size_t NI, NJ, NK;
    double SRCI, SRCJ, SRCK;
    size_t BSIZE_I;
    size_t BSIZE_J;
    size_t BSIZE_K;
    char *LSM_UNLOCKED;
    int OMP_MAX_THREADS;
    double EPS, EPS_MULT, EIK3D_Time;
    const char *IMP_NAME;
    size_t SRCidx_i, i, j, k;
    size_t *SRCidx, SRCidx_NI, SRCidx_NJ;

    if(argc < 6){
        fprintf(stderr,"Error: invalid command line parameters\n");
        app_info(argv[0],1);
        errcode = OPENST_ERR_PARAM_INVALID;
        goto EXIT;
    }

    app_info(argv[0],0);

    FILE_IN = argv[1];
    FILE_OUT = argv[2];

    SRCI = atof(argv[3]);
    SRCJ = atof(argv[4]);
    SRCK = atof(argv[5]);

    printf("SRCI = %e; SRCJ = %e; SRCK = %e\n",SRCI,SRCJ,SRCK);

    if(argc > 6){
        V_DATASET_NAME = argv[6];
    } else {
        V_DATASET_NAME = DEFAULT_V_DATASET_NAME;
    }

    if(argc > 7){
        EPS_MULT = atof(argv[7]);
    } else {
        EPS_MULT = DEFAULT_EPS_MULT;
    }

    if(argc > 8){
        max_iter = atoi(argv[8]);
    } else {
        max_iter = DEFAULT_MAX_ITER;
    }

    hdf5_read_model(FILE_IN,V_DATASET_NAME,&V,&HI,&HJ,&HK,&NI,&NJ,&NK);

#ifdef _MSC_VER
    printf("NI = %Iu; NJ = %Iu; NK = %Iu\n",NI,NJ,HK);
#else
    printf("NI = %zu; NJ = %zu; NK = %zu\n",NI,NJ,NK);
#endif
    printf("HI = %e; HJ = %e; HK = %e\n",HI,HJ,HK);
    printf("V[0,0,0]: %e\n",V[OPENST_MEMADR_3D(0,0,0,NI,NJ,NK)]);
    printf("V[0,0,1]: %e\n",V[OPENST_MEMADR_3D(0,0,1,NI,NJ,NK)]);
    printf("V[0,1,0]: %e\n",V[OPENST_MEMADR_3D(0,1,0,NI,NJ,NK)]);
    printf("V[1,0,0]: %e\n",V[OPENST_MEMADR_3D(1,0,0,NI,NJ,NK)]);

    if(OpenST_CRS_IsPointNotWithinBounds(SRCI,SRCJ,SRCK,NI,NJ,NK,HI,HJ,HK)){
        fprintf(stderr,"Error: SRC is not within domain bounds\n");
        errcode = OPENST_ERR_PARAM_INVALID;
        goto EXIT;
    }

    if(argc > 9){
        BSIZE_I = (size_t) atoi(argv[9]);
        BSIZE_J = (size_t) atoi(argv[10]);
        BSIZE_K = (size_t) atoi(argv[11]);
    } else {
        OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);
    }

    U = (double *)malloc(NI * NJ * NK * sizeof(double));
    assert(U);
    LSM_UNLOCKED = (char *)malloc(NI * NJ * NK * sizeof(char));
    assert(LSM_UNLOCKED);

    OMP_MAX_THREADS = omp_get_max_threads();

    if(EPS_MULT != 0.0){
        OpenST_AOP_GetArrStats(V, NI * NJ * NK, &vmin, &vmax, &vmean);
        printf("Vmin = %e; Vmean = %e; Vmax = %e\n",vmin,vmean,vmax);
        EPS = EPS_MULT * OpenST_BRT3D_SuggestTSTEP(vmax, HI, HJ, HK);
    } else {
        EPS = 0.0;
    }

    t1 = omp_get_wtime();
#ifndef TEST_FSM
    IMP_NAME = OPENST_LSM3D_IMP_NAME;
    errcode = OpenST_LSM3D_Init_2(U,LSM_UNLOCKED,V,
                        NI,NJ,NK,
                        HI,HJ,HK,
                        SRCI,SRCJ,SRCK,
                        &SRCidx,&SRCidx_NI,&SRCidx_NJ,
                        OPENST_FSM3D_INIT_DEFAULT);
#else
    IMP_NAME = OPENST_FSM3D_IMP_NAME;
    errcode = OpenST_FSM3D_Init_2(U,V,
                        NI,NJ,NK,
                        HI,HJ,HK,
                        SRCI,SRCJ,SRCK,
                        &SRCidx,&SRCidx_NI,&SRCidx_NJ,
                        OPENST_FSM3D_INIT_DEFAULT);
#endif
    t2 = omp_get_wtime();

    if(errcode != OPENST_ERR_SUCCESS){
        fprintf(stderr,"Error: Initialization errcode = %i\n",errcode);
        goto EXIT;
    }

    printf("Initialization time: %e sec\n", t2 - t1);
    printf("U initialized in:\n");
    for(SRCidx_i = 0; SRCidx_i < SRCidx_NI; ++SRCidx_i){
        i = SRCidx[OPENST_MEMADR_2D(SRCidx_i,0,SRCidx_NI,SRCidx_NJ)];
        j = SRCidx[OPENST_MEMADR_2D(SRCidx_i,1,SRCidx_NI,SRCidx_NJ)];
        k = SRCidx[OPENST_MEMADR_2D(SRCidx_i,2,SRCidx_NI,SRCidx_NJ)];
#ifdef _MSC_VER
        printf("U[%Iu,%Iu,%Iu] = %e\n",i,j,k,U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)]);
#else
        printf("U[%zu,%zu,%zu] = %e\n",i,j,k,U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)]);
#endif
    }

    t1 = omp_get_wtime();
#ifndef TEST_FSM
    it = OpenST_LSM3D_Compute(U,LSM_UNLOCKED,V,
                              NI,NJ,NK,
                              HI,HJ,HK,
                              max_iter,&converged,
                              BSIZE_I,BSIZE_J,BSIZE_K,EPS);
#else
    it = OpenST_FSM3D_Compute(U,V,
                              NI,NJ,NK,
                              HI,HJ,HK,
                              max_iter,&converged,
                              BSIZE_I,BSIZE_J,BSIZE_K,EPS);
#endif
    t2 = omp_get_wtime();
    EIK3D_Time = t2 - t1;

    hdf5_write_time(FILE_OUT,"T",U,&HI,&HJ,&HK,NI,NJ,NK);

    printf("\n====================================================\n");
    printf("TEST_ID,METHOD_ID,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K," \
           "EPS,max_iter,it,converged,sec," \
           "NI,NJ,NK,SRCI,SRCJ,SRCK,OPENST_BUILDINFO_LINK_TYPE\n");

#ifdef _MSC_VER
    printf("%s,%s,%i,%Iu,%Iu,%Iu,%e,%i,%i,%i,%e," \
           "%Iu,%Iu,%Iu,%e,%e,%e,%s\n",
#else
    printf("%s,%s,%i,%zu,%zu,%zu,%e,%i,%i,%i,%e," \
           "%zu,%zu,%zu,%e,%e,%e,%s\n",
#endif
           TEST_ID,IMP_NAME,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,
           EPS,max_iter,it,converged,EIK3D_Time,
           NI,NJ,NK,SRCI,SRCJ,SRCK,OPENST_BUILDINFO_LINK_TYPE_STATIC
           );

EXIT:
    return errcode;
}
