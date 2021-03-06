#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"

#define TEST_ID "EIKONAL_EX1"

#define DEFAULT_DOMAIN_SIZE OPENST_FLOAT_1_0
#define DEFAULT_SRC OPENST_FLOAT_0_5
#define DEFAULT_GRID_SIZE 50u
#define DEFAULT_MAX_ITER 10
#define DEFAULT_EPS_MULT OPENST_FLOAT_1_0
#define DEFAULT_V OPENST_FLOAT_1_0


void ex_init(OPENST_FLOAT *V, size_t NI, size_t NJ, size_t NK){
    size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                V[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = DEFAULT_V;
            }
        }
    }
}


void ex_check(OPENST_FLOAT *U,
              size_t NI, size_t NJ, size_t NK,
              OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
              OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
              OPENST_FLOAT *L1, OPENST_FLOAT *L2, OPENST_FLOAT *Linf,
              OPENST_FLOAT *Umin, OPENST_FLOAT *Umax, OPENST_FLOAT *Umean){

    size_t i, j, k, NN;
    OPENST_FLOAT uval, umin, umax, umean;
    OPENST_FLOAT di, dj, dk, dist, diff, l1, l2, linf;
    l1 = OPENST_FLOAT_0_0;
    l2 = OPENST_FLOAT_0_0;
    linf = -OPENST_FLOAT_INF;
    umin = OPENST_FLOAT_INF;
    umax = -OPENST_FLOAT_INF;
    umean = OPENST_FLOAT_0_0;
    NN = NI * NJ * NK;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                uval = U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
                umin = OPENST_FLOAT_FMIN(umin,uval);
                umax = OPENST_FLOAT_FMAX(umax,uval);
                umean += uval;
                di = SRCI - (OPENST_FLOAT)i * HI;
                dj = SRCJ - (OPENST_FLOAT)j * HJ;
                dk = SRCK - (OPENST_FLOAT)k * HK;
                dist = OPENST_FLOAT_SQRT(di * di + dj * dj + dk * dk);
                diff = uval - dist / DEFAULT_V;
                l1 += OPENST_FLOAT_FABS(diff);
                l2 += (diff * diff);
                linf = OPENST_FLOAT_FMAX(linf,OPENST_FLOAT_FABS(diff));
            }
        }
    }
    *L1 = l1/(OPENST_FLOAT)NN;
    *L2 = l2/(OPENST_FLOAT)NN;
    *Linf = linf;
    *Umin = umin;
    *Umax = umax;
    *Umean = umean/(OPENST_FLOAT)NN;
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
    OPENST_FLOAT *U, *V, HI, HJ, HK;
    int max_iter;
    int it, converged;
    double t1,t2;
    OPENST_FLOAT L1, L2, Linf, Umin, Umax, Umean;
    size_t NI, NJ, NK;
    OPENST_FLOAT SRCI, SRCJ, SRCK;
    size_t BSIZE_I;
    size_t BSIZE_J;
    size_t BSIZE_K;
    char *LSM_UNLOCKED;
    int OMP_MAX_THREADS;
    OPENST_FLOAT vmin, vmax, vmean;
    OPENST_FLOAT EPS, EPS_MULT, EIK3D_Time, BRT3D_TSTEP;
    const char *IMP_NAME, *IMP_BLOCKSERIAL_NAME;
    int usage_flag;
    size_t SRCidx_i, i, j, k;
    size_t *SRCidx, SRCidx_NI, SRCidx_NJ;

    if(argc > 1){
        usage_flag = 0;
        NI = (size_t) atoi(argv[1]);
        NJ = (size_t) atoi(argv[2]);
        NK = (size_t) atoi(argv[3]);
    } else {
        usage_flag = 1;
        NI = DEFAULT_GRID_SIZE;
        NJ = DEFAULT_GRID_SIZE;
        NK = DEFAULT_GRID_SIZE;
    }

    if(argc > 4){
        /* This is only a test program, ignore safe handling of inputs */
        BSIZE_I = (size_t) atoi(argv[4]);
        BSIZE_J = (size_t) atoi(argv[5]);
        BSIZE_K = (size_t) atoi(argv[6]);
    } else {
        OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);
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

    app_info(argv[0],usage_flag);

    U = (OPENST_FLOAT *)malloc(NI * NJ * NK * sizeof(OPENST_FLOAT));
    assert(U);
    V = (OPENST_FLOAT *)malloc(NI * NJ * NK * sizeof(OPENST_FLOAT));
    assert(V);
    LSM_UNLOCKED = (char *)malloc(NI * NJ * NK * sizeof(char));
    assert(LSM_UNLOCKED);

    SRCI = DEFAULT_SRC;
    SRCJ = DEFAULT_SRC;
    SRCK = DEFAULT_SRC;
    HI = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NI - 1);
    HJ = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NJ - 1);
    HK = DEFAULT_DOMAIN_SIZE / (OPENST_FLOAT)(NK - 1);

    printf("V = %e\n",DEFAULT_V);
    printf("HI = %e; HJ = %e; HK = %e\n",HI,HJ,HK);

    OMP_MAX_THREADS = omp_get_max_threads();

    ex_init(V, NI, NJ, NK);

    if(EPS_MULT != OPENST_FLOAT_0_0){
        OpenST_AOP_GetArrStats(V, NI * NJ * NK, &vmin, &vmax, &vmean);
        printf("Vmin = %e; Vmean = %e; Vmax = %e\n",vmin,vmean,vmax);
        BRT3D_TSTEP = OpenST_BRT3D_SuggestTSTEP(vmax, HI, HJ, HK);
        EPS = BRT3D_TSTEP * EPS_MULT;
    } else {
        EPS = OPENST_FLOAT_0_0;
    }

    t1 = omp_get_wtime();
#ifndef TEST_FSM
    IMP_NAME = OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME;
    IMP_BLOCKSERIAL_NAME = OPENST_LSM3D_BLOCKSERIAL_IMP_NAME;
    errcode = OpenST_LSM3D_Init_2(U,LSM_UNLOCKED,V,
                        NI,NJ,NK,
                        HI,HJ,HK,
                        SRCI,SRCJ,SRCK,
                        &SRCidx,&SRCidx_NI,&SRCidx_NJ,
                        OPENST_FSM3D_INIT_DEFAULT);
#else
    IMP_NAME = OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME;
    IMP_BLOCKSERIAL_NAME = OPENST_FSM3D_BLOCKSERIAL_IMP_NAME;
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

    ex_check(U, NI, NJ, NK, HI, HJ, HK,
             SRCI, SRCJ, SRCK, &L1, &L2, &Linf, &Umin, &Umax, &Umean);

    printf("\n====================================================\n");
    printf("TEST_ID,IMP_ID,BLOCK_IMP_ID,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K," \
           "EPS,max_iter,it,converged,sec,L1(MAE),L2(MSE),L_inf," \
           "Umin,Umean,Umax," \
           "NI,NJ,NK,SRCI,SRCJ,SRCK,OPENST_BUILDINFO_LINK_TYPE\n");

#ifdef _MSC_VER
    printf("%s,%s,%s,%i,%Iu,%Iu,%Iu,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e," \
           "%Iu,%Iu,%Iu,%e,%e,%e,%s\n",
#else
    printf("%s,%s,%s,%i,%zu,%zu,%zu,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e," \
           "%zu,%zu,%zu,%e,%e,%e,%s\n",
#endif
           TEST_ID,IMP_NAME,IMP_BLOCKSERIAL_NAME,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,
           EPS,max_iter,it,converged,EIK3D_Time,L1,L2,Linf,Umin,Umean,Umax,
           NI,NJ,NK,SRCI,SRCJ,SRCK,OPENST_BUILDINFO_LINK_TYPE_STATIC
           );

EXIT:
    if(converged && (it == 9) && (errcode == OPENST_ERR_SUCCESS)){
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}
