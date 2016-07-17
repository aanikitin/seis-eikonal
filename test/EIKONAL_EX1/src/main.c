#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "openst.h"

#define TEST_ID "EIKONAL_EX1"

#define DEFAULT_SIZE 100u
#define DEFAULT_MAX_ITER 10
#define DEFAULT_EPS_MULT 0.0


void ex_init(double *V, size_t NI, size_t NJ, size_t NK){
    size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                V[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = 1.0;
            }
        }
    }
}


void ex_check(double *U, double H, size_t NI, size_t NJ, size_t NK,
              size_t SRCI, size_t SRCJ, size_t SRCK,
              double *L1, double *L2, double *Linf,
              double *Umin, double *Umax, double *Umean
              ){
    size_t i, j, k, NN;
    double uval, umin, umax, umean;
    double di, dj, dk, dist, diff, l1, l2, linf;
    l1 = 0.0;
    l2 = 0.0;
    linf = -INFINITY;
    umin = INFINITY;
    umax = -INFINITY;
    umean = 0.0;
    NN = NI * NJ * NK;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                uval = U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
                umin = fmin(umin,uval);
                umax = fmax(umax,uval);
                umean += uval;
                di = ((double)SRCI - (double)i) * H;
                dj = ((double)SRCJ - (double)j) * H;
                dk = ((double)SRCK - (double)k) * H;
                dist = sqrt(di * di + dj * dj + dk * dk);
                diff = uval - dist;
                l1 += fabs(diff);
                l2 += (diff * diff);
                linf = fmax(linf,fabs(diff));
            }
        }
    }
    *L1 = l1/NN;
    *L2 = l2/NN;
    *Linf = linf;
    *Umin = umin;
    *Umax = umax;
    *Umean = umean/NN;
}


void app_info(char *BIN_NAME,int usage){
    printf("TEST_ID: %s\n",TEST_ID);
    if(usage){
        printf("Usage: %s [N] [BSIZE_I BSIZE_J BSIZE_K] [EPS_MULT] [MAX_ITER]\n" \
               "\t[N] - grid size [N x N x N]\n" \
               "\t[BSIZE_I x BSIZE_J x BSIZE_K] - task size\n" \
               "\t[EPS_MULT] - EPS convergence parameter multiplier\n" \
               "\t[MAX_ITER] - maximum number of iterations\n"
               "Running using default values...\n\n",BIN_NAME);
    }

    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);
}


int main(int argc, char *argv[]){
    double *U, *V, H;
    int max_iter;
    int it, converged;
    double t1,t2;
    double L1, L2, Linf, Umin, Umax, Umean;
    size_t NI, NJ, NK;
    size_t SRCI, SRCJ, SRCK;
    size_t BSIZE_I;
    size_t BSIZE_J;
    size_t BSIZE_K;
    char *LSM_UNLOCKED;
    int OMP_MAX_THREADS;
    double vmin, vmax, vmean;
    double EPS, EPS_MULT, EIK3D_Time, BRT3D_TSTEP;
    const char *IMP_NAME;
    int usage_flag, errcode;

    NI = DEFAULT_SIZE;
    EPS_MULT = DEFAULT_EPS_MULT;
    max_iter = DEFAULT_MAX_ITER;
    usage_flag = 1;

    if(argc > 1){
        NI = atoi(argv[1]);
        usage_flag = 0;
    }
    
    NJ = NI;
    NK = NI;

    if(argc >= 5){
        BSIZE_I = atoi(argv[2]);
        BSIZE_J = atoi(argv[3]);
        BSIZE_K = atoi(argv[4]);
    } else {
        OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);
    }

    if(argc > 5){
        EPS_MULT = atof(argv[5]);
    }

    if(argc > 6){
        max_iter = atoi(argv[6]);
    }

    app_info(argv[0],usage_flag);

    U = malloc(NI * NJ * NK * sizeof(double));
    assert(U);
    V = malloc(NI * NJ * NK * sizeof(double));
    assert(V);
    LSM_UNLOCKED = malloc(NI * NJ * NK * sizeof(char));
    assert(LSM_UNLOCKED);

    SRCI = NI/2;
    SRCJ = NJ/2;
    SRCK = NK/2;
    H = 1.0/(NI - 1);

    OMP_MAX_THREADS = omp_get_max_threads();

    ex_init(V, NI, NJ, NK);

    if(EPS_MULT != 0.0){
        OpenST_AOP_GetArrStats(V, NI * NJ * NK, &vmin, &vmax, &vmean);
        printf("V[min;max;mean]: [%e;%e;%e]\n",vmin,vmax,vmean);
        BRT3D_TSTEP = OpenST_BRT3D_SuggestTSTEP(vmax, H, H, H);
        EPS = BRT3D_TSTEP * EPS_MULT;
    } else {
        EPS = 0.0;
    }

#ifndef TEST_FSM
    IMP_NAME = OPENST_LSM3D_IMP_NAME;
    OpenST_LSM3D_Init(U,LSM_UNLOCKED, NI, NJ, NK, SRCI,SRCJ,SRCK);
    t1 = omp_get_wtime();
    it = OpenST_LSM3D_Compute(U,LSM_UNLOCKED,V,H,NI,NJ,NK,SRCI,SRCJ,SRCK,
                              max_iter,&converged,
                              BSIZE_I,BSIZE_J,BSIZE_K,EPS_MULT);
    t2 = omp_get_wtime();
    EIK3D_Time = t2 - t1;
#else
    IMP_NAME = OPENST_FSM3D_IMP_NAME;
    OpenST_FSM3D_Init(U,NI,NJ,NK,SRCI,SRCJ,SRCK);
    t1 = omp_get_wtime();
    it = OpenST_FSM3D_Compute(U,V,H,NI,NJ,NK,SRCI,SRCJ,SRCK,
                              max_iter,&converged,
                              BSIZE_I,BSIZE_J,BSIZE_K,EPS);
    t2 = omp_get_wtime();
    EIK3D_Time = t2 - t1;
#endif

    ex_check(U, H, NI, NJ, NK, SRCI,SRCJ,SRCK, &L1, &L2, &Linf, &Umin, &Umax, &Umean);

    printf("====================================================\n");
    printf("TEST_ID,METHOD_ID,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K," \
           "EPS,max_iter,it,converged,sec,L1(MAE),L2(MSE),L_inf," \
           "Umin,Umean,Umax," \
           "NI,NJ,NK,SRCI,SRCJ,SRCK,OPENST_BUILDINFO_LINK_TYPE\n");

#ifdef _WIN32
    printf("%s,%s,%i,%Iu,%Iu,%Iu,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e," \
	   "%Iu,%Iu,%Iu,%Iu,%Iu,%Iu,%s\n",
#else
    printf("%s,%s,%i,%zu,%zu,%zu,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e," \
           "%zu,%zu,%zu,%zu,%zu,%zu,%s\n",
#endif
           TEST_ID,IMP_NAME,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,
           EPS,max_iter,it,converged,EIK3D_Time,L1,L2,Linf,Umin,Umean,Umax,
           NI,NJ,NK,SRCI,SRCJ,SRCK,OPENST_BUILDINFO_LINK_TYPE_STATIC
           );

    if(converged && it == 9){
        errcode = EXIT_SUCCESS;
    } else {
        errcode = EXIT_FAILURE;
    }
    return errcode;
}
