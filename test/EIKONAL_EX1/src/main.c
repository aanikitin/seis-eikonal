#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

//#define OPENST_LINK_SHARED
#include "openst.h"
#include "openst.h"

#define TEST_NAME "EIKONAL_EX1"

#define BSIZE_DEF 8u

void ex1_init(double *F, size_t NI, size_t NJ, size_t NK){
    size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                F[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = 1.0;
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
    linf = -INFINITY;
    NN = NI * NJ * NK;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                di = ((double)SRCI - (double)i) * H;
                dj = ((double)SRCJ - (double)j) * H;
                dk = ((double)SRCK - (double)k) * H;
                dist = sqrt(di * di + dj * dj + dk * dk);
                diff = U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] - dist;
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

void app_info(char *BIN_NAME,int usage){
    if(usage){
        printf("\nUsage: %s NI NJ NK [BSIZE_I BSIZE_J BSIZE_K]\n\t[NI x NJ x NK] - grid size\n\t[BSIZE_I x BSIZE_J x BSIZE_K] - task size\n",BIN_NAME);
    }
}

int main(int argc, char *argv[]){
    double *U, *V, H;
    double t1,t2;
    double L1, L2, Linf;
    int it, converged;
    int max_iter;
    size_t NI, NJ, NK;
    size_t SRCI, SRCJ, SRCK;
    size_t BSIZE_I;
    size_t BSIZE_J;
    size_t BSIZE_K;
    char *LSM_UNLOCKED;
    double *RAY;
    size_t RAY_NI, RAY_NJ;
    size_t RCVI, RCVJ, RCVK;
    int OMP_MAX_THREADS;
    double vmin, vmax, vmean;
    double EPS, LSM3D_Time, BRT3D_Time, BRT3D_TSTEP;
    int BRT3D_ERR;
    const char *IMP_NAME;
    
    printf("%s\n",OPENST_VERSION_STR_FULL_STATIC);
    printf("%s\n",OPENST_BUILDINFO_STR_FULL_STATIC);

    if(argc >= 4){
        NI = atoi(argv[1]);
        NJ = atoi(argv[2]);
        NK = atoi(argv[3]);
        app_info(argv[0],0);
    } else {
        app_info(argv[0],1);
        return EXIT_FAILURE;
    }

    if(argc >= 7){
        BSIZE_I = atoi(argv[4]);
        BSIZE_J = atoi(argv[5]);
        BSIZE_K = atoi(argv[6]);
    } else {
        OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);
    }

    U = malloc(NI * NJ * NK * sizeof(double));
    assert(U);
    V = malloc(NI * NJ * NK * sizeof(double));
    assert(V);
    LSM_UNLOCKED = malloc(NI * NJ * NK * sizeof(char));
    assert(LSM_UNLOCKED);

    SRCI = NI/2;
    SRCJ = NJ/2;
    SRCK = NK/2;
    max_iter = 10;
    H = 2.0/NI;
    
    RCVI = 0;
    RCVJ = 0;
    RCVK = 0;

    OMP_MAX_THREADS = omp_get_max_threads();
    printf("OMP_MAX_THREADS: %i\n",OMP_MAX_THREADS);

#if defined(_WIN32) || defined(_WIN64)
	printf("SRC: [%Iu,%Iu,%Iu]\n",SRCI,SRCJ,SRCK);
#else
    printf("SRC: [%zu,%zu,%zu]\n",SRCI,SRCJ,SRCK);
#endif

    ex1_init(V, NI, NJ, NK);
    
    OpenST_AOP_GetArrStats(V, NI * NJ * NK, &vmin, &vmax, &vmean);
    
    BRT3D_TSTEP = OpenST_BRT3D_SuggestTSTEP(vmax, H, H, H);
    EPS = BRT3D_TSTEP;
    
    printf("V[min;max;mean]: [%e;%e;%e]\n",vmin,vmax,vmean);

#ifdef LSM3D_IMP
    IMP_NAME = OPENST_LSM3D_IMP_NAME;
    OpenST_LSM3D_Init(U,LSM_UNLOCKED, NI, NJ, NK, SRCI,SRCJ,SRCK);
    t1 = omp_get_wtime();
    it = OpenST_LSM3D_Compute(U,LSM_UNLOCKED,V,H,NI,NJ,NK,SRCI,SRCJ,SRCK,max_iter,&converged,BSIZE_I,BSIZE_J,BSIZE_K,EPS);
    t2 = omp_get_wtime();
    LSM3D_Time = t2 - t1;
#else
    IMP_NAME = OPENST_FSM3D_IMP_NAME;
    OpenST_FSM3D_Init(U,NI,NJ,NK,SRCI,SRCJ,SRCK);
    t1 = omp_get_wtime();
    it = OpenST_FSM3D_Compute(U,V,H,NI,NJ,NK,SRCI,SRCJ,SRCK,max_iter,&converged,BSIZE_I,BSIZE_J,BSIZE_K,EPS);
    t2 = omp_get_wtime();
    LSM3D_Time = t2 - t1;
#endif
   
    ex1_check(U, H, NI, NJ, NK, SRCI,SRCJ,SRCK, &L1, &L2, &Linf);

    t1 = omp_get_wtime();
    BRT3D_ERR = OpenST_BRT3D_Trace(U, V, NI, NJ, NK, H, H, H, BRT3D_TSTEP,
                        (double)RCVI*H, (double)RCVJ*H, (double)RCVK*H,
                        (double)SRCI*H, (double)SRCJ*H, (double)SRCK*H,
                        &RAY, &RAY_NI, &RAY_NJ);
    t2 = omp_get_wtime();
    BRT3D_Time = t2 - t1;

    printf("====================================================\n");
    printf("TEST_ID,METHOD_ID,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,EPS,max_iter,it,converged,sec,L1(MAE),L2(MSE),L_inf,NI,NJ,NK,SRCI,SRCJ,SRCK,RCVI,RCVJ,RCVK,BRT3D_TSTEP,BRT3D_sec,BRT3D_ERR,OPENST_BUILDINFO_LINK_TYPE\n");

#ifdef _WIN32
	printf("%s,%s,%i,%Iu,%Iu,%Iu,%e,%i,%i,%i,%e,%e,%e,%e,%Iu,%Iu,%Iu,%Iu,%Iu,%Iu,%Iu,%Iu,%Iu,%e,%e,%i,%s\n",
#else
    printf("%s,%s,%i,%zu,%zu,%zu,%e,%i,%i,%i,%e,%e,%e,%e,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%e,%e,%i,%s\n",
#endif
      TEST_NAME,IMP_NAME,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,EPS,max_iter,it,converged,LSM3D_Time,L1,L2,Linf,NI,NJ,NK,SRCI,SRCJ,SRCK,RCVI,RCVJ,RCVK,BRT3D_TSTEP,BRT3D_Time,BRT3D_ERR,OPENST_BUILDINFO_LINK_TYPE_STATIC
    );
    
    return EXIT_SUCCESS;
}
