#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#define OPENST_LINK_SHARED
#include "openst.h"
#include "fio_hdf5.h"

size_t BSIZE_I;
size_t BSIZE_J;
size_t BSIZE_K;
char *LSM_UNLOCKED;


void app_info(char *BIN_NAME, int usage){
    if(usage){
        printf("\nUsage: %s FIN FOUT [SRCI SRCJ SRCK [BSIZE_I BSIZE_J BSIZE_K]]\n\t[SRCI SRCJ SRCK] "\
               "- SRC indices\n\t[BSIZE_I x BSIZE_J x BSIZE_K] - task size\n",BIN_NAME);
    }
}


int main(int argc, char *argv[]){
    char *filename, *filename_out;
    double *U, *V, H;
    double t1,t2;
    int it, converged;
    int max_iter;
    size_t NI, NJ, NK;
    size_t SRCI, SRCJ, SRCK;
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

    if(argc >= 3){
        filename = argv[1];
        filename_out = argv[2];
        app_info(argv[0],0);
    } else {
        app_info(argv[0],1);
        return EXIT_FAILURE;
    }

    hdf5_read_model(filename, &V, &H, &NI, &NJ, &NK);

    if(argc >= 6){
        SRCI = atoi(argv[3]);
        SRCJ = atoi(argv[4]);
        SRCK = atoi(argv[5]);
    } else {
        SRCI = NI/2;
        SRCJ = NJ/2;
        SRCK = NK/2;
    }
    
    printf("H: %f\n",H);
    printf("V[0]: %f\n",V[0 * NJ * NK + 0 * NK + 0]);
    printf("V[i]: %f\n",V[1 * NJ * NK + 0 * NK + 0]);
    printf("V[j]: %f\n",V[0 * NJ * NK + 1 * NK + 0]);
    printf("V[k]: %f\n",V[0 * NJ * NK + 0 * NK + 1]);
    printf("V[end]: %f\n",V[(NI - 1) * NJ * NK + (NJ - 1) * NK + (NK - 1)]);
    printf("V[i]: %f\n",V[(NI - 2) * NJ * NK + (NJ - 1) * NK + (NK - 1)]);
    printf("V[j]: %f\n",V[(NI - 1) * NJ * NK + (NJ - 2) * NK + (NK - 1)]);
    printf("V[k]: %f\n",V[(NI - 1) * NJ * NK + (NJ - 1) * NK + (NK - 2)]);
    
    if(argc >= 9){
        BSIZE_I = atoi(argv[6]);
        BSIZE_J = atoi(argv[7]);
        BSIZE_K = atoi(argv[8]);
    } else {
        OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);
    }

    U = malloc(NI * NJ * NK * sizeof(double));
    assert(U);
    LSM_UNLOCKED = malloc(NI * NJ * NK * sizeof(char));
    assert(LSM_UNLOCKED);

    max_iter = 1000;
    RCVI = 0;
    RCVJ = 0;
    RCVK = 0;

    OMP_MAX_THREADS = omp_get_max_threads();

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

    //hdf5_write_time(filename_out, U, &H, NI, NJ, NK);

    t1 = omp_get_wtime();
    BRT3D_ERR = OpenST_BRT3D_Trace(U, V, NI, NJ, NK, H, H, H, BRT3D_TSTEP,
                                   (double)RCVI*H, (double)RCVJ*H, (double)RCVK*H,
                                   (double)SRCI*H, (double)SRCJ*H, (double)SRCK*H,
                                   &RAY, &RAY_NI, &RAY_NJ);
    t2 = omp_get_wtime();
    BRT3D_Time = t2 - t1;

    printf("====================================================\n");
    printf("TEST_ID,METHOD_ID,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,EPS,max_iter,it,converged,sec,L1(MAE),L2(MSE),L_inf,NI,NJ,NK,SRCI,SRCJ,SRCK,RCVI,RCVJ,RCVK,BRT3D_TSTEP,BRT3D_sec,BRT3D_ERR,OPENST_BUILDINFO_LINK_TYPE\n");

    printf("%s,%s,%i,%zu,%zu,%zu,%e,%i,%i,%i,%e,%e,%e,%e,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%e,%e,%i,%s\n",
           filename,IMP_NAME,OMP_MAX_THREADS,BSIZE_I,BSIZE_J,BSIZE_K,EPS,max_iter,
           it,converged,LSM3D_Time,0.0,0.0,0.0,NI,NJ,NK,SRCI,SRCJ,SRCK,RCVI,
           RCVJ,RCVK,BRT3D_TSTEP,BRT3D_Time,BRT3D_ERR,
           OPENST_BUILDINFO_LINK_TYPE_STATIC
           );

    return EXIT_SUCCESS;
}
