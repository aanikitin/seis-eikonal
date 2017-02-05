#include "openst/eikonal/fsm.h"


OPENST_ERR OpenST_FSM3D(double *U, double *V,
                        size_t NI, size_t NJ, size_t NK,
                        double HI, double HJ, double HK,
                        double SRCI, double SRCJ, double SRCK,
                        double EPS, int max_iter,
                        int *it, int *converged){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t BSIZE_I, BSIZE_J, BSIZE_K;

    OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);

    if((errcode = OpenST_FSM3D_Init(U,V,
                      NI,NJ,NK,
                      HI,HJ,HK,
                      SRCI,SRCJ,SRCK)) != OPENST_ERR_SUCCESS){
        goto EXIT;
    }

    *it = OpenST_FSM3D_Compute(U,V,
                              NI,NJ,NK,
                              HI,HJ,HK,
                              max_iter,converged,
                              BSIZE_I,BSIZE_J,BSIZE_K,EPS);

EXIT:
    return errcode;
}


void OpenST_FSM3D_SuggestBlockSize(size_t NI, size_t NJ, size_t NK, 
                                   size_t *BSIZE_I, size_t *BSIZE_J,
                                   size_t *BSIZE_K){
    /* TODO: heuristic to determine optimal block size
     * (based on profiling maybe) */
    *BSIZE_I = 1;
    *BSIZE_J = 10;
    *BSIZE_K = NK;
}


inline void OpenST_FSM3D_GetSweepOrder(int it, int *REVI, int *REVJ, int *REVK){
    int order = it % 8;
    if(order & 4){
        *REVI = 1;
    } else {
        *REVI = 0;
    }
    if(order & 2){
        *REVJ = 1;
    } else {
        *REVJ = 0;
    }
    if(order & 1){
        *REVK = 1;
    } else {
        *REVK = 0;
    }
}


int OpenST_FSM3D_Compute(double *U, double *V,
                         size_t NI, size_t NJ, size_t NK,
                         double HI, double HJ, double HK,
                         int max_iter, int *converged,
                         size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                         double EPS){
    return OpenST_FSM3D_ComputePartial(U, V,
                                       NI, NJ, NK,
                                       HI, HJ, HK,
                                       0, max_iter, converged,
                                       BSIZE_I, BSIZE_J, BSIZE_K, EPS);
}
