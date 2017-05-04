#include "openst/eikonal/lsm.h"


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


OPENST_ERR OpenST_LSM3D(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                        size_t NI, size_t NJ, size_t NK,
                        OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                        OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                        OPENST_FLOAT EPS, int max_iter,
                        int *it, int *converged){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t BSIZE_I, BSIZE_J, BSIZE_K;

    OpenST_FSM3D_SuggestBlockSize(NI,NJ,NK,&BSIZE_I,&BSIZE_J,&BSIZE_K);

    if((errcode = OpenST_LSM3D_Init(U,LSM_UNLOCKED,V,
                      NI,NJ,NK,
                      HI,HJ,HK,
                      SRCI,SRCJ,SRCK)) != OPENST_ERR_SUCCESS){
        goto EXIT;
    }

    *it = OpenST_LSM3D_Compute(U,LSM_UNLOCKED,V,
                              NI,NJ,NK,
                              HI,HJ,HK,
                              max_iter,converged,
                              BSIZE_I,BSIZE_J,BSIZE_K,EPS);

EXIT:
    return errcode;
}


int OpenST_LSM3D_Compute(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                         size_t NI, size_t NJ, size_t NK,
                         OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                         int max_iter, int *converged,
                         size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                         OPENST_FLOAT EPS){
    return OpenST_LSM3D_ComputePartial(U, LSM_UNLOCKED, V,
                                       NI, NJ, NK,
                                       HI, HJ, HK,
                                       0, max_iter, converged,
                                       BSIZE_I, BSIZE_J, BSIZE_K, EPS);
}
