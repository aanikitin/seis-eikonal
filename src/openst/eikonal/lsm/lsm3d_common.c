#include "openst/eikonal/lsm.h"


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


OPENST_ERR OpenST_LSM3D(double *U, char *LSM_UNLOCKED, double *V,
                        size_t NI, size_t NJ, size_t NK,
                        double HI, double HJ, double HK,
                        double SRCI, double SRCJ, double SRCK,
                        double EPS, int max_iter,
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


int OpenST_LSM3D_Compute(double *U, char *LSM_UNLOCKED, double *V,
                         size_t NI, size_t NJ, size_t NK,
                         double HI, double HJ, double HK,
                         int max_iter, int *converged,
                         size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                         double EPS){
    return OpenST_LSM3D_ComputePartial(U, LSM_UNLOCKED, V,
                                       NI, NJ, NK,
                                       HI, HJ, HK,
                                       0, max_iter, converged,
                                       BSIZE_I, BSIZE_J, BSIZE_K, EPS);
}
