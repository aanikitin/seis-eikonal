#include "openst/eikonal/lsm.h"

#define M_LSM3D_IMP_NAME "LSM"


const char OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME[] = M_LSM3D_IMP_NAME;
const size_t OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME_LENGTH = sizeof(M_LSM3D_IMP_NAME);


int OpenST_LSM3D_ComputePartial(double *U, char *LSM_UNLOCKED, double *V,
                                size_t NI, size_t NJ, size_t NK,
                                double HI, double HJ, double HK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconverged;
    int REVI, REVJ, REVK;

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
    notconverged = 0;

    for(it = start_iter; it < max_iter; ++it){

        ++total_it;

        OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

        notconverged = OpenST_LSM3D_BlockSerial(U, LSM_UNLOCKED, V,
                                                NI, NJ, NK,
                                                HI, HJ, HK,
                                                REVI, REVJ, REVK,
                                                0, 0, 0,
                                                NI, NJ, NK, EPS);

        if(!notconverged){
            break;
        }

    }

    *converged = !notconverged;

    return total_it;
}
