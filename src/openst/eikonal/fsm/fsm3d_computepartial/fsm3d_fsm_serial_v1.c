#include "openst/eikonal/fsm.h"

#define M_FSM3D_IMP_NAME "fsm3d_fsm_serial_v1.c"


const char OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME[] = M_FSM3D_IMP_NAME;
const size_t OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME_LENGTH = sizeof(M_FSM3D_IMP_NAME);


int OpenST_FSM3D_ComputePartial(double *U, double *V,
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

        notconverged = OpenST_FSM3D_BlockSerial(U, V,
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
