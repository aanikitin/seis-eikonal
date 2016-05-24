#include "openst/eikonal/fsm.h"

#define M_FSM3D_IMP_NAME "fsm3d_fsm_serial_v1.c"


const char OPENST_FSM3D_IMP_NAME[] = M_FSM3D_IMP_NAME;
const size_t OPENST_FSM3D_IMP_NAME_LENGTH = sizeof(M_FSM3D_IMP_NAME);


int OpenST_FSM3D_ComputePartial(double *U, double *V, double H,
                          size_t NI, size_t NJ, size_t NK,
                          size_t SRCI, size_t SRCJ, size_t SRCK,
                          int start_iter, int max_iter, int *converged,
                          size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                          double EPS){

    int total_it, it, convergedl;
    int REVI, REVJ, REVK;

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
    convergedl = 0;

    for(it = start_iter; it < max_iter; ++it){
        ++total_it;

        OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

        convergedl = !(OpenST_FSM3D_BlockSerial(U, V, H, NI, NJ, NK,
                                          SRCI, SRCJ, SRCK,
                                          REVI, REVJ, REVK,
                                          0, 0, 0,
                                          NI, NJ, NK, EPS));
        
        if(convergedl){
            break;
        }
    }
    *converged = convergedl;
    return total_it;
}


int OpenST_FSM3D_Compute(double *U, double *V, double H,
                  size_t NI, size_t NJ, size_t NK,
                  size_t SRCI, size_t SRCJ, size_t SRCK,
                  int max_iter, int *converged,
                  size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                  double EPS){
    return OpenST_FSM3D_ComputePartial(U, V, H,
                                        NI, NJ, NK,
                                        SRCI, SRCJ, SRCK,
                                        0, max_iter, converged,
                                        BSIZE_I, BSIZE_J, BSIZE_K, EPS);
}
