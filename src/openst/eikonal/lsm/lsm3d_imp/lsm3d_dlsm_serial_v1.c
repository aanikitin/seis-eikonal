#include "openst/eikonal/lsm.h"

#define M_LSM3D_IMP_NAME "lsm3d_dlsm_serial_v1.c"


const char OPENST_LSM3D_IMP_NAME[] = M_LSM3D_IMP_NAME;
const size_t OPENST_LSM3D_IMP_NAME_LENGTH = sizeof(M_LSM3D_IMP_NAME);


int OpenST_LSM3D_ComputePartial(double *U, char *LSM_UNLOCKED, double *V, double H,
                                size_t NI, size_t NJ, size_t NK,
                                size_t SRCI, size_t SRCJ, size_t SRCK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){


    int total_it, it, notconvergedl;
    int REVI, REVJ, REVK;
    size_t levelr, K1, K2, kr, level, I1, I2, ir, jr;

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
    notconvergedl = 0;

    for(it = start_iter; it < max_iter; ++it){

        ++total_it;
        notconvergedl = 0;

        OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

        for(levelr = 0; levelr < NI + NJ + NK - 2; ++levelr){

            K1 = (NI + NJ - 2 < levelr) ?
                        (levelr - NI - NJ + 2) : 0;
            K2 = (NK - 1 > levelr) ? levelr : NK - 1;

            for(kr = K1; kr <= K2; ++kr){
                level = levelr - kr;

                I1 = (NJ - 1 < level) ? (level - NJ + 1) : 0;
                I2 = (NI - 1 > level) ? level : NI - 1;

                for(ir = I1; ir <= I2; ++ir){
                    jr = level - ir;

                    if(OpenST_LSM3D_NodeUpdate(U, LSM_UNLOCKED, V, H,
                                               NI, NJ, NK,
                                               SRCI, SRCJ, SRCK,
                                               REVI, REVJ, REVK,
                                               ir, jr, kr, EPS)){
                        notconvergedl = 1;
                    }
                }
            }
        }
        if(!notconvergedl){
            break;
        }
    }
    *converged = (notconvergedl == 0);
    return total_it;
}


int OpenST_LSM3D_Compute(double *U, char *LSM_UNLOCKED, double *V, double H,
                         size_t NI, size_t NJ, size_t NK,
                         size_t SRCI, size_t SRCJ, size_t SRCK, int max_iter,
                         int *converged, size_t BSIZE_I, size_t BSIZE_J,
                         size_t BSIZE_K, double EPS){
    return OpenST_LSM3D_ComputePartial(U, LSM_UNLOCKED, V, H,
                                       NI, NJ, NK,
                                       SRCI, SRCJ, SRCK,
                                       0, max_iter, converged,
                                       BSIZE_I, BSIZE_J, BSIZE_K, EPS);
}
