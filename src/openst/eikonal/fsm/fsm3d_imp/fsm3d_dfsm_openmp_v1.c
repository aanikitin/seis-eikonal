#include "openst/eikonal/fsm.h"

#define M_FSM3D_IMP_NAME "fsm3d_dfsm_openmp_v1.c"


const char OPENST_FSM3D_IMP_NAME[] = M_FSM3D_IMP_NAME;
const size_t OPENST_FSM3D_IMP_NAME_LENGTH = sizeof(M_FSM3D_IMP_NAME);


int OpenST_FSM3D_ComputePartial(double *U, double *V, double H,
                                size_t NI, size_t NJ, size_t NK,
                                size_t SRCI, size_t SRCJ, size_t SRCK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconvergedl, notconvergedt;
    int REVI, REVJ, REVK;
    size_t levelr, K1, K2, kr, level, I1, I2, ir, jr;

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
    notconvergedl = 0;

#pragma omp parallel default(none) \
    shared(total_it, notconvergedl, NI, NJ, NK, SRCI, SRCJ, SRCK, \
    U, V, H, max_iter, EPS, start_iter) \
    private(it, notconvergedt, \
    REVI, REVJ, REVK, levelr, K1, K2, kr, level, I1, I2, ir, jr)
    {
        for(it = start_iter; it < max_iter; ++it){
#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;
            }

            notconvergedt = 0;

            OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

            for(levelr = 0; levelr < NI + NJ + NK - 2; ++levelr){

                K1 = (NI + NJ - 2 < levelr) ?
                            (levelr - NI - NJ + 2) : 0;
                K2 = (NK - 1 > levelr) ? levelr : NK - 1;

                for(kr = K1; kr <= K2; ++kr){
                    level = levelr - kr;

                    I1 = (NJ - 1 < level) ? (level - NJ + 1) : 0;
                    I2 = (NI - 1 > level) ? level : NI - 1;

#pragma omp for nowait
                    for(ir = I1; ir <= I2; ++ir){
                        jr = level - ir;

                        if(OpenST_FSM3D_NodeUpdate(U, V, H, NI, NJ, NK,
                                                          SRCI, SRCJ, SRCK,
                                                          REVI, REVJ, REVK,
                                                          ir, jr, kr, EPS)){
                            notconvergedt = 1;
                        }
                    }
                }
#pragma omp barrier
            }
#pragma omp atomic
            notconvergedl += notconvergedt;
#pragma omp barrier
#pragma omp flush (notconvergedl)
            if(!notconvergedl){
                break;
            }
#pragma omp barrier
        }
    }
    *converged = (notconvergedl == 0);
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
