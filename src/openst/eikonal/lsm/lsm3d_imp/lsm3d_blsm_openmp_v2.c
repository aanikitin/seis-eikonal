#include "openst/eikonal/lsm.h"

#define M_LSM3D_IMP_NAME "lsm3d_blsm_openmp_v2.c"


const char OPENST_LSM3D_IMP_NAME[] = M_LSM3D_IMP_NAME;
const size_t OPENST_LSM3D_IMP_NAME_LENGTH = sizeof(M_LSM3D_IMP_NAME);


int OpenST_LSM3D_ComputePartial(double *U, char *LSM_UNLOCKED, double *V, double H,
                                size_t NI, size_t NJ, size_t NK,
                                size_t SRCI, size_t SRCJ, size_t SRCK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconvergedl, notconvergedt;
    int REVI, REVJ, REVK;
    size_t NBI, NBJ, NBK;

#if (_OPENMP > 200203)
    size_t levelr, K1, K2, kr, level, I1, I2, ir, jr;
#else
#warning size_t downcast to ptrdiff_t due to missing OpenMP 3.0 support
    ptrdiff_t levelr, K1, K2, kr, level, I1, I2, ir, jr;
#endif

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
    notconvergedl = 0;
    
    NBI = NI/BSIZE_I + (NI % BSIZE_I > 0);
    NBJ = NJ/BSIZE_J + (NJ % BSIZE_J > 0);
    NBK = NK/BSIZE_K + (NK % BSIZE_K > 0);

#pragma omp parallel default(none) \
    shared(BSIZE_I, BSIZE_J, BSIZE_K, NBI, NBJ, NBK, \
    start_iter, total_it, notconvergedl, NI, NJ, NK, SRCI, SRCJ, SRCK, \
    U, LSM_UNLOCKED, V, H, max_iter, EPS) \
    private(it, REVI, REVJ, REVK, notconvergedt, \
    levelr, K1, K2, level, I1, I2, ir, jr, kr)
    {
        for(it = start_iter; it < max_iter; ++it){
#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;
            }

            notconvergedt = 0;

            OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

            for(levelr = 0; levelr < NBI + NBJ + NBK - 2; ++levelr){

                K1 = (NBI + NBJ - 2 < levelr) ?
                            (levelr - NBI - NBJ + 2) : 0;
                K2 = (NBK - 1 > levelr) ? levelr : NBK - 1;

                for(kr = K1; kr <= K2; ++kr){
                    level = levelr - kr;

                    I1 = (NBJ - 1 < level) ? (level - NBJ + 1) : 0;
                    I2 = (NBI - 1 > level) ? level : NBI - 1;

#pragma omp for nowait schedule(dynamic,1)
                    for(ir = I1; ir <= I2; ++ir){
                        jr = level - ir;

                        if(OpenST_LSM3D_BlockSerial(U, LSM_UNLOCKED, V, H,
                                                    NI, NJ, NK,
                                                    SRCI, SRCJ, SRCK,
                                                    REVI, REVJ, REVK,
                                                    ir * BSIZE_I, jr * BSIZE_J,
                                                    kr * BSIZE_K,
                                                    BSIZE_I, BSIZE_J, BSIZE_K,
                                                    EPS)){
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


int OpenST_LSM3D_Compute(double *U, char *LSM_UNLOCKED, double *V, double H,
                         size_t NI, size_t NJ, size_t NK,
                         size_t SRCI, size_t SRCJ, size_t SRCK,
                         int max_iter, int *converged,
                         size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                         double EPS){
    return OpenST_LSM3D_ComputePartial(U, LSM_UNLOCKED, V, H,
                                       NI, NJ, NK,
                                       SRCI, SRCJ, SRCK,
                                       0, max_iter, converged,
                                       BSIZE_I, BSIZE_J, BSIZE_K, EPS);
}
