#include "openst/eikonal/fsm.h"

#define M_FSM3D_IMP_NAME "BFSMv1"


const char OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME[] = M_FSM3D_IMP_NAME;
const size_t OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME_LENGTH = sizeof(M_FSM3D_IMP_NAME);


int OpenST_FSM3D_ComputePartial(double *U, double *V,
                                size_t NI, size_t NJ, size_t NK,
                                double HI, double HJ, double HK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconvergedl, notconvergedt;
    int REVI, REVJ, REVK;
    size_t NBI, NBJ, NBK;

#if (_OPENMP > 200203)
    size_t levelr, K1, K2, kr, level, I1, I2, ir, jr;
#else
	#pragma message("WARNING: size_t to ptrdiff_t cast enabled")
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
    start_iter, total_it, notconvergedl, NI, NJ, NK, \
    U, V, HI, HJ, HK, max_iter, EPS) \
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

                        if(OpenST_FSM3D_BlockSerial(U, V,
                                                    NI, NJ, NK,
                                                    HI, HJ, HK,
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
