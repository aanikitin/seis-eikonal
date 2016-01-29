/** @file fsm_openmp_v1.c
 *  @brief FSM3D OpenMP implementation of the [Detrixhe et al, 2013] parallel
 *         algorithm.
 *
 *  References:
 *  1) Detrixhe, Miles, Frédéric Gibou, and Chohong Min. "A parallel fast
 *  sweeping method for the Eikonal equation." Journal of Computational Physics
 *  237 (2013): 46-55.
 *
 *  @author Alexandr Nikitin
 */

#include "fsm.h"

int FSM3D(double *U, double *F, double H, size_t NI, size_t NJ, size_t NK,
          size_t SRCI, size_t SRCJ, size_t SRCK, int max_iter, int *converged){

    int total_it, it, notconvergedl, notconvergedt, order;
    int REVI, REVJ, REVK;
    size_t levelr, K1, K2, kr, level, I1, I2, ir, jr;

    total_it = 0;

#pragma omp parallel default(none) \
    shared(total_it, notconvergedl, NI, NJ, NK, SRCI, SRCJ, SRCK, \
    U, F, H, max_iter) \
    private(it, notconvergedt, order, \
    REVI, REVJ, REVK, levelr, K1, K2, kr, level, I1, I2, ir, jr)
    {
        for(it = 0; it < max_iter; ++it){
#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;
            }

            order = it % 8;
            notconvergedt = 0;

            if(order & 4){
                REVI = 1;
            } else {
                REVI = 0;
            }
            if(order & 2){
                REVJ = 1;
            } else {
                REVJ = 0;
            }
            if(order & 1){
                REVK = 1;
            } else {
                REVK = 0;
            }

            for(levelr = 0; levelr < NI + NJ + NK - 2; ++levelr){

                K1 = (NI + NJ - 2 < levelr) ?
                            (levelr - NI - NJ + 2) : 0;
                K2 = (NK - 1 > levelr) ? levelr : NK - 1;

                for(kr = K1; kr <= K2; ++kr){
                    level = levelr - kr;

                    I1 = (NJ - 1 < level) ? (level - NJ + 1) : 0;
                    I2 = (NI - 1 > level) ? level : NI - 1;

#pragma omp for nowait schedule(static,1)
                    for(ir = I1; ir <= I2; ++ir){
                        jr = level - ir;

                        notconvergedt = FSM3D_node_update(U, F, H, NI, NJ, NK,
                                                          SRCI, SRCJ, SRCK, REVI, REVJ, REVK, ir, jr, kr);
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

