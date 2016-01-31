/** @file fsm_openmp_v5.c
 *  @brief FSM3D OpenMP version similar to fsm_openmp_v1.c but with blocks
 *  in OpenMP tasks.
 *
 *  @author Alexandr Nikitin
 */


#include "fsm.h"
#include "stdlib.h"

#define DEBUG_LOG 0

extern size_t BSIZE_I;
extern size_t BSIZE_J;
extern size_t BSIZE_K;


int FSM3D(double *U, double *F, double H, size_t NI, size_t NJ, size_t NK,
          size_t SRCI, size_t SRCJ, size_t SRCK, int max_iter, int *converged){

    int total_it, it, notconvergedl, order;
    int REVI, REVJ, REVK;
    size_t levelr, K1, K2, kr, level, I1, I2, ir, jr;
    size_t NBI, NBJ, NBK;

    total_it = 0;
    notconvergedl = 0;

    NBI = NI/BSIZE_I + (NI % BSIZE_I > 0);
    NBJ = NJ/BSIZE_J + (NJ % BSIZE_J > 0);
    NBK = NK/BSIZE_K + (NK % BSIZE_K > 0);
    
    int *notconvergedt = malloc(sizeof(int) * NBI * NBJ * NBK);

#pragma omp parallel default(none) \
    shared(BSIZE_I, BSIZE_J, BSIZE_K, NBI, NBJ, NBK, total_it, notconvergedl, \
    NI, NJ, NK, SRCI, SRCJ, SRCK, levelr, K1, K2, ir, jr, kr, level, I1, I2, \
    U, F, H, max_iter, order, REVI, REVJ, REVK, notconvergedt) \
    private(it)
    {

        for(it = 0; it < max_iter; ++it){

#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;

                order = it % 8;
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

                for(levelr = 0; levelr < NBI + NBJ + NBK - 2; ++levelr){

                    K1 = (NBI + NBJ - 2 < levelr) ?
                                (levelr - NBI - NBJ + 2) : 0;
                    K2 = (NBK - 1 > levelr) ? levelr : NBK - 1;

                    for(kr = K1; kr <= K2; ++kr){
                        level = levelr - kr;

                        I1 = (NBJ - 1 < level) ? (level - NBJ + 1) : 0;
                        I2 = (NBI - 1 > level) ? level : NBI - 1;

                        for(ir = I1; ir <= I2; ++ir){
                            jr = level - ir;

#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK)
                            {
                                notconvergedt[ir * NBJ * NBK + jr * NBK + kr] = FSM3D_serial(U, F, H, NI, NJ, NK,
                                                              SRCI, SRCJ, SRCK,
                                                              REVI, REVJ, REVK,
                                                              ir * BSIZE_I, jr * BSIZE_J,
                                                              kr * BSIZE_K,
                                                              BSIZE_I, BSIZE_J, BSIZE_K);
                            }
                        }
                    }
#pragma omp taskwait
                }
            }

#pragma omp barrier
#pragma omp for reduction(+:notconvergedl)
	    for(ir = 0; ir < NBI * NBJ * NBK; ++ir){
		notconvergedl += notconvergedt[ir];
	    }
	    #pragma omp barrier
	    #pragma omp flush(notconvergedl)
            if(!notconvergedl){
                break;
            }
#pragma omp barrier
        }
    }

    *converged = (notconvergedl == 0);
    return total_it;
}
