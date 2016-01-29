/** @file fsm_openmp_v2.c
 *  @brief FSM3D OpenMP first proposed parallel implementation using
 *  static scheduling.
 *
 *  @author Alexandr Nikitin
 */

#include "fsm.h"

#define DEBUG_LOG 0

extern size_t BSIZE_I;
extern size_t BSIZE_J;
extern size_t BSIZE_K;


int FSM3D(double *U, double *F, double H, size_t NI, size_t NJ, size_t NK,
          size_t SRCI, size_t SRCJ, size_t SRCK, int max_iter, int *converged){

    int total_it, it, notconvergedl, order;
    int REVI, REVJ, REVK;
    int nth, tid, notconvergedt;

    size_t NLBI, NBLS, NTS, NBI, NBJ, NBK;
    int tid_last, done, nth_max, nth_limit;
    size_t bls_rem, bi, bj, bk, ts;

    total_it = 0;
    notconvergedl = 0;

    NBI = NI/BSIZE_I + (NI % BSIZE_I > 0);
    NBJ = NJ/BSIZE_J + (NJ % BSIZE_J > 0);
    NBK = NK/BSIZE_K + (NK % BSIZE_K > 0);

    nth_max = omp_get_max_threads();
    if(nth_max > NBI){
        nth_limit = NBI;
    } else {
        nth_limit = nth_max;
    }

#pragma omp parallel num_threads(nth_limit) default(none) \
    shared(BSIZE_I, BSIZE_J, BSIZE_K, nth, NLBI, NBI, NBJ, NBK, total_it, notconvergedl, \
    NI, NJ, NK, SRCI, SRCJ, SRCK, \
    U, F, H, max_iter, NBLS, NTS, tid_last) \
    private(tid, it, notconvergedt, order, \
    REVI, REVJ, REVK, ts, bls_rem, \
    bi, bj, bk, done)
    {

#pragma omp single
        {
            nth = omp_get_num_threads();
            NLBI = NBI/nth + (NBI % nth > 0);
            tid_last = (NBI % nth > 0) ? NBI % nth - 1 : nth - 1;
            NBLS = (NBK * NBJ > nth) ? 0 : (nth - NBK * NBJ);
            NTS = tid_last + NLBI * NBK * NBJ + (NLBI - 1) * NBLS;
        }

        tid = omp_get_thread_num();

        for(it = 0; it < max_iter; ++it){

#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;
            }

            notconvergedt = 0;

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

            bi = tid;
            bj = 0;
            bk = 0;
            bls_rem = tid;
            done = 0;

            for(ts = 0; ts < NTS; ++ts){

                if(!done){
                    if(bls_rem){
                        --bls_rem;
                    } else {

                        if(FSM3D_serial(U, F, H, NI, NJ, NK,
                                        SRCI, SRCJ, SRCK,
                                        REVI, REVJ, REVK,
                                        bi * BSIZE_I, bj * BSIZE_J, bk * BSIZE_K,
                                        BSIZE_I, BSIZE_J, BSIZE_K)){
                            notconvergedt = 1;
                        }

                        //update block index
                        ++bk;
                        if(bk == NBK){
                            bk = 0;
                            ++bj;
                            if(bj == NBJ){
                                bj = 0;
                                bls_rem = NBLS;
                                bi += nth;
                                if(bi >= NBI){
                                    done = 1;
                                }
                            }
                        }

                    }
                }

#pragma omp barrier
            }

#pragma omp atomic
            notconvergedl += notconvergedt;
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
