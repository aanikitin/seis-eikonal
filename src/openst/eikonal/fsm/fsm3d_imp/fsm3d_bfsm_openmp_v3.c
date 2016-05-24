#include "openst/eikonal/fsm.h"

#define M_FSM3D_IMP_NAME "fsm3d_bfsm_openmp_v3.c"


const char OPENST_FSM3D_IMP_NAME[] = M_FSM3D_IMP_NAME;
const size_t OPENST_FSM3D_IMP_NAME_LENGTH = sizeof(M_FSM3D_IMP_NAME);


int OpenST_FSM3D_ComputePartial(double *U, double *V, double H,
                                size_t NI, size_t NJ, size_t NK,
                                size_t SRCI, size_t SRCJ, size_t SRCK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconvergedl;
    int REVI, REVJ, REVK;
    int nth, tid, notconvergedt;

    size_t NLBI, NBLS, NTS, NBI, NBJ, NBK;
    int tid_last, done, nth_max, nth_limit;
    size_t bls_rem, bi, bj, bk, ts;

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
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
    U, V, H, max_iter, NBLS, NTS, tid_last, EPS, start_iter) \
    private(tid, it, notconvergedt, \
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

        for(it = start_iter; it < max_iter; ++it){

#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;
            }

            notconvergedt = 0;

            OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

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
                        
                        if(OpenST_FSM3D_BlockSerial(U, V, H, NI, NJ, NK,
                                                    SRCI, SRCJ, SRCK,
                                                    REVI, REVJ, REVK,
                                                    bi * BSIZE_I, bj * BSIZE_J,
                                                    bk * BSIZE_K,
                                                    BSIZE_I, BSIZE_J, BSIZE_K,
                                                    EPS)){
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

