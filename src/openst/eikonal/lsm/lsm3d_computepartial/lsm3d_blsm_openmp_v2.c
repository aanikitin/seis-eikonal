#include "openst/eikonal/lsm.h"

#define M_LSM3D_IMP_NAME "BLSMv2"


const char OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME[] = M_LSM3D_IMP_NAME;
const size_t OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME_LENGTH = sizeof(M_LSM3D_IMP_NAME);

#include <pthread.h>
#include <semaphore.h>
#include <errno.h>

sem_t *sem;

int OpenST_LSM3D_ComputePartial(double *U, char *LSM_UNLOCKED, double *V,
                                size_t NI, size_t NJ, size_t NK,
                                double HI, double HJ, double HK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconvergedl;
    int REVI, REVJ, REVK;
    int nth, tid, tid_prev, notconvergedt;

    size_t NBI, NBJ, NBK;
    int nth_max, nth_limit;
    size_t bi, bj, bk;

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
    shared(sem, BSIZE_I, BSIZE_J, BSIZE_K, nth, NBI, NBJ, NBK, total_it, notconvergedl, \
    NI, NJ, NK, \
    U, LSM_UNLOCKED, V, HI, HJ, HK, max_iter, EPS, start_iter) \
    private(tid, tid_prev, it, notconvergedt, \
    REVI, REVJ, REVK, \
    bi, bj, bk)
    {

#pragma omp single
        {
            nth = omp_get_num_threads();
            sem = (sem_t *) malloc(sizeof(sem_t) * (nth));
            for(int i = 0; i < nth; ++i){
              sem_init(&sem[i], 0, 0);
            }
        }

        tid = omp_get_thread_num();
        tid_prev = (tid + nth - 1) % nth;

        for(it = start_iter; it < max_iter; ++it){
#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;
            }

            notconvergedt = 0;

            OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

            for(bi = tid; bi < NBI; bi += nth){
                for(bj = 0; bj < NBJ; ++bj){
                    for(bk = 0; bk < NBK; ++bk){
                        if(bi != 0){
                            sem_wait(&sem[tid_prev]);
                        }
                        if(OpenST_LSM3D_BlockSerial(U, LSM_UNLOCKED, V,
                                                    NI, NJ, NK,
                                                    HI, HJ, HK,
                                                    REVI, REVJ, REVK,
                                                    bi * BSIZE_I, bj * BSIZE_J,
                                                    bk * BSIZE_K,
                                                    BSIZE_I, BSIZE_J, BSIZE_K,
                                                    EPS)){
                            notconvergedt = 1;
                        }
                        if(bi != (NBI - 1)){
                            sem_post(&sem[tid]);
                        }
                    }
                }
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

#pragma omp single
        {
            for(int i = 0; i < nth; ++i){
              sem_destroy(&sem[i]);
            }
        }
    }

    *converged = (notconvergedl == 0);
    return total_it;
}
