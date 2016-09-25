#include "openst/eikonal/fsm.h"

#define M_FSM3D_IMP_NAME "fsm3d_bfsm_openmp_v1.c"


const char OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME[] = M_FSM3D_IMP_NAME;
const size_t OPENST_FSM3D_COMPUTEPARTIAL_IMP_NAME_LENGTH = sizeof(M_FSM3D_IMP_NAME);


int OpenST_FSM3D_ComputePartial(double *U, double *V,
                                size_t NI, size_t NJ, size_t NK,
                                double HI, double HJ, double HK,
                                int start_iter, int max_iter, int *converged,
                                size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                double EPS){

    int total_it, it, notconvergedl;
    int REVI, REVJ, REVK;
    size_t ir, jr, kr;
    size_t NBI, NBJ, NBK;
    int *notconvergedt;
    int ***U3d;

    if(start_iter >= max_iter){
        return max_iter;
    }

    total_it = start_iter;
    notconvergedl = 0;

    NBI = NI/BSIZE_I + (NI % BSIZE_I > 0);
    NBJ = NJ/BSIZE_J + (NJ % BSIZE_J > 0);
    NBK = NK/BSIZE_K + (NK % BSIZE_K > 0);

    notconvergedt = (int *)malloc(sizeof(int) * NBI * NBJ * NBK);
    U3d = (int ***)malloc(sizeof(int **) * NBI);
    for(ir = 0; ir < NBI; ++ir){
        U3d[ir] = (int **)malloc(sizeof(int *) * NBJ);
        for(jr = 0; jr < NBJ; ++jr){
            U3d[ir][jr] = &notconvergedt[ir * NBJ * NBK + jr * NBK];
        }
    }

#pragma omp parallel default(none) \
    shared(BSIZE_I, BSIZE_J, BSIZE_K, NBI, NBJ, NBK, total_it, notconvergedl, \
    NI, NJ, NK, ir, jr, kr, \
    U, V, HI, HJ, HK, start_iter, max_iter, notconvergedt, REVI, REVJ, REVK, U3d, \
    EPS) \
    private(it)
    {

        for(it = start_iter; it < max_iter; ++it){

#pragma omp single nowait
            {
                ++total_it;
                notconvergedl = 0;

                OpenST_FSM3D_GetSweepOrder(it, &REVI, &REVJ, &REVK);

#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(out: U3d[0:1][0:1][0:1])
                notconvergedt[0] =
                        OpenST_FSM3D_BlockSerial(U, V,
                                                 NI, NJ, NK,
                                                 HI, HJ, HK,
                                                 REVI, REVJ, REVK,
                                                 0, 0, 0,
                                                 BSIZE_I, BSIZE_J, BSIZE_K, EPS);


                for(ir = 1; ir < NBI; ++ir){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][0 : 1][0: 1]) \
    depend(out: U3d[ir : 1][0 : 1][0 : 1])
                    notconvergedt[ir * NBJ * NBK] =
                            OpenST_FSM3D_BlockSerial(U, V,
                                                     NI, NJ, NK,
                                                     HI, HJ, HK,
                                                     REVI, REVJ, REVK,
                                                     ir * BSIZE_I, 0, 0,
                                                     BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                }

                for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[0 : 1][0 : 1][(kr - 1) : 1]) \
    depend(out: U3d[0 : 1][0 : 1][kr : 1])
                    notconvergedt[kr] =
                            OpenST_FSM3D_BlockSerial(U, V,
                                                     NI, NJ, NK,
                                                     HI, HJ, HK,
                                                     REVI, REVJ, REVK,
                                                     0, 0, kr * BSIZE_K,
                                                     BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                }


                for(ir = 1; ir < NBI; ++ir){
                    for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][0 : 1][kr : 1]) \
    depend(in: U3d[ir : 1][0 : 1][(kr - 1) : 1]) \
    depend(out: U3d[ir : 1][0 : 1][kr : 1])
                        notconvergedt[ir * NBJ * NBK + kr] =
                                OpenST_FSM3D_BlockSerial(U, V,
                                                         NI, NJ, NK,
                                                         HI, HJ, HK,
                                                         REVI, REVJ, REVK,
                                                         ir * BSIZE_I, 0,
                                                         kr * BSIZE_K,
                                                         BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                    }
                }

                for(jr = 1; jr < NBJ; ++jr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[0 : 1][(jr - 1) : 1][0 : 1]) \
    depend(out: U3d[0 : 1][jr  : 1][0 : 1])
                    notconvergedt[jr * NBK] =
                            OpenST_FSM3D_BlockSerial(U, V,
                                                     NI, NJ, NK,
                                                     HI, HJ, HK,
                                                     REVI, REVJ, REVK,
                                                     0, jr * BSIZE_J, 0,
                                                     BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                }

                for(jr = 1; jr < NBJ; ++jr){
                    for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[0 : 1][(jr - 1) : 1][kr : 1]) \
    depend(in: U3d[0 : 1][jr : 1][(kr - 1) : 1]) \
    depend(out: U3d[0 : 1][jr : 1][kr : 1])
                        notconvergedt[jr * NBK + kr] =
                                OpenST_FSM3D_BlockSerial(U, V,
                                                         NI, NJ, NK,
                                                         HI, HJ, HK,
                                                         REVI, REVJ, REVK,
                                                         0, jr * BSIZE_J,
                                                         kr * BSIZE_K,
                                                         BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                    }
                }

                for(ir = 1; ir < NBI; ++ir){
                    for(jr = 1; jr < NBJ; ++jr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][jr : 1][0 : 1]) \
    depend(in: U3d[ir : 1][(jr - 1) : 1][0 : 1]) \
    depend(out: U3d[ir : 1][jr : 1][0 : 1])
                        notconvergedt[ir * NBJ * NBK + jr * NBK] =
                                OpenST_FSM3D_BlockSerial(U, V,
                                                         NI, NJ, NK,
                                                         HI, HJ, HK,
                                                         REVI, REVJ, REVK,
                                                         ir * BSIZE_I, jr * BSIZE_J,
                                                         0,
                                                         BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                    }
                }

                for(ir = 1; ir < NBI; ++ir){
                    for(jr = 1; jr < NBJ; ++jr){
                        for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][jr : 1][kr : 1]) \
    depend(in: U3d[ir : 1][(jr - 1) : 1][kr : 1]) \
    depend(in: U3d[ir : 1][jr : 1][(kr - 1) : 1]) \
    depend(out: U3d[ir : 1][jr : 1][kr : 1])
                            notconvergedt[ir * NBJ * NBK + jr * NBK + kr] =
                                    OpenST_FSM3D_BlockSerial(U, V,
                                                             NI, NJ, NK,
                                                             HI, HJ, HK,
                                                             REVI, REVJ, REVK,
                                                             ir * BSIZE_I, jr * BSIZE_J,
                                                             kr * BSIZE_K,
                                                             BSIZE_I, BSIZE_J, BSIZE_K, EPS);
                        }
                    }
                }
            }

#pragma omp taskwait
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

    for(ir = 0; ir < NBI; ++ir){
        free(U3d[ir]);
    }
    free(U3d);
    free(notconvergedt);

    return total_it;
}
