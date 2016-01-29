/** @file fsm_openmp_v3.c
 *  @brief FSM3D OpenMP second proposed parallel implementation using
 *  OpenMP 4.0 task depend clause.
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
    size_t ir, jr, kr;
    size_t NBI, NBJ, NBK;

    total_it = 0;
    notconvergedl = 0;

    NBI = NI/BSIZE_I + (NI % BSIZE_I > 0);
    NBJ = NJ/BSIZE_J + (NJ % BSIZE_J > 0);
    NBK = NK/BSIZE_K + (NK % BSIZE_K > 0);

    int *notconvergedt = malloc(sizeof(int) * NBI * NBJ * NBK);
    int ***U3d = (int ***)malloc(sizeof(int **) * NBI);
    for(ir = 0; ir < NBI; ++ir){
        U3d[ir] = (int **)malloc(sizeof(int *) * NBJ);
        for(jr = 0; jr < NBJ; ++jr){
            U3d[ir][jr] = &notconvergedt[ir * NBJ * NBK + jr * NBK];
        }
    }

#pragma omp parallel default(none) \
    shared(BSIZE_I, BSIZE_J, BSIZE_K, NBI, NBJ, NBK, total_it, notconvergedl, \
    NI, NJ, NK, SRCI, SRCJ, SRCK, ir, jr, kr, \
    U, F, H, max_iter, notconvergedt, order, REVI, REVJ, REVK, U3d) \
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

#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(out: U3d[0:1][0:1][0:1])
                notconvergedt[0] =
                        FSM3D_serial(U, F, H, NI, NJ, NK,
                                     SRCI, SRCJ, SRCK,
                                     REVI, REVJ, REVK,
                                     0, 0, 0,
                                     BSIZE_I, BSIZE_J, BSIZE_K);


                for(ir = 1; ir < NBI; ++ir){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][0 : 1][0: 1]) \
    depend(out: U3d[ir : 1][0 : 1][0 : 1])
                    notconvergedt[ir * NBJ * NBK] =
                            FSM3D_serial(U, F, H, NI, NJ, NK,
                                         SRCI, SRCJ, SRCK,
                                         REVI, REVJ, REVK,
                                         ir * BSIZE_I, 0, 0,
                                         BSIZE_I, BSIZE_J, BSIZE_K);
                }

                for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[0 : 1][0 : 1][(kr - 1) : 1]) \
    depend(out: U3d[0 : 1][0 : 1][kr : 1])
                    notconvergedt[kr] =
                            FSM3D_serial(U, F, H, NI, NJ, NK,
                                         SRCI, SRCJ, SRCK,
                                         REVI, REVJ, REVK,
                                         0, 0, kr * BSIZE_K,
                                         BSIZE_I, BSIZE_J, BSIZE_K);
                }


                for(ir = 1; ir < NBI; ++ir){
                    for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][0 : 1][kr : 1]) \
    depend(in: U3d[ir : 1][0 : 1][(kr - 1) : 1]) \
    depend(out: U3d[ir : 1][0 : 1][kr : 1])
                        notconvergedt[ir * NBJ * NBK + kr] =
                                FSM3D_serial(U, F, H, NI, NJ, NK,
                                             SRCI, SRCJ, SRCK,
                                             REVI, REVJ, REVK,
                                             ir * BSIZE_I, 0,
                                             kr * BSIZE_K,
                                             BSIZE_I, BSIZE_J, BSIZE_K);
                    }
                }

                for(jr = 1; jr < NBJ; ++jr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[0 : 1][(jr - 1) : 1][0 : 1]) \
    depend(out: U3d[0 : 1][jr  : 1][0 : 1])
                    notconvergedt[jr * NBK] =
                            FSM3D_serial(U, F, H, NI, NJ, NK,
                                         SRCI, SRCJ, SRCK,
                                         REVI, REVJ, REVK,
                                         0, jr * BSIZE_J, 0,
                                         BSIZE_I, BSIZE_J, BSIZE_K);
                }

                for(jr = 1; jr < NBJ; ++jr){
                    for(kr = 1; kr < NBK; ++kr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[0 : 1][(jr - 1) : 1][kr : 1]) \
    depend(in: U3d[0 : 1][jr : 1][(kr - 1) : 1]) \
    depend(out: U3d[0 : 1][jr : 1][kr : 1])
                        notconvergedt[jr * NBK + kr] =
                                FSM3D_serial(U, F, H, NI, NJ, NK,
                                             SRCI, SRCJ, SRCK,
                                             REVI, REVJ, REVK,
                                             0, jr * BSIZE_J,
                                             kr * BSIZE_K,
                                             BSIZE_I, BSIZE_J, BSIZE_K);
                    }
                }

                for(ir = 1; ir < NBI; ++ir){
                    for(jr = 1; jr < NBJ; ++jr){
#pragma omp task default(shared) firstprivate(ir, jr, kr, REVI, REVJ, REVK) \
    depend(in: U3d[(ir - 1) : 1][jr : 1][0 : 1]) \
    depend(in: U3d[ir : 1][(jr - 1) : 1][0 : 1]) \
    depend(out: U3d[ir : 1][jr : 1][0 : 1])
                        notconvergedt[ir * NBJ * NBK + jr * NBK] =
                                FSM3D_serial(U, F, H, NI, NJ, NK,
                                             SRCI, SRCJ, SRCK,
                                             REVI, REVJ, REVK,
                                             ir * BSIZE_I, jr * BSIZE_J,
                                             0,
                                             BSIZE_I, BSIZE_J, BSIZE_K);
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
                                    FSM3D_serial(U, F, H, NI, NJ, NK,
                                                 SRCI, SRCJ, SRCK,
                                                 REVI, REVJ, REVK,
                                                 ir * BSIZE_I, jr * BSIZE_J,
                                                 kr * BSIZE_K,
                                                 BSIZE_I, BSIZE_J, BSIZE_K);
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
    return total_it;
}
