/** @file fsm_openmp_v1.c
 *  @brief FSM3D serial implementation of the [Zhao, 2005] original FSM
 *  algorithm.
 *
 *  References:
 *  1) Zhao, Hongkai. "A fast sweeping method for eikonal equations."
 *  Mathematics of computation 74.250 (2005): 603-627.
 *
 *  @author Alexandr Nikitin
 */

#include "fsm.h"

int FSM3D(double *U, double *F, double H, size_t NI, size_t NJ, size_t NK,
    size_t SRCI, size_t SRCJ, size_t SRCK, int max_iter, int *converged){

    int total_it, it, convergedl, order;
    int REVI, REVJ, REVK;

    total_it = 0;
    for(it = 0; it < max_iter; ++it){
        ++total_it;
        convergedl = 1;
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

        convergedl = !(FSM3D_serial(U, F, H, NI, NJ, NK, 
                                     SRCI, SRCJ, SRCK,
                                     REVI, REVJ, REVK,
                                     0, 0, 0,
                                     NI, NJ, NK));
        
        if(convergedl){
            break;
        }
    }
    *converged = convergedl;
    return total_it;
}
