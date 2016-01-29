#include "fsm.h"
#include <float.h>


void FSM3DInit(double *U, size_t NI, size_t NJ, size_t NK,
               size_t SRCI, size_t SRCJ, size_t SRCK){
  size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                if(i == SRCI && j == SRCJ && k == SRCK){
                    U[M_MEM_ADR_3D(i,j,k,NI,NJ,NK)] = 0.0;
                } else {
                    U[M_MEM_ADR_3D(i,j,k,NI,NJ,NK)] = DBL_MAX;
                }
            }
        }
    }
}


int FSM3D_node_update(double *U, double *F, double H,
                      size_t NI, size_t NJ, size_t NK,
                      size_t SRCI, size_t SRCJ, size_t SRCK,
                      int REVI, int REVJ, int REVK,
                      size_t ir, size_t jr, size_t kr){

    size_t i,j,k;
    size_t mem_cur, mem_il, mem_ir, mem_jl, mem_jr, mem_kl, mem_kr;
    double uold, uxmin, uymin, uzmin, a1, a2, a3, t1, t2, t3, unew;
    int notconverged;

    notconverged = 0;

    if(REVI){
        i = NI - 1 - ir;
    } else {
        i = ir;
    }
    if(REVJ){
        j = NJ - 1 - jr;
    } else {
        j = jr;
    }
    if(REVK){
        k = NK - 1 - kr;
    } else {
        k = kr;
    }

    if(i == SRCI && j == SRCJ && k == SRCK){
        return notconverged;
    }

    mem_cur = M_MEM_ADR_3D(i,j,k,NI,NJ,NK);
    mem_il = M_MEM_ADR_3D(i - 1,j,k,NI,NJ,NK);
    mem_ir = M_MEM_ADR_3D(i + 1,j,k,NI,NJ,NK);
    mem_jl = M_MEM_ADR_3D(i,j - 1,k,NI,NJ,NK);
    mem_jr = M_MEM_ADR_3D(i,j + 1,k,NI,NJ,NK);
    mem_kl = M_MEM_ADR_3D(i,j,k - 1,NI,NJ,NK);
    mem_kr = M_MEM_ADR_3D(i,j,k + 1,NI,NJ,NK);

    uold = U[mem_cur];

    if(j == 0){
        uxmin = U[mem_jr];
    } else if (j == NJ - 1){
        uxmin = U[mem_jl];
    } else {
        uxmin = fmin(U[mem_jl], U[mem_jr]);
    }

    if(i == 0){
        uymin = U[mem_ir];
    } else if (i == NI - 1){
        uymin = U[mem_il];
    } else {
        uymin = fmin(U[mem_il], U[mem_ir]);
    }

    if(k == 0){
        uzmin = U[mem_kr];
    } else if (k == NK - 1){
        uzmin = U[mem_kl];
    } else {
        uzmin = fmin(U[mem_kl], U[mem_kr]);
    }

    if (uxmin <= uymin && uxmin <= uzmin){
        a1 = uxmin;
        if (uymin <= uzmin){
            a2 = uymin;
            a3 = uzmin;
        } else {
            a2 = uzmin;
            a3 = uymin;
        }
    } else if (uymin <= uxmin && uymin <= uzmin){
        a1 = uymin;
        if (uxmin <= uzmin){
            a2 = uxmin;
            a3 = uzmin;
        } else {
            a2 = uzmin;
            a3 = uxmin;
        }
    } else {
        a1 = uzmin;
        if (uxmin <= uymin){
            a2 = uxmin;
            a3 = uymin;
        } else {
            a2 = uymin;
            a3 = uxmin;
        }
    }

    t1 = a1 - a2;
    t2 = F[mem_cur] * H;
    unew = a1 + t2;
    if (unew > a2){
        unew = (a1 + a2 + sqrt(2 * t2 * t2 - t1 * t1))/2;
        if (unew > a3){
            t3 = a1 + a2 + a3;
            unew = (t3 + sqrt(t3 * t3  - 3 *
                              (a1 * a1 + a2 * a2 + a3 * a3
                               - t2 * t2)))/3;
        }
    }

    if(unew < uold){
        notconverged = 1;
        U[mem_cur] = unew;
    }

    return notconverged;

}


int FSM3D_serial(double *U, double *F, double H,
                 size_t NI, size_t NJ, size_t NK,
                 size_t SRCI, size_t SRCJ, size_t SRCK,
                 int REVI, int REVJ, int REVK,
                 size_t istart, size_t jstart, size_t kstart,
                 size_t isize, size_t jsize, size_t ksize){

    int notconverged;
    size_t ir, jr, kr;

    notconverged = 0;
    for(ir = istart; ir < istart + isize && ir < NI; ++ir){
        for(jr = jstart; jr < jstart + jsize && jr < NJ; ++jr){
            for(kr = kstart; kr < kstart + ksize && kr < NK; ++kr){
                notconverged = FSM3D_node_update(U, F, H, NI, NJ, NK,
                    SRCI, SRCJ, SRCK, REVI, REVJ, REVK, ir, jr, kr);
            }
        }
    }

    return notconverged;
}
