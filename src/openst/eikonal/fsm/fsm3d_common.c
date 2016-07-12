#include "openst/eikonal/fsm.h"


void OpenST_FSM3D_SuggestBlockSize(size_t NI, size_t NJ, size_t NK, 
                              size_t *BSIZE_I, size_t *BSIZE_J,
                              size_t *BSIZE_K){
    /* TODO: heuristic to determine optimal block size
     * (based on profiling maybe) */
    *BSIZE_I = 1;
    *BSIZE_J = 10;
    *BSIZE_K = NK;
}


void OpenST_FSM3D_GetSweepOrder(int it, int *REVI, int *REVJ, int *REVK){
    int order = it % 8;
    if(order & 4){
        *REVI = 1;
    } else {
        *REVI = 0;
    }
    if(order & 2){
        *REVJ = 1;
    } else {
        *REVJ = 0;
    }
    if(order & 1){
        *REVK = 1;
    } else {
        *REVK = 0;
    }
}


void OpenST_FSM3D_Init(double *U, size_t NI, size_t NJ, size_t NK,
                size_t SRCI, size_t SRCJ, size_t SRCK){
    size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = DBL_MAX;
            }
        }
    }
    U[OPENST_MEMADR_3D(SRCI,SRCJ,SRCK,NI,NJ,NK)] = 0.0;
}


int OpenST_FSM3D_NodeUpdate(double *U, double *V, double H,
                      size_t NI, size_t NJ, size_t NK,
                      size_t SRCI, size_t SRCJ, size_t SRCK,
                      int REVI, int REVJ, int REVK,
                      size_t ir, size_t jr, size_t kr, double EPS){

    size_t i,j,k;
    size_t mem_cur, mem_il, mem_ir, mem_jl, mem_jr, mem_kl, mem_kr;
    double uold, uxmin, uymin, uzmin, a1, a2, a3, t, t1, t2, t3, unew;
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

    mem_cur = OPENST_MEMADR_3D(i,j,k,NI,NJ,NK);
    mem_il = OPENST_MEMADR_3D(i - 1,j,k,NI,NJ,NK);
    mem_ir = OPENST_MEMADR_3D(i + 1,j,k,NI,NJ,NK);
    mem_jl = OPENST_MEMADR_3D(i,j - 1,k,NI,NJ,NK);
    mem_jr = OPENST_MEMADR_3D(i,j + 1,k,NI,NJ,NK);
    mem_kl = OPENST_MEMADR_3D(i,j,k - 1,NI,NJ,NK);
    mem_kr = OPENST_MEMADR_3D(i,j,k + 1,NI,NJ,NK);

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

    a1 = uxmin;
    a2 = uymin;
    a3 = uzmin;
    
    if(a1 > a3){
      t = a3;
      a1 = a3;
      a3 = t;
    }
    if(a1 > a2){
      t = a1;
      a1 = a2;
      a2 = t;
    }
    if(a2 > a3){
      t = a3;
      a3 = a2;
      a2 = t;
    }

    t1 = a1 - a2;
    t2 = H / V[mem_cur];
    unew = a1 + t2;
    if (unew > a2){
        unew = (a1 + a2 + sqrt(2.0 * t2 * t2 - t1 * t1))/2.0;
        if (unew > a3){
            t3 = a1 + a2 + a3;
            unew = (t3 + sqrt(t3 * t3  - 3.0 *
                              (a1 * a1 + a2 * a2 + a3 * a3
                               - t2 * t2)))/3.0;
        }
    }

    if(unew < uold){

        if(fabs(uold - unew) > EPS){
            notconverged = 1;
        }

        U[mem_cur] = unew;
    }

    return notconverged;

}


int OpenST_FSM3D_BlockSerial(double *U, double *V, double H,
                       size_t NI, size_t NJ, size_t NK,
                       size_t SRCI, size_t SRCJ, size_t SRCK,
                       int REVI, int REVJ, int REVK,
                       size_t istart, size_t jstart, size_t kstart,
                       size_t isize, size_t jsize, size_t ksize, double EPS){

    int notconverged;
    size_t ir, jr, kr;

    notconverged = 0;
    for(ir = istart; ir < istart + isize && ir < NI; ++ir){
        for(jr = jstart; jr < jstart + jsize && jr < NJ; ++jr){
            for(kr = kstart; kr < kstart + ksize && kr < NK; ++kr){
                if(OpenST_FSM3D_NodeUpdate(U, V, H, NI, NJ, NK,
                                     SRCI, SRCJ, SRCK,
                                     REVI, REVJ, REVK,
                                     ir, jr, kr, EPS)){
                    notconverged = 1;
                }
            }
        }
    }

    return notconverged;
}
