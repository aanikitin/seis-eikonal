#include "openst/eikonal/lsm.h"


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


void OpenST_LSM3D_Init(double *U, char *LSM_UNLOCKED, size_t NI, size_t NJ, size_t NK,
                       size_t SRCI, size_t SRCJ, size_t SRCK){
    size_t i, j, k;
    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = DBL_MAX;
                LSM_UNLOCKED[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = 0;
            }
        }
    }
    U[OPENST_MEMADR_3D(SRCI,SRCJ,SRCK,NI,NJ,NK)] = 0.0;
    if(SRCI > 0){
        LSM_UNLOCKED[OPENST_MEMADR_3D(SRCI - 1,SRCJ,SRCK,NI,NJ,NK)] = 1;
    }
    if(SRCI < NI - 1){
        LSM_UNLOCKED[OPENST_MEMADR_3D(SRCI + 1,SRCJ,SRCK,NI,NJ,NK)] = 1;
    }
    if(SRCJ > 0){
        LSM_UNLOCKED[OPENST_MEMADR_3D(SRCI,SRCJ - 1,SRCK,NI,NJ,NK)] = 1;
    }
    if(SRCJ < NJ - 1){
        LSM_UNLOCKED[OPENST_MEMADR_3D(SRCI,SRCJ + 1,SRCK,NI,NJ,NK)] = 1;
    }
    if(SRCK > 0){
        LSM_UNLOCKED[OPENST_MEMADR_3D(SRCI,SRCJ,SRCK - 1,NI,NJ,NK)] = 1;
    }
    if(SRCK < NK - 1){
        LSM_UNLOCKED[OPENST_MEMADR_3D(SRCI,SRCJ,SRCK + 1,NI,NJ,NK)] = 1;
    }
}


int OpenST_LSM3D_NodeUpdate(double *U, char *LSM_UNLOCKED, double *V, double H,
                            size_t NI, size_t NJ, size_t NK,
                            size_t SRCI, size_t SRCJ, size_t SRCK,
                            int REVI, int REVJ, int REVK,
                            size_t ir, size_t jr, size_t kr, double EPS){

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
    
    mem_cur = OPENST_MEMADR_3D(i,j,k,NI,NJ,NK);
    
    if(LSM_UNLOCKED[mem_cur]){

        if(i == SRCI && j == SRCJ && k == SRCK){
            return notconverged;
        }

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

        /* in LSM tests this is actually slightly faster than sorting network on swaps */
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

            if(i > 0){
                if(!LSM_UNLOCKED[mem_il]){
                    if(U[mem_il] > U[mem_cur]){
                        LSM_UNLOCKED[mem_il] = 1;
                    }
                }
            }

            if(i < NI - 1){
                if(!LSM_UNLOCKED[mem_ir]){
                    if(U[mem_ir] > U[mem_cur]){
                        LSM_UNLOCKED[mem_ir] = 1;
                    }
                }
            }

            if(j > 0){
                if(!LSM_UNLOCKED[mem_jl]){
                    if(U[mem_jl] > U[mem_cur]){
                        LSM_UNLOCKED[mem_jl] = 1;
                    }
                }
            }

            if(j < NJ - 1){
                if(!LSM_UNLOCKED[mem_jr]){
                    if(U[mem_jr] > U[mem_cur]){
                        LSM_UNLOCKED[mem_jr] = 1;
                    }
                }
            }

            if(k > 0){
                if(!LSM_UNLOCKED[mem_kl]){
                    if(U[mem_kl] > U[mem_cur]){
                        LSM_UNLOCKED[mem_kl] = 1;
                    }
                }
            }

            if(k < NK - 1){
                if(!LSM_UNLOCKED[mem_kr]){
                    if(U[mem_kr] > U[mem_cur]){
                        LSM_UNLOCKED[mem_kr] = 1;
                    }
                }
            }

        }

        LSM_UNLOCKED[mem_cur] = 0;

    }

    return notconverged;
}


int OpenST_LSM3D_BlockSerial(double *U, char *LSM_UNLOCKED, double *V, double H,
                             size_t NI, size_t NJ, size_t NK,
                             size_t SRCI, size_t SRCJ, size_t SRCK,
                             int REVI, int REVJ, int REVK,
                             size_t istart, size_t jstart, size_t kstart,
                             size_t isize, size_t jsize, size_t ksize,
                             double EPS){

    int notconverged;
    size_t i, j, ir, jr, kr, iend, jend, kend;

    iend = MIN(istart + isize, NI);
    jend = MIN(jstart + jsize, NJ);
    kend = MIN(kstart + ksize, NK);

    notconverged = 0;
    for(ir = istart; ir < iend; ir += 2u){
        for(jr = jstart; jr < jend; jr += 2u){
            for(i = ir; (i < (ir + 2u)) && (i < iend); ++i){
                for(j = jr; (j < (jr + 2u)) && (j < jend); ++j){
                    for(kr = kstart; kr < kend; ++kr){
                        if(OpenST_LSM3D_NodeUpdate(U, LSM_UNLOCKED, V, H,
                                                   NI, NJ, NK,
                                                   SRCI, SRCJ, SRCK,
                                                   REVI, REVJ, REVK,
                                                   i, j, kr, EPS)){
                            notconverged = 1;
                        }
                    }
                }
            }
        }
    }

    return notconverged;
}
