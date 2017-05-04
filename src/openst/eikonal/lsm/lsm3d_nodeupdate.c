#include "openst/eikonal/lsm.h"


int OpenST_LSM3D_NodeUpdate_1H(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                               size_t NI, size_t NJ, size_t NK,
                               OPENST_FLOAT H,
                               int REVI, int REVJ, int REVK,
                               size_t ir, size_t jr, size_t kr, OPENST_FLOAT EPS){

    size_t i,j,k;
    size_t mem_cur, mem_il, mem_ir, mem_jl, mem_jr, mem_kl, mem_kr;
    OPENST_FLOAT uold, uxmin, uymin, uzmin, a1, a2, a3, t1, t2, t3, unew;
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
            uxmin = OPENST_FLOAT_FMIN(U[mem_jl], U[mem_jr]);
        }

        if(i == 0){
            uymin = U[mem_ir];
        } else if (i == NI - 1){
            uymin = U[mem_il];
        } else {
            uymin = OPENST_FLOAT_FMIN(U[mem_il], U[mem_ir]);
        }

        if(k == 0){
            uzmin = U[mem_kr];
        } else if (k == NK - 1){
            uzmin = U[mem_kl];
        } else {
            uzmin = OPENST_FLOAT_FMIN(U[mem_kl], U[mem_kr]);
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
        t2 = H / V[mem_cur];
        unew = a1 + t2;
        if (unew > a2){
            unew = (a1 + a2 + OPENST_FLOAT_SQRT(OPENST_FLOAT_2_0 * t2 * t2 - t1 * t1))/OPENST_FLOAT_2_0;
            if (unew > a3){
                t3 = a1 + a2 + a3;
                unew = (t3 + OPENST_FLOAT_SQRT(t3 * t3  - OPENST_FLOAT_3_0 *
                                  (a1 * a1 + a2 * a2 + a3 * a3
                                   - t2 * t2)))/OPENST_FLOAT_3_0;
            }
        }

        if(unew < uold){

            if(OPENST_FLOAT_FABS(uold - unew) > EPS){
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


int OpenST_LSM3D_NodeUpdate(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                            size_t NI, size_t NJ, size_t NK,
                            OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                            int REVI, int REVJ, int REVK,
                            size_t ir, size_t jr, size_t kr, OPENST_FLOAT EPS){

    size_t i,j,k;
    size_t mem_cur, mem_il, mem_ir, mem_jl, mem_jr, mem_kl, mem_kr;
    OPENST_FLOAT uold, v, uimin, ujmin, ukmin;
    OPENST_FLOAT a1, a2, a3;
    OPENST_FLOAT h1, h2, h3;
    OPENST_FLOAT m1, m2, m3, mf;
    OPENST_FLOAT a, bp, b, c, D;
    OPENST_FLOAT unew;
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

        mem_il = OPENST_MEMADR_3D(i - 1,j,k,NI,NJ,NK);
        mem_ir = OPENST_MEMADR_3D(i + 1,j,k,NI,NJ,NK);
        mem_jl = OPENST_MEMADR_3D(i,j - 1,k,NI,NJ,NK);
        mem_jr = OPENST_MEMADR_3D(i,j + 1,k,NI,NJ,NK);
        mem_kl = OPENST_MEMADR_3D(i,j,k - 1,NI,NJ,NK);
        mem_kr = OPENST_MEMADR_3D(i,j,k + 1,NI,NJ,NK);

        uold = U[mem_cur];
        v = V[mem_cur];

        if(j == 0){
            ujmin = U[mem_jr];
        } else if (j == NJ - 1){
            ujmin = U[mem_jl];
        } else {
            ujmin = OPENST_FLOAT_FMIN(U[mem_jl], U[mem_jr]);
        }

        if(i == 0){
            uimin = U[mem_ir];
        } else if (i == NI - 1){
            uimin = U[mem_il];
        } else {
            uimin = OPENST_FLOAT_FMIN(U[mem_il], U[mem_ir]);
        }

        if(k == 0){
            ukmin = U[mem_kr];
        } else if (k == NK - 1){
            ukmin = U[mem_kl];
        } else {
            ukmin = OPENST_FLOAT_FMIN(U[mem_kl], U[mem_kr]);
        }

        if (ujmin <= uimin && ujmin <= ukmin){
            a1 = ujmin;
            h1 = HJ;
            if (uimin <= ukmin){
                a2 = uimin;
                h2 = HI;
                a3 = ukmin;
                h3 = HK;
            } else {
                a2 = ukmin;
                h2 = HK;
                a3 = uimin;
                h3 = HI;
            }
        } else if (uimin <= ujmin && uimin <= ukmin){
            a1 = uimin;
            h1 = HI;
            if (ujmin <= ukmin){
                a2 = ujmin;
                h2 = HJ;
                a3 = ukmin;
                h3 = HK;
            } else {
                a2 = ukmin;
                h2 = HK;
                a3 = ujmin;
                h3 = HJ;
            }
        } else {
            a1 = ukmin;
            h1 = HK;
            if (ujmin <= uimin){
                a2 = ujmin;
                h2 = HJ;
                a3 = uimin;
                h3 = HI;
            } else {
                a2 = uimin;
                h2 = HI;
                a3 = ujmin;
                h3 = HJ;
            }
        }

        unew = a1 + h1 / v;

        if (unew > a2){
            m1 = h2 * h2 * h3 * h3;
            m2 = h1 * h1 * h3 * h3;
            mf = h1 * h1 * h2 * h2 * h3 * h3;

            a = m1 + m2;
            bp = m1 * a1 + m2 * a2;
            c = m1 * a1 * a1 + m2 * a2 * a2 - mf / (v * v);

            b = -OPENST_FLOAT_2_0 * bp;
            D = b * b - OPENST_FLOAT_4_0 * a * c;

            unew = (-b + OPENST_FLOAT_SQRT(D)) / (OPENST_FLOAT_2_0 * a);

            if (unew > a3){
                m3 = h1 * h1 * h2 * h2;
                a += m3;
                bp += m3 * a3;
                c += m3 * a3 * a3;

                b = -OPENST_FLOAT_2_0 * bp;
                D = b * b - OPENST_FLOAT_4_0 * a * c;

                unew = (-b + OPENST_FLOAT_SQRT(D)) / (OPENST_FLOAT_2_0 * a);
            }
        }

        if(unew < uold){

            if(OPENST_FLOAT_FABS(uold - unew) > EPS){
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
