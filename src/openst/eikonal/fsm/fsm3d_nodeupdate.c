/* TODO: this is significantly slower than 1H version, try to optimize?*/
#include "openst/eikonal/fsm.h"


int OpenST_FSM3D_NodeUpdate(double *U, double *V,
                            size_t NI, size_t NJ, size_t NK,
                            double HI, double HJ, double HK,
                            int REVI, int REVJ, int REVK,
                            size_t ir, size_t jr, size_t kr,
                            double EPS){

    size_t i,j,k;
    size_t mem_cur, mem_il, mem_ir, mem_jl, mem_jr, mem_kl, mem_kr;
    double uold, v, uimin, ujmin, ukmin;
    double a1, a2, a3;
    double h1, h2, h3;
    double m1, m2, m3, mf;
    double a, bp, b, c, D;
    double unew;
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
        ujmin = fmin(U[mem_jl], U[mem_jr]);
    }

    if(i == 0){
        uimin = U[mem_ir];
    } else if (i == NI - 1){
        uimin = U[mem_il];
    } else {
        uimin = fmin(U[mem_il], U[mem_ir]);
    }

    if(k == 0){
        ukmin = U[mem_kr];
    } else if (k == NK - 1){
        ukmin = U[mem_kl];
    } else {
        ukmin = fmin(U[mem_kl], U[mem_kr]);
    }

    /* rolled back sorting algorithm because of loss of accuracy when sorting
     * using sorting by swapping */
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

        b = -2.0 * bp;
        D = b * b - 4.0 * a * c;

        unew = (-b + sqrt(D)) / (2.0 * a);

        if (unew > a3){
            m3 = h1 * h1 * h2 * h2;
            a += m3;
            bp += m3 * a3;
            c += m3 * a3 * a3;

            b = -2.0 * bp;
            D = b * b - 4.0 * a * c;

            unew = (-b + sqrt(D)) / (2.0 * a);
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
