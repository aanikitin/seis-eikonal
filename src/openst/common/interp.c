#include "openst/common/interp.h"


void OpenST_INTERP_Trilinear(double *A, size_t NI, size_t NJ, size_t NK,
                             double HI, double HJ, double HK,
                             double PI, double PJ, double PK,
                             double *VAL){

    size_t ii[2], ji[2], ki[2];
    double ic[2], jc[2], kc[2];
    int interp_i, interp_j, interp_k, interp_dims;

    OpenST_INTERP_Trilinear_Neighboors(
        HI, HJ, HK,
        PI, PJ, PK,
        ii, ji, ki,
        ic, jc, kc,
        &interp_i, &interp_j, &interp_k, &interp_dims);

    OpenST_INTERP_Trilinear_Compute(
        A,
        NI, NJ, NK,
        PI, PJ, PK,
        ii, ji, ki,
        ic, jc, kc,
        interp_i, interp_j, interp_k, interp_dims, VAL);

}


void OpenST_INTERP_Linear_Formula(double *f,
                                 double i,
                                 double f0, double f1,
                                 double i0, double i1){
    double div;

    div = i1 - i0;

    *f = f0 + (f1 - f0) * (i - i0) / (div);
}


void OpenST_INTERP_Bilinear_Formula(double *f,
                                   double i, double j,
                                   double f00, double f01,
                                   double f10, double f11,
                                   double i0, double j0,
                                   double i1, double j1){
    double divi, divj, m0, m1, fij0, fij1;

    divi = i1 - i0;

    m0 = (i1 - i) / divi;
    m1 = (i - i0) / divi;

    fij0 = m0 * f00 + m1 * f10;
    fij1 = m0 * f01 + m1 * f11;

    divj = j1 - j0;

    m0 = (j1 - j) / divj;
    m1 = (j - j0) / divj;

    *f = m0 * fij0 + m1 * fij1;
}

/* checked */
void OpenST_INTERP_Trilinear_Formula(double *f,
                                    double i, double j, double k,
                                    double f000, double f001,
                                    double f010, double f011,
                                    double f100, double f101,
                                    double f110, double f111,
                                    double i0, double j0, double k0,
                                    double i1, double j1, double k1){
    double divi, divj, divk, mi, mj, mk;
    double f00, f01, f10, f11, f0, f1;

    divi = i1 - i0;
    mi = (i - i0) / divi;

    f00 = f000 * (1.0 - mi) + f100 * mi;
    f01 = f001 * (1.0 - mi) + f101 * mi;
    f10 = f010 * (1.0 - mi) + f110 * mi;
    f11 = f011 * (1.0 - mi) + f111 * mi;

    divj = j1 - j0;
    mj = (j - j0) / divj;

    f0 = f00 * (1 - mj) + f10 * mj;
    f1 = f01 * (1 - mj) + f11 * mj;

    divk = k1 - k0;
    mk = (k - k0) / divk;

    *f = f0 * (1.0 - mk) + f1 * mk;
}


OPENST_API void OpenST_INTERP_Trilinear_Neighboors(
    double HI, double HJ, double HK,
    double PI, double PJ, double PK,
    size_t ii[static 2], size_t ji[static 2], size_t ki[static 2],
    double ic[static 2], double jc[static 2], double kc[static 2],
    int *interp_i, int *interp_j, int *interp_k, int *interp_dims){

    ii[0] = (size_t) floor(PI / HI);
    ic[0] = (double) ii[0] * HI;

    ii[1] = (size_t) ceil(PI / HI);
    ic[1] = (double) ii[1] * HI;

    ji[0] = (size_t) floor(PJ / HJ);
    jc[0] = (double) ji[0] * HJ;

    ji[1] = (size_t) ceil(PJ / HJ);
    jc[1] = (double) ji[1] * HJ;

    ki[0] = (size_t) floor(PK / HK);
    kc[0] = (double) ki[0] * HK;

    ki[1] = (size_t) ceil(PK / HK);
    kc[1] = (double) ki[1] * HK;

    *interp_dims = 0;
    if((*interp_i = (ii[0] != ii[1]))){
        ++(*interp_dims);
    }
    if((*interp_j = (ji[0] != ji[1]))){
        ++(*interp_dims);
    }
    if((*interp_k = (ki[0] != ki[1]))){
        ++(*interp_dims);
    }

}


OPENST_API void OpenST_INTERP_Trilinear_Compute(
    double *A,
    size_t NI, size_t NJ, size_t NK,
    double PI, double PJ, double PK,
    size_t ii[static 2], size_t ji[static 2], size_t ki[static 2],
    double ic[static 2], double jc[static 2], double kc[static 2],
    int interp_i, int interp_j, int interp_k, int interp_dims,
    double *VAL){

    if(interp_dims == 0){

        *VAL = A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)];

    } else {

        if (interp_dims == 3) {

            OpenST_INTERP_Trilinear_Formula(VAL, PI, PJ, PK,
                    A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[0], ji[1], ki[1], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[1], ji[0], ki[1], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[1], ji[1], ki[0], NI, NJ, NK)],
                    A[OPENST_MEMADR_3D(ii[1], ji[1], ki[1], NI, NJ, NK)],
                    ic[0], jc[0], kc[0],
                    ic[1], jc[1], kc[1]);

        } else if (interp_dims == 2) {

            if (!interp_i) {

                OpenST_INTERP_Bilinear_Formula(VAL, PJ, PK,
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[1], ki[1], NI, NJ, NK)],
                        jc[0], kc[0],
                        jc[1], kc[1]);

            } else if (!interp_j) {

                OpenST_INTERP_Bilinear_Formula(VAL, PI, PK,
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[1], ji[0], ki[1], NI, NJ, NK)],
                        ic[0], kc[0],
                        ic[1], kc[1]);

            } else if (!interp_k) {

                OpenST_INTERP_Bilinear_Formula(VAL, PI, PJ,
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[1], ji[1], ki[0], NI, NJ, NK)],
                        ic[0], jc[0],
                        ic[1], jc[1]);

            }

        } else if (interp_dims == 1) {

            if (interp_i) {

                OpenST_INTERP_Linear_Formula(VAL, PI,
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                        ic[0], ic[1]);

            } else if (interp_j) {

                OpenST_INTERP_Linear_Formula(VAL, PJ,
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                        jc[0], jc[1]);

            } else if (interp_k) {

                OpenST_INTERP_Linear_Formula(VAL, PK,
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        A[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                        kc[0], kc[1]);

            }

        }

    }

}
