#include "openst/common/interp.h"


OPENST_ERR OpenST_INTERP_3D(OPENST_FLOAT *A, size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
                             OPENST_FLOAT *VAL){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;

    if ( OpenST_CRS_IsPointNotWithinBounds(PI, PJ, PK,
                                        NI, NJ, NK,
                                        HI, HJ, HK) ) {
        errcode = OPENST_ERR_PARAM_INVALID;
        goto EXIT;
    }

    OpenST_INTERP_Trilinear(A, NI, NJ, NK,
                            HI, HJ, HK,
                            PI, PJ, PK,
                            VAL);

EXIT:
    return errcode;
}


void OpenST_INTERP_Trilinear(OPENST_FLOAT *A, size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
                             OPENST_FLOAT *VAL){

    size_t ii[2], ji[2], ki[2];
    OPENST_FLOAT ic[2], jc[2], kc[2];
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


void OpenST_INTERP_Linear_Formula(OPENST_FLOAT *f,
                                 OPENST_FLOAT i,
                                 OPENST_FLOAT f0, OPENST_FLOAT f1,
                                 OPENST_FLOAT i0, OPENST_FLOAT i1){
    OPENST_FLOAT div;

    div = i1 - i0;

    *f = f0 + (f1 - f0) * (i - i0) / (div);
}


void OpenST_INTERP_Bilinear_Formula(OPENST_FLOAT *f,
                                   OPENST_FLOAT i, OPENST_FLOAT j,
                                   OPENST_FLOAT f00, OPENST_FLOAT f01,
                                   OPENST_FLOAT f10, OPENST_FLOAT f11,
                                   OPENST_FLOAT i0, OPENST_FLOAT j0,
                                   OPENST_FLOAT i1, OPENST_FLOAT j1){
    OPENST_FLOAT divi, divj, m0, m1, fij0, fij1;

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


void OpenST_INTERP_Trilinear_Formula(OPENST_FLOAT *f,
                                    OPENST_FLOAT i, OPENST_FLOAT j, OPENST_FLOAT k,
                                    OPENST_FLOAT f000, OPENST_FLOAT f001,
                                    OPENST_FLOAT f010, OPENST_FLOAT f011,
                                    OPENST_FLOAT f100, OPENST_FLOAT f101,
                                    OPENST_FLOAT f110, OPENST_FLOAT f111,
                                    OPENST_FLOAT i0, OPENST_FLOAT j0, OPENST_FLOAT k0,
                                    OPENST_FLOAT i1, OPENST_FLOAT j1, OPENST_FLOAT k1){
    OPENST_FLOAT divi, divj, divk, mi, mj, mk;
    OPENST_FLOAT f00, f01, f10, f11, f0, f1;

    divi = i1 - i0;
    mi = (i - i0) / divi;

    f00 = f000 * (OPENST_FLOAT_1_0 - mi) + f100 * mi;
    f01 = f001 * (OPENST_FLOAT_1_0 - mi) + f101 * mi;
    f10 = f010 * (OPENST_FLOAT_1_0 - mi) + f110 * mi;
    f11 = f011 * (OPENST_FLOAT_1_0 - mi) + f111 * mi;

    divj = j1 - j0;
    mj = (j - j0) / divj;

    f0 = f00 * (1 - mj) + f10 * mj;
    f1 = f01 * (1 - mj) + f11 * mj;

    divk = k1 - k0;
    mk = (k - k0) / divk;

    *f = f0 * (OPENST_FLOAT_1_0 - mk) + f1 * mk;
}


OPENST_API void OpenST_INTERP_Trilinear_Neighboors(
    OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
    OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
    size_t ii[2], size_t ji[2], size_t ki[2],
    OPENST_FLOAT ic[2], OPENST_FLOAT jc[2], OPENST_FLOAT kc[2],
    int *interp_i, int *interp_j, int *interp_k, int *interp_dims){

    *interp_i = OpenST_FLOAT_GetNeighboorSizeT(PI / HI, &ii[0], &ii[1]);
    *interp_j = OpenST_FLOAT_GetNeighboorSizeT(PJ / HJ, &ji[0], &ji[1]);
    *interp_k = OpenST_FLOAT_GetNeighboorSizeT(PK / HK, &ki[0], &ki[1]);

    *interp_dims = 0;

    if (*interp_i) {
        ++(*interp_dims);
        ic[0] = (OPENST_FLOAT) ii[0] * HI;
        ic[1] = (OPENST_FLOAT) ii[1] * HI;
    } else {
        ic[0] = PI;
        ic[1] = PI;
    }

    if (*interp_j) {
        ++(*interp_dims);
        jc[0] = (OPENST_FLOAT) ji[0] * HJ;
        jc[1] = (OPENST_FLOAT) ji[1] * HJ;
    } else {
        jc[0] = PJ;
        jc[1] = PJ;
    }

    if (*interp_k) {
        ++(*interp_dims);
        kc[0] = (OPENST_FLOAT) ki[0] * HK;
        kc[1] = (OPENST_FLOAT) ki[1] * HK;
    } else {
        kc[0] = PK;
        kc[1] = PK;
    }

}


OPENST_API void OpenST_INTERP_Trilinear_Compute(
    OPENST_FLOAT *A,
    size_t NI, size_t NJ, size_t NK,
    OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
    size_t ii[2], size_t ji[2], size_t ki[2],
    OPENST_FLOAT ic[2], OPENST_FLOAT jc[2], OPENST_FLOAT kc[2],
    int interp_i, int interp_j, int interp_k, int interp_dims,
    OPENST_FLOAT *VAL){

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
