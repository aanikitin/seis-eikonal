#include "openst/common/grad.h"


void OpenST_INTERP_Linear(double *f, double i,
                                    double f0, double f1,
                                    double i0, double i1){
    double div;

    div = i1 - i0;

    *f = f0 + (f1 - f0) * (i - i0) / (div);
}


void OpenST_INTERP_Bilinear(double *f, double i, double j,
                                      double f00, double f01,
                                      double f10, double f11,
                                      double i0,
                                      double j0, double i1, double j1){
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


void OpenST_INTERP_Trilinear(double *f,
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

    f00 = f000 * (1 - mi) + f100 * mi;
    f01 = f001 * (1 - mi) + f101 * mi;
    f10 = f010 * (1 - mi) + f110 * mi;
    f11 = f011 * (1 - mi) + f111 * mi;

    divj = j1 - j0;
    mj = (j - j0) / divj;

    f0 = f00 * (1 - mj) + f10 * mj;
    f1 = f01 * (1 - mj) + f11 * mj;

    divk = k1 - k0;
    mk = (k - k0) / divk;

    *f = f0 * (1 - mk) + f1 * mk;
}
