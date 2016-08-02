#include "openst/common/grad.h"


void OpenST_GRAD_Grad_Kernel(double *left, double *center, double *right,
                             double H, double *grad){
    double f1, f2;
    double divh;

    divh = H;
    if(left == NULL){
        f1 = *center;
        f2 = *right;
    } else if (right == NULL) {
        f1 = *left;
        f2 = *center;
    } else {
        f1 = *left;
        f2 = *right;
        divh *= 2.0;
    }
    *grad = (f2 - f1)/divh;

}


void OpenST_GRAD_Grad3D(double *A, size_t NI, size_t NJ, size_t NK,
            double HI, double HJ, double HK,
            size_t i, size_t j, size_t k,
            double *gradi, double *gradj, double *gradk){

    double f1, f2;
    double divh;

    divh = HI;
    if(i == 0){
        f1 = A[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i+1,j,k,NI,NJ,NK)];
    } else if (i == NI - 1) {
        f1 = A[OPENST_MEMADR_3D(i-1,j,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
    } else {
        f1 = A[OPENST_MEMADR_3D(i-1,j,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i+1,j,k,NI,NJ,NK)];
        divh *= 2.0;
    }
    *gradi = (f2 - f1)/divh;

    divh = HJ;
    if(j == 0){
        f1 = A[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j+1,k,NI,NJ,NK)];
    } else if (j == NJ - 1) {
        f1 = A[OPENST_MEMADR_3D(i,j-1,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
    } else {
        f1 = A[OPENST_MEMADR_3D(i,j-1,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j+1,k,NI,NJ,NK)];
        divh *= 2.0;
    }
    *gradj = (f2 - f1)/divh;

    divh = HK;
    if(k == 0){
        f1 = A[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j,k+1,NI,NJ,NK)];
    } else if (k == NK - 1) {
        f1 = A[OPENST_MEMADR_3D(i,j,k-1,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)];
    } else {
        f1 = A[OPENST_MEMADR_3D(i,j,k-1,NI,NJ,NK)];
        f2 = A[OPENST_MEMADR_3D(i,j,k+1,NI,NJ,NK)];
        divh *= 2.0;
    }
    *gradk = (f2 - f1)/divh;
}
