#include "openst/common/grad.h"


void OpenST_GRAD_Grad_Kernel(OPENST_FLOAT *left, OPENST_FLOAT *center, OPENST_FLOAT *right,
                             OPENST_FLOAT H, OPENST_FLOAT *grad){
    OPENST_FLOAT f1, f2;
    OPENST_FLOAT divh;

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
        divh *= OPENST_FLOAT_2_0;
    }
    *grad = (f2 - f1)/divh;

}


void OpenST_GRAD_Grad3D(OPENST_FLOAT *A, size_t NI, size_t NJ, size_t NK,
            OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
            size_t i, size_t j, size_t k,
            OPENST_FLOAT *gradi, OPENST_FLOAT *gradj, OPENST_FLOAT *gradk){

    OPENST_FLOAT f1, f2;
    OPENST_FLOAT divh;

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
        divh *= OPENST_FLOAT_2_0;
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
        divh *= OPENST_FLOAT_2_0;
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
        divh *= OPENST_FLOAT_2_0;
    }
    *gradk = (f2 - f1)/divh;
}
