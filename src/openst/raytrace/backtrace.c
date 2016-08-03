//TODO: recheck and cleanup
//TODO: improve accuracy
#include "openst/raytrace/backtrace.h"

#define TSTEP_DEFAULT_MULT 0.999
#define DEBUG_LOG 0

#include <float.h>

#if DEBUG_LOG
#include <stdio.h>
#endif


double OpenST_BRT3D_SuggestTSTEP(double vmax, double HI, double HJ, double HK) {
    double tstep = fmin(HI, (fmin(HJ, HK))) / vmax * TSTEP_DEFAULT_MULT;
    return tstep;
}


OPENST_ERR OpenST_BRT3D_Step(double *T, double *V,
    size_t NI, size_t NJ, size_t NK,
    double HI, double HJ, double HK,
    double TSTEP,
    double CURI, double CURJ, double CURK,
    double *DSTI, double *DSTJ, double *DSTK) {

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    double TC, TIL, TIR, TJL, TJR, TKL, TKR;
    double *TILa, *TIRa, *TJLa, *TJRa, *TKLa, *TKRa;
    double gradi, gradj, gradk, grad_length;
    double vel, vali, valj, valk;

    /*
     * TODO: check can fail on MAX domain edge on certain grid sizes and steps
     * store DOMAIN_MAX and compare with it in the future release,
     * do not recalculate
     */
    if((CURI < 0) || (CURJ < 0) || (CURK < 0) ||
            (CURI > ((double)(NI - 1) * HI)) ||
			(CURJ > ((double)(NJ - 1) * HJ)) ||
            (CURK > ((double)(NK - 1) * HK)) ){
        errcode = OPENST_ERR_ALG_DOMAIN;
        goto EXIT;
    }

    OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK,&TC);

    TILa = &TIL;
    if(CURI >= HI/2.0){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI - HI/2.0,CURJ,CURK,TILa);
    } else {
        TILa = NULL;
    }

    TIRa = &TIR;
    if( (CURI + HI/2.0) <= ((double)(NI - 1) * HI) ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI + HI/2.0,CURJ,CURK,TIRa);
    } else {
        TIRa = NULL;
    }

    TJLa = &TJL;
    if(CURJ >= HJ/2.0){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ - HJ/2.0,CURK,TJLa);
    } else {
        TJLa = NULL;
    }

    TJRa = &TJR;
    if( (CURJ + HJ / 2.0) <= ((double)(NJ - 1) * HJ) ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ + HJ/2.0,CURK,TJRa);
    } else {
        TJRa = NULL;
    }

    TKLa = &TKL;
    if(CURK >= HK/2.0){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK - HK/2.0,TKLa);
    } else {
        TKLa = NULL;
    }

    TKRa = &TKR;
    if( (CURK + HK / 2.0) <= ((double)(NK - 1) * HK) ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK + HK/2.0,TKRa);
    } else {
        TKRa = NULL;
    }

    OpenST_GRAD_Grad_Kernel(TILa,&TC,TIRa,HI,&gradi);
    OpenST_GRAD_Grad_Kernel(TJLa,&TC,TJRa,HJ,&gradj);
    OpenST_GRAD_Grad_Kernel(TKLa,&TC,TKRa,HK,&gradk);

    grad_length = sqrt(pow(gradi, 2.0) + pow(gradj, 2.0) + pow(gradk, 2.0));

    OpenST_INTERP_Trilinear(V,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK,&vel);

#if DEBUG_LOG
    printf("GRAD: [%f %f %f]/%f V %f\n", gradi, gradj, gradk, grad_length, vel);
#endif

    if (grad_length > 0) {
        vali = CURI - gradi / grad_length * TSTEP * vel;
        valj = CURJ - gradj / grad_length * TSTEP * vel;
        valk = CURK - gradk / grad_length * TSTEP * vel;
    }
    else {
        errcode = OPENST_ERR_DIV_BY_ZERO;
        goto EXIT;
    }

    *DSTI = vali;
    *DSTJ = valj;
    *DSTK = valk;

EXIT:
    return errcode;
}


OPENST_ERR OpenST_BRT3D_Trace(double *T, double *V,
    size_t NI, size_t NJ, size_t NK,
    double HI, double HJ, double HK, double TSTEP,
    double RCVI, double RCVJ, double RCVK,
    double SRCI, double SRCJ, double SRCK,
    size_t MAX_SEG,
    double **RAY, size_t *RAY_NI, size_t *RAY_NJ) {

    OPENST_ERR errcode = OPENST_ERR_ALG_EARLY_TERM;
    double CUR[3];
    double DST[3];
    double SRCI_L, SRCI_H, SRCJ_L, SRCJ_H, SRCK_L, SRCK_H;
    size_t numseg;

    struct OpenST_DYNARR arr;
    struct OpenST_DYNARR *arrptr = &arr;

    CUR[0] = RCVI;
    CUR[1] = RCVJ;
    CUR[2] = RCVK;

    SRCI_L = SRCI - HI / 2.0;
    SRCI_H = SRCI + HI / 2.0;
    SRCJ_L = SRCJ - HJ / 2.0;
    SRCJ_H = SRCJ + HJ / 2.0;
    SRCK_L = SRCK - HK / 2.0;
    SRCK_H = SRCK + HK / 2.0;

    if (OpenST_DYNARR_Init(arrptr, NI + NJ + NK, sizeof(double) * 3) == NULL) {
        errcode = OPENST_ERR_MEMORY_ALLOC;
        goto EXIT;
    }

    numseg = 0;
    while (numseg < MAX_SEG) {

#if DEBUG_LOG
        printf("[%e;%e;%e]\n", CUR[0], CUR[1], CUR[2]);
        fflush(stdout);
#endif

        if (OpenST_DYNARR_Pushback(arrptr, CUR) == NULL) {
            errcode = OPENST_ERR_MEMORY_ALLOC;
            goto EXIT;
        }

        ++numseg;

        if ((CUR[0] > SRCI_L) && (CUR[0] < SRCI_H) &&
            (CUR[1] > SRCJ_L) && (CUR[1] < SRCJ_H) &&
            (CUR[2] > SRCK_L) && (CUR[2] < SRCK_H)) {
            errcode = OPENST_ERR_SUCCESS;
#if DEBUG_LOG
            printf("SOURCE REACHED\n");
            fflush(stdout);
#endif
            break;
        }

        if ((errcode = OpenST_BRT3D_Step(T, V, NI, NJ, NK, HI, HJ, HK, TSTEP,
            CUR[0], CUR[1], CUR[2],
            &DST[0], &DST[1], &DST[2]))) {
            goto EXIT;
        }

        CUR[0] = DST[0];
        CUR[1] = DST[1];
        CUR[2] = DST[2];
    }

    if (OpenST_DYNARR_Shrink(arrptr) == NULL) {
        errcode = OPENST_ERR_MEMORY_ALLOC;
        goto EXIT;
    }

    *RAY = arr.data;
    *RAY_NJ = 3;
    *RAY_NI = arr.num;

EXIT:
    return errcode;
}
