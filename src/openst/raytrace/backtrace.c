//TODO: integrate by ds
//TODO: improve accuracy (Runge-Kutta 4th order, improve eikonal accuracy)
//TODO: improve handling of near source termination
#include "openst/raytrace/backtrace.h"

#define TSTEP_DEFAULT_MULT OPENST_FLOAT_SFX(0.999)
#define DEBUG_LOG 0

#if DEBUG_LOG
#include <stdio.h>
#endif


OPENST_FLOAT OpenST_BRT3D_SuggestTSTEP(OPENST_FLOAT vmax, OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK) {
    OPENST_FLOAT tstep = OPENST_FLOAT_FMIN(HI, (OPENST_FLOAT_FMIN(HJ, HK))) / vmax * TSTEP_DEFAULT_MULT;
    return tstep;
}


OPENST_ERR OpenST_BRT3D_Step(OPENST_FLOAT *T, OPENST_FLOAT *V,
    size_t NI, size_t NJ, size_t NK,
    OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
    OPENST_FLOAT TSTEP,
    OPENST_FLOAT CURI, OPENST_FLOAT CURJ, OPENST_FLOAT CURK,
    OPENST_FLOAT *DSTI, OPENST_FLOAT *DSTJ, OPENST_FLOAT *DSTK) {

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    OPENST_FLOAT TC, TIL, TIR, TJL, TJR, TKL, TKR;
    OPENST_FLOAT *TILa, *TIRa, *TJLa, *TJRa, *TKLa, *TKRa;
    OPENST_FLOAT gradi, gradj, gradk, grad_length;
    OPENST_FLOAT vel, vali, valj, valk;

    if( OpenST_CRS_IsPointNotWithinBounds(CURI, CURJ, CURK,
                                       NI, NJ, NK,
                                       HI, HJ, HK) ){
        errcode = OPENST_ERR_ALG_OUT_OF_BOUNDS;
        goto EXIT;
    }

    OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK,&TC);

    TILa = &TIL;
    if(CURI >= HI){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI - HI,CURJ,CURK,TILa);
    } else {
        TILa = NULL;
    }

    TIRa = &TIR;
    if( (CURI + HI) <= ((OPENST_FLOAT)(NI - 1) * HI) ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI + HI,CURJ,CURK,TIRa);
    } else {
        TIRa = NULL;
    }

    TJLa = &TJL;
    if(CURJ >= HJ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ - HJ,CURK,TJLa);
    } else {
        TJLa = NULL;
    }

    TJRa = &TJR;
    if( (CURJ + HJ) <= ((OPENST_FLOAT)(NJ - 1) * HJ) ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ + HJ,CURK,TJRa);
    } else {
        TJRa = NULL;
    }

    TKLa = &TKL;
    if(CURK >= HK){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK - HK,TKLa);
    } else {
        TKLa = NULL;
    }

    TKRa = &TKR;
    if( (CURK + HK) <= ((OPENST_FLOAT)(NK - 1) * HK) ){
        OpenST_INTERP_Trilinear(T,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK + HK,TKRa);
    } else {
        TKRa = NULL;
    }

    OpenST_GRAD_Grad_Kernel(TILa,&TC,TIRa,HI,&gradi);
    OpenST_GRAD_Grad_Kernel(TJLa,&TC,TJRa,HJ,&gradj);
    OpenST_GRAD_Grad_Kernel(TKLa,&TC,TKRa,HK,&gradk);

    grad_length = OPENST_FLOAT_SQRT(OPENST_FLOAT_POW(gradi, OPENST_FLOAT_2_0) + OPENST_FLOAT_POW(gradj, OPENST_FLOAT_2_0) + OPENST_FLOAT_POW(gradk, OPENST_FLOAT_2_0));

    OpenST_INTERP_Trilinear(V,NI,NJ,NK,HI,HJ,HK,CURI,CURJ,CURK,&vel);

#if DEBUG_LOG
    printf("GRAD: [%f %f %f]/%f V %f\n", gradi, gradj, gradk, grad_length, vel);
#endif

    if (grad_length > 0) {
        vali = CURI - gradi / grad_length * vel * TSTEP;
        valj = CURJ - gradj / grad_length * vel * TSTEP;
        valk = CURK - gradk / grad_length * vel * TSTEP;
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


OPENST_ERR OpenST_BRT3D_Trace(OPENST_FLOAT *T, OPENST_FLOAT *V,
    size_t NI, size_t NJ, size_t NK,
    OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK, OPENST_FLOAT TSTEP,
    OPENST_FLOAT RCVI, OPENST_FLOAT RCVJ, OPENST_FLOAT RCVK,
    OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
    size_t MAX_SEG,
    OPENST_FLOAT **RAY, size_t *RAY_NI, size_t *RAY_NJ) {

    OPENST_ERR errcode = OPENST_ERR_ALG_TERMINATED;
    OPENST_FLOAT CUR[3];
    OPENST_FLOAT DST[3];
    OPENST_FLOAT SRCI_L, SRCI_H, SRCJ_L, SRCJ_H, SRCK_L, SRCK_H;
    size_t numseg;

    struct OpenST_DYNARR arr;
    struct OpenST_DYNARR *arrptr = &arr;

    CUR[0] = RCVI;
    CUR[1] = RCVJ;
    CUR[2] = RCVK;

    SRCI_L = SRCI - HI / OPENST_FLOAT_2_0;
    SRCI_H = SRCI + HI / OPENST_FLOAT_2_0;
    SRCJ_L = SRCJ - HJ / OPENST_FLOAT_2_0;
    SRCJ_H = SRCJ + HJ / OPENST_FLOAT_2_0;
    SRCK_L = SRCK - HK / OPENST_FLOAT_2_0;
    SRCK_H = SRCK + HK / OPENST_FLOAT_2_0;

    if (OpenST_DYNARR_Init(arrptr, NI + NJ + NK, sizeof(OPENST_FLOAT) * 3) == NULL) {
        errcode = OPENST_ERR_MEMORY_ALLOC;
        goto EXIT;
    }

    numseg = 0;
    while (numseg < MAX_SEG) {

#if DEBUG_LOG
        printf("[%.20e;%.20e;%.20e]\n", CUR[0], CUR[1], CUR[2]);
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
            goto FINISH;
        }

        CUR[0] = DST[0];
        CUR[1] = DST[1];
        CUR[2] = DST[2];
    }

FINISH:
    if (OpenST_DYNARR_Shrink(arrptr) == NULL) {
        errcode = OPENST_ERR_MEMORY_ALLOC;
        goto EXIT;
    }

    *RAY = (OPENST_FLOAT *) arr.data;
    *RAY_NJ = 3;
    *RAY_NI = arr.num;

EXIT:
    return errcode;
}
