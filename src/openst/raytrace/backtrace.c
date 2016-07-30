#include "openst/raytrace/backtrace.h"

#define DEBUG_LOG 0

#include <float.h>

#if DEBUG_LOG
#include <stdio.h>
#endif


double OpenST_BRT3D_SuggestTSTEP(double vmax, double HI, double HJ, double HK){
    double tstep = fmin(HI,(fmin(HJ,HK))) / vmax;
    return tstep;
}


OPENST_ERR OpenST_BRT3D_Step(double *T, double *V,
                     size_t NI, size_t NJ, size_t NK,
                     double HI, double HJ, double HK,
                     double TSTEP,
                     double CURI, double CURJ, double CURK,
                     size_t ind_cur_i, size_t ind_cur_j, size_t ind_cur_k,
                     double *DSTI, double *DSTJ, double *DSTK){

    double gradi, gradj, gradk, grad_length;
    double vel, vali, valj, valk;

    OpenST_GRAD_Grad3D(T, NI, NJ, NK, HI, HJ, HK, ind_cur_i, ind_cur_j, ind_cur_k,
           &gradi, &gradj, &gradk);

    grad_length = sqrt(pow(gradi,2.0) + pow(gradj,2.0) + pow(gradk,2.0));

    vel = V[OPENST_MEMADR_3D(ind_cur_i,ind_cur_j,ind_cur_k,NI,NJ,NK)];

#if DEBUG_LOG
    printf("GRAD: [%f %f %f]/%f V %f\n", gradi, gradj, gradk, grad_length, vel);
#endif

    if(grad_length > 0){
        vali = CURI - gradi/grad_length * TSTEP * vel;
        valj = CURJ - gradj/grad_length * TSTEP * vel;
        valk = CURK - gradk/grad_length * TSTEP * vel;
    } else {
        return OPENST_ERR_DIV_BY_ZERO;
    }

    *DSTI = vali;
    *DSTJ = valj;
    *DSTK = valk;

    return OPENST_ERR_SUCCESS;
}


OPENST_ERR OpenST_BRT3D_Trace(double *T, double *V,
                size_t NI, size_t NJ, size_t NK,
                double HI, double HJ, double HK, double TSTEP,
                double RCVI, double RCVJ, double RCVK,
                double SRCI, double SRCJ, double SRCK,
                double **RAY, size_t *RAY_NI, size_t *RAY_NJ){

    OPENST_ERR errcode;
    size_t ind_src_i, ind_src_j, ind_src_k;
    size_t ind_cur_i, ind_cur_j, ind_cur_k;
    double CUR[3];
    double DST[3];
    struct OpenST_DYNARR arr;
    struct OpenST_DYNARR *arrptr = &arr;

    if((errcode = OpenST_CRS_Cart2Ind(SRCI, HI, &ind_src_i))){
        goto EXIT;
    }
    if((errcode = OpenST_CRS_Cart2Ind(SRCJ, HJ, &ind_src_j))){
        goto EXIT;
    }
    if((errcode = OpenST_CRS_Cart2Ind(SRCK, HK, &ind_src_k))){
        goto EXIT;
    }

    CUR[0] = RCVI;
    CUR[1] = RCVJ;
    CUR[2] = RCVK;

    if(OpenST_DYNARR_Init(arrptr, NI + NJ + NK, sizeof(double) * 3) == NULL){
        errcode = OPENST_ERR_MEMORY_ALLOC;
        goto EXIT;
    }

    while(1){
        if(CUR[0] < 0 || CUR[1] < 0 || CUR[2] < 0){
            errcode = OPENST_ERR_ALGORITHM;
            break;
        }

        if((errcode = OpenST_CRS_Cart2Ind(CUR[0], HI, &ind_cur_i))){
            goto EXIT;
        }
        if((errcode = OpenST_CRS_Cart2Ind(CUR[1], HJ, &ind_cur_j))){
            goto EXIT;
        }
        if((errcode = OpenST_CRS_Cart2Ind(CUR[2], HK, &ind_cur_k))){
            goto EXIT;
        }

#if DEBUG_LOG
        printf("[%zu;%zu;%zu] - [%e;%e;%e]\n",ind_cur_i,ind_cur_j,ind_cur_k,CUR[0],CUR[1],CUR[2]);
        fflush(stdout);
#endif

        if(ind_cur_i == ind_src_i
            && ind_cur_j == ind_src_j && ind_cur_k == ind_src_k){
            errcode = OPENST_ERR_SUCCESS;
#if DEBUG_LOG
        printf("SOURCE REACHED\n");
        fflush(stdout);
#endif
            break;
        }

        if(ind_cur_i >= NI || ind_cur_j >= NJ || ind_cur_k >= NK){
            errcode = OPENST_ERR_ALGORITHM;
            break;
        }

        if(OpenST_DYNARR_Pushback(arrptr, CUR) == NULL){
            errcode = OPENST_ERR_MEMORY_ALLOC;
            goto EXIT;
        }

        if((errcode = OpenST_BRT3D_Step(T, V, NI, NJ, NK, HI, HJ, HK, TSTEP,
                         CUR[0], CUR[1], CUR[2],
                         ind_cur_i, ind_cur_j, ind_cur_k,
                         &DST[0], &DST[1], &DST[2]))){
            goto EXIT;
        }

        CUR[0] = DST[0];
        CUR[1] = DST[1];
        CUR[2] = DST[2];
    }

    if(OpenST_DYNARR_Shrink(arrptr) == NULL){
        errcode = OPENST_ERR_MEMORY_ALLOC;
        goto EXIT;
    }

    *RAY = arr.data;
    *RAY_NJ = 3;
    *RAY_NI = arr.num;

EXIT:
    return errcode;
}
