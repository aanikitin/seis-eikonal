/* TODO: add more accurate initialization schemes */
#include "openst/eikonal/fsm.h"


OPENST_ERR OpenST_FSM3D_InitSRC_Point(OPENST_FLOAT *U, OPENST_FLOAT *V,
                                      size_t NI, size_t NJ, size_t NK,
                                      OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                      OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                                      size_t **SRCidx, size_t *SRCidx_NI,
                                      size_t *SRCidx_NJ){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t ii, ji, ki;
    OPENST_FLOAT ic, jc, kc;
    OPENST_FLOAT dist;
    int output_srcidx;
    size_t *SRCidx_loc = NULL;
    size_t SRCidx_NI_loc;
    size_t SRCidx_NJ_loc;

    if ( OpenST_CRS_IsPointNotWithinBounds(SRCI, SRCJ, SRCK,
                                        NI, NJ, NK,
                                        HI, HJ, HK) ) {
        errcode = OPENST_ERR_PARAM_INVALID;
        goto EXIT;
    }

    if((SRCidx != NULL) && (SRCidx_NI != NULL) && (SRCidx_NJ != NULL)){
        output_srcidx = 1;
        SRCidx_NI_loc = 1;
        SRCidx_NJ_loc = 3;
        SRCidx_loc = (size_t *)malloc(sizeof(size_t) * SRCidx_NI_loc
                                      * SRCidx_NJ_loc);
        if(SRCidx_loc == NULL){
            errcode = OPENST_ERR_MEMORY_ALLOC;
            goto EXIT;
        }
    } else {
        output_srcidx = 0;
    }

    ii = (size_t) OPENST_FLOAT_ROUND(SRCI / HI);
    ic = (OPENST_FLOAT)ii * HI;
    
    ji = (size_t) OPENST_FLOAT_ROUND(SRCJ / HJ);
    jc = (OPENST_FLOAT)ji * HJ;
    
    ki = (size_t) OPENST_FLOAT_ROUND(SRCK / HK);
    kc = (OPENST_FLOAT)ki * HK;

    dist = OpenST_CRS_Distance3D(SRCI, SRCJ, SRCK, ic, jc, kc);

    U[OPENST_MEMADR_3D(ii, ji, ki, NI, NJ, NK)] =
            dist / V[OPENST_MEMADR_3D(ii, ji, ki, NI, NJ, NK)];

    if(output_srcidx){
        SRCidx_loc[0] = ii;
        SRCidx_loc[1] = ji;
        SRCidx_loc[2] = ki;

        *SRCidx = SRCidx_loc;
        *SRCidx_NI = SRCidx_NI_loc;
        *SRCidx_NJ = SRCidx_NJ_loc;
    }

EXIT:
    return errcode;
}


OPENST_ERR OpenST_FSM3D_InitSRC_Linear(OPENST_FLOAT *U, OPENST_FLOAT *V,
                                       size_t NI, size_t NJ, size_t NK,
                                       OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                       OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                                       size_t **SRCidx,
                                       size_t *SRCidx_NI, size_t *SRCidx_NJ){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    OPENST_FLOAT SRCV;
    OPENST_FLOAT dist;
    size_t ii[2], ji[2], ki[2];
    OPENST_FLOAT ic[2], jc[2], kc[2];
    int i, j, k, interp_dims, interp_i, interp_j, interp_k;
    int output_srcidx;
    size_t *SRCidx_loc = NULL;
    size_t SRCidx_NI_loc;
    size_t SRCidx_NJ_loc;
    size_t SRCidx_ind;

    if ( OpenST_CRS_IsPointNotWithinBounds(SRCI, SRCJ, SRCK,
                                        NI, NJ, NK,
                                        HI, HJ, HK) ) {
        errcode = OPENST_ERR_PARAM_INVALID;
        goto EXIT;
    }

    if((SRCidx != NULL) && (SRCidx_NI != NULL) && (SRCidx_NJ != NULL)){
        output_srcidx = 1;
    } else {
        output_srcidx = 0;
    }

    OpenST_INTERP_Trilinear_Neighboors(
                HI, HJ, HK,
                SRCI, SRCJ, SRCK,
                ii, ji, ki,
                ic, jc, kc,
                &interp_i, &interp_j, &interp_k, &interp_dims);

    if(output_srcidx){
        SRCidx_NI_loc = 1;
        SRCidx_NI_loc <<= interp_dims;
        SRCidx_NJ_loc = 3;
        SRCidx_loc = (size_t *)malloc(sizeof(size_t) * SRCidx_NI_loc
                                      * SRCidx_NJ_loc);
        if(SRCidx_loc == NULL){
            errcode = OPENST_ERR_MEMORY_ALLOC;
            goto EXIT;
        }
    }

    if(interp_dims == 0){

        U[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)] = OPENST_FLOAT_0_0;

        if(output_srcidx){
            SRCidx_loc[0] = ii[0];
            SRCidx_loc[1] = ji[0];
            SRCidx_loc[2] = ki[0];
        }

    } else {

        OpenST_INTERP_Trilinear_Compute(
                    V,
                    NI, NJ, NK,
                    SRCI, SRCJ, SRCK,
                    ii, ji, ki,
                    ic, jc, kc,
                    interp_i, interp_j, interp_k, interp_dims, &SRCV);

        SRCidx_ind = 0;
        for(i = 0; i < interp_i + 1; ++i){
            for(j = 0; j < interp_j + 1; ++j){
                for(k = 0; k < interp_k + 1; ++k){

                    dist = OpenST_CRS_Distance3D(SRCI, SRCJ, SRCK,
                                                 ic[i], jc[j], kc[k]);

                    U[OPENST_MEMADR_3D(ii[i], ji[j], ki[k], NI, NJ, NK)] =
                            dist / SRCV;

                    if(output_srcidx){
                        SRCidx_loc[OPENST_MEMADR_2D(SRCidx_ind,0,SRCidx_NI_loc,
                                                    SRCidx_NJ_loc)] = ii[i];
                        SRCidx_loc[OPENST_MEMADR_2D(SRCidx_ind,1,SRCidx_NI_loc,
                                                    SRCidx_NJ_loc)] = ji[j];
                        SRCidx_loc[OPENST_MEMADR_2D(SRCidx_ind,2,SRCidx_NI_loc,
                                                    SRCidx_NJ_loc)] = ki[k];
                        ++SRCidx_ind;
                    }
                }
            }
        }

    }

    if(output_srcidx){
        *SRCidx = SRCidx_loc;
        *SRCidx_NI = SRCidx_NI_loc;
        *SRCidx_NJ = SRCidx_NJ_loc;
    }

EXIT:
    return errcode;
}


OPENST_ERR OpenST_FSM3D_Init_2(OPENST_FLOAT *U, OPENST_FLOAT *V,
                               size_t NI, size_t NJ, size_t NK,
                               OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                               OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                               size_t **SRCidx, size_t *SRCidx_NI,
                               size_t *SRCidx_NJ,
                               OPENST_FSM3D_INIT_METHOD method){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t i, j, k;

    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = OPENST_FLOAT_INF;
            }
        }
    }

    switch(method){
    case OPENST_FSM3D_INIT_LINEAR_INTERP:
        errcode = OpenST_FSM3D_InitSRC_Linear(U, V, NI, NJ, NK,
                                              HI, HJ, HK,
                                              SRCI, SRCJ, SRCK,
                                              SRCidx, SRCidx_NI, SRCidx_NJ);
        break;
    case OPENST_FSM3D_INIT_POINT:
        errcode = OpenST_FSM3D_InitSRC_Point(U, V, NI, NJ, NK,
                                             HI, HJ, HK,
                                             SRCI, SRCJ, SRCK,
                                             SRCidx, SRCidx_NI, SRCidx_NJ);
        break;
    default:
        errcode = OPENST_ERR_PARAM_INVALID;
    }

    return errcode;
}


OPENST_ERR OpenST_FSM3D_Init(OPENST_FLOAT *U, OPENST_FLOAT *V,
                             size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK){
    return OpenST_FSM3D_Init_2(U,V,
                               NI,NJ,NK,
                               HI,HJ,HK,
                               SRCI,SRCJ,SRCK,
                               NULL,NULL,NULL,
                               OPENST_FSM3D_INIT_DEFAULT);
}
