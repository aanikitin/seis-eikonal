/* TODO: add more accurate initialization schemes */
#include "openst/eikonal/fsm.h"


OPENST_ERR OpenST_FSM3D_InitSRC_Point(double *U, double *V,
                                      size_t NI, size_t NJ, size_t NK,
                                      double HI, double HJ, double HK,
                                      double SRCI, double SRCJ, double SRCK,
                                      size_t **SRCidx, size_t *SRCidx_NI,
                                      size_t *SRCidx_NJ){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t ii, ji, ki;
    double ic, jc, kc;
    double dist;
    int output_srcidx;
    size_t *SRCidx_loc = NULL;
    size_t SRCidx_NI_loc;
    size_t SRCidx_NJ_loc;

    if((SRCidx != NULL) && (SRCidx_NI != NULL) && (SRCidx_NJ != NULL)){
        output_srcidx = 1;
        SRCidx_NI_loc = 1;
        SRCidx_NJ_loc = 3;
        SRCidx_loc = (size_t *)malloc(sizeof(size_t) * SRCidx_NI_loc
                                      * SRCidx_NJ_loc);
        if(SRCidx_loc == NULL){
            errcode = OPENST_ERR_MEMORY;
            goto EXIT;
        }
    } else {
        output_srcidx = 0;
    }

    ii = (size_t) round(SRCI / HI);
    ic = (double)ii * HI;
    
    ji = (size_t) round(SRCJ / HJ);
    jc = (double)ji * HJ;
    
    ki = (size_t) round(SRCK / HK);
    kc = (double)ki * HK;

    dist = OpenST_CRS_Distance3D(SRCI, SRCJ, SRCK, ic, jc, kc);

    U[OPENST_MEMADR_3D(ii, ji, ki, NI, NJ, NK)] =
            dist * V[OPENST_MEMADR_3D(ii, ji, ki, NI, NJ, NK)];

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


OPENST_ERR OpenST_FSM3D_InitSRC_Linear(double *U, double *V,
                                       size_t NI, size_t NJ, size_t NK,
                                       double HI, double HJ, double HK,
                                       double SRCI, double SRCJ, double SRCK,
                                       size_t **SRCidx,
                                       size_t *SRCidx_NI, size_t *SRCidx_NJ){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t i, j, k;
    double SRCV;
    double dist;
    size_t ii[2], ji[2], ki[2];
    double ic[2], jc[2], kc[2];
    unsigned int interp_dims, interp_i, interp_j, interp_k;
    int output_srcidx;
    size_t *SRCidx_loc = NULL;
    size_t SRCidx_NI_loc;
    size_t SRCidx_NJ_loc;
    size_t SRCidx_ind;

    if((SRCidx != NULL) && (SRCidx_NI != NULL) && (SRCidx_NJ != NULL)){
        output_srcidx = 1;
    } else {
        output_srcidx = 0;
    }

    ii[0] = (size_t) floor(SRCI / HI);
    ic[0] = (double) ii[0] * HI;

    ii[1] = (size_t) ceil(SRCI / HI);
    ic[1] = (double) ii[1] * HI;

    ji[0] = (size_t) floor(SRCJ / HJ);
    jc[0] = (double) ji[0] * HJ;

    ji[1] = (size_t) ceil(SRCJ / HJ);
    jc[1] = (double) ji[1] * HJ;

    ki[0] = (size_t) floor(SRCK / HK);
    kc[0] = (double) ki[0] * HK;

    ki[1] = (size_t) ceil(SRCK / HK);
    kc[1] = (double) ki[1] * HK;

    interp_dims = 0;
    if((interp_i = (ii[0] != ii[1]))){
        ++interp_dims;
    }
    if((interp_j = (ji[0] != ji[1]))){
        ++interp_dims;
    }
    if((interp_k = (ki[0] != ki[1]))){
        ++interp_dims;
    }

    if(output_srcidx){
        SRCidx_NI_loc = 1;
        SRCidx_NI_loc <<= interp_dims;
        SRCidx_NJ_loc = 3;
        SRCidx_loc = (size_t *)malloc(sizeof(size_t) * SRCidx_NI_loc
                                      * SRCidx_NJ_loc);
        if(SRCidx_loc == NULL){
            errcode = OPENST_ERR_MEMORY;
            goto EXIT;
        }
    }

    if(interp_dims == 0){

        U[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)] = 0.0;

        if(output_srcidx){
            SRCidx_loc[0] = ii[0];
            SRCidx_loc[1] = ji[0];
            SRCidx_loc[2] = ki[0];
        }

    } else {

        if (interp_dims == 3) {

            OpenST_INTERP_Trilinear(&SRCV, SRCI, SRCJ, SRCK,
                                    V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[0], ji[1], ki[1], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[1], ji[0], ki[1], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[1], ji[1], ki[0], NI, NJ, NK)],
                    V[OPENST_MEMADR_3D(ii[1], ji[1], ki[1], NI, NJ, NK)],
                    ic[0], jc[0], kc[0],
                    ic[1], jc[1], kc[1]);

        } else if (interp_dims == 2) {

            if (!interp_i) {
                OpenST_INTERP_Bilinear(&SRCV, SRCJ, SRCK,
                                       V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[1], ki[1], NI, NJ, NK)],
                        jc[0], kc[0],
                        jc[1], kc[1]);
            } else if (!interp_j) {
                OpenST_INTERP_Bilinear(&SRCV, SRCI, SRCK,
                                       V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[1], ji[0], ki[1], NI, NJ, NK)],
                        ic[0], kc[0],
                        ic[1], kc[1]);
            } else if (!interp_k) {
                OpenST_INTERP_Bilinear(&SRCV, SRCI, SRCK,
                                       V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[1], ji[1], ki[0], NI, NJ, NK)],
                        ic[0], jc[0],
                        ic[1], jc[1]);
            }

        } else if (interp_dims == 1) {

            if (interp_i) {
                OpenST_INTERP_Linear(&SRCV, SRCI,
                                     V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[1], ji[0], ki[0], NI, NJ, NK)],
                        ic[0], ic[1]);
            } else if (interp_j) {
                OpenST_INTERP_Linear(&SRCV, SRCJ,
                                     V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[1], ki[0], NI, NJ, NK)],
                        jc[0], jc[1]);
            } else if (interp_k) {
                OpenST_INTERP_Linear(&SRCV, SRCK,
                                     V[OPENST_MEMADR_3D(ii[0], ji[0], ki[0], NI, NJ, NK)],
                        V[OPENST_MEMADR_3D(ii[0], ji[0], ki[1], NI, NJ, NK)],
                        kc[0], kc[1]);
            }

        }

        SRCidx_ind = 0;
        for(i = 0; i < interp_i + 1; ++i){
            for(j = 0; j < interp_j + 1; ++j){
                for(k = 0; k < interp_k + 1; ++k){

                    dist = OpenST_CRS_Distance3D(SRCI, SRCJ, SRCK,
                                                 ic[i], jc[j], kc[k]);

                    U[OPENST_MEMADR_3D(ii[i], ji[j], ki[k], NI, NJ, NK)] =
                            dist * SRCV;

                    if(output_srcidx){
                        SRCidx_loc[OPENST_MEMADR_2D(SRCidx_ind,0u,SRCidx_NI_loc,
                                                    SRCidx_NJ_loc)] = ii[i];
                        SRCidx_loc[OPENST_MEMADR_2D(SRCidx_ind,1u,SRCidx_NI_loc,
                                                    SRCidx_NJ_loc)] = ji[j];
                        SRCidx_loc[OPENST_MEMADR_2D(SRCidx_ind,2u,SRCidx_NI_loc,
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


OPENST_ERR OpenST_FSM3D_Init_2(double *U, double *V,
                             size_t NI, size_t NJ, size_t NK,
                             double HI, double HJ, double HK,
                             double SRCI, double SRCJ, double SRCK,
                             size_t **SRCidx, size_t *SRCidx_NI,
                             size_t *SRCidx_NJ,
                             OPENST_FSM3D_INIT_METHOD method){

    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t i, j, k;

    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                U[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = INFINITY;
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


OPENST_ERR OpenST_FSM3D_Init(double *U, double *V,
                             size_t NI, size_t NJ, size_t NK,
                             double HI, double HJ, double HK,
                             double SRCI, double SRCJ, double SRCK){
    return OpenST_FSM3D_Init_2(U,V,
                             NI,NJ,NK,
                             HI,HJ,HK,
                             SRCI,SRCJ,SRCK,
                             NULL,NULL,NULL,
                             OPENST_FSM3D_INIT_DEFAULT);
}
