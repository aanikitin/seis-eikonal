#include "openst/eikonal/lsm.h"


OPENST_ERR OpenST_LSM3D_Init_2(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                             size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                             size_t **SRCidx, size_t *SRCidx_NI,
                             size_t *SRCidx_NJ,
                             OPENST_FSM3D_INIT_METHOD method){
    OPENST_ERR errcode = OPENST_ERR_SUCCESS;
    size_t srcind, i, j, k;
    size_t *SRCidx_loc, SRCidx_NI_loc, SRCidx_NJ_loc;
    int output_srcidx;
    
    if((SRCidx != NULL) && (SRCidx_NI != NULL) && (SRCidx_NJ != NULL)){
        output_srcidx = 1;
    } else {
        output_srcidx = 0;
    }

    errcode = OpenST_FSM3D_Init_2(U, V, NI, NJ, NK, HI, HJ, HK,
                                SRCI, SRCJ, SRCK,
                                &SRCidx_loc, &SRCidx_NI_loc, &SRCidx_NJ_loc,
                                method);

    for(i = 0; i < NI; ++i){
        for(j = 0; j < NJ; ++j){
            for(k = 0; k < NK; ++k){
                LSM_UNLOCKED[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = 0;
            }
        }
    }

    for(srcind = 0; srcind < SRCidx_NI_loc; ++srcind){
        i = SRCidx_loc[OPENST_MEMADR_2D(srcind,0u,SRCidx_NI_loc,SRCidx_NJ_loc)];
        j = SRCidx_loc[OPENST_MEMADR_2D(srcind,1u,SRCidx_NI_loc,SRCidx_NJ_loc)];
        k = SRCidx_loc[OPENST_MEMADR_2D(srcind,2u,SRCidx_NI_loc,SRCidx_NJ_loc)];

        if(i > 0){
            LSM_UNLOCKED[OPENST_MEMADR_3D(i - 1,j,k,NI,NJ,NK)] = 1;
        }
        if(i < NI - 1){
            LSM_UNLOCKED[OPENST_MEMADR_3D(i + 1,j,k,NI,NJ,NK)] = 1;
        }
        if(j > 0){
            LSM_UNLOCKED[OPENST_MEMADR_3D(i,j - 1,k,NI,NJ,NK)] = 1;
        }
        if(j < NJ - 1){
            LSM_UNLOCKED[OPENST_MEMADR_3D(i,j + 1,k,NI,NJ,NK)] = 1;
        }
        if(k > 0){
            LSM_UNLOCKED[OPENST_MEMADR_3D(i,j,k - 1,NI,NJ,NK)] = 1;
        }
        if(k < NK - 1){
            LSM_UNLOCKED[OPENST_MEMADR_3D(i,j,k + 1,NI,NJ,NK)] = 1;
        }
    }

    for(srcind = 0; srcind < SRCidx_NI_loc; ++srcind){
        i = SRCidx_loc[OPENST_MEMADR_2D(srcind,0u,SRCidx_NI_loc,SRCidx_NJ_loc)];
        j = SRCidx_loc[OPENST_MEMADR_2D(srcind,1u,SRCidx_NI_loc,SRCidx_NJ_loc)];
        k = SRCidx_loc[OPENST_MEMADR_2D(srcind,2u,SRCidx_NI_loc,SRCidx_NJ_loc)];

        LSM_UNLOCKED[OPENST_MEMADR_3D(i,j,k,NI,NJ,NK)] = 0;
    }

    if(output_srcidx){
        *SRCidx = SRCidx_loc;
        *SRCidx_NI = SRCidx_NI_loc;
        *SRCidx_NJ = SRCidx_NJ_loc;
    }

    return errcode;
}


OPENST_ERR OpenST_LSM3D_Init(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                             size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK){
    return OpenST_LSM3D_Init_2(U,LSM_UNLOCKED,V,
                             NI,NJ,NK,
                             HI,HJ,HK,
                             SRCI,SRCJ,SRCK,
                             NULL,NULL,NULL,
                             OPENST_FSM3D_INIT_DEFAULT);
}
