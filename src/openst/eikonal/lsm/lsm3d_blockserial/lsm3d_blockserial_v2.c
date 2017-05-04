#include "openst/eikonal/lsm.h"

#define M_LSM3D_BLOCKSERIAL_IMP_NAME "lsm3d_blockserial_v2.c"


const char OPENST_LSM3D_BLOCKSERIAL_IMP_NAME[] = M_LSM3D_BLOCKSERIAL_IMP_NAME;
const size_t OPENST_LSM3D_BLOCKSERIAL_IMP_NAME_LENGTH = sizeof(OPENST_LSM3D_BLOCKSERIAL_IMP_NAME);


int OpenST_LSM3D_BlockSerial_1H(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                size_t NI, size_t NJ, size_t NK,
                                OPENST_FLOAT H,
                                int REVI, int REVJ, int REVK,
                                size_t istart, size_t jstart, size_t kstart,
                                size_t isize, size_t jsize, size_t ksize,
                                OPENST_FLOAT EPS){

    int notconverged;
    size_t i, j, ir, jr, kr, iend, jend, kend;

    iend = OPENST_HACK_MIN(istart + isize, NI);
    jend = OPENST_HACK_MIN(jstart + jsize, NJ);
    kend = OPENST_HACK_MIN(kstart + ksize, NK);

    notconverged = 0;
    for(ir = istart; ir < iend; ir += 2u){
        for(jr = jstart; jr < jend; jr += 2u){
            for(i = ir; (i < (ir + 2u)) && (i < iend); ++i){
                for(j = jr; (j < (jr + 2u)) && (j < jend); ++j){
                    for(kr = kstart; kr < kend; ++kr){
                        if(OpenST_LSM3D_NodeUpdate_1H(U, LSM_UNLOCKED, V,
                                                      NI, NJ, NK,
                                                      H,
                                                      REVI, REVJ, REVK,
                                                      i, j, kr, EPS)){
                            notconverged = 1;
                        }
                    }
                }
            }
        }
    }

    return notconverged;
}


int OpenST_LSM3D_BlockSerial(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                             size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             int REVI, int REVJ, int REVK,
                             size_t istart, size_t jstart, size_t kstart,
                             size_t isize, size_t jsize, size_t ksize,
                             OPENST_FLOAT EPS){

    int notconverged;
    size_t i, j, ir, jr, kr, iend, jend, kend;

    if((HI == HJ) && (HI == HK)){
        return OpenST_LSM3D_BlockSerial_1H(U, LSM_UNLOCKED, V,
                                           NI, NJ, NK,
                                           HI,
                                           REVI, REVJ, REVK,
                                           istart, jstart, kstart,
                                           isize, jsize, ksize, EPS);
    }

    iend = OPENST_HACK_MIN(istart + isize, NI);
    jend = OPENST_HACK_MIN(jstart + jsize, NJ);
    kend = OPENST_HACK_MIN(kstart + ksize, NK);

    notconverged = 0;
    for(ir = istart; ir < iend; ir += 2u){
        for(jr = jstart; jr < jend; jr += 2u){
            for(i = ir; (i < (ir + 2u)) && (i < iend); ++i){
                for(j = jr; (j < (jr + 2u)) && (j < jend); ++j){
                    for(kr = kstart; kr < kend; ++kr){
                        if(OpenST_LSM3D_NodeUpdate(U, LSM_UNLOCKED, V,
                                                   NI, NJ, NK,
                                                   HI, HJ, HK,
                                                   REVI, REVJ, REVK,
                                                   i, j, kr, EPS)){
                            notconverged = 1;
                        }
                    }
                }
            }
        }
    }

    return notconverged;
}
