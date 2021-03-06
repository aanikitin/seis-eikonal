// WARNING: Will suffer significant performance drop if link-time optimization 
// is not enabled
#ifndef OPENST_EIKONAL_LSM_H
#define OPENST_EIKONAL_LSM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <omp.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <float.h>

#include "openst/common/float.h"
#include "openst/common/memadr.h"
#include "openst/common/macros.h"
#include "openst/common/error.h"
#include "openst/eikonal/fsm.h"

OPENST_API extern const char OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME[];
OPENST_API extern const size_t OPENST_LSM3D_COMPUTEPARTIAL_IMP_NAME_LENGTH;

OPENST_API extern const char OPENST_LSM3D_BLOCKSERIAL_IMP_NAME[];
OPENST_API extern const size_t OPENST_LSM3D_BLOCKSERIAL_IMP_NAME_LENGTH;

OPENST_API OPENST_ERR OpenST_LSM3D(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                   size_t NI, size_t NJ, size_t NK,
                                   OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                   OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                                   OPENST_FLOAT EPS, int max_iter,
                                   int *it, int *converged);

OPENST_API OPENST_ERR OpenST_LSM3D_Init(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                        size_t NI, size_t NJ, size_t NK,
                                        OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                        OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK);

OPENST_API OPENST_ERR OpenST_LSM3D_Init_2(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                          size_t NI, size_t NJ, size_t NK,
                                          OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                          OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK,
                                          size_t **SRCidx, size_t *SRCidx_NI,
                                          size_t *SRCidx_NJ,
                                          OPENST_FSM3D_INIT_METHOD method);

OPENST_API int OpenST_LSM3D_Compute(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                    size_t NI, size_t NJ, size_t NK,
                                    OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                    int max_iter, int *converged,
                                    size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                    OPENST_FLOAT EPS);

OPENST_API int OpenST_LSM3D_ComputePartial(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                           size_t NI, size_t NJ, size_t NK,
                                           OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                           int start_iter, int max_iter, int *converged,
                                           size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                           OPENST_FLOAT EPS);

OPENST_API int OpenST_LSM3D_BlockSerial(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                        size_t NI, size_t NJ, size_t NK,
                                        OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                        int REVI, int REVJ, int REVK,
                                        size_t istart, size_t jstart, size_t kstart,
                                        size_t isize, size_t jsize, size_t ksize,
                                        OPENST_FLOAT EPS);

OPENST_API int OpenST_LSM3D_BlockSerial_1H(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                           size_t NI, size_t NJ, size_t NK,
                                           OPENST_FLOAT H,
                                           int REVI, int REVJ, int REVK,
                                           size_t istart, size_t jstart, size_t kstart,
                                           size_t isize, size_t jsize, size_t ksize,
                                           OPENST_FLOAT EPS);

OPENST_API int OpenST_LSM3D_NodeUpdate(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                       size_t NI, size_t NJ, size_t NK,
                                       OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                                       int REVI, int REVJ, int REVK,
                                       size_t ir, size_t jr, size_t kr,
                                       OPENST_FLOAT EPS);

OPENST_API int OpenST_LSM3D_NodeUpdate_1H(OPENST_FLOAT *U, char *LSM_UNLOCKED, OPENST_FLOAT *V,
                                          size_t NI, size_t NJ, size_t NK,
                                          OPENST_FLOAT H,
                                          int REVI, int REVJ, int REVK,
                                          size_t ir, size_t jr, size_t kr, OPENST_FLOAT EPS);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_EIKONAL_LSM_H */
