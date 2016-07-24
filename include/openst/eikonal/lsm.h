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

#include "openst/common/memadr.h"
#include "openst/common/macros.h"
#include "openst/common/error.h"
#include "openst/eikonal/fsm.h"

OPENST_API extern const char OPENST_LSM3D_IMP_NAME[];
OPENST_API extern const size_t OPENST_LSM3D_IMP_NAME_LENGTH;

OPENST_API OPENST_ERR OpenST_LSM3D_Init(double *U, char *LSM_UNLOCKED, double *V,
                                        size_t NI, size_t NJ, size_t NK,
                                        double HI, double HJ, double HK,
                                        double SRCI, double SRCJ, double SRCK,
                                        size_t **SRCidx,
                                        size_t *SRCidx_NI, size_t *SRCidx_NJ,
                                        OPENST_FSM3D_INIT_METHOD method);

OPENST_API int OpenST_LSM3D_Compute(double *U, char *LSM_UNLOCKED, double *V,
                                    size_t NI, size_t NJ, size_t NK,
                                    double HI, double HJ, double HK,
                                    int max_iter, int *converged,
                                    size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                    double EPS);

OPENST_API int OpenST_LSM3D_ComputePartial(double *U, char *LSM_UNLOCKED, double *V,
                                           size_t NI, size_t NJ, size_t NK,
                                           double HI, double HJ, double HK,
                                           int start_iter, int max_iter, int *converged,
                                           size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                           double EPS);

OPENST_API int OpenST_LSM3D_BlockSerial(double *U, char *LSM_UNLOCKED, double *V,
                                        size_t NI, size_t NJ, size_t NK,
                                        double HI, double HJ, double HK,
                                        int REVI, int REVJ, int REVK,
                                        size_t istart, size_t jstart, size_t kstart,
                                        size_t isize, size_t jsize, size_t ksize,
                                        double EPS);

OPENST_API int OpenST_LSM3D_NodeUpdate(double *U, char *LSM_UNLOCKED, double *V,
                                       size_t NI, size_t NJ, size_t NK,
                                       double HI, double HJ, double HK,
                                       int REVI, int REVJ, int REVK,
                                       size_t ir, size_t jr, size_t kr,
                                       double EPS);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_EIKONAL_LSM_H */
