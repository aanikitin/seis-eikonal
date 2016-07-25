#ifndef OPENST_EIKONAL_FSM_H
#define OPENST_EIKONAL_FSM_H

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
#include "openst/common/coordsys.h"
#include "openst/common/interp.h"

OPENST_API extern const char OPENST_FSM3D_IMP_NAME[];
OPENST_API extern const size_t OPENST_FSM3D_IMP_NAME_LENGTH;

typedef enum OPENST_FSM3D_INIT_METHOD_enum{
    OPENST_FSM3D_INIT_POINT,
    OPENST_FSM3D_INIT_LINEAR_INTERP,
    OPENST_FSM3D_INIT_DEFAULT = OPENST_FSM3D_INIT_LINEAR_INTERP
} OPENST_FSM3D_INIT_METHOD;

OPENST_API OPENST_ERR OpenST_FSM3D(double *U, double *V,
                                   size_t NI, size_t NJ, size_t NK,
                                   double HI, double HJ, double HK,
                                   double SRCI, double SRCJ, double SRCK,
                                   double EPS, int max_iter,
                                   int *it, int *converged);

OPENST_API OPENST_ERR OpenST_FSM3D_Init(double *U, double *V,
                                        size_t NI, size_t NJ, size_t NK,
                                        double HI, double HJ, double HK,
                                        double SRCI, double SRCJ, double SRCK);

OPENST_API OPENST_ERR OpenST_FSM3D_Init_2(double *U, double *V,
                                          size_t NI, size_t NJ, size_t NK,
                                          double HI, double HJ, double HK,
                                          double SRCI, double SRCJ, double SRCK,
                                          size_t **SRCidx, size_t *SRCidx_NI,
                                          size_t *SRCidx_NJ,
                                          OPENST_FSM3D_INIT_METHOD method);

OPENST_API void OpenST_FSM3D_SuggestBlockSize(size_t NI, size_t NJ, size_t NK,
                                              size_t *BSIZE_I, size_t *BSIZE_J,
                                              size_t *BSIZE_K);

OPENST_API int OpenST_FSM3D_Compute(double *U, double *V,
                                    size_t NI, size_t NJ, size_t NK,
                                    double HI, double HJ, double HK,
                                    int max_iter,
                                    int *converged, size_t BSIZE_I, size_t BSIZE_J,
                                    size_t BSIZE_K, double EPS);

OPENST_API int OpenST_FSM3D_ComputePartial(double *U, double *V,
                                           size_t NI, size_t NJ, size_t NK,
                                           double HI, double HJ, double HK,
                                           int start_iter, int max_iter, int *converged,
                                           size_t BSIZE_I, size_t BSIZE_J, size_t BSIZE_K,
                                           double EPS);

OPENST_API int OpenST_FSM3D_BlockSerial(double *U, double *V,
                                        size_t NI, size_t NJ, size_t NK,
                                        double HI, double HJ, double HK,
                                        int REVI, int REVJ, int REVK,
                                        size_t istart, size_t jstart, size_t kstart,
                                        size_t isize, size_t jsize, size_t ksize,
                                        double EPS);

OPENST_API int OpenST_FSM3D_NodeUpdate(double *U, double *V,
                                       size_t NI, size_t NJ, size_t NK,
                                       double HI, double HJ, double HK,
                                       int REVI, int REVJ, int REVK,
                                       size_t ir, size_t jr, size_t kr, double EPS);

OPENST_API void OpenST_FSM3D_GetSweepOrder(int it, int *REVI, int *REVJ, int *REVK);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_EIKONAL_FSM_H */
