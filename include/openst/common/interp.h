#ifndef OPENST_COMMON_INTERP_H
#define OPENST_COMMON_INTERP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#include "openst/common/float.h"
#include "openst/common/macros.h"
#include "openst/common/memadr.h"
#include "openst/common/error.h"
#include "openst/common/coordsys.h"

typedef enum OPENST_INTERP_METHOD_enum{
    OPENST_INTERP_LINEAR,
    OPENST_INTERP_DEFAULT = OPENST_INTERP_LINEAR
} OPENST_INTERP_METHOD;

OPENST_API OPENST_ERR OpenST_INTERP_3D(OPENST_FLOAT *A,
                              size_t NI, size_t NJ, size_t NK,
                              OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                              OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
                              OPENST_FLOAT *VAL);

OPENST_API void OpenST_INTERP_Trilinear(OPENST_FLOAT *A,
                             size_t NI, size_t NJ, size_t NK,
                             OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                             OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
                             OPENST_FLOAT *VAL);

OPENST_API void OpenST_INTERP_Linear_Formula(OPENST_FLOAT *f,
                                    OPENST_FLOAT i,
                                    OPENST_FLOAT f0, OPENST_FLOAT f1,
                                    OPENST_FLOAT i0, OPENST_FLOAT i1);

OPENST_API void OpenST_INTERP_Bilinear_Formula(OPENST_FLOAT *f,
                                      OPENST_FLOAT i, OPENST_FLOAT j,
                                      OPENST_FLOAT f00, OPENST_FLOAT f01,
                                      OPENST_FLOAT f10, OPENST_FLOAT f11,
                                      OPENST_FLOAT i0, OPENST_FLOAT j0,
                                      OPENST_FLOAT i1, OPENST_FLOAT j1);

OPENST_API void OpenST_INTERP_Trilinear_Formula(OPENST_FLOAT *f,
                                       OPENST_FLOAT i, OPENST_FLOAT j, OPENST_FLOAT k,
                                       OPENST_FLOAT f000, OPENST_FLOAT f001,
                                       OPENST_FLOAT f010, OPENST_FLOAT f011,
                                       OPENST_FLOAT f100, OPENST_FLOAT f101,
                                       OPENST_FLOAT f110, OPENST_FLOAT f111,
                                       OPENST_FLOAT i0, OPENST_FLOAT j0, OPENST_FLOAT k0,
                                       OPENST_FLOAT i1, OPENST_FLOAT j1, OPENST_FLOAT k1);

OPENST_API void OpenST_INTERP_Trilinear_Neighboors(
    OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
    OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
    size_t ii[2], size_t ji[2], size_t ki[2],
    OPENST_FLOAT ic[2], OPENST_FLOAT jc[2], OPENST_FLOAT kc[2],
    int *interp_i, int *interp_j, int *interp_k, int *interp_dims);

OPENST_API void OpenST_INTERP_Trilinear_Compute(
    OPENST_FLOAT *A,
    size_t NI, size_t NJ, size_t NK,
    OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
    size_t ii[2], size_t ji[2], size_t ki[2],
    OPENST_FLOAT ic[2], OPENST_FLOAT jc[2], OPENST_FLOAT kc[2],
    int interp_i, int interp_j, int interp_k, int interp_dims,
    OPENST_FLOAT *VAL);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_INTERP_H */
