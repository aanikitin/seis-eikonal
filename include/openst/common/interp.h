#ifndef OPENST_COMMON_INTERP_H
#define OPENST_COMMON_INTERP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#include "openst/common/macros.h"
#include "openst/common/memadr.h"
#include "openst/common/error.h"
#include "openst/common/float.h"

typedef enum OPENST_INTERP_METHOD_enum{
    OPENST_INTERP_LINEAR,
    OPENST_INTERP_DEFAULT = OPENST_INTERP_LINEAR
} OPENST_INTERP_METHOD;

OPENST_API void OpenST_INTERP_Trilinear(double *A,
                             size_t NI, size_t NJ, size_t NK,
                             double HI, double HJ, double HK,
                             double PI, double PJ, double PK,
                             double *VAL);

OPENST_API void OpenST_INTERP_Linear_Formula(double *f,
                                    double i,
                                    double f0, double f1,
                                    double i0, double i1);

OPENST_API void OpenST_INTERP_Bilinear_Formula(double *f,
                                      double i, double j,
                                      double f00, double f01,
                                      double f10, double f11,
                                      double i0, double j0,
                                      double i1, double j1);

OPENST_API void OpenST_INTERP_Trilinear_Formula(double *f,
                                       double i, double j, double k,
                                       double f000, double f001,
                                       double f010, double f011,
                                       double f100, double f101,
                                       double f110, double f111,
                                       double i0, double j0, double k0,
                                       double i1, double j1, double k1);

OPENST_API void OpenST_INTERP_Trilinear_Neighboors(
    double HI, double HJ, double HK,
    double PI, double PJ, double PK,
    size_t ii[2], size_t ji[2], size_t ki[2],
    double ic[2], double jc[2], double kc[2],
    int *interp_i, int *interp_j, int *interp_k, int *interp_dims);

OPENST_API void OpenST_INTERP_Trilinear_Compute(
    double *A,
    size_t NI, size_t NJ, size_t NK,
    double PI, double PJ, double PK,
    size_t ii[2], size_t ji[2], size_t ki[2],
    double ic[2], double jc[2], double kc[2],
    int interp_i, int interp_j, int interp_k, int interp_dims,
    double *VAL);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_INTERP_H */
