#ifndef OPENST_COMMON_INTERP_H
#define OPENST_COMMON_INTERP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

#include "openst/common/macros.h"

OPENST_API void OpenST_INTERP_Linear(double *f, double i,
                                    double f0, double f1,
                                    double i0, double i1);

OPENST_API void OpenST_INTERP_Bilinear(double *f, double i, double j,
                                      double f00, double f01,
                                      double f10, double f11,
                                      double i0, double j0,
                                      double i1, double j1);

OPENST_API void OpenST_INTERP_Trilinear(double *f,
                                       double i, double j, double k,
                                       double f000, double f001,
                                       double f010, double f011,
                                       double f100, double f101,
                                       double f110, double f111,
                                       double i0, double j0, double k0,
                                       double i1, double j1, double k1);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_INTERP_H */
