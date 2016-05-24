#ifndef OPENST_COMMON_ARRAYOPS_H
#define OPENST_COMMON_ARRAYOPS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <float.h>
#include <stddef.h>
#include <math.h>

#include "openst/common/macros.h"

OPENST_API void OpenST_AOP_GetArrStats(double *A, size_t numel, double *min,
                        double *max, double *mean);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_ARRAYOPS_H */
