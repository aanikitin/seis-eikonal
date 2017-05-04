#ifndef OPENST_COMMON_ARRAYOPS_H
#define OPENST_COMMON_ARRAYOPS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <float.h>
#include <stddef.h>
#include <math.h>

#include "openst/common/float.h"
#include "openst/common/hacks.h"
#include "openst/common/macros.h"

OPENST_API void OpenST_AOP_GetArrStats(OPENST_FLOAT *A, size_t numel, OPENST_FLOAT *min,
                        OPENST_FLOAT *max, OPENST_FLOAT *mean);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_ARRAYOPS_H */
