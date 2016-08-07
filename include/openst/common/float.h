#ifndef OPENST_COMMON_FLOAT_H
#define OPENST_COMMON_FLOAT_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include <math.h>
#include <stddef.h>
#include <float.h>

#include "openst/common/macros.h"

OPENST_API int OpenST_FLOAT_ApproximatelyEqual(double a, double b, double epsilon);

OPENST_API int OpenST_FLOAT_EssentiallyEqual(double a, double b, double epsilon);

OPENST_API int OpenST_FLOAT_DefinitelyGreater(double a, double b, double epsilon);

OPENST_API int OpenST_FLOAT_DefinitelyLess(double a, double b, double epsilon);

OPENST_API int OpenST_FLOAT_GetNeighboorSizeT(double real, size_t *left, size_t *right);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_FLOAT_H */
