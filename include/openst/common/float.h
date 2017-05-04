#ifndef OPENST_COMMON_FLOAT_H
#define OPENST_COMMON_FLOAT_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include <math.h>
#include <stddef.h>
#include <float.h>

#include "openst/common/macros.h"

#ifdef OPENST_FLOAT_PRECISION_SINGLE
typedef float OPENST_FLOAT;
#define OPENST_FLOAT_SFX(s) s##f
/* tgmath.h is not supported in older VS */
#define OPENST_FLOAT_FABS(a) fabsf((a))
#define OPENST_FLOAT_FMIN(a,b) fminf((a),(b))
#define OPENST_FLOAT_FMAX(a,b) fmaxf((a),(b))
#define OPENST_FLOAT_FLOOR(a) floorf((a))
#define OPENST_FLOAT_CEIL(a) ceilf((a))
#define OPENST_FLOAT_ROUND(a) roundf((a))
#define OPENST_FLOAT_POW(a,b) powf((a),(b))
#define OPENST_FLOAT_SQRT(a) sqrtf((a))
/* constants */
#define OPENST_FLOAT_INF INFINITY
#define OPENST_FLOAT_0_0 0.0f
#define OPENST_FLOAT_0_5 0.5f
#define OPENST_FLOAT_1_0 1.0f
#define OPENST_FLOAT_2_0 2.0f
#define OPENST_FLOAT_3_0 3.0f
#define OPENST_FLOAT_4_0 4.0f
#else
typedef double OPENST_FLOAT;
#define OPENST_FLOAT_SFX(s) s
/* tgmath.h is not supported in older VS */
#define OPENST_FLOAT_FABS(a) fabs((a))
#define OPENST_FLOAT_FMIN(a,b) fmin((a),(b))
#define OPENST_FLOAT_FMAX(a,b) fmax((a),(b))
#define OPENST_FLOAT_FLOOR(a) floor((a))
#define OPENST_FLOAT_CEIL(a) ceil((a))
#define OPENST_FLOAT_ROUND(a) round((a))
#define OPENST_FLOAT_POW(a,b) pow((a),(b))
#define OPENST_FLOAT_SQRT(a) sqrt((a))
/* constants */
#define OPENST_FLOAT_INF INFINITY
#define OPENST_FLOAT_0_0 0.0
#define OPENST_FLOAT_0_5 0.5
#define OPENST_FLOAT_1_0 1.0
#define OPENST_FLOAT_2_0 2.0
#define OPENST_FLOAT_3_0 3.0
#define OPENST_FLOAT_4_0 4.0
#endif

OPENST_API int OpenST_FLOAT_ApproximatelyEqual(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon);

OPENST_API int OpenST_FLOAT_EssentiallyEqual(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon);

OPENST_API int OpenST_FLOAT_DefinitelyGreater(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon);

OPENST_API int OpenST_FLOAT_DefinitelyLess(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon);

OPENST_API int OpenST_FLOAT_GetNeighboorSizeT(OPENST_FLOAT real, size_t *left, size_t *right);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_FLOAT_H */
