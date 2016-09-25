#ifndef OPENST_COMMON_HACKS_H
#define OPENST_COMMON_HACKS_H

#ifdef __cplusplus
extern "C" {
#endif

#define OPENST_HACK_MIN(a,b) (((a)<(b))?(a):(b))
#define OPENST_HACK_MAX(a,b) (((a)>(b))?(a):(b))

/* Enable dirty hacks for missing C99 support in old Visual Studio versions */
#ifdef _MSC_VER
#if (_MSC_VER < 1900)
#define OPENST_HACK 1
#pragma message("!!! WARNING !!! OPENST_HACK is enabled! Upgrade to VS 2015.")
#else
#define OPENST_HACK 0
#endif
#endif

#if OPENST_HACK

#include <float.h>
#include <math.h>

#define OPENST_HACK_INFINITY 1
#define OPENST_HACK_FLOAT_FUNC 1

/* Missing C99's INFINITY in older versions of Visual Studio */
#if OPENST_HACK_INFINITY
#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif
#endif

#if OPENST_HACK_FLOAT_FUNC
#define fmin(a,b) ((a) < (b) ? (a) : (b))
#define fmax(a,b) ((b) > (a) ? (b) : (a))
#define round(a) (((floor((a)) + 0.5) < ceil((a))) ? floor((a)) : ceil((a)))
#endif

#endif /* if OPENST_HACK */

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_HACKS_H */
