#ifndef OPENST_MEX_H
#define OPENST_MEX_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mex.h"
#include "matrix.h"
#include <math.h>

#define OPENST_MEX_ERRID __FUNCTION__##":Error"

#define OPENST_MEX_CHECK(x,MSG) \
    do { \
    OPENST_MEX_ERRCODE errcode = (x); \
    if (errcode != 0) { \
    OPENST_MEX_PrintError(errcode, OPENST_MEX_ERRID, \
    __FUNCTION__, __LINE__, MSG); \
} \
} while (0)

typedef enum OPENST_MEX_ERRCODE_enum
{
    OPENST_MEX_ERRCODE_SUCCESS = 0,
    OPENST_MEX_ERRCODE_MEMORY,
    OPENST_MEX_ERRCODE_NOT_REAL,
    OPENST_MEX_ERRCODE_WRONG_TYPE,
    OPENST_MEX_ERRCODE_WRONG_DIMS
} OPENST_MEX_ERRCODE;

void OPENST_MEX_PrintError(OPENST_MEX_ERRCODE errcode, char* errid, 
                           char *file, int line, char *MSG);

mwSize OPENST_MEX_GetNumel(mwSize ndims, mwSize *dims);

int OPENST_MEX_ValueIsInteger(OPENST_FLOAT value);

OPENST_MEX_ERRCODE OPENST_MEX_GetDoubleArray(mxArray *INP, 
                                             OPENST_FLOAT **ptr, mwSize *ndims,
                                             mwSize **dims);

OPENST_MEX_ERRCODE OPENST_MEX_GetDoubleScalar(mxArray *INP, OPENST_FLOAT *out);

int OPENST_MEX_DoubleIsInteger(OPENST_FLOAT);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_MEX_H */
