/* Copyright 2014-2017 Alexandr Nikitin. */

#include "mex.h"
#include "matrix.h"
#include "OpenST_MEX.h"
#include "openst.h"

#include <Windows.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    OPENST_FLOAT *V;
    mwSize V_ndims;
    mwSize *V_dims;
    OPENST_FLOAT *SRC;
    mwSize SRC_ndims;
    mwSize *SRC_dims;
    OPENST_FLOAT *H;
    mwSize H_ndims;
    mwSize *H_dims;
    OPENST_FLOAT EPS, MAX_ITER;
    OPENST_FLOAT *U;
    char *LSM_UNLOCKED;
    int max_iter_int, it, converged;
    size_t NI, NJ, NK;
    OPENST_FLOAT SRCI, SRCJ, SRCK, HI, HJ, HK;
    
    LARGE_INTEGER frequency;
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    OPENST_FLOAT elapsedSeconds;
    QueryPerformanceFrequency(&frequency);
    
    OPENST_MEX_CHECK((nrhs < 5), "Not enough input arguments.");
    OPENST_MEX_CHECK((nrhs > 5), "Too many input arguments.");
    OPENST_MEX_CHECK((nlhs < 4), "Not enough output arguments.");
    OPENST_MEX_CHECK((nlhs > 4), "Too many output arguments.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[0], &V,
            &V_ndims, &V_dims)), NULL);
    OPENST_MEX_CHECK((V_ndims != 3), "V must be 3D.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[1], &SRC, &SRC_ndims,
            &SRC_dims)), NULL);
    OPENST_MEX_CHECK((OPENST_MEX_GetNumel(SRC_ndims, SRC_dims) != 3),
            "numel(SRC) must be 3.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[2], &H, &H_ndims,
            &H_dims)), NULL);
    OPENST_MEX_CHECK((OPENST_MEX_GetNumel(H_ndims, H_dims) != 3),
            "numel(H) must be 3.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleScalar(prhs[3], &EPS)),
            NULL);
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleScalar(prhs[4], &MAX_ITER)),
            NULL);
    OPENST_MEX_CHECK(!OPENST_MEX_ValueIsInteger(MAX_ITER),
            "MAX_ITER must be an integer value.");
    
    NI = V_dims[2];
    NJ = V_dims[1];
    NK = V_dims[0];
    HI = H[0];
    HJ = H[1];
    HK = H[2];
    SRCI = SRC[0];
    SRCJ = SRC[1];
    SRCK = SRC[2];
    
    OPENST_MEX_CHECK(((plhs[0] =
            mxCreateNumericArray(3, V_dims, mxDOUBLE_CLASS,
            mxREAL)) == NULL), NULL);
    OPENST_MEX_CHECK(((U = mxGetPr(plhs[0])) == NULL), NULL);
    
    OPENST_MEX_CHECK(((LSM_UNLOCKED = mxMalloc(NI * NJ * NK *
            sizeof(OPENST_FLOAT))) == NULL), NULL);
    
    QueryPerformanceCounter(&start);
    OPENST_MEX_CHECK(((OpenST_LSM3D(U, LSM_UNLOCKED, V,
            NI, NJ, NK,
            HI, HJ, HK,
            SRCI, SRCJ, SRCK,
            EPS, (int) MAX_ITER,
            &it, &converged)
            != OPENST_ERR_SUCCESS)), "OpenST_LSM3D Error");
    QueryPerformanceCounter(&end);
    elapsedSeconds = (end.QuadPart - start.QuadPart) / (OPENST_FLOAT)frequency.QuadPart;
    
    OPENST_MEX_CHECK(((plhs[1] =
            mxCreateDoubleScalar((OPENST_FLOAT) converged)) == NULL), NULL);
    
    OPENST_MEX_CHECK(((plhs[2] =
            mxCreateDoubleScalar((OPENST_FLOAT) it)) == NULL), NULL);
    
    OPENST_MEX_CHECK(((plhs[3] =
            mxCreateDoubleScalar((OPENST_FLOAT) elapsedSeconds)) == NULL), NULL);
    
    mxFree(LSM_UNLOCKED);
    
}
