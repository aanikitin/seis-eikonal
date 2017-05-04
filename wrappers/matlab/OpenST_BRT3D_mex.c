/* Copyright 2014-2017 Alexandr Nikitin. */

#include "mex.h"
#include "matrix.h"
#include "OpenST_MEX.h"
#include "openst.h"

#include <stdlib.h>
#include <Windows.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    OPENST_FLOAT *T;
    mwSize T_ndims;
    mwSize *T_dims;
    OPENST_FLOAT *V;
    mwSize V_ndims;
    mwSize *V_dims;
    OPENST_FLOAT *H;
    mwSize H_ndims;
    mwSize *H_dims;
    OPENST_FLOAT TSTEP;
    OPENST_FLOAT *RCV;
    mwSize RCV_ndims;
    mwSize *RCV_dims;
    OPENST_FLOAT *SRC;
    mwSize SRC_ndims;
    mwSize *SRC_dims;
    size_t NI, NJ, NK;
    OPENST_FLOAT HI, HJ, HK;
    OPENST_FLOAT RCVI, RCVJ, RCVK, SRCI, SRCJ, SRCK;
    OPENST_FLOAT MAX_SEG;
    OPENST_FLOAT *RAY, *RAYout;
    size_t RAY_NI, RAY_NJ;
    mwSize RAY_dims[2];
    
    LARGE_INTEGER frequency;
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    OPENST_FLOAT elapsedSeconds;
    QueryPerformanceFrequency(&frequency);
    
    OPENST_MEX_CHECK((nrhs < 7), "Not enough input arguments.");
    OPENST_MEX_CHECK((nrhs > 7), "Too many input arguments.");
    OPENST_MEX_CHECK((nlhs < 2), "Not enough output arguments.");
    OPENST_MEX_CHECK((nlhs > 2), "Too many output arguments.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[0], &T,
            &T_ndims, &T_dims)), NULL);
    OPENST_MEX_CHECK((T_ndims != 3), "T must be 3D.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[1], &V,
            &V_ndims, &V_dims)), NULL);
    OPENST_MEX_CHECK((V_ndims != 3), "V must be 3D.");
    
    OPENST_MEX_CHECK(( (T_dims[0] != V_dims[0]) ||
            (T_dims[1] != V_dims[1]) ||
            (T_dims[2] != V_dims[2]) ),
            "T and V dimensions must agree.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[2], &H, &H_ndims,
            &H_dims)), NULL);
    OPENST_MEX_CHECK((OPENST_MEX_GetNumel(H_ndims, H_dims) != 3),
            "numel(H) must be 3.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleScalar(prhs[3], &TSTEP)),
            NULL);
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[4], &RCV, &RCV_ndims,
            &RCV_dims)), NULL);
    OPENST_MEX_CHECK((OPENST_MEX_GetNumel(RCV_ndims, RCV_dims) != 3),
            "numel(RCV) must be 3.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleArray(prhs[5], &SRC, &SRC_ndims,
            &SRC_dims)), NULL);
    OPENST_MEX_CHECK((OPENST_MEX_GetNumel(SRC_ndims, SRC_dims) != 3),
            "numel(SRC) must be 3.");
    
    OPENST_MEX_CHECK((OPENST_MEX_GetDoubleScalar(prhs[6], &MAX_SEG)),
            NULL);
    OPENST_MEX_CHECK((!OPENST_MEX_ValueIsInteger(MAX_SEG)),
            "MAX_SEG must be an integer value.");
    
    NI = V_dims[2];
    NJ = V_dims[1];
    NK = V_dims[0];
    HI = H[0];
    HJ = H[1];
    HK = H[2];
    SRCI = SRC[0];
    SRCJ = SRC[1];
    SRCK = SRC[2];
    RCVI = RCV[0];
    RCVJ = RCV[1];
    RCVK = RCV[2];
    
    QueryPerformanceCounter(&start);
    OPENST_MEX_CHECK((OpenST_BRT3D_Trace(T, V,
            NI, NJ, NK, HI, HJ, HK, TSTEP,
            RCVI, RCVJ, RCVK, SRCI, SRCJ, SRCK, (size_t) MAX_SEG,
            &RAY, &RAY_NI, &RAY_NJ)), "OpenST_BRT3D Error");
    QueryPerformanceCounter(&end);
    elapsedSeconds = (end.QuadPart - start.QuadPart) / (OPENST_FLOAT)frequency.QuadPart;
    
    RAY_dims[0] = RAY_NJ;
    RAY_dims[1] = RAY_NI;
    
    OPENST_MEX_CHECK(((plhs[0] =
            mxCreateNumericArray(2, RAY_dims, mxDOUBLE_CLASS,
            mxREAL)) == NULL), NULL);
    OPENST_MEX_CHECK(((RAYout = mxGetPr(plhs[0])) == NULL), NULL);
    
    OPENST_MEX_CHECK(((plhs[1] =
            mxCreateDoubleScalar((OPENST_FLOAT) elapsedSeconds)) == NULL), NULL);
    
    memcpy(RAYout, RAY, RAY_NI * RAY_NJ * sizeof(OPENST_FLOAT));
    
}
