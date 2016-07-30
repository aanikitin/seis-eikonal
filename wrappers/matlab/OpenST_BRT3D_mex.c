#include "mex.h"
#include "matrix.h"
#include "OpenST_MEX.h"
#include "openst.h"
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *T;
    mwSize T_ndims;
    mwSize *T_dims;
    double *V;
    mwSize V_ndims;
    mwSize *V_dims;
    double *H;
    mwSize H_ndims;
    mwSize *H_dims;
    double TSTEP;
    double *RCV;
    mwSize RCV_ndims;
    mwSize *RCV_dims;
    double *SRC;
    mwSize SRC_ndims;
    mwSize *SRC_dims;
    size_t NI, NJ, NK;
    double HI, HJ, HK;
    double RCVI, RCVJ, RCVK, SRCI, SRCJ, SRCK;
    double *RAY, *RAYout;
    size_t RAY_NI, RAY_NJ;
    mwSize RAY_dims[2];
    
    OPENST_MEX_CHECK((nrhs < 6), "Not enough input arguments.");
    OPENST_MEX_CHECK((nrhs > 6), "Too many input arguments.");
    OPENST_MEX_CHECK((nlhs < 1), "Not enough output arguments.");
    OPENST_MEX_CHECK((nlhs > 1), "Too many output arguments.");
    
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
    
    OPENST_MEX_CHECK((OpenST_BRT3D_Trace(T, V,
            NI, NJ, NK, HI, HJ, HK, TSTEP,
            RCVI, RCVJ, RCVK, SRCI, SRCJ, SRCK,
            &RAY, &RAY_NI, &RAY_NJ)), "OpenST_BRT3D Error");
    
    RAY_dims[0] = RAY_NJ;
    RAY_dims[1] = RAY_NI;
    
    OPENST_MEX_CHECK(((plhs[0] =
            mxCreateNumericArray(2, RAY_dims, mxDOUBLE_CLASS,
            mxREAL)) == NULL), NULL);
    OPENST_MEX_CHECK(((RAYout = mxGetPr(plhs[0])) == NULL), NULL);
    
    memcpy(RAYout, RAY, RAY_NI * RAY_NJ * sizeof(double));
    
}
