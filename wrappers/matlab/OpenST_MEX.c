#include "OpenST_MEX.h"


void OPENST_MEX_PrintError(OPENST_MEX_ERRCODE errcode, char* errid,
        char *file, int line, char *msg){
    if(msg == NULL){
        mexErrMsgIdAndTxt(errid,"MEX Error in [%s,%i]: code %i.", \
                file,line,errcode,msg);
    } else {
        mexErrMsgIdAndTxt(errid,"MEX Error in [%s,%i]: code %i: %s", \
                file,line,errcode,msg);
    }
}


OPENST_MEX_ERRCODE OPENST_MEX_GetDoubleArray(mxArray *INP,
        double **ptr, mwSize *ndims,
        mwSize **dims){
    
    if( !mxIsDouble(INP) ){
        return OPENST_MEX_ERRCODE_WRONG_TYPE;
    }
    
    *ndims = mxGetNumberOfDimensions(INP);
    
    *dims = mxGetDimensions(INP);
    
    *ptr = mxGetPr(INP);
    
    if( *ptr == NULL) {
        return OPENST_MEX_ERRCODE_NOT_REAL;
    }
    
    return OPENST_MEX_ERRCODE_SUCCESS;
}


mwSize OPENST_MEX_GetNumel(mwSize ndims, mwSize *dims){
    mwSize numel, i;
    numel = 1;
    for(i = 0; i < ndims; ++i){
        numel *= dims[i];
    }
    return numel;
}


int OPENST_MEX_ValueIsInteger(double value){
    return (ceil(value) == floor(value));
}


OPENST_MEX_ERRCODE OPENST_MEX_GetDoubleScalar(mxArray *INP, double *out){
    
    if( !mxIsDouble(INP) ){
        return OPENST_MEX_ERRCODE_WRONG_TYPE;
    }
    if( mxIsComplex(INP) ){
        return OPENST_MEX_ERRCODE_NOT_REAL;
    }
    if( mxGetNumberOfElements(INP) != 1 ){
        return OPENST_MEX_ERRCODE_WRONG_DIMS;
    }
    *out = mxGetScalar(INP);
    return OPENST_MEX_ERRCODE_SUCCESS;
}
