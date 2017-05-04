#include "openst/common/arrayops.h"


void OpenST_AOP_GetArrStats(OPENST_FLOAT *A, size_t numel, OPENST_FLOAT *min,
                        OPENST_FLOAT *max, OPENST_FLOAT *mean){
    size_t i;

    OPENST_FLOAT tmin = OPENST_FLOAT_INF;
    OPENST_FLOAT tmax = -OPENST_FLOAT_INF;
    OPENST_FLOAT tmean = OPENST_FLOAT_0_0;
    OPENST_FLOAT cur;

    if(min == NULL && max == NULL && mean == NULL){
        return;
    }

    for(i = 0; i < numel; ++i){
        cur = A[i];
        if(min != NULL){
            tmin = OPENST_FLOAT_FMIN(tmin,cur);
        }
        if(max != NULL){
            tmax = OPENST_FLOAT_FMAX(tmax,cur);
        }
        if(mean != NULL){
            tmean += cur;
        }
    }

    if(mean != NULL){
        tmean /= (OPENST_FLOAT)numel;
		*mean = tmean;
    }

	if (min != NULL) {
		*min = tmin;
	}

	if (max != NULL) {
		*max = tmax;
	}

}
