#include "openst/common/arrayops.h"


void OpenST_AOP_GetArrStats(double *A, size_t numel, double *min,
                        double *max, double *mean){
    size_t i;

    double tmin = INFINITY;
    double tmax = -INFINITY;
    double tmean = 0.0;
    double cur;

    if(min == NULL && max == NULL && mean == NULL){
        return;
    }

    for(i = 0; i < numel; ++i){
        cur = A[i];
        if(min != NULL){
            tmin = fmin(tmin,cur);
        }
        if(max != NULL){
            tmax = fmax(tmax,cur);
        }
        if(mean != NULL){
            tmean += cur;
        }
    }

    if(mean != NULL){
        tmean /= numel;
    }

    *min = tmin;
    *max = tmax;
    *mean = tmean;

}
