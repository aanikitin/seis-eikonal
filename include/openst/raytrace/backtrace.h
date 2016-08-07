#ifndef OPENST_RAYTRACE_BACKTRACE_H
#define OPENST_RAYTRACE_BACKTRACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <math.h>

#include "openst/common/memadr.h"
#include "openst/common/dynarr.h"
#include "openst/common/coordsys.h"
#include "openst/common/grad.h"
#include "openst/common/macros.h"
#include "openst/common/error.h"
#include "openst/common/interp.h"
#include "openst/common/float.h"

OPENST_API double OpenST_BRT3D_SuggestTSTEP(double vmax,
                               double HI, double HJ, double HK);

OPENST_API OPENST_ERR OpenST_BRT3D_Trace(double *T, double *V,
                size_t NI, size_t NJ, size_t NK,
                double HI, double HJ, double HK, double TSTEP,
                double RCVI, double RCVJ, double RCVK,
                double SRCI, double SRCJ, double SRCK, size_t MAX_SEG,
                double **RAY, size_t *RAY_NI, size_t *RAY_NJ);

OPENST_API OPENST_ERR OpenST_BRT3D_Step(double *T, double *V,
                     size_t NI, size_t NJ, size_t NK,
                     double HI, double HJ, double HK,
                     double TSTEP,
                     double CURI, double CURJ, double CURK,
                     double *DSTI, double *DSTJ, double *DSTK);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_RAYTRACE_BACKTRACE_H */
