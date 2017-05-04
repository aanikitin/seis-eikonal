#ifndef OPENST_RAYTRACE_BACKTRACE_H
#define OPENST_RAYTRACE_BACKTRACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <math.h>

#include "openst/common/float.h"
#include "openst/common/hacks.h"
#include "openst/common/memadr.h"
#include "openst/common/dynarr.h"
#include "openst/common/coordsys.h"
#include "openst/common/grad.h"
#include "openst/common/macros.h"
#include "openst/common/error.h"
#include "openst/common/interp.h"

OPENST_API OPENST_FLOAT OpenST_BRT3D_SuggestTSTEP(OPENST_FLOAT vmax,
                               OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK);

OPENST_API OPENST_ERR OpenST_BRT3D_Trace(OPENST_FLOAT *T, OPENST_FLOAT *V,
                size_t NI, size_t NJ, size_t NK,
                OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK, OPENST_FLOAT TSTEP,
                OPENST_FLOAT RCVI, OPENST_FLOAT RCVJ, OPENST_FLOAT RCVK,
                OPENST_FLOAT SRCI, OPENST_FLOAT SRCJ, OPENST_FLOAT SRCK, size_t MAX_SEG,
                OPENST_FLOAT **RAY, size_t *RAY_NI, size_t *RAY_NJ);

OPENST_API OPENST_ERR OpenST_BRT3D_Step(OPENST_FLOAT *T, OPENST_FLOAT *V,
                     size_t NI, size_t NJ, size_t NK,
                     OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
                     OPENST_FLOAT TSTEP,
                     OPENST_FLOAT CURI, OPENST_FLOAT CURJ, OPENST_FLOAT CURK,
                     OPENST_FLOAT *DSTI, OPENST_FLOAT *DSTJ, OPENST_FLOAT *DSTK);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_RAYTRACE_BACKTRACE_H */
