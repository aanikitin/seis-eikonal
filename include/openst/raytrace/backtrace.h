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

#define OPENST_RAYTRACE_SUCCESS 0
#define OPENST_RAYTRACE_ERR_MEM -1
#define OPENST_RAYTRACE_ERR_COORD -2
#define OPENST_RAYTRACE_ERR_BORDER_PASSED -3
#define OPENST_RAYTRACE_ERR_GRAD_ZERO -4

OPENST_API double OpenST_BRT3D_SuggestTSTEP(double vmax,
                               double HI, double HJ, double HK);

OPENST_API int OpenST_BRT3D_Step(double *T, double *V,
                     size_t NI, size_t NJ, size_t NK,
                     double HI, double HJ, double HK,
                     double TSTEP,
                     double CURI, double CURJ, double CURK,
                     size_t ind_cur_i, size_t ind_cur_j, size_t ind_cur_k,
                     double *DSTI, double *DSTJ, double *DSTK);

OPENST_API int OpenST_BRT3D_Trace(double *T, double *V,
                size_t NI, size_t NJ, size_t NK,
                double HI, double HJ, double HK, double TSTEP,
                double RCVI, double RCVJ, double RCVK,
                double SRCI, double SRCJ, double SRCK,
                double **RAY, size_t *RAY_NI, size_t *RAY_NJ);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_RAYTRACE_BACKTRACE_H */
