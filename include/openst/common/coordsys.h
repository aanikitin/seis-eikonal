#ifndef OPENST_COMMON_COORDSYS_H
#define OPENST_COMMON_COORDSYS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#include "openst/common/macros.h"
#include "openst/common/error.h"
#include "openst/common/float.h"

OPENST_API double OpenST_CRS_Distance3D(double i0, double j0, double k0,
                                        double i1, double j1, double k1);

OPENST_API OPENST_ERR OpenST_CRS_Cart2IndFloor(double coord,
                                               double cell_size, size_t *ind);

OPENST_API OPENST_ERR OpenST_CRS_Cart2IndRound(double coord,
                                               double cell_size, size_t *ind);

OPENST_API OPENST_ERR OpenST_CRS_Cart2Ind(double coord,
                                          double cell_size, size_t *ind);

OPENST_API int OpenST_CRS_IsPointWithinBounds(double PI, double PJ, double PK,
                                       size_t NI, size_t NJ, size_t NK,
                                       double HI, double HJ, double HK);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_COORDSYS_H */
