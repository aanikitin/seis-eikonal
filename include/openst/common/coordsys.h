#ifndef OPENST_COMMON_COORDSYS_H
#define OPENST_COMMON_COORDSYS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#include "openst/common/float.h"
#include "openst/common/macros.h"
#include "openst/common/error.h"

OPENST_API OPENST_FLOAT OpenST_CRS_Distance3D(OPENST_FLOAT i0, OPENST_FLOAT j0, OPENST_FLOAT k0,
                                        OPENST_FLOAT i1, OPENST_FLOAT j1, OPENST_FLOAT k1);

OPENST_API OPENST_ERR OpenST_CRS_Cart2IndFloor(OPENST_FLOAT coord,
                                               OPENST_FLOAT cell_size, size_t *ind);

OPENST_API OPENST_ERR OpenST_CRS_Cart2IndRound(OPENST_FLOAT coord,
                                               OPENST_FLOAT cell_size, size_t *ind);

OPENST_API OPENST_ERR OpenST_CRS_Cart2Ind(OPENST_FLOAT coord,
                                          OPENST_FLOAT cell_size, size_t *ind);

OPENST_API int OpenST_CRS_IsPointNotWithinBounds(OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
                                       size_t NI, size_t NJ, size_t NK,
                                       OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_COORDSYS_H */
