#ifndef OPENST_COMMON_COORDSYS_H
#define OPENST_COMMON_COORDSYS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#include "openst/common/macros.h"

OPENST_API double OpenST_CRS_Distance3D(double i0, double j0, double k0,
                             double i1, double j1, double k1);

OPENST_API int OpenST_CRS_Cart2IndFloor(double coord, double cell_size, size_t *ind);

OPENST_API int OpenST_CRS_Cart2IndRound(double coord, double cell_size, size_t *ind);

OPENST_API int OpenST_CRS_Cart2Ind(double coord, double cell_size, size_t *ind);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_COORDSYS_H */
