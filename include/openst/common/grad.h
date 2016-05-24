#ifndef OPENST_COMMON_GRAD_H
#define OPENST_COMMON_GRAD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

#include "openst/common/memadr.h"
#include "openst/common/macros.h"

OPENST_API void OpenST_GRAD_Grad3D(double *A, size_t NI, size_t NJ, size_t NK,
            double HI, double HJ, double HK,
            size_t i, size_t j, size_t k,
            double *gradi, double *gradj, double *gradk);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_GRAD_H */
