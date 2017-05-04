#ifndef OPENST_COMMON_GRAD_H
#define OPENST_COMMON_GRAD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

#include "openst/common/float.h"
#include "openst/common/memadr.h"
#include "openst/common/macros.h"

OPENST_API void OpenST_GRAD_Grad_Kernel(OPENST_FLOAT *left, OPENST_FLOAT *center,
                                        OPENST_FLOAT *right, OPENST_FLOAT H, OPENST_FLOAT *grad);

OPENST_API void OpenST_GRAD_Grad3D(OPENST_FLOAT *A, size_t NI, size_t NJ, size_t NK,
            OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK,
            size_t i, size_t j, size_t k,
            OPENST_FLOAT *gradi, OPENST_FLOAT *gradj, OPENST_FLOAT *gradk);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_GRAD_H */
