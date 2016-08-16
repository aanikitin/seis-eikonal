#ifndef OPENST_COMMON_INTEROP_H
#define OPENST_COMMON_INTEROP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

#include "openst/common/macros.h"

OPENST_API void OpenST_INTEROP_FreePointer(void *ptr);

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_INTEROP_H */
