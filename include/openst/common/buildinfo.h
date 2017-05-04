#ifndef OPENST_COMMON_BUILDINFO_H
#define OPENST_COMMON_BUILDINFO_H

#ifdef __cplusplus
extern "C" {
#endif

#include "openst/common/float.h"
#include "openst/common/macros.h"
#include "openst/common/version.h"

OPENST_API extern const char OPENST_BUILDINFO_C_COMPILER_STATIC[];
OPENST_API extern const char OPENST_BUILDINFO_C_COMPILER_ID_STATIC[];
OPENST_API extern const char OPENST_BUILDINFO_C_COMPILER_VERSION_STATIC[];
OPENST_API extern const char OPENST_BUILDINFO_C_FLAGS_STATIC[];
OPENST_API extern const char OPENST_BUILDINFO_LINK_TYPE_STATIC[];
OPENST_API extern const char OPENST_BUILDINFO_FLOAT_PRECISION_STATIC[];
OPENST_API extern const char OPENST_BUILDINFO_STR_FULL_STATIC[];

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_BUILDINFO_H */
