#ifndef OPENST_COMMON_VERSION_H
#define OPENST_COMMON_VERSION_H

#ifdef __cplusplus
extern "C" {
#endif

#include "openst/common/macros.h"

#define OPENST_VERSION_STABLE 0
#define OPENST_VERSION_MAJOR 0
#define OPENST_VERSION_MINOR 2
#define OPENST_VERSION_PATCH 0

OPENST_API extern const int OPENST_VERSION_STABLE_STATIC;
OPENST_API extern const int OPENST_VERSION_MAJOR_STATIC;
OPENST_API extern const int OPENST_VERSION_MINOR_STATIC;
OPENST_API extern const int OPENST_VERSION_PATCH_STATIC;

OPENST_API extern const char OPENST_VERSION_STR_STATIC[];
OPENST_API extern const char OPENST_VERSION_STR_FULL_STATIC[];

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_VERSION_H */
