#include "openst/common/version.h"

#define M_OPENST_VERSION_RELEASE_STR OPENST_M_XSTR(OPENST_VERSION_STABLE)
#define M_OPENST_VERSION_MAJOR_STR OPENST_M_XSTR(OPENST_VERSION_MAJOR)
#define M_OPENST_VERSION_MINOR_STR OPENST_M_XSTR(OPENST_VERSION_MINOR)
#define M_OPENST_VERSION_PATCH_STR OPENST_M_XSTR(OPENST_VERSION_PATCH)

#if OPENST_VERSION_STABLE
    #define M_OPENST_VERSION_RELEASE_STR_FULL "Stable Version"
#else
    #define M_OPENST_VERSION_RELEASE_STR_FULL "Experimental Version"
#endif

OPENST_API const int OPENST_VERSION_STABLE_STATIC = OPENST_VERSION_STABLE;
OPENST_API const int OPENST_VERSION_MAJOR_STATIC = OPENST_VERSION_MAJOR;
OPENST_API const int OPENST_VERSION_MINOR_STATIC = OPENST_VERSION_MINOR;
OPENST_API const int OPENST_VERSION_PATCH_STATIC = OPENST_VERSION_PATCH;

#define M_OPENST_VERSION_STR M_OPENST_VERSION_MAJOR_STR"."\
    M_OPENST_VERSION_MINOR_STR"."\
    M_OPENST_VERSION_PATCH_STR

#define M_OPENST_VERSION_STR_FULL "OpenST "M_OPENST_VERSION_RELEASE_STR_FULL" "\
    M_OPENST_VERSION_STR


OPENST_API const char OPENST_VERSION_STR_STATIC[] = M_OPENST_VERSION_STR;
OPENST_API const char OPENST_VERSION_STR_FULL_STATIC[] = M_OPENST_VERSION_STR_FULL;
