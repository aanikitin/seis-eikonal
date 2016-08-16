#include "openst/common/interop.h"


void OpenST_INTEROP_FreePointer(void *ptr){
    free(ptr);
}
