#include "openst/common/float.h"


int OpenST_FLOAT_ApproximatelyEqual(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon)
{
    return OPENST_FLOAT_FABS(a - b) <= ( (OPENST_FLOAT_FABS(a) < OPENST_FLOAT_FABS(b) ? OPENST_FLOAT_FABS(b) : OPENST_FLOAT_FABS(a)) * epsilon);
}


int OpenST_FLOAT_EssentiallyEqual(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon)
{
    return OPENST_FLOAT_FABS(a - b) <= ( (OPENST_FLOAT_FABS(a) > OPENST_FLOAT_FABS(b) ? OPENST_FLOAT_FABS(b) : OPENST_FLOAT_FABS(a)) * epsilon);
}


int OpenST_FLOAT_DefinitelyGreater(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon)
{
    return (a - b) > ( (OPENST_FLOAT_FABS(a) < OPENST_FLOAT_FABS(b) ? OPENST_FLOAT_FABS(b) : OPENST_FLOAT_FABS(a)) * epsilon);
}


int OpenST_FLOAT_DefinitelyLess(OPENST_FLOAT a, OPENST_FLOAT b, OPENST_FLOAT epsilon)
{
    return (b - a) > ( (OPENST_FLOAT_FABS(a) < OPENST_FLOAT_FABS(b) ? OPENST_FLOAT_FABS(b) : OPENST_FLOAT_FABS(a)) * epsilon);
}


int OpenST_FLOAT_GetNeighboorSizeT(OPENST_FLOAT real, size_t *left, size_t *right){
    int is_inbetween;
    OPENST_FLOAT ifloor, iceil;

    ifloor = OPENST_FLOAT_FLOOR(real);
    iceil = OPENST_FLOAT_CEIL(real);

    if (OpenST_FLOAT_EssentiallyEqual(real, ifloor, DBL_EPSILON)) {
        *left = (size_t) ifloor;
        *right = *left;
        is_inbetween = 0;
    } else if (OpenST_FLOAT_EssentiallyEqual(real, iceil, DBL_EPSILON)) {
        *left = (size_t) iceil;
        *right = *left;
        is_inbetween = 0;
    } else {
        *left = (size_t) ifloor;
        *right = (size_t) iceil;
        is_inbetween = 1;
    }

    return is_inbetween;
}
