#include "openst/common/float.h"


int OpenST_FLOAT_ApproximatelyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


int OpenST_FLOAT_EssentiallyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


int OpenST_FLOAT_DefinitelyGreater(double a, double b, double epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


int OpenST_FLOAT_DefinitelyLess(double a, double b, double epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


int OpenST_FLOAT_GetNeighboorSizeT(double real, size_t *left, size_t *right){
    int is_inbetween;
    double ifloor, iceil;

    ifloor = floor(real);
    iceil = ceil(real);

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
