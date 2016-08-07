#include "openst/common/coordsys.h"


double OpenST_CRS_Distance3D(double i0, double j0, double k0,
                             double i1, double j1, double k1){
    return sqrt(pow(i1 - i0, 2.0) + pow(j1 - j0, 2.0) + pow(k1 - k0, 2.0));
}


OPENST_ERR OpenST_CRS_Cart2IndFloor(double coord, double cell_size, size_t *ind)
{
    double d_ind;
    if(coord < 0) return OPENST_ERR_TYPE_CONVERSION;
    d_ind = coord/cell_size;
    if(d_ind > SIZE_MAX) return OPENST_ERR_TYPE_CONVERSION;
    *ind = (size_t)(d_ind);
    return OPENST_ERR_SUCCESS;
}


OPENST_ERR OpenST_CRS_Cart2IndRound(double coord, double cell_size, size_t *ind)
{
    double d_ind;
    if(coord < 0) return OPENST_ERR_TYPE_CONVERSION;
    d_ind = coord/cell_size;
    if(d_ind > SIZE_MAX) return OPENST_ERR_TYPE_CONVERSION;
    *ind = (size_t)(d_ind + 0.5);
    return OPENST_ERR_SUCCESS;
}


OPENST_ERR OpenST_CRS_Cart2Ind(double coord, double cell_size, size_t *ind)
{
    return OpenST_CRS_Cart2IndRound(coord, cell_size, ind);
}


//TODO: change to stored domain bounds in the future release
int OpenST_CRS_IsPointWithinBounds(double PI, double PJ, double PK,
                                       size_t NI, size_t NJ, size_t NK,
                                       double HI, double HJ, double HK){
    int il, ih, jl, jh, kl, kh;

    il = OpenST_FLOAT_DefinitelyLess(PI, 0, DBL_EPSILON);
    ih = OpenST_FLOAT_DefinitelyGreater(PI, ((double)(NI - 1) * HI), DBL_EPSILON);

    jl = OpenST_FLOAT_DefinitelyLess(PJ, 0, DBL_EPSILON);
    jh = OpenST_FLOAT_DefinitelyGreater(PJ, ((double)(NJ - 1) * HJ), DBL_EPSILON);

    kl = OpenST_FLOAT_DefinitelyLess(PK, 0, DBL_EPSILON);
    kh = OpenST_FLOAT_DefinitelyGreater(PK, ((double)(NK - 1) * HK), DBL_EPSILON);

    return (il || ih || jl || jh || kl || kh);
}
