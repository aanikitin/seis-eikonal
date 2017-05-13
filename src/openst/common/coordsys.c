#include "openst/common/coordsys.h"


OPENST_FLOAT OpenST_CRS_Distance3D(OPENST_FLOAT i0, OPENST_FLOAT j0, OPENST_FLOAT k0,
                             OPENST_FLOAT i1, OPENST_FLOAT j1, OPENST_FLOAT k1){
    return OPENST_FLOAT_SQRT(OPENST_FLOAT_POW(i1 - i0, OPENST_FLOAT_2_0) + OPENST_FLOAT_POW(j1 - j0, OPENST_FLOAT_2_0) + OPENST_FLOAT_POW(k1 - k0, OPENST_FLOAT_2_0));
}


OPENST_ERR OpenST_CRS_Cart2IndFloor(OPENST_FLOAT coord, OPENST_FLOAT cell_size, size_t *ind)
{
    OPENST_FLOAT d_ind;
    if(coord < 0) return OPENST_ERR_TYPE_CONVERSION;
    d_ind = coord/cell_size;
    if(d_ind > SIZE_MAX) return OPENST_ERR_TYPE_CONVERSION;
    *ind = (size_t)(d_ind);
    return OPENST_ERR_SUCCESS;
}


OPENST_ERR OpenST_CRS_Cart2IndRound(OPENST_FLOAT coord, OPENST_FLOAT cell_size, size_t *ind)
{
    OPENST_FLOAT d_ind;
    if(coord < 0) return OPENST_ERR_TYPE_CONVERSION;
    d_ind = coord/cell_size;
    if(d_ind > SIZE_MAX) return OPENST_ERR_TYPE_CONVERSION;
    *ind = (size_t)(d_ind + OPENST_FLOAT_0_5);
    return OPENST_ERR_SUCCESS;
}


OPENST_ERR OpenST_CRS_Cart2Ind(OPENST_FLOAT coord, OPENST_FLOAT cell_size, size_t *ind)
{
    return OpenST_CRS_Cart2IndRound(coord, cell_size, ind);
}


//TODO: change to stored domain bounds in the future release
int OpenST_CRS_IsPointNotWithinBounds_3D(OPENST_FLOAT PI, OPENST_FLOAT PJ, OPENST_FLOAT PK,
                                       size_t NI, size_t NJ, size_t NK,
                                       OPENST_FLOAT HI, OPENST_FLOAT HJ, OPENST_FLOAT HK){
    int il, ih, jl, jh, kl, kh;

    il = OpenST_FLOAT_DefinitelyLess(PI, 0, DBL_EPSILON);
    ih = OpenST_FLOAT_DefinitelyGreater(PI, ((OPENST_FLOAT)(NI - 1) * HI), DBL_EPSILON);

    jl = OpenST_FLOAT_DefinitelyLess(PJ, 0, DBL_EPSILON);
    jh = OpenST_FLOAT_DefinitelyGreater(PJ, ((OPENST_FLOAT)(NJ - 1) * HJ), DBL_EPSILON);

    kl = OpenST_FLOAT_DefinitelyLess(PK, 0, DBL_EPSILON);
    kh = OpenST_FLOAT_DefinitelyGreater(PK, ((OPENST_FLOAT)(NK - 1) * HK), DBL_EPSILON);

    return (il || ih || jl || jh || kl || kh);
}


//TODO: change to stored domain bounds in the future release
int OpenST_CRS_IsPointNotWithinBounds_2D(OPENST_FLOAT PI, OPENST_FLOAT PJ,
                                       size_t NI, size_t NJ,
                                       OPENST_FLOAT HI, OPENST_FLOAT HJ){
    int il, ih, jl, jh;

    il = OpenST_FLOAT_DefinitelyLess(PI, 0, DBL_EPSILON);
    ih = OpenST_FLOAT_DefinitelyGreater(PI, ((OPENST_FLOAT)(NI - 1) * HI), DBL_EPSILON);

    jl = OpenST_FLOAT_DefinitelyLess(PJ, 0, DBL_EPSILON);
    jh = OpenST_FLOAT_DefinitelyGreater(PJ, ((OPENST_FLOAT)(NJ - 1) * HJ), DBL_EPSILON);

    return (il || ih || jl || jh);
}
