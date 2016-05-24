#include "openst/common/coordsys.h"


int OpenST_CRS_Cart2IndFloor(double coord, double cell_size, size_t *ind)
{
  double d_ind;
  if(coord < 0) return -1;
  d_ind = coord/cell_size;
  if(d_ind > SIZE_MAX) return -2;
  *ind = (size_t)(d_ind);
  return 0;
}


int OpenST_CRS_Cart2IndRound(double coord, double cell_size, size_t *ind)
{
  double d_ind;
  if(coord < 0) return -1;
  d_ind = coord/cell_size;
  if(d_ind > SIZE_MAX) return -2;
  *ind = (size_t)(d_ind + 0.5);
  return 0;
}


int OpenST_CRS_Cart2Ind(double coord, double cell_size, size_t *ind)
{
  return OpenST_CRS_Cart2IndRound(coord, cell_size, ind);
}
