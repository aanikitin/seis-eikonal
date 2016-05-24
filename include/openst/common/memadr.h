#ifndef OPENST_COMMON_MEMADR_H
#define OPENST_COMMON_MEMADR_H

#ifdef __cplusplus
extern "C" {
#endif

#define OPENST_MEMADR_2D_RM(i,j,NI,NJ) ((i) * (NJ) + (j))

#define OPENST_MEMADR_2D_CM(i,j,NI,NJ) ((i) + (NI) * (j))

#define OPENST_MEMADR_2D(i,j,NI,NJ) OPENST_MEMADR_2D_RM(i,j,NI,NJ)

#define OPENST_MEMADR_3D_RM(i,j,k,NI,NJ,NK) \
    ((i) * (NJ) * (NK) + (j) * (NK) + (k))

#define OPENST_MEMADR_3D_CM(i,j,k,NI,NJ,NK) \
    ((i) + (NI) * (j) + (NI) * (NJ) * (k))

#define OPENST_MEMADR_3D(i,j,k,NI,NJ,NK) OPENST_MEMADR_3D_RM(i,j,k,NI,NJ,NK)

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_MEMADR_H */
