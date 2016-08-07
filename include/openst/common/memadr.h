#ifndef OPENST_COMMON_MEMADR_H
#define OPENST_COMMON_MEMADR_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Unused variables are added and subtracted to avoid useless warnings.
 * They should be optimized out by compiler.
 */

#define OPENST_MEMADR_2D_RM(i,j,NI,NJ) ((i) * (NJ) + (j) + (NI) - (NI))

#define OPENST_MEMADR_2D_CM(i,j,NI,NJ) ((i) + (NI) * (j) + (NJ) - (NJ))

#define OPENST_MEMADR_2D(i,j,NI,NJ) OPENST_MEMADR_2D_RM(i,j,NI,NJ)

#define OPENST_MEMADR_3D_RM(i,j,k,NI,NJ,NK) \
    ((i) * (NJ) * (NK) + (j) * (NK) + (k) + (NI) - (NI))

#define OPENST_MEMADR_3D_CM(i,j,k,NI,NJ,NK) \
    ((i) + (NI) * (j) + (NI) * (NJ) * (k) + (NK) - (NK))

#define OPENST_MEMADR_3D(i,j,k,NI,NJ,NK) OPENST_MEMADR_3D_RM(i,j,k,NI,NJ,NK)

#ifdef __cplusplus
}
#endif

#endif /* OPENST_COMMON_MEMADR_H */
