#ifndef FIO_HDF5_H
#define FIO_HDF5_H

#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>


#include "fio_hdf5.h"


#define M_H5ERRCHK(H5CALL) do { \
    if ((H5CALL) < 0){ \
    fprintf(stderr, "HDF5 ERROR: %s:%i\n", __FILE__, __LINE__);\
    exit(EXIT_FAILURE);\
    } \
    } while (0)


#define M_MEMCHK(CALL) do { \
    if ((CALL) == NULL){ \
    fprintf(stderr, "MEMERROR: %s:%i\n", __FILE__, __LINE__);\
    exit(EXIT_FAILURE);\
    } \
    } while (0)


int hdf5_read_array_double(hid_t file, char *dataset_name, double **arr,
                           size_t *ndims, size_t **dims);


int hdf5_read_array_double_3d(hid_t file, char *dataset_name, double **arr,
                              size_t *NI, size_t *NJ, size_t *NK);


int hdf5_read_model(char *filename, double **V, double *H, 
			      size_t *NI, size_t *NJ, size_t *NK);

int hdf5_write_time(char *filename, double *U, double *H, 
			      size_t NI, size_t NJ, size_t NK);

#endif
