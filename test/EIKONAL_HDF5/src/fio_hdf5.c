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
                           size_t *ndims, size_t **dims){
    
    hid_t dataset;
    hid_t filespace;
    herr_t status;
    int i;
    int lndims;
    hsize_t *ldims;
    double *larr;
    size_t total_size;
    
    M_H5ERRCHK(dataset = H5Dopen(file, dataset_name, H5P_DEFAULT));
    M_H5ERRCHK(filespace = H5Dget_space(dataset));
    M_H5ERRCHK(lndims = H5Sget_simple_extent_ndims(filespace));

    M_MEMCHK(ldims = (hsize_t *)malloc(sizeof(hsize_t) * (hsize_t)lndims));

    M_H5ERRCHK(status = H5Sget_simple_extent_dims(filespace, ldims, NULL));

    total_size = sizeof(double);
    for(i = 0; i < lndims; ++i){
        total_size *= ldims[i];
    }

    M_MEMCHK(larr = (double *)malloc(total_size));

    M_H5ERRCHK(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, larr));

    M_H5ERRCHK(H5Dclose(dataset));
    M_H5ERRCHK(H5Sclose(filespace));

    M_MEMCHK(*dims = (size_t *)malloc(sizeof(size_t) * (size_t)lndims));

    for(i = 0; i < lndims; ++i){
        (*dims)[i] = ldims[i];
    }
    free(ldims);

    *arr = larr;
    *ndims = (size_t)lndims;
    for(i = 0; i < lndims; ++i){
        total_size *= ldims[i];
    }

    return 0;
}


int hdf5_read_array_double_3d(hid_t file, char *dataset_name, double **arr,
                              size_t *NI, size_t *NJ, size_t *NK){
    size_t ndims;
    size_t *dims;
    int errcode = 0;

    errcode = hdf5_read_array_double(file, dataset_name, arr, &ndims, &dims);

    *NI = dims[0];
    *NJ = dims[1];
    *NK = dims[2];

    return errcode;
}


int hdf5_read_model_double(char *filename, char *dataset,
                    double **V,
                    double *HI, double *HJ, double *HK,
                    size_t *NI, size_t *NJ, size_t *NK){
    int errcode = 0;
    hid_t file;
    
    M_H5ERRCHK(file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT));
    
    errcode = hdf5_read_array_double_3d(file, dataset, V, NI, NJ, NK);
    
    M_H5ERRCHK(H5LTread_dataset_double(file, "HI", HI));
    M_H5ERRCHK(H5LTread_dataset_double(file, "HJ", HJ));
    M_H5ERRCHK(H5LTread_dataset_double(file, "HK", HK));
    
    M_H5ERRCHK(H5Fclose(file));
    return errcode;
}


int hdf5_write_time_double(char *filename, char *dataset,
                    double *U,
                    double *HI, double *HJ, double *HK,
                    size_t NI, size_t NJ, size_t NK){
    int ndims = 3;
    hsize_t dims[3];
    int errcode = 0;
    hid_t file;
    
    M_H5ERRCHK(file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    
    dims[0] = NI;
    dims[1] = NJ;
    dims[2] = NK;
    
    M_H5ERRCHK(H5LTmake_dataset_double(file, dataset, ndims, dims, U));
    
    ndims = 1;
    dims[0] = 1;
    M_H5ERRCHK(H5LTmake_dataset_double(file, "HI", ndims, dims, HI));
    M_H5ERRCHK(H5LTmake_dataset_double(file, "HJ", ndims, dims, HJ));
    M_H5ERRCHK(H5LTmake_dataset_double(file, "HK", ndims, dims, HK));
    
    M_H5ERRCHK(H5Fclose(file));
    return errcode;
}


int hdf5_read_array_float(hid_t file, char *dataset_name, float **arr,
                           size_t *ndims, size_t **dims){

    hid_t dataset;
    hid_t filespace;
    herr_t status;
    int i;
    int lndims;
    hsize_t *ldims;
    float *larr;
    size_t total_size;

    M_H5ERRCHK(dataset = H5Dopen(file, dataset_name, H5P_DEFAULT));
    M_H5ERRCHK(filespace = H5Dget_space(dataset));
    M_H5ERRCHK(lndims = H5Sget_simple_extent_ndims(filespace));

    M_MEMCHK(ldims = (hsize_t *)malloc(sizeof(hsize_t) * (hsize_t)lndims));

    M_H5ERRCHK(status = H5Sget_simple_extent_dims(filespace, ldims, NULL));

    total_size = sizeof(float);
    for(i = 0; i < lndims; ++i){
        total_size *= ldims[i];
    }

    M_MEMCHK(larr = (float *)malloc(total_size));

    M_H5ERRCHK(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, larr));

    M_H5ERRCHK(H5Dclose(dataset));
    M_H5ERRCHK(H5Sclose(filespace));

    M_MEMCHK(*dims = (size_t *)malloc(sizeof(size_t) * (size_t)lndims));

    for(i = 0; i < lndims; ++i){
        (*dims)[i] = ldims[i];
    }
    free(ldims);

    *arr = larr;
    *ndims = (size_t)lndims;
    for(i = 0; i < lndims; ++i){
        total_size *= ldims[i];
    }

    return 0;
}


int hdf5_read_array_float_3d(hid_t file, char *dataset_name, float **arr,
                              size_t *NI, size_t *NJ, size_t *NK){
    size_t ndims;
    size_t *dims;
    int errcode = 0;

    errcode = hdf5_read_array_float(file, dataset_name, arr, &ndims, &dims);

    *NI = dims[0];
    *NJ = dims[1];
    *NK = dims[2];

    return errcode;
}


int hdf5_read_model_float(char *filename, char *dataset,
                    float **V,
                    float *HI, float *HJ, float *HK,
                    size_t *NI, size_t *NJ, size_t *NK){
    int errcode = 0;
    hid_t file;

    M_H5ERRCHK(file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT));

    errcode = hdf5_read_array_float_3d(file, dataset, V, NI, NJ, NK);

    M_H5ERRCHK(H5LTread_dataset_float(file, "HI", HI));
    M_H5ERRCHK(H5LTread_dataset_float(file, "HJ", HJ));
    M_H5ERRCHK(H5LTread_dataset_float(file, "HK", HK));

    M_H5ERRCHK(H5Fclose(file));
    return errcode;
}


int hdf5_write_time_float(char *filename, char *dataset,
                    float *U,
                    float *HI, float *HJ, float *HK,
                    size_t NI, size_t NJ, size_t NK){
    int ndims = 3;
    hsize_t dims[3];
    int errcode = 0;
    hid_t file;

    M_H5ERRCHK(file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));

    dims[0] = NI;
    dims[1] = NJ;
    dims[2] = NK;

    M_H5ERRCHK(H5LTmake_dataset_float(file, dataset, ndims, dims, U));

    ndims = 1;
    dims[0] = 1;
    M_H5ERRCHK(H5LTmake_dataset_float(file, "HI", ndims, dims, HI));
    M_H5ERRCHK(H5LTmake_dataset_float(file, "HJ", ndims, dims, HJ));
    M_H5ERRCHK(H5LTmake_dataset_float(file, "HK", ndims, dims, HK));

    M_H5ERRCHK(H5Fclose(file));
    return errcode;
}
