#include "hdf5.h"

#include "hdf5ReaderWriter.h"

void createAndSetUpH5Env(hid_t *file_id, hid_t *dataset_id, hsize_t chunk_dims[], const char* restrict fileNameHDF5){
    *file_id = H5Fcreate(fileNameHDF5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Initialize dataspace with 0 rows
    const hsize_t n_cols = chunk_dims[1];
    hsize_t dims[2] = {0, n_cols};
    hsize_t max_dims[2] = {H5S_UNLIMITED, n_cols};
    hid_t dataspace_id = H5Screate_simple(2, dims, max_dims);

    //Create dataset creation property list
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    //hsize_t chunk_dims[2] = {n_rows_slab, n_cols};
    H5Pset_chunk(plist, 2, chunk_dims);

    //Create dataset
    *dataset_id = H5Dcreate2(*file_id, "/dset", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, plist, H5P_DEFAULT);

    H5Pclose(plist);
    H5Sclose(dataspace_id);
}

void extendDatasetH5(hid_t *dataset_id, hsize_t chunk_dims[], int t){
    hsize_t dims[2];
    dims[0]= t*chunk_dims[0];
    dims[1]=chunk_dims[1];
    H5Dset_extent(*dataset_id, dims);
}

void appendBufferToDatasetH5(double* buffer, hid_t *dataset_id, hid_t mem_space, hsize_t start[], hsize_t count[], int t, int n_particles){
    hid_t dataspace_id = H5Dget_space(*dataset_id);
    start[0] = (t-1)*n_particles;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(*dataset_id, H5T_NATIVE_DOUBLE, mem_space, dataspace_id, H5P_DEFAULT, buffer);
    H5Sclose(dataspace_id);
}

void closeFilesH5(hid_t *mem_space, hid_t *dataset_id, hid_t *file_id){
    H5Sclose(*mem_space);
    H5Dclose(*dataset_id);
    H5Fclose(*file_id);
}
