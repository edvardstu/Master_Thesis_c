#ifndef HDF5READERWRITER_H
#define HDF5READERWRITER_H

#include "hdf5.h"

void createAndSetUpH5Env(hid_t *file_id, hid_t *dataset_id, hsize_t chunk_dims[], const char* restrict fileNameHDF5);

void extendDatasetH5(hid_t *dataset_id, hsize_t chunk_dims[], int t);

void appendBufferToDatasetH5(double* buffer, hid_t *dataset_id, hid_t mem_space, hsize_t start[], hsize_t count[], int t, int n_particles);

void closeFilesH5(hid_t *mem_space, hid_t *dataset_id, hid_t *file_id);

#endif
