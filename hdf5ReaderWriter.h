#ifndef HDF5READERWRITER_H
#define HDF5READERWRITER_H

#include "hdf5.h"

void createAndSetUpH5Env(hid_t *file_id, hid_t *dataset_id, hsize_t chunk_dims[], const char* restrict fileNameHDF5);

void extendDatasetH5(hid_t *dataset_id, hsize_t chunk_dims[], int writeNumber);

hsize_t n_cols;

void appendBufferToDatasetH5(double buffer[][n_cols], hid_t *dataset_id, hid_t mem_space, hsize_t start[], hsize_t count[], int writeNumber, int n_particles);

void closeFilesH5(hid_t *mem_space, hid_t *dataset_id, hid_t *file_id);

void transposeSimulationData(double buffer[][n_cols], int n_particles, double time, double* x, double* y, double* theta, double* vx, double* vy);

#endif
