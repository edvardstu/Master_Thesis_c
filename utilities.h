#ifndef UTILITIES_H
#define UTILITIES_H

#include <gsl/gsl_rng.h>

extern int errno;

void testFunc();

unsigned long int random_seed();

void setUpRNG(const gsl_rng_type **T, gsl_rng **r, bool rndSeed);

double randDouble(double min, double max, gsl_rng ** r);

void openFile(const char * restrict fileName, FILE** fp);

void closeFile(const char * restrict fileName, FILE** fp);

void sunflower(double x[], double y[], int n_particles, int alpha, double R);

double radius(int index, int n_particles, int n_boundary, double R);

void sunflower_fixed_boundary(double x[], double y[], int n_particles, int alpha, double R, int n_fixed, double R_reduced);

void uniformRectangle(double x[], double y[], int n_particles, double L, double H, gsl_rng** r);

void uniformFunnel(double x[], double y[], int n_particles, double L, double H, double H_FUNNEL, gsl_rng** r);

void fixedBoundaryFunnel(double x[], double y[], double theta[], int n_particles, int n_fixed, double L, double H, double H_FUNNEL);

double walltime();

const char* restrict createFileNameBase(const char* restrict fileNameBase, bool overwrite);

const char* restrict createFileName(const char* restrict fileNameBase);

const char* restrict createFileNameH5(const char* restrict fileNameBase);

const char* restrict createFileNameSucs(const char* restrict fileNameBase);

void writeSimulationParameters(const char* restrict fileNameBase, double r, double l, double h, double h_funnel, double r_particle, unsigned int n_particles, unsigned int n_fixed_particles, double u_0, double D_r, unsigned int n_steps, double dt, double gamma_t, double gamma_r, double lambda_har, double kappa_har, double gamma_pp, double lambda_pp, double r_cut_off_force_2, double r_cut_off_torque_2, double sigma_pp);

void writeFinalState(const char* restrict fileNameBase, int n_particles, double time, double x[], double y[], double theta[], double vx[], double vy[], double d_r);

void swapPointers(double **Y_i, double **Y_i_prev);

const char* restrict createFileNamePrevious(const char* restrict fileNameBase, bool overwrite);

void readInitialState(const char* restrict fileName, int n_particles, double *time, double x[], double y[], double theta[], double vx[], double vy[], double *d_r);

#endif
