#ifndef SYSTEMS_H
#define SYSTEMS_H

#include <gsl/gsl_rng.h>

//void corePeriodicTube(double t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng * r, double N_PARTICLES, double N_FIXED_PARTICLES, double DT, double D_R, double a, double U_0, double L, double H, double LAMBDA_HAR, double KAPPA_HAR);

void corePeriodic(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double gamma_pp, double lambda_pp, double r_cut_off_force_2, double r_cut_off_torque_2,  double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n);

void corePeriodicTube(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double r_cut_off_torque_2, double gamma_pp, double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n);

void corePeriodicFunnel(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double h_funnel, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double r_cut_off_torque_2, double gamma_pp, double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n);

void corePeriodicBM(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double gamma_pp, double lambda_pp, double r_cut_off_force_2_matrix[n_particles][n_particles], double r_cut_off_torque_2,  double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n, int n_A);

void upadteIntegrationParameters(double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int t, double theta, double fx_b, double fy_b, double torque_b, double fx_n, double fy_n, double torque_n, double number_n, double fs_scale, double u_0);

void updateParticleParameters(double* x, double*y, double* theta, double*vx, double* vy, double Y_x, double Y_y, double Y_th, double Y_x_prev, double Y_y_prev, double Y_th_prev, double f_AB1, double f_AB2, double dt, double d_r, double a, gsl_rng **r);

void updateParticleParametersPeriodicX(double* x, double*y, double* theta, double*vx, double* vy, double Y_x, double Y_y, double Y_th, double Y_x_prev, double Y_y_prev, double Y_th_prev, double f_AB1, double f_AB2, double dt, double d_r, double a, gsl_rng **r, double l);

void updateParticleParametersPeriodicXY(double* x, double*y, double* theta, double*vx, double* vy, double Y_x, double Y_y, double Y_th, double Y_x_prev, double Y_y_prev, double Y_th_prev, double f_AB1, double f_AB2, double dt, double d_r, double a, gsl_rng **r, double l, double h);

#endif
