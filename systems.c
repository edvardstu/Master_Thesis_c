#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//include <errno.h>
//#include <string.h>
#include <gsl/gsl_rng.h>
//#include <sys/time.h>
//#include <unistd.h>
#include <stdbool.h>

#include "systems.h"
#include "utilities.h"
#include "interactions.h"

//#include "str_builder.h"

void corePeriodic(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double r_cut_off_torque_2, double gamma_pp, double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n){

    for (index_p = 0; index_p < n_particles; index_p++){
        fx_b = 0;
        fy_b = 0;
        torque_b = 0;
        double r_cut_off_force_2 = 1.21;

        for (index_n = index_p + 1; index_n < n_particles; index_n++){
            delta_x = x[index_p]-x[index_n];
            if (fabs(delta_x) > (l-sqrt(r_cut_off_force_2))){
                delta_x = delta_x - l*(1.0-2.0*signbit(delta_x));
            }

            delta_y = y[index_p]-y[index_n];
            if (fabs(delta_y) > (h-sqrt(r_cut_off_force_2))){
                delta_y = delta_y - h*(1.0-2.0*signbit(delta_y));
            }

            r_pn_2 = delta_x*delta_x + delta_y*delta_y;

            //if (r_pn_2 < 1.08683418){
            //if (r_pn_2 < 2*sigma_pp*sigma_pp) {
            //if (r_pn_2 < 4.0) {
            //if (r_pn_2 < r_cut_off_torque_2){

            if (r_pn_2 < r_cut_off_force_2){

                    //if (r_pn_2 < 1.0){

                torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], gamma_pp, r_pn_2);
                torque_n[index_p] += temp_torque_n;
                torque_n[index_n] -= temp_torque_n;
                number_n[index_n]++;
                number_n[index_p]++;
                //}

                //forceWeeksChandlerAndersen(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                //forceLennardJonesShifted(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                //forceLennardJonesRepAndExpRep(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                forceHarmonicPP(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, sqrt(r_cut_off_force_2), 10.0);
                //forceOneOverRQuad(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                //forceOneOverRQuadSig(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, sigma_pp);

                fx_n[index_p] += temp_fx_n;
                fy_n[index_p] += temp_fy_n;
                fx_n[index_n] -= temp_fx_n;
                fy_n[index_n] -= temp_fy_n;
                deformation_n[index_p] += r_pn_2;
                deformation_n[index_n] += r_pn_2;
                //}
            }


        } //For index_n


        //Update Adams-Bashforth helping parameters
        upadteIntegrationParameters(&Y_x[index_p], &Y_y[index_p], &Y_th[index_p], &Y_x_prev[index_p], &Y_y_prev[index_p], &Y_th_prev[index_p], t, theta[index_p], fx_b, fy_b, torque_b, fx_n[index_p], fy_n[index_p], torque_n[index_p], number_n[index_p], fs_scale, u_0);

        //Update particle parameters
        if (index_p >= n_fixed_particles){
            updateParticleParametersPeriodicXY(&x[index_p], &y[index_p], &theta[index_p], &vx[index_p], &vy[index_p], Y_x[index_p], Y_y[index_p], Y_th[index_p], Y_x_prev[index_p], Y_y_prev[index_p], Y_th_prev[index_p], f_AB1, f_AB2, dt, d_r, a, r, l, h);
        }
    }
}



void corePeriodicTube(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double r_cut_off_torque_2, double gamma_pp, double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n){

    for (index_p = 0; index_p < n_particles; index_p++){
        fx_b = 0;
        fy_b = 0;
        torque_b = 0;
        if (fabs(y[index_p])> h/2){
            forceHarmonicOneD(&fy_b, y[index_p], (1.0-2.0*signbit(y[index_p]))*h/2, lambda_har);
            torqueHarmonicOneD(&torque_b, theta[index_p], M_PI/2, lambda_har, kappa_har);
        }

        //
        for (index_n = index_p + 1; index_n < n_particles; index_n++){
            delta_x = x[index_p]-x[index_n];
            if (fabs(delta_x) > (l-sqrt(2)*sigma_pp)){
                delta_x = delta_x - l*(1.0-2.0*signbit(delta_x));
            }
            delta_y = y[index_p]-y[index_n];
            r_pn_2 = delta_x*delta_x + delta_y*delta_y;

            if (r_pn_2 < 2*sigma_pp*sigma_pp) {
                if (r_pn_2 < r_cut_off_torque_2){
                    //if (r_pn_2 < R_CUT_OFF_FORCE*R_CUT_OFF_FORCE){
                    //if (r_pn_2 < 1.0){

                    torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], gamma_pp, r_pn_2);
                    torque_n[index_p] += temp_torque_n;
                    torque_n[index_n] -= temp_torque_n;
                    number_n[index_n]++;
                    number_n[index_p]++;
                }
                //forceHarmonicPP(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, R_CUT_OFF_FORCE, LAMBDA_PP);
                //forceOneOverRQuad(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                forceOneOverRQuadSig(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, sigma_pp);
                fx_n[index_p] += temp_fx_n;
                fy_n[index_p] += temp_fy_n;
                fx_n[index_n] -= temp_fx_n;
                fy_n[index_n] -= temp_fy_n;
                deformation_n[index_p] += r_pn_2;
                deformation_n[index_n] += r_pn_2;
            }


        } //For index_n



        //Update Adams-Bashforth helping parameters
        upadteIntegrationParameters(&Y_x[index_p], &Y_y[index_p], &Y_th[index_p], &Y_x_prev[index_p], &Y_y_prev[index_p], &Y_th_prev[index_p], t, theta[index_p], fx_b, fy_b, torque_b, fx_n[index_p], fy_n[index_p], torque_n[index_p], number_n[index_p], fs_scale, u_0);

        //Update particle parameters
        if (index_p >= n_fixed_particles){
            updateParticleParametersPeriodicX(&x[index_p], &y[index_p], &theta[index_p], &vx[index_p], &vy[index_p], Y_x[index_p], Y_y[index_p], Y_th[index_p], Y_x_prev[index_p], Y_y_prev[index_p], Y_th_prev[index_p], f_AB1, f_AB2, dt, d_r, a, r, l);
        }
    }
}

void corePeriodicFunnel(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double h_funnel, double lambda_har, double kappa_har, double delta_x, double delta_y, double r_pn_2, double r_cut_off_torque_2, double gamma_pp, double sigma_pp, double temp_fx_n, double temp_fy_n, double temp_torque_n){

    for (index_p = 0; index_p < n_particles; index_p++){
        fx_b = 0;
        fy_b = 0;
        torque_b = 0;
        if (fabs(y[index_p])> (fabs(x[index_p])*(h-h_funnel)/l)+h_funnel/2){
            forceHarmonicFunnel(&fx_b, &fy_b, x[index_p], y[index_p], l,  h, h_funnel, lambda_har);
            //NB, angle is wrong here...
            torqueHarmonicOneD(&torque_b, theta[index_p], M_PI/2, lambda_har, kappa_har);
        }

        //
        for (index_n = index_p + 1; index_n < n_particles; index_n++){
            delta_x = x[index_p]-x[index_n];
            if (fabs(delta_x) > (l-sqrt(2)*sigma_pp)){
                delta_x = delta_x - l*(1.0-2.0*signbit(delta_x));
            }
            delta_y = y[index_p]-y[index_n];
            r_pn_2 = delta_x*delta_x + delta_y*delta_y;

            if (r_pn_2 < 2*sigma_pp*sigma_pp) {
                if (r_pn_2 < r_cut_off_torque_2){
                    //if (r_pn_2 < R_CUT_OFF_FORCE*R_CUT_OFF_FORCE){
                    //if (r_pn_2 < 1.0){

                    torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], gamma_pp, r_pn_2);
                    torque_n[index_p] += temp_torque_n;
                    torque_n[index_n] -= temp_torque_n;
                    number_n[index_n]++;
                    number_n[index_p]++;
                }
                //forceHarmonicPP(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, R_CUT_OFF_FORCE, LAMBDA_PP);
                //forceOneOverRQuad(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                forceOneOverRQuadSig(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, sigma_pp);
                fx_n[index_p] += temp_fx_n;
                fy_n[index_p] += temp_fy_n;
                fx_n[index_n] -= temp_fx_n;
                fy_n[index_n] -= temp_fy_n;
                deformation_n[index_p] += r_pn_2;
                deformation_n[index_n] += r_pn_2;
            }


        } //For index_n



        //Update Adams-Bashforth helping parameters
        upadteIntegrationParameters(&Y_x[index_p], &Y_y[index_p], &Y_th[index_p], &Y_x_prev[index_p], &Y_y_prev[index_p], &Y_th_prev[index_p], t, theta[index_p], fx_b, fy_b, torque_b, fx_n[index_p], fy_n[index_p], torque_n[index_p], number_n[index_p], fs_scale, u_0);

        //Update particle parameters
        if (index_p >= n_fixed_particles){
            updateParticleParametersPeriodicX(&x[index_p], &y[index_p], &theta[index_p], &vx[index_p], &vy[index_p], Y_x[index_p], Y_y[index_p], Y_th[index_p], Y_x_prev[index_p], Y_y_prev[index_p], Y_th_prev[index_p], f_AB1, f_AB2, dt, d_r, a, r, l);
        }
    }
}

void upadteIntegrationParameters(double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int t, double theta, double fx_b, double fy_b, double torque_b, double fx_n, double fy_n, double torque_n, double number_n, double fs_scale, double u_0){
    *Y_x = u_0*cos(theta) + (fx_b + fx_n)*fs_scale;
    *Y_y = u_0*sin(theta) + (fy_b + fy_n)*fs_scale;
    if (number_n > 0){
        *Y_th = (torque_b + torque_n/number_n)*fs_scale;
    } else {
        *Y_th = torque_b*fs_scale;
    }

    if (t==1){
        Y_x_prev = Y_x;
        Y_y_prev = Y_y;
        Y_th_prev = Y_th;
    }
}

void updateParticleParameters(double* x, double*y, double* theta, double*vx, double* vy, double Y_x, double Y_y, double Y_th, double Y_x_prev, double Y_y_prev, double Y_th_prev, double f_AB1, double f_AB2, double dt, double d_r, double a, gsl_rng **r){
    *x = *x + (f_AB1*Y_x - f_AB2*Y_x_prev)*dt;
    *y = *y + (f_AB1*Y_y - f_AB2*Y_y_prev)*dt;
    *theta = *theta + (f_AB1*Y_th - f_AB2*Y_th_prev)*dt + sqrt(2*d_r*dt)*randDouble(-a, a, r);//&r);
    *vx = f_AB1*Y_x - f_AB2*Y_x_prev;
    *vy = f_AB1*Y_y - f_AB2*Y_y_prev;
}

void updateParticleParametersPeriodicX(double* x, double*y, double* theta, double*vx, double* vy, double Y_x, double Y_y, double Y_th, double Y_x_prev, double Y_y_prev, double Y_th_prev, double f_AB1, double f_AB2, double dt, double d_r, double a, gsl_rng **r, double l){
    updateParticleParameters(x, y, theta, vx, vy, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, f_AB1, f_AB2, dt, d_r, a, r);
    if (*x <-l/2){
        *x +=l;
    } else if (*x >l/2){
        *x -=l;
    }
}

void updateParticleParametersPeriodicXY(double* x, double*y, double* theta, double*vx, double* vy, double Y_x, double Y_y, double Y_th, double Y_x_prev, double Y_y_prev, double Y_th_prev, double f_AB1, double f_AB2, double dt, double d_r, double a, gsl_rng **r, double l, double h){
    updateParticleParameters(x, y, theta, vx, vy, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, f_AB1, f_AB2, dt, d_r, a, r);
    if (*x <-l/2){
        *x +=l;
    } else if (*x >l/2){
        *x -=l;
    }

    if (*y < -h/2){
        *y +=h;
    } else if (*y > h/2){
        *y -=h;
    }
}
