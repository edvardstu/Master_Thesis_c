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

void corePeriodicTube(int t, double* x, double* y, double* theta, double* vx, double* vy, double fx_b, double fy_b, double torque_b, double* fx_n, double* fy_n, double* torque_n, double* number_n, double* deformation_n, double* Y_x, double* Y_y, double* Y_th, double* Y_x_prev, double* Y_y_prev, double* Y_th_prev, int index_p, int index_n, double fs_scale, double f_AB1, double f_AB2, gsl_rng** r, int n_particles, int n_fixed_particles, double dt, double d_r, double a, double u_0, double l, double h, double lambda_har, double kappa_har){

    for (index_p = 0; index_p < n_particles; index_p++){
        fx_b = 0;
        fy_b = 0;
        torque_b = 0;
        if (fabs(y[index_p])> h/2){
            forceHarmonicOneD(&fy_b, y[index_p], signbit(y[index_p])*h/2, lambda_har);
            torqueHarmonicOneD(&torque_b, theta[index_p], M_PI/2, lambda_har, kappa_har);
        }

        // for (index_n = index_p + 1; index_n < n_particles; index_n++){
        //     delta_x = x[index_p]-x[index_n];
        //     if (fabs(delta_x) > (l-sqrt(2)*SIGMA_PP)){
        //         dx = dx - l*signbit(dx);
        //     }
        //     delta_y = y[index_p]-y[index_n];
        //     r_pn_2 = delta_x*delta_x + delta_y*delta_y;
        //
        //     if (r_pn_2 < 2*SIGMA_PP*SIGMA_PP) {
        //         if (r_pn_2 < R_CUT_OFF_TORQUE_2){
        //             //if (r_pn_2 < R_CUT_OFF_FORCE*R_CUT_OFF_FORCE){
        //             //if (r_pn_2 < 1.0){
        //
        //             torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], GAMMA_PP, r_pn_2);
        //             torque_n[index_p] += temp_torque_n;
        //             torque_n[index_n] -= temp_torque_n;
        //             number_n[index_n]++;
        //             number_n[index_p]++;
        //         }
        //         //forceHarmonicPP(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, R_CUT_OFF_FORCE, LAMBDA_PP);
        //         //forceOneOverRQuad(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
        //         forceOneOverRQuadSig(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, SIGMA_PP);
        //         fx_n[index_p] += temp_fx_n;
        //         fy_n[index_p] += temp_fy_n;
        //         fx_n[index_n] -= temp_fx_n;
        //         fy_n[index_n] -= temp_fy_n;
        //         deformation_n[index_p] += r_pn_2;
        //         deformation_n[index_n] += r_pn_2;
        //     }
        //
        //
        // } //For index_n



        //Update Adams-Bashforth helping parameters
        Y_x[index_p] = u_0*cos(theta[index_p]) + (fx_b + fx_n[index_p])*fs_scale;
        Y_y[index_p] = u_0*sin(theta[index_p]) + (fy_b + fy_n[index_p])*fs_scale;
        if (number_n[index_p] > 0){
            Y_th[index_p] = (torque_b + torque_n[index_p]/number_n[index_p])*fs_scale;
        } else {
            Y_th[index_p] = torque_b*fs_scale;
        }

        if (t==1){
            Y_x_prev[index_p] = Y_x[index_p];
            Y_y_prev[index_p] = Y_y[index_p];
            Y_th_prev[index_p] = Y_th[index_p];
        }

        //Update particle parameters
        if (index_p >= n_fixed_particles){
            x[index_p] = x[index_p] + (f_AB1*Y_x[index_p] - f_AB2*Y_x_prev[index_p])*dt;
            if (x[index_p]<-l/2){
                x[index_p]+=l;
            } else if (x[index_p]>l/2){
                x[index_p]-=l;
            }
            y[index_p] = y[index_p] + (f_AB1*Y_y[index_p] - f_AB2*Y_y_prev[index_p])*dt;
            theta[index_p] = theta[index_p] + (f_AB1*Y_th[index_p] - f_AB2*Y_th_prev[index_p])*dt + sqrt(2*d_r*dt)*randDouble(-a, a, r);//&r);
            vx[index_p] = f_AB1*Y_x[index_p] - f_AB2*Y_x_prev[index_p];
            vy[index_p] = f_AB1*Y_y[index_p] - f_AB2*Y_y_prev[index_p];
        }
    }
}
