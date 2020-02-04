#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "interactions.h"

//In all the following functions it is assumed that the particles are all ready within range of the potentials
//This is to improve computational time

void forceHarmonicCircular(double *fx_b, double *fy_b, double r_coord, double x, double y, double r_boundary, double lambda_harmonic){
    //Assume the circle is placed in (x,y)=(0,0)
    //double beta = atan2(y, x);
    double f = lambda_harmonic*(r_coord-r_boundary);

    *fx_b = -f*x/r_coord;
    *fy_b = -f*y/r_coord;
}

void forceHarmonicOneD(double *fx_b, double x, double x_boundary, double lambda_harmonic){
    *fx_b = -lambda_harmonic*(x-x_boundary);
}


void forceHarmonicInfWell(double *fx_b, double *fy_b, double x, double y, double L, double lambda_harmonic){
    if (x<0.0) {
        forceHarmonicOneD(fx_b, x, 0.0, lambda_harmonic);
    } else if (x>L) {
        forceHarmonicOneD(fx_b, x, L, lambda_harmonic);
    }

    if (y<0.0){
        forceHarmonicOneD(fy_b, y, 0.0, lambda_harmonic);
    }
}

void torqueHarmonicCircular(double *torque_b, double r_coord, double x, double y, double theta, double lambda_harmonic, double kappa_harmonic){
    double beta = atan2(y, x);
    *torque_b = lambda_harmonic*kappa_harmonic*sin(2*(theta-beta));
}

void torqueHarmonicInfWell(double *torque_b, double x, double y, double theta, double L, double lambda_harmonic, double kappa_harmonic){
    if ((x<0.0) || (x>L)){
        *torque_b = lambda_harmonic*kappa_harmonic*sin(2*theta);
    } else if (y<0.0){
        *torque_b = lambda_harmonic*kappa_harmonic*sin(2*theta-M_PI);
    }
}

void forceWeeksChandlerAndersen(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y){
    //Here it is assumed that sigma=1/2^(1/6) which gives a particle radius of 1
    double r_pn_6 = r_pn_2*r_pn_2*r_pn_2;
    double f = 12*(1.0/r_pn_6-1.0)/(r_pn_6*r_pn_2);
    *fx_n = f*delta_x;
    *fy_n = f*delta_y;
}

void forceHarmonicPP(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double r_cut_off_force, double lambda_pp){
    double f = -lambda_pp*(sqrt(r_pn_2)-r_cut_off_force);
    *fx_n = f*delta_x/sqrt(r_pn_2);
    *fy_n = f*delta_y/sqrt(r_pn_2);
}

void forceOneOverRSquared(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double sigma_pp){
    double r_pn = sqrt(r_pn_2);
    double r_pn_3 = r_pn_2*r_pn;
    double f = 4*(2.0*sigma_pp*sigma_pp/r_pn-sigma_pp)/r_pn_3;
    *fx_n = f*delta_x;
    *fy_n = f*delta_y;
}

void forceOneOverRQuad(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y){
    //Here it is assumed that sigma=1/2^(1/6) which gives a particle radius of 1
    double f = 12*(1.0/r_pn_2-1.0)/r_pn_2;
    *fx_n = f*delta_x;
    *fy_n = f*delta_y;
}

void forceOneOverRQuadSig(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double sigma_pp){
    //Here it is assumed that sigma=1/2^(1/6) which gives a particle radius of 1
    double f = 12*((2.0*sigma_pp*sigma_pp)/r_pn_2-1.0)/r_pn_2;
    *fx_n = f*delta_x;
    *fy_n = f*delta_y;
}

/*void torqueWeeksChandlerAndersen(double *torque_n, double theta_p, double theta_n, double gamma_pp){
    *torque_n = gamma_pp*sin(theta_n-theta_p);
}*/

void torqueWeeksChandlerAndersen(double *torque_n, double theta_p, double theta_n, double gamma_pp, double r_pn_2){
    *torque_n = (1/r_pn_2)*gamma_pp*sin(theta_n-theta_p);
    //*torque_n = (1/r_pn_2)*sin(theta_n-theta_p);
}
