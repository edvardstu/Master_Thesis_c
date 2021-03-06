#ifndef INTERACTIONS_H
#define INTERACTIONS_H

void forceHarmonicCircular(double *fx_p, double *fy_p, double r_coord, double x, double y, double r_boundary, double lambda_harmonic);

void forceHarmonicInfWell(double *fx_b, double *fy_b, double x, double y, double L, double lambda_harmonic);

void forceHarmonicOneD(double *fx_b, double x, double x_boundary, double lambda_harmonic);

void forceHarmonicFunnel(double *fx_b, double *fy_b, double x, double y, double l, double h, double h_funnel, double lambda_harmonic);

void torqueHarmonicCircular(double *torque_b, double r_coord, double x, double y, double theta, double lambda_harmonic, double kappa_harmonic);

void torqueHarmonicOneD(double *torque_b, double theta, double beta, double lambda_harmonic, double kappa_harmonic);

void torqueHarmonicInfWell(double *torque_b, double x, double y, double theta, double L, double lambda_harmonic, double kappa_harmonic);

void forceWeeksChandlerAndersen(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y);

void forceLennardJonesShifted(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y);

void forceLennardJonesRepAndExpRep(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y);

void forceHarmonicPP(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double r_cut_off_force, double lambda_pp);

void forceOneOverRSquared(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double sigma_pp);

void forceOneOverRQuad(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y);

void forceOneOverRQuadSig(double *fx_n, double *fy_n, double r_pn_2, double delta_x, double delta_y, double sigma_pp);

/*void torqueWeeksChandlerAndersen(double *torque_n, double theta_p, double theta_n, double gamma_pp);*/

void torqueWeeksChandlerAndersen(double *torque_n, double theta_p, double theta_n, double gamma_pp, double r_pn_2);

#endif
