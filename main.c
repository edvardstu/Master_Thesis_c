#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include <omp.h>
#include <stdbool.h>

#include "hdf5.h"
#include "utilities.h"
#include "interactions.h"
#include "str_builder.h"
#include "systems.h"
#include "hdf5ReaderWriter.h"



//System and particle parameters
#define R 30.0//17.0
#define L 25.0
#define H 25.0
#define H_FUNNEL 9.3
#define R_PARTICLE 0.5
#define N_PARTICLES_I 1000
#define N_FIXED_PARTICLES 0
#define N_A 500
#define U_0_I 1.0
#define D_R_C 0.001 //0.2

#define N_STEPS 1000
#define DT 0.005 //0.0005 for U_0=10.0

//Diffusive parameters
#define GAMMA_T 1
#define GAMMA_R 1

//Boundary interatction
#define LAMBDA_HAR 10.0//20.0//200.0 //FS
#define KAPPA_HAR 10.0  //GS

//Particle particle interaction
#define GAMMA_PP 0.0//1.0 //GAM
#define LAMBDA_PP 20.0
#define R_CUT_OFF_FORCE_2 1.1025  //
#define R_CUT_OFF_TORQUE_2 1.1025 //9.0
const double SIGMA_PP = 1/sqrt(2);//1.5;//sqrt(2.0);//

double R_CUT_OFF_FORCE_2_ARRAY[] = {0.7225, 1.1025, 1.5625}; //AA, AB, BB

double R_CUT_OFF_FORCE_2_MARTIX[N_PARTICLES_I][N_PARTICLES_I];

//double** R_CUT_OFF_FORCE_2_MARTIX = (double **)malloc(N_PARTICLES_I * sizeof(double*));
//for (int i = 0; i < N_PARTICLES_I; i++) R_CUT_OFF_FORCE_2_MARTIX[i] = (double *)malloc(N_PARTICLES_I * sizeof(double));
//double (*R_CUT_OFF_FORCE_2_MARTIX)[N_PARTICLES_I] = malloc(N_PARTICLES_I * sizeof( *R_CUT_OFF_FORCE_2_MARTIX));

const double a = sqrt(3);

//Infinite well variables
//#define L 50.0
#define U_S 0.5

//Type of boundary for simulation
enum barrier {Circular, Periodic, PeriodicTube, PeriodicFunnel, PeriodicBM};
//Parameter to sweep
enum sweep {None, ParticleNumber, DiffusionCoefficient, SelfPropulsion};


int main(int argc, char **argv) {



    #ifdef _OPENMP
        printf("OpenMP is defined\n");
    #else
        printf("OpenMP is not defined\n");
    #endif



    /*
    //#pragma omp parallel for
    for (int i=0; i<10; i++){
        if (i==0){
            int threads = omp_get_num_threads();
            printf("The number of threads is: %d\n", threads);
        }
        int p = omp_get_thread_num();
        printf("Thread number: %i\n", p);
    }*/

    double time_start = walltime();
    bool useAB = false;
    bool rndSeed = true;
    const bool overwrite = false;
    bool continueFromPrev = false;
    enum barrier simulationBarrier = PeriodicBM;
    enum sweep simulationSweep = SelfPropulsion;
    int writeInterval = 50;
    //const char * restrict fileNameBaseInit = "results/periodic_2D/phaseDiagram/gamma_pp_0_1/test";
    //const char * restrict fileNameBaseInit = "results/periodic_2D/phaseDiagram/gamma_pp_0_0/01/sweepPropulsion";
    //const char * restrict fileNameBaseInit = "/media/edvardst/My Book/NTNU/Programming/Master_Thesis/Periodic_2D/phaseDiagram/gamma_pp_0_0/highDensity/sweepPropulsion";
    const char * restrict fileNameBaseInit = "/home/edvardst/Documents/NTNU/Programming/Master_Thesis/Master_Thesis_c/results/paperTest/dummy";







    double f_AB1 = 1.5;
    double f_AB2 = 0.5;
    if (!useAB){
        f_AB1 = 1.0;
        f_AB2 = 0.0;
    }



    int numberOfTimeframes;
    //double U_0_array[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};
    double U_0_array[] = {1.0};
    //double U_0_array[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0};
    //Timeframe setup
    if (simulationSweep == SelfPropulsion){
        numberOfTimeframes = sizeof U_0_array / sizeof *U_0_array;
    } else {
        numberOfTimeframes = 1;
    }


    for (int k=0; k<numberOfTimeframes; k++){
    printf("Timeframe %d of %d\n", k+1, numberOfTimeframes);

    ///////////////////Setting up sweep variables /////////////////////
    double D_R = D_R_C;
    double D_R_I;
    int N_PARTICLES = N_PARTICLES_I;
    double U_0 = U_0_I;
    double time=0.0;

    switch (simulationSweep) {
        case None:
            //Nothing
            break;
        case DiffusionCoefficient:
            //Should be something here...
            break;
        case ParticleNumber:
            N_PARTICLES = N_PARTICLES_I*(k+1);
            break;
        case SelfPropulsion:
            U_0 = U_0_array[k];
            break;
    }
    ///////////////////////////////////////////////////////////////////

    //Calculate density. Assume particle radius is 0.5.
    //Implies that p-p-WCA-interaction has a cutoff distance of 1.0
    //Assume random close packing of 0.82 (https://www.sciencedirect.com/science/article/abs/pii/0021979771903389)

    /*
    //Check the number of particles compared to the size of the system
    double r_particle = (R_CUT_OFF_FORCE-U_0/LAMBDA_PP)/2;
    printf("The particle radius is %f\n", r_particle);
    int max_p = floor((R*R*0.9069)/(R_PARTICLE*R_PARTICLE));
    if (N_PARTICLES>max_p){
        printf("Too many particles\n");
        printf("Maximum number of particles is %d\n", max_p);
        //exit(-1);
    } else {
        printf("The system occupancy is %f\n", (double)N_PARTICLES/max_p);
    }
    */


    //////////////////////Setting up variables ////////////////////////
    FILE *fp;//, *fp2;
    double x[N_PARTICLES], y[N_PARTICLES], theta[N_PARTICLES], vx[N_PARTICLES], vy[N_PARTICLES];
    double fs_scale;
    //Loop parameters
    unsigned int t, index_p, index_n, i;
    //Parameters for boundary
    double r_coord, fx_b, fy_b, torque_b;
    //Parameters for particle particle interaction
    double delta_x, delta_y, temp_fx_n, temp_fy_n, temp_torque_n, r_pn_2;
    //Parameters for four point susceptibility
    double omega;
    //Helping variables for Adams_Bashforth
    double *Y_x = malloc(N_PARTICLES*sizeof(double));
    double *Y_x_prev = malloc(N_PARTICLES*sizeof(double));
    double *Y_y = malloc(N_PARTICLES*sizeof(double));
    double *Y_y_prev = malloc(N_PARTICLES*sizeof(double));
    double *Y_th = malloc(N_PARTICLES*sizeof(double));
    double *Y_th_prev = malloc(N_PARTICLES*sizeof(double));

    //Parameters for loading previos Results
    const char * restrict fileNameBase;
    const char * restrict fileName;
    const char * restrict fileNameHDF5;
    const char * restrict fileNamePrevious;
    // const char * restrict fileNameSusc;

    //bool continueFromPrev = true;

    ///////////////////////////////////////////////////////////////////


    //////////////////////// System fractions /////////////////////////
    double particleArea = M_PI*0.25;
    double areaFaction = particleArea*N_PARTICLES/(H*L);
    double packingFraction = areaFaction/0.82;
    printf("The system has an area fraction of: %f, and a packing fraction of: %f\n",areaFaction, packingFraction);
    double particleAreaA = M_PI_4*sqrt(R_CUT_OFF_FORCE_2_ARRAY[0]);
    double particleAreaB = M_PI_4*sqrt(R_CUT_OFF_FORCE_2_ARRAY[2]);
    areaFaction = (particleAreaA*N_A + particleAreaB*(N_PARTICLES-N_A))/(H*L);
    packingFraction = areaFaction/0.82;
    printf("The binary mixture has an area fraction of: %f, and a packing fraction of: %f\n",areaFaction, packingFraction);
    ///////////////////////////////////////////////////////////////////

    ////////////////////// Creating file names ////////////////////////
    fileNamePrevious = createFileNamePrevious(fileNameBaseInit, overwrite);
    fileNameBase = createFileNameBase(fileNameBaseInit, overwrite);
    fileName = createFileName(fileNameBase);
    fileNameHDF5 = createFileNameH5(fileNameBase);
    //fileNameSusc = createFileNameSucs(fileNameBase);
    ///////////////////////////////////////////////////////////////////

    //////////////////////Setting up GSL RNG ////////////////////////
    const gsl_rng_type * T;
    gsl_rng * r;
    setUpRNG(&T, &r, rndSeed);
    /////////////////////////////////////////////////////////////////

    //////////Setting up particle position and angle ////////////////
    if (!continueFromPrev){
        switch (simulationBarrier) {
            case Circular:
                sunflower(x, y, N_PARTICLES, 0, R);
                //sunflower_fixed_boundary(x, y, N_PARTICLES, 0, R+1.0, N_FIXED_PARTICLES, R*0.95);
                break;
            case Periodic:
                uniformRectangle(x, y, N_PARTICLES, L, H, &r);
                break;
            case PeriodicTube:
                uniformRectangle(x, y, N_PARTICLES, L, H, &r);
                break;
            case PeriodicFunnel:
                uniformFunnel(x, y, N_PARTICLES, L, H, H_FUNNEL, &r);
                fixedBoundaryFunnel(x, y, theta, N_PARTICLES, N_FIXED_PARTICLES, L, H, H_FUNNEL);
                break;
            case PeriodicBM:
                uniformRectangle(x, y, N_PARTICLES, L, H, &r);
                for (i=0; i<N_PARTICLES; i++){
                    for (int j=0; j<N_PARTICLES; j++){
                        if ((i < N_A) && (j < N_A)){
                            R_CUT_OFF_FORCE_2_MARTIX[i][j] = R_CUT_OFF_FORCE_2_ARRAY[0];
                        } else if ((i >= N_A) && (j >= N_A)){
                            R_CUT_OFF_FORCE_2_MARTIX[i][j] = R_CUT_OFF_FORCE_2_ARRAY[2];
                        } else {
                            R_CUT_OFF_FORCE_2_MARTIX[i][j] = R_CUT_OFF_FORCE_2_ARRAY[1];
                        }
                    }
                }
                break;
        }
        for (i=0; i<N_FIXED_PARTICLES; i++){
            //theta[i] = atan2(y[i], x[i]) + M_PI/2;
            //theta[i]=randDouble(-M_PI, M_PI, &r);
            vx[i] = 0.1*U_0*cos(theta[i]);
            vy[i] = 0.1*U_0*sin(theta[i]);
        }
        for (i=N_FIXED_PARTICLES;i<N_PARTICLES;i++){
            theta[i]=randDouble(-M_PI, M_PI, &r);
            vx[i] = U_0*cos(theta[i]);
            vy[i] = U_0*sin(theta[i]);
        }
    } else {
        readInitialState(fileNamePrevious, N_PARTICLES, &time, x, y, theta, vx, vy, &D_R_I);
        D_R = D_R_I;
        if (D_R_C != D_R_I){
            printf("D_r was initilized to %.3f, but change to %.3f\n", D_R_C, D_R_I);
        }
    }
    /////////////////////////////////////////////////////////////////

    ///////////////////Opening file for results ///////////////////////
    writeSimulationParameters(fileNameBase, R, L, H, H_FUNNEL, R_PARTICLE, N_PARTICLES, N_FIXED_PARTICLES, U_0, D_R, N_STEPS, DT, GAMMA_T, GAMMA_R, LAMBDA_HAR, KAPPA_HAR, GAMMA_PP, LAMBDA_PP, R_CUT_OFF_FORCE_2, R_CUT_OFF_TORQUE_2, SIGMA_PP);

    //openFile(fileName, &fp);
    // openFile(fileNameSusc, &fp2);
    //for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n", i, time, x[i], y[i], theta[i], vx[i], vy[i], D_R, 0.0);
    /////////////////////////////////////////////////////////////////


    //////////////// Parameters for HDF5 writer ////////////////////////////////
    hid_t file_id, dataset_id;
    n_cols = 7;
    const hsize_t n_rows_slab = N_PARTICLES;
    hsize_t chunk_dims[2] = {n_rows_slab, n_cols};
    hid_t mem_space = H5Screate_simple(2, chunk_dims, NULL);
    herr_t status;
    hsize_t dims[2];

    double buffer[n_rows_slab][n_cols];
    hsize_t count[2] = {n_rows_slab, n_cols};
    hsize_t start[2] = {0, 0};
    ////////////////////////////////////////////////////////////////////////////

    //////////////////////// Create HDF5 file //////////////////////////////////
    createAndSetUpH5Env(&file_id, &dataset_id, chunk_dims, fileNameHDF5);
    ////////////////////////////////////////////////////////////////////////////

    /////////////////////// Write inital state /////////////////////////////////
    transposeSimulationData(buffer, N_PARTICLES, time, x, y, theta, vx, vy);
    int writeNumber = 0;
    extendDatasetH5(&dataset_id, chunk_dims, writeNumber);
    appendBufferToDatasetH5(buffer, &dataset_id, mem_space, start, count, writeNumber, n_rows_slab);
    ////////////////////////////////////////////////////////////////////////////

    fs_scale = 0.0;
    //For each time step
    double n_scale_steps = 1000; //50000;

    double fx_n[N_PARTICLES], fy_n[N_PARTICLES], torque_n[N_PARTICLES], number_n[N_PARTICLES], deformation_n[N_PARTICLES];


    //#pragma omp parallel default(shared) private(index_p, index_n, fx_b, fy_b, torque_b, r_coord, delta_x, delta_y, r_pn_2, temp_torque_n, temp_fx_n, temp_fy_n)
    //#pragma omp parallel default(shared) num_threads(1)
    for (t = 1; t <= N_STEPS; t++){
        int thread_n = omp_get_thread_num();
        if (t % 10000 ==0 ) printf("Thread %d, step %d\n",thread_n, t);

        if (thread_n == 0){
        //if (t % 10000 ==0 ) printf("%d\n",t);


        if (t <= n_scale_steps){
            if (!continueFromPrev){
                fs_scale=t/n_scale_steps;
            } else {
                fs_scale=1.0;
                D_R = D_R_I + 0.2*(t/n_scale_steps);
            }
        }
        //Should not be here, only for solver testing
        //fs_scale=1.0;

        for (int j = 0; j<N_PARTICLES; j++){
            fx_n[j] = 0;
            fy_n[j] = 0;
            torque_n[j] = 0;
            number_n[j] = 0;
            deformation_n[j] = 0;
        }
        // double fx_n[N_PARTICLES] = {0};
        // double fy_n[N_PARTICLES] = {0};
        // double torque_n[N_PARTICLES] = {0};
        // double number_n[N_PARTICLES] = {0};
        // double deformation_n[N_PARTICLES] = {0};

        }
        //#pragma omp barrier

        //#pragma omp for private(index_p, index_n, fx_b, fy_b, torque_b, r_coord, delta_x, delta_y, r_pn_2, temp_torque_n, temp_fx_n, temp_fy_n)
        //#pragma omp for default(shared) private(index_p, index_n, fx_b, fy_b, torque_b, r_coord, delta_x, delta_y, r_pn_2, temp_torque_n, temp_fx_n, temp_fy_n)
        switch (simulationBarrier) {
        case Circular:
        for (index_p = 0; index_p < N_PARTICLES; index_p++){
            //Find forces and torque from wall
            //Assume circular potentail has a centre in (0,0)
            fx_b = 0;
            fy_b = 0;
            torque_b = 0;

            /* Infinite well
            forceHarmonicInfWell(&fx_b, &fy_b, x[index_p], y[index_p], L, LAMBDA_HAR);
            fy_b -= U_S;
            torqueHarmonicInfWell(&torque_b, x[index_p], y[index_p], theta[index_p], L, LAMBDA_HAR, KAPPA_HAR);
            */

            r_coord = sqrt(x[index_p]*x[index_p]+y[index_p]*y[index_p]);
            if (r_coord > R){
                forceHarmonicCircular(&fx_b, &fy_b, r_coord, x[index_p], y[index_p], R, LAMBDA_HAR);
                torqueHarmonicCircular(&torque_b, r_coord, x[index_p], y[index_p], theta[index_p], LAMBDA_HAR, KAPPA_HAR);
            }

            for (index_n = index_p + 1; index_n < N_PARTICLES; index_n++){
                delta_x = x[index_p]-x[index_n];
                delta_y = y[index_p]-y[index_n];
                r_pn_2 = delta_x*delta_x + delta_y*delta_y;

                if (r_pn_2 < 2*SIGMA_PP*SIGMA_PP) {
                    if (r_pn_2 < R_CUT_OFF_TORQUE_2){
                        //if (r_pn_2 < R_CUT_OFF_FORCE*R_CUT_OFF_FORCE){
                        //if (r_pn_2 < 1.0){

                        torqueWeeksChandlerAndersen(&temp_torque_n, theta[index_p], theta[index_n], GAMMA_PP, r_pn_2);
                        torque_n[index_p] += temp_torque_n;
                        torque_n[index_n] -= temp_torque_n;
                        number_n[index_n]++;
                        number_n[index_p]++;
                    }
                    //forceHarmonicPP(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, R_CUT_OFF_FORCE, LAMBDA_PP);
                    //forceOneOverRQuad(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y);
                    forceOneOverRQuadSig(&temp_fx_n, &temp_fy_n, r_pn_2, delta_x, delta_y, SIGMA_PP);
                    fx_n[index_p] += temp_fx_n;
                    fy_n[index_p] += temp_fy_n;
                    fx_n[index_n] -= temp_fx_n;
                    fy_n[index_n] -= temp_fy_n;
                    deformation_n[index_p] += r_pn_2;
                    deformation_n[index_n] += r_pn_2;
                }


            } //For index_n

            //Update Adams-Bashforth helping parameters
            Y_x[index_p] = U_0*cos(theta[index_p]) + (fx_b + fx_n[index_p])*fs_scale;
            Y_y[index_p] = U_0*sin(theta[index_p]) + (fy_b + fy_n[index_p])*fs_scale;
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
            if (index_p >= N_FIXED_PARTICLES){
                x[index_p] = x[index_p] + (f_AB1*Y_x[index_p] - f_AB2*Y_x_prev[index_p])*DT;
                y[index_p] = y[index_p] + (f_AB1*Y_y[index_p] - f_AB2*Y_y_prev[index_p])*DT;
                theta[index_p] = theta[index_p] + (f_AB1*Y_th[index_p] - f_AB2*Y_th_prev[index_p])*DT + sqrt(2*D_R*DT)*randDouble(-a, a, &r);

                vx[index_p] = f_AB1*Y_x[index_p] - f_AB2*Y_x_prev[index_p];
                vy[index_p] = f_AB1*Y_y[index_p] - f_AB2*Y_y_prev[index_p];
            }

        } //For index_p
            break;

            case Periodic:
                //corePeriodicTube(t, x, y, theta, vx, vy, fx_b, fy_b, torque_b, fx_n, fy_n, torque_n, number_n, deformation_n, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, index_p, index_n, fs_scale, f_AB1, f_AB2, &r, N_PARTICLES, N_FIXED_PARTICLES, DT, D_R, a, U_0, L, H, LAMBDA_HAR, KAPPA_HAR);
                corePeriodic(t, x, y, theta, vx, vy, fx_b, fy_b, torque_b, fx_n, fy_n, torque_n, number_n, deformation_n, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, index_p, index_n, fs_scale, f_AB1, f_AB2, &r, N_PARTICLES, N_FIXED_PARTICLES, DT, D_R, a, U_0, L, H, LAMBDA_HAR, KAPPA_HAR, delta_x, delta_y, r_pn_2, GAMMA_PP, LAMBDA_PP, R_CUT_OFF_FORCE_2, R_CUT_OFF_TORQUE_2, SIGMA_PP, temp_fx_n, temp_fy_n, temp_torque_n);
                break;

            case PeriodicTube:
                corePeriodicTube(t, x, y, theta, vx, vy, fx_b, fy_b, torque_b, fx_n, fy_n, torque_n, number_n, deformation_n, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, index_p, index_n, fs_scale, f_AB1, f_AB2, &r, N_PARTICLES, N_FIXED_PARTICLES, DT, D_R, a, U_0, L, H, LAMBDA_HAR, KAPPA_HAR, delta_x, delta_y, r_pn_2, R_CUT_OFF_TORQUE_2, GAMMA_PP, SIGMA_PP, temp_fx_n, temp_fy_n, temp_torque_n);
                break;

            case PeriodicFunnel:
                corePeriodicFunnel(t, x, y, theta, vx, vy, fx_b, fy_b, torque_b, fx_n, fy_n, torque_n, number_n, deformation_n, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, index_p, index_n, fs_scale, f_AB1, f_AB2, &r, N_PARTICLES, N_FIXED_PARTICLES, DT, D_R, a, U_0, L, H, H_FUNNEL, LAMBDA_HAR, KAPPA_HAR, delta_x, delta_y, r_pn_2, R_CUT_OFF_TORQUE_2, GAMMA_PP, SIGMA_PP, temp_fx_n, temp_fy_n, temp_torque_n);
                break;

            case PeriodicBM:
                corePeriodicBM(t, x, y, theta, vx, vy, fx_b, fy_b, torque_b, fx_n, fy_n, torque_n, number_n, deformation_n, Y_x, Y_y, Y_th, Y_x_prev, Y_y_prev, Y_th_prev, index_p, index_n, fs_scale, f_AB1, f_AB2, &r, N_PARTICLES, N_FIXED_PARTICLES, DT, D_R, a, U_0, L, H, LAMBDA_HAR, KAPPA_HAR, delta_x, delta_y, r_pn_2, GAMMA_PP, LAMBDA_PP, R_CUT_OFF_FORCE_2_MARTIX, R_CUT_OFF_TORQUE_2, SIGMA_PP, temp_fx_n, temp_fy_n, temp_torque_n, N_A);
                break;
        } // switch simulationBarrier

        //#pragma omp barrier
        if (thread_n == 0){
        time += DT;
        if (t % writeInterval == 0){
            // for (i=0;i<N_PARTICLES;i++) fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n", i, time, x[i], y[i], theta[i], vx[i], vy[i], D_R, deformation_n[i]);
            transposeSimulationData(buffer, N_PARTICLES, time, x, y, theta, vx, vy);
            int writeNumber = t/writeInterval;
            extendDatasetH5(&dataset_id, chunk_dims, writeNumber);
            appendBufferToDatasetH5(buffer, &dataset_id, mem_space, start, count, writeNumber, n_rows_slab);
        }


        // omega = 0.0;
        // for (i=0; i<N_PARTICLES; i++){
        //     delta_x = vx[i];
        //     delta_y = vy[i];
        //     if ((delta_x*delta_x+delta_y*delta_y)<(0.25*U_0*U_0)){
        //         omega+=1.0;
        //     }
        // }
        // omega/=N_PARTICLES;
        // fprintf(fp2, "%lf %lf\n", time, omega);

        //Swapping pointer of help parameters for AB
        swapPointers(&Y_x, &Y_x_prev);
        swapPointers(&Y_y, &Y_y_prev);
        swapPointers(&Y_th, &Y_th_prev);
        }
        //#pragma omp barrier
    } //End of for t
    writeFinalState(fileNameBase, N_PARTICLES, time, x, y, theta, vx, vy, D_R);



    ///////////////////Closing file with results //////////////////////
    //closeFile(fileName, &fp);
    // closeFile(fileNameSusc, &fp2);
    ///////////////////////////////////////////////////////////////////

    closeFilesH5(&mem_space, &dataset_id, &file_id);

    //Freeing RNG
    gsl_rng_free (r);

    double time_end = walltime();
    printf("Simulation time: %7.3f s\n",(time_end-time_start));
    switch (simulationSweep) {
        case DiffusionCoefficient:
            continueFromPrev=true;
        default:
            continueFromPrev=false;
    }

} //End of for k, timeFrames
//} //End of solver tester
return 0;
}
