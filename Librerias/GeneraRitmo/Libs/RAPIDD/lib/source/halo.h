/*Please include all definitions of the functions and variables used in halo.c
  that you can use in other parts of the code (needed for correct linking)*/

#pragma once
#include <gsl/gsl_spline.h>

int f_shm(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval);
double shm_halo (double vmin, double vesc, double v0, double ve, double beta, int i);
double shm_halo_analitic(double velmin, double ve, double vesc, double v0, double beta);
double norm_shm(double vesc, double v0, double beta);
double Norm_fornasa (double vesc, double vt, double vc, double k);
double fornasa_halo (double vmin, double vesc, double vt, double vc, double ve, double k, int i);
double test_fornasa (double v, double vesc, double vt, double vc, double k, int i);

double lisanti_halo (double vmin, double vesc, double v0, double ve, double k, int i);
double Nk_uncert (double vesc, double v0, double k);
double halo (double vmin, int i);
double halo_w (double vmin, gsl_interp_accel *ga, gsl_spline * gs);
double halo_f (char * profile, double vmin, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i);

void define_halo (char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i);

void define_and_write_halo(char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i);
void define_and_write_halo_path(char* path, char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i);
void read_halo(char* path);
int access_check();
int access_check_time();
//double time_ve(double ve, double ve0, double t0, double T, int t);
void define_and_write_halo_time(char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i, double ve0, double t0, double T);
double time_ve(double ve, double ve0, double t0, double T, int t);
