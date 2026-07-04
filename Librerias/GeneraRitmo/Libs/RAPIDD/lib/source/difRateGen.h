#pragma once
#include <gsl/gsl_spline.h>

struct difcros_params { char * target; double Er; double mchi; double jchi; double v; char * Nucleon; double exposure; };
struct halo_params { char * profile; double vesc; double v0; double beta; double vt; double vc; double ve; double k; int i; double ve0; double t0; double T; char * halo_path; };

double difrate_SI_standard(double rho, int A, int Z, double fp, double fn, double mchi, double Er);
double total_difrate_isotope_dEr(int A, int Z, double rhochi, void * input_difcros);
double difrate_dER(double rhochi, void * input_difcros, double logenergy, char* model);
double difrate_dER_2(double rhochi, void * input_difcros, double energy, char* model);
double difrate_isotope_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
double difrate_isotope_v2_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs);
double difrate_isotope_v0_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs);

double difrate_dER_w(double rhochi, void * input_difcros, double logenergy, char* model, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
double total_difrate_isotope_dEr_w(int A, int Z, double rhochi, void * input_difcros, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
double difrate_dER_python(double rhochi, double mass, double logenergy, char* model, char* Target, char*ISO_switch);
double difrate_dER_python_2(double rhochi, double mass, double energy, char* model, char* Target, char*ISO_switch);
