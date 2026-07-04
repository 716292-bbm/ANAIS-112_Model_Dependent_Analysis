#pragma once

double counts_bin( char* halo_path, double rhochi, void * input_difcros, double exposure, double E1, double E2, char* model);
double counts_eff_bin( double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model);
double counts_eff_bin_python( char* halo_path, double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target);
double counts_bin_python( char* halo_path, double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target, char*basis);
double counts_bin_smeared( char* halo_path, double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model, double res);

double counts_bin_noread(       double rhochi, void * input_difcros, double exposure, double E1, double E2, char* model);
double counts_bin_python_noread(double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target);


double counts_bin_w( double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
