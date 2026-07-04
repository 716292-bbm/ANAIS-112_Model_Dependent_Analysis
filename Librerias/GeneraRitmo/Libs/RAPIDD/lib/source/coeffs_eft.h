#pragma once

#define OPERATOR_NUMBER 16

double Cp(int i);
double Cn(int i);
double Cp_trace(int i);
double Cn_trace(int i);
double give_mass();
double give_med_mass();
double give_rho();
double give_vesc();
double give_v0();
double give_k();
double give_vc();
void set_coeffs();
void print_coeffs();
void print_SI_coeffs();
void read_halo_params(const char* fname);

//void read_coeffs();
void read_coeffs(const char* fname);
void read_coeffs_low_MED(const char* fname);
void read_coeffs_iso(const char* fname);
void write_coeffs(char* output_path);
void set_coeffs_matching_Standard( double mchi, double sigma, double sigmaSD);
void set_coeffs_matching_Standard_all( double mchi, double sigma_p,double sigma_n, double sigmaSD_p, double sigmaSD_n);
void set_any_coeffs(double C, int i);
void set_any_Ncoeff(double input, int i, char* nucleon);
void set_med_mass(double med_mass);
void set_any_coeffs_trace(double input, int i, char*nucleon);
void set_coeffs_anapole(double mchi, double ANAPOLE);
void set_coeffs_magmom(double mchi, double mu_x);
void set_coeffs_elecmom(double mchi, double mu_x);
void set_coeffs_charrad(double mchi, double b_x);
void coeffs_set_master(char* coeffs, char* model);
void copy_trace();
