/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
Binning_general.c 
##########################################*/

#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include "../source/difRateGen.h"
#include "../source/halo.h"
#include "../source/coeffs_eft.h"
#include "../source/phys_consts.h"
#include "efficiency_curve.h"

#include "binning_general.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))

 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

struct binning_params      { double rhochi; double exposure; void * input_difcros; double E1; double E2; char* model;};
struct binning_params_w      { double rhochi; double exposure; void * input_difcros; double E1; double E2; char* model; gsl_interp_accel *ga; gsl_spline * gs; gsl_interp_accel *ga2; gsl_spline * gs2;};

struct binning_eff_params_w      { double rhochi; double exposure; void * input_difcros; double E1; double E2; char* model; gsl_interp_accel *ga; gsl_spline * gs; gsl_interp_accel *ga2; gsl_spline * gs2; gsl_interp_accel *gaeff; gsl_spline * gseff;};

struct binning_params_smear { double rhochi; double exposure; void * input_difcros; double E1; double E2; char* model; double res;};


double difrate_integrand(double x, void * p){

	struct binning_params * params = (struct binning_params *)p;
	double isotopes[10];
	double prefact[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	double E1 = (params->E1);
	double E2 = (params->E2);
	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
    double logenergy = log10(x);
    char* model = (params -> model);

    double rate = exposure*difrate_dER(rhochi, val_difcros, log10(x), model);

	return rate;
}

double difrate_integrand_w(double x, void * p){

        /*printf("INTEGRAND\n");*/
	struct binning_params_w * params = (struct binning_params_w *)p;
	/*double isotopes[10];*/
	/*double prefact[10];*/
	/*int atomic_numbers[10]; */
	/*int num_isos, znum;*/
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	/*double E1 = (params->E1);*/
	/*double E2 = (params->E2);*/
        gsl_interp_accel * ga = (params->ga);
        gsl_spline       * gs = (params->gs);
        gsl_interp_accel * ga2 = (params->ga2);
        gsl_spline       * gs2 = (params->gs2);
	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
    char* model = (params -> model);

    double rate = exposure*difrate_dER_w(rhochi, val_difcros, x, model, ga, gs, ga2, gs2);

	return rate;
}

double difrate_eff_integrand(double x, void * p){

        /*printf("INTEGRAND\n");*/
	struct binning_params * params = (struct binning_params *)p;
	/*double isotopes[10];*/
	/*double prefact[10];*/
	/*int atomic_numbers[10]; */
	/*int num_isos, znum;*/
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	/*double E1 = (params->E1);*/
	/*double E2 = (params->E2);*/

	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
    char* model = (params -> model);

    double rate = exposure*efficiency(x)*difrate_dER(rhochi, val_difcros, log10(x), model);

	return rate;
}


double counts_bin( char* halo_path, double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model){

	double result, abserr;
        read_halo(halo_path);
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	
	gsl_function F;
	struct binning_params params = {rhochi, exposure, val_difcros, E1, E2, model};
	
	F.function = &difrate_integrand;
	F.params = &params; 
	double Eleft = E1;
	double Eright = E2;
	size_t nevals;
	gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
	gsl_integration_cquad (&F, Eleft, Eright, 1.e-10, 1.e-10, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
	
    return result;
}


double counts_eff_bin( double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model){
	double result, abserr;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	
	gsl_function F;
	struct binning_params params = {rhochi, exposure, val_difcros, E1, E2, model};
	
	F.function = &difrate_eff_integrand;
	F.params = &params; 
	double Eleft = E1;
	double Eright = E2;
	size_t nevals;
	gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
	gsl_integration_cquad (&F, Eleft, Eright, 1.e-1, 1.e-1, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
	
    return result;
}


double counts_bin_w(double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2){

	double result, abserr;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	
	gsl_function F;
	struct binning_params_w params = {rhochi, exposure, val_difcros, E1, E2, model, ga, gs, ga2, gs2};
	
	F.function = &difrate_integrand_w;
	F.params = &params; 
	double Eleft = E1;
	double Eright = E2;
	size_t nevals;
	gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
	gsl_integration_cquad (&F, Eleft, Eright, 1.e-1, 1.e-1, v, &result, &abserr, &nevals);
        gsl_integration_cquad_workspace_free(v);
	
    return result;
}

double counts_bin_python( char* halo_path, double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target, char*basis){
	struct difcros_params struct_difcros_Xe = {Target, 1., mass, 0.5, 220, basis, 0.0};
	double counts = counts_bin( halo_path, rhochi, &struct_difcros_Xe, exposure, E1,  E2, model);
	return counts ;
}

double counts_eff_bin_python( char* halo_path, double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target){
	struct difcros_params struct_difcros_Xe = {Target, 1., mass, 0.5, 220, "All", 0.0};
	double counts = counts_eff_bin( rhochi, &struct_difcros_Xe, exposure, E1,  E2, model);
	return counts ;
}

double counts_bin_noread(double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model)
{
  double result, abserr;
  struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
  gsl_function F;
  struct binning_params params = {rhochi, exposure, val_difcros, E1, E2, model};

  F.function = &difrate_integrand;
  F.params = &params; 
  double Eleft = E1;
  double Eright = E2;
  size_t nevals;
  gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
  gsl_integration_cquad (&F, Eleft, Eright, 1.e-1, 1.e-1, v, &result, &abserr, &nevals);
  gsl_integration_cquad_workspace_free(v);
	
   return result;
}

double counts_bin_python_noread(double rhochi, double mass, double exposure, double E1, double E2, char*model, char*Target)
{
	struct difcros_params struct_difcros_Xe = {Target, 1., mass, 0.5, 220, "All", 0.0};
	double counts = counts_bin_noread(rhochi, &struct_difcros_Xe, exposure, E1,  E2, model);
	return counts ;
}

double difrate_integrand_smearing(double x, void * p){

	struct binning_params_smear * params = (struct binning_params_smear *)p;
	double isotopes[10];
	double prefact[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	double E1 = (params->E1);
	double E2 = (params->E2);
	double res = (params->res);
	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
    double logenergy = log10(x);
    char* model = (params -> model);
	double rate;
	if ( x > 0.1){
		double p_smear = 0.5*(erfc(E2-x)/(sqrt(2)*res) - erfc(E1-x)/(sqrt(2)*res));
    	rate = exposure*difrate_dER(rhochi, val_difcros, log10(x), model);
	}
	else
	{
		rate = 0.0;
	}
	
	
	
	return rate;
}

double counts_bin_smeared( char* halo_path, double rhochi, void * input_difcros, double exposure, double E1, double E2, char*model, double res){

	double result, abserr;
	read_halo(halo_path);
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	
	gsl_function F;
	struct binning_params_smear params = {rhochi, exposure, val_difcros, E1, E2, model, res};
	
	F.function = &difrate_integrand;
	F.params = &params; 
	double Eleft = E1;
	double Eright = E2;
	size_t nevals;
	gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
	gsl_integration_cquad (&F, Eleft, Eright, 1.e-2, 1.e-2, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
	
    return result;
}
