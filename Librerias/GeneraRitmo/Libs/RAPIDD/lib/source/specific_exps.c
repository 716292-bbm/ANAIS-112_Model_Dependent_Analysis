/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             

##########################################*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include <string.h>
#include <unistd.h> 

#include "../source/difRateGen.h"
#include "../source/coeffs_eft.h"
#include "../source/phys_consts.h"
#include "../source/binning_general.h"
#include "efficiency_curve.h"
#include<gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>



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

int leff_length;
float leff_data[1000];
float leff_er[1000];
char leff_path_true[26];

struct binning_params      { double rhochi; double exposure; void * input_difcros; double E1; double E2; char* model;};


// change the name!
double bin_Xenon1T( double rhochi, double mass, double E1, double E2, char* model, char* basis){
  double exposure = 0.9*1e3 * 278.8 * 0.475;

  struct difcros_params struct_difcros_Xe = {"Xe", 1., mass, 0.5, 220, basis, exposure};
  
  double integrated = counts_eff_bin( rhochi, & struct_difcros_Xe, exposure, E1, E2, model);
  
  return integrated;

}

double bin_DS50( double rhochi, double mass, double E1, double E2, char* model, char* basis){
  double exposure = 16660;

  struct difcros_params struct_difcros_Ar = {"Ar", 1., mass, 0.5, 220, basis, exposure};
  
  double integrated = counts_eff_bin( rhochi, & struct_difcros_Ar, exposure, E1, E2, model);
  
  return integrated;
}

//Linhard emperically determined for lux  k = 0.1735 in https://arxiv.org/pdf/1712.05696.pdf

double lindhard_lux(int A, int Z, double ER){
  double k = 0.1735;
  double eps = 11.5*pow(Z,-7.0/3.0)*ER ;
  double g = 3.0*pow(eps, 0.15) + 0.7*pow(eps, 0.6) + eps;
  double L = k*g/(1+k*g);
  return L;
}

//Xenon1T and XENONnT

double res_fn_xenon(double ERee){
  double a = 0.310;
  double b = 0.0037;
  return a*sqrt(ERee) + b * ERee; 
}

double res_fn_xenon_nr(int A, int Z, double ERnr){
  double ERee = ERnr*lindhard_lux(A, Z, ERnr);
  double a = 0.310;
  double b = 0.0037;
  return a/sqrt(ERee) + b ; 
}

double difrate_effres_integrand_Xenon1T(double x, void * p){

        /*printf("INTEGRAND\n");*/
	struct binning_params * params = (struct binning_params *)p;
	/*double isotopes[10];*/
	/*double prefact[10];*/
	/*int atomic_numbers[10]; */
	/*int num_isos, znum;*/
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	double E1 = (params->E1);
	double E2 = (params->E2);
  double res = res_fn_xenon_nr(132, 74, x)*x;

	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
  char* model = (params -> model);
  double p_smear = 0.5*(erf((E2-x)/(sqrt(2)*res)) - erf((E1-x)/(sqrt(2)*res)));
  //printf("p_smear %f res %f  Er %f \n", p_smear, res, x);
  double rate = exposure*p_smear*efficiency(x)*difrate_dER(rhochi, val_difcros, log10(x), model);

	return rate;
}

double counts_effres_bin_Xenon1T( double rhochi, double mass, double E1, double E2, char*model, char*basis){
	double result, abserr;
	struct difcros_params struct_difcros_Xe = {"Xe", 1., mass, 0.5, 220, basis, 0.0};
	double exposure = 0.9*1e3 * 278.8 * 0.475;
	gsl_function F;
	struct binning_params params = {rhochi, exposure, &struct_difcros_Xe, E1, E2, model};
	
	F.function = &difrate_effres_integrand_Xenon1T;
	F.params = &params; 
  // Now integral is not over the bins, its over whole range, for now lets put Eth 1-50 keV
  // Also set Emin as 2.0 so you can't have counts under 2.0 keV 
  double Emin = 2.0;
  if(E2 <= Emin){
    result = 0.0;
  }
  else{
    double Eth = 0.1;
    double Emax = 50.0;
    size_t nevals;
    gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
    gsl_integration_cquad (&F, Eth, Emax, 1.e-4, 1.e-4, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
  }
	
  return result;
}

// LUX and LZ related 

double k_lindhard(int A, int Z){
  return 0.133*pow(Z, 2.0/3.0)*(1.0/sqrt(A));

}


double res_fn_lux(double ERee){
  double a = 0.33; 
  return a/sqrt(ERee);
}

double res_fn_lux_nr(int A, int Z, double ERnr){
  double ERee = ERnr*lindhard_lux(A, Z, ERnr);
  double a = 0.33;
  return a/sqrt(ERee); 
}

double bin_LZ( double rhochi, double mass, double E1, double E2, char* model, char* basis){
  double exposure = 5.6*1.0e3 * 1.0e3;

  struct difcros_params struct_difcros_Xe = {"Xe", 1., mass, 0.5, 220, basis, exposure};
  
  double integrated = counts_eff_bin( rhochi, & struct_difcros_Xe, exposure, E1, E2, model);
  
  return integrated;

}

/// TESTING resolution effects 

double difrate_effres_integrand_LZ(double x, void * p){

        /*printf("INTEGRAND\n");*/
	struct binning_params * params = (struct binning_params *)p;
	/*double isotopes[10];*/
	/*double prefact[10];*/
	/*int atomic_numbers[10]; */
	/*int num_isos, znum;*/
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	double E1 = (params->E1);
	double E2 = (params->E2);
  double res = res_fn_lux_nr(132, 74, x)*x;

	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
  char* model = (params -> model);
  double p_smear = 0.5*(erf((E2-x)/(sqrt(2)*res)) - erf((E1-x)/(sqrt(2)*res)));
  //printf("p_smear %f res %f  Er %f \n", p_smear, res, x);
  double rate = exposure*p_smear*efficiency(x)*difrate_dER(rhochi, val_difcros, log10(x), model);

	return rate;
}

double counts_effres_bin_LZ( double rhochi, double mass, double E1, double E2, char*model, char*basis){
	double result, abserr;
	struct difcros_params struct_difcros_Xe = {"Xe", 1., mass, 0.5, 220, basis, 0.0};
	double exposure = 5.6*1.0e3 * 1.0e3;
	gsl_function F;
	struct binning_params params = {rhochi, exposure, &struct_difcros_Xe, E1, E2, model};
	
	F.function = &difrate_effres_integrand_LZ;
	F.params = &params; 
  // Now integral is not over the bins, its over whole range, for now lets put Eth 1-50 keV
  double Emin = 1.0;
  if(E2<=Emin){
    result=0.0;
  }
  else{
    double Eth = 0.8;
    double Emax = 50.0;
    size_t nevals;
    gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
    gsl_integration_cquad (&F, Eth, Emax, 1.e-2, 1.e-2, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
  }
	
  return result;
}



double bin_DS20k( double rhochi, double mass, double E1, double E2, char* model, char* basis){
  double exposure = 200*1.0e3 * 365;

  struct difcros_params struct_difcros_Ar = {"Ar", 1., mass, 0.5, 220, basis, exposure};
  
  double integrated = counts_eff_bin( rhochi, & struct_difcros_Ar, exposure, E1, E2, model);
  
  return integrated;

}

double bin_DEAP3600( double rhochi, double mass, double E1, double E2, char* model, char* basis){
  double exposure = 758*1e3 ; //200*1.0e3 * 365;

  struct difcros_params struct_difcros_Ar = {"Ar", 1., mass, 0.5, 220, basis, exposure};
  
  double integrated = counts_eff_bin( rhochi, & struct_difcros_Ar, exposure, E1, E2, model);
  
  return integrated;

}

// TODO!!!!! 


//======================================================================
//  Leff in fig21 of 
// https://iopscience.iop.org/article/10.1088/1748-0221/12/10/P10015/pdf
//======================================================================

void read_DS50_LEFF(char* path){
  FILE * table;
  int i;
  int j;
  int check_length;
  leff_length = countlines(path);
  //printf("%s\n", path );
  if (access(path, F_OK) == -1){
    //printf("Incorrect path to efficiency, taking as 1\n");
    snprintf(leff_path_true, sizeof(leff_path_true), "False");
    //printf("%s\n", leff_path_true );
  }
  else{
    snprintf(leff_path_true, sizeof(leff_path_true), "True");
    //printf("Reading pre-calculated table from file %s...\n", path);
    table = fopen(path, "r");

    //printf("Getting table format...\n");
    //printf("Done! %i\n", leff_length);

    //printf("Reading data table...\n");
      for (j = 0; j < leff_length; j++){
        fscanf(table, "%E %E", &leff_er[j], &leff_data[j]);

      }

      fclose(table);
      //printf("printing here %i\n",leff_length );

  //printf("Done!\n");
}
}

double DS50_LEFF ( double E_r){
  if ( strncmp(leff_path_true, "True", 10) == 0){
    int m;
    double x[leff_length], y[leff_length];
    double result;

    for (m = 0; m < leff_length; m++){
      x[m] = leff_er[m];
      y[m] = leff_data[m];
    }
    if (E_r >= x[0] && E_r <= x[leff_length-1]){
        /*printf("%f\n", x[length-1]);*/
        gsl_interp_accel *leff = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, leff_length);

        gsl_spline_init (spline, x, y, leff_length);

        result = gsl_spline_eval (spline,E_r, leff);

        gsl_spline_free (spline);
        gsl_interp_accel_free (leff);
        if (result < 0.0) {
          return 0.0;
        }
        else{
          return result;
        }
      }
      if (E_r < x[0]) {  // eff(E_r) data out of the range
        return y[0];
      }
      if (E_r > x[leff_length-1]) {  // eff(E_r) data out of the range
        return y[leff_length-1];
      }

    }
  else{
    return 0.0;
  }
}

double lindhard_DS(int A, int Z, double ER){
  double k = k_lindhard(A, Z);
  double eps = 11.5*pow(Z,-7.0/3.0)*ER ;
  double g = 3.0*pow(eps, 0.15) + 0.7*pow(eps, 0.6) + eps;
  double Leff = k*g/(1+k*g);
  return Leff;
}

double res_fn_DS20(double ERee){
  double a = 0.452; 
  return a/sqrt(ERee);
}

double res_fn_DS20_nr(int A, int Z, double ERnr){
  double ERee = ERnr*DS50_LEFF(ERnr);
  //double ERee = ERnr*lindhard_DS(A, Z, ERnr);
  double a = 0.452;
  return a/sqrt(ERee); 
}


/// TESTING resolution effects 

double difrate_effres_integrand_DS20(double x, void * p){

        /*printf("INTEGRAND\n");*/
	struct binning_params * params = (struct binning_params *)p;
	/*double isotopes[10];*/
	/*double prefact[10];*/
	/*int atomic_numbers[10]; */
	/*int num_isos, znum;*/
	double rhochi = (params->rhochi);
	double exposure = (params->exposure);
	double E1 = (params->E1);
	double E2 = (params->E2);
  double res = res_fn_DS20_nr(40, 18, x)*x;

	struct difcros_params * val_difcros = (params-> input_difcros);
	val_difcros->Er = x; 
  char* model = (params -> model);
  double p_smear = 0.5*(erf((E2-x)/(sqrt(2)*res)) - erf((E1-x)/(sqrt(2)*res)));
  //printf("p_smear %f res %f  Er %f \n", p_smear, res, x);
  double rate = exposure*p_smear*efficiency(x)*difrate_dER(rhochi, val_difcros, log10(x), model);

	return rate;
}

double counts_effres_bin_DS50( double rhochi, double mass, double E1, double E2, char*model, char*basis){
	double result, abserr;
  double exposure = 16660;

  struct difcros_params struct_difcros_Ar = {"Ar", 1., mass, 0.5, 220, basis, exposure};
  
	gsl_function F;
	struct binning_params params = {rhochi, exposure, &struct_difcros_Ar, E1, E2, model};
	
	F.function = &difrate_effres_integrand_DS20;
	F.params = &params; 

  double Emin = 47.0;
  if(E2<=Emin){
    result=0.0;
  }
  else{
  // Now integral is not over the bins, its over whole range, for now lets put Eth 1-50 keV
    double Eth = 45.0;
    double Emax = 200.0;
    size_t nevals;
    gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
    gsl_integration_cquad (&F, Eth, Emax, 1.e-4, 1.e-4, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
  }
	
  return result;
}

double counts_effres_bin_DS20k( double rhochi, double mass, double E1, double E2, char*model, char*basis){
	double result, abserr;
  double exposure = 200*1.0e3 * 365;

  struct difcros_params struct_difcros_Ar = {"Ar", 1., mass, 0.5, 220, basis, exposure};
  
	gsl_function F;
	struct binning_params params = {rhochi, exposure, &struct_difcros_Ar, E1, E2, model};
	
	F.function = &difrate_effres_integrand_DS20;
	F.params = &params; 

  double Emin = 20.0;
  if(E1<Emin){
    result=0.0;
  }
  else{
  // Now integral is not over the bins, its over whole range, for now lets put Eth 1-50 keV
    double Eth = 17.0;
    double Emax = 200.0;
    size_t nevals;
    gsl_integration_cquad_workspace * v = gsl_integration_cquad_workspace_alloc (100); 
    gsl_integration_cquad (&F, Eth, Emax, 1.e-4, 1.e-4, v, &result, &abserr, &nevals);
    gsl_integration_cquad_workspace_free(v);
  }
	
  return result;
}
//====================================================================
//    Modified from eq 5.15 of 
//  https://iopscience.iop.org/article/10.1088/1748-0221/12/10/P10015 
//====================================================================



// double res_fn_DS_nr(int A, int Z, double ERnr){
//   double ERee = ERnr*lindhard_DS50(A, Z, ERnr);
//   double a = 0.33;
//   return a/sqrt(ERee); 
// }
