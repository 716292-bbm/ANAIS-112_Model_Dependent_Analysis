/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
difRateGen.c 
Calculates differential rates
early versions of this code were developed by M. Peiró
and E. Gerstmayr. 
##########################################*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "difRateGen.h"
#include "difCrosSec.h"
#include "halo.h"
#include "phys_consts.h"
#include "coeffs_eft.h"
#include "effFormFact.h"

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef C_NORM
#define C_NORM 1816.9117020201415
#endif


#define VERBOSE 0

/*########################################################################################
Combining results from the velocity distribution integrals and the differential cross section
The return value is given in events/(keV*kg*day)
########################################################################################*/

/*#####################################################################################
First we define the spin independent differential rate in the canonical way
#####################################################################################*/

double difrate_SI_standard(double rho, int A, int Z, double fp, double fn, double mchi, double Er){
	double mN = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mN, mchi);

	double velmin = c*sqrt(1. / (2.*mN*Er*1e-6))*mN*Er*1e-6 / muN;
	double rate = (2.*rho / (mchi*M_PI))*pow((Z*fp + (A - Z)*fn), 2.)*pow(fsi(A, Er), 2.)*halo(velmin, 0)*4.36e+5;
	//printf("rate %.5E \n", rate);
	return rate;
}  

/*#####################################################################################
Below are the functions for the EFT. 
#####################################################################################*/

double difrate_isotope_v0_dEr(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;

	double Er = (val_difcros->Er);
	double mchi = (val_difcros->mchi);
	double jchi = (val_difcros->jchi);
	double v = (val_difcros->v);
	char * Nucleon = (val_difcros->Nucleon);
        double pbGeVfactor = 2.67e-9;
	v = c; //so effectively v=1


	double mtarget = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mchi, mtarget);

	double vmin = c*sqrt((mtarget*Er*1.e-6) / 2.) / muN;
	//printf("halo %.5e\n", halo(vmin, 0));
	return (rhochi / mchi)*difcros_isotope_v0_dEr(A, Z, Er, mchi, jchi, v, Nucleon, F_i, F_j)*halo(vmin, 0)*4.36e+5*(1./pbGeVfactor);
	//return (rhochi / mchi)*difcros_isotope_v0_dEr(A, Z, Er, mchi, jchi, v, Nucleon, F_i, F_j)*shm_halo_analitic(vmin,232,544,220, 0)*4.36e+5*(1./pbGeVfactor);
	
	// shm_halo_analitic(double velmin, double ve, double vesc, double v0, double beta)
}

double difrate_isotope_v0_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs){
                                                                                                    
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;

	double Er = (val_difcros->Er);
	double mchi = (val_difcros->mchi);
	double jchi = (val_difcros->jchi);
	double v = (val_difcros->v);
	char * Nucleon = (val_difcros->Nucleon);
        double pbGeVfactor = 2.67e-9;
	v = c; //so effectively v=1


	double mtarget = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mchi, mtarget);

	double vmin = c*sqrt((mtarget*Er*1.e-6) / 2.) / muN;
        /*printf("vmin: %f\n", vmin);*/
        /*printf("halo %.5e\n", halo(vmin, 0));*/
	return (rhochi / mchi)*difcros_isotope_v0_dEr(A, Z, Er, mchi, jchi, v, Nucleon, F_i, F_j)*halo_w(vmin, ga, gs)*4.36e+5*(1./pbGeVfactor);
}

double difrate_isotope_v2_dEr(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;

	double Er = (val_difcros->Er);
	double mchi = (val_difcros->mchi);
	double jchi = (val_difcros->jchi);
	double v = (val_difcros->v);
	char * Nucleon = (val_difcros->Nucleon);
        const double pbGeVfactor = 2.67e-9;
	v = c; //so effectively v=1


	double mtarget = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mchi, mtarget);

	double vmin = c*sqrt((mtarget*Er*1.e-6) / 2.) / muN;

	//return (rhochi / mchi)*difcros_33_v0_isotope_dEr(A, Z, Er, mchi, jchi, v, Nucleon, Cp, Cn)*shm_halo_analitic(vmin, ve, vesc, v0, beta)*4.36e+5;
	/*return (rhochi / mchi)*difcros_isotope_v2_dEr(A, Z, Er, mchi, jchi, v, Nucleon, F_i, F_j)*halo(vmin, 2)*(1/pbGeVfactor)*4.36e+5/pow(c,2.);*/
	return (rhochi / mchi)*difcros_isotope_v2_dEr(A, Z, Er, mchi, jchi, v, Nucleon, F_i, F_j)*halo(vmin, 2) * C_NORM;
}

double difrate_isotope_v2_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs){
  struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
  double Er = (val_difcros->Er);
  double mchi = (val_difcros->mchi);
  double jchi = (val_difcros->jchi);
  double v = (val_difcros->v);
  char * Nucleon = (val_difcros->Nucleon);
  double pbGeVfactor = 2.67e-9;
  v = c; //so effectively v=1
  double mtarget = approx_mass_nucleus(A, Z);
  double muN = reduced_mass(mchi, mtarget);
  double vmin = c*sqrt((mtarget*Er*1.e-6) / 2.) / muN;
  return (rhochi / mchi)*difcros_isotope_v2_dEr(A, Z, Er, mchi, jchi, v, Nucleon, F_i, F_j) * halo_w(vmin, ga, gs) * C_NORM;
}


double difrate_isotope_dEr(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;

	double mchi = (val_difcros->mchi);

	//return pow(4.*mchi*mtarget,2.)*(difrate_isotope_v0_dEr(A, Z, rhochi, input_difcros, F_i, F_j) + difrate_isotope_v2_dEr(A, Z, rhochi, input_difcros, F_i, F_j));
    return (16*mchi*mchi)*(1./(higgs_vev*higgs_vev*higgs_vev*higgs_vev))*(difrate_isotope_v0_dEr(A, Z, rhochi, input_difcros, F_i, F_j) + difrate_isotope_v2_dEr(A, Z, rhochi, input_difcros, F_i, F_j));

}

double difrate_isotope_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga0, gsl_spline * gs0, gsl_interp_accel *ga2, gsl_spline * gs2)
{
  struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
  double mchi = (val_difcros->mchi);
  return (16*mchi*mchi)*(1./(higgs_vev*higgs_vev*higgs_vev*higgs_vev))*(difrate_isotope_v0_dEr_w(A, Z, rhochi, input_difcros, F_i, F_j, ga0, gs0) + difrate_isotope_v2_dEr_w(A, Z, rhochi, input_difcros, F_i, F_j, ga2, gs2));

}

double total_difrate_isotope_dEr(int A, int Z, double rhochi, void * input_difcros){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	double rate = 0.0;
	int i;
	for ( i = 1; i < 16; i++) // FIXME what about Cp(0) et al?
        {
          if (Cp(i) != 0. || Cn(i) != 0.)
          {
		//printf("CN %.5E CP %.5E \r\n", Cp(i), Cn(i));
		//printf("Op number %i \n", i);
		rate += difrate_isotope_dEr(A, Z, rhochi, val_difcros, i, i);
		//printf("here %.5E \n", rate);
          }
        }

	//INTEREFERENCE TERMS
	if ((Cp(1) != 0. || Cn(1) != 0.) & (Cp(3) != 0. || Cn(3) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 1, 3) + difrate_isotope_dEr(A, Z, rhochi, input_difcros, 3, 1);
	}

	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(5) != 0. || Cn(5) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 4, 5) + difrate_isotope_dEr(A, Z, rhochi, input_difcros, 5, 4);
	}

	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(6) != 0. || Cn(6) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 4, 6) + difrate_isotope_dEr(A, Z, rhochi, input_difcros, 6, 4);
	}
	if ((Cp(8) != 0. || Cn(8) != 0.) & (Cp(9) != 0. || Cn(9) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 8, 9) + difrate_isotope_dEr(A, Z, rhochi, input_difcros, 9, 8);

	}

        // FIXME fabs?
	if (rate >= 0.0){

		return rate;
        }
	if (rate < 0.0){
		return -rate;
	}
 //printf("here %.5E \n", rate);
}

double total_difrate_isotope_dEr_w(int A, int Z, double rhochi, void * input_difcros, gsl_interp_accel *ga0, gsl_spline * gs0, gsl_interp_accel *ga2, gsl_spline * gs2){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;

	double rate = 0.0;
	
        int i;
	for ( i = 1; i < 16; i++)
        {
          if (Cp(i) != 0. || Cn(i) != 0.)
          {
		//printf("CN %.5E CP %.5E \r\n", Cp(i), Cn(i));
		//printf("Op number %i \n", i);
		rate += difrate_isotope_dEr_w(A, Z, rhochi, val_difcros, i, i, ga0,gs0,ga2,gs2);
		//printf("here %.5E \n", rate);
          }
	}

	//INTEREFERENCE TERMS
	if ((Cp(1) != 0. || Cn(1) != 0.) & (Cp(3) != 0. || Cn(3) != 0.)){
		rate += difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 1, 3, ga0,gs0,ga2,gs2) + difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 3, 1, ga0,gs0,ga2,gs2);
	}

	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(5) != 0. || Cn(5) != 0.)){
		rate += difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 4, 5, ga0,gs0,ga2,gs2) + difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 5, 4, ga0,gs0,ga2,gs2);
	}

	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(6) != 0. || Cn(6) != 0.)){
		rate += difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 4, 6, ga0,gs0,ga2,gs2) + difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 6, 4, ga0,gs0,ga2,gs2);
	}
	if ((Cp(8) != 0. || Cn(8) != 0.) & (Cp(9) != 0. || Cn(9) != 0.)){
		rate += difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 8, 9, ga0,gs0,ga2,gs2) + difrate_isotope_dEr_w(A, Z, rhochi, input_difcros, 9, 8, ga0,gs0,ga2,gs2);

	}

        return fabs(rate);
 //printf("here %.5E \n", rate);
}




double difrate_dER(double rhochi, void * input_difcros, double logenergy, char* model){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double znumarr[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	char result[256];
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 74;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Xe[l];	
			atomic_numbers[l]=atomic_numbers_Xe[l];	
		}
	}
	if (strncmp(target, "Ge", 10) == 0){
		num_isos = 5;
		znum = 32; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ge[l];	
			atomic_numbers[l]=atomic_numbers_Ge[l];	
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
		num_isos = 1;
		znum = 18; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ar[l];	
			atomic_numbers[l]=atomic_numbers_Ar[l];	
		}
	}
	if (strncmp(target, "F", 10) == 0){
		num_isos = 1;
		znum = 9; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_F[l];	
			atomic_numbers[l]=atomic_numbers_F[l];	
		}
	}
	if (strncmp(target, "CaWO4", 10) == 0){
		num_isos = 2;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=Z_numbers_CaWO4[l];
			isotopes[l]=isotopes_CaWO4[l];	
			atomic_numbers[l]=atomic_numbers_CaWO4[l];
			printf("Z %i A %i number %f\n", Z_numbers_CaWO4[l],atomic_numbers[l],isotopes[l]);

		}
	}
	if (strncmp(target, "Na", 10) == 0){
		num_isos = 1;
		znum = 11; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=1.0;
			atomic_numbers[l]=23;
		}
	}

		if (strncmp(target, "I", 10) == 0){
		num_isos = 1;
		znum = 53; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=1.0;
			atomic_numbers[l]=127;
		}
	}
	double counts=0.0;
	double energy = pow(10., logenergy);
	val_difcros->Er = energy;
	for (int l = 0; l< num_isos; l++){
		//printf("check %i, %lf , %i, %lf\n", l, isotopes[l], atomic_numbers[l],total_difrate_isotope_dEr(atomic_numbers[l],znum, rhochi, val_difcros));
		if (strncmp(model, "Light_Med", 10)==0)
		{
			double med_mass = give_med_mass();
			printf("Here we are, mass =%lf\n", med_mass );
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "ElecMom", 10)==0)
		{	
			double C11p = Cp(11);
			set_coeffs();
			set_any_Ncoeff(C11p, 11, "p");
			double med_mass = 0.0;
			printf("Here we are, mass =%lf\n", med_mass );
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "MagMom", 10)==0){
			double C1p, C4p, C4n, C5p, C6p, C6n;
			C1p = Cp(1);
			C4p = Cp(4);
			C4n = Cn(4);
			C5p = Cp(5);
			C6p = Cp(6);
			C6n = Cn(6);
			set_coeffs();
			set_any_Ncoeff(C1p, 1, "p");
			set_any_Ncoeff(C4p, 4, "p");
			set_any_Ncoeff(C4n, 4, "n");
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
			
			set_coeffs();
			set_any_Ncoeff(C5p, 5, "p");
			set_any_Ncoeff(C6p, 6, "p");
			set_any_Ncoeff(C6n, 6, "n");
			double med_mass = 0.0;
			//printf("Here we are, mass =%lf\n", med_mass );
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else {
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
			//printf("%.5E \n", counts);
		}
	}
	return counts;
}


double difrate_dER_w(double rhochi, void * input_difcros, double linenergy, char* model, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double znumarr[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	/*char result[256];*/
        num_isos = 7;
        znum = 74;
        for (int l = 0; l < num_isos; l++){
            znumarr[l]=znum;
            isotopes[l]=isotopes_Xe[l];	
            atomic_numbers[l]=atomic_numbers_Xe[l];	
        }
	double counts=0.0;
	/*double energy = pow(10., logenergy);*/
	val_difcros->Er = linenergy; // difrate_integrand already sets this FIXME
	for (int l = 0; l< num_isos; l++)
        {
            counts += isotopes[l] * total_difrate_isotope_dEr_w(atomic_numbers[l],znumarr[l], rhochi, val_difcros, ga, gs, ga2, gs2);
	}
	return counts;
}



double difrate_dER_python(double rhochi, double mass, double logenergy, char* model, char* Target, char*ISO_switch){
  struct difcros_params struct_difcros_T= {Target, 1., mass, 0.5, 220, ISO_switch, 0.0};
  
  double rate = difrate_dER(rhochi, &struct_difcros_T, logenergy, model);
  return rate;
}

double difrate_dER_python_2(double rhochi, double mass, double energy, char* model, char* Target, char*ISO_switch){
  struct difcros_params struct_difcros_T= {Target, 1., mass, 0.5, 220, ISO_switch, 0.0};
  
  double rate = difrate_dER_2(rhochi, &struct_difcros_T, energy, model);
  return rate;
}

double difrate_dER_2(double rhochi, void * input_difcros, double energy, char* model){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double znumarr[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	char result[256];
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 74;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Xe[l];	
			atomic_numbers[l]=atomic_numbers_Xe[l];	
		}
	}
	if (strncmp(target, "Ge", 10) == 0){
		num_isos = 5;
		znum = 32; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ge[l];	
			atomic_numbers[l]=atomic_numbers_Ge[l];	
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
		num_isos = 1;
		znum = 18; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ar[l];	
			atomic_numbers[l]=atomic_numbers_Ar[l];	
		}
	}
	if (strncmp(target, "F", 10) == 0){
		num_isos = 1;
		znum = 9; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_F[l];	
			atomic_numbers[l]=atomic_numbers_F[l];	
		}
	}
	if (strncmp(target, "CaWO4", 10) == 0){
		num_isos = 2;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=Z_numbers_CaWO4[l];
			isotopes[l]=isotopes_CaWO4[l];	
			atomic_numbers[l]=atomic_numbers_CaWO4[l];
			printf("Z %i A %i number %f\n", Z_numbers_CaWO4[l],atomic_numbers[l],isotopes[l]);

		}
	}
	double counts=0.0;
	//double energy = pow(10., logenergy);
	val_difcros->Er = energy;
	for (int l = 0; l< num_isos; l++){
		//printf("check %i, %lf , %i, %lf\n", l, isotopes[l], atomic_numbers[l],total_difrate_isotope_dEr(atomic_numbers[l],znum, rhochi, val_difcros));
		if (strncmp(model, "Light_Med", 10)==0)
		{
			double med_mass = give_med_mass();
			printf("Here we are, mass =%lf\n", med_mass );
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "ElecMom", 10)==0)
		{	
			double C11p = Cp(11);
			set_coeffs();
			set_any_Ncoeff(C11p, 11, "p");
			double med_mass = 0.0;
			printf("Here we are, mass =%lf\n", med_mass );
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "MagMom", 10)==0){
			double C1p, C4p, C4n, C5p, C6p, C6n;
			C1p = Cp(1);
			C4p = Cp(4);
			C4n = Cn(4);
			C5p = Cp(5);
			C6p = Cp(6);
			C6n = Cn(6);
			set_coeffs();
			set_any_Ncoeff(C1p, 1, "p");
			set_any_Ncoeff(C4p, 4, "p");
			set_any_Ncoeff(C4n, 4, "n");
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
			
			set_coeffs();
			set_any_Ncoeff(C5p, 5, "p");
			set_any_Ncoeff(C6p, 6, "p");
			set_any_Ncoeff(C6n, 6, "n");
			double med_mass = 0.0;
			//printf("Here we are, mass =%lf\n", med_mass );
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else {
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
			//printf("%.5E \n", counts);
		}
	}
	return counts;
}
