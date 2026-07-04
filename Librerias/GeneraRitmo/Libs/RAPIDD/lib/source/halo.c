/*########################################################################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
file : halo.c
Effective operators for DM direct detection 
It uses 2D integrations to change from galactic to Earth (lab) coordinates.

early versions of this code were developed by M. Peiró
and E. Gerstmayr. 

########################################################################################*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>//Contains gamma and beta functions.
#include <unistd.h> //for access check
#include <time.h> //for elapsed time
#include <sys/stat.h> //for mkdir


#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))


#include "halo.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../cubature.h"

#define VERBOSE 0
#define length 100 //Number of divisions between 0 and (vesc+ve) to perform the interpolation.
#define power 10 //Maximum predefined power of the velocity in the halo integral

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

float vel[length];
float fv[power][length];

#define Narray 50
double pQ[Narray+1];
double pQa[Narray+1];


double Normalization_shm;
double Normalization_fornasa;
double Normalization_lisanti;

/*########################################################################################
Standard Halo distribution (cutoff and smooth distributions controlled by beta)
See: 1005.0579 for more info.
In the following the integer "i" will control the power in velocity included in
calculation. For i=0 (usual calculation); i=1 (operator proportional to v);
i=2 (operator proportional v^2) and so on.
########################################################################################*/

double norm_shm(double vesc, double v0, double beta)//SMH function normalization
{
double xesc=vesc/v0;
  return pow(M_PI,1.5)*pow(v0,3.)*(gsl_sf_erf(xesc)-4/sqrt(M_PI)
  *exp(-xesc*xesc)*(xesc/2.+beta*pow(xesc,3.)/3.));
}

struct Nshm_params { double vesc_Nshm; double v0_Nshm; double beta_Nshm;};

double N_shm(double v,void * p) {
    struct Nshm_params * params
    = (struct Nshm_params *)p;
    double vesc_Nshm = (params->vesc_Nshm);
    double v0_Nshm = (params->v0_Nshm);
    double beta_Nshm = (params->beta_Nshm);
    return 4.*M_PI*v*v*(exp(-(v*v)/(v0_Nshm*v0_Nshm))
     - beta_Nshm*exp(-vesc_Nshm*vesc_Nshm/(v0_Nshm*v0_Nshm)));
}

double norm_shm_num(double vesc, double v0, double beta){
    double result, error;
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    struct Nshm_params params = {vesc,v0,beta};

    F.function = &N_shm;
    F.params = &params;

    gsl_integration_qags (&F, 0, vesc, 1e-8, 1e-8, 1000 ,w, &result, &error);
    gsl_integration_workspace_free (w);

    return result;

}



double shm_halo_analitic(double velmin, double ve, double vesc, double v0, double beta)
{
double xesc=vesc/v0;
double xmin=velmin/v0;
double xe=ve/v0;

if((xe+xmin)<xesc){
  return pow(M_PI,1.5)*pow(v0,2.)/(2.*xe*norm_shm(vesc,v0,beta))
  *(gsl_sf_erf(xmin+xe)-gsl_sf_erf(xmin-xe)-4.*xe/sqrt(M_PI)*exp(-xesc*xesc)
  *(1.+beta*(pow(xesc,2.)-pow(xe,2.)/3.-pow(xmin,2.))));
}
if(xmin>fabs(xesc-xe) && (xe+xesc)>xmin){
  return pow(M_PI,1.5)*pow(v0,2.)/(2.*xe*norm_shm(vesc,v0,beta))
  *(gsl_sf_erf(xesc)+gsl_sf_erf(xe-xmin)-2./sqrt(M_PI)*exp(-xesc*xesc)
  *(xesc+xe-xmin-1./3.*beta*(xe-2.*xesc-xmin)*pow(xesc+xe-xmin,2.)));
}
if(xe>(xmin+xesc)){
  return 1./(xe*v0);
}
 if(xmin>(xe+xesc)){
  return 0.;
}
 return 0.;

}

struct shm_params { double vesc_param;
	double v0_param; double ve_param; double beta_param; int i_param;};

int sign(double x){
    if (x > 0) return 1;
    if (x < 0) return -1;

return 0;

}

int f_shm(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){

    struct shm_params * params
     = (struct shm_params *)p;
    double vesc_param = (params->vesc_param);
    double v0_param = (params->v0_param);
    double ve_param = (params->ve_param);
    double beta_param = (params->beta_param);
    int i_param = (params->i_param);


    fval[0] = 2.*M_PI*k[1] * pow(k[1], i_param*1.0)/*(1./N)*/*
    (exp(-(k[1] * k[1] + ve_param*ve_param - 2.*k[1] * ve_param*k[0]) / (v0_param*v0_param))
     - beta_param*exp(-vesc_param*vesc_param/(v0_param*v0_param)));

    /*Edited by Andrew originally 2.*M_PI*k[1] * pow(k[1], i_param*1.0)
    (exp(-(k[1] * k[1] + ve_param*ve_param + 2.*k[1] * ve_param*k[0]) / (v0_param*v0_param))
     - beta_param*exp(-vesc_param*vesc_param/(v0_param*v0_param)))*/

return 0;

}

/* This is added to use the smooth cut-off expression used in arXiv:1112.0524v2 added on 12/02/2016 */

int f_shm_ann(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){

    struct shm_params * params
    = (struct shm_params *)p;
    double vesc_param = (params->vesc_param);
    double v0_param = (params->v0_param);
    double ve_param = (params->ve_param);
    double beta_param = (params->beta_param);
    int i_param = (params->i_param);

    fval[0] = 2.*M_PI*k[1] * pow(k[1], i_param*1.0)/*(1./N)*/*
    (exp((-3/2)*(k[1] * k[1] + ve_param*ve_param - 2.*k[1] * ve_param*k[0]+beta_param*(-vesc_param)*vesc_param)/(v0_param*v0_param)));


    return 0;

}


/* The following is the 2D integration in angle and velocity */

double shm_halo (double vmin, double vesc, double v0, double ve, double beta, int i){

    double val, err;
    struct shm_params params = {vesc,v0,ve,beta,i};

    //limits of integration
    double xl[2]={-1, vmin};
	double xu[2] = {1, vesc + ve};

if(vmin>vesc+ve) return 0.;


else {
    hcubature(1, &f_shm, &params,
	      2, xl, xu,
	      0, 0, 1e-13, ERROR_INDIVIDUAL, &val, &err);

return val;}
}

/*########################################################################################
Extracted phase-space densities from N-body simulations
See: 1010.4300 for more info.
In the following the integer "i" will control the power in velocity included in
calculation. For i=0 (usual calculation); i=1 (operator proportional to v);
i=2 (operator proportional v^2) and so on.

Use equation 15 in 1208.6426
########################################################################################*/

//we calculate the norm of the distribution

//we calculate the norm of the distribution

struct Nk_params { double vesc_Nk; double v0_Nk; double k_Nk;};

double Nk_intd(double v,void * p) {
    struct Nk_params * params
    = (struct Nk_params *)p;
    double vesc_Nk = (params->vesc_Nk);
    double v0_Nk = (params->v0_Nk);
    double k_Nk = (params->k_Nk);
    double x=pow(v/v0_Nk,1);
    double xesc=pow(vesc_Nk/v0_Nk,1);
    if(k_Nk>=0.1) return 4*M_PI*v*v*pow(exp(-(v*v)/(k_Nk*v0_Nk*v0_Nk))-exp(-(vesc_Nk*vesc_Nk)/(v0_Nk*v0_Nk*k_Nk)),k_Nk);
    else return  v*v*exp(-x);
}

double Nk_uncert (double vesc, double v0, double k){
    double result, error;

    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    struct Nk_params params = {vesc,v0,k};

    F.function = &Nk_intd;
    F.params = &params;

    gsl_integration_qags (&F, 0, vesc, 1e-6, 1e-6, 1000 ,w, &result, &error);
    //Fast integration --> see GSL documentation for more info
    //size_t neval; gsl_integration_qng (&F, 0, vesc, 1e-1, 1e-1, &result, &error, &neval);
    gsl_integration_workspace_free (w);

    return result;
}

struct lisanti_params { double vesc_lis; double v0_lis; double ve_lis; double k_lis; int i_lis;};


int fint(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){

    struct lisanti_params * params
    = (struct lisanti_params *)p;
    double vesc_lis = (params->vesc_lis);
    double v0_lis = (params->v0_lis);
    double ve_lis = (params->ve_lis);
    double k_lis = (params->k_lis);
    int i_lis = (params->i_lis);
    //printf("\n i %i \n", i_lis);
    double A = exp(-(k[1]*k[1]+ve_lis*ve_lis+2.*k[1]*ve_lis*k[0])/(k_lis*v0_lis*v0_lis));
    double B = exp(-(vesc_lis*vesc_lis)/(k_lis*v0_lis*v0_lis));
    if ( A-B > 0.0 ){
	fval[0] = 2*M_PI*k[1]*pow(k[1],(i_lis*1.0))*pow(A-B,k_lis);

        //printf("\n halo %.5e | velocity %.5e | Cos(angle) %.5e | k %.5e| v0 %.5e | 1st term %.5e| 2nd %.5e | both %.5e  \n", fval[0], k[1], k[0], k_lis, v0_lis, A, B, pow(A-B,k_lis) );
    }
    else {
	fval[0] = 0;
    }

    return 0;

}

//INTEGRATION OF THE UNCERTAINTIES DISTRIBUTION

double lisanti_halo (double vmin, double vesc, double v0, double ve, double k, int i)

{

    double  val, err;
    struct lisanti_params params = {vesc,v0,ve,k,i};
    //printf("ve + vesc %lf\n", ve+vesc);
    //limits of integration
    double xl[2]={-1,vmin};
    double xu[2]={1,vesc+ve};

    hcubature(1, &fint, &params,
              2, xl, xu,
              0, 0, 1e-6, ERROR_INDIVIDUAL, &val, &err);

    //printf("Computed integral = %0.10g +/- %g\n", val, err);


    return val;

}



/*############ NORMALIZATIONS ##########################################################*/

void normalization(char * profile, double vesc, double v0, double beta, double vt, double vc, double k){
if (strncmp (profile,"SHM",10) == 0){
Normalization_shm = norm_shm_num(vesc, v0, beta);
}
if (strncmp (profile,"Lisanti",10) == 0){
Normalization_lisanti = Nk_uncert(vesc, v0, k);
}
else {
Normalization_fornasa = 0.0;
printf("Need to have a table we have\n");
}
}

/*****************************************************************************************
Velocity bins
*****************************************************************************************/
int velocity(double vesc, double ve){
   int i;
   for (i=0 ; i<length ; i++){
      vel[i] = (vesc+ve)*i/length;

   }
   return 0;
}

/*****************************************************************************************
 int f(v)*v^(i+1) bins
*****************************************************************************************/
int fv_v(char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){
velocity(vesc,ve);
int m;
for(m=0 ; m<i+1 ; m++){
   if (strncmp (profile,"Lisanti",10) == 0){
   	int j;
    for (j=0 ; j<length ; j++){
      fv[m][j] = lisanti_halo (vel[j], vesc, v0, ve, k, m)/Normalization_lisanti;
   }}
	if (strncmp (profile,"SHM",10) == 0){
   	int j;
    for (j=0 ; j<length ; j++){
      fv[m][j] = shm_halo (vel[j], vesc, v0, ve, beta,m)/Normalization_shm;
   }}
}

   return 0;
}

void define_halo (char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){
if(i>power-1){ printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n",power-1);}
else{
if (strncmp (profile,"SHM",10) == 0 || strncmp (profile,"Lisanti",10) == 0){

printf("Calculating halo integrals for %s up to order %d in v...\n",profile,i);
normalization(profile,vesc,v0,beta,vt,vc,k);
fv_v(profile, vesc, v0, beta, vt, vc, ve, k, i);
printf("Done!\n");

}

else { printf("The profile does NOT exist! please use: Lisanti or SHM\n");}
}
}

/*#################################################
               Interpolation
##################################################*/

double halo (double vmin, int i){
   int m;
   double x[length], y[length];
   double result;


   for (m = 0; m < length; m++)
   {
      x[m] = vel[m];
      y[m] = fv[i][m];

   }


   if (vmin >= x[0] && vmin <= x[length-1]){
       /*printf("%f\n", x[length-1]);*/
       gsl_interp_accel *fofv = gsl_interp_accel_alloc ();
       gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, length);

       gsl_spline_init (spline, x, y, length);

       result = gsl_spline_eval (spline,vmin, fofv);

       gsl_spline_free (spline);
       gsl_interp_accel_free (fofv);

       return result;

   }
   else {  // f(v) data out of the range
      return 0.;
   }
}

double halo_w (double vmin, gsl_interp_accel *ga, gsl_spline * gs ){
   if (vmin >= vel[0] && vmin <= vel[length-1])
   {
       return gsl_spline_eval (gs, vmin, ga);
   }
   else {  // f(v) data out of the range
      return 0.;
   }
}

double halo_f (char * profile, double vmin, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){
if (strncmp (profile,"SHM",10) == 0){
double N = norm_shm_num(vesc, v0, beta);
return shm_halo (vmin, vesc, v0, ve, beta,i)/N;
}
if (strncmp (profile,"Lisanti",10) == 0){
double N = Nk_uncert(vesc, v0, k);
return lisanti_halo (vmin, vesc, v0, ve, k, i)/N;
}
else{
printf("We do not have this halow implemented\n");
return 0.0;
}
}


/*########################################################################

				Modified functions for writing and reading
									30-07-2015

###########################################################################*/

int write_fv_v(char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){
	fv_v(profile, vesc, v0, beta, vt, vc, ve, k, i);
	FILE * table;

	if (i == 0){

		table = fopen("halo_table.dat", "w+");
		fprintf(table, "%d %d \r\n", length, i);
		int j;
		for (j = 0; j < length; j++){
			printf("\n i %i \n", i);
			fprintf(table, "%.5E %.5E \r\n", vel[j], fv[0][j]);
		}
		fclose(table);

	}
	if (i == 1){

		table = fopen("halo_table.dat", "w+");
		fprintf(table, "%d %d \r\n", length, i);
		int j;
		for (j = 0; j < length; j++){
			printf("\n i %i \n", i);
			fprintf(table, "%.5E %.5E %.5E \r\n", vel[j], fv[0][j], fv[1][j]);
		}
		fclose(table);

	}
	if (i == 2){

		table = fopen("halo_table.dat", "w+");
		fprintf(table, "%d %d \r\n", length, i);
		int j;
		for (j = 0; j < length; j++){

			fprintf(table, "%.5E %.5E %.5E %.5E \r\n", vel[j], fv[0][j], fv[1][j], fv[2][j]);
		}
		fclose(table);

	}
	else printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1);

	return 0;
}

int write_fv_v_varpath(char * path, char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){
	fv_v(profile, vesc, v0, beta, vt, vc, ve, k, i);
	FILE * table;

	if (i == 0){

		table = fopen(path, "w+");
		fprintf(table, "%d %d \r\n", length, i);
		int j;
		for (j = 0; j < length; j++){
			fprintf(table, "%.5E %.5E \r\n", vel[j], fv[0][j]);
		}
		fclose(table);

	}
	if (i == 1){

		table = fopen(path, "w+");
		fprintf(table, "%d %d \r\n", length, i);
		int j;
		for (j = 0; j < length; j++){
			fprintf(table, "%.5E %.5E %.5E \r\n", vel[j], fv[0][j], fv[1][j]);
		}
		fclose(table);

	}
	if (i == 2){

		table = fopen(path, "w+");
		fprintf(table, "%d %d \r\n", length, i);
		int j;
		for (j = 0; j < length; j++){
			fprintf(table, "%.5E %.5E %.5E %.5E \r\n", vel[j], fv[0][j], fv[1][j], fv[2][j]);
		}
		fclose(table);

	}
	else printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1);

	return 0;
}

void define_and_write_halo(char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){


	if (i>power - 1){ printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1); }
	else{
		if (strncmp(profile, "NFW", 10) == 0 || strncmp(profile, "Einasto", 10) == 0 ||
			strncmp(profile, "Burkert", 10) == 0 || strncmp(profile, "SHM", 10) == 0
			|| strncmp(profile, "Lisanti", 10) == 0){

			printf("Calculating halo integrals for %s up to order %d in v...\n", profile, i);
			normalization(profile, vesc, v0, beta, vt, vc, k);
			write_fv_v(profile, vesc, v0, beta, vt, vc, ve, k, i);
			printf("Done!\n");

		}

		else { printf("The profile does NOT exist! please use: NFW, Einasto, Burkert, Lisanti or SHM\n"); }
	}
}

void define_and_write_halo_path(char* path, char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i){


	if (i>power - 1){ printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1); }
	else{
		if (strncmp(profile, "NFW", 10) == 0 || strncmp(profile, "Einasto", 10) == 0 ||
			strncmp(profile, "Burkert", 10) == 0 || strncmp(profile, "SHM", 10) == 0
			|| strncmp(profile, "Lisanti", 10) == 0){

			//printf("Calculating halo integrals for %s up to order %d in v...\n", profile, i);
			normalization(profile, vesc, v0, beta, vt, vc, k);
			write_fv_v_varpath(path, profile, vesc, v0, beta, vt, vc, ve, k, i);
			//printf("Done!\n");

		}

		else { printf("The profile does NOT exist! please use: NFW, Einasto, Burkert, Lisanti or SHM\n"); }
	}
}

void read_halo(char* path){

	FILE * table;
	int i;
	int j;
	int check_length;

	//printf("Reading pre-calculated halo table from file %s...\n", path);
	table = fopen(path, "r");

	//printf("Getting table format...\n");
	fscanf(table, "%d %d", &check_length, &i);
	//printf("Done!\n");

	//printf("Reading data table...\n");
		for (j = 0; j < length; j++){

			if (i == 0) fscanf(table, "%E %E", &vel[j], &fv[0][j]);
			if (i == 1) fscanf(table, "%E %E %E", &vel[j], &fv[0][j], &fv[1][j]);
			if (i == 2) fscanf(table, "%E %E %E %E", &vel[j], &fv[0][j], &fv[1][j], &fv[2][j]);

		}

	fclose(table);
	//printf("Done!\n");
}

int access_check(char* path){

	printf("Checking if halo table exists already...\n");

	if (access(path, F_OK) != -1){
		printf("Halo table exists already. If new calculation is desired, please delete or rename the file (%s)...\n", path);
		return 1;
	}
	else{
		printf("No file %s found. New table will be calculated...\n", path);
		return 0;
	}

}

int access_check_time(){

	printf("Checking if annual modulated halo table exists already...\n");
	if (access("halo_table/halo_table_0.dat", F_OK) != -1){
		printf("Halo table exists already. If new calculation is desired, please delete or rename the file (halo_table_0.dat)...\n");
		return 1;
	}
	else{
		printf("No file halo_table_0.dat found. New table will be calculated...\n");
		mkdir("halo_table", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		return 0;
	}

}

double time_ve(double ve, double ve0, double t0, double T, int t){

	return ve + ve0*cos(2*M_PI*(t-t0)/T);


}
void define_and_write_halo_time(char * profile, double vesc, double v0, double beta, double vt, double vc, double ve, double k, int i, double ve0, double t0, double T){
  double vsun = sqrt(vc*vc + vc*12.24 + 11.10*11.10 + 12.24*12.24 + 7.25*7.25);
	if (i>power - 1){ printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1); }
	else{
		if (strncmp(profile, "NFW", 10) == 0 || strncmp(profile, "Einasto", 10) == 0 ||
			strncmp(profile, "Burkert", 10) == 0 || strncmp(profile, "SHM", 10) == 0
			|| strncmp(profile, "Lisanti", 10) == 0){

			printf("Calculating halo integrals for %s up to order %d in v...\n", profile, i);
			normalization(profile, vesc, v0, beta, vt, vsun, k);

			clock_t start, end;
			double cpu_time_used;

			start = clock();

			int t;
			for (t = 0; t < T; t++){
				char path[32];
				snprintf(path, sizeof(char) * 32, "halo_table/halo_table_%i.dat", t);
				printf("Calculating for t=%i...\n", t);
				write_fv_v_varpath(path, profile, vesc, v0, beta, vt, vsun, time_ve(ve, ve0, t0, T, t), k, i);
			}

			printf("Done!\n");
			end = clock();
			cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
			printf("Elapsed time: %f seconds\n", cpu_time_used);

		}

		else { printf("The profile does NOT exist! please use: NFW, Einasto, Burkert, Lisanti or SHM\n"); }
	}
}
