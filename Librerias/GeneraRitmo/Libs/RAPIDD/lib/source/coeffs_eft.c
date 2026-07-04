/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
			coeffs_eft.c
Set EFT coefficients for calculations here
##########################################*/

#include "coeffs_eft.h"
#include "phys_consts.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double C[OPERATOR_NUMBER];
double CP[OPERATOR_NUMBER];
double CN[OPERATOR_NUMBER];
double CP_trace[OPERATOR_NUMBER];
double CN_trace[OPERATOR_NUMBER];
double C0[OPERATOR_NUMBER];
double C1[OPERATOR_NUMBER];
double C0_gen_SI[5];
double C1_gen_SI[5];
double MASS;
double MED_MASS;
double rho;
double vesc;
double v0;
double k;
double vc;



void set_coeffs(){
	//Indices in papers run from 1 to 15, so define first element of array to be zero, out of convenience//
	CP[0] = 0;
	C[0] = 0;
	CN[0] = 0;
	//Operator 1
	C[1] = 0.0;
    CP[1] = 0.0;
    CN[1] = 0.0;
	//Operator 2
	C[2] = 0.0;
	CP[2] = 0.0;
	CN[2] = 0.0;
	//Operator 3
	C[3] = 0.0;
	CP[3] = 0.0;
	CN[3] = 0.0;
	//Operator 4
	C[4] = 0.0;
	CP[4] = 0.0;
	CN[4] = 0.0;
	//Operator 5
	C[5] = 0.0;
	CP[5] = 0.0;
	CN[5] = 0.0;
	//Operator 6
	C[6] = 0.0;
	CP[6] = 0.0;
	CN[6] = 0.0;
	//Operator 7
	C[7] = 0.0;
	CP[7] = 0.0;
	CN[7] = 0.0;
	//Operator 8
	C[8] = 0.0;
	CP[8] = 0.0;
	CN[8] = 0.0;
	//Operator 9
	C[9] = 0.0;
	CP[9] = 0.0;
	CN[9] = 0.0;
	//Operator 10
    C[10] = 0.0;
    CP[10] = 0.0;
	CN[10] = 0.0;
	//Operator 11
	C[11] = 0.0;
	CP[11] = 0.0;
	CN[11] = 0.0;
	//Operator 12
	C[12] = 0.0;
	CP[12] = 0.0;
	CN[12] = 0.0;
	//Operator 13
	C[13] = 0.0;
	CP[13] = 0.0;
	CN[13] = 0.0;
	//Operator 14
	C[14] = 0.0;
	CP[14] = 0.0;
	CN[14] = 0.0;
	//Operator 15
	C[15] = 0.0;
	CP[15] = 0.0;
	CN[15] = 0.0;
	//Mass
	MASS = 0.0;
  //MED_MASS = 0.0;
}


double Cp(int i){
	return CP[i];
}

double Cn(int i){
	return CN[i];
}

double Cp_trace(int i){
	return CP_trace[i];
}

double Cn_trace(int i){
	return CN_trace[i];
}

double give_mass(){
  return MASS;
}

double give_med_mass(){
  return MED_MASS;
}

double give_rho(){
  return rho;
}

double give_vesc(){
  return vesc;
}

double give_v0(){
  return v0;
}

double give_k(){
  return k;
}

double give_vc(){
  return vc;
}

void print_coeffs() {
  for (int i = 0; i < OPERATOR_NUMBER; i++) {
    if (fabs(CP[i])>0 || fabs(CN[i]) >0 )
      printf("Set parameters for operator %i: %e (proton) %e (neutron)\n",i, CP[i], CN[i]);
  }
}

void print_SI_coeffs() {
  for (int i = 0; i < 5; i++) {
    if (fabs(C0_gen_SI[i])>0 || fabs(C1_gen_SI[i]) >0 )
      printf("Set parameters for operator %i: %e (proton) %e (neutron)\n",i, C0_gen_SI[i], C1_gen_SI[i]);
  }
}

/// Read in Parameters from file fname
void read_coeffs(const char* fname){
  FILE * coeffs_table;

  CP[0] = 0;
  CN[0] = 0;

  printf("Reading input parameters from %s\n", fname);
  coeffs_table = fopen(fname, "r");
  for (int i = 0; i < OPERATOR_NUMBER; i++) {
    if (i==0) {
      fscanf(coeffs_table, "%lf", &MASS);
    }
    else {
      fscanf(coeffs_table, "%lf %lf", &CP[i], &CN[i]);
    }
  }
  fclose(coeffs_table);
}

void read_coeffs_low_MED(const char* fname){
  FILE * coeffs_table;

  CP[0] = 0;
  CN[0] = 0;

  printf("Reading input parameters from %s\n", fname);
  coeffs_table = fopen(fname, "r");
  for (int i = 0; i < OPERATOR_NUMBER; i++) {
    if (i==0) {
      fscanf(coeffs_table, "%lf %lf", &MASS, &MED_MASS);
    }
    else {
      fscanf(coeffs_table, "%lf %lf", &CP[i], &CN[i]);
    }
  }
  fclose(coeffs_table);
}

/// Read in Parameters in isoscalar and isovector form
void read_coeffs_iso(const char* fname){
  FILE * coeffs_table;

  CP[0] = 0;
  CN[0] = 0;



  printf("Reading input parameters from %s\n", fname);
  coeffs_table = fopen(fname, "r");
  for (int i = 0; i < OPERATOR_NUMBER; i++) {
    if (i==0) {
      fscanf(coeffs_table, "%lf", &MASS);
    }
    else {
      fscanf(coeffs_table, "%lf %lf", &C0[i], &C1[i]);
      CP[i] = C0[i] + C1[i];
      CN[i] = C0[i] - C1[i];
    }
  }
  fclose(coeffs_table);
}

void read_halo_params(const char* fname){
  FILE * coeffs_table;

  CP[0] = 0;
  CN[0] = 0;

  printf("Reading halo parameters from %s\n", fname);
  coeffs_table = fopen(fname, "r");
  for (int i = 0; i < 5; i++) {
    if ( i == 0){
      fscanf(coeffs_table, "%lf", &rho);
      //printf("\n rho %f\n",rho);
    }
    if ( i == 1){
      fscanf(coeffs_table, "%lf", &vesc);
    }
    if ( i == 2){
      fscanf(coeffs_table, "%lf", &v0);
    }
    if ( i == 3){
      fscanf(coeffs_table, "%lf", &k);
    }
    if (i == 4){
      //printf("no way\n" );
      fscanf(coeffs_table, "%lf", &vc);
    }
  }
  fclose(coeffs_table);
}

void write_coeffs(char * output_path){

	FILE * coeffs_table;

	int i;

	coeffs_table = fopen(output_path, "w+");

	for (i = 0; i < OPERATOR_NUMBER; i++){
		if (i==0){
			fprintf(coeffs_table, "%.5E \r\n", give_mass());
		}
		else {
			fprintf(coeffs_table, "%.5E %.5E \r\n", Cp(i), Cn(i));
		}
	}
	fclose(coeffs_table);

}

void set_coeffs_matching_Standard( double mchi, double sigma, double sigmaSD){

    int i;
    double pbGeVfactor = 2.67e-9;
    for (i=0; i < OPERATOR_NUMBER; i++){
        if(i != 1){
            CP[i] = 0;
            CN[i] = 0;
        }
        if ( i == 1){
            double fn = sigma_to_f_11(sigma, mchi);
            double fp = sigma_to_f_11(sigma, mchi);
            //CP[i] = 8* mproton*mchi*fp*sqrt(pbGeVfactor);
            //CN[i] = 8* mneutron*mchi*fn*sqrt(pbGeVfactor);
            CP[i] = 2*higgs_vev*higgs_vev*fp*sqrt(pbGeVfactor);
            CN[i] = 2*higgs_vev*higgs_vev*fn*sqrt(pbGeVfactor);
            //printf("\n coeffs %.5E \n", CP[i]);

        }
	else if ( i == 4){
	    double a_44 = sigma_to_a_44(sigmaSD, mchi);
	    CP[i]=8*higgs_vev*higgs_vev*sqrt(2)*G_F*a_44*sqrt(pbGeVfactor);
   	    CN[i]=8*higgs_vev*higgs_vev*sqrt(2)*G_F*a_44*sqrt(pbGeVfactor);
         }
    }
}

void set_coeffs_matching_Standard_all( double mchi, double sigma_p,double sigma_n, double sigmaSD_p, double sigmaSD_n){

    int i;
    double pbGeVfactor = 2.67e-9;
    for (i=0; i < OPERATOR_NUMBER; i++){
        if (i == 0){
          MASS = mchi;
        }
        else if ( i == 1){
          double fn = sigma_to_f_11(sigma_n, mchi);
          double fp = sigma_to_f_11(sigma_p, mchi);
          //CP[i] = 8* mproton*mchi*fp*sqrt(pbGeVfactor);
          //CN[i] = 8* mneutron*mchi*fn*sqrt(pbGeVfactor);
          CP[i] = 2*higgs_vev*higgs_vev*fp*sqrt(pbGeVfactor);
          CN[i] = 2*higgs_vev*higgs_vev*fn*sqrt(pbGeVfactor);
          //printf("\n coeffs %.5E \n", CP[i]);

        }
      	else if ( i == 4){
	        double ap_44 = sigma_to_a_44(sigmaSD_p, mchi);
			    double an_44 = sigma_to_a_44(sigmaSD_n, mchi);
	        CP[i]=8*higgs_vev*higgs_vev*sqrt(2)*G_F*ap_44*sqrt(pbGeVfactor);
   	      CN[i]=8*higgs_vev*higgs_vev*sqrt(2)*G_F*an_44*sqrt(pbGeVfactor);
        }
        else{
          CP[i] = 0.0;
          CN[i] = 0.0;
        }   
    }
}

void set_any_coeffs(double C, int i){

    if ( i == 0){

    MASS = C;
    }

    if( i > 0 && i < OPERATOR_NUMBER){

    CP[i] = C;
    CN[i] = C;
    }

}

void set_any_Ncoeff(double input, int i, char* nucleon){

	if ( strncmp (nucleon,"n",10) == 0){
		CN[i] = input;
		//printf("%s %i %lf", nucleon, i, CN[i]);
	}
	if ( strncmp (nucleon,"p",10) == 0){
		CP[i] = input;
		//printf("%s %i %lf", nucleon, i, CP[i]);
	}
}

void set_med_mass(double med_mass){
  MED_MASS = med_mass;
}
void set_any_coeffs_trace(double input, int i, char*nucleon){


	if ( strncmp (nucleon,"n",10) == 0){
		CN_trace[i] = input;
		//printf("%s %i %lf", nucleon, i, CN[i]);
	}
	if ( strncmp (nucleon,"p",10) == 0){
		CP_trace[i] = input;
		//printf("%s %i %lf", nucleon, i, CP[i]);
	}

}

void set_coeffs_anapole(double mchi, double ANAPOLE){
    int i;
		MASS = mchi;
    double pbGeVfactor = 2.67e-9;
    double mup = reduced_mass(mchi, mproton);
    for (i=0; i < OPERATOR_NUMBER; i++){
        if(i != 8 && i != 9){
            CP[i] = 0;
            CN[i] = 0;
        }
        if (i == 8){
            CN[i] = 0;
            CP[i] = ANAPOLE*higgs_vev*higgs_vev*2*e;
        }
        if (i == 9){
            CP[i] = -ANAPOLE*higgs_vev*higgs_vev*e*g_p;
            CN[i] = -ANAPOLE*higgs_vev*higgs_vev*e*g_n;
        }
    }
}

void set_coeffs_magmom(double mchi, double mu_x){
	int i;
	MASS = mchi;
	double mup = reduced_mass(mchi, mproton);
	for (i=0; i < OPERATOR_NUMBER; i++){
			if(i != 1 && i != 4 && i !=5 && i != 6){
				CP[i] = 0;
				CN[i] = 0;
			}
			else if(i==1){
				CN[i] = 0.0;
				CP[i] = mu_x*higgs_vev*higgs_vev*e/(2*mchi);
			}
			else if(i==4){
				CN[i] = (g_n *mu_x*higgs_vev*higgs_vev *e)/(mneutron);
				CP[i] = (g_p *mu_x*higgs_vev*higgs_vev*e)/(mproton);
			}
			else if(i==5){
				CN[i] = 0.0;
				CP[i] = 2*mu_x*higgs_vev*higgs_vev*e*mproton;
			}
			else if(i==6){
				CN[i] = -(g_n*mu_x*higgs_vev*higgs_vev*e*mneutron);
				CP[i] = -(g_p*mu_x*higgs_vev*higgs_vev*e*mproton);
			}
		}
}

void set_coeffs_elecmom(double mchi, double d_x){
	int i;
	MASS = mchi;
  MED_MASS = 0.0;
	double mup = reduced_mass(mchi, mproton);
	for (i=0; i < OPERATOR_NUMBER; i++){
			if(i != 11){
				CP[i] = 0;
				CN[i] = 0;
			}
			else if(i==11){
				CN[i] = 0.0;
				CP[i] = 2*d_x*e*higgs_vev*higgs_vev*mproton;
			}
		}
}

void set_coeffs_charrad(double mchi, double b_x){
	int i;
	MASS = mchi;
	double mup = reduced_mass(mchi, mproton);
	for (i=0; i < OPERATOR_NUMBER; i++){
			if(i != 1 ){
				CP[i] = 0;
				CN[i] = 0;
			}
			else if(i==1){
				CN[i] = 0.0;
				CP[i] = b_x*higgs_vev*higgs_vev*e;
			}
		}
}

void coeffs_set_master(char* coeffs, char* model){
  if ( strncmp(model, "None", 10) == 0){
    	set_coeffs();
    	read_coeffs(coeffs);
  }
  else if( strncmp(model, "Light_Med", 10) == 0){
    	set_coeffs();
      printf("light mediator\n");
    	read_coeffs_low_MED(coeffs);
  }
  else {
  	printf("Implement this model %s\n", model);
  }
  printf("leaving\n");
}

void copy_trace(){
  for (int i = 0; i < OPERATOR_NUMBER; i++){
    CP[i] = Cp_trace(i);
    CN[i] = Cn_trace(i);
  }
}