/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
difCrosSec.c
difCrosSEffective operators for DM direct detectionec.c
early versions of this code were developed by M. Peiró
and E. Gerstmayr. 

Subroutine calculating equation (53) and (55)
          of arXiv: 1203.3542  

##########################################*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "effFormFact.h"
#include "phys_consts.h"
#include "coeffs_eft.h"
#include "MenenFF.h"

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define VERBOSE 0

/*########################################################################################
Since we take the nuclear responses from effFormFact.c, we have divided the differential
cross section into two pieces, a velocity independent and a v^2 pieces.  
########################################################################################*/

double difcros_isotope_v0_dEr(int A, int Z, double Er, double mchi, double jchi, double v, char * Nucleon, int F_i, int F_j){

	double mN = approx_mass_nucleus(A, Z);
	/*double E = 0.5*mchi*pow(v / c, 2.);*/
        double E = 0.5*mchi*v*v/c/c;
	/*double r = 4.*mchi*mN / pow(mchi + mN, 2.);*/
        double mchimN_2 =(mchi + mN)*(mchi + mN);
	double r = 4.*mchi*mN /mchimN_2;
	double Jacobian = 2. / (r*E);
	//printf("here is v %f \n", v);
	//double prefactor = pow(mN, 2.) / (32.*M_PI*pow(mchi + mN, 2.)*pow(mproton, 2.))*Jacobian;
        /*double prefactor = pow(mN, 1.) / (32.*M_PI*pow(mchi + mN, 2.))*Jacobian;*/
        double prefactor = mN / (32.*M_PI*mchimN_2)*Jacobian;


	if (strncmp(Nucleon, "pp", 10) == 0){
		return prefactor*Cp(F_i)*Cp(F_j)*FormFact_v0(A, Z, F_i, F_j, "p", "p", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "nn", 10) == 0){
		return prefactor*Cn(F_i)*Cn(F_j)*FormFact_v0(A, Z, F_i, F_j, "n", "n", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "pn", 10) == 0){
		return prefactor*Cp(F_i)*Cn(F_j)*FormFact_v0(A,Z, F_i, F_j, "p", "n", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "np", 10) == 0){
		return prefactor*Cp(F_i)*Cn(F_j)*FormFact_v0(A,Z, F_i, F_j, "n", "p", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "All", 10) == 0){
		return prefactor*(Cp(F_i)*Cp(F_j)*FormFact_v0(A, Z, F_i,F_j, "p", "p", Er, mchi, jchi) + Cn(F_i)*Cn(F_j)*FormFact_v0(A, Z, F_i,F_j, "n", "n", Er, mchi, jchi) +
			Cp(F_i)*Cn(F_j)*FormFact_v0(A, Z, F_i,F_j, "p", "n", Er, mchi, jchi) + Cn(F_i)*Cp(F_j)*FormFact_v0(A, Z, F_i,F_j, "n", "p", Er, mchi, jchi));
	}
	if (strncmp(Nucleon, "ISO", 10) == 0){
		return prefactor*(Cp(F_i)*Cp(F_j)*FormFactMen_v0(A, Z, F_i,F_j, "+", "+", Er, mchi, jchi) + Cn(F_i)*Cn(F_j)*FormFactMen_v0(A, Z, F_i,F_j, "-", "-", Er, mchi, jchi) +
			Cp(F_i)*Cn(F_j)*FormFactMen_v0(A, Z, F_i,F_j, "+", "-", Er, mchi, jchi) + Cn(F_i)*Cp(F_j)*FormFactMen_v0(A, Z, F_i,F_j, "-", "+", Er, mchi, jchi));
	}
	else return 0.;

}

double difcros_isotope_v2_dEr(int A, int Z, double Er, double mchi, double jchi, double v, char * Nucleon, int F_i, int F_j){

	double mN = approx_mass_nucleus(A, Z);
	/*double E = 0.5*mchi*pow(v / c, 2.);*/
        double E = 0.5*mchi*v*v/c/c;
	/*double r = 4.*mchi*mN / pow(mchi + mN, 2.);*/
        double mchimN_2 =(mchi + mN)*(mchi + mN);
	double r = 4.*mchi*mN /mchimN_2;
	double Jacobian = 2. / (r*E);
	//double prefactor = pow(mN, 2.) / (32.*M_PI*pow(mchi + mN, 2.)*pow(mproton, 2.))*Jacobian;
        /*double prefactor = pow(mN, 1.) / (32.*M_PI*pow(mchi + mN, 2.))*Jacobian;*/
        double prefactor = mN / (32.*M_PI*mchimN_2)*Jacobian;


	if (strncmp(Nucleon, "pp", 10) == 0){
		return prefactor*Cp(F_i) * Cp(F_j) * FormFact_v2(A, Z, F_i,F_j, "p", "p", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "nn", 10) == 0){
		return prefactor*Cn(F_i) * Cn(F_j) * FormFact_v2(A, Z, F_i,F_j, "n", "n", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "pn", 10) == 0){
		return prefactor*Cp(F_i) * Cn(F_j) * FormFact_v2(A, Z, F_i,F_j, "p", "n", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "np", 10) == 0){
		return prefactor*Cp(F_i) * Cn(F_j) * FormFact_v2(A, Z, F_i,F_j, "n", "p", Er, mchi, jchi);
	}
	if (strncmp(Nucleon, "All", 10) == 0){
		return prefactor*(Cp(F_i) * Cp(F_j) * FormFact_v2(A, Z, F_i,F_j, "p", "p", Er, mchi, jchi) + Cn(F_i) * Cn(F_j) * FormFact_v2(A, Z, F_i,F_j, "n", "n", Er, mchi, jchi) +
			Cp(F_i) * Cn(F_j) * FormFact_v2(A, Z, F_i,F_j, "p", "n", Er, mchi, jchi) + Cn(F_i) * Cp(F_j) * FormFact_v2(A, Z, F_i,F_j, "n", "p", Er, mchi, jchi));
	}

	if (strncmp(Nucleon, "ISO", 10) == 0){
		return prefactor*(Cp(F_i)*Cp(F_j)*FormFactMen_v2(A, Z, F_i,F_j, "+", "+", Er, mchi, jchi) + Cn(F_i)*Cn(F_j)*FormFactMen_v2(A, Z, F_i,F_j, "-", "-", Er, mchi, jchi) +
			Cp(F_i)*Cn(F_j)*FormFactMen_v2(A, Z, F_i,F_j, "+", "-", Er, mchi, jchi) + Cn(F_i)*Cp(F_j)*FormFactMen_v2(A, Z, F_i,F_j, "-", "+", Er, mchi, jchi));
	}
	else return 0.;

}


