/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
Physical constants
##########################################*/
#include <math.h>
#include<stdio.h>
#include<stdlib.h>

#include "phys_consts.h"


/////Frequently used equations///////////

double reduced_mass(double mass_1, double mass_2){

	return  mass_1*mass_2 / (mass_1 + mass_2);

}

////////Function to calculate the approximate mass of a nucleus by adding the masses of neutrons and protons together////////////
double approx_mass_nucleus(double A, double Z){

	return Z*mproton + (A - Z)*mneutron;

}

////////Function to calculate couplings f (f_n, f_p) for SI interaction by using the cross section sigma//////////////////////////
double sigma_to_f_11(double sigma, double mchi){

	double mun = reduced_mass(mchi, mproton); //reduced mass DM particle and nucleon
	return sqrt(sigma*M_PI) / (2.*mun);

}

////////Function to calculate couplings a (a_n, a_p) for SD interaction by using the cross section sigma//////////////////////////
double sigma_to_a_44(double sigma, double mchi){

	double mun = reduced_mass(mchi, mproton); //reduced mass DM particle and nucleon
	return sqrt(sigma*M_PI) / (sqrt(24.)*mun*G_F);

}


int countlines(char *filename){
  // count the number of lines in the file called filename
  FILE *fp = fopen(filename,"r");
  int ch=0;
  int lines=0;

  if (fp == NULL){
  return 0;
  }


  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }

  fclose(fp);
	//printf("attempt %i\n", lines );
  return lines;
}
