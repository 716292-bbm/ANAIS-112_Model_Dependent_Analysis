#ifndef BLABLA
#define BLABLA
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

#pragma once

const double static mproton = 0.938272046; // proton mass
const double static mneutron = 0.939565560; // neutron mass
const double static c = 2.99792458e+5; //Speed of light in km/s
const double static G_F = 1.1664e-5; //Fermi coupling constant in GeV^-2
const double static higgs_vev = 246.2 ; //GeV 
const double static e = 0.3; // e charge
const double static s_w = 0.4895916665957459;
const double static g_n = -3.83;
const double static g_p = 5.59;

const double static isotopes_Ar[1] = {1.0};
const int static atomic_numbers_Ar[1] = {40};

const double static isotopes_F[1] = {1.0};
const int static atomic_numbers_F[1] = {19};

const double static isotopes_Ge[5] = {0.2123, 0.2766, 0.0773, 0.3594, 0.0774};
const int static atomic_numbers_Ge[5] = {70,72,73,74,76};

const double static isotopes_Xe[7] = {0.0192, 0.2644, 0.0408, 0.2118, 0.2689,0.1044,0.0887};
const int static atomic_numbers_Xe[7] = {128,129,130,131,132,134,136};

const double static isotopes_CaWO4[2] = {0.616, 0.0};
//const double static isotopes_CaWO4[2] = {0.616, 0.384};
const int static atomic_numbers_CaWO4[2] = {16,40};
const int static Z_numbers_CaWO4[2] = {8,20};
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double reduced_mass(double mass_1, double mass_2);
double approx_mass_nucleus(double A, double Z);

double sigma_to_f_11(double sigma, double mchi);
double sigma_to_a_44(double sigma, double mchi);


//double isotopes_Ar[1];
//int atomic_numbers_Ar[1];

//double isotopes_F[1];
//int atomic_numbers_F[1];

//double isotopes_Ge[5];
//int atomic_numbers_Ge[5];

//double isotopes_Xe[7];
//int atomic_numbers_Xe[7];


//double isotopes_CaWO4[2];
//int atomic_numbers_CaWO4[2];
//int Z_numbers_CaWO4[2];

int countlines(char *filename);

#endif
