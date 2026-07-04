/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
##########################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>


#include "cubature.h"
#include "source/phys_consts.h"


/*From the source folder we include the following files*/
#include "source/halo.h"

#include "source/difRateGen.h"
#include "source/coeffs_eft.h"
#include "source/WriteData.h"

#define VERBOSE 0

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


//=====================================================================================
// MAIN CODE
//=====================================================================================

int main() {
  /*########################################################################


                                  Astrophysical and Particle Physics Input


  ###########################################################################*/

  



  //read_coeffs("coeffs_table.dat");
  //double mass;
  //mass = give_mass();
  double mass = 10.0;
  printf("%.5e", mass);
  
  // struct difcros_params struct_difcros_Xe = {"Xe", 1., mass, 0.5, 220, "All", 0.0};
  // struct difcros_params struct_difcros_Ge = {"Ge", 1., mass,0.5,220.,"All", 0.0};
  // struct difcros_params struct_difcros_Ar = {"Ar", 1., mass, 0.5, 220, "All", 0.0};
  // struct difcros_params struct_difcros_F = {"F", 1., mass, 0.5, 220, "All", 0.0};

  // //{char * target; double Er; double mchi; double jchi; double v; char * Nucleon; double exposure;};

  // /*########################################################################


  //             Calculate halo files or retrieve data from table (halo_table.dat)


  //  ###########################################################################*/


  // double rhochi = 0.4;
  // struct halo_params struct_halo = { "SHM", 544, 220., 0., 90., 245., 220., 0., 2, 15., 151., 365., halo_file_path};
  // printf("ve %lf \n", struct_halo.ve);
  // define_and_write_halo(struct_halo.profile, struct_halo.vesc, struct_halo.v0, struct_halo.beta, struct_halo.vt, struct_halo.vc, struct_halo.ve, struct_halo.k, struct_halo.i);


  // if (access_check(struct_halo.halo_path) != 1){
  //   define_and_write_halo(struct_halo.profile, struct_halo.vesc, struct_halo.v0, struct_halo.beta, struct_halo.vt, struct_halo.vc, struct_halo.ve, struct_halo.k, struct_halo.i);
  // }
  // if (access_check_time() != 1) {
  //   define_and_write_halo_time(struct_halo.profile, struct_halo.vesc, struct_halo.v0, struct_halo.beta, struct_halo.vt, struct_halo.vc, struct_halo.ve, struct_halo.k, struct_halo.i, struct_halo.ve0, struct_halo.t0, struct_halo.T);
  // }


  //write_binned_data(coeffs_table, experiment_tables , results_path, struct_halo.halo_path, rhochi, &struct_difcros_Xe, "None");
  //write_binned_data_and_each_operator(argv[1], argv[2], argv[3], struct_halo.halo_path, rhochi, &struct_difcros_Ge, "None");
  //write_binned_data_and_each_operator(argv[1], argv[2], argv[3], struct_halo.halo_path, rhochi, &struct_difcros_Ar, "None");
  //write_binned_data_and_each_operator(argv[1], argv[2], argv[3], struct_halo.halo_path, rhochi, &struct_difcros_F, "None");
  //write_difrate_total(argv[1], argv[2], struct_halo.halo_path, rhochi, &struct_difcros_Xe, "None");
  //printf("here is the 1st isotope %i, with has abundance %lf\n", atomic_numbers_Ge[0], isotopes_Ge[0]);




  return 0;


}
