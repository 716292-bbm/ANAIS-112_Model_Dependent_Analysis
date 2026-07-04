//========================================================
// EFFICIENCY table reader
//========================================================
// A.Cheek 11/10/2016

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 
#include<math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_spline.h>

#include "phys_consts.h"

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

int length;
float data[1000];
float er[1000];
char eff_path_true[26];
/*
const int datacount = 100;
data = malloc(sizeof(int) * datacount);
er = malloc(sizeof(int) * datacount);
if (!data) {
  perror("Error allocating memory");
  abort();
}
if (!er) {
  perror("Error allocating memory");
  abort();
}*/

void read_efficiency(char* path){
  FILE * table;
  int i;
  int j;
  int check_length;
  length = countlines(path);
  //printf("%s\n", path );
  if (access(path, F_OK) == -1){
    //printf("Incorrect path to efficiency, taking as 1\n");
    snprintf(eff_path_true, sizeof(eff_path_true), "False");
    //printf("%s\n", eff_path_true );
  }
  else{
    snprintf(eff_path_true, sizeof(eff_path_true), "True");
    //printf("Reading pre-calculated table from file %s...\n", path);
    table = fopen(path, "r");

    //printf("Getting table format...\n");
    //printf("Done! %i\n", length);

    //printf("Reading data table...\n");
      for (j = 0; j < length; j++){
        fscanf(table, "%E %E", &er[j], &data[j]);

      }

      fclose(table);
      //printf("printing here %i\n",length );

  //printf("Done!\n");
}
}

double efficiency ( double E_r){
  if ( strncmp(eff_path_true, "True", 10) == 0){
    int m;
    double x[length], y[length];
    double result;

    for (m = 0; m < length; m++){
      x[m] = er[m];
      y[m] = data[m];
    }
    if (E_r >= x[0] && E_r <= x[length-1]){
        /*printf("%f\n", x[length-1]);*/
        gsl_interp_accel *eff = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, length);

        gsl_spline_init (spline, x, y, length);

        result = gsl_spline_eval (spline,E_r, eff);

        gsl_spline_free (spline);
        gsl_interp_accel_free (eff);
        if (result < 0.0) {
          return 0.0;
        }
        else{
          return result;
        }
      }
      else {  // eff(E_r) data out of the range
        return 0.;
      }
    }
  else{
    return 1.0;
  }
}

// double efficiency_w( double E_r, gsl_interp_accel *ga, gsl_spline * gs){
//     if (E_r >= x[0] && E_r <= x[length-1]){

//         result = gsl_spline_eval (ga ,E_r, ga);

//         if (result < 0.0) {
//           return 0.0;
//         }
//         else{
//           return result;
//         }
//       }
//       else {  // eff(E_r) data out of the range
//         return 0.;
//       }
//   }


void write_efficiency ( char* path, char* eff_path){
  read_efficiency(eff_path);
  int j;
  FILE* resultsfile;

  if (!(resultsfile = fopen(path, "w+"))){
    //printf("opening file failed...\n");
  }

  for (j = 0; j<100; j++){

    fprintf(resultsfile, "%.5E %.5E \r\n", j*1.0, efficiency(j*1.0) );

  }
  fclose(resultsfile);

}
