/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
WriteData.c 
##########################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "effFormFact.h"
#include "halo.h"
#include "difRateGen.h"
#include "phys_consts.h"
#include "binning_general.h"
#include "coeffs_eft.h"

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define VERBOSE 0

/*########################################################################################
Functions to write (differential) rates and form factors into .dat files
########################################################################################*/

/*########################################################################


Write Results into .dat files


###########################################################################*/

const char* write_header(char* header){
	for (int j=1; j<16; j++){
	//printf("\nhere\n");
		if( Cp(j) != 0.0 || Cn(j) != 0.0) {
			printf("\n in if\n");
			char element[10];
			set_any_coeffs_trace(Cn(j),j, "n");
			set_any_coeffs_trace(Cp(j),j, "p");
			printf("\n %i Couplings %.5E %.5E\n", j, Cn_trace(j), Cp_trace(j));
			snprintf(element, sizeof element, " c%i_c%i", j,j);
			strcat( header, element);
		}
		else{
			set_any_coeffs_trace(0.0,j, "n");
			set_any_coeffs_trace(0.0,j, "p");
		}
	}
	if ((Cp(1) != 0. || Cn(1) != 0.) & (Cp(3) != 0. || Cn(3) != 0.)){
		char element[10];
		snprintf(element, sizeof element, " c1_c3");
		strcat( header, element);

	}
	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(5) != 0. || Cn(5) != 0.)){
		char element[10];
		snprintf(element, sizeof element, " c4_c5");
		strcat( header, element);
	}
	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(6) != 0. || Cn(6) != 0.)){
		char element[10];
		snprintf(element, sizeof element, " c4_c6");
		strcat( header, element);
	}
	if ((Cp(8) != 0. || Cn(8) != 0.) & (Cp(9) != 0. || Cn(9) != 0.)){
		char element[10];
		snprintf(element, sizeof element, " c8_c9");
		strcat( header, element);
	}
	printf("header %s \n", header);
	return header;
}

void write_difrate_line(char* result_file, double rhochi, void * input_difcros, double line, char* model){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double prefact[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	char result[256];
	double energy = pow(10., line);
	val_difcros->Er = energy;
	copy_trace();
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 74;
		for (int l = 0; l < num_isos; l++){
			isotopes[l]=isotopes_Xe[l];	
			atomic_numbers[l]=atomic_numbers_Xe[l];	
			if (strncmp(model, "Light_Med", 10)==0){
				double med_mass = give_med_mass();
				printf("Here we are, mass =%lf\n", med_mass );
				prefact[l] = 1./(2*approx_mass_nucleus(atomic_numbers[l],znum)*energy*1.e-6 + med_mass*med_mass);
			}
			else{
				prefact[l] = 1.0;
			}
		}
	}
	if (strncmp(target, "Ge", 10) == 0){
		num_isos = 5;
		znum = 32; 
		for (int l = 0; l < num_isos; l++){
			isotopes[l]=isotopes_Ge[l];	
			atomic_numbers[l]=atomic_numbers_Ge[l];	
			if (strncmp(model, "Light_Med", 10)==0){
				double med_mass = give_med_mass();
				printf("Here we are, mass =%lf\n", med_mass );
				prefact[l] = 1./(2*approx_mass_nucleus(atomic_numbers[l],znum)*energy*1.e-6 + med_mass*med_mass);
			}
			else{
				prefact[l] = 1.0;
			}
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
		num_isos = 1;
		znum = 18; 
		for (int l = 0; l < num_isos; l++){
			isotopes[l]=isotopes_Ar[l];	
			atomic_numbers[l]=atomic_numbers_Ar[l];	
			if (strncmp(model, "Light_Med", 10)==0){
				double med_mass = give_med_mass();
				printf("Here we are, mass =%lf\n", med_mass );
				prefact[l] = 1./(2*approx_mass_nucleus(atomic_numbers[l],znum)*energy*1.e-6 + med_mass*med_mass);
			}
			else{
				prefact[l] = 1.0;
			}
		}
	}
	if (strncmp(target, "F", 10) == 0){
		num_isos = 1;
		znum = 9; 
		for (int l = 0; l < num_isos; l++){
			isotopes[l]=isotopes_F[l];	
			atomic_numbers[l]=atomic_numbers_F[l];	
			if (strncmp(model, "Light_Med", 10)==0){
				double med_mass = give_med_mass();
				printf("Here we are, mass =%lf\n", med_mass );
				prefact[l] = 1./(2*approx_mass_nucleus(atomic_numbers[l],znum)*energy*1.e-6 + med_mass*med_mass);
			}
			else{
				prefact[l] = 1.0;
			}
		}
	}

	// Output file name
	double counts[16];

	counts[0]=0.0;
	char result_tot[256];
	for (int l = 0; l< num_isos; l++){
		printf("check %i, %lf , %i, %lf\n", l, isotopes[l], atomic_numbers[l],total_difrate_isotope_dEr(atomic_numbers[l],znum, rhochi, val_difcros));
		counts[0] +=  prefact[l]*prefact[l]*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znum, rhochi, val_difcros);
	}
	//printf("the result %lf %lf\n", energy, counts[0]);
	snprintf(result, sizeof result, "%.5E %.5E", energy, counts[0]);
	for (int k = 1; k<16; k++){
		if(Cp_trace(k) == 0.0 && Cn_trace(k)==0.0){
			//printf("\n this should happen\n");
			continue;
		}
		if(Cp_trace(k) != 0. || Cn_trace(k) !=0.){
			//printf("\n why here %i %f %f", k, CN[k], CP[k] );
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(k), k, "p");
			set_any_Ncoeff(Cn_trace(k), k, "n");
			counts[k]=0.0;
			for (int l = 0; l< num_isos; l++){
				counts[k] +=  prefact[l]*prefact[l]*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],32, rhochi,val_difcros);
			}
			snprintf(reselement, sizeof reselement, " %.5E", counts[k]);
			printf("\n test %s %s \n", result, reselement);
			strcat(result, reselement);
			printf("result %s\n", result);
		}
	}
	if ((Cp_trace(1) != 0. || Cn_trace(1) != 0.) & (Cp_trace(3) != 0. || Cn_trace(3) != 0.)){
		char reselement[256];
		printf("\n I should not appear\n");
		set_coeffs();
		set_any_Ncoeff(Cp_trace(1), 1, "p");
		set_any_Ncoeff(Cn_trace(1), 1, "n");
		set_any_Ncoeff(Cp_trace(3), 3, "p");
		set_any_Ncoeff(Cn_trace(3), 3, "n");
		double counts13=0.0;
		for (int l = 0; l< num_isos; l++){
			counts13 +=  prefact[l]*prefact[l]*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],32, rhochi, val_difcros);
		}
		double int13 = counts13 - counts[1] - counts[3];
		snprintf(reselement, sizeof reselement, " %.5E", int13);
		strcat(result, reselement);
	}
	if ((Cp_trace(4) != 0. || Cn_trace(4) != 0.) & (Cp_trace(5) != 0. || Cn_trace(5) != 0.)){
		char reselement[256];
		//printf("\n I should not appear\n");
		set_coeffs();
		set_any_Ncoeff(Cp_trace(4), 4, "p");
		set_any_Ncoeff(Cn_trace(4), 4, "n");
		set_any_Ncoeff(Cp_trace(5), 5, "p");
		set_any_Ncoeff(Cn_trace(5), 5, "n");
		double counts45=0.0;
		for (int l = 0; l< num_isos; l++){
			counts45 +=  prefact[l]*prefact[l]*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],32, rhochi, val_difcros);
		}
		double int45 = counts45 - counts[4] - counts[5];
		snprintf(reselement, sizeof reselement, " %.5E", int45);
		strcat(result, reselement);
	}
	if ((Cp_trace(4) != 0. || Cn_trace(4) != 0.) & (Cp_trace(6) != 0. || Cn_trace(6) != 0.)){
		//printf("\n I should not appear\n");
		char reselement[256];
		set_coeffs();
		set_any_Ncoeff(Cp_trace(4), 4, "p");
		set_any_Ncoeff(Cn_trace(4), 4, "n");
		set_any_Ncoeff(Cp_trace(6), 6, "p");
		set_any_Ncoeff(Cn_trace(6), 6, "n");
		double counts46=0.0;
		for (int l = 0; l< num_isos; l++){
			counts46 +=  prefact[l]*prefact[l]*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],32, rhochi, val_difcros);
		}
		double int46 = counts46 - counts[4] - counts[6];
		snprintf(reselement, sizeof reselement, " %.5E", int46);
		strcat(result, reselement);
	}
	if ((Cp_trace(8) != 0. || Cn_trace(8) != 0.) & (Cp_trace(9) != 0. || Cn_trace(9) != 0.)){
		//printf("\n I should appear\n");
		char reselement[256];
		set_coeffs();
		set_any_Ncoeff(Cp_trace(8), 8, "p");
		set_any_Ncoeff(Cn_trace(8), 8, "n");
		set_any_Ncoeff(Cp_trace(9), 9, "p");
		set_any_Ncoeff(Cn_trace(9), 9, "n");
		double counts89=0.0;
		for (int l = 0; l< num_isos; l++){
			counts89 +=  prefact[l]*prefact[l]*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],32, rhochi, val_difcros);
		}
		double int89 = counts89 - counts[8] - counts[9];
		snprintf(reselement, sizeof reselement, " %.5E", int89);
		strcat(result, reselement);
	}

	//printf("here is the result 2 %s\n", result);
	FILE* resultsfile;
	resultsfile = fopen(result_file, "a");
	fprintf(resultsfile,"%s\r\n", result);
	fclose(resultsfile);
}

void write_difrate (char* coefficients, char* path, char* halo_path, double rhochi, void* input_difcros, char* model){
	coeffs_set_master( coefficients, model);
	printf("before reading the halo\n");
	read_halo(halo_path);
	printf("Yes here \n");
	char binned_data[256], header[256];
	//char* result[256];
	snprintf(header, sizeof header, "# Er counts") ;

	//double dE = (E_max - E_th)/(N_bins*1.);

	double CN[16], CP[16], counts[16], i;
	char Op[]= {};
	double endpoint = 3.0;
	double startpoint = -2.0;
	int steps = 100;
	double interval = (endpoint - startpoint)/steps;
	int j;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);

	char header_fin[256];
	snprintf( header_fin, sizeof(header_fin), "%s", write_header(header));

	if (strncmp(target, "Ge", 10) == 0){
		printf("here3\n");
		set_coeffs();
		snprintf(binned_data, sizeof binned_data, "%s/GERMANIUMspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(i=startpoint;i<endpoint;i+=interval){
			write_difrate_line(binned_data, rhochi, val_difcros, i, model);
		}
	}

	if (strncmp(target, "Xe", 10) == 0){
		set_coeffs();
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/XENONspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(i=startpoint;i<endpoint;i+=interval){
			write_difrate_line(binned_data, rhochi, val_difcros, i, model);
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
		set_coeffs();
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/ARGONspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(i=startpoint;i<endpoint;i+=interval){
			write_difrate_line(binned_data, rhochi, val_difcros, i, model);
		}
	}
	if (strncmp(target, "F", 10) == 0){
		set_coeffs();
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/FLUORINEspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(i=startpoint;i<endpoint;i+=interval){
			write_difrate_line(binned_data, rhochi, val_difcros, i, model);
		}
	}

}

void write_difrate_total(char* coefficients, char* path, char* halo_path, double rhochi, void* input_difcros, char* model){
	coeffs_set_master( coefficients, model);
	read_halo(halo_path);

	char binned_data[256], header[256];
	//char* result[256];
	snprintf(header, sizeof header, "# Er counts") ;

	//double dE = (E_max - E_th)/(N_bins*1.);

	double i;
	double endpoint = 3.0;
	double startpoint = -2.0;
	int steps = 100;
	double interval = (endpoint - startpoint)/steps;
	int j;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);

	char header_fin[256];
	snprintf( header_fin, sizeof(header_fin), "Er counts");

	if (strncmp(target, "Ge", 10) == 0){
		snprintf(binned_data, sizeof binned_data, "%s/GERMANIUMspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		for(i=startpoint;i<endpoint;i+=interval){
			fprintf(resultsfile,"%.5E %.5E\r\n",
				pow(10.0,i), difrate_dER(rhochi, val_difcros, i, model));
		}
	}

	if (strncmp(target, "Xe", 10) == 0){
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/XENONspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		for(i=startpoint;i<endpoint;i+=interval){
			fprintf(resultsfile,"%.5E %.5E\r\n",
				pow(10.0,i), difrate_dER( rhochi, val_difcros, i, model));
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/ARGONspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		for(i=startpoint;i<endpoint;i+=interval){
			fprintf(resultsfile,"%.5E %.5E\r\n",
				pow(10.0,i), difrate_dER(rhochi, val_difcros, i, model));
		}
	}
	if (strncmp(target, "F", 10) == 0){
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/FLUORINEspec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		for(i=startpoint;i<endpoint;i+=interval){
			fprintf(resultsfile,"%.5E %.5E\r\n",
				pow(10.0,i), difrate_dER(rhochi, val_difcros, i, model));
		}
	}

	if (strncmp(target, "CaWO4", 10) == 0){
        printf("Yes \n");
		snprintf(binned_data, sizeof binned_data, "%s/CaWO4spec.dat", path);
		FILE* resultsfile;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		for(i=startpoint;i<endpoint;i+=interval){
			fprintf(resultsfile,"%.5E %.5E\r\n",
				pow(10.0,i), difrate_dER(rhochi, val_difcros, i, model));
		}
	}
}

void write_binned_data( char* coefficients, char* exp_path, char* path, char * halo_path, double rhochi, void * input_difcros, char* model){

	coeffs_set_master( coefficients, model);
	read_halo(halo_path);

	char binned_data[256], header[256];
	char experimental_card[256];
	//double dE = (E_max - E_th)/(N_bins*1.);

	double exposure, E1, E2, E_th, E_max;
	int N_bins;
	double counts;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	snprintf(header, sizeof header, "#EL ER counts") ;

	if (strncmp(target, "Ge", 10) == 0){
	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/GERMANIUM.dat", exp_path);
	    FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
		printf("\nopening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
          // Output file name
            snprintf(binned_data, sizeof binned_data, "%s/GERMANIUM.dat", path);
            FILE* resultsfile;
            int j;
            resultsfile = fopen(binned_data, "w+");
			fprintf(resultsfile, "%s \r\n", header);
            for(j=0;j<N_bins;j++){
                    E1 = E_th + j*dE;
                    E2 = E_th + (j+1)*dE;
                    //printf("Yes %.5E \n", E1);
                    counts = counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
                    fprintf(resultsfile, "%.5E %.5E %.5E \r\n", E1, E2, counts);

            }
            fclose(resultsfile);
            printf("Output written to %s\n", binned_data);

	}
	else if (strncmp(target, "Xe", 10) == 0){

	    	//read experiment card
	    	snprintf(experimental_card, sizeof experimental_card, "%s/XENON.dat", exp_path);
	   	FILE* exp_card;
	    	if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    	}
	    	fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    	double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/XENON.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			counts = counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			fprintf(resultsfile, "%.5E %.5E %.5E \r\n", E1, E2, counts);

		}
		fclose(resultsfile);
		printf("Output written to %s\n", binned_data);
	}

	else if (strncmp(target, "Ar", 10) == 0){

	    	//read experiment card
	    	snprintf(experimental_card, sizeof experimental_card, "%s/ARGON.dat", exp_path);
	   	FILE* exp_card;
	    	if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    	}
	    	fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    	double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/ARGON.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			counts = counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			fprintf(resultsfile, "%.5E %.5E %.5E \r\n", E1, E2, counts);

		}
		fclose(resultsfile);
		printf("Output written to %s\n", binned_data);
	}

	else if (strncmp(target, "CaWO4", 10) == 0){

	    	//read experiment card
	    	snprintf(experimental_card, sizeof experimental_card, "%s/CaWO4.dat", exp_path);
	   	FILE* exp_card;
	    	if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    	}
	    	fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    	double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/CaWO4.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			counts = counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			fprintf(resultsfile, "%.5E %.5E %.5E \r\n", E1, E2, counts);

		}
		fclose(resultsfile);
		printf("Output written to %s\n", binned_data);
	}

	else if (strncmp(target, "Xe", 10) != 0 && strncmp(target, "Ge", 10) != 0 && strncmp(target, "Ar", 10) != 0){
		printf("I haven't written the code for this detector yet!\n");
	}

}

void write_binned_rate_line(char* result_file, double rhochi, char* halo_path, void * input_difcros, double E1, double E2, double exposure, char* model){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	char result[256];
	copy_trace();

	// Output file name
	double counts[16];

	counts[0]=0.0;
	char result_tot[256];
	counts[0] +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
	//printf("the result %lf %lf\n",  E1, E2, counts[0]);
	snprintf(result, sizeof result, "%.5E %.5E %.5E", E1, E2, counts[0]);
	for (int k = 1; k<16; k++){
		if(Cp_trace(k) == 0.0 && Cn_trace(k)==0.0){
			//printf("\n this should happen\n");
			continue;
		}
		if(Cp_trace(k) != 0. || Cn_trace(k) !=0.){
			//printf("\n why here %i %f %f", k, CN[k], CP[k] );
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(k), k, "p");
			set_any_Ncoeff(Cn_trace(k), k, "n");
			counts[k]=0.0;
			counts[k] +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			snprintf(reselement, sizeof reselement, " %.5E", counts[k]);
			printf("\n test %s %s \n", result, reselement);
			strcat(result, reselement);
			printf("result %s\n", result);
		}
	}
	if ((Cp_trace(1) != 0. || Cn_trace(1) != 0.) & (Cp_trace(3) != 0. || Cn_trace(3) != 0.)){
		char reselement[256];
		printf("\n I should not appear\n");
		set_coeffs();
		set_any_Ncoeff(Cp_trace(1), 1, "p");
		set_any_Ncoeff(Cn_trace(1), 1, "n");
		set_any_Ncoeff(Cp_trace(3), 3, "p");
		set_any_Ncoeff(Cn_trace(3), 3, "n");
		double counts13=0.0;
		counts13 += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
		double int13 = counts13 - counts[1] - counts[3];
		snprintf(reselement, sizeof reselement, " %.5E", int13);
		strcat(result, reselement);
	}
	if ((Cp_trace(4) != 0. || Cn_trace(4) != 0.) & (Cp_trace(5) != 0. || Cn_trace(5) != 0.)){
		char reselement[256];
		//printf("\n I should not appear\n");
		set_coeffs();
		set_any_Ncoeff(Cp_trace(4), 4, "p");
		set_any_Ncoeff(Cn_trace(4), 4, "n");
		set_any_Ncoeff(Cp_trace(5), 5, "p");
		set_any_Ncoeff(Cn_trace(5), 5, "n");
		double counts45=0.0;
		counts45 +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
		double int45 = counts45 - counts[4] - counts[5];
		snprintf(reselement, sizeof reselement, " %.5E", int45);
		strcat(result, reselement);
	}
	if ((Cp_trace(4) != 0. || Cn_trace(4) != 0.) & (Cp_trace(6) != 0. || Cn_trace(6) != 0.)){
		//printf("\n I should not appear\n");
		char reselement[256];
		set_coeffs();
		set_any_Ncoeff(Cp_trace(4), 4, "p");
		set_any_Ncoeff(Cn_trace(4), 4, "n");
		set_any_Ncoeff(Cp_trace(6), 6, "p");
		set_any_Ncoeff(Cn_trace(6), 6, "n");
		double counts46=0.0;
		counts46 +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
		double int46 = counts46 - counts[4] - counts[6];
		snprintf(reselement, sizeof reselement, " %.5E", int46);
		strcat(result, reselement);
	}
	if ((Cp_trace(8) != 0. || Cn_trace(8) != 0.) & (Cp_trace(9) != 0. || Cn_trace(9) != 0.)){
		//printf("\n I should appear\n");
		char reselement[256];
		set_coeffs();
		set_any_Ncoeff(Cp_trace(8), 8, "p");
		set_any_Ncoeff(Cn_trace(8), 8, "n");
		set_any_Ncoeff(Cp_trace(9), 9, "p");
		set_any_Ncoeff(Cn_trace(9), 9, "n");
		double counts89=0.0;
		counts89 +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
		double int89 = counts89 - counts[8] - counts[9];
		snprintf(reselement, sizeof reselement, " %.5E", int89);
		strcat(result, reselement);
	}

	//printf("here is the result 2 %s\n", result);
	FILE* resultsfile;
	resultsfile = fopen(result_file, "a");
	fprintf(resultsfile,"%s\r\n", result);
	fclose(resultsfile);
}

void write_binned_data_each( char* coefficients, char* exp_path, char* path, char * halo_path, double rhochi, void * input_difcros, char* model){

	coeffs_set_master( coefficients, model);
	read_halo(halo_path);

	char binned_data[256], header[256];
	char experimental_card[256];
	//double dE = (E_max - E_th)/(N_bins*1.);

	double exposure, E1, E2, E_th, E_max;
	int N_bins;
	double counts;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	snprintf(header, sizeof header, "#EL ER counts") ;
	char header_fin[256];
	snprintf( header_fin, sizeof(header_fin), "%s", write_header(header));

	if (strncmp(target, "Ge", 10) == 0){
	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/GERMANIUM.dat", exp_path);
	    FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\nopening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
        // Output file name
        snprintf(binned_data, sizeof binned_data, "%s/GERMANIUM.dat", path);
        FILE* resultsfile;
        int j;
        resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
        for(j=0;j<N_bins;j++){
            E1 = E_th + j*dE;
            E2 = E_th + (j+1)*dE;
            //printf("Yes %.5E \n", E1);
            write_binned_rate_line(binned_data, rhochi, halo_path, input_difcros, E1, E2, exposure, model);
        }
        printf("Output written to %s\n", binned_data);

	}
	else if (strncmp(target, "Xe", 10) == 0){

	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/XENON.dat", exp_path);
	   	FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/XENON.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			write_binned_rate_line(binned_data, rhochi, halo_path, input_difcros, E1, E2, exposure, model);
		}
		printf("Output written to %s\n", binned_data);
	}

	else if (strncmp(target, "Ar", 10) == 0){

	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/ARGON.dat", exp_path);
	   	FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/ARGON.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			write_binned_rate_line(binned_data, rhochi, halo_path, input_difcros, E1, E2, exposure, model);

		}
		printf("Output written to %s\n", binned_data);
	}
	else if (strncmp(target, "Xe", 10) != 0 && strncmp(target, "Ge", 10) != 0 && strncmp(target, "Ar", 10) != 0){
		printf("I haven't written the code for this detector yet!\n");
	}

}

const char* write_header_NP(char* header){
	for (int j=1; j<16; j++){
	//printf("\nhere\n");
		if( Cp(j) != 0.0 || Cn(j) != 0.0) {
			printf("\n in if\n");
			char element[20];
			set_any_coeffs_trace(Cn(j),j, "n");
			set_any_coeffs_trace(Cp(j),j, "p");
			printf("\n %i Couplings %.5E %.5E\n", j, Cn_trace(j), Cp_trace(j));
			if ( Cp(j) != 0.0){
				snprintf(element, sizeof element, " c%ip_c%ip", j,j);
				strcat( header, element);
			}
			if ( Cn(j) != 0.0){
				snprintf(element, sizeof element, " c%in_c%in", j,j);
				strcat( header, element);
			}
			if ( Cn(j) != 0.0 && Cp(j)!=0.0){
				snprintf(element, sizeof element, " c%ip_c%in", j,j);
				strcat( header, element);
			}
		}
		else{
			set_any_coeffs_trace(0.0,j, "n");
			set_any_coeffs_trace(0.0,j, "p");
		}
	}
	if ((Cp(1) != 0. || Cn(1) != 0.) & (Cp(3) != 0. || Cn(3) != 0.)){
		if (Cp(1)!=0. && Cp(3) != 0.){
			char element[10];
			snprintf(element, sizeof element, " c1p_c3p");
			strcat( header, element);
		}
		if (Cn(1) != 0. && Cn(3) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c1n_c3n");
			strcat( header, element);
		}
		if (Cp(1)!=0. && Cn(3) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c1p_c3n");
			strcat( header, element);
		}

		if (Cn(1)!=0. && Cp(3) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c1n_c3p");
			strcat( header, element);
		}

	}
	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(5) != 0. || Cn(5) != 0.)){
		if (Cp(4)!=0. && Cp(5) != 0.){
			char element[10];
			snprintf(element, sizeof element, " c4p_c5p");
			strcat( header, element);
		}
		if (Cn(4) != 0. && Cn(5) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c4n_c5n");
			strcat( header, element);
		}
		if (Cp(4)!=0. && Cn(5) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c4p_c5n");
			strcat( header, element);
		}

		if (Cn(4)!=0. && Cp(5) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c4n_c5p");
			strcat( header, element);
		}
	}
	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(6) != 0. || Cn(6) != 0.)){
		if (Cp(4)!=0. && Cp(6) != 0.){
			char element[10];
			snprintf(element, sizeof element, " c4p_c6p");
			strcat( header, element);
		}
		if (Cn(4) != 0. && Cn(6) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c4n_c6n");
			strcat( header, element);
		}
		if (Cp(4)!=0. && Cn(6) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c4p_c6n");
			strcat( header, element);
		}

		if (Cn(4)!=0. && Cp(6) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c4n_c6p");
			strcat( header, element);
		}
	}

	if ((Cp(8) != 0. || Cn(8) != 0.) & (Cp(9) != 0. || Cn(9) != 0.)){
		if (Cp(8)!=0. && Cp(9) != 0.){
			char element[10];
			snprintf(element, sizeof element, " c8p_c9p");
			strcat( header, element);
		}
		if (Cn(8) != 0. && Cn(9) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c8n_c9n");
			strcat( header, element);
		}
		if (Cp(8)!=0. && Cn(9) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c8p_c9n");
			strcat( header, element);
		}

		if (Cn(8)!=0. && Cp(9) != 0.) {
			char element[10];
			snprintf(element, sizeof element, " c8n_c9p");
			strcat( header, element);
		}
	}
	printf("header %s \n", header);
	return header;
}

void write_binned_rate_line_NP(char* result_file, double rhochi, char* halo_path, void * input_difcros, double E1, double E2, double exposure, char* model){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	char result[256];
	copy_trace();

	// Output file name
	double counts[45];

	counts[0]=0.0;
	char result_tot[256];
	counts[0] +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
	//printf("the result %lf %lf %lf\n",  E1, E2, counts[0]);
	snprintf(result, sizeof result, "%.5E %.5E %.5E", E1, E2, counts[0]);
	for (int k = 1; k<16; k++){
		if(Cp_trace(k) == 0.0 && Cn_trace(k)==0.0){
			//printf("\n this should happen\n");
			continue;
		}
		if(Cp_trace(k) != 0. || Cn_trace(k) !=0.){
			//printf("\n why here %i %f %f %i\n", k, Cp_trace(k), Cn_trace(k));

			if ( Cp_trace(k) != 0.0){
				char reselement[256];
				set_coeffs();
				set_any_Ncoeff(Cp_trace(k), k, "p");
				counts[1+3*(k-1)]=0.0;
				counts[1+3*(k-1)] +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
				snprintf(reselement, sizeof reselement, " %.5E", counts[1+3*(k-1)]);
				//printf("\n test %s %s \n", result, reselement);
				strcat(result, reselement);
				//printf("result %s\n", result);
			}
			if ( Cn_trace(k) != 0.0){
				char reselement[256];
				set_coeffs();
				set_any_Ncoeff(Cn_trace(k), k, "n");
				counts[2+3*(k-1)]=0.0;
				counts[2+3*(k-1)] +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
				snprintf(reselement, sizeof reselement, " %.5E", counts[2+3*(k-1)]);
				//printf("\n test %s %s \n", result, reselement);
				strcat(result, reselement);
				//printf("result %s\n", result);
			}
			if ( Cn_trace(k) != 0.0 && Cp_trace(k)!=0.0){
				char reselement[256];
				set_coeffs();
				set_any_Ncoeff(Cn_trace(k), k, "n");
				set_any_Ncoeff(Cp_trace(k), k, "p");
				counts[3+3*(k-1)]=0.0;
				counts[3+3*(k-1)] +=  counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model) - counts[1+3*(k-1)]- counts[2+3*(k-1)];
				snprintf(reselement, sizeof reselement, " %.5E", counts[3+3*(k-1)]);
				//printf("\n test %s %s \n", result, reselement);
				strcat(result, reselement);
				//printf("result %s\n", result);
			}
		}
	}
	if ((Cp_trace(1) != 0. || Cn_trace(1) != 0.) & (Cp_trace(3) != 0. || Cn_trace(3) != 0.)){
		//printf("\n I should only appear with interference of Ops 1 and 3\n");
		if(Cp_trace(1) != 0. & Cp_trace(3) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(1), 1, "p");
			set_any_Ncoeff(Cp_trace(3), 3, "p");
			double counts1p3p=0.0;
			counts1p3p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int1p3p = counts1p3p - counts[1] - counts[7];
			snprintf(reselement, sizeof reselement, " %.5E", int1p3p);
			strcat(result, reselement);
		}
		if(Cn_trace(1) != 0. & Cn_trace(3) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(1), 1, "n");
			set_any_Ncoeff(Cn_trace(3), 3, "n");
			double counts1n3n=0.0;
			counts1n3n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int1n3n = counts1n3n - counts[2] - counts[8];
			snprintf(reselement, sizeof reselement, " %.5E", int1n3n);
			strcat(result, reselement);
		}
		if(Cn_trace(1) != 0. & Cp_trace(3) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(1), 1, "n");
			set_any_Ncoeff(Cp_trace(3), 3, "p");
			double counts1n3p=0.0;
			counts1n3p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int1n3p = counts1n3p - counts[2] - counts[7];
			snprintf(reselement, sizeof reselement, " %.5E", int1n3p);
			strcat(result, reselement);
		}
		if(Cp_trace(1) != 0. & Cn_trace(3) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(1), 1, "p");
			set_any_Ncoeff(Cn_trace(3), 3, "n");
			double counts1p3n=0.0;
			counts1p3n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int1p3n = counts1p3n - counts[1] - counts[8];
			snprintf(reselement, sizeof reselement, " %.5E", int1p3n);
			strcat(result, reselement);
		}
	}

	if ((Cp_trace(4) != 0. || Cn_trace(4) != 0.) & (Cp_trace(5) != 0. || Cn_trace(5) != 0.)){
		//printf("\n I should only appear with interference of Ops 4 and 5\n");
		if(Cp_trace(4) != 0. & Cp_trace(5) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(4), 4, "p");
			set_any_Ncoeff(Cp_trace(5), 5, "p");
			double counts4p5p=0.0;
			counts4p5p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4p5p = counts4p5p - counts[10] - counts[13];
			//snprintf(reselement, sizeof reselement, " %.5E", int4p5p);
			strcat(result, reselement);
		}
		if(Cn_trace(4) != 0. & Cn_trace(5) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(4), 4, "n");
			set_any_Ncoeff(Cn_trace(5), 5, "n");
			double counts4n5n=0.0;
			counts4n5n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4n5n = counts4n5n - counts[11] - counts[14];
			snprintf(reselement, sizeof reselement, " %.5E", int4n5n);
			strcat(result, reselement);
		}
		if(Cn_trace(4) != 0. & Cp_trace(5) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(4), 4, "n");
			set_any_Ncoeff(Cp_trace(5), 5, "p");
			double counts4n5p=0.0;
			counts4n5p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4n5p = counts4n5p - counts[11] - counts[13];
			snprintf(reselement, sizeof reselement, " %.5E", int4n5p);
			strcat(result, reselement);
		}
		if(Cp_trace(4) != 0. & Cn_trace(5) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(4), 4, "p");
			set_any_Ncoeff(Cn_trace(5), 5, "n");
			double counts4p5n=0.0;
			counts4p5n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4p5n = counts4p5n - counts[10] - counts[14];
			snprintf(reselement, sizeof reselement, " %.5E", int4p5n);
			strcat(result, reselement);
		}
	}
	if ((Cp_trace(4) != 0. || Cn_trace(4) != 0.) & (Cp_trace(6) != 0. || Cn_trace(6) != 0.)){
		//printf("\n I should only appear with interference of Ops 4 and 6\n");
		if(Cp_trace(4) != 0. & Cp_trace(6) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(4), 4, "p");
			set_any_Ncoeff(Cp_trace(6), 6, "p");
			double counts4p6p=0.0;
			counts4p6p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4p6p = counts4p6p - counts[10] - counts[16];
			snprintf(reselement, sizeof reselement, " %.5E", int4p6p);
			strcat(result, reselement);
		}
		if(Cn_trace(4) != 0. & Cn_trace(6) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(4), 4, "n");
			set_any_Ncoeff(Cn_trace(6), 6, "n");
			double counts4n6n=0.0;
			counts4n6n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4n6n = counts4n6n - counts[11] - counts[17];
			snprintf(reselement, sizeof reselement, " %.5E", int4n6n);
			strcat(result, reselement);
		}
		if(Cn_trace(4) != 0. & Cp_trace(6) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(4), 4, "n");
			set_any_Ncoeff(Cp_trace(6), 6, "p");
			double counts4n6p=0.0;
			counts4n6p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4n6p = counts4n6p - counts[11] - counts[16];
			snprintf(reselement, sizeof reselement, " %.5E", int4n6p);
			strcat(result, reselement);
		}
		if(Cp_trace(4) != 0. & Cn_trace(6) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(4), 4, "p");
			set_any_Ncoeff(Cn_trace(6), 6, "n");
			double counts4p6n=0.0;
			counts4p6n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int4p6n = counts4p6n - counts[10] - counts[17];
			snprintf(reselement, sizeof reselement, " %.5E", int4p6n);
			strcat(result, reselement);
		}
	}

	if ((Cp_trace(8) != 0. || Cn_trace(8) != 0.) & (Cp_trace(9) != 0. || Cn_trace(9) != 0.)){
		//printf("\n I should only appear with interference of Ops 8 and 9\n");
		if(Cp_trace(8) != 0. & Cp_trace(9) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(8), 8, "p");
			set_any_Ncoeff(Cp_trace(9), 9, "p");
			double counts8p9p=0.0;
			counts8p9p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int8p9p = counts8p9p - counts[22] - counts[25];
			snprintf(reselement, sizeof reselement, " %.5E", int8p9p);
			strcat(result, reselement);
		}
		if(Cn_trace(8) != 0. & Cn_trace(9) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(8), 8, "n");
			set_any_Ncoeff(Cn_trace(9), 9, "n");
			double counts8n9n=0.0;
			counts8n9n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int8n9n = counts8n9n - counts[23] - counts[26];
			snprintf(reselement, sizeof reselement, " %.5E", int8n9n);
			strcat(result, reselement);
		}
		if(Cn_trace(8) != 0. & Cp_trace(9) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cn_trace(8), 8, "n");
			set_any_Ncoeff(Cp_trace(9), 9, "p");
			double counts8n9p=0.0;
			counts8n9p += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int8n9p = counts8n9p - counts[23] - counts[25];
			snprintf(reselement, sizeof reselement, " %.5E", int8n9p);
			strcat(result, reselement);
		}
		if(Cp_trace(8) != 0. & Cn_trace(9) != 0.){
			char reselement[256];
			set_coeffs();
			set_any_Ncoeff(Cp_trace(8), 8, "p");
			set_any_Ncoeff(Cn_trace(9), 9, "n");
			double counts8p9n=0.0;
			counts8p9n += counts_bin(halo_path, rhochi, input_difcros, exposure, E1, E2, model);
			double int8p9n = counts8p9n - counts[22] - counts[26];
			snprintf(reselement, sizeof reselement, " %.5E", int8p9n);
			strcat(result, reselement);
		}
	}

	//printf("here is the result 2 %s\n", result);
	FILE* resultsfile;
	resultsfile = fopen(result_file, "a");
	fprintf(resultsfile,"%s\r\n", result);
	fclose(resultsfile);
}

void write_binned_data_NP( char* coefficients, char* exp_path, char* path, char * halo_path, double rhochi, void * input_difcros, char* model){

	coeffs_set_master( coefficients, model);
	read_halo(halo_path);

	char binned_data[256], header[256];
	char experimental_card[256];
	//double dE = (E_max - E_th)/(N_bins*1.);

	double exposure, E1, E2, E_th, E_max;
	int N_bins;
	double counts;
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	snprintf(header, sizeof header, "#EL ER counts") ;
	char header_fin[256];
	snprintf( header_fin, sizeof(header_fin), "%s", write_header_NP(header));

	if (strncmp(target, "Ge", 10) == 0){
	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/GERMANIUM.dat", exp_path);
	    FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\nopening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
        // Output file name
        snprintf(binned_data, sizeof binned_data, "%s/GERMANIUM.dat", path);
        FILE* resultsfile;
        int j;
        resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
        for(j=0;j<N_bins;j++){
            E1 = E_th + j*dE;
            E2 = E_th + (j+1)*dE;
            //printf("Yes %.5E \n", E1);
            write_binned_rate_line_NP(binned_data, rhochi, halo_path, input_difcros, E1, E2, exposure, model);
        }
        printf("Output written to %s\n", binned_data);

	}
	else if (strncmp(target, "Xe", 10) == 0){

	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/XENON.dat", exp_path);
	   	FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/XENON.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			write_binned_rate_line_NP(binned_data, rhochi, halo_path, input_difcros, E1, E2, exposure, model);
		}
		printf("Output written to %s\n", binned_data);
	}

	else if (strncmp(target, "Ar", 10) == 0){

	    //read experiment card
	    snprintf(experimental_card, sizeof experimental_card, "%s/ARGON.dat", exp_path);
	   	FILE* exp_card;
	    if (!(exp_card = fopen(experimental_card, "r"))){
			printf("\n opening experiment card failed... make one!\n");
	    }
	    fscanf(exp_card, "%lf %lf %lf %i", &exposure, &E_th, &E_max, &N_bins);
	    double dE = (E_max - E_th)/(N_bins*1.);
		snprintf(binned_data, sizeof binned_data, "%s/ARGON.dat", path);

		FILE* resultsfile;
		int j;
		resultsfile = fopen(binned_data, "w+");
		fprintf(resultsfile, "%s \r\n", header_fin);
		fclose(resultsfile);
		for(j=0;j<N_bins;j++){
			E1 = E_th + j*dE;
			E2 = E_th + (j+1)*dE;
			//printf("Yes %.5E \n", E1);
			write_binned_rate_line_NP(binned_data, rhochi, halo_path, input_difcros, E1, E2, exposure, model);

		}
		printf("Output written to %s\n", binned_data);
	}
	else if (strncmp(target, "Xe", 10) != 0 && strncmp(target, "Ge", 10) != 0 && strncmp(target, "Ar", 10) != 0){
		printf("I haven't written the code for this detector yet!\n");
	}

}