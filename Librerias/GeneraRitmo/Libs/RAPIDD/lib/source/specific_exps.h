#pragma once
double bin_Xenon1T( double rhochi, double mass, double E1, double E2, char* model, char* basis);
double bin_DS50( double rhochi, double mass, double E1, double E2, char* model, char* basis);
double bin_LZ( double rhochi, double mass, double E1, double E2, char* model, char* basis);
double bin_DS20k( double rhochi, double mass, double E1, double E2, char* model, char* basis);
double res_fn_lux_nr(int A, int Z, double ERnr);
double counts_effres_bin_Xenon1T( double rhochi, double mass, double E1, double E2, char*model, char*basis);
double counts_effres_bin_LZ( double rhochi, double mass, double E1, double E2, char*model, char*basis);
void read_DS50_LEFF(char* path);
double DS50_LEFF ( double E_r);
double res_fn_DS20_nr(double ERnr);
double counts_effres_bin_DS50( double rhochi, double mass, double E1, double E2, char*model, char*basis);
double counts_effres_bin_DS20k( double rhochi, double mass, double E1, double E2, char*model, char*basis);