#pragma once

double f_Men( double Er, int A, int Z, double a0,
             double a1, double a2, double a3, double a4, double a5);
double FF40_Men( char* Symbol, char * N, double Er);
double FF128_Men( char * Symbol, char * N, double Er);
double FF129_Men( char * Symbol, char * N, double Er);
double FF130_Men( char * Symbol, char * N, double Er);
double FF131_Men( char * Symbol, char * N, double Er);
double FF132_Men( char * Symbol, char * N, double Er);
double FF134_Men( char * Symbol, char * N, double Er);
double FF136_Men( char * Symbol, char * N, double Er);

double FFMen_40_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_128_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_129_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_130_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_131_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_132_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_134_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FFMen_136_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);

double FormFactMen_v0(int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFactMen_v2(int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
