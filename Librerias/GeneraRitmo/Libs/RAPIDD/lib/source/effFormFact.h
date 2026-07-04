#pragma once

double myFF131(int procid, int pnid, double Er);
static size_t MEGA=0;
double fsi(int A, double Er);
double f(double Er, int A, int Z, double a0, double a1, double a2, double a3, double a4, 
         double a5, double a6, double a7, double a8, double a9, double a10);

double fholger(double Er, int A, int Z, double const * coeffs);
double FF19( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF23( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF28( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF70( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF72( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF73( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF74( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF76( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF127( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF128( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF129( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF130( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF131( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF132( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF134( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);
double FF136( char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);

double FormFact_19_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_19_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);

double FormFact_23_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_23_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_28_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_28_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);

double FormFact_70_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_70_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_72_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_72_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_73_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_73_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_74_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_74_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_76_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_76_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);

double FormFact_127_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_127_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);

double FormFact_128_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_128_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_129_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_129_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_130_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_130_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_131_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_131_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_132_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_132_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_132_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_134_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_136_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_136_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);

double FormFact_v0(int A,int Z,int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_v2(int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_gen_SI(int A, char* symbol1, char* symbol2, char * N1, char * N2, double Er);
