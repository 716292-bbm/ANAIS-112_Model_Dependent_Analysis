/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
##########################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define VERBOSE 0


/*########################################################################################
Dictionary: Used for "Symbol" and "Symbol2"
- M = M
- Sigma' = Sp
- Sigma''= Spp
- Delta = D
- Phi'' = Ppp
########################################################################################*/


double f_Men( double Er, int A, int Z, double a0, double a1, double a2, double a3, double a4,
			 double a5){
    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=Z*mproton+(A-Z)*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    //double b = bfm/0.197; //oscillator param.
    double b =  1./(sqrt(938.919*1.0e-6*(45*pow(A,-1.0/3.0) - 25*pow(A,-2./3.))));
    double y = pow(q*b,2.)/2.;
	double test = exp(-y/2.0)*(a0+a1*y+a2*pow(y,2.)+a3*pow(y,3.)+a4*pow(y,4.)+a5*pow(y,5.));
	//printf("%f, %f, %f, %f\n", exp(-y/2.0), a0, y, test);

return test;
}

double f_WIMPY( double Er, int A, int Z, double a0, double a1, double a2, double a3, double a4,
			 double a5, double a6, double a7, double a8, double a9, double a10){
    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=Z*mproton+(A-Z)*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    //double b = bfm/0.197; //oscillator param.
    double b =  1./(sqrt(938.919*1.0e-6*(45*pow(A,-1.0/3.0) - 25*pow(A,-2./3.))));
    double y = pow(q*b,2.)/4.;
	double test = exp(-y*2.0)*(a0+a1*y+a2*pow(y,2.)+a3*pow(y,3.)+a4*pow(y,4.)+a5*pow(y,5.) + a6*pow(y,6) + a7*pow(y,7) + a8*pow(y,8) + a9*pow(y,9) + a10*pow(y, 10));
	//printf("%f, %f, %f, %f\n", exp(-y/2.0), a0, y, test);

return test;
}

/*########################################################################################
40 ARGON
########################################################################################*/

double FF40_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 40, 18, 40.,-20.9778,2.41486,-0.0368597,0.0, 0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 40, 18, -4.,3.42422,-0.618209,0.0268957,0.0,0.0);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 40, 18, -4.79093,1.4068,-0.0683192,0.0,0.0, 0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 40, 18, 0.326509,-0.452519,0.0589909,0.0,0.0,0.0);}
	else { return 0.;}
}


else return 0.;

}

/*########################################################################################
128 Xenon
########################################################################################*/

double FF128_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 128, 54, 128.,-126.455,35.82,-3.66991,0.125062,-5.63731E-4 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 128, 54, -20.,29.0588,-11.7104,1.68447,-0.0820044,6.65781E-4);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 128, 54, -25.211, 17.592, -3.46466,0.224722,-0.00353316,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 128, 54, 3.89629,-4.73163,1.48489,-0.140203,0.00344765,0.0);}
	else { return 0.;}
}

else return 0.;

}

/*########################################################################################
129 Xenon
########################################################################################*/

double FF129_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 129, 54, 129.,-128.09,36.4367,-3.75317,0.129553,-6.55816E-4 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 129, 54, -21.,30.6854,-12.3687,1.77928,-0.0868754,7.39474E-4);}
	else { return 0.;}
}
if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 129, 54, -26.1264, 18.4401, -3.64669,0.239379,-0.00399779,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 129, 54, 5.47022,-5.96963,1.7533,-0.160094,0.00387983,0.0);}
	else { return 0.;}
}
else return 0.;
}

double FF129_WIMPY(char * Symbol1, char* Symbol2, char* N1, char* N2, double Er){
    double prefact =  16*M_PI/(2*(1/2) + 1);
    if ( strncmp (Symbol1, "D", 10) == 0 && strncmp (Symbol2, Symbol1, 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.0011820311074370736, - 0.004394673084553479, 0.004939677080500773, - 0.0016591951043697303, 0.00030252674687763364, - 
            0.000058711097141676065, 7.537054290128716e-6, - 5.317774154567649e-7, 6.839858687345737e-8, 0.0, 0.0);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.0018401150865632236, - 0.002472619780561839, 0.0029954093078786064, - 0.0017485667587378354,  0.000811856665941505, -0.00015793665651111697, -1.4445982959217275e-6, 1.793241995856393e-6, 6.839864920594645e-8, 0.0, 0.0);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N2, "+", 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.0014748129622369048, - 0.003732479311881731, 0.0032428502221175724, - 0.002132479031795405,  0.0005643668780295333, -0.000059101554569479804, 6.832830281777335e-6, - 6.307317604952918e-7, -6.839861803969482e-8, 0.0, 0.0);
        }

        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N2, "-", 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.0014748129622369048, - 0.003732479311881731, 0.0032428502221175724, - 0.002132479031795405,  0.0005643668780295333, -0.000059101554569479804, 6.832830281777335e-6, - 6.307317604952918e-7, -6.839861803969482e-8, 0.0, 0.0);
        }

        else {return 0.0;};
        
    }
    if ( strncmp (Symbol1, "Sp", 10) == 0 && strncmp (Symbol2, Symbol1, 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.020808872289191243, -0.13646274750687623, 0.36821805666884777, -0.5255360956521645,  0.426380506625309, -0.19888283928562725, 0.05256093049853532, -0.007362921124572906, 0.0004447576508398376, -3.9111225196870595e-6,9.324872435292655e-9);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.018513885042581883, - 0.12378141258326486, 0.33151263238856177, -0.4595787230968264, 0.35830376879188996, -0.1610290190383875, 0.041462915554454936, -0.005751180220212278, 0.00035166268082252914, -3.4610263745635837e-6, 9.324880933165722e-9);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N2, "+", 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, -0.019627864617116512, 0.1299735688442808, -0.3493496139848393, 0.4916039983714224, -0.39101463905664985, 0.17899534758174646, -0.046682820272984135, 0.00650716109272763, -0.00039549447332930254, 3.6860745496702217e-6, -9.32487668422822e-9);
        }

        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N2, "-", 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, -0.019627864617116512, 0.1299735688442808, -0.3493496139848393, 0.4916039983714224, -0.39101463905664985, 0.17899534758174646, -0.046682820272984135, 0.00650716109272763, -0.00039549447332930254, 3.6860745496702217e-6, -9.32487668422822e-9);
        }

    }

    if ( strncmp (Symbol1, "Spp", 10) == 0 && strncmp (Symbol2, Symbol1, 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, 0.010404436144595613,-0.036544366808683726,0.06435315117845676,-0.06636061394953534,0.04319773517899248, -0.01706358482436624, 0.00405027170634069, -0.0005404221227330367, 0.000032952872479730084, 1.2849348515446646e-7, 1.295121171568424e-10);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74,0.009256942521290935,-0.030494848130468506, 0.051234027256647774, -0.05089830170682882, 0.03238601981523851, -0.012741768322147035, 0.0030733393742958747, -0.0004235714049895079,0.0000273160570250448, 1.1691210783186218e-7, 1.2951223518285722e-10);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N2, "+", 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, -0.00981393230855825, 0.03340000636472022, -0.057450420574708144, 0.058128093534203534, -0.03740107149254411, 0.014743897053985473, -0.0035285222323281308, 0.00047863089552092176, -0.000030005008237028564, 1.2270279913173943e-7, -1.2951217616983635e-10);
        }

        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N2, "-", 10) == 0 ){
            return prefact*f_WIMPY(Er, 129, 74, -0.00981393230855825, 0.03340000636472022, -0.057450420574708144, 0.058128093534203534, -0.03740107149254411, 0.014743897053985473, -0.0035285222323281308, 0.00047863089552092176, -0.000030005008237028564, 1.2270279913173943e-7, -1.2951217616983635e-10);
        }
    }

    if (strncmp ( Symbol1, "Sp", 10) == 0 && strncmp (Symbol2, "D", 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10 ) == 0 ){
            return prefact * f_WIMPY(Er, 129, 74, -0.004959509487490609, 
             0.02548149208769477, -0.04924247373152658, 0.044204217456558045, 
             -0.018907732881088244, 0.004161062710306827, -0.0005713348788979224,
             0.00006870842906480167, -5.394486037493884e-6, 
             2.5254862845584265e-8, 0.0);
        }
        
        if (strncmp ( N1, "-", 10) == 0 && strncmp (N1, N2, 10) == 0 ){
            return prefact * f_WIMPY(Er, 129, 74, 0.005836752451299621, 
             -0.023433394287691554, 0.036186086760354196, -0.03191890688184716, 
              0.018400710623717113, -0.00595820810978778,  
              0.0008768656769386072, - 0.00002255412440479619, 
             -4.355746488114697e-6, 2.5254885860657158e-8, 0.0);
        }

        if (strncmp (N1, "-", 10) == 0 && strncmp(N2, "+", 10) == 0 ){
            return prefact * f_WIMPY( Er, 129, 74, -0.0007573308995633562,
            0.0000626295604509357, 0.003882782425347979, - 0.0030477405118422753, -0.0014400560697404767, 0.001152942372012883, 0.000025920840621028315, - 0.00009748538710262024, - 0.000025956434470508036, 
            + 0.00003616093400834701, - 0.000013518274847402887 );
        }

        if (strncmp (N1, "+", 10) == 0 && strncmp(N2, "-", 10) == 0){
            return prefact * f_WIMPY(Er, 129, 74, -0.006187949566189774, 
             0.02444748348692596, - 0.03875562702832236, 0.03455930974221026,
             - 0.020260108044808798, 0.006707850234682113, -0.000997842266798669,
             0.000024945423288469897, 4.96525424878263e-6, -2.5254874353118093e-8, 0.0);
        }
    }
else return 0.0;

}

/*########################################################################################
130 Xenon
########################################################################################*/
double FF130_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 130, 54, 130.,-129.753,37.2381,-3.89291,0.139778,-9.30032E-4 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 130, 54, -22.,32.2019,-13.1152,1.90775,-0.0948184,8.47975E-4);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 130, 54, -27.7106,19.7108, -3.85805,0.252667,-0.00444209,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 130, 54, 6.28519,-6.63842,1.85406,-0.166079,0.00413453,0.0);}
	else { return 0.;}
}

else return 0.;
}


/*########################################################################################
131 Xenon
########################################################################################*/

double FF131_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 131, 54, 131.,-131.26,37.8232,-3.97171,0.142995,-9.12955E-4 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 131, 54, -23.,33.7021,-13.7433,2.00031,-0.0991364,8.60686E-4);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 131, 54, -28.0443, 20.0888, -3.94934, 0.260624, -0.00468846,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 131, 54, 6.90542, -7.17962, 1.97217, -0.175248, 0.00437613, 0.0);}
	else { return 0.;}
}

else return 0.;
}


double FF131_WIMPY(char * Symbol1, char* Symbol2, char* N1, char* N2, double Er){
    double prefact =  16*M_PI/(2*(3/2) + 1);
    if ( strncmp (Symbol1, "D", 10) == 0 && strncmp (Symbol2, Symbol1, 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, 0.06411286031263148, - 0.1629545548505601, 0.19709458160535825, - 0.1325758204184924,  0.05128514880137004, - 0.010332337385168932, 0.0008906474066539766, 0.000012720818128338462, 5.019144903677493e-8, 0.0, 0.0);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, 0.028508020457524867, - 0.07447644460876486, 0.09838919993393962, - 0.07292392714216149, 0.03155489196208029, - 0.007004702220659644,0.0006579409769396572, - 
            0.000010874867430961905,  5.0191080015468504e-8, 0.0, 0.0);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N2, "+", 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, -0.04275196759662569, 0.11017515921130767, - 0.13938377026013277,  0.09855265772071016, -0.04028701628009529, 0.008513050482966399,- 0.0007657792399722185, 0.000011797839364032412, -5.019126452577812e-8, 0.0, 0.0);
        }

        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N2, "-", 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, -0.04275196759662569, 0.11017515921130767, - 0.13938377026013277,  0.09855265772071016, -0.04028701628009529, 0.008513050482966399,- 0.0007657792399722185, 0.000011797839364032412, -5.019126452577812e-8, 0.0, 0.0);
        }

        else {return 0.0;};
    }

    if ( strncmp (Symbol1, "Sp", 10) == 0 && strncmp (Symbol2, Symbol1, 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, 0.014707821388042745, -0.13732537753381133, 0.47309062668703133, - 0.7781456806574102, 0.7126868955838259, - 0.37244667575862184, 0.10845443880671564, - 0.016137538638583637, 0.0009330377654740855, 5.08149166686281e-6, 7.015771540521167e-9 );
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, 0.013272091746950735, -0.12505298374885687, 0.44292909076165154, -0.7661782395006832, 0.7373288715284319,-0.4002496665646047, 0.11947776284330879, -0.0180035776961701, 0.0010450837366378113, 5.380611474527931e-6, 7.015720311671177e-9);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N2, "+", 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, -0.013971526575856616, 0.13104695249126264, -0.45776368441557164, 0.7724000950817769, -0.7251181061534898, 0.386159363063301, -0.11383865586358587,0.01704522062905004,-0.0009874486522898913,-5.231052109747236e-6,-7.015745926048718e-9);
        }

        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N2, "-", 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, -0.013971526575856616, 0.13104695249126264, -0.45776368441557164, 0.7724000950817769, -0.7251181061534898, 0.386159363063301, -0.11383865586358587,0.01704522062905004,-0.0009874486522898913,-5.231052109747236e-6,-7.015745926048718e-9);
        }
    }

    if ( strncmp (Symbol1, "Spp", 10) == 0 && strncmp (Symbol2, Symbol1, 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, 0.007353910694021366,0.024013554202436335,-0.02014650487125114, -0.037184108285885614, 0.10207295648693686, -0.08818794415698591, 0.03593113383888596, -0.007020155839698539, 0.0005402993168705856, -6.167594109512812e-7, 1.787035062316364e-10);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N1, N2, 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, 0.0066360458734753686, 0.0230091769076141, -0.0167703486318344, -0.036938300365609425, 0.09593468456183336, -0.08255616478254374, 0.03381938601451815, -0.006662512003029585, 0.0005185483277899781, -6.047848694091502e-7, 1.7870236298441898e-10);
        }
        if (strncmp ( N1, "-", 10) == 0 && strncmp ( N2, "+", 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, -0.006985763287928304, -0.02351657691879264, 0.018478701186953755, 0.03696179519459378, -0.09886450319704412, 0.08528791992678614, -0.03485107661177335, 0.006838186074602294, -0.0005292745155866782, 6.107721235197574e-7, 1.787029346070819e-10);
        }

        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N2, "-", 10) == 0 ){
            return prefact*f_WIMPY(Er, 131, 74, -0.006985763287928304, -0.02351657691879264, 0.018478701186953755, 0.03696179519459378, -0.09886450319704412, 0.08528791992678614, -0.03485107661177335, 0.006838186074602294, -0.0005292745155866782, 6.107721235197574e-7, 1.787029346070819e-10);;
        }
    }

    if (strncmp ( Symbol1, "Sp", 10) == 0 && strncmp (Symbol2, "D", 10) == 0 ){
        if (strncmp ( N1, "+", 10) == 0 && strncmp ( N1, N2, 10 ) == 0 ){
            return prefact * f_WIMPY(Er, 131, 74, 0.03070766188029818,
             - 0.18238166875057954, 0.3638294378190443,
             - 0.3790471186386412, 0.2250586156968798,
             - 0.07512192860735903, 0.012931718514154415,
             - 0.0009148005098435476, 4.362060624958545e-6,
              1.8757261217738304e-8, 0.0);
        }
        
        if (strncmp ( N1, "-", 10) == 0 && strncmp (N1, N2, 10) == 0 ){
            return prefact * f_WIMPY(Er, 131, 74, 0.01945150541825024,
             - 0.11704675619987992, 0.24537358064222317,
             - 0.2736722082374217, 0.17541123903932163,
             - 0.06270119580834876, 0.0113877960442417,
             - 0.0008427934001863716, 5.105039770518882e-6,
              1.8757123768142482e-8, 0.0);
        }

        if (strncmp (N1, "-", 10) == 0 && strncmp(N2, "+", 10) == 0 ){
            return prefact * f_WIMPY( Er, 131, 74, -0.02478334389448439,
             0.14797579402377817, - 0.30204350108235556,
             0.3236531604737859, - 0.19797413282902399,
             0.06896764198720394, - 0.0149928673789674,
             0.004490422206553065, - 0.00249909264032691, 
             0.0011070951729064053, - 0.00032913985114829176 );
        }

        if (strncmp (N1, "+", 10) == 0 && strncmp(N2, "-", 10) == 0){
            return prefact * f_WIMPY(Er, 131, 74,-0.020476593296150224, 
             0.12234117376193052, - 0.248891931123951,
             0.2712732536833822, - 0.1703408873962704,
             0.05994062216364239, - 0.01078825392947655, 
             0.0007959353218256179, - 4.706373145642153e-6,
             -1.875719249281273e-8, 0.0);
        }
    }
else return 0.0;

}
/*########################################################################################
132 Xenon
########################################################################################*/

double FF132_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 132, 54, 132.,-132.835,38.4665,-4.06999,0.149636,-0.00111463 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 132, 54, -24.,35.253,-14.4437,2.11305,-0.105689,9.61344E-4);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 132, 54, -28.7972, 20.7751, -4.0995, 0.272865, -0.00507527,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 132, 54, 7.93145, -8.01086, 2.12817, -0.186148, 0.00469887, 0.0);}
	else { return 0.;}
}

else return 0.;
}

/*########################################################################################
134 Xenon
########################################################################################*/

double FF134_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 134, 54, 134.,-135.861,39.6872,-4.24713,0.159053,-0.00125724 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 134, 54, -26.,38.2701,-15.773,2.32061,-0.116557,0.00106693);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 134, 54, -29.5095, 21.5578, -4.27308, 0.287393, -0.00555437,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 134, 54, 9.3351, -9.20279, 2.35489, -0.202364, 0.00519463, 0.0);}
	else { return 0.;}
}

else return 0.;
}

/*########################################################################################
136 Xenon
########################################################################################*/
double FF136_Men( char * Symbol, char * N, double Er){
if (strncmp (Symbol,"M",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 136, 54, 136.,-138.787,40.9048,-4.41984,0.165388,-0.00109211 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 136, 54, -28.,41.2081,-17.0848,2.52635,-0.12686,0.00110965);}
	else { return 0.;}
}

if (strncmp (Symbol,"Ppp",10) == 0){
	if (strncmp (N,"+",10) == 0){
		return f_Men(Er, 136, 54, -29.8571, 22.0402, -4.37033, 0.296134, -0.0059684,0.0 );}
	if (strncmp (N,"-",10) == 0){
		return f_Men(Er, 136, 54, 10.1433, -9.96123, 2.48784, -0.212062, 0.00559688, 0.0);}
	else { return 0.;}
}

else return 0.;
}

/*########################################################################################
40 Argon
########################################################################################*/
double FFMen_40_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=22*mproton+18*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF40_Men( "M", N1, Er) * FF40_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
     return pow(q,4.)/(4.*mN*mN*mN*mN)*FF40_Men( "Ppp", N1, Er) * FF40_Men( "Ppp", N2, Er) ;
    }

    if (i == 5 && j == i){
     return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF40_Men( "M", N1, Er)*FF40_Men( "M", N2, Er));
    }


    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF40_Men( "M", N1, Er)*FF40_Men( "M", N2, Er));
    }



    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF40_Men( "M", N1, Er)*FF40_Men( "M", N2, Er);
    }

    if (i == 12 && j == i){
     return Cj/16.*pow(q,2.)/(mN*mN) *FF40_Men( "Ppp", N1, Er)*FF40_Men( "Ppp", N2, Er);
    }
    
    if (i == 15 && j == i){
     printf("O15\n");
     return Cj/16.*pow(q,6.)/(mN*mN*mN*mN*mN*mN)*FF40_Men( "Ppp", N1, Er)*FF40_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     return pow(q,2.)/(2.*mN*mN)*FF40_Men( "M", N1, Er)*FF40_Men( "Ppp", N2, Er);
    }
    
    if (i == 11 && j == 12){
     return Cj/4.*pow(q,2.)/(2.*mN*mN)*FF40_Men( "M", N1, Er)*FF40_Men( "Ppp", N2, Er);
    }
    
    if (i == 11 && j == 15){
     return Cj/4.*(-pow(q,4.)/(2.*mN*mN*mN*mN)*FF40_Men( "M", N1, Er)*FF40_Men( "Ppp", N2, Er));
    }
    
    if (i == 15 && j == 12) {
      printf("O15 O12 interference\n");
      return Cj/16.*(-2) *pow(q,4.)/(mN*mN*mN*mN)* FF40_Men( "Ppp", N1, Er)*FF40_Men( "Ppp", N2, Er);
    }
    else return 0.0;
}

 double FFMen_40_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

 double mproton=0.938; //Proton mass in GeV
 double mneutron=0.940; //Neutron mass in GeV
 double mn=18*mproton+22*mneutron;
 double mN = mproton;
 double q = sqrt(2.*mn*Er*1.e-6);//GeV
 double Cj = 4.*jchi*(jchi+1.)/3.;


 if (i == 5 && j == i){

 return Cj/4.*pow(q,2.)/(mN*mN)*FF40_Men( "M", N1, Er)*FF40_Men( "M", N2, Er);

 }

 if (i == 8 && j == i){
 return Cj/4.*FF40_Men( "M", N1, Er)*FF40_Men( "M", N2, Er);

 }

 else return 0.0;
 }



/*########################################################################################
128 Xenon
########################################################################################*/

double FFMen_128_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+74*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF128_Men( "M", N1, Er) * FF128_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
        return pow(q,4.)/(4.*mN*mN*mN*mN)*FF128_Men( "Ppp", N1, Er) * FF128_Men( "Ppp", N2, Er) ;
    }


    if (i == 5 && j == i){
        return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF128_Men( "M", N1, Er)*FF128_Men( "M", N2, Er));
    }


    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF128_Men( "M", N1, Er)*FF128_Men( "M", N2, Er));
    }



    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF128_Men( "M", N1, Er)*FF128_Men( "M", N2, Er);
    }

    if (i == 15 && j == i){
     return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*FF128_Men( "Ppp", N1, Er)*FF128_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     //return pow(q,2.)/(2.*mN*mN)*FF128_Men( "M", N1, Er)*FF128_Men( "Ppp", N2, Er);
    }
else return 0.0;
}

double FFMen_128_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+74*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF128_Men( "M", N1, Er)*FF128_Men( "M", N2, Er);

}

if (i == 8 && j == i){
return Cj/4.*FF128_Men( "M", N1, Er)*FF128_Men( "M", N2, Er);

}

else return 0.0;
}


/*########################################################################################
129 Xenon
########################################################################################*/

double FFMen_129_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+75*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double q2 = 2.*mn*Er*1.e-6;//GeV
    double q4=q2*q2;
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF129_Men( "M", N1, Er) * FF129_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
        return pow(q,4.)/(4.*mN*mN*mN*mN)*FF129_Men( "Ppp", N1, Er) * FF129_Men( "Ppp", N2, Er)-q4/(4.*muT*muT*mN*mN)*FF129_WIMPY( "Sp", "Sp", N1, N2, Er) ;
    }

    if (i == 4 && j == i){
    return Cj/16.*(FF129_WIMPY( "Spp", "Spp", N1, N2, Er)+FF129_WIMPY( "Sp", "Sp", N1, N2, Er));

    }

    if (i == 5 && j == i){
        return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF129_Men( "M", N1, Er)*FF129_Men( "M", N2, Er)+q4/(mN*mN*mN*mN)*FF129_WIMPY( "D", "D", N1, N2, Er));
    }
    if (i == 6 && j == i){
    return Cj/16.*q4/pow(mN,4.)*FF129_WIMPY( "Spp", "Spp", N1, N2, Er);

    }

    if (i == 7 && j == i){
    return 1./8.*(-q2/(4.*muT*muT)*FF129_WIMPY( "Sp", "Sp", N1, N2, Er));

    }

    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF129_Men( "M", N1, Er)*FF129_Men( "M", N2, Er)+q2/(mN*mN)*FF129_WIMPY( "D", "D", N1, N2, Er));
    }

    if (i == 9 && j == i){
    return Cj/16.*q2/(mN*mN)*FF129_WIMPY( "Sp", "Sp", N1, N2, Er);

    }

    if (i == 10 && j == i){
    return 1./4.*q2/(mN*mN)*FF129_WIMPY( "Spp", "Spp", N1, N2, Er);

    }


    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF129_Men( "M", N1, Er)*FF129_Men( "M", N2, Er);
    }

    if (i == 15 && j == i){
     return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*FF129_Men( "Ppp", N1, Er)*FF129_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     return pow(q,2.)/(2.*mN*mN)*FF129_Men( "M", N1, Er)*FF129_Men( "Ppp", N2, Er);
    }

    if (i == 8 && j == 9){
    return Cj * q2/(8.*mN*mN)*FF129_WIMPY( "Sp", "D", N1, N2, Er);

    }
else return 0.0;
}

double FFMen_129_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+75*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 3 && j == i){
return pow(q,2.)/(mN*mN)*FF129_WIMPY( "Sp", "Sp", N1, N2, Er);

}

if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF129_Men( "M", N1, Er)*FF129_Men( "M", N2, Er);

}

if (i == 7 && j == i){
return 1./8.*FF129_WIMPY( "Sp", "Sp", N1, N2, Er);

}
if (i == 8 && j == i){
return Cj/4.*FF129_Men( "M", N1, Er)*FF129_Men( "M", N2, Er);

}
if (i == 12 && j == i){
return Cj/16.*(FF129_WIMPY( "Spp", "Spp", N1, N2, Er)+0.5*FF129_WIMPY( "Sp", "Sp", N1, N2, Er));

}

if (i == 13 && j == i){
return Cj/16.*pow(q,2.)/(mN*mN)*FF129_WIMPY( "Spp", "Spp", N1, N2, Er);

}

if (i == 14 && j == i){
return Cj/16.*pow(q,2.)/(2*mN*mN)*FF129_WIMPY( "Sp", "Sp", N1, N2, Er);

}

if (i == 15 && j == i){
return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*0.5*FF129_WIMPY( "Sp", "Sp", N1, N2, Er);
}

else return 0.0;
}



/*########################################################################################
130 Xenon
########################################################################################*/

double FFMen_130_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+76*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF130_Men( "M", N1, Er) * FF130_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
        return pow(q,4.)/(4.*mN*mN*mN*mN)*FF130_Men( "Ppp", N1, Er) * FF130_Men( "Ppp", N2, Er) ;
    }

    if (i == 5 && j == i){
        return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF130_Men( "M", N1, Er)*FF130_Men( "M", N2, Er));
    }


    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF130_Men( "M", N1, Er)*FF130_Men( "M", N2, Er));
    }



    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF130_Men( "M", N1, Er)*FF130_Men( "M", N2, Er);
    }

    if (i == 15 && j == i){
     return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*FF130_Men( "Ppp", N1, Er)*FF130_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     return pow(q,2.)/(2.*mN*mN)*FF130_Men( "M", N1, Er)*FF130_Men( "Ppp", N2, Er);
    }
else return 0.0;
}

double FFMen_130_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+76*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;


if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF130_Men( "M", N1, Er)*FF130_Men( "M", N2, Er);

}

if (i == 8 && j == i){
return Cj/4.*FF130_Men( "M", N1, Er)*FF130_Men( "M", N2, Er);

}

else return 0.0;
}


/*########################################################################################
131 Xenon
########################################################################################*/

double FFMen_131_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+77*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double q2 = 2.*mn*Er*1.e-6;//GeV
    double q4=q2*q2;
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF131_Men( "M", N1, Er) * FF131_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
        return pow(q,4.)/(4.*mN*mN*mN*mN)*FF131_Men( "Ppp", N1, Er) * FF131_Men( "Ppp", N2, Er) -q4/(4.*muT*muT*mN*mN)*FF131_WIMPY( "Sp", "Sp", N1, N2, Er) ;
    }

    if (i == 4 && j == i){
    return Cj/16.*(FF131_WIMPY( "Spp", "Spp", N1, N2, Er)+FF131_WIMPY( "Sp", "Sp", N1, N2, Er));

    }

    if (i == 5 && j == i){
        return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF131_Men( "M", N1, Er)*FF131_Men( "M", N2, Er)+q4/(mN*mN*mN*mN)*FF131_WIMPY( "D", "D", N1, N2, Er));
    }

    if (i == 6 && j == i){
    return Cj/16.*q4/pow(mN,4.)*FF131_WIMPY( "Spp", "Spp", N1, N2, Er);

    }

    if (i == 7 && j == i){
        return 1./8.*(-q2/(4.*muT*muT)*FF131_WIMPY( "Sp", "Sp", N1, N2, Er));

    }

    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF131_Men( "M", N1, Er)*FF131_Men( "M", N2, Er) + q2/(mN*mN)*FF131_WIMPY( "D", "D", N1, N2, Er));
    }

    if (i == 9 && j == i){
    return Cj/16.*q2/(mN*mN)*FF131_WIMPY( "Sp", "Sp", N1, N2, Er);

    }

    if (i == 10 && j == i){
    return 1./4.*q2/(mN*mN)*FF131_WIMPY( "Spp", "Spp", N1, N2, Er);

    }

    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF131_Men( "M", N1, Er)*FF131_Men( "M", N2, Er);
    }


    if (i == 15 && j == i){
     return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*FF131_Men( "Ppp", N1, Er)*FF131_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     return pow(q,2.)/(2.*mN*mN)*FF131_Men( "M", N1, Er)*FF131_Men( "Ppp", N2, Er);
    }

    if (i == 8 && j == 9){
    return Cj * q2/(8.*mN*mN)*FF131_WIMPY( "Sp", "D", N1, N2, Er);

    }
else return 0.0;
}

double FFMen_131_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+77*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 3 && j == i){
return pow(q,2.)/(mN*mN)*FF131_WIMPY( "Sp", "Sp", N1, N2, Er);

}

if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF131_Men( "M", N1, Er)*FF131_Men( "M", N2, Er);

}

if (i == 7 && j == i){
return 1./8.*FF131_WIMPY("Sp", "Sp", N1, N2, Er);

}

if (i == 8 && j == i){
return Cj/4.*FF131_Men( "M", N1, Er)*FF131_Men( "M", N2, Er);

}

if (i == 12 && j == i){
return Cj/16.*(FF131_WIMPY( "Spp", "Spp", N1, N2, Er)+0.5*FF131_WIMPY( "Sp", "Sp", N1, N2, Er));

}
if (i == 13 && j == i){
return Cj/16.*pow(q,2.)/(mN*mN)*FF131_WIMPY( "Spp", "Spp", N1, N2, Er);

}
if (i == 14 && j == i){
return Cj/16.*pow(q,2.)/(2*mN*mN)*FF131_WIMPY( "Sp", "Sp", N1, N2, Er);

}

if (i == 15 && j == i){
return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*0.5*FF131_WIMPY( "Sp", "Sp", N1, N2, Er);
}

else return 0.0;
}

/*########################################################################################
132 Xenon
########################################################################################*/

double FFMen_132_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+78*mneutron;
    double q2 = 2.*mn*Er*1.e-6;//GeV
    double q4=q2*q2;
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 1 && j == i){
return FF132_Men( "M", N1, Er)*FF132_Men( "M", N2, Er);
}
if (i == 3 && j == i){
return q4/(4.*mN*mN*mN*mN)*FF132_Men( "Ppp", N1, Er)*FF132_Men( "Ppp", N1, Er);

}

if (i == 5 && j == i){
return Cj/4.*(-q4/(4.*muT*muT*mN*mN)*FF132_Men( "M", N1, Er)*FF132_Men( "M", N2, Er));

}

if (i == 8 && j == i){
return Cj/4.*(-q2/(4.*muT*muT)*FF132_Men( "M", N1, Er)*FF132_Men( "M", N2, Er));

}


if (i == 11 && j == i){
return Cj/4.*q2/(mN*mN)*FF132_Men( "M", N1, Er)*FF132_Men( "M", N2, Er);

}

if (i == 15 && j == i){
return Cj/16.*q4/(mN*mN*mN*mN)*FF132_Men( "Ppp", N1, Er)*FF132_Men( "Ppp", N1, Er);

}

if (i == 1 && j == 3){
return q2/(2.*mN*mN)*FF132_Men( "M", N1, Er)*FF132_Men( "Ppp", N1, Er);

}

else return 0.0;
}

double FFMen_132_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+78*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF132_Men( "M", N1, Er)*FF132_Men( "M", N2, Er);

}

if (i == 8 && j == i){
return Cj/4.*FF132_Men( "M", N1, Er)*FF132_Men( "M", N2, Er);

}

else return 0.0;
}

/*########################################################################################
134 Xenon
########################################################################################*/

double FFMen_134_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+80*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF134_Men( "M", N1, Er) * FF134_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
        return pow(q,4.)/(4.*mN*mN*mN*mN)*FF134_Men( "Ppp", N1, Er) * FF134_Men( "Ppp", N2, Er) ;
    }

    if (i == 5 && j == i){
        return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF134_Men( "M", N1, Er)*FF134_Men( "M", N2, Er));
    }


    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF134_Men( "M", N1, Er)*FF134_Men( "M", N2, Er));
    }



    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF134_Men( "M", N1, Er)*FF134_Men( "M", N2, Er);
    }

    if (i == 15 && j == i){
     return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*FF134_Men( "Ppp", N1, Er)*FF134_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     return pow(q,2.)/(2.*mN*mN)*FF134_Men( "M", N1, Er)*FF134_Men( "Ppp", N2, Er);
    }
else return 0.0;
}

double FFMen_134_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+80*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF134_Men( "M", N1, Er)*FF134_Men( "M", N2, Er);

}

if (i == 8 && j == i){
return Cj/4.*FF134_Men( "M", N1, Er)*FF134_Men( "M", N2, Er);

}

else return 0.0;
}


/*########################################################################################
136 Xenon
########################################################################################*/

double FFMen_136_v0(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+82*mneutron;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double mN = mproton;
    double muT = mn*mchi/(mn+mchi);
    double Cj = 4.*jchi*(jchi+1.)/3.;

    if (i == 1 && j == i){
        return (FF136_Men( "M", N1, Er) * FF136_Men( "M", N2, Er));
    }

    if (i == 3 && j == i){
        return pow(q,4.)/(4.*mN*mN*mN*mN)*FF136_Men( "Ppp", N1, Er) * FF136_Men( "Ppp", N2, Er) ;
    }

    if (i == 5 && j == i){
        return Cj/4.*(-pow(q,4.)/(4.*muT*muT*mN*mN)*FF136_Men( "M", N1, Er)*FF136_Men( "M", N2, Er));
    }


    if (i == 8 && j == i){
     return Cj/4.*(-pow(q,2.)/(4.*muT*muT)*FF136_Men( "M", N1, Er)*FF136_Men( "M", N2, Er));
    }



    if (i == 11 && j == i){
     return Cj/4.*pow(q,2.)/(mN*mN)*FF136_Men( "M", N1, Er)*FF136_Men( "M", N2, Er);
    }

    if (i == 15 && j == i){
     return Cj/16.*pow(q,4.)/(mN*mN*mN*mN)*FF136_Men( "Ppp", N1, Er)*FF136_Men( "Ppp", N2, Er);
    }

    if (i == 1 && j == 3){
     return pow(q,2.)/(2.*mN*mN)*FF136_Men( "M", N1, Er)*FF136_Men( "Ppp", N2, Er);
    }
else return 0.0;
}

double FFMen_136_v2(int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

    double mproton=0.938; //Proton mass in GeV
    double mneutron=0.940; //Neutron mass in GeV
    double mn=54*mproton+82*mneutron;
	double mN = mproton;
    double q = sqrt(2.*mn*Er*1.e-6);//GeV
    double Cj = 4.*jchi*(jchi+1.)/3.;

if (i == 5 && j == i){
return Cj/4.*pow(q,2.)/(mN*mN)*FF136_Men( "M", N1, Er)*FF136_Men( "M", N2, Er);

}


if (i == 8 && j == i){
return Cj/4.*FF136_Men( "M", N1, Er)*FF136_Men( "M", N2, Er);

}

else return 0.0;
}


double FormFactMen_v0(int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

if ( A==40){
    return FFMen_40_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==128){
return FFMen_128_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==129){
return FFMen_129_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==130){
return FFMen_130_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==131){
return FFMen_131_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==132){
return FFMen_132_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==134){
return FFMen_134_v0(i, j, N1, N2, Er, mchi, jchi);
}

if (A==136){
return FFMen_136_v0(i, j, N1, N2, Er, mchi, jchi);
}

else {
printf("The isotope, %i chosen is not implemented\n", A);
return 0.;
}

}

double FormFactMen_v2(int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi){

if ( A==40){
    return FFMen_40_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==128){
return FFMen_128_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==129){
return FFMen_129_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==130){
return FFMen_130_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==131){
return FFMen_131_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==132){
return FFMen_132_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==134){
return FFMen_134_v2(i, j, N1, N2, Er, mchi, jchi);
}

if (A==136){
return FFMen_136_v2(i, j, N1, N2, Er, mchi, jchi);
}

else {
printf("The isotope, %i chosen is not implemented\n", A);
return 0.;
}

}
