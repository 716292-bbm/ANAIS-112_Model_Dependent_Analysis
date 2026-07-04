#include "TString.h"
#include "TGaxis.h"
#include "TImage.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TEventList.h"
#include "TRandom.h"
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <RooFitResult.h>
#include <TMatrixDSym.h>

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooChi2Var.h"
#include "RooHistPdf.h"
#include "RooLinkedList.h"
#include "RooMinimizer.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimWSTool.h"
#include "RooSimultaneous.h"
#include "RooFormulaVar.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TText.h>

#include <ADB.h>
#include <DMRatePlotHandle.hh>

using namespace RooFit;
using namespace std;

#define ARCHIVO

//#define MIGDAL

////////////////////////////////////
////////////////////////////////////
// MARIA: GLOBALS, Eee
double minEne = 0;
double maxEne = 100; // MARIA extend to 30 keVee (100 keV NR for QFNa=0.3)
double binEne = 0.1;
////////////////////////////////////
// GLOBALS: QF
// if TGraph != null, use it, otherwise, use cte
TGraph *gQNa = 0;
TGraph *gQI = 0;
double QNa = 0;
double QI = 0;

int Conv2(double e1, double e2, double ebin, double p1, double ps,
          double *S0, TH1 *ritmo_sr)
{
  double PRECISION = 0.0001; // Maximum precision
  double minEne = 0.001;
  int ind;
  double ei, edif, en;                                 // loop variables
  double ef, estep, en1, en2, ens;                     // interval variables
  double S01, SM1, PHI1, NORM1, S02, SM2, PHI2, NORM2; // Integral var.
  double weight;
  int indArray;

  // Number of bins in the auxiliar array
  int arrayDim = 400;
  int nSig = 5; // Number of +-sigma of the convolution interval

  // Compute auxiliar array with the rates. FFm e1-nSig*sigma (or 0.1) to e2+nSig*sigma
  // double auxsig = (ResolutionK<0 ? -nSig*ResolutionK : nSig*ResolutionK*sqrt(e1));
  double auxsig = nSig * fabs(p1 + ps * sqrt(e1));
  double arrayEi = (e1 - auxsig < minEne ? minEne : e1 - auxsig);
  double arrayEf = e2;
  // if (ResolutionK<0) arrayEf -= nSig * ResolutionK;
  // else arrayEf += nSig * ResolutionK * sqrt (e2);
  arrayEf += nSig * fabs(p1 + ps * sqrt(e2));
  double arrayEbin = (arrayEf - arrayEi) / arrayDim;

  // Auxiliar arrays with the rates (Unmodulated, modulated, phase)
  double *array_S0 = new double[arrayDim];

  ei = arrayEi;
  for (ind = 0; ind < arrayDim; ind++)
  {
    ei += arrayEbin;
    array_S0[ind] = ritmo_sr->GetBinContent(ritmo_sr->FindBin(ei));
  }

  // Loop ei in the number of output bins
  ind = 0;
  for (ei = e1; ei < e2; ei += ebin)
  {
    // PRECISION CONTROL
    if (fabs(e2 - ei) < PRECISION)
      break;

    ef = ei + ebin;
    estep = (ef - ei) / 20.;

    // initialize integrating variables
    S01 = 0;
    NORM1 = 0;

    // Loop edif: integral in the bin ei-ef
    for (edif = ei; edif <= ef; edif += estep)
    {
      // auxsig = (ResolutionK<0 ? -nSig*ResolutionK : nSig*ResolutionK*sqrt(edif));
      auxsig = nSig * fabs(p1 + ps * sqrt(edif));
      en1 = (edif - auxsig < minEne ? minEne : edif - auxsig);
      en2 = edif + auxsig;
      ens = (en2 - en1) / 100.;
      // Initialize integrating variables
      S02 = 0;
      NORM2 = 0;

      // Control of 0 resolution
      if (auxsig == 0)
      {
        indArray = (int)((edif - arrayEi) / arrayEbin);
        S02 = array_S0[indArray];
        NORM2 = 1;
        en2 = 0;
      }

      // Loop en: Gaussian convolution
      for (en = en1; en <= en2; en += ens)
      {
        // auxsig = (ResolutionK<0 ? -ResolutionK : ResolutionK*sqrt(en));
        auxsig = fabs(p1 + ps * sqrt(en));
        weight = exp(-(edif - en) * (edif - en) / 2. / auxsig / auxsig) / auxsig;

        // Look for the index corresponding to en in the rate array
        indArray = (int)((en - arrayEi) / arrayEbin);

        // Integrate
        // MARIA: Mas adelante se puede hacer una integracion trapezoidal...
        S02 += weight * array_S0[indArray];
        NORM2 += weight;

      } // end of en loop

      // Integrate
      S01 += S02 / NORM2;
      NORM1++;

    } // end of edif loop

    // Fill the output arrays
    S0[ind] = S01 / NORM1;
    ind++;

  } // End of ei loop

  delete[] array_S0;

  return 0;
}

// TODO
TH1F *DMModelGetRate(double mw, double sigma, int qfModel, int SpinModel) // Devuelve un puntero a un histograma de ROOT
{
#ifdef ARCHIVO

  TFile *f = new TFile("rate_DMAnalysis_SDn_Al.root", "UPDATE");
  // f->mkdir("SI_rates")->cd();
  TString file_name = "rate_SDn";

#endif

  DMRate *rate = new DMRate();                    // Se crea un objeto tipo DMRate
  int err = rate->Initialize((char *)"rate.dat"); // Inicializa con los valores del archivo rate.dat
  if (err != 0)
    return 0;
  rate->SetMW(mw); // Fija la masa del WIMP

  if (SpinModel == 0) // Spin-Independent
  {
    rate->SetSigSI(sigma); // Fija la seccion eficaz SI
    rate->SetSigSD(0);     // Fija la seccion eficaz SD
    rate->SetTheta(0);     // Fija theta en 0
  }

  if (SpinModel == 1) // Spin-Dependent Proton
  {
    rate->SetSigSI(0);     // Fija la seccion eficaz SI
    rate->SetSigSD(sigma); // Fija la seccion eficaz SD
    for (int iel = 0; iel < rate->GetNElements(); iel++)
    {
      rate->GetElement(iel)->SetTheta(0); // Fija theta en 0
    }
  }

  if (SpinModel == 2) // Spin-Dependent Neutron
  {
    rate->SetSigSI(0);     // Fija la seccion eficaz SI
    rate->SetSigSD(sigma); // Fija la seccion eficaz SD
    for (int iel = 0; iel < rate->GetNElements(); iel++)
    {
      rate->GetElement(iel)->SetTheta(1.57079632679); // Fija theta en pi/2
    }
  }

  if (gQNa)
    rate->GetElement(0)->SetREF(gQNa);
  else
    rate->GetElement(0)->SetREF(QNa);

  DMRatePlotHandle plothl(rate);

  TH1F *hrate = plothl.GetRate(1000); // Calcula los valores del histograma

#ifdef ARCHIVO
  TString histName = file_name + Form("_mw%0.1f", mw);
  TString histTitle = Form("DMA Rate SDn - M_{W} = %0.1f GeV;Energy (keV);Rate (events/kg/day/keV)", mw);

  hrate->SetName(histName);
  hrate->SetTitle(histTitle);

  // Escribir el histograma en el archivo
  hrate->Write();
  hrate->SetDirectory(0); // Desasocia el histograma del archivo

  f->cd();
  f->Close();
#endif
  return hrate;
}

// FORM WIMPY AND RAPIDD
TH1D *DMModelGetRateEeeFromFile(std::string fileName, std::string nameNa, std::string nameI)
{
  TFile *f = TFile::Open(fileName.c_str(), "READ");
  if (!f || f->IsZombie() || !f->IsOpen())
  {
    std::cerr << "[DMModelGetRateEeeFromFile] ERROR: no se pudo abrir '"
              << fileName.c_str() << "'\n";
    return nullptr;
  }

  TH1D *hNa = static_cast<TH1D *>(f->Get(nameNa.c_str()));
  hNa->SetDirectory(0);
  if (!hNa || hNa->IsZombie())
  {
    std::cerr << "[DMModelGetRateEeeFromFile] ERROR: histograma '"
              << nameNa.c_str() << "' no encontrado\n";
    f->Close();
    delete f;
    return nullptr;
  }
  TH1D *hI = static_cast<TH1D *>(f->Get(nameI.c_str()));
  hI->SetDirectory(0);
  if (!hI || hI->IsZombie())
  {
    std::cerr << "[DMModelGetRateEeeFromFile] ERROR: histograma '"
              << nameI.c_str() << "' no encontrado\n";
    f->Close();
    delete f;
    return nullptr;
  }

  // create new histogram till EMax
  int nBins = (maxEne - minEne) / binEne;
  TH1D *hWimp = new TH1D(Form("%s_%s", nameNa.c_str(), nameI.c_str()), "", nBins, minEne, maxEne);
  hWimp->SetDirectory(0);
  for (int ii = 1; ii <= nBins; ii++)
  {
    double Eee = hWimp->GetBinCenter(ii);
    // Na
    double qNa = QNa;
    if (gQNa)
      qNa = gQNa->Eval(Eee);
    double ENR_Na = Eee / qNa;
    double rateNa = hNa->GetBinContent(hNa->FindBin(ENR_Na)) / qNa;

    // I
    double qI = QI;
    if (gQI)
      qI = gQI->Eval(Eee);
    double ENR_I = Eee / qI;
    double rateI = hI->GetBinContent(hI->FindBin(ENR_I)) / qI;

    hWimp->SetBinContent(ii, (23 * rateNa + 127 * rateI) / 150.);
  }

  f->Close();
  delete f;

  return hWimp;
}

// arg1: eneIni
// arg2: eneEnd
// arg3: (optional) include ALE population (identified with ANOD) (default false)
int main(int argc, char **argv)
{
  // std::cout << " argc" << argc << std::endl;
  // std::cout << " argv" << argv << std::endl;

  if (argc < 5)
  {
    std::cout << "Usage fitSimulMakeExclusion cl eneIni eneEnd spinModel[0-SI, 1-SDp, 2-SD-neutron] qf [optional, default Tamara, else use 1 for DAMA and 2 for cte 0.2 0.06] includeALE[optional, default: false]" << std::endl;
    return 1;
  }

  // Se lee desde el string que se le pasa a la funcion los valores:
  double cl = atof(argv[1]);  // Establece el Confidence Level
  double min = atof(argv[2]); // Establece el nivel de energia minimo
  double max = atof(argv[3]); // Establece el nivel de energia maximo
  // int SpinModel = atof(argv[4]); // Lee el modelo de Spin (0-SI 1-SD-proton 2-SD-neutron)

  // BORJA Esto es para cambiar entre SPIN-DEPENDENT y SPIN-INDEPENDENT
  // BORJA TODO - GENERALIZAR A OPERADORES TEORIA EFECTIVA

  int SpinModel = 2;

  int resolution_p = 2; // 2 tiene en cuenta la resolución (convolucion) , 1 tiene en cuenta la resolucion (Generacion de eventos segun gausiana), 0 no la tiene en cuenta

  double param_1, param_2;

  if (resolution_p == 2)
  {

    const char *lowResFile = "fitsResolution.root";
    TFile *fres = TFile::Open(lowResFile, "READ");
    TF1 *ffres = (TF1 *)fres->Get("fresD0");

    param_1 = ffres->GetParameter(0);
    param_2 = ffres->GetParameter(1);
  }

  // int thmodel = 1; // 0-DMAnalysis, 1-Python, 2-RAPIDD, 3-WIMPYDD 4-DMAnalysis(Por Archivo) 5-Migdall

  int thmodel = atof(argv[4]);

  std::cout << " cl" << cl << std::endl;
  std::cout << " min" << min << std::endl;
  std::cout << " max" << max << std::endl;
  std::cout << " SpinModel" << SpinModel << std::endl;

  if (SpinModel < 0 || SpinModel > 2)
  {
    std::cout << " Spin Model not Valid!!" << std::endl;
    exit(0);
  }

  int qfModel = 1; // Establece el modelo de quenching factor
  std::cout << " SpinModel: (1-DAMA, 2-ANAIS CTE, 3-TAMARA)" << SpinModel << std::endl;
  if (argc > 5)
    qfModel = atoi(argv[5]); // Si se le pasa por los parametros de funcion cambia el modelo de QF

  // ANOD SÍ O NO
  bool ANOD = 1;

  if (argc > 6)
    ANOD = 1; // Si se le pasa por parametros incluye ALE population

  // MARIA 100326. Set QF MODE
  if (qfModel == 1) // DAMA
  {
    QNa = 0.3;
    QI = 0.09;
  }
  else if (qfModel == 2) // ANAIS CTE
  {
    // QNa = 0.2;
    // QI = 0.06;
    QNa = 1.0;
    QI = 1.0;
  }
  else // TAMARA
  {
    TFile *f = new TFile("QFTamara.root", "READ");
    if (!f)
    {
      std::cout << " CANNOT READ QUENCHING FACTOR FILE!!" << std::endl;
      exit(0);
    }
    gQNa = (TGraph *)f->Get("gNa");
    if (!gQNa)
    {
      std::cout << " CANNOT FIND gNa IN QUENCHING FACTOR FILE!!" << std::endl;
      exit(0);
    }
    gQI = (TGraph *)f->Get("gI");
    if (!gQI)
    {
      std::cout << " CANNOT FIND gI IN QUENCHING FACTOR FILE!!" << std::endl;
      exit(0);
    }
  }

  // RooWorkspace* w = new RooWorkspace("w", "wimpFitSimul Workspace");
  // Crea una variable real de RooFit llamada energia:

  RooRealVar energy("energy", "energy", minEne, maxEne, "keV");

  // w->import(energy);
  // CHANGE HERE FOR DIFFERENT YEARS - Seleccionamos la exposicion
  // CHECK COSINE : 3 years, 50kg
  std::vector<int> years = {1, 2, 3, 4, 5, 6};
  // std::vector<int> years = {1,2,3};
  int nyears = years.size(); // Define el numero de anios
  std::vector<int> detectors = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  // std::vector<int> detectors = {2,4,6,7,8}; // use the cleanest 5
  int ndet = detectors.size(); // Define el numero de detectores

  RooCategory detCat("det", "det"); // detector category for the simultaneuos fit
  for (int det = 0; det < ndet; det++)
    detCat.defineType(Form("det%d", detectors[det]));
  // w->import(detCat);
  detCat.Print("V");

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // READ DATA
  //////////////////////////////////////////////////////////////////////////

  // READ ANAIS EXPOSURE
  ADBTime dt;
  // livetime is different for every detector
  // TODO, by now copy it from Ivan mail for 6 years

  std::vector<double> exposure; // Genera un vector dinamico para las exposiciones

  // Live time D0: 2031.38 days
  // Live time D1: 2033.20 days
  // Live time D2: 2029.52 days
  // Live time D3: 2022.55 days
  // Live time D4: 2033.01 days
  // Live time D5: 2030.18 days
  // Live time D6: 2032.27 days
  // Live time D7: 2031.02 days
  // Live time D8: 2020.29 days
  //  CHECK COSINE : 3 years, 60kg

  // double scaleFactor= 8./6.; <---- CON ESTO SE ESCALA LA EXPOSICION

  double scaleFactor = 1;

  // Agrega al vector dinamico exposure los siguientes valores:
  exposure.push_back(2031.38 * 12.5 * scaleFactor); // D0
  exposure.push_back(2033.20 * 12.5 * scaleFactor); // D1
  exposure.push_back(2029.52 * 12.5 * scaleFactor); // D2
  exposure.push_back(2022.55 * 12.5 * scaleFactor); // D3
  exposure.push_back(2033.01 * 12.5 * scaleFactor); // D4
  exposure.push_back(2030.18 * 12.5 * scaleFactor); // D5
  exposure.push_back(2032.27 * 12.5 * scaleFactor); // D6
  exposure.push_back(2031.02 * 12.5 * scaleFactor); // D7
  exposure.push_back(2020.29 * 12.5 * scaleFactor); // D8
  // for(int i=0;i<nyears;i++) exposure.push_back(dt.GetExposure(-1,years[i],1)/ndet); // exposure in kgxday

  // READ ANAIS BKG

  // Genera el string del nombre del archivo en base a los anios seleccionados anteriormente
  std::string name = "../../data/BEhistos_year";
  for (int i = 0; i < nyears; i++)
    name += std::to_string(years[i]);
  name += ".root";

  // Abre el archivo root con el nombre del string anterior

  TFile *file = new TFile(name.c_str(), "read");

  // MAP OF HISTOGRAMS FOR THE SIMULTANEOUS FIT - Carga los histogramas del archivo ROOT
  std::map<std::string, TH1 *> hdataMap;
  for (int det = 0; det < ndet; det++)
  {
    std::string histname = "hbea_";
    for (int i = 0; i < nyears; i++)
      histname += std::to_string(years[i]);
    histname += "y_D";
    histname += std::to_string(detectors[det]);
    TH1F *h = (TH1F *)file->Get(histname.c_str());
    if (!h)
    {
      std::cout << " file " << histname.c_str() << " not found in " << name.c_str() << " file " << std::endl;
      continue;
    }
    h->Scale(h->GetBinWidth(1) * exposure[det]); // change to counts for the extended fit
    // RooDataHist* dataDet = new RooDataHist(Form("dataD%d",detectors[det]),"",RooArgSet(energy),h);
    // w->import(*dataDet);
    hdataMap[Form("det%d", detectors[det])] = h;
  }

  // create RooDataHist from the map - Toma cada histograma en hdataMap y lo convierte en dataset de ROOFIT
  RooDataHist *dataData = new RooDataHist("data", "data", RooArgSet(energy), detCat, hdataMap);
  // w->import(*dataData);

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // READ BACKGROUND MODEL
  ///////////////////////////////////////////////////////////////////////////
  // single hits
  std::string bkgname = "../../backgroundModel/backgroundModel_single_y123456";
  // for (int i=0; i<nyears;i++) bkgname+=std::to_string(years[i]);
  bkgname += (ANOD ? "_conANOD.root" : ".root");    // Con esto carga si tenemos el ALE o no
  TFile *fbkg = new TFile(bkgname.c_str(), "read"); // Abre el archivo que tiene el modelo de fondo

  // Scale factor (cross-section)
  // here it is a normalization factor
  // -> large range, positive or negative in order to not bias the fit
  // common for all detectors
  RooRealVar *nNorm = new RooRealVar("nNorm", "", 1, -1e10, 1e10);
  nNorm->setRange(-1e6, 1e6);
  nNorm->setVal(1e4);
  nNorm->setConstant(kFALSE); // Marca la variable como libre en el ajuste

  ////////////////////////////////
  // LOOP IN WIMP MASSES

#ifdef MIGDAL
                              // PRUEBA CARMEN MIGDAL
  std::vector<double> mw = {
      0.25, 0.29, 0.33, 0.37, 0.43, 0.49, 0.56, 0.64,
      0.73, 0.84, 0.96, 1.1, 1.3, 1.4, 1.6, 1.9,
      2.1, 2.4, 2.8, 3.2, 3.7, 4.2, 4.8, 5.5,
      6.3, 7.1, 8.2, 9.3, 11.0, 12.0, 14.0, 16.0,
      18.0, 21.0, 24.0, 27.0, 31.0, 36.0, 41.0, 47.0,
      53.0, 61.0, 70.0, 80.0, 91.0, 100.0, 120.0, 140.0,
      160.0, 180.0, 200.0, 230.0, 270.0, 310.0, 350.0, 400.0,
      460.0, 520.0, 600.0, 680.0, 780.0, 890.0, 1000.0, 1200.0,
      1300.0, 1500.0, 1700.0, 2000.0, 2300.0, 2600.0, 3000.0, 3400.0,
      3900.0, 4500.0, 5100.0, 5800.0, 6700.0, 7600.0, 8700.0, 10000.0};

#else

  // BORJA - Funcion que genera el array de masas
  std::vector<double> mw;

  for (int imw = 2; imw <= 9; ++imw)
  {
    for (int j = 0; j < 4; ++j)
    {
      double a = imw + 0.2 * j;
      mw.push_back(a);
    }
  }

  for (int imw = 1; imw <= 9; ++imw)
  {
    for (int j = 0; j < 4; ++j)
    {
      double a = imw + 0.2 * j;
      mw.push_back(a * 10);
    }
  }

  for (int imw = 1; imw <= 9; ++imw)
  {
    for (int j = 0; j < 4; ++j)
    {
      double a = imw + 0.2 * j;
      mw.push_back(a * 100);
    }
  }

  for (int imw = 1; imw <= 10; ++imw)
  {
    for (int j = 0; j < 4; ++j)
    {
      double a = imw + 0.2 * j;
      mw.push_back(a * 1000);
    }
  }

#endif

  // BORJA - Genera dos arrays con los nombres de los histogramas

  std::vector<std::string> nombres_Na;
  std::vector<std::string> nombres_I;
  nombres_Na.reserve(mw.size());
  nombres_I.reserve(mw.size());

#ifdef MIGDAL
  // CARMEN: MIGDAL
  std::ostringstream oss;
  std::vector<std::string> mw_str = {
      "0.25", "0.29", "0.33", "0.37", "0.43", "0.49", "0.56", "0.64",
      "0.73", "0.84", "0.96", "1.1", "1.3", "1.4", "1.6", "1.9",
      "2.1", "2.4", "2.8", "3.2", "3.7", "4.2", "4.8", "5.5",
      "6.3", "7.1", "8.2", "9.3", "11.0", "12.0", "14.0", "16.0",
      "18.0", "21.0", "24.0", "27.0", "31.0", "36.0", "41.0", "47.0",
      "53.0", "61.0", "70.0", "80.0", "91.0", "100.0", "120.0", "140.0",
      "160.0", "180.0", "200.0", "230.0", "270.0", "310.0", "350.0", "400.0",
      "460.0", "520.0", "600.0", "680.0", "780.0", "890.0", "1000.0", "1200.0",
      "1300.0", "1500.0", "1700.0", "2000.0", "2300.0", "2600.0", "3000.0", "3400.0",
      "3900.0", "4500.0", "5100.0", "5800.0", "6700.0", "7600.0", "8700.0", "10000.0"};

  for (size_t i = 0; i < mw.size(); ++i)
  {
    nombres_Na.push_back("hist_Na_mw_" + mw_str[i]);
    nombres_I.push_back("hist_I_mw_" + mw_str[i]);
  }

#else
  std::ostringstream oss;

  for (double m : mw)
  {
    oss.str("");
    oss.clear();
    oss << std::fixed << std::setprecision(1) << m;
    std::string mass_str = oss.str();

    nombres_Na.push_back("hist_Na_mw_" + mass_str);
    nombres_I.push_back("hist_I_mw_" + mass_str);
  }

#endif

  // MARIA: CHOSE HERE
  std::string fileName;

  if (thmodel == 0)
  {
    fileName = "RAPIDD_SI_TH1D.root";
  }
  if (thmodel == 1)
  {
    fileName = "ANAIS_SI_TH1D.root";
  }
  if (thmodel == 2)
  {
    fileName = "RAPIDD_SI_TH1D.root";
  }
  if (thmodel == 3)
  {
    if (SpinModel == 0)
      fileName = "WIMPYDD_SI_TH1D.root";
    if (SpinModel == 1)
      fileName = "WIMPYDD_SDp_TH1D.root";
    if (SpinModel == 2)
      fileName = "WIMPYDD_SDn_TH1D.root";
  }
  if (thmodel == 4)
  {
    fileName = "DMA_SI_TH1D.root";
  }

  if (thmodel == 5)
  {
    fileName = "MIGDAL.root";
  }

  std::vector<double> sigma;
  for (int imw = 0; imw < (int)mw.size(); imw++)
  {
    /////////////////////////////////////////
    // simultaneous fit
    RooSimultaneous *simWimpHyp = new RooSimultaneous("simWimpHyp", "simultaneuos pdf WimpHyp", detCat);

    // Read histogram with DM rate for mW and cross_section=1 - CALCULA EL ESPECTRO TEORICO
    TH1 *hWimpNoRes;
    if (thmodel == 0)
    {
      hWimpNoRes = DMModelGetRate(mw[imw], 1, qfModel, SpinModel);
    }
    else
    {
      hWimpNoRes = DMModelGetRateEeeFromFile(fileName, nombres_Na[imw], nombres_I[imw]);
    }
    // TH1F *hWimp = (TH1F *)DMModelGetRate(mw[imw], 1, qfModel, SpinModel);
    // TH1D *hWimp = (TH1D *)DMModelGetRate_table_clone("file.root", hist_mw[imw]); // histogram in NR energy

    // hWimp in electron-equivalent

    // Apply resolution in hWimp
    int nBins = (maxEne - minEne) / binEne;
    const int nBins_2 = 1000;
    double rate_res[nBins_2];
    double integral_hWimpNoRes = hWimpNoRes->Integral(1, nBins);

    TH1 *hWimp = (TH1 *)hWimpNoRes->Clone(Form("%s_res", hWimpNoRes->GetName()));

    if (resolution_p == 1)
    {
      hWimp->Reset();
      TRandom ran;
      int nEvRes = 10000000;

      // CARMEN 26/03/2026: GET fRes
      const char *lowResFile = "/media/storage2/tamara/resMartaByDet/fitsResolution.root";
      TFile *fres = TFile::Open(lowResFile, "READ");
      TF1 *ffres = (TF1 *)fres->Get("fresD0");

      for (int ii = 0; ii < nEvRes; ii++)
      {
        double ee = hWimpNoRes->GetRandom();
        // double sigma = fres->Eval(ee);
        double sigma = ffres->GetParameter(0) + ffres->GetParameter(1) * sqrt(ee);
        double eeRes = ran.Gaus(ee, sigma);
        hWimp->Fill(eeRes);
      }

      double integral_hWimp = hWimp->Integral(1, nBins);
      double factor_norm = integral_hWimpNoRes / integral_hWimp;
      hWimp->Scale(factor_norm);
    }

    if (resolution_p == 2)
    {
      hWimp->Reset();

      int asd = Conv2(0, 100, 0.1, param_1, param_2, rate_res, hWimpNoRes);

      for (int ii = 1; ii <= nBins; ii++)
      {
        hWimp->SetBinContent(ii, rate_res[ii - 1]);
      }
    }

    RooDataHist *rdhWimp = new RooDataHist("rateWimp", "", RooArgSet(energy), hWimp);
    RooHistPdf *pdf_Wimp = new RooHistPdf("pdf_Wimp", "", energy, energy, *rdhWimp, 1);

    // Construct the model for every detector
    // nbkg* fbk + nNorm*fwimp
    for (int det = 0; det < ndet; det++)
    {
      //  BACKGROUND MODEL
      // The model for every detector is : nNorm*pdf_Wimp + nbkg*pdf_bkg

      // Lee el histograma del fondo del detector
      TH1F *hbkg = (TH1F *)fbkg->Get(Form("hD%d", detectors[det]));

      hbkg->Smooth(); // smootheamos para quitar fluctuaciones (opcional)

      // Comandos ROOFIT
      RooDataHist *rdhBkg = new RooDataHist(Form("rdhBkg_D%d", detectors[det]), "", RooArgSet(energy), hbkg);
      RooHistPdf *pdf_bkg = new RooHistPdf(Form("pdf_bkg_D%d", detectors[det]), "", energy, energy, *rdhBkg, 1);

      // Convierte a tasa de numero de cuentas por energia (Exposicion por Tasa de Ritmo)
      double fscale = exposure[detectors[det]] * hbkg->GetBinWidth(1); // Normalization factor : kgxdayxbinWidth

      // Integra el fondo en el rango de energia del fit
      double intebkg = hbkg->Integral(hbkg->FindBin(min), hbkg->FindBin(max) - 1) * fscale; // in counts

      // Creamos variable del numero de cuentas total del fondo nbkg
      RooRealVar *nbkg = new RooRealVar(Form("nbkg_D%d", detectors[det]), "", 1000, intebkg / 2, intebkg * 2);
      nbkg->setRange(intebkg / 2, intebkg * 2); // Establecemos rango
      nbkg->setVal(intebkg);                    // inicializamos
      nbkg->setConstant(kTRUE);                 // FIJAMOS el parametro (el fit no lo toca)

      // Hace el fit de n_norm*pdf_Wimp + n_bkg*pdf_bkg

      RooAddPdf *model_det = new RooAddPdf(Form("model_det%d", detectors[det]), "WIMP + background",
                                           RooArgList(*pdf_Wimp, *pdf_bkg),
                                           RooArgList(*nNorm, *nbkg));
      simWimpHyp->addPdf(*model_det, Form("det%d", detectors[det]));
    }

    ///////////////////////////////////////
    ///////////////////////////////////////
    ///////////////////////////////////////
    // FIT
    ///////////////////////////////////////
    double sig; // Aqui se guarda la seccion eficaz limite de la masa del WIMP
    if (min > 0)
    {
      nNorm->setVal(1e4);                    // Parametro libre
      energy.setRange("fitRange", min, max); // Establecemos rango de Energia para el FIT

      // Hace el ajuste
      simWimpHyp->fitTo(*dataData, SumCoefRange("fitRange"), RooFit::Range("fitRange"), SumW2Error(false), Save(true), Extended());

      // uncomment to print results
      // RooFitResult *results = simWimpHyp->fitTo(*dataData, SumCoefRange("fitRange"), RooFit::Range("fitRange"), SumW2Error(false), Save(true), Extended());
      // RooFitResult* results= simWimpHyp->fitTo(*dataData, Save(true), Extended());
      // results->Print();

      ////////// Get sigma
      double nNorm_val = nNorm->getVal();
      double nNorm_err = nNorm->getError();

      // Aplica el Confidence Level al valor calculado
      double valcl = (cl == 95 ? nNorm_val + nNorm_err * 1.64 : nNorm_val + nNorm_err * 1.28);
      std::cout << " counts at 90\% CL " << valcl << std::endl;

      // Calcula sigma limite dividiendo el valor con el confidence level entre el valor teorico del RITMO obtenido con la simulacion del fondo
      sig = valcl / hWimp->Integral(hWimp->FindBin(min), hWimp->FindBin(max), "width");
    }
    else
    {
      // MARIA: ESTO NO FUNCIONA
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
      sig = 1e300;
      double minInt = 4; // 2 keV
      double step = 1;   // step 0.5 keV
      for (double eneIni = 1; eneIni <= 6 - minInt; eneIni += step)
      {
        for (double eneEnd = 1 + minInt; eneEnd <= 6; eneEnd += step)
        {
          nNorm->setVal(1e4);
          energy.setRange("fitRange", eneIni, eneEnd);
          simWimpHyp->fitTo(*dataData, SumCoefRange("fitRange"), RooFit::Range("fitRange"), SumW2Error(false), Save(true), Extended());
          double nNorm_val = nNorm->getVal();
          double nNorm_err = nNorm->getError();
          double valcl = (cl == 95 ? nNorm_val + nNorm_err * 1.64 : nNorm_val + nNorm_err * 1.28);
          double integral = hWimp->Integral(hWimp->FindBin(eneIni), hWimp->FindBin(eneEnd), "width");
          std::cout << " ENE: " << eneIni << " - " << eneEnd << " counts at 90\% CL " << valcl << " integral " << integral << std::endl;
          if (valcl > 0 && valcl / integral < sig)
            sig = valcl / integral;
        }
      }
    }

    double totalExposure = 0;
    for (int det = 0; det < ndet; det++)
      totalExposure += exposure[det];
    std::cout << " total exposure " << totalExposure << " kgxday" << std::endl;

    // Pasamos de pb a cm^2 y al dividir por totalExposure estamos convirtiendo el Ritmo de la simulacion en numero de cuentas.
    sig *= 1e-36 / totalExposure;

    // std::cout << " mw: "<< mw[imw] << " sigma: " << sig << " cm2 " << std::endl;

    // Anadimos sig al array de sigmas
    sigma.push_back(sig);

    delete simWimpHyp;
    delete rdhWimp;
    delete pdf_Wimp;
    delete hWimp;
  }
  // end loop in mw

  // Muestra los valores de mW y de sigma por pantalla
  for (int imw = 0; imw < (int)mw.size(); imw++)
    std::cout << " mw: " << mw[imw] << " sigma: " << sigma[imw] << " cm2 " << std::endl;

  // TGraph representa graficos XY: TGraph(Numero de puntos, Array x, Array y)
  TGraph *gex = new TGraph(mw.size(), mw.data(), sigma.data());

  // Abre el archivo de root results.root
  //  TFile * f = new TFile("plots/results.root","recreate");
  TFile *f = new TFile("plots/SI_varios_2.root", "update");

  // Creamos una variable para el nombre del objeto con los datos
  std::string gname = "gA112_6y_";
  // std::string gname = "gCOSINE_3y_";
  // std::string gname = "gCOSINE_1y_";
  gname += (cl == 90 ? "90_" : "95_");
  if (min > 0)
  {
    gname += std::to_string(min);
    gname += "_";
    gname += std::to_string(max);
  }
  else
    gname += "opt";
  gname += "_QF";
  if (qfModel == 1)
    gname += "dama";
  else if (qfModel == 2)
    gname += "cte";
  else
    gname += "tamara";
  if (ANOD)
    gname += "_ANOD";

  if (thmodel == 0)
  {
    gname += "_DM_SG";
  }
  if (thmodel == 1)
  {
    gname += "_PY";
  }
  if (thmodel == 2)
  {
    gname += "_RA";
  }
  if (thmodel == 3)
  {
    gname += "_WI";
  }

  if (thmodel == 4)
  {
    gname += "_DMA";
  }

  if (thmodel == 4)
  {
    gname += "_MIGDALL";
  }

  double pb2cm = 1e-36;
  double mL = 1;
  double mH = 1e4;
  double sL = 5e-50;
  double sH = 1e-37;

  if (SpinModel == 0) // Spin-Independent
  {
    gname += "_SI";
  }

  if (SpinModel == 1) // Spin-Dependent Proton
  {
    gname += "_SDproton";
    sL = 5e-40;
    sH = 1e-33;
  }

  if (SpinModel == 2) // Spin-Dependent Neutron
  {
    gname += "_SDneutron";
    sL = 5e-40;
    sH = 1e-33;
  }

  gex->Write(gname.c_str());
  f->Close();

  TCanvas *c = new TCanvas("c", "", 1200, 900);

  // Genera el lienzo para hacer la grafica
  // mL y mH son los calores x_min y x_max
  // sL y sH son los calores y_min y y_max
  TH1F *frame = gPad->DrawFrame(mL, sL, mH, sH);

  // Ponemos titulo a los ejes
  frame->GetYaxis()->SetTitle("#sigma_{SI} (cm2)");
  frame->GetXaxis()->SetTitle("Wimp mass (GeV)");

  // Generamos un puntero a la imagen SI_mw_em1_e4_sdp_5em50_em37.JPG

  TImage *img = TImage::Open("plots/SI_mw_em1_e4_sdp_5em50_em37.JPG");

  // Poner los numeros en notacion cientifica
  frame->GetXaxis()->SetNoExponent(0);

  // Ejes logaritmicos
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  // Permite que se estire la imagen para llenar el area del pad
  img->SetConstRatio(kFALSE);
  // Dibuja la imagen en el pad
  img->Draw();

  // Genera un eje secundario en pb en vez de cm^2

  TGaxis *axis = new TGaxis(mH, sL, mH, sH, sL / pb2cm, sH / pb2cm, 50510, "+LG");
  axis->SetLabelOffset(0.01);         // Distancias de las etiquetas respecto al eje
  axis->SetLabelFont(42);             // Tipo de letra de las etiquetas
  axis->SetTitle("#sigma_{SI} (pb)"); // Titulo del eje
  axis->SetTitleFont(42);             // Tipo de letra del titulo
  axis->SetLabelSize(0.05);           // Tamano de las etiquetas
  axis->SetTitleSize(0.06);           // Tamano del titulo
  axis->SetTitleOffset(1.0);          // Distancia del titulo al eje
  axis->Draw();                       // Dibujamos el eje

  // Actualiza el lienzo
  gPad->Update();

  // Estilo de la linea principal, border (black, thicker)
  gex->SetLineColor(kBlack); // color negro
  gex->SetLineWidth(4);      // grosor de línea
  gex->SetLineStyle(9);      // estilo de línea (punteada, dash-dotted, etc.)
  gex->Draw("lsame");        // "l" = line, "same" = dibujar sobre el mismo pad

  // main line (yellow, thinner)
  TGraph *gex2 = (TGraph *)gex->Clone();
  gex2->SetLineColor(kRed); // color rojo
  gex2->SetLineWidth(2);    // grosor más delgado
  gex->SetLineStyle(9);     // mantiene estilo del primer gráfico
  gex2->Draw("L SAME");     // dibuja la línea continua encima

  std::string gname2 = "plots/";
  gname2 += gname;
  gname2 += ".png";
  c->SaveAs(gname2.c_str());
}
