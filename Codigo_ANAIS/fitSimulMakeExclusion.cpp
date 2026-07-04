// ===========================================================================
//  fitSimulMakeExclusion
// ---------------------------------------------------------------------------
//  Calcula la curva de exclusion (seccion eficaz limite frente a la masa del
//  WIMP) para ANAIS mediante un ajuste simultaneo (RooFit) de los 9 detectores.
//
//  Para cada masa de WIMP:
//    1. Se obtiene el espectro teorico (ritmo) de la senal de materia oscura.
//    2. Se le aplica la resolucion energetica del detector.
//    3. Se ajustan datos = nNorm*senal + nbkg*fondo simultaneamente en todos
//       los detectores, dejando libre nNorm (proporcional a la seccion eficaz).
//    4. Del valor y error de nNorm se deriva la seccion eficaz limite al CL
//       pedido.
//  Finalmente se dibuja y guarda la curva sigma(mW).
// ===========================================================================

// --- ROOT: histogramas, ficheros, graficos, ajustes, dibujo ---------------
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
#include <TMatrixDSym.h>

// --- C++ estandar ---------------------------------------------------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

// --- RooFit: ajuste estadistico (pdf, datasets, fit simultaneo) -----------
#include <RooFitResult.h>
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

// --- Librerias propias de ANAIS -------------------------------------------
#include <ADB.h>                 // Acceso a la base de datos (exposiciones, etc.)
#include <DMRatePlotHandle.hh>   // Calculo del ritmo teorico de materia oscura

using namespace RooFit;
using namespace std;

// Macro que activa la escritura/lectura del ritmo teorico en un fichero ROOT
// dentro de DMModelGetRate (ver mas abajo).
// #define ARCHIVO

// ===========================================================================
//  VARIABLES GLOBALES
// ===========================================================================

// --- Rango y binado del eje de energia (en keV, escala electron-equivalente) ---
// Se usan tanto para construir los histogramas de senal como el rango del fit.
double minEne = 0;
double maxEne = 100; // MARIA: extendido a ~30 keVee (100 keV NR para QFNa=0.3)
double binEne = 0.1; // Anchura de bin en keV

// --- Quenching Factor (QF) del Na y del I ---------------------------------
// El QF relaciona energia de retroceso nuclear (NR) con energia electron-
// equivalente (Eee): Eee = QF * ENR.
// Convencion: si el TGraph correspondiente no es nulo, se usa el QF dependiente
// de la energia (gQNa / gQI); si es nulo, se usa el valor constante (QNa / QI).
TGraph *gQNa = 0; // QF del sodio dependiente de energia (o nullptr)
TGraph *gQI = 0;  // QF del yodo dependiente de energia (o nullptr)
double QNa = 0;   // QF del sodio constante
double QI = 0;    // QF del yodo constante

// ===========================================================================
//  Conv2: convolucion del espectro con la resolucion gaussiana del detector
// ---------------------------------------------------------------------------
//  Convoluciona el ritmo de entrada (ritmo_sr) con una gaussiana de anchura
//  dependiente de la energia, sigma(E) = |p1 + ps*sqrt(E)|, y devuelve el
//  espectro resuelto muestreado en bins de anchura 'ebin' entre e1 y e2.
//
//  Parametros:
//    e1, e2   : rango de energia de salida (keV)
//    ebin     : anchura de bin de salida (keV)
//    p1, ps   : parametros de la resolucion sigma(E) = |p1 + ps*sqrt(E)|
//    S0       : array de salida con el ritmo convolucionado (lo rellena)
//    ritmo_sr : histograma de entrada (ritmo sin resolucion)
//
//  Devuelve 0 siempre (el resultado util se escribe en S0).
// ===========================================================================
int Conv2(double e1, double e2, double ebin, double p1, double ps,
          double *S0, TH1 *ritmo_sr)
{
  double PRECISION = 0.0001; // Precision para cortar el bucle en e2
  double minEne = 0.001;     // Energia minima permitida (evita E<=0 en sqrt)
  int ind;
  double ei, edif, en;                                 // variables de bucle
  double ef, estep, en1, en2, ens;                     // variables de intervalo
  double S01, SM1, PHI1, NORM1, S02, SM2, PHI2, NORM2; // acumuladores de integral
  double weight;
  int indArray;

  // Array auxiliar que precalcula el ritmo de entrada finamente muestreado,
  // para no llamar a GetBinContent dentro de los bucles internos.
  int arrayDim = 400;
  int nSig = 5; // Numero de +-sigma que abarca el intervalo de convolucion

  // Rango del array auxiliar: desde e1-nSig*sigma (o minEne) hasta e2+nSig*sigma
  // double auxsig = (ResolutionK<0 ? -nSig*ResolutionK : nSig*ResolutionK*sqrt(e1));
  double auxsig = nSig * fabs(p1 + ps * sqrt(e1));
  double arrayEi = (e1 - auxsig < minEne ? minEne : e1 - auxsig);
  double arrayEf = e2;
  arrayEf += nSig * fabs(p1 + ps * sqrt(e2));
  double arrayEbin = (arrayEf - arrayEi) / arrayDim; // anchura de bin del array aux.

  // Rellena el array auxiliar con el ritmo de entrada muestreado
  double *array_S0 = new double[arrayDim];

  ei = arrayEi;
  for (ind = 0; ind < arrayDim; ind++)
  {
    ei += arrayEbin;
    array_S0[ind] = ritmo_sr->GetBinContent(ritmo_sr->FindBin(ei));
  }

  // --- Bucle externo: recorre cada bin de energia de salida [ei, ef] ------
  ind = 0;
  for (ei = e1; ei < e2; ei += ebin)
  {
    // Corte de seguridad al llegar al final del rango
    if (fabs(e2 - ei) < PRECISION)
      break;

    ef = ei + ebin;
    estep = (ef - ei) / 20.; // 20 sub-pasos para integrar dentro del bin

    // Acumuladores del promedio dentro del bin de salida
    S01 = 0;
    NORM1 = 0;

    // --- Integral sobre el bin de salida (recorre edif de ei a ef) --------
    for (edif = ei; edif <= ef; edif += estep)
    {
      // Anchura gaussiana en esta energia y ventana de +-nSig*sigma
      auxsig = nSig * fabs(p1 + ps * sqrt(edif));
      en1 = (edif - auxsig < minEne ? minEne : edif - auxsig);
      en2 = edif + auxsig;
      ens = (en2 - en1) / 100.; // 100 sub-pasos para la convolucion gaussiana

      // Acumuladores de la convolucion en este punto edif
      S02 = 0;
      NORM2 = 0;

      // Caso resolucion nula: se toma directamente el valor del array
      if (auxsig == 0)
      {
        indArray = (int)((edif - arrayEi) / arrayEbin);
        S02 = array_S0[indArray];
        NORM2 = 1;
        en2 = 0; // evita entrar en el bucle gaussiano de abajo
      }

      // --- Convolucion gaussiana: integra ritmo(en) * gauss(edif-en) ------
      for (en = en1; en <= en2; en += ens)
      {
        // Peso gaussiano (sin normalizar; se normaliza con NORM2 al final)
        auxsig = fabs(p1 + ps * sqrt(en));
        weight = exp(-(edif - en) * (edif - en) / 2. / auxsig / auxsig) / auxsig;

        // Indice del array auxiliar correspondiente a la energia 'en'
        indArray = (int)((en - arrayEi) / arrayEbin);

        // Acumula numerador y normalizacion
        // MARIA: mas adelante se podria hacer integracion trapezoidal...
        S02 += weight * array_S0[indArray];
        NORM2 += weight;

      } // fin bucle en (convolucion)

      // Anade el valor convolucionado en edif al promedio del bin de salida
      S01 += S02 / NORM2;
      NORM1++;

    } // fin bucle edif (integral en el bin)

    // Valor final del bin de salida (promedio sobre el bin)
    S0[ind] = S01 / NORM1;
    ind++;

  } // fin bucle ei (bins de salida)

  delete[] array_S0;

  return 0;
}

// ===========================================================================
//  DMModelGetRate: calcula el ritmo teorico de senal WIMP con DMRate
// ---------------------------------------------------------------------------
//  Genera el espectro esperado de materia oscura (eventos/kg/dia/keV) para una
//  masa mw y seccion eficaz sigma dadas, segun el modelo de spin elegido.
//  Devuelve un histograma con el ritmo (en energia de retroceso nuclear).
//
//  Parametros:
//    mw        : masa del WIMP (GeV)
//    sigma     : seccion eficaz (aqui suele pasarse 1, sirve de normalizacion)
//    qfModel   : modelo de quenching factor (no usado directamente aqui)
//    SpinModel : 0-SI (independiente de spin), 1-SDp (proton), 2-SDn (neutron)
// ===========================================================================
TH1F *DMModelGetRate(double mw, double sigma, int qfModel, int SpinModel)
{
#ifdef ARCHIVO
  // Fichero donde se guardara el histograma calculado (modo ARCHIVO)
  TFile *f = new TFile("rate_DMAnalysis.root", "UPDATE");
  // f->mkdir("SI_rates")->cd();
  TString file_name = "rate";
#endif

  DMRate *rate = new DMRate();                    // Objeto que calcula el ritmo
  int err = rate->Initialize((char *)"rate.dat"); // Inicializa desde rate.dat
  if (err != 0)
    return 0;         // Si falla la inicializacion, aborta
  rate->SetMW(mw);    // Fija la masa del WIMP

  if (SpinModel == 0) // Spin-Independent
  {
    rate->SetSigSI(sigma); // Seccion eficaz SI = sigma
    rate->SetSigSD(0);     // Seccion eficaz SD = 0
    rate->SetTheta(0);     // theta = 0
  }

  if (SpinModel == 1) // Spin-Dependent Proton
  {
    rate->SetSigSI(0);     // SI = 0
    rate->SetSigSD(sigma); // SD = sigma
    for (int iel = 0; iel < rate->GetNElements(); iel++)
    {
      rate->GetElement(iel)->SetTheta(0); // theta = 0 (acoplo a proton)
    }
  }

  if (SpinModel == 2) // Spin-Dependent Neutron
  {
    rate->SetSigSI(0);     // SI = 0
    rate->SetSigSD(sigma); // SD = sigma
    for (int iel = 0; iel < rate->GetNElements(); iel++)
    {
      rate->GetElement(iel)->SetTheta(1.57079632679); // theta = pi/2 (neutron)
    }
  }

  // Asigna el quenching factor del Na al primer elemento (grafico o constante)
  if (gQNa)
    rate->GetElement(0)->SetREF(gQNa);
  else
    rate->GetElement(0)->SetREF(QNa);

  if (gQI)
    rate->GetElement(1)->SetREF(gQI);
  else
    rate->GetElement(1)->SetREF(QI);


  DMRatePlotHandle plothl(rate);

  TH1F *hrate = plothl.GetRate(1000); // Calcula el histograma con 1000 bins

#ifdef ARCHIVO
  // Nombre y titulo del histograma, y lo escribe en el fichero
  TString histName = file_name + Form("_mw%0.1f", mw);
  TString histTitle = Form("DMA Rate - M_{W} = %0.1f GeV;Energy (keV);Rate (events/kg/day/keV)", mw);

  hrate->SetName(histName);
  hrate->SetTitle(histTitle);

  hrate->Write();         // Guarda el histograma en el fichero
  hrate->SetDirectory(0); // Lo desasocia del fichero (para poder cerrarlo)

  f->cd();
  f->Close();
#endif
  return hrate;
}

// ===========================================================================
//  DMModelGetRateEeeFromFile: lee el ritmo teorico desde un fichero externo
// ---------------------------------------------------------------------------
//  Para modelos generados fuera de este programa (WIMPYDD, RAPIDD, Python,
//  Migdal...). Lee los histogramas de ritmo en energia de retroceso nuclear
//  del Na y del I, los pasa a energia electron-equivalente aplicando el QF, y
//  los combina segun la abundancia estequiometrica del NaI (23 Na + 127 I).
//
//  Parametros:
//    fileName : fichero ROOT con los histogramas
//    nameNa   : nombre del histograma del sodio dentro del fichero
//    nameI    : nombre del histograma del yodo dentro del fichero
//  Devuelve el histograma combinado en energia electron-equivalente (o nullptr).
// ===========================================================================
TH1D *DMModelGetRateEeeFromFile(std::string fileName, std::string nameNa, std::string nameI)
{
  // Abre el fichero en modo lectura y verifica que sea valido
  TFile *f = TFile::Open(fileName.c_str(), "READ");
  if (!f || f->IsZombie() || !f->IsOpen())
  {
    std::cerr << "[DMModelGetRateEeeFromFile] ERROR: no se pudo abrir '"
              << fileName.c_str() << "'\n";
    return nullptr;
  }

  // Lee el histograma del sodio y lo desasocia del fichero
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
  // Lee el histograma del yodo y lo desasocia del fichero
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

  // Histograma de salida en energia electron-equivalente (Eee)
  int nBins = (maxEne - minEne) / binEne;
  TH1D *hWimp = new TH1D(Form("%s_%s", nameNa.c_str(), nameI.c_str()), "", nBins, minEne, maxEne);
  hWimp->SetDirectory(0);

  // Para cada bin en Eee, convierte a energia NR con el QF y suma Na + I
  for (int ii = 1; ii <= nBins; ii++)
  {
    double Eee = hWimp->GetBinCenter(ii);

    // --- Sodio: Eee -> ENR con el QF del Na; el ritmo se escala por 1/QF ---
    double qNa = QNa;
    if (gQNa)
      qNa = gQNa->Eval(Eee); // QF dependiente de energia si hay grafico
    double ENR_Na = Eee / qNa;
    double rateNa = hNa->GetBinContent(hNa->FindBin(ENR_Na)) / qNa;

    // --- Yodo: idem con el QF del I ----------------------------------------
    double qI = QI;
    if (gQI)
      qI = gQI->Eval(Eee);
    double ENR_I = Eee / qI;
    double rateI = hI->GetBinContent(hI->FindBin(ENR_I)) / qI;

    // Combinacion por abundancia molar del NaI: 23 (Na) + 127 (I) = 150
    hWimp->SetBinContent(ii, (23 * rateNa + 127 * rateI) / 150.);
  }

  f->Close();
  delete f;

  return hWimp;
}

// ===========================================================================
//  main: lee configuracion por linea de comandos y calcula la exclusion
// ===========================================================================
int main(int argc, char **argv)
{
  // --- Comprobacion del numero de argumentos --------------------------------
  // Todos los modos se controlan por argumentos: NO hay que recompilar para
  // cambiar de modelo, spin, QF, ANOD o tratamiento de la resolucion.
  if (argc < 9)
  {
    std::cout << "Usage: fitSimulMakeExclusion cl eneIni eneEnd thmodel spinModel qfModel ANOD resolution_p" << std::endl;
    std::cout << "  cl           : Confidence Level (90, 95)" << std::endl;
    std::cout << "  eneIni       : Energia minima del fit (keV)" << std::endl;
    std::cout << "  eneEnd       : Energia maxima del fit (keV)" << std::endl;
    std::cout << "  thmodel      : 0-DMAnalysis, 1-Python, 2-RAPIDD, 3-WIMPYDD, 4-DMAnalysis(archivo), 5-Migdal" << std::endl;
    std::cout << "  spinModel    : 0-SI, 1-SD-proton, 2-SD-neutron" << std::endl;
    std::cout << "  qfModel      : 1-DAMA, 2-ANAIS CTE, 3-TAMARA" << std::endl;
    std::cout << "  ANOD         : incluir poblacion ALE (0-no, 1-si)" << std::endl;
    std::cout << "  resolution_p : 0-sin resolucion, 1-gausiana, 2-convolucion" << std::endl;
    return 1;
  }

  // --- Lectura de argumentos ------------------------------------------------
  double cl = atof(argv[1]);         // Confidence Level (90 o 95)
  double min = atof(argv[2]);        // Energia minima del fit (keV)
  double max = atof(argv[3]);        // Energia maxima del fit (keV)
  int thmodel = atoi(argv[4]);       // Modelo teorico (0-5, ver Usage)
  int SpinModel = atoi(argv[5]);     // 0-SI, 1-SD-proton, 2-SD-neutron
  int qfModel = atoi(argv[6]);       // 1-DAMA, 2-ANAIS CTE, 3-TAMARA
  bool ANOD = atoi(argv[7]);         // Incluir poblacion ALE (0/1)
  int resolution_p = atoi(argv[8]);  // 0-sin resolucion, 1-gausiana, 2-convolucion

  // --- Parametros de la resolucion energetica -------------------------------
  // Solo se cargan si se usa la convolucion (resolution_p == 2).
  // sigma(E) = param_1 + param_2*sqrt(E), leidos de un ajuste previo.
  double param_1 = 0, param_2 = 0;

  if (resolution_p == 2)
  {
    const char *lowResFile = "fitsResolution.root";
    TFile *fres = TFile::Open(lowResFile, "READ");
    TF1 *ffres = (TF1 *)fres->Get("fresD0");

    param_1 = ffres->GetParameter(0);
    param_2 = ffres->GetParameter(1);
  }

  // --- Resumen de la configuracion por pantalla -----------------------------
  std::cout << " cl           " << cl << std::endl;
  std::cout << " min          " << min << std::endl;
  std::cout << " max          " << max << std::endl;
  std::cout << " thmodel      " << thmodel << std::endl;
  std::cout << " SpinModel    " << SpinModel << " (0-SI, 1-SDp, 2-SDn)" << std::endl;
  std::cout << " qfModel      " << qfModel << " (1-DAMA, 2-ANAIS CTE, 3-TAMARA)" << std::endl;
  std::cout << " ANOD         " << ANOD << std::endl;
  std::cout << " resolution_p " << resolution_p << std::endl;

  if (SpinModel < 0 || SpinModel > 2)
  {
    std::cout << " Spin Model not Valid!!" << std::endl;
    exit(0);
  }

  // --- Configuracion del Quenching Factor segun qfModel ---------------------
  // MARIA 100326. Set QF MODE
  if (qfModel == 1) // DAMA: QF constantes de DAMA
  {
    QNa = 0.3;
    QI = 0.09;
  }
  else if (qfModel == 2) // ANAIS CTE: aqui puestos a 1 (sin quenching)
  {
    // QNa = 0.2;
    // QI = 0.06;
    QNa = 1.0;
    QI = 1.0;
  }
  else // TAMARA: QF dependiente de energia leido de graficos en fichero
  {
    TFile *f = new TFile("QFTamara.root", "READ");
    if (!f)
    {
      std::cout << " CANNOT READ QUENCHING FACTOR FILE!!" << std::endl;
      exit(0);
    }
    gQNa = (TGraph *)f->Get("gNa"); // QF(E) del sodio
    if (!gQNa)
    {
      std::cout << " CANNOT FIND gNa IN QUENCHING FACTOR FILE!!" << std::endl;
      exit(0);
    }
    gQI = (TGraph *)f->Get("gI"); // QF(E) del yodo
    if (!gQI)
    {
      std::cout << " CANNOT FIND gI IN QUENCHING FACTOR FILE!!" << std::endl;
      exit(0);
    }
  }

  // --- Variable de energia de RooFit (eje del ajuste) -----------------------
  RooRealVar energy("energy", "energy", minEne, maxEne, "keV");

  // --- Seleccion de anios y detectores --------------------------------------
  // CHANGE HERE FOR DIFFERENT YEARS - Seleccionamos la exposicion
  std::vector<int> years = {1, 2, 3, 4, 5, 6};
  // std::vector<int> years = {1,2,3};
  int nyears = years.size(); // Numero de anios
  std::vector<int> detectors = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  // std::vector<int> detectors = {2,4,6,7,8}; // usar los 5 mas limpios
  int ndet = detectors.size(); // Numero de detectores

  // Categoria de RooFit: un "tipo" por detector, necesario para el fit simultaneo
  RooCategory detCat("det", "det");
  for (int det = 0; det < ndet; det++)
    detCat.defineType(Form("det%d", detectors[det]));
  detCat.Print("V");

  ///////////////////////////////////////////////////////////////////////////
  // LECTURA DE DATOS
  ///////////////////////////////////////////////////////////////////////////

  // --- Exposicion (kg x dia) de cada detector -------------------------------
  ADBTime dt;
  // El livetime es distinto para cada detector (valores de 6 anios).

  std::vector<double> exposure; // Vector con la exposicion de cada detector

  // Live time D0..D8 (dias):
  //   D0:2031.38 D1:2033.20 D2:2029.52 D3:2022.55 D4:2033.01
  //   D5:2030.18 D6:2032.27 D7:2031.02 D8:2020.29

  // Con scaleFactor se puede escalar la exposicion (p.ej. 8./6. para 8 anios)
  double scaleFactor = 1;

  // Exposicion = livetime(dias) * masa(12.5 kg) * scaleFactor
  exposure.push_back(2031.38 * 12.5 * scaleFactor); // D0
  exposure.push_back(2033.20 * 12.5 * scaleFactor); // D1
  exposure.push_back(2029.52 * 12.5 * scaleFactor); // D2
  exposure.push_back(2022.55 * 12.5 * scaleFactor); // D3
  exposure.push_back(2033.01 * 12.5 * scaleFactor); // D4
  exposure.push_back(2030.18 * 12.5 * scaleFactor); // D5
  exposure.push_back(2032.27 * 12.5 * scaleFactor); // D6
  exposure.push_back(2031.02 * 12.5 * scaleFactor); // D7
  exposure.push_back(2020.29 * 12.5 * scaleFactor); // D8
  // for(int i=0;i<nyears;i++) exposure.push_back(dt.GetExposure(-1,years[i],1)/ndet);

  // --- Lectura de los datos experimentales (fondo medido) -------------------
  // Construye el nombre del fichero de datos a partir de los anios elegidos.
  std::string name = "../../data/BEhistos_year";
  for (int i = 0; i < nyears; i++)
    name += std::to_string(years[i]);
  name += ".root";

  TFile *file = new TFile(name.c_str(), "read"); // Fichero con los histogramas de datos

  // Mapa {detector -> histograma de datos} para el fit simultaneo
  std::map<std::string, TH1 *> hdataMap;
  for (int det = 0; det < ndet; det++)
  {
    // Nombre del histograma: hbea_<anios>y_D<detector>
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
    // Pasa de ritmo (cuentas/kg/dia/keV) a numero de cuentas para el fit extendido
    h->Scale(h->GetBinWidth(1) * exposure[det]);
    hdataMap[Form("det%d", detectors[det])] = h;
  }

  // Combina el mapa de histogramas en un unico dataset de RooFit por categoria
  RooDataHist *dataData = new RooDataHist("data", "data", RooArgSet(energy), detCat, hdataMap);

  ///////////////////////////////////////////////////////////////////////////
  // MODELO DE FONDO
  ///////////////////////////////////////////////////////////////////////////
  // single hits: fichero con el modelo de fondo (con o sin ALE segun ANOD)
  std::string bkgname = "../../backgroundModel/backgroundModel_single_y123456";
  // for (int i=0; i<nyears;i++) bkgname+=std::to_string(years[i]);
  bkgname += (ANOD ? "_conANOD.root" : ".root"); // incluye ALE si ANOD=1
  TFile *fbkg = new TFile(bkgname.c_str(), "read");

  // --- Factor de normalizacion de la senal (proporcional a la seccion eficaz) ---
  // Es comun a todos los detectores y se deja LIBRE en el fit. Rango amplio y
  // que admite valores negativos para no sesgar el ajuste.
  RooRealVar *nNorm = new RooRealVar("nNorm", "", 1, -1e10, 1e10);
  nNorm->setRange(-1e6, 1e6);
  nNorm->setVal(1e4);
  nNorm->setConstant(kFALSE); // parametro libre del ajuste

  ///////////////////////////////////////////////////////////////////////////
  // ARRAY DE MASAS DE WIMP
  ///////////////////////////////////////////////////////////////////////////
  // El array de masas depende de thmodel (5 = Migdal) en tiempo de ejecucion,
  // no de una macro de compilacion.
  std::vector<double> mw;

  if (thmodel == 5) // Migdal: lista fija de masas (CARMEN)
  {
    mw = {
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
  }
  else
  {
    // Resto de modelos: rejilla logaritmica generada, con 4 puntos por decada
    // (imw + 0.2*j) en las decadas de las unidades, decenas, centenas y miles.
    for (int imw = 2; imw <= 9; ++imw)
      for (int j = 0; j < 4; ++j)
        mw.push_back(imw + 0.2 * j);

    for (int imw = 1; imw <= 9; ++imw)
      for (int j = 0; j < 4; ++j)
        mw.push_back((imw + 0.2 * j) * 10);

    for (int imw = 1; imw <= 9; ++imw)
      for (int j = 0; j < 4; ++j)
        mw.push_back((imw + 0.2 * j) * 100);

    for (int imw = 1; imw <= 10; ++imw)
      for (int j = 0; j < 4; ++j)
        mw.push_back((imw + 0.2 * j) * 1000);
  }

  ///////////////////////////////////////////////////////////////////////////
  // NOMBRES DE LOS HISTOGRAMAS DE SENAL (uno por masa, para Na y para I)
  ///////////////////////////////////////////////////////////////////////////
  std::vector<std::string> nombres_Na;
  std::vector<std::string> nombres_I;
  nombres_Na.reserve(mw.size());
  nombres_I.reserve(mw.size());

  if (thmodel == 5) // Migdal: nombres con las cadenas exactas del fichero
  {
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
  }
  else
  {
    // Resto de modelos: formatea cada masa con 1 decimal (ej. "12.0")
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
  }

  ///////////////////////////////////////////////////////////////////////////
  // FICHERO DE SENAL segun el modelo teorico elegido
  ///////////////////////////////////////////////////////////////////////////
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
  if (thmodel == 3) // WIMPYDD: el fichero depende del modelo de spin
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

  ///////////////////////////////////////////////////////////////////////////
  // BUCLE PRINCIPAL SOBRE LAS MASAS DE WIMP
  ///////////////////////////////////////////////////////////////////////////
  std::vector<double> sigma; // Seccion eficaz limite para cada masa
  for (int imw = 0; imw < (int)mw.size(); imw++)
  {
    // pdf simultaneo (uno por detector) para esta masa
    RooSimultaneous *simWimpHyp = new RooSimultaneous("simWimpHyp", "simultaneuos pdf WimpHyp", detCat);

    // --- Espectro teorico de senal (sin resolucion) para esta masa ---------
    // thmodel==0 lo calcula con DMRate; el resto lo lee de fichero.
    TH1 *hWimpNoRes;
    if (thmodel == 0)
    {
      hWimpNoRes = DMModelGetRate(mw[imw], 1, qfModel, SpinModel);
    }
    else
    {
      hWimpNoRes = DMModelGetRateEeeFromFile(fileName, nombres_Na[imw], nombres_I[imw]);
    }

    // --- Aplicacion de la resolucion energetica ----------------------------
    int nBins = (maxEne - minEne) / binEne;
    const int nBins_2 = 1000;
    double rate_res[nBins_2]; // buffer de salida para la convolucion (Conv2)
    double integral_hWimpNoRes = hWimpNoRes->Integral(1, nBins);

    // hWimp sera el espectro final (con resolucion aplicada)
    TH1 *hWimp = (TH1 *)hWimpNoRes->Clone(Form("%s_res", hWimpNoRes->GetName()));

    if (resolution_p == 1)
    {
      // Metodo 1: Monte Carlo. Se generan eventos segun el espectro y se
      // desplazan con una gaussiana de anchura sigma(E), reconstruyendo el
      // espectro resuelto por muestreo.
      hWimp->Reset();
      TRandom ran;
      int nEvRes = 10000000; // numero de eventos generados

      // CARMEN 26/03/2026: carga la funcion de resolucion sigma(E)
      const char *lowResFile = "/media/storage2/tamara/resMartaByDet/fitsResolution.root";
      TFile *fres = TFile::Open(lowResFile, "READ");
      TF1 *ffres = (TF1 *)fres->Get("fresD0");

      for (int ii = 0; ii < nEvRes; ii++)
      {
        double ee = hWimpNoRes->GetRandom();
        double sigma = ffres->GetParameter(0) + ffres->GetParameter(1) * sqrt(ee);
        double eeRes = ran.Gaus(ee, sigma); // aplica la resolucion
        hWimp->Fill(eeRes);
      }

      // Renormaliza para conservar la integral del espectro original
      double integral_hWimp = hWimp->Integral(1, nBins);
      double factor_norm = integral_hWimpNoRes / integral_hWimp;
      hWimp->Scale(factor_norm);
    }

    if (resolution_p == 2)
    {
      // Metodo 2: convolucion analitica con Conv2 (mas rapido y sin ruido MC)
      hWimp->Reset();

      int asd = Conv2(0, 100, 0.1, param_1, param_2, rate_res, hWimpNoRes);

      // Vuelca el resultado de la convolucion al histograma
      for (int ii = 1; ii <= nBins; ii++)
      {
        hWimp->SetBinContent(ii, rate_res[ii - 1]);
      }
    }
    // (resolution_p == 0: no se aplica resolucion, hWimp queda igual al original)

    // --- Convierte el espectro de senal en pdf de RooFit -------------------
    RooDataHist *rdhWimp = new RooDataHist("rateWimp", "", RooArgSet(energy), hWimp);
    RooHistPdf *pdf_Wimp = new RooHistPdf("pdf_Wimp", "", energy, energy, *rdhWimp, 1);

    // Construct the model for every detector
    // --- Construye el modelo de cada detector: nNorm*senal + nbkg*fondo -----
    for (int det = 0; det < ndet; det++)
    {
      // Lee el histograma del fondo de este detector
      TH1F *hbkg = (TH1F *)fbkg->Get(Form("hD%d", detectors[det]));

      hbkg->Smooth(); // suaviza para quitar fluctuaciones (opcional)

      // pdf del fondo a partir del histograma
      RooDataHist *rdhBkg = new RooDataHist(Form("rdhBkg_D%d", detectors[det]), "", RooArgSet(energy), hbkg);
      RooHistPdf *pdf_bkg = new RooHistPdf(Form("pdf_bkg_D%d", detectors[det]), "", energy, energy, *rdhBkg, 1);

      // Factor de escala a numero de cuentas: exposicion * anchura de bin
      double fscale = exposure[detectors[det]] * hbkg->GetBinWidth(1); // kg*dia*binWidth

      // Numero de cuentas de fondo integrado en el rango del fit
      double intebkg = hbkg->Integral(hbkg->FindBin(min), hbkg->FindBin(max) - 1) * fscale;

      // Numero total de cuentas de fondo: FIJO en el ajuste (no se toca)
      RooRealVar *nbkg = new RooRealVar(Form("nbkg_D%d", detectors[det]), "", 1000, intebkg / 2, intebkg * 2);
      nbkg->setRange(intebkg / 2, intebkg * 2);
      nbkg->setVal(intebkg);
      nbkg->setConstant(kTRUE); // fondo fijado; solo nNorm (senal) es libre

      // Modelo del detector = nNorm*pdf_Wimp + nbkg*pdf_bkg
      RooAddPdf *model_det = new RooAddPdf(Form("model_det%d", detectors[det]), "WIMP + background",
                                           RooArgList(*pdf_Wimp, *pdf_bkg),
                                           RooArgList(*nNorm, *nbkg));
      simWimpHyp->addPdf(*model_det, Form("det%d", detectors[det]));
    }

    ///////////////////////////////////////
    // AJUSTE (FIT)
    ///////////////////////////////////////
    double sig; // seccion eficaz limite para esta masa
    if (min > 0)
    {
      // --- Caso normal: rango de fit fijo [min, max] -----------------------
      nNorm->setVal(1e4);                    // valor inicial del parametro libre
      energy.setRange("fitRange", min, max); // rango de energia del ajuste

      // Ajuste extendido simultaneo de todos los detectores
      simWimpHyp->fitTo(*dataData, SumCoefRange("fitRange"), RooFit::Range("fitRange"), SumW2Error(false), Save(true), Extended());

      // Para imprimir el resultado del fit, descomentar:
      // RooFitResult *results = simWimpHyp->fitTo(*dataData, SumCoefRange("fitRange"), RooFit::Range("fitRange"), SumW2Error(false), Save(true), Extended());
      // results->Print();

      // Valor y error del numero de cuentas de senal ajustado
      double nNorm_val = nNorm->getVal();
      double nNorm_err = nNorm->getError();

      // Limite superior al CL pedido (1.64 sigma para 95%, 1.28 para 90%)
      double valcl = (cl == 95 ? nNorm_val + nNorm_err * 1.64 : nNorm_val + nNorm_err * 1.28);
      std::cout << " counts at 90\% CL " << valcl << std::endl;

      // sigma limite = cuentas limite / integral del ritmo teorico en el rango
      sig = valcl / hWimp->Integral(hWimp->FindBin(min), hWimp->FindBin(max), "width");
    }
    else
    {
      // --- Caso "optimo": explora rangos y toma el mejor limite ------------
      // MARIA: ESTO NO FUNCIONA (pendiente de revisar)
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
      sig = 1e300;
      double minInt = 4; // anchura minima del intervalo
      double step = 1;   // paso de exploracion
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
          // Se queda con el rango que da el limite mas restrictivo (menor sig)
          if (valcl > 0 && valcl / integral < sig)
            sig = valcl / integral;
        }
      }
    }

    // --- Normalizacion final de la seccion eficaz --------------------------
    // Exposicion total (kg*dia) sumando todos los detectores
    double totalExposure = 0;
    for (int det = 0; det < ndet; det++)
      totalExposure += exposure[det];
    std::cout << " total exposure " << totalExposure << " kgxday" << std::endl;

    // Convierte a cm^2 (1 pb = 1e-36 cm^2) y normaliza por la exposicion total
    sig *= 1e-36 / totalExposure;

    // Guarda la seccion eficaz de esta masa
    sigma.push_back(sig);

    // Libera memoria de los objetos de esta iteracion
    delete simWimpHyp;
    delete rdhWimp;
    delete pdf_Wimp;
    delete hWimp;
  }
  // fin del bucle en masas

  ///////////////////////////////////////////////////////////////////////////
  // RESULTADOS Y GRAFICA
  ///////////////////////////////////////////////////////////////////////////

  // Imprime la tabla final mW -> sigma
  for (int imw = 0; imw < (int)mw.size(); imw++)
    std::cout << " mw: " << mw[imw] << " sigma: " << sigma[imw] << " cm2 " << std::endl;

  // Grafico XY con la curva de exclusion: sigma frente a mW
  TGraph *gex = new TGraph(mw.size(), mw.data(), sigma.data());

  // Fichero ROOT donde se guarda el grafico
  //  TFile * f = new TFile("plots/results.root","recreate");
  TFile *f = new TFile("plots/SI_varios_2.root", "update");

  // --- Construye el nombre del objeto codificando la configuracion usada ----
  std::string gname = "gA112_6y_";
  // std::string gname = "gCOSINE_3y_";
  gname += (cl == 90 ? "90_" : "95_"); // CL
  if (min > 0)                          // rango de energia (o "opt")
  {
    gname += std::to_string(min);
    gname += "_";
    gname += std::to_string(max);
  }
  else
    gname += "opt";
  gname += "_QF"; // modelo de quenching factor
  if (qfModel == 1)
    gname += "dama";
  else if (qfModel == 2)
    gname += "cte";
  else
    gname += "tamara";
  if (ANOD)
    gname += "_ANOD"; // sufijo si se incluyo la poblacion ALE

  // Sufijo segun el modelo teorico
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

  if (thmodel == 5)
  {
    gname += "_MIGDALL";
  }

  if (resolution_p == 0)
  {
    gname += "_rp0";
  }

   if (resolution_p == 1)
  {
    gname += "_rp1";
  }

     if (resolution_p == 2)
  {
    gname += "_rp2";
  }

  // --- Limites de los ejes del plot -----------------------------------------
  double pb2cm = 1e-36;  // factor de conversion pb -> cm^2
  double mL = 1;         // masa minima (eje X)
  double mH = 1e4;       // masa maxima (eje X)
  double sL = 5e-50;     // sigma minima (eje Y) para SI
  double sH = 1e-37;     // sigma maxima (eje Y) para SI

  // Sufijo por modelo de spin y, para SD, ajuste del rango del eje Y
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


  
  if (thmodel == 5)
  {
    mL = 0.01; 
    sL = 5e-45;
    sH = 1e-30;
  }

  // Guarda la curva en el fichero ROOT con el nombre construido
  gex->Write(gname.c_str());
  f->Close();

  // --- Dibujo de la curva de exclusion --------------------------------------
  TCanvas *c = new TCanvas("c", "", 1200, 900);

  // Marco del plot: (mL,sL) esquina inferior-izquierda, (mH,sH) superior-derecha
  TH1F *frame = gPad->DrawFrame(mL, sL, mH, sH);

  // Titulos de los ejes
  frame->GetYaxis()->SetTitle("#sigma_{SI} (cm2)");
  frame->GetXaxis()->SetTitle("Wimp mass (GeV)");

  // Imagen de fondo (curvas de referencia de otros experimentos)
  TImage *img = TImage::Open("plots/SI_mw_em1_e4_sdp_5em50_em37.JPG");

  // Notacion cientifica en el eje X
  frame->GetXaxis()->SetNoExponent(0);

  // Ejes logaritmicos en X e Y
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  // Permite estirar la imagen para rellenar el pad y la dibuja de fondo
  img->SetConstRatio(kFALSE);
  img->Draw();

  // --- Eje secundario en pb (a la derecha) ----------------------------------
  TGaxis *axis = new TGaxis(mH, sL, mH, sH, sL / pb2cm, sH / pb2cm, 50510, "+LG");
  axis->SetLabelOffset(0.01);         // separacion etiquetas-eje
  axis->SetLabelFont(42);             // fuente de las etiquetas
  axis->SetTitle("#sigma_{SI} (pb)"); // titulo del eje
  axis->SetTitleFont(42);             // fuente del titulo
  axis->SetLabelSize(0.05);           // tamano de las etiquetas
  axis->SetTitleSize(0.06);           // tamano del titulo
  axis->SetTitleOffset(1.0);          // separacion titulo-eje
  axis->Draw();

  gPad->Update();

  // --- Dibujo de la curva (doble trazo: contorno negro + linea roja) --------
  // Contorno (negro, grueso)
  gex->SetLineColor(kBlack);
  gex->SetLineWidth(4);
  gex->SetLineStyle(9);
  gex->Draw("lsame"); // "l" = linea, "same" = sobre el mismo pad

  // Linea principal (roja, mas fina) por encima
  TGraph *gex2 = (TGraph *)gex->Clone();
  gex2->SetLineColor(kRed);
  gex2->SetLineWidth(2);
  gex->SetLineStyle(9);
  gex2->Draw("L SAME");

  // --- Guarda la figura en PNG con el mismo nombre que el grafico -----------
  std::string gname2 = "plots/";
  gname2 += gname;
  gname2 += ".png";
  c->SaveAs(gname2.c_str());
}