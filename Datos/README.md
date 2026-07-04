# `Datos/`

Datos experimentales de entrada del análisis.

## Contenido

- **`BEhistos_year123456.root`** — Histogramas de bajo fondo (*Bulk Events*) medidos por
  ANAIS-112, correspondientes a los **seis primeros años** de toma de datos y promediados
  sobre los nueve detectores. Contiene el espectro en energía en la región de baja energía
  (region of interest), que es la entrada experimental frente a la que se ajusta el modelo
  de fondo más señal para construir las curvas de exclusión.

## Uso

Este archivo lo leen los notebooks del análisis (por ejemplo `TFG_Calculo_N_Exp.ipynb`) y
el programa `Codigo_ANAIS/fitSimulMakeExclusion.cpp`, que accede a los histogramas por
detector con el patrón de nombre `hbea_<años>y_D<n>`.

> Formato ROOT: puede abrirse desde Python con `uproot` (sin necesidad de tener ROOT
> instalado) o directamente con ROOT.
