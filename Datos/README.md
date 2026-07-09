# `Datos/`

Datos experimentales de entrada del análisis.

## Contenido

- **`BEhistos_year123456.root`** — Histogramas de bajo fondo (*Bulk Events*) medidos por
  ANAIS-112, correspondientes a los **seis primeros años** de toma de datos y promediados
  sobre los nueve detectores. Contiene el espectro en energía en la región de baja energía
  (region of interest), que es la entrada experimental frente a la que se ajusta el modelo
  de fondo más señal para construir las curvas de exclusión.
- **`backgroundModel_single_y123456.root`** — Modelo de fondo de ANAIS-112 obtenido por
  simulación Monte Carlo (single hits), con un histograma por detector (`hD0`–`hD8`). Es la
  componente de fondo del ajuste `datos = señal + fondo`.
- **`backgroundModel_single_y123456_conANOD.root`** — Misma versión del modelo de fondo que
  incluye además la población de eventos anómalos (*Anomalous Light Events*, ALE). El programa
  de exclusión lo usa en lugar del anterior cuando se activa la opción `ANOD=1`.

## Uso

Este archivo lo leen los notebooks del análisis (por ejemplo `TFG_Calculo_N_Exp.ipynb`) y
el programa `Codigo_ANAIS/fitSimulMakeExclusion.cpp`, que accede a los histogramas por
detector con el patrón de nombre `hbea_<años>y_D<n>`.

> Formato ROOT: puede abrirse desde Python con `uproot` (sin necesidad de tener ROOT
> instalado) o directamente con ROOT.