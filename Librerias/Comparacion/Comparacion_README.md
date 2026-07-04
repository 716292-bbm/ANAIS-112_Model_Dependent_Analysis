# `Librerias/Comparacion/`

Notebooks y datos para **comparar** los ritmos y las curvas de exclusión obtenidos con los
distintos códigos de cálculo (código propio de ANAIS, RAPIDD, WIMPyDD y el `DMAnalysis`
interno), y para estudios auxiliares del análisis. Es la carpeta donde se realiza la
validación cruzada que respalda los resultados del trabajo.

## Notebooks

- **`Comparacion_Ritmo.ipynb`** — Compara el ritmo diferencial esperado calculado con cada
  código, para comprobar que el código propio reproduce los de referencia.
- **`Comparacion_entre_DMAnalysis.ipynb`** — Comparación específica con el código interno
  `DMAnalysis`.
- **`Comparacion_Plots_Con_Fondo_1.ipynb`**, **`Comparacion_Plots_Con_Fondo_2.ipynb`** —
  Curvas de exclusión incorporando el modelo de fondo, comparadas entre métodos.
- **`Comparacion_Plots_Exclusion_Sin_Fondo.ipynb`** — Curvas de exclusión con el método
  simplificado (sin modelo de fondo).
- **`Oprimizar_Intervalo_SI.ipynb`** — Estudio del intervalo de integración óptimo en energía
  para la interacción SI.
- **`GenerarLista_DMAnalysis.ipynb`** — Genera las listas de ritmo con el código `DMAnalysis`.

## Datos y resultados

- **`Archivos_ROOT/`** — Ritmos `rateDMAnalysis_SI_mw*.root` generados para cada masa de WIMP.
- **`Results/`** — Contiene `DMA_SI_TH1D.root`, el histograma de ritmo de `DMAnalysis` en el
  formato que consume el programa de exclusión.
- **`SI_varios*.root`** — Colecciones de curvas de exclusión SI bajo distintas condiciones
  (intervalos, resolución, correcciones).
- **`DAMA1.txt`, `DAMA2.txt`** — Contornos de la región compatible con la señal de DAMA/LIBRA,
  usados como referencia en las comparaciones.
- **`SI_ANAIS_*.txt`, `list_energies*.txt`, `datos_guardados_3.txt`** — Listas de puntos
  $(\sigma, m_\chi)$ y de energías empleadas en los cálculos y representaciones.
