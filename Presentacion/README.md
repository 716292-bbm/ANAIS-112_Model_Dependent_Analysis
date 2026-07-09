# `Presentacion/`

Notebooks, datos y figuras organizados según las **etapas sucesivas del análisis** descritas
en la memoria del TFG. Cada subcarpeta corresponde a una mejora del método de construcción de
las curvas de exclusión, desde el planteamiento más simple hasta el resultado final. La mayor
parte de las figuras se guardan en formato `.svg`.

## Subcarpetas

### `Preliminar/`
Método simplificado: se atribuye la **totalidad** de los eventos registrados a la interacción
con WIMPs, sin modelo de fondo. Es el punto de partida del análisis (`SI_preliminar.ipynb`).

### `Fondo/`
Incorporación del **modelo de fondo** de ANAIS-112 obtenido por simulación Monte Carlo. El
espectro se ajusta como fondo más señal, lo que mejora la sensibilidad en varios órdenes de
magnitud (`SI_fondo.ipynb`).

### `ModeloFondo/`
Representación del **espectro de baja energía** de ANAIS-112 comparando los datos medidos, el
modelo de fondo obtenido por simulación Monte Carlo y el fondo más la población de ALEs
(`ModeloFondo.ipynb`, figura `Modelosdefondo.pdf`).

### `Intervalo/`
Optimización del **intervalo de integración** en energía: para cada masa de WIMP se busca el
subintervalo que proporciona el límite más restrictivo (`SI_Intervalo.ipynb`).

### `ALE/`
Inclusión de la población de **eventos anómalos** (*Anomalous Light Events*), que mejora el
modelo de fondo en la región de 1–2 keVee (`SI_ALE.ipynb`).

### `FactorQ/`
Estudio del impacto del **factor de *quenching*** sobre las curvas de exclusión (`SI_QF.ipynb`).

### `QF_DAMA/` y `QF_Tamara/`
Resultados completos (SI, SDp, SDn) empleando cada uno de los dos modelos de *quenching*: el
**valor constante de DAMA/LIBRA** (`QF_DAMA/`) y el **dependiente de la energía** medido para
ANAIS (`QF_Tamara/`). Cada carpeta contiene los notebooks por interacción, las listas de
energías y los contornos de DAMA y COSINE para comparar.

### `Final/`
**Curvas de exclusión finales** del trabajo para las tres interacciones, con el intervalo
óptimo y comparando ambos modelos de *quenching* (`SI_FINAL.ipynb`, `SDp_FINAL.ipynb`,
`SDn_FINAL.ipynb`).

### `Migdal/`
Estudio del **efecto Migdal**: se compara el ritmo del retroceso nuclear (NR) con el de la
señal Migdal para distintas masas de WIMP (`migdal.ipynb`, figura `rate_migdal.pdf`, datos en
los `.csv` locales). Los límites resultantes se recogen en los ficheros `Final/*_Migdal_*`.

### `Exposición/`
Estudio de la dependencia de la sensibilidad con la **exposición** del experimento
(`Plot_Exposición.ipynb`).

### `RitmoPorDetector/`
Representación del **ritmo diferencial medido por cada uno** de los nueve detectores
(`Ritmo_Por_Detector.ipynb`).

### `PlanoGalactico/`
Figura ilustrativa del **"viento de WIMPs"**: representa el esquema de velocidades que da lugar a la modulación anual de la señal
(junio/diciembre). Es material gráfico de apoyo (`wimp_wind.ipynb`, `wimp_wind_2.pdf`).