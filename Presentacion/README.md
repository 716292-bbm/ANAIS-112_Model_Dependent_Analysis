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

### `Exposición/`
Estudio de la dependencia de la sensibilidad con la **exposición** del experimento
(`Plot_Exposición.ipynb`).

### `RitmoPorDetector/`
Representación del **ritmo diferencial medido por cada uno** de los nueve detectores
(`Ritmo_Por_Detector.ipynb`).

### `Comparacion(ohare)/`
Comparación de los resultados con el **compendio público de límites** de otros experimentos
de detección directa (basado en el repositorio de datos de Ciaran O'Hare). Contiene una gran
cantidad de archivos de datos de terceros; **no es desarrollo original de este trabajo** y
podría eliminarse sin afectar al análisis propio.
