# ANAIS-112 · Análisis dependiente del modelo (TFG)

Código, datos y notebooks empleados en el Trabajo de Fin de Grado de Física
**«Búsquedas de materia oscura con ANAIS-112»** (Universidad de Zaragoza, junio de 2026).

El objetivo del trabajo es obtener **curvas de exclusión** sobre la sección eficaz de
interacción WIMP–nucleón a partir de los datos del experimento ANAIS-112, mediante un
análisis dependiente del modelo. Se implementa numéricamente el modelo de interacción
WIMP–núcleo para cristales de NaI(Tl), distinguiendo las contribuciones independiente (SI)
y dependiente del espín (SD, en protón y neutrón), y se estudia el impacto del factor de
*quenching* y del modelo de fondo sobre los límites obtenidos.

- **Autor:** Borja Jimeno Altelarrea
- **Directoras:** Dra. María Lucía Martínez Pérez · Carmen Seoane Herce
- **Departamento de Física Teórica — Facultad de Ciencias — Universidad de Zaragoza**

---

## Requisitos

El núcleo del análisis está escrito en **Python 3** con Jupyter. Dependencias principales:

- `numpy`, `scipy`, `matplotlib`
- `uproot` (lectura de archivos ROOT sin necesidad de ROOT instalado)

Instalación rápida:

```bash
pip install numpy scipy matplotlib uproot jupyter
```

Las librerías externas de validación (`WIMPyDD` y `RAPIDD`, incluidas en `Librerias/`)
tienen sus propias dependencias y, en el caso de RAPIDD, requieren compilar una parte en
C con CMake. Consulta los README propios dentro de cada carpeta.

---

## Estructura del repositorio

### Raíz — flujo principal del análisis

| Archivo | Descripción |
|---|---|
| `funciones_ritmo_teorico.py` | Módulo central del **ritmo teórico**. Implementa el modelo de halo estándar (SHM) y la función de velocidad inversa media ($\eta$), el factor de forma de Helm (`FF`), el ritmo diferencial WIMP–núcleo (`rate`, `RateNaI`), la conversión a energía equivalente en electrones mediante el factor de *quenching* (`getQFNa`, `getQFI`, `rate_ee`) y el cálculo del número de cuentas esperado (`numero_cuentas_teo`). Incluye las constantes físicas y los parámetros del modelo de halo. |
| `funciones_ritmo_exp.py` | Módulo del **ritmo experimental**. Carga los histogramas medidos por ANAIS-112, integra el ritmo en la región de interés, gestiona las exposiciones por detector y calcula el número de cuentas observado y su límite superior a un nivel de confianza dado (`numero_cuentas_exp_CL`). |
| `TFG_Calculo_N_Teo.ipynb` | Notebook de cálculo del **número de eventos teórico** esperado para cada masa de WIMP e interacción. |
| `TFG_Calculo_N_Exp.ipynb` | Notebook de cálculo del **número de eventos experimental** a partir de los datos medidos. |
| `TFG_Plot_Exclusion.ipynb` | Notebook que combina ambos para construir y representar las **curvas de exclusión** en el plano (σ, mχ). |

### `Datos/`
Histogramas de bajo fondo de ANAIS-112 en formato ROOT (`BEhistos_year123456.root`),
correspondientes a los primeros seis años de toma de datos, promediados sobre los nueve
detectores. Es la entrada experimental del análisis.

### `Codigo_ANAIS/`
`fitSimulMakeExclusion.cpp` — rutina en C++ (ROOT) empleada para el ajuste del espectro
medido con el modelo de fondo más la señal de materia oscura, base del método que
incorpora el modelo de fondo.

### `Librerias/`
Herramientas de generación de ritmos y de validación cruzada.

- **`GeneraRitmo/`** — Generación de listas de ritmo esperado (SI y SD) y comparación con
  códigos de referencia. Contiene la subcarpeta `Libs/` con:
  - `ANAIS/` — funciones propias de ritmo teórico y experimental (versión de librería).
  - `WIMPYDD/` — copia de la librería **WIMPyDD** (teoría efectiva de campos no relativista),
    usada para validar los códigos propios.
  - `RAPIDD/` — copia de la librería **RAPIDD**, usada igualmente como validación independiente.
- **`Comparacion/`** — Notebooks y datos para comparar los resultados propios con WIMPyDD/RAPIDD
  y con la señal de DAMA/LIBRA, y para optimizar el intervalo de integración.

### `Presentacion/`
Notebooks organizados según las **etapas sucesivas del análisis** descritas en la memoria,
cada uno con sus datos y figuras (`.svg`):

- `Preliminar/` — método simplificado (toda la señal atribuida a WIMPs).
- `Fondo/` — incorporación del modelo de fondo.
- `Intervalo/` — optimización del intervalo de integración.
- `ALE/` — inclusión de la población de eventos anómalos (*Anomalous Light Events*).
- `FactorQ/`, `QF_DAMA/`, `QF_Tamara/` — estudio de los distintos modelos de factor de
  *quenching* (constante de DAMA vs. dependiente de la energía).
- `Final/` — curvas de exclusión finales para SI, SDp y SDn.
- `Exposición/`, `RitmoPorDetector/` — estudios de exposición y ritmo por detector.
- `Comparacion(ohare)/` — comparación con el compendio de límites de otros experimentos.

---

## Uso rápido

1. Instala las dependencias (ver arriba).
2. Abre los notebooks de la raíz en orden: primero `TFG_Calculo_N_Teo.ipynb` y
   `TFG_Calculo_N_Exp.ipynb`, y después `TFG_Plot_Exclusion.ipynb`.
3. Los módulos `funciones_ritmo_teorico.py` y `funciones_ritmo_exp.py` se importan
   automáticamente desde los notebooks (`from funciones_ritmo_teorico import *`).
4. Para reproducir una etapa concreta del análisis, entra en la subcarpeta correspondiente
   de `Presentacion/` y ejecuta su notebook.



---

## Licencia

Este repositorio se distribuye bajo la **GNU Affero General Public License v3.0**
(ver el archivo [`LICENSE`](LICENSE)). Permite usar, modificar y redistribuir el código,
también con fines comerciales, siempre que se mantenga la atribución y que cualquier obra
derivada se publique bajo la misma licencia, con el código fuente disponible.

