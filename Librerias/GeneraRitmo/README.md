# `Librerias/GeneraRitmo/`

Núcleo de la **generación del ritmo teórico** de materia oscura. Aquí se calcula el espectro
esperado en ANAIS-112 para cada masa de WIMP y tipo de interacción (SI, SDp, SDn), usando
tanto el código propio de ANAIS como las librerías externas de referencia, y se exporta al
formato ROOT que consume el ajuste de exclusión.

**Conexión con el resto del proyecto:** los archivos `*_TH1D.root` que se generan aquí (en
`Results/`) son exactamente la entrada teórica que lee
`Codigo_ANAIS/fitSimulMakeExclusion.cpp` a través de su variable `thmodel`. Esta carpeta es,
por tanto, el paso previo a la construcción de las curvas de exclusión.

## Notebooks de generación

- **`GenerarLista_SI.ipynb`** — Genera las listas de ritmo para la interacción independiente
  del espín (SI). `GenerarLista_SI(prueba).ipynb` es una versión de pruebas.
- **`GenerarLista_SD.ipynb`** — Genera las listas de ritmo para las interacciones dependientes
  del espín (SDp y SDn).
- **`Genera_Lista_WIMPYDD_SD.ipynb`** — Genera las listas SD específicamente con WIMPyDD.
- **`WIMPYDD.ipynb`, `WIMPYDD_SD.ipynb`, `WIMPYDD_SD_NREFT.ipynb`** — Cálculo del ritmo con la
  librería WIMPyDD para SI, SD y en el marco de la teoría efectiva no relativista (NREFT).
- **`NREFT_example.ipynb`** — Ejemplo de uso de la teoría efectiva de campos no relativista.
- **`Abrir_Root.ipynb`** — Utilidad para inspeccionar los archivos ROOT generados.
- **`RepresentarResultados SI.ipynb`** — Representación de los ritmos/resultados SI.

## Subcarpetas y archivos

- **`Libs/`** — Copias de las librerías empleadas:
  - `ANAIS/` — código propio (`funciones_ritmo_teorico.py`, `funciones_ritmo_exp.py`) en
    forma de librería importable.
  - `WIMPYDD/` — librería WIMPyDD (software de terceros, con su propio README).
  - `RAPIDD/` — librería RAPIDD (software de terceros; requiere compilar con CMake, ver su
    README).
- **`Patch/package.py`** — Parche aplicado sobre una de las librerías externas para adaptarla
  al uso de este trabajo.
- **`Results/`** — Histogramas de ritmo finales en formato ROOT, listos para el ajuste:
  `ANAIS_SI_TH1D.root`, `RAPIDD_SI_TH1D.root`, `WIMPYDD_SI_TH1D.root`,
  `WIMPYDD_SDp_TH1D.root`, `WIMPYDD_SDn_TH1D.root`.
- **`rate_DMAnalysis_*.root`** — Ritmos intermedios generados con el código `DMAnalysis`.

## Requisitos

Además de `numpy`/`scipy`/`matplotlib`, el uso de las librerías externas requiere instalarlas
según sus propias instrucciones (WIMPyDD se instala como paquete de Python; RAPIDD necesita
`cmake`, las librerías `gsl` y `pkg-config`).
