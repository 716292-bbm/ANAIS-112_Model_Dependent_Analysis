# `Librerias/`

Herramientas para **generar el ritmo teórico** esperado de materia oscura y para **validar**
los códigos propios contra herramientas de referencia de la comunidad (WIMPyDD y RAPIDD).

El flujo general es: aquí se calculan los espectros teóricos de ritmo para cada masa de WIMP
y tipo de interacción, se exportan a archivos ROOT, y esos archivos son los que después
consume el ajuste de exclusión en `Codigo_ANAIS/fitSimulMakeExclusion.cpp`.

## Subcarpetas

### `GeneraRitmo/`
Núcleo de la generación de ritmos. Contiene los notebooks que calculan el espectro esperado
(SI, SDp, SDn) con el código propio de ANAIS y con las librerías externas, y los exportan al
formato de histograma ROOT (`*_TH1D.root`) que espera el programa de exclusión. Incluye la
subcarpeta `Libs/` con las copias de las librerías **ANAIS** (código propio), **WIMPyDD** y
**RAPIDD**. Ver su README propio para el detalle.

### `Comparacion/`
Notebooks y datos para **comparar** los resultados obtenidos con los distintos códigos de
ritmo entre sí y con la señal de DAMA/LIBRA, y para estudios auxiliares como la optimización
del intervalo de integración. Ver su README propio.

## Nota

Las librerías externas incluidas (**WIMPyDD** y **RAPIDD**, dentro de `GeneraRitmo/Libs/`)
son software de terceros con sus propias licencias y READMEs. Se incluyen aquí para poder
reproducir la validación cruzada de los códigos propios; no son desarrollo original de este
trabajo.
