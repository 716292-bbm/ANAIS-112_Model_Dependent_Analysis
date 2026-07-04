# `fitSimulMakeExclusion.cpp`

Programa en **C++ / ROOT** que construye las curvas de exclusión de ANAIS-112 mediante un
**ajuste simultáneo** (*simultaneous fit*) del espectro medido con RooFit. Para cada masa
de WIMP $m_\chi$ recorre el espectro de bajo fondo de los nueve detectores, ajusta un modelo
de fondo más una señal teórica de materia oscura, y despeja la sección eficaz límite a un
nivel de confianza dado. El resultado es un `TGraph` en el plano $(\sigma, m_\chi)$.

El código admite la interacción independiente del espín (SI) y las dependientes del espín
en protón (SDp) y neutrón (SDn), varios modelos de factor de *quenching*, la inclusión
opcional de la población ALE, y —lo más relevante para este trabajo— la posibilidad de
tomar el ritmo teórico desde distintos códigos externos (**el código propio de ANAIS,
RAPIDD, WIMPyDD** o  **DMAnalysis** ).

# Partes del código desarrolladas en este trabajo:

---

## 1. Carga del ritmo teórico desde código externo (parte desarrollada)

El espectro esperado de materia oscura no se calcula dentro de este programa, sino que se
lee **ya generado** desde archivos ROOT producidos por cada código de cálculo de ritmo.
Esto permite validar el código propio contra herramientas de referencia sin cambiar el
resto del análisis y en un futuro introducir nuevos análisis (por ejemplo otros operadores de teoría efectiva)

### Selección del código: la variable `thmodel`

El modelo teórico se elige con el argumento de línea de comandos `argv[4]`, que se guarda
en la variable entera `thmodel`. Su valor selecciona el nombre del archivo ROOT del que se
leerá el ritmo:

| `thmodel` | Código de origen | Archivo ROOT |
|:---:|---|---|
| `0` | `DMAnalysis` interno (vía `DMModelGetRate`) | `RAPIDD_SI_TH1D.root` |
| `1` | **Código desarrollado en este trabajo** (Python) | `ANAIS_SI_TH1D.root` |
| `2` | **RAPIDD** | `RAPIDD_SI_TH1D.root` |
| `3` | **WIMPyDD** | `WIMPYDD_SI_TH1D.root`, `WIMPYDD_SDp_TH1D.root` o `WIMPYDD_SDn_TH1D.root` |
| `4` | `DMAnalysis` por archivo | `DMA_SI_TH1D.root` |
| `5` | Migdal | `MIGDAL.root` |

En el caso de WIMPyDD (`thmodel == 3`) el nombre del archivo depende además del tipo de
interacción (`SpinModel`), de modo que se cargan histogramas distintos para SI, SDp y SDn.
Esta cadena de `if` que asigna `fileName` es el punto donde se decide, en la práctica, qué
código de ritmo se está usando.

El sufijo del nombre del `TGraph` de salida también refleja el código empleado
(`_PY` para Python/ANAIS, `_RA` para RAPIDD, `_WI` para WIMPyDD, `_DM_SG`/`_DMA` para
DMAnalysis), lo que facilita comparar las curvas de exclusión obtenidas con cada uno.

### Lectura y reconstrucción del espectro: `DMModelGetRateEeeFromFile`

Salvo el caso `thmodel == 0` (que llama a `DMModelGetRate` y calcula el ritmo con la
librería `DMRate` interna), todos los modelos externos se cargan con la función:

```cpp
TH1D *DMModelGetRateEeeFromFile(std::string fileName,
                                std::string nameNa,
                                std::string nameI)
```

Esta función es la interfaz común con los códigos externos. Su lógica es:

1. **Abre el archivo ROOT** indicado por `fileName` y recupera dos histogramas: uno con el
   ritmo del **sodio** (`nameNa`) y otro con el del **yodo** (`nameI`), cada uno en función de la energía de retroceso nuclear. Los nombres de estos histogramas se construyen antes, en el bucle principal, con el patrón `hist_Na_mw_<masa>` / `hist_I_mw_<masa>`, de forma que cada masa de WIMP tiene su par de histogramas.

2. **Convierte la energía de retroceso nuclear a energía equivalente en electrones**
   ($E_{ee}$), que es la magnitud que mide realmente el detector, aplicando el factor de *quenching* $Q$ elemento por elemento.E l *quenching* puede ser constante (`QNa`, `QI`) o dependiente de la energía si se ha cargado como `TGraph` (`gQNa`, `gQI`), en cuyo caso se evalúa en cada
   punto con `gQNa->Eval(Eee)`.

3. **Combina sodio y yodo** en el ritmo del cristal de NaI, ponderando por el número másico
   de cada especie. Devuelve un único histograma `hWimp` con el ritmo del NaI en $E_{ee}$, listo para
   convolucionar con la resolución.

De este modo, añadir un nuevo código de cálculo de ritmo se reduce a generar un archivo
ROOT con los histogramas de Na e I por masa y añadir una entrada en la selección de
`fileName`; el resto del ajuste no cambia.

---

## 2. Convolución gaussiana con la resolución del detector (parte desarrollada)

El ritmo leído del archivo hay que **convolucionarlo con una gaussiana** (cuya anchura depende de la energía) antes de ajustarlo a los datos.

### Selección del método: `resolution_p`

La variable `resolution_p` controla cómo se aplica la resolución:

- `0` — no se aplica resolución.
- `1` — resolución por **muestreo Monte Carlo**: se sortean $10^7$ eventos del espectro con
  `GetRandom()` y cada uno se desplaza con `TRandom::Gaus(ee, sigma)`, rellenando el
  histograma resultante. Después se renormaliza para conservar la integral.
- `2` — resolución por **convolución numérica directa** (el método usado), implementado en
  la función `Conv2`.

En ambos casos la anchura de la gaussiana sigue la forma
$\sigma_E = p_1 + p_s\sqrt{E}$, con los parámetros $p_1$ y $p_s$ (`param_1`, `param_2`)
leídos de un ajuste de resolución almacenado en `fitsResolution.root` (función `fresD0`).

### La función `Conv2`

```cpp
int Conv2(double e1, double e2, double ebin, double p1, double ps,
          double *S0, TH1 *ritmo_sr)
```

Realiza la convolución del espectro `ritmo_sr` con la gaussiana de resolución, devolviendo
en el array `S0` el ritmo emborronado en cada bin de salida. Los pasos:

1. **Array auxiliar de ritmos.** Se muestrea el espectro de entrada en un array denso de
   `arrayDim = 400` puntos, extendido más allá del rango de interés en $\pm 5\sigma$
   (`nSig = 5`) para que las colas de la gaussiana no se trunquen. La anchura local se
   estima con $\sigma = |p_1 + p_s\sqrt{E}|$.

2. **Triple bucle de integración.** Para cada bin de salida (bucle en `ei`):
   - se subdivide el bin en 20 pasos (bucle en `edif`) para integrar dentro del propio bin;
   - en cada punto se integra la gaussiana sobre un intervalo de $\pm 5\sigma$ (bucle en `en`), con peso multiplicado por el valor del ritmo interpolado del array auxiliar.
   - El valor convolucionado se **normaliza** dividiendo por la suma de pesos (`NORM2`),
     lo que garantiza que la convolución conserva el número de cuentas.

3. **Caso $\sigma = 0$.** Si la anchura es nula se toma directamente el valor del array sin
   convolucionar, evitando la división por cero.

El resultado (`rate_res`) se vuelca en el histograma `hWimp`, que ya representa el espectro
teórico tal como lo vería el detector. Este es el histograma que se convierte en
`RooHistPdf` (`pdf_Wimp`) y entra en el ajuste.

---

## 3. Resto del programa (resumen)

### Configuración y argumentos de entrada
`main` recibe por línea de comandos: nivel de confianza (`cl`), energías mínima y máxima
del intervalo de ajuste (`min`, `max`), el modelo teórico (`thmodel`), y opcionalmente el
modelo de *quenching* (`qf`) y la inclusión de ALE. El tipo de interacción está fijado en
el código con `SpinModel` (0 = SI, 1 = SDp, 2 = SDn).

### Factor de *quenching*
Según `qfModel` se establece: `1` → valores constantes de DAMA ($Q_{Na}=0.3$, $Q_I=0.09$);
`2` → constante de ANAIS; en otro caso se cargan las curvas dependientes de la energía
(`gNa`, `gI`) desde `QFTamara.root`. Estas curvas son las que usa la conversión a $E_{ee}$
descrita en la sección 1.

### Lectura de datos y exposición
Se fijan los tiempos vivos de los nueve detectores y se calcula la exposición
(tiempo vivo × 12.5 kg). El espectro medido se lee de `BEhistos_year123456.root` y se
carga como `RooDataHist` con una categoría por detector (`RooCategory`), preparándolo para
el ajuste simultáneo.

### Modelo de fondo
Se lee el modelo de fondo por detector (con o sin ALE según el sufijo `_conANOD`) y se
convierte en `RooHistPdf`. Su normalización se fija a la integral del fondo simulado y se
mantiene **constante** en el ajuste, de modo que el único parámetro libre es la
normalización de la señal WIMP (`nNorm`).

### Bucle en masas y ajuste
Para cada $m_\chi$ (el vector `mw` cubre de ~2 GeV a $10^4$ GeV) se construye, por detector,
el modelo `nNorm·pdf_Wimp + nbkg·pdf_bkg` y se combinan en un `RooSimultaneous`. El ajuste
se hace con `fitTo(...)` de forma extendida sobre el rango de energía elegido.

### Cálculo de la sección eficaz límite
Del ajuste se obtiene `nNorm` y su error; se aplica el factor del nivel de confianza
(1.28 para 90 %, 1.64 para 95 %) y se divide por la integral del espectro teórico para
despejar la sección eficaz. Finalmente se convierte de pb a cm² ($\times 10^{-36}$) y se
normaliza por la exposición total. Existe además una rama (con `min <= 0`) que busca el
**intervalo de integración óptimo** barriendo subintervalos entre 1 y 6 keVee.

### Salida
Los pares $(m_\chi, \sigma)$ se guardan en un `TGraph` cuyo nombre codifica todas las
opciones del análisis (años, CL, intervalo, modelo de QF, ALE, código teórico y tipo de
interacción), se escribe en un archivo ROOT y se representa en un lienzo con ejes
logarítmicos y un eje secundario en pb.

---

## Uso

```bash
fitSimulMakeExclusion <cl> <eneIni> <eneEnd> <thmodel> [qf] [includeALE]
```

- `cl` — nivel de confianza (p. ej. `90` o `95`).
- `eneIni`, `eneEnd` — límites del intervalo de ajuste en keVee (usar `eneIni <= 0` para
  buscar el intervalo óptimo).
- `thmodel` — código del ritmo teórico: `1` ANAIS, `2` RAPIDD, `3` WIMPyDD, `4` DMAnalysis,
  `5` Migdal.
- `qf` *(opcional)* — `1` DAMA, `2` constante ANAIS, cualquier otro valor → dependiente de
  la energía.
- `includeALE` *(opcional)* — incluye la población ALE.

