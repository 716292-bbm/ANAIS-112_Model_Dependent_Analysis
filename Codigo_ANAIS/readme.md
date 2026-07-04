# fitSimulMakeExclusion — Documentación

Programa para calcular la **curva de exclusión** de materia oscura (sección eficaz WIMP-nucleón límite frente a la masa del WIMP) del experimento **ANAIS**, mediante un ajuste simultáneo con RooFit de los 9 detectores.

---

## 1. Qué hace, en resumen

Para cada masa de WIMP de una lista, el programa:

1. Obtiene el **espectro teórico** (ritmo esperado de eventos de materia oscura), ya sea calculándolo con `DMRate` o leyéndolo de un fichero externo.
2. Le aplica la **resolución energética** del detector (por convolución, por Monte Carlo, o nada).
3. Realiza un **ajuste simultáneo** de los datos experimentales en los 9 detectores con el modelo `datos = nNorm·señal + nbkg·fondo`, dejando libre únicamente `nNorm` (proporcional a la sección eficaz).
4. Del valor y el error de `nNorm` deriva la **sección eficaz límite** al nivel de confianza pedido.

Al final construye la curva σ(mW), la guarda en un fichero ROOT y la dibuja en un PNG.

---

## 2. Inputs (argumentos de línea de comandos)

El programa se ejecuta con **8 argumentos obligatorios**, todos posicionales y en este orden:

```
./fitSimulMakeExclusion cl eneIni eneEnd thmodel spinModel qfModel ANOD resolution_p
```

| # | Argumento | Tipo | Significado | Valores |
|---|-----------|------|-------------|---------|
| 1 | `cl` | double | Nivel de confianza del límite | `90` o `95` |
| 2 | `eneIni` | double | Energía mínima del ajuste (keV) | p.ej. `1.0` |
| 3 | `eneEnd` | double | Energía máxima del ajuste (keV) | p.ej. `6.0` |
| 4 | `thmodel` | int | Modelo teórico / origen del espectro de señal | `0`-`5` (ver abajo) |
| 5 | `spinModel` | int | Tipo de acoplo | `0`-SI, `1`-SD protón, `2`-SD neutrón |
| 6 | `qfModel` | int | Modelo de *quenching factor* | `1`-DAMA, `2`-ANAIS cte, `3`-TAMARA |
| 7 | `ANOD` | int (bool) | Incluir población ALE en el fondo | `0`-no, `1`-sí |
| 8 | `resolution_p` | int | Tratamiento de la resolución | `0`-ninguna, `1`-gaussiana (MC), `2`-convolución |

Si se pasan **menos de 8 argumentos**, el programa imprime la ayuda de uso y termina con código `1`.

### Detalle del argumento `thmodel`

| Valor | Modelo | Fichero de señal usado |
|-------|--------|------------------------|
| `0` | DMAnalysis (calculado en vivo con `DMRate`) | — (no lee fichero) |
| `1` | Python | `ANAIS_SI_TH1D.root` |
| `2` | RAPIDD | `RAPIDD_SI_TH1D.root` |
| `3` | WIMPYDD | `WIMPYDD_SI/SDp/SDn_TH1D.root` (según spin) |
| `4` | DMAnalysis (por archivo) | `DMA_SI_TH1D.root` |
| `5` | Migdal | `MIGDAL.root` |

> Nota: solo `thmodel == 0` calcula el espectro internamente con `DMRate`. El resto lo **lee de un fichero** con `DMModelGetRateEeeFromFile`.

### Detalle del argumento `qfModel`

El *quenching factor* (QF) relaciona la energía de retroceso nuclear con la energía electrón-equivalente medida: 

| Valor | Modelo | QF Na | QF I |
|-------|--------|-------|------|
| `1` | DAMA | 0.30 (cte) | 0.09 (cte) |
| `2` | ANAIS cte | 1.0 (cte) | 1.0 (cte) |
| `3` | TAMARA | QF(E) leído de `QFTamara.root` (`gNa`) | QF(E) leído de `QFTamara.root` (`gI`) |

### Detalle del argumento `resolution_p`

| Valor | Método | Cómo funciona |
|-------|--------|---------------|
| `0` | Sin resolución | El espectro se usa tal cual |
| `1` | Gaussiana (Monte Carlo) | Genera 10⁷ eventos y los desplaza con σ(E); reconstruye el espectro |
| `2` | Convolución | Convoluciona el espectro con la gaussiana σ(E) mediante `Conv2` (más rápido y sin ruido) |

En los modos 1 y 2, la anchura de la resolución es `σ(E) = p1 + p2·√E`, con `p1`, `p2` leídos de un ajuste previo (`fitsResolution.root` para el modo 2).

---

## 3. Ficheros de entrada requeridos

Además de los argumentos, el programa espera encontrar (rutas relativas al directorio de ejecución):

- **Datos experimentales:** `../../data/BEhistos_yearN.root` — histogramas de fondo medido, uno por detector (`hbea_...y_D0..D8`). El nombre incorpora los años seleccionados.
- **Modelo de fondo:** `../../backgroundModel/backgroundModel_single_y123456.root` (o `..._conANOD.root` si `ANOD=1`) — histogramas `hD0..hD8`.
- **Resolución:** `fitsResolution.root` (con `resolution_p=2`) o la ruta a `.../resMartaByDet/fitsResolution.root` (con `resolution_p=1`), objeto `fresD0`.
- **Quenching factor:** `QFTamara.root` (solo con `qfModel=3`), con los `TGraph` `gNa` y `gI`.
- **Espectros de señal:** el `.root` correspondiente al `thmodel` (ver tabla), salvo `thmodel=0`.
- **Para `DMRate` (`thmodel=0`):** `rate.dat`.
- **Fondo del plot:** `plots/SI_mw_em1_e4_sdp_5em50_em37.JPG` (imagen de referencia).

---

## 4. Outputs

- **Fichero ROOT:** `plots/SI_varios_2.root` (abierto en modo `update`). Dentro se guarda un `TGraph` con la curva σ(mW). El **nombre del objeto** codifica toda la configuración usada, por ejemplo:

  ```
  gA112_6y_90_1.000000_6.000000_QFdama_ANOD_DM_SG_rp2_SI
  ```

  El nombre se compone de: prefijo + CL + rango de energía + modelo de QF + (ANOD) + modelo teórico + modo de resolución (`_rp0/_rp1/_rp2`) + modelo de spin. Esto permite guardar muchas curvas en el mismo fichero sin que se pisen.

- **Figura PNG:** `plots/<mismo_nombre>.png` — la curva de exclusión dibujada sobre la imagen de referencia, con ejes logarítmicos y un eje secundario en pb.

- **Salida por pantalla:** resumen de la configuración, cuentas al CL para cada masa, exposición total y la tabla final `mW → sigma`.

---

## 5. Funciones del programa

### `int Conv2(e1, e2, ebin, p1, ps, S0, ritmo_sr)`

Convoluciona un espectro con la resolución gaussiana del detector.

- **Entrada:** rango `[e1, e2]` y binado `ebin` de salida; parámetros de resolución `p1`, `ps` (σ(E) = |p1 + ps·√E|); histograma de entrada `ritmo_sr`.
- **Salida:** rellena el array `S0` con el ritmo convolucionado. Devuelve `0` siempre.
- **Cómo funciona:** precalcula el ritmo de entrada en un array fino (para no llamar a `GetBinContent` en los bucles internos), y para cada bin de salida integra el producto `ritmo(en)·gauss(edif−en)` sobre una ventana de ±5σ, normalizando por el peso gaussiano acumulado. Es el método usado cuando `resolution_p == 2`.

### `TH1F* DMModelGetRate(mw, sigma, qfModel, SpinModel)`

Calcula el espectro teórico de señal WIMP **en vivo** con la librería `DMRate` (camino de `thmodel == 0`).

- **Entrada:** masa del WIMP `mw` (GeV); `sigma` (normalización, se suele pasar 1); `qfModel` (no se usa directamente aquí); `SpinModel`.
- **Salida:** histograma `TH1F` con el ritmo (1000 bins), en energía de retroceso nuclear.
- **Cómo funciona:** inicializa `DMRate` desde `rate.dat`, fija la masa y, según `SpinModel`, activa el acoplo SI (θ=0), SD-protón (θ=0) o SD-neutrón (θ=π/2). Asigna el QF del Na al elemento 0 y el del I al elemento 1 (constante o gráfico según corresponda). Con la macro `ARCHIVO` activada, además guarda el histograma en un `.root`.

### `TH1D* DMModelGetRateEeeFromFile(fileName, nameNa, nameI)`

Lee el espectro de señal desde un **fichero externo** (usado por todos los `thmodel` salvo el 0).

- **Entrada:** `fileName` (ROOT con los histogramas); `nameNa` y `nameI` (nombres de los histogramas de Na y I dentro del fichero).
- **Salida:** histograma `TH1D` combinado, en energía electrón-equivalente. Devuelve `nullptr` si algo falla (con mensaje de error).
- **Cómo funciona:** lee los histogramas de Na y de I (que están en energía de retroceso nuclear), y para cada bin en energía electrón-equivalente convierte con el QF (`ENR = Eee/QF`, y el ritmo se escala por `1/QF`). Combina ambos según la abundancia molar del NaI: `(23·Na + 127·I) / 150`.

### `int main(int argc, char** argv)`

Orquesta todo el cálculo. Sus bloques principales:

1. **Argumentos:** comprueba que hay 8 y los lee.
2. **Resolución:** si `resolution_p == 2`, carga `p1`, `p2` de `fitsResolution.root`.
3. **Quenching factor:** fija QF constantes o carga los gráficos según `qfModel`.
4. **Datos y categorías:** define la variable de energía de RooFit, la lista de años y detectores, y una `RooCategory` con un tipo por detector (necesaria para el fit simultáneo).
5. **Lectura de datos:** carga los histogramas de datos, los escala a número de cuentas y los mete en un `RooDataHist` por categoría.
6. **Modelo de fondo:** abre el fichero de fondo (con o sin ANOD) y crea el parámetro libre `nNorm`.
7. **Lista de masas y nombres:** construye el array de masas (lista fija para Migdal, rejilla logarítmica para el resto) y los nombres de los histogramas de señal. Selecciona el fichero de señal según `thmodel`.
8. **Bucle sobre masas:** para cada masa obtiene el espectro, le aplica la resolución, construye el modelo `nNorm·señal + nbkg·fondo` por detector (con `nbkg` **fijo** y `nNorm` **libre**), hace el ajuste extendido simultáneo, y calcula la sección eficaz límite (factor 1.64 para 95 %, 1.28 para 90 %), normalizada por la exposición total.
9. **Resultados y gráfica:** imprime la tabla, construye el nombre codificado, guarda el `TGraph` en el `.root` y dibuja el PNG.

---

## 6. Cómo funciona `runExclusion.sh`

Es un script auxiliar que **automatiza el barrido** de muchas ejecuciones del programa sin tener que escribirlas a mano. Ejecuta `fitSimulMakeExclusion` sobre **todas las combinaciones** de los parámetros que se le indiquen.

### Estructura

Toda la configuración está en listas al principio del script (sección `CONFIG`), que son lo único que hay que editar:

```bash
CL=90                              # Nivel de confianza (fijo)

INTERVALS=(                        # Intervalos de energia [min max]
  "1.0 6.0"
  "2.0 6.0"
)

THMODELS=(0 1 2 3 4 5)             # Modelos teoricos a recorrer
SPINMODELS=(0 1 2)                 # Modelos de spin
QFMODELS=(1 2 3)                   # Modelos de quenching factor
ANODS=(0 1)                        # Poblacion ALE (no / si)
RESOLUTIONS=(0 1 2)               # Tratamiento de la resolucion

EXE=./fitSimulMakeExclusion       # Ejecutable a llamar
```

### Qué hace

El script anida bucles (`for`) sobre cada lista y, para **cada combinación**, lanza:

```bash
./fitSimulMakeExclusion $CL $eneIni $eneEnd $thmodel $spin $qf $anod $res
```

Es decir, recorre **intervalos × thmodel × spin × qf × ANOD × resolución**. Con las listas de ejemplo son `2 × 6 × 3 × 3 × 2 × 3 = 648` ejecuciones.

Antes de cada ejecución imprime una cabecera indicando qué combinación se está corriendo, y comprueba el código de salida para llevar la cuenta de cuántas terminaron bien y cuántas fallaron. Al final imprime un resumen:

```
Terminado. Total: 648   OK: 645   Fallidas: 3
```

### Cómo usarlo

```bash
bash runExclusion.sh
```

