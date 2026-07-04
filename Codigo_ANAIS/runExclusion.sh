#!/bin/bash
# ===========================================================================
#  runExclusion.sh
# ---------------------------------------------------------------------------
#  Ejecuta fitSimulMakeExclusion sobre TODAS las combinaciones de:
#    intervalos de energia  x  thmodel  x  spinModel  x  qfModel  x  ANOD  x  resolution
#
#  Orden de los argumentos del ejecutable:
#    ./fitSimulMakeExclusion cl eneIni eneEnd thmodel spinModel qfModel ANOD resolution_p
#
#  Para cambiar que se recorre, edita solo las listas de la seccion CONFIG.
# ===========================================================================

# --- Detiene el script si un comando falla (quita '-e' si prefieres continuar) ---
set -u

# ===========================================================================
#  CONFIG - edita aqui lo que quieras barrer
# ===========================================================================

# Confidence Level (fijo)
CL=90

# Intervalos de energia [min max]: una entrada "min max" por linea
INTERVALS=(
  "1.0 6.0"
  "2.0 6.0"
)

# Modelos teoricos: 0-DMAnalysis 1-Python 2-RAPIDD 3-WIMPYDD 4-DMAnalysis(archivo) 5-Migdal
THMODELS=(0 1 2 3 4 5)

# Modelos de spin: 0-SI 1-SD-proton 2-SD-neutron
SPINMODELS=(0 1 2)

# Modelos de quenching factor: 1-DAMA 2-ANAIS CTE 3-TAMARA
QFMODELS=(1 2 3)

# Poblacion ALE (ANOD): 0-no 1-si
ANODS=(0 1)

# Tratamiento de la resolucion: 0-sin resolucion 1-gausiana 2-convolucion
RESOLUTIONS=(0 1 2)

# Ejecutable a llamar
EXE=./fitSimulMakeExclusion

# ===========================================================================
#  BUCLE - no suele hacer falta tocar nada de aqui para abajo
# ===========================================================================

# Contadores para el resumen final
total=0
ok=0
fail=0

for interval in "${INTERVALS[@]}"; do
  # Separa "min max" en dos variables
  read -r eneIni eneEnd <<< "$interval"

  for thmodel in "${THMODELS[@]}"; do
    for spin in "${SPINMODELS[@]}"; do
      for qf in "${QFMODELS[@]}"; do
        for anod in "${ANODS[@]}"; do
          for res in "${RESOLUTIONS[@]}"; do

            total=$((total + 1))

            echo "=================================================================="
            echo " [$total] cl=$CL ene=[$eneIni,$eneEnd] th=$thmodel spin=$spin qf=$qf anod=$anod res=$res"
            echo "=================================================================="

            "$EXE" "$CL" "$eneIni" "$eneEnd" "$thmodel" "$spin" "$qf" "$anod" "$res"

            # Comprueba si la ejecucion termino bien
            if [ $? -eq 0 ]; then
              ok=$((ok + 1))
            else
              fail=$((fail + 1))
              echo " >>> AVISO: fallo en la combinacion [$total]"
            fi

          done
        done
      done
    done
  done
done

echo "=================================================================="
echo " Terminado. Total: $total   OK: $ok   Fallidas: $fail"
echo "=================================================================="
