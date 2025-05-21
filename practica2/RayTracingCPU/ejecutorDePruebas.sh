#!/bin/bash

RESOLUCIONES=("80" "640" "1200")
NS_VALUES=("10" "50")
MPI_PROCS=("2" "4" "10")
OMP_THREADS=("2" "8" "16")

for res in "${RESOLUCIONES[@]}"; do
  for ns in "${NS_VALUES[@]}"; do
    for mpi in "${MPI_PROCS[@]}"; do
      for omp_threads in "${OMP_THREADS[@]}"; do

        echo "Ejecutando configuración: ${res}x${res}, ${mpi} procesos MPI, ${omp_threads} hilos OMP, ${ns} muestras"

        CONFIG_DIR="resultados/${res}x${res}_${mpi}p_${omp_threads}h_${ns}ns"
        mkdir -p "$CONFIG_DIR"

        RESULT_FILE="${CONFIG_DIR}/tiempos.txt"
        echo "Tiempos de ejecución para ${res}x${res}, ${mpi} procesos MPI, ${omp_threads} hilos OMP, ${ns} ns" > "$RESULT_FILE"
        echo "Iteración | Serie | OMP Columnas | OMP Filas | OMP Bloques | MPI Columnas | MPI Filas | MPI Bloques | MPI+OMP Columnas | MPI+OMP Filas | MPI+OMP Bloques" >> "$RESULT_FILE"

        for iter in {1..10}; do
          echo "  → Iteración $iter"
          export OMP_NUM_THREADS=$omp_threads

          # Ejecutar y capturar la salida del programa
          OUTPUT=$(mpirun -np $mpi ./rayTracerCompleto $res $res $ns) #MODIFICAR AQUÍ SI SE LE HA DADO OTRO NOMBRE AL COMPILAR 

          # Mover todas las imágenes .bmp a la carpeta de la configuración (se sobreescriben)
          for img in *.bmp; do
            if [ -f "$img" ]; then
              mv "$img" "$CONFIG_DIR/"
            fi
          done

          # Extraer tiempos de cada estrategia
          serie=$(echo "$OUTPUT" | grep "repartición en serie" | grep -oE "[0-9]+\.[0-9]+")
          omp_col=$(echo "$OUTPUT" | grep "columnas con OMP" | grep -oE "[0-9]+\.[0-9]+")
          omp_fil=$(echo "$OUTPUT" | grep "filas con OMP:" | grep -oE "[0-9]+\.[0-9]+")
          omp_blo=$(echo "$OUTPUT" | grep "bloques con OMP:" | grep -oE "[0-9]+\.[0-9]+")
          mpi_col=$(echo "$OUTPUT" | grep "columnas con MPI:" | grep -oE "[0-9]+\.[0-9]+")
          mpi_fil=$(echo "$OUTPUT" | grep "filas con MPI:" | grep -oE "[0-9]+\.[0-9]+")
          mpi_blo=$(echo "$OUTPUT" | grep "bloques con MPI:" | grep -oE "[0-9]+\.[0-9]+")
          mpimp_col=$(echo "$OUTPUT" | grep "columnas con MPI y OMP:" | grep -oE "[0-9]+\.[0-9]+")
          mpimp_fil=$(echo "$OUTPUT" | grep "filas con MPI y OMP:" | grep -oE "[0-9]+\.[0-9]+")
          mpimp_blo=$(echo "$OUTPUT" | grep "bloques con MPI y OMP:" | grep -oE "[0-9]+\.[0-9]+")

          # Escribir resultados en el fichero por configuración
          echo "$iter | $serie | $omp_col | $omp_fil | $omp_blo | $mpi_col | $mpi_fil | $mpi_blo | $mpimp_col | $mpimp_fil | $mpimp_blo" >> "$RESULT_FILE"
        done 
      done
    done
  done
done

