#!/bin/bash

RESOLUCIONES=("80" "640" "1200")
NS_VALUES=("10" "50")

# === 1. SOLO MPI ===
MPI_PROCS=("2" "4" "10")

for res in "${RESOLUCIONES[@]}"; do
  for ns in "${NS_VALUES[@]}"; do
    for mpi in "${MPI_PROCS[@]}"; do
      echo "→ MPI: ${res}x${res}, ${mpi} procesos, ${ns} muestras"
      CONFIG_DIR="resultados_mpi/${res}x${res}_${mpi}p_${ns}ns"
      mkdir -p "$CONFIG_DIR"
      RESULT_FILE="${CONFIG_DIR}/tiempos.txt"
      echo "Tiempos de ejecución (MPI) para ${res}x${res}, ${mpi} procesos, ${ns} ns" > "$RESULT_FILE"
      echo "Iteración | MPI Columnas | MPI Filas | MPI Bloques" >> "$RESULT_FILE"

      for iter in {1..10}; do
        OUTPUT=$(mpirun -np $mpi ./rayTracerMPI $res $res $ns)
        col=$(echo "$OUTPUT" | grep "columnas con MPI" | grep -oE "[0-9]+\.[0-9]+")
        fil=$(echo "$OUTPUT" | grep "filas con MPI" | grep -oE "[0-9]+\.[0-9]+")
        blo=$(echo "$OUTPUT" | grep "bloques con MPI" | grep -oE "[0-9]+\.[0-9]+")
        echo "$iter | $col | $fil | $blo" >> "$RESULT_FILE"
      done

      echo "" >> "$RESULT_FILE"
    done
  done
done

# === 2. SOLO OMP ===
OMP_THREADS=("2" "10" "16")

for res in "${RESOLUCIONES[@]}"; do
  for ns in "${NS_VALUES[@]}"; do
    for omp in "${OMP_THREADS[@]}"; do
      echo "→ OMP: ${res}x${res}, ${omp} hilos, ${ns} muestras"
      CONFIG_DIR="resultados_omp/${res}x${res}_${omp}t_${ns}ns"
      mkdir -p "$CONFIG_DIR"
      RESULT_FILE="${CONFIG_DIR}/tiempos.txt"
      echo "Tiempos de ejecución (OMP) para ${res}x${res}, ${omp} hilos, ${ns} ns" > "$RESULT_FILE"
      echo "Iteración | OMP Columnas | OMP Filas | OMP Bloques" >> "$RESULT_FILE"

      for iter in {1..10}; do
        export OMP_NUM_THREADS=$omp
        OUTPUT=$(./rayTracerOMP $res $res $ns)
        col=$(echo "$OUTPUT" | grep "columnas con OMP" | grep -oE "[0-9]+\.[0-9]+")
        fil=$(echo "$OUTPUT" | grep "filas con OMP" | grep -oE "[0-9]+\.[0-9]+")
        blo=$(echo "$OUTPUT" | grep "bloques con OMP" | grep -oE "[0-9]+\.[0-9]+")
        echo "$iter | $col | $fil | $blo" >> "$RESULT_FILE"
      done

      echo "" >> "$RESULT_FILE"
    done
  done
done

# === 3. MPI + OMP ===
MPI_VALUES=("2" "2" "10" "10" "4")
OMP_VALUES=("2" "16" "2" "16" "8")

for res in "${RESOLUCIONES[@]}"; do
  for ns in "${NS_VALUES[@]}"; do
    for i in {0..4}; do
      mpi=${MPI_VALUES[$i]}
      omp=${OMP_VALUES[$i]}
      echo "→ MPI+OMP: ${res}x${res}, ${mpi}p, ${omp}t, ${ns} ns"
      CONFIG_DIR="resultados_mpiomp/${res}x${res}_${mpi}p_${omp}t_${ns}ns"
      mkdir -p "$CONFIG_DIR"
      RESULT_FILE="${CONFIG_DIR}/tiempos.txt"
      echo "Tiempos de ejecución (MPI+OMP) para ${res}x${res}, ${mpi}p, ${omp}t, ${ns} ns" > "$RESULT_FILE"
      echo "Iteración | MPI+OMP Columnas | MPI+OMP Filas | MPI+OMP Bloques" >> "$RESULT_FILE"

      for iter in {1..10}; do
        export OMP_NUM_THREADS=$omp
        OUTPUT=$(mpirun -np $mpi ./rayTracerMPIOMP $res $res $ns)
        col=$(echo "$OUTPUT" | grep "columnas con MPI y OMP" | grep -oE "[0-9]+\.[0-9]+")
        fil=$(echo "$OUTPUT" | grep "filas con MPI y OMP" | grep -oE "[0-9]+\.[0-9]+")
        blo=$(echo "$OUTPUT" | grep "bloques con MPI y OMP" | grep -oE "[0-9]+\.[0-9]+")
        echo "$iter | $col | $fil | $blo" >> "$RESULT_FILE"
      done

      echo "" >> "$RESULT_FILE"
    done
  done
done