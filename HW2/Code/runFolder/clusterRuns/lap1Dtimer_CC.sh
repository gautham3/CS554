#!/bin/bash

module load openmpi
echo "MPI module loaded"

mpicc ../../Drivers/lap1D_par.c -lm -o lap1D_par.exe
echo "lap1D_par.c Compiled with mpicc"

sbatch ./lap1D_par_CC.sbatch 
echo "Job submitted via sbatch"
