#!/bin/bash

N=50	## Size of 1D data array
ntimer=1000	## No. of iterations for each solve
p=5 	## No. of MPI ranks to use

mpicc ../../Drivers/lap1D_par.c -lm -o lap1D_par.exe 
echo $'lap1D_par.c Compiled with mpicc \n'

mpirun --oversubscribe -np $p lap1D_par.exe $N $ntimer
echo $'Job Run using mpirun - completed \n'

rm lap1D_par.exe