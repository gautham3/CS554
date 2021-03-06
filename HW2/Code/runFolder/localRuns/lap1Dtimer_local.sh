#!/bin/bash

N=10	## Size of max 1D data array=2^N
ntimer=1000	## no. of 1D Lap solves done to avg time
p=4 	## No. of MPI ranks to use

mpicc ../../Drivers/lap1D_par.c -lm -o lap1D_par.exe 
echo $'lap1D_par.c Compiled with mpicc \n'

mpirun --oversubscribe -np $p lap1D_par.exe $N $ntimer
echo $'Job Run using mpirun - completed \n'

rm lap1D_par.exe
