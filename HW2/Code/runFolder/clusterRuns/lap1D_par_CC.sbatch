#!/bin/bash

#####################################################
##        Important parameters of your job         ##
##             are specified here 		   ##
#####################################################

#SBATCH -J Lap_1D_Timer				## job name
#SBATCH -t 00:30:00					## total computing time
#SBATCH -N 1						## number of nodes 
#SBATCH -n 4						## total number of MPI tasks
####SBATCH --mem-per-cpu=7200			## total memory required (if not specified, will get the entire node memory)
#SBATCH --partition=cs				## queue
#SBATCH --output=../Results/1Dtimer.out%j		## file that will receive output from execution
#SBATCH --error=../Results/1Dtimer.err%j			## file that will receive any error messages

#SBATCH --mail-user=gautham3@illinois.edu
#SBATCH --mail-type=begin 			## email me when the job starts
#SBATCH --mail-type=end 			## email me when the job finishes


########## Run your executable ######################
n=10			## max size of full 1-D problem, timed till here in sizes powers of 2
ntimer=1000		## no. of 1D Lap solves done to avg time

module load openmpi
srun ./lap1D_par.exe $n $ntimer
