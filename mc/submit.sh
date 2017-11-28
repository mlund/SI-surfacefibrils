#!/bin/bash
#SBATCH -p snic
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A snic2017-1-48
#
# job time, change for what your job requires 
#SBATCH -t 20:00:00
# 
# job name
#SBATCH -J fibril
#
# filenames stdout and stderr - customise, include %j
#SBATCH -o fibrils.out
#SBATCH -e fibrils.err

module purge
module add GCC/6.2.0-2.27
module add CMake

../../cc > out