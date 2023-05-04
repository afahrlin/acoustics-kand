#!/bin/bash -l
#SBATCH -A naiss2023-22-290
#SBATCH -t 4-00:00:00
#SBATCH -J wallsource
#SBATCH -n 10
#SBATCH -p core

module load matlab

matlab -nojvm -batch "main; exit"