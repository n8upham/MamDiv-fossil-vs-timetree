#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=PSR-%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=25000
#SBATCH -p bigmem

module load R/3.6.1-foss-2018b-X11-20180604

Rscript ./getPulledSpeciation_100trees_forHPC.R
