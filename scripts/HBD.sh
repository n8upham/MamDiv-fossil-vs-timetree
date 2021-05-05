#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathan.upham@yale.edu
#SBATCH --output=HBD-%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=20G
#SBATCH -p bigmem

module load R/3.6.1-foss-2018b-X11-20180604

Rscript ./estHBD_fixedMu_100trees_all6_forHPC.R
