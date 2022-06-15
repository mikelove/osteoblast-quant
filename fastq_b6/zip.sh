#!/bin/bash
#
#SBATCH --job-name=zip
#SBATCH --time=360
#SBATCH --mem=10000
#SBATCH -n 8
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load pigz
for f in `ls B6*`; do echo $f; pigz $f; done
