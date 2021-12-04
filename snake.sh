#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=120
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load python/3.8.8
module load multiqc
snakemake -j 1 --latency-wait 30 --cluster "sbatch -n 1 --mem=10000 --time=120"
