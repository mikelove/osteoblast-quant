#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1200
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load python
module load samtools
snakemake -j 9 --latency-wait 30 --cluster "sbatch -n 12 --mem=15000 --time=360"
