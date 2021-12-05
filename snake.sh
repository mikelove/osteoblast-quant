#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=240
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load python
module load samtools
snakemake -j 1 --latency-wait 30 --cluster "sbatch -n 1 --mem=15000 --time=120"
