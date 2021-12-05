#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1200
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load python
module load wasp
module load samtools
snakemake -j 2 --latency-wait 30 --cluster "sbatch -n 1 --mem=15000 --time=1200"
