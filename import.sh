#!/bin/bash
#
#SBATCH --job-name=import_allelic
#SBATCH --time=60
#SBATCH --mem=15000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

module load r/4.2.1
R CMD BATCH --no-save --no-restore import_allelic.R
R CMD BATCH --no-save --no-restore export_allelic.R
