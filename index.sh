#!/bin/bash
#
#SBATCH --job-name=index
#SBATCH --time=120
#SBATCH --mem=10000
#SBATCH --ntasks=12
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

/proj/milovelab/bin/salmon-1.5.2_linux_x86_64/bin/salmon index -t /proj/milovelab/anno/Mus_musculus.GRCm38.v102.fa.gz -i /proj/milovelab/anno/Mus_musculus.GRCm38.v102-salmon_1.5.2
