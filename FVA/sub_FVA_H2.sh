#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=2:00:00
#SBATCH --account=cdm8

echo " "
echo "JOB started on $(hostname -s) at $(date)"
module load gams
gams FVA_H2.gms
echo "job ended at $(date)"
