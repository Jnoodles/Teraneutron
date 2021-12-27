#!/bin/sh
#SBATCH --partition=AllNodes
#SBATCH --job-name=tetraneutron
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=output0-5
ls /local/storage
srun singularity exec -B /local/storage:/local/storage /home/jwzhu/openmp.sif ./command.sh $SLURM_ARRAY_TASK_ID
