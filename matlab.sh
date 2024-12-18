#!/bin/bash -l
#SBATCH --array=1-7
#SBATCH --nodes 1 
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 72
#SBATCH --time 72:00:00
#SBATCH --qos=serial

#SBATCH --partition=bigmem
#SBATCH --mail-user=ZHE.ZHANG@epfl.ch
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --error=matlab_%A_%a.err
#SBATCH --output=matlab_%A_%a.out
set -e

ulimit -n 131072
echo STARTING AT `date`
echo "Slurm Job Id SLURM_ARRAY_JOB_ID is ${SLURM_ARRAY_JOB_ID}"

module purge
module load matlab/R2024a
module list

matlab -nosplash -nodisplay -nodesktop -r "Main($SLURM_ARRAY_TASK_ID)", exit

$ ln -s /scratch/zhang2/ServiceHost /home/zhang2/.MathWorks/ServiceHost
echo FINISHED at `date`
