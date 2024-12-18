#!/bin/bash -l
#SBATCH --chdir /home/zhang2/MF_CC_2024/Array_test/
#SBATCH --nodes 4 
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 72
#SBATCH --time 72:00:00
#SBATCH --qos=parallel
#SBATCH --array=1-4

#SBATCH --partition=bigmem
#SBATCH --mail-user=ZHE.ZHANG@epfl.ch
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --error=matlab_%A_%a.err
#SBATCH --output=matlab_%A_%a.out
ulimit -u 16384
echo STARTING AT `date`

module purge 
module load matlab

srun matlab -nosplash -nodisplay -nodesktop -r "Main_test($SLURM_ARRAY_TASK_ID)"

echo "I am array task number" $SLURM_ARRAY_TASK_ID
echo FINISHED at `date`
