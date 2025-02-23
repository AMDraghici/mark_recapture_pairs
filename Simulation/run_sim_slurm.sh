#!/bin/bash
#SBATCH --array=1-648
#SBATCH --time=90:00:00
#SBATCH --mem=5G
#SBATCH --error=Logs/rep_%A_%a_long.log
#SBATCH --output=Logs/rep_%A_%a_long.log

module load r/4.1.2

Rscript --vanilla run_simulation.R $SLURM_ARRAY_TASK_ID &> Logs/long_run_${SLURM_ARRAY_TASK_ID}.Rout
