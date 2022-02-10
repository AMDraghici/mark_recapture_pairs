#!/bin/bash
#SBATCH --array=1-100
#SBATCH --time=84:00:00
#SBATCH --mem=16G
#SBATCH --error=Logs/rep_%A_%a_long.log
#SBATCH --output=Logs/rep_%A_%a_long.log

module load jags
module load r/4.1.2

Rscript --vanilla run_simulation.R $SLURM_ARRAY_TASK_ID param_grid.rds &> Logs/long_run_${SLURM_ARRAY_TASK_ID}.Rout
