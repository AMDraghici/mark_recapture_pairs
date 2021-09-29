#!/bin/bash
#SBATCH --array=10-12
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --error=Logs/rep_%A_%a_long.log
#SBATCH --output=Logs/rep_%A_%a_long.log

module load r

Rscript --vanilla run_test.R $SLURM_ARRAY_TASK_ID param_grid.rds &> Logs/long_run_${SLURM_ARRAY_TASK_ID}.Rout