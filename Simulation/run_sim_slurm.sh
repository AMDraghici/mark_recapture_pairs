#!/bin/bash
#SBATCH --array=483,556,667,682,694,756,812,815,960,961,962,963,964,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,
#SBATCH --time=30:00:00
#SBATCH --mem=4G
#SBATCH --error=Logs/rep_%A_%a_long.log
#SBATCH --output=Logs/rep_%A_%a_long.log

module load r/4.1.2

Rscript --vanilla run_simulation.R $SLURM_ARRAY_TASK_ID &> Logs/long_run_${SLURM_ARRAY_TASK_ID}.Rout
