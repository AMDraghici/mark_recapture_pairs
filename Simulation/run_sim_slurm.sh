#!/bin/bash
#SBATCH --array=4,5,6,7,12,13,14,15,22,23,24,25,31,32,50,64,65,90,107,108,118,132,133,186,187,235,236,237,238,239,240,241,242,243,298,315,316,317,384,394,395,408,409,415,416,417,418,428,429,430,431,434,444,454,455,464,465,466,467,468,488,489,490,491,496,497,498,499,503,504,505,506,512,513,514,515,516,523,524,525,526,527,533,534,535,536,537,545,546,547,548,549,553,554,555,556,557,561,562,563,564,572,592,612,624,634,642
#SBATCH --time=25:00:00
#SBATCH --mem=4G
#SBATCH --error=Logs/rep_%A_%a_long.log
#SBATCH --output=Logs/rep_%A_%a_long.log

module load r/4.1.2

Rscript --vanilla run_simulation.R $SLURM_ARRAY_TASK_ID &> Logs/long_run_${SLURM_ARRAY_TASK_ID}.Rout
