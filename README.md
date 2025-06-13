# Modelling Dependent fates of long-term mates in mark-recapture studies using conditional data methods

The `R` and `C++` code used to conduct and analyze the simulation study from our article "Modelling Dependent fates of long-term mates in mark-recapture studies using conditional data methods" can be found in this repository. For those unfamiliar with Github or Git, a Zip folder containing the files in this repository can be downloaded to your local machine by pressing the big green button labelled Code at the top of the page.

A link to the article associated to this work can be found at the bottom of this document. 

## Details on R Code

### Scripts

Contains all the custom functions `R` written by the authors of this manuscript. 

- `fn_correlation_estimators.R`: contains the novel methods outlined in the article. Namely the functions here can be used to calculate the recapture and survival correlation estimate between mated pairs and to execute the bootstrap algorithm to estimate the standard deviation of the correlation estimates. 

- `fn_generic.R`: general purpose functions for running the simulation studies, computing existing mark-recapture metrics, and data processing. 

### Src

-  `generate_pair_data.cpp`: `C++` code used to generate mark-recapture data in which mated pairs are correlated, intended to be used in the simulation studies discussed in the manuscript. `C++` was utilized here because executing the bootstrap iterations for each simulation replicate can lead to several hundreds of thousands of iterations for just a few hundred scenarios. The speed-up here reduced computational overhead drastically. 

### Simulation

- `run_simulation.R`: code used to execute simulation study. Sources in methods from Scripts and Src. Code to build a scenario grid of parameters for each replicate exists in `fn_generic.R` and can be easily updated by the end user wishing to experiment with different parameter settings. 

- `run_sim_slurm.sh`: bash script to execute multiple simulation scenarios at once in parallel on an HPC cluster. Sharcnet was used for this work: https://www.sharcnet.ca/my/front/

- `investigate_simulation.R`: Code used to process the results of multiple simulations into clean dataframes. 

Article Link: TBD
Reference: TBD
