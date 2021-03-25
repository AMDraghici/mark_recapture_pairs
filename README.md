# Pair-Swap Joint Fates CJS Model

This repository is where I am tracking my work for this project. The scripts in the repo are briefly documented below. 

## 00_fn_sim_pair_data.R

Functions for generating data from the mark-recapture pair-swap joint fates CJS model. Broken up into several modular functions which allows quick trackable changes.

This script contains an outdated method for producing datasets. Refer to 00_fn_sim_pair_data_rework.R (This will be deleted soon)

## 00_fn_sim_pair_data_rework.R

Current method for generating model data. Also includes function to convert generated data into standard CJS data for model comparison. 

## 01_fn_process_live_data.R

Script containing functions for preparing the harlequin duck data for modelling with the pair-swap joint fates model. Mostly just data munging. 

## 02_fn_model_code.R

Higher level functions for users. Only contains two atm - 

- One for generating several datasets from 00 using parallel and
- One for running JAGS in parallel (this is actually just a general purpose function for the most part).

## 10_mod_pair_swap_time.R

Current iteration of the pair-swap joint fates model with time-variable parameters (difficult to get estimates for these -> needs lots of data).

TODO:
- The pair-swap mechanism needs work - it is incorrect in its current state. 
- The way the model handles the first time point is also really sloppy could use some work. Nimble might help since I can build out functions which can at least organize the code better. 
- Method for dealing with the likelihood calculation is correct (I think) but really inefficent wrt to memory/speed once population starts getting larger. There's probably a bunch of calculations that can be streamlined. 
- Priors are way too flat and will slow things down quite a bit (low priority thing to fix right now). 

## 11_mod_pair_swap_notime.R

See above but without time variable parameters (only one survival rate, etc). 
FYI: Recruitment is still time-variable and the model is fine with estimating this as it stands. 
List of TODOs is the same as 10. 

## 11_cjs_mod_nimble_WIP.R

Attempt at NIMBLIZING model 11. Not currently working. Oddly enough the code compiles in C like any other nimble model. When I try to run it Rstudio crashes with message:

================== C Stack trace ========================================================

(No symbol) [0x0x3b9e3b2]
(No symbol) [0x0x2185a2f589a0]
RtlAllocateHeap [0x0x7ffdb06fb86b+3787]
(No symbol) [0x0x64]
(No symbol) [0x0x4a577f9740]
(No symbol [0x0x1]

and then the R session encounters a fatal error. The Eigen methods in C did show some warnings during compilation about multiplication but they were not clear. Might have be a simple matter of how some variable/operation is defined in the syntax. 

Appeal of using NIMBLE: 

- Can write own MCMC sampler 
- Can write own LL function (might help with efficency issues)
- Can break down sampler into own functions and custom distributions (of which I have many) which will clean the syntax

## 14_cjs_mod.R

Vanilla CJS model coded up in JAGS. For comparison studies. 

## 14_cjs_mod_nimble.R

Nimblized code for vanilla CJS model. Works well but code is messy. Mostly just a learning exercise at the time being. 

## 20_call_sim_study.R

Code for sourcing in the functions from the above scripts. Contains some notes to self at the bottom. Likely to be replaced by this document. 

