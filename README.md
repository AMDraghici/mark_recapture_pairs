# Pair-Swap Joint Fates CJS Model

This repository is where I am tracking my work for this project. The scripts in the repo are briefly documented below. 

## 00_fn_sim_pair_data.R

Functions for generating data from the mark-recapture pair-swap joint fates CJS model. Broken up into several modular functions which allows quick traceable changes.  Also includes function to convert generated data into standard CJS data for model comparison. 

## 02_fn_model_code.R

Higher level functions for users. 

- One for generating several datasets from 00 using parallel and
- One for running JAGS in parallel (this is actually just a general purpose function for the most part).
- Some functions for plotting results

## 03_fn_process_hduck_data.R

Script containing functions for preparing the harlequin duck data for modeling with the pair-swap joint fates model. Mostly just data munging. 

## 10_mod_pair_swap_notime.R

Current iteration of the pair-swap joint fates model without time-variable parameters (difficult to get estimates for these -> needs lots of data).

TODO:
- The pair-swap mechanism needs work - it is incorrect in its current state. 
- Priors are way too flat and will slow things down quite a bit (low priority thing to fix right now). 

## 11_cjs_mod_standard.R

Vanilla CJS model coded up in JAGS. For comparison studies. 

## 12_cjs_mod_nimble.R

Nimblized code for vanilla CJS model. Works well but code is messy. Mostly just a learning exercise at the time being. 

Appeal of using NIMBLE: 

- Can write own MCMC sampler 
- Can write own LL function (might help with efficency issues)
- Can break down sampler into own functions and custom distributions (of which I have many) which will clean the syntax

## 20_call_sim_study.R

Code for sourcing in the functions from the above scripts. Contains some notes to self at the bottom.
