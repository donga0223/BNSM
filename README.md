# BNSM

This repository contains code for the paper `Bayesian Resolution of Discrepant
Self-Reported Network Ties`.

[R_function](https://github.com/donga0223/BNSM/tree/main/R_function) contains code for simulation and application chapter 
 - [simulation](https://github.com/donga0223/BNSM/tree/main/R_function/simulation) contains code for simulation study in simulation section
    - `thetaprior_function_forpaper.R` contain all function to generate data and model estimation 
    - `BNSM_simulation_cluster.R` fits a BNSM model to simulated data for a single simulation replicate
    - `BNSM_simulation_cluster_sh.R` creates cluster jobs to fit a BNSM model; each job calls `BNSM_simulation_cluster.R`
    - `BNSM_simulation_missing_cluster.R` fits a BNSM model to simulated data with missingness for a single simulation replicate
    - `BNSM_simulation_missing_cluster_sh.R` creates cluster jobs to fit a BNSM model; each job calls `BNSM_simulation_missing_cluster.R`
    - Run `BNSM_simulation_summary.R` to produce summary plot.
 - [data_matched_sim](https://github.com/donga0223/BNSM/tree/main/R_function/data_matched_sim) contains code for data matched simulation in application section
