# BNSM

This repository contains code for the paper `Bayesian Resolution of Discrepant
Self-Reported Network Ties`.

[R_function](https://github.com/donga0223/BNSM/tree/main/R_function) contains code for simulation and application chapter 
 - `thetaprior_function_forpaper.R` It contains all the functions to generate data and estimate the model.
 - `BNSM_summary_function.R` This code is needed to generate a summary plot for the application.
 - `Butts_PPC.R` This code is for conducting a posterior predictive check when using Butts method
 - `Butts_PPC_plots.R` To perform the posterior predictive check after fitting the real data to the Butts method..
 
 - [simulation](https://github.com/donga0223/BNSM/tree/main/R_function/simulation) contains code for simulation study in the simulation section
    - `BNSM_simulation_cluster.R` fits a BNSM model to simulated data for a single simulation replicate
    - `BNSM_simulation_cluster_sh.R` creates cluster jobs to fit a BNSM model; each job calls `BNSM_simulation_cluster.R`
    - `BNSM_simulation_missing_cluster.R` fits a BNSM model to simulated data with missingness for a single simulation replicate
    - `BNSM_simulation_missing_cluster_sh.R` creates cluster jobs to fit a BNSM model; each job calls `BNSM_simulation_missing_cluster.R`
    - Run `BNSM_simulation_summary.R` to produce summary plot.
    
 - [application](https://github.com/donga0223/BNSM/tree/main/R_function/application) contains code for fitting the model to [real_school_data](https://github.com/donga0223/BNSM/tree/main/R_function/application/real_school_data) and conducting data matched simulation studies in the application section
    - [real_school_data](https://github.com/donga0223/BNSM/tree/main/R_function/application/real_school_data) contains code for fitting the model to real data 
      - `BNSM_application.R` fits a BNSM model to real data
      - `BNSM_application_sh.R` creates cluster jobs to fit a BNSM model; each job calls `BNSM_application.R`
      - `BNSM_application_summary.R` to produce summary plot 
    - [data_matched_sim](https://github.com/donga0223/BNSM/tree/main/R_function/application/data_matched_sim) contains code for fitting the model to data matched simulation 
      - `BNSM_data_matched_sim.R` fits a BNSM model to data matched simulation
      - `BNSM_data_matched_sim_sh.R` creates cluster jobs to fit a BNSM model; each job calls `BNSM_data_matched_sim.R`
      - `BNSM_data_matched_sim_summary.R` to produce a summary plot, it produces two plots: One for summarizing and comparing three models, and the other for comparing butts model with our proposed model. 
 