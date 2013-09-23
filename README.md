prediction_legacy
=================
This git repo contains codes for computing Bayesian predictions using MCMC results from MCMC_legacy folder. 
The codes have 2 versions: one runs the prediction procedure locally; and the other one is for the Amazon cloud map-reduce
procedure.

local version:
1. step 1: invcor_generate_local.R: generate cholesky decomposition matrices for MCMC covariance matrices;
2. step 2: predict_MCMC_meansd_chol_local.R: given MCMC results, prediction grids with covariates attached to make predictions locally.

cloud version:
divide the prediction data into small files, so that each file can be used independently
1. step 1: 1.invcor_generate_foreachtif_SOC.R:  get the index of clusters needed for each small file (i.e. clusters which have data correlated with locations in small files);
2. step 2: 2.invcor_generate.R and 3.seperate_chollist.R: generate cholesky decomposition matrices for MCMC covariance matrices for each small file;
3. step 3: send data into cloud;
4. step 4: run predict_MCMC_meansd_onenode_chol_SOC.R

Note: the data names need to be changed for predicting different soil property.

Author: Jiehua Chen <jc3288@columbia.edu>