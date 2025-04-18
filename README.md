# Propensity Score Matching under Multiple Treatments in Causal Inference

## Overview

This repository contains all the codes and analyzes of propensity score matching in multi-treatments settings. This is also the final project of STAT840. The main methods implemented in this repo is the Vector Matching VM and two of its extensions called VM2 and VM_MD. The main finding is that, although Vector Matching has satisfied performance in balancing the covariates and estimating the true treatment effects, it is computational intensive especially when the number of covariates and treatment increases.

All the data are simulated and stored in the folder `simulation` which also contains the codes about how to generate the simulation data. The output paper is `psm_multi-treatment.pdf`, which is knitted by `psm_multi-treatment.rmd`.

## File Structure

The repo is structured as:

- `literatures` contains all the recent researches cited by this paper.
- `simulation` contains all the simulated data with corresponding codes.
- `psm_multi-treatment.pdf` is the paper.
- `psm_multi-treatment.rmd` is the raw rmd file of the paper.
- The rest `.R` files are all supplementry codes.


 
