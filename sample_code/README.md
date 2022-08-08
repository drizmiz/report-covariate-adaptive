# R functions for model-robust inference for clinical trials using stratified randomization and covariate adjustment

## Overview

This repo stores R functions that can caculate point estimate and variance of the population average treatment effect under stratified permuated block randomization (or the biased-coin covariate-adaptive design) with adjustment for baseline variables. 

The estimators covered in this repo are the ANCOVA estimator for continuous outcomes, standardized logistic regression estimator for binary outcomes, the doubly-robust weighted-least-squared (DR-WLS) estimator for continuous outcomes with missing outcomes, the Kaplan-Meier esitmator of survival functions for time-to-event outcomes, and the augmented inverse probability weighted (AIPW) estiamtor of the restricted mean survival time for time-to-event outcomes.

For details, please see the paper [here](https://doi.org/10.1080/01621459.2021.1981338).

## Functions

- __R/ICAD.R__ Inference under Covariate-Adaptive Design (ICAD) with continuous or binary outcomes using the ANCOVA estimator, standardized logistic regression estimator, and the DR-WLS estimator depending on the outcome time and where missing outcomes present.

- __R/ICAD-km.R__  Inference of survival functions under Covariate-Adaptive Design with time-to-event outcomes using the Kaplan-Meier estimator.

- __R/ICAD-AIPW.R__ Inference of the restricted mean survival time under Covariate-Adaptive Design (ICAD) with time-to-event outcomes using the AIPW estimator.

## Code of data application

The R scripts under the "data-analysis" folder contain the code of data application in the [paper](https://arxiv.org/abs/1910.13954).
