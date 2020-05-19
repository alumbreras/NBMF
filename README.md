# Bernoulli Non-Negative Matrix Factorization

Code for our paper "Bayesian Mean-parameterized Nonnegative Binary Matrix Factorization"

## Code structure

* R: algorithms in R. The main file is BernoulliNMF.R, that contains the function BernoulliNMF, that performs the inference of the latent matrices W and H. The inference method is chosen through the parameters. Both MCMC and VB methods are mostly writenn in Rcpp ('src' folder).

* src: implementation of MCMC and VB algorithms in Rcpp. These methods are called from the ones in the R folder.

* man: automatically generated folder with documentation of the R functions.

* data: datasets ready to be used.

* experiments: scripts to run the experiments.

## Experiments

