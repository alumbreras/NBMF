# Bernoulli Non-Negative Matrix Factorization

Code for our paper "Bayesian Mean-parameterized Nonnegative Binary Matrix Factorization"

## Project structure

* R: algorithms in R. The main file is BernoulliNMF.R, that contains the function BernoulliNMF, that performs the inference of the latent matrices W and H. The inference method is chosen through the parameters. Both MCMC and VB methods are mostly writenn in Rcpp ('src' folder).

* `src`: implementation of MCMC and VB algorithms in Rcpp. These methods are called from the ones in the R folder.

* `man`: automatically generated folder with documentation of the R functions.

* `data`: datasets ready to be used.

* `experiments`: scripts to run the experiments.

## Experiments

Experiments and figures are reproducible with the scripts in the `experiments` folder.

First, in order to get acquainted with the general workflow and functions, we suggest you go trough the file
`xp_paper_basic_example.R`

The other files are:

* `xp_paper_fig2_synthetic`: generates synthetic data from the different generative models (Figure 2).

* `xp_paper_fig5_originals`: plots the datasets used in the experiments (Figure 3).

* `xp_paper_fig6_reconstruct`: infers the latent factors and plot the expectation of V (Figure 6). 

* `xp_paper_fig7-8-9_dictionaries`: infers and plots the expectation of the dictionaries W (Figures 7 to 9). 

* `xp_paper_fig10_predictive_likelihood`: infers the latent factors and posterior of V, and plots the log-likelihood on unobserved parts of the matrix (Figure 10).

* `xp_paper_fig11_VB_convergences`: runs VB algorithms and plots their predictive log-likelihood to compare speed of convergence (Figure 11). 