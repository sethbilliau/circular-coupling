# Circularly Coupled MCMC (WIP)
Illustrating Circularly Coupled MCMC from Radford Neal's “Circularly-Coupled Markov Chain Sampling.” 

Authors: Seth Billiau, Christine Cai, Michael Yin with invaluable contributions from Prof. Pierre Jacob and his course on Couplings and Monte Carlo

### Introduction
Circular coupling is a method of coupling that circularly reuses common random numbers. In Neal's words, circular coupling is in fact a "refinement" of the common random numbers coupling described in a paper by Johnson (1995). Circular coupling is designed specifically to avoid discarding burn-in times, which can introduce extra bias: instead of running the Markov Chain process and discarding the first few terms, we can just take a "wrapped around" chain because it will be approximately distributed as the Markov chain's stationary distribution. Neal's proposed diagnostic for convergence, auxiliary chains, also allow for more rigorous checks of convergence speed than just using one sequence of CRN like in Johnson, all at the same starting point. However, circular coupling has its shortcomings, and quite a few at that. Most notably, it does not necessarily guarantee convergence to the stationary distribution, and can only guarantee a rough upper bound of total variation distance from the stationary distribution. This upper bound is also near-impossible to verify computationally in the first place, which is why one must approximate with diagnostic methods (that are sometimes quite computationally expensive). The coupling is also only successful if the wrap around chain converges by some fixed time N, but it is difficult to determine if it will even converge within a fixed time N -- and it is also difficult to determine how much to increase N if it does not coalesce, so many bounds can seem somewhat arbitrary.

`report.pdf` contains our full project report. 

### Repository Description

- `randomGridMetropolis.Rmd` is a self-contained notebook that illustrates coupled random-grid Metrpolis sampling. This code produced the figures and results in section 2.1 of the report. 
- `coupledGibbs.Rmd` is a self-contained notebook that illustrates coupled Gibbs sampling. This code produced the figures and results in section 2.2 of the report. 
- `logisticRegression.R` contains the code used to implement the Bayesian Logistic regression problem in section 3. `helpers.R` contains helper functions for `logisticRegression.R`. `seed29.csv` and `seed29_secondrun.csv` contain the results of the original and wrap around chains for the logistic regression example. 

Original Paper link: https://arxiv.org/abs/1711.04399. 
