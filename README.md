# Bayesian Inference for Partial Orders from Rank Data with Ties

This repository contains the code, simulation results, and figures for my
summer project on Bayesian inference for partial orders from ranking data
with ties.

The project develops an observation model for tied ranking data and explores
Bayesian inference for an underlying partial order using MCMC. The methods
are investigated through simulation experiments and an application to
primate intelligence ranking data.

## Repository structure

### `Codebase/`
Contains the functions used throughout the project, including the
implementation of the MCMC algorithms, likelihood and prior models,
data-generation procedures, and other utilities used in the simulations.

### `Code and simulation results/`
Contains the scripts used to run the simulation experiments together with
the corresponding simulation outputs. These scripts are the ones used to
produce the results reported in the project.

### `Plots/`
Contains the figures generated from the simulation experiments and used in
the project report.

## Project report

The LaTeX source for the project report is available on Overleaf:

https://www.overleaf.com/read/qvtxsfcgsvkv#9fed78

The report describes the statistical models, theoretical results, MCMC
methodology, simulation experiments, and application in detail.

## Reproducing the results

The scripts in `Code and simulation results/` contain the exact code used
to run the simulations and generate the corresponding plots.

The supporting functions required by these scripts can be found in
`Codebase/`.

Note that some of the MCMC experiments use a large number of iterations
and may therefore take a substantial amount of time to reproduce.

## Contents

In particular, the repository contains code for:

- generating partial orders and ranking data;
- uniform and frontier-softmax observation models;
- ranking models with and without noise;
- generating tied rankings (bucket orders);
- MCMC inference for latent partial orders;
- simulation studies and posterior diagnostics;
- consensus partial-order estimation; and
- analysis of the primate ranking data.

## Author

Ho Man Yiu

Summer Project, University of Oxford
