## Overview
Implementation of a Bayesian Hyperbolic Latent Space Clustering Model (HLSCM).

This repository includes:
- MCMC estimation of latent positions and model parameters
- Hyperbolic latent position clustering
- Case-control approximate likelihood for faster inference
- Hyperbolic Procrustes alignment
- Hyperbolic Barycenter computation
- Poincaré disk visualisation for posterior summaries such as barycenters and posterior clouds
- Posterior predictive checks and diagnostics

## Dependencies
- igraph
- MASS
- mvtnorm
- clue
- hydra

## Status
Research Code under development. Further work planned for ease of use, particularly in hyperparameter selection.
