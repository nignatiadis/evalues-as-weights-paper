# evalues-as-weights-paper

Reproduction code for the paper:

 > E-values as unnormalized weights in multiple testing
 
 
This repository contains the following:

## `EValueWeightsPaper` R package

A helper package with functions implementing the different methods evaluated as well as functions generating synthetic datasets for simulation experiments. It may be installed as follows:

```{r}
devtools::install_github("nignatiadis/evalues-as-weights-paper")
```

## Simulation scripts, and figures

Code to reproduce the simulations is available in the folder `simulations`. The simulation scripts depend upon the `EValueWeightsPaper` package.

## Reproduction of data analysis

See the `data_analysis` folder. The main analysis has been implemented in `bottomly.qmd` as a [quarto](notebook).
