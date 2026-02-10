# Trawl process: Misspecification testing

This repository contains R code for the simulation and analysis of trawl
processes described in Section 5.1 of the main article. 


## Overview

The code demonstrates:
- Simulation of Poisson and Gaussian trawl processes with piecewise exponential trawl functions
- Nonparametric trawl function estimation
- Data-driven breakpoint detection via hypothesis testing
- Hybrid estimation combining nonparametric and parametric approaches
- Comprehensive hypothesis testing with Bonferroni correction

## Requirements

Required R packages:
```r
install.packages(c("ambit", "ggplot2", "latex2exp", "reshape2", "boot", 
                   "tidyr", "dplyr", "minpack.lm", "forecast", "gridExtra", "grid"))
```

## Usage

### Step 1: Simulate data
```r
source("Step1_Simulation.R")
```

This creates:
- `output/data/path_poisson.csv` - Simulated Poisson trawl process
- `output/data/path_gaussian.csv` - Simulated Gaussian trawl process
- Supporting files with trawl functions and theoretical values

### Step 2: Run analysis
```r
source("Step2_Analysis.R")
```

This performs comprehensive analysis and creates plots in:
- `output/plots/poisson/` - All Poisson analysis plots
- `output/plots/gaussian/` - All Gaussian analysis plots

## Key features

### Hybrid estimator
The hybrid estimator combines nonparametric estimates for short lags with a parametric exponential tail:

- Uses nonparametric estimates up to detected breakpoint
- Fits exponential decay after breakpoint via log-linear regression
- Ensures perfect continuity at the breakpoint
- Robust to noisy tail estimates

### Data-driven breakpoint detection
- Identifies structural breaks via hypothesis test statistics
- Finds lag with maximum |T(t)|
- Comprehensive testing at lags 0-60 with Bonferroni correction

## Output Files

### Data Files (output/data/)
- Simulated trawl process paths (CSV)
- Hypothesis test results (CSV)
- Forecast accuracy comparisons (CSV)

### Plots (output/plots/{poisson,gaussian}/)
- Simulated paths
- ACF with theoretical overlay
- Nonparametric trawl estimates with confidence intervals
- Parametric fits comparison
- Hybrid estimation plots
- Breakpoint identification diagnostics
- Hypothesis test results (test statistics, p-values, differences)
- Forecast accuracy comparison

All plots saved in both PNG and EPS formats.


