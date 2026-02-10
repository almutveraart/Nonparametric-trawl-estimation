# Empirical study of trawl processes

This repository contains R code for conducting
the empirical study in Section 5.2 of the paper, 
comparing different forecasting methods for high-frequency financial data.

## Overview

The code implements and compares four forecasting approaches:
1. Nonparametric trawl estimation
2. ACF-based slice estimation
3. Exponential trawl fit
4. Long Memory (LM) trawl fit

## Requirements

### R packages
The following R packages are required:
- `ggplot2` - for visualization
- `forecast` - for forecasting methods
- `reshape2` - for data reshaping
- `tidyr` - for data tidying
- `dplyr` - for data manipulation
- `ambit` - for trawl process functions
- `trawl` - for parametric trawl fitting

### Installation
```r
# Install required packages
required_packages <- c("ggplot2", "forecast", "reshape2", "tidyr", "dplyr", "ambit", "trawl")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

## Input data format

The code expects data files in the following format:
- File naming: `{TICKER}_5s_data.txt` (e.g., `A_5s_data.txt`)
- Format: Text file with semicolon-separated values (`;`)
- Structure: Matrix where each row represents a trading day and columns represent 5-second intervals

Default tickers included in the analysis:
- A (Agilent Technologies)
- DFS (Discover Financial Services)
- WAT (Waters Corporation)
- WM (Waste Management)

## Usage

### Running the full analysis
```r
# Source the main script
source("EmpiricalStudy-Final.R")
```

### Testing mode
By default, the code runs the full analysis using all available days of data. For faster testing with a subset of data:
1. Open `EmpiricalStudy-Final.R`
2. Uncomment line 26: `used_days <- 10`
3. Comment out line 28: `# used_days <- ndays`

## Output files

The script generates the following outputs for each ticker:

### Plots (EPS format)
- `{TICKER}_Ratio_Naive.eps` - MSE ratio vs naive forecast
- `{TICKER}_Ratio_Exp.eps` - MSE ratio vs exponential trawl
- `{TICKER}_Ratio_LM.eps` - MSE ratio vs LM trawl
- `{TICKER}_Ratio_acf.eps` - MSE ratio vs ACF-based method
- `{TICKER}_DM.eps` - Diebold-Mariano test p-values
- `{TICKER}_AdjPvalues.eps` - Histogram of adjusted p-values
- `{TICKER}_TrawlFct_Boxplot.eps` - Trawl function estimates
- `{TICKER}_WeightComp_Intersec.eps` - Weight comparison for set intersection
- `{TICKER}_WeightComp_Setdiff.eps` - Weight comparison for set difference

### Text output
- `output.txt` - Summary statistics and percentages for each ticker

## Files in this repository

- `EmpiricalStudy-Final.R` - Main analysis script
- `AdditionalRFunctions.R` - Plotting functions
- `README.md` - This file

## Parameters

Key parameters (defined in `EmpiricalStudy-Final.R`):
- `my_Delta = 5/60` - Sample grid width (5 seconds = 5/60 minutes)
- `h_vector = 1:10` - Forecast horizons
- `is_length = floor(80/100*n)` - In-sample length (80% of data)

## Statistical tests

The code performs Diebold-Mariano tests to compare forecast accuracy between methods:
- Uses Benjamini-Hochberg (BH) adjustment for multiple testing
- Reports percentage of significant results at α = 0.05 level


## Notes

- The code uses the `ambit` package for core trawl process functions (including `my_mse()`)
- All messages from `melt()` operations are suppressed for cleaner output
- The script automatically creates an output file and handles cleanup via `on.exit()`
