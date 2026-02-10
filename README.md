# Nonparametric estimation of trawl processes

This repository contains R code implementing the methods described in:

**O. Sauri and A. E. D. Veraart (2026).** *Nonparametric estimation of trawl processes: Theory and Applications.* arXiv:2209.05894. [https://arxiv.org/abs/2209.05894](https://arxiv.org/abs/2209.05894)

## Repository Structure

```
.
├── simulation/          # Section 5.1: Misspecification testing
│   ├── Step1_Simulation.R
│   ├── Step2_Analysis.R
│   └── README.md
│
└── empirical-study/     # Section 5.2: Spread data forecasting
    ├── EmpiricalStudy-Final.R
    ├── AdditionalRFunctions.R
    └── README.md
```

## Quick Start

### Installation
```r
# Core packages required across both projects
install.packages(c("ambit", "trawl", "ggplot2", "forecast", 
                   "dplyr", "tidyr", "reshape2"))
```

### Running the code

**Simulation study (Section 5.1):**
```r
source("simulation/Step1_Simulation.R")
source("simulation/Step2_Analysis.R")
```

**Empirical study (Section 5.2):**
```r
source("empirical-study/EmpiricalStudy-Final.R")
```

## Key methods implemented

### Nonparametric estimation
- Trawl function estimation from discrete observations
- Asymptotic confidence intervals
- Hybrid estimator combining nonparametric and parametric fits

### Misspecification testing
- Data-driven breakpoint detection via hypothesis testing
- Bonferroni-corrected multiple testing procedures
- Comparison with parametric exponential and long-memory trawls

### Forecasting applications
- Four forecasting methods: nonparametric, ACF-based, exponential, and LM trawl
- Diebold-Mariano tests with Benjamini-Hochberg adjustment
- MSE ratio analysis across multiple forecast horizons

## Data requirements

**Simulation:** Self-contained (generates synthetic data)

**Empirical study:** Requires high-frequency financial data in format:
- Files: `{TICKER}_5s_data.txt` (semicolon-separated)
- Structure: Rows = trading days, Columns = 5-second intervals
- Default tickers: A, DFS, WAT, WM

## Output

Both projects generate publication-ready plots (EPS format) and detailed statistical summaries. See individual README files for complete output descriptions.

## Citation

If you use this code, please cite:

```bibtex
@article{sauriveraart2026,
  title={Nonparametric estimation of trawl processes: Theory and Applications},
  author={Sauri, Orimar and Veraart, Almut E. D.},
  journal={arXiv preprint arXiv:2209.05894},
  year={2026}
}
```

## Acknowledgments

Claude.ai (Sonnet 4.5) was used to structure the code and repositories, 
improve the documentation, make the code more efficient, improve the graphics
and provide first drafts of the README files.

