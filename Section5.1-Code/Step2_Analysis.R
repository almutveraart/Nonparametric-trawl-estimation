################################################################################
# PART 2: Trawl Process Analysis with Hypothesis Testing
#
# This script analyses previously simulated trawl processes with comprehensive
# hypothesis testing at all lags from 1-60.
#
# Run Step1_Simulation.R first to generate the data.
################################################################################

# Load required libraries
library(ambit)
library(ggplot2)
library(latex2exp)
library(reshape2)
library(boot)
library(tidyr)
library(dplyr)
library(minpack.lm)
library(forecast)
library(gridExtra)  # For multi-panel plots
library(grid)       # For textGrob function

# Source all analysis functions
source("Step2_Functions.R")

# Suppress harmless warnings during analysis
# (log10 of p-values near zero, and EPS transparency warnings)
options(warn = -1)
################################################################################
# 1. LOAD SIMULATED DATA
################################################################################

cat("================================================================================\n")
cat("LOADING SIMULATED DATA\n")
cat("================================================================================\n")

# Check if data files exist
if (!file.exists("output/data/path_poisson.csv")) {
  stop("ERROR: Simulated data not found. Please run Step1_Simulation.R first.")
}

# Load paths
path_Poi <- read.csv("output/data/path_poisson.csv")$value
path_Gauss <- read.csv("output/data/path_gaussian.csv")$value

# Load simulation parameters
simulation_info <- readRDS("output/data/simulation_info.rds")
n <- simulation_info$n
Delta <- simulation_info$Delta

# Load trawl functions
a <- readRDS("output/data/trawl_function_a.rds")
a2 <- readRDS("output/data/trawl_function_a2.rds")

# Load theoretical values
theoretical_values <- readRDS("output/data/theoretical_values.rds")
rho <- readRDS("output/data/rho_function.rds")

cat("✓ Data loaded successfully\n")
cat("  Poisson path: ", length(path_Poi), " observations\n", sep = "")
cat("  Gaussian path: ", length(path_Gauss), " observations\n", sep = "")
cat("  n = ", n, ", Delta = ", Delta, "\n", sep = "")
cat("\n")

################################################################################
# 2. CREATE OUTPUT DIRECTORIES
################################################################################

if (!dir.exists("output/plots")) {
  dir.create("output/plots")
}

################################################################################
# 3. PLOT THEORETICAL TRAWL FUNCTIONS
################################################################################

cat("Plotting theoretical trawl functions and ACF...\n")

# Plot 1: Piecewise exponential trawl function
# PNG version
png("output/plots/theoretical_trawl_piecewise.png", width = 800, height = 800, res = 150)
par(cex.lab = 1.5, cex.axis = 1.3, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$a_vals, type = 'l', lwd = 2, col = 'blue',
     main = "", xlab = "x", ylab = "a(x)")
abline(v = 2, lty = 2, col = "red")
grid(col = "gray80", lty = 1)
dev.off()

# EPS version
setEPS()
postscript("output/plots/theoretical_trawl_piecewise.eps", width = 4, height = 4)
par(cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$a_vals, type = 'l', lwd = 2, col = 'blue',
     main = "", xlab = "x", ylab = "a(x)")
abline(v = 2, lty = 2, col = "red")
grid(col = "gray80", lty = 1)
dev.off()

# Plot 2: Theoretical ACF for piecewise exponential
# PNG version
png("output/plots/theoretical_acf_piecewise.png", width = 800, height = 800, res = 150)
par(cex.lab = 1.5, cex.axis = 1.3, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$rho_vals, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = expression(rho(x)),
     main = "")
grid(col = "gray80", lty = 1)
dev.off()

# EPS version
setEPS()
postscript("output/plots/theoretical_acf_piecewise.eps", width = 4, height = 4)
par(cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$rho_vals, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = expression(rho(x)),
     main = "")
grid(col = "gray80", lty = 1)
dev.off()

# Plot 3: Simple exponential trawl function
# PNG version
png("output/plots/theoretical_trawl_exponential.png", width = 800, height = 800, res = 150)
par(cex.lab = 1.5, cex.axis = 1.3, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$a2_vals, type = 'l', lwd = 2, col = 'blue',
     main = "", xlab = "x", ylab = "a(x)")
grid(col = "gray80", lty = 1)
dev.off()

# EPS version
setEPS()
postscript("output/plots/theoretical_trawl_exponential.eps", width = 4, height = 4)
par(cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$a2_vals, type = 'l', lwd = 2, col = 'blue',
     main = "", xlab = "x", ylab = "a(x)")
grid(col = "gray80", lty = 1)
dev.off()

# Plot 4: Theoretical ACF for simple exponential
# PNG version
png("output/plots/theoretical_acf_exponential.png", width = 800, height = 800, res = 150)
par(cex.lab = 1.5, cex.axis = 1.3, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$rho2_vals, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = expression(rho(x)),
     main = "")
grid(col = "gray80", lty = 1)
dev.off()

# EPS version
setEPS()
postscript("output/plots/theoretical_acf_exponential.eps", width = 4, height = 4)
par(cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5, mar = c(5, 5, 2, 2), las = 1)
plot(theoretical_values$x_vals, theoretical_values$rho2_vals, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = expression(rho(x)),
     main = "")
grid(col = "gray80", lty = 1)
dev.off()

################################################################################
# 4. RUN ANALYSIS FOR BOTH POISSON AND GAUSSIAN
################################################################################

cat("\n\n")
cat("================================================================================\n")
cat("STARTING COMPREHENSIVE TRAWL PROCESS ANALYSIS (ALL LAGS 1-60)\n")
cat("================================================================================\n")
cat("NOTE: This version tests at ALL lags 1-60 for comprehensive assessment.\n")
cat("      Bonferroni correction will be more conservative (α/61 for 61 tests).\n")
cat("\n")

# Analyze Poisson process
results_poisson <- analyse_trawl_process(path_Poi, "Poisson", Delta, n)

# Analyze Gaussian process
results_gaussian <- analyse_trawl_process(path_Gauss, "Gaussian", Delta, n)

################
# Restore default warning level
options(warn = 0)

################################################################################
# END OF ANALYSIS SCRIPT
################################################################################
