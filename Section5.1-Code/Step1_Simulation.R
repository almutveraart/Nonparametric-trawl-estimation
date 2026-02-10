################################################################################
# PART 1: Trawl Process Simulation
#
# This script simulates the Poisson and Gaussian trawl processes
# used in Section 5.1 of the article Sauri and Veraart (2026)
################################################################################

# Load required libraries
library(ambit)

# Set global seed for reproducibility
set.seed(123)

################################################################################
# 1. SIMULATION PARAMETERS
################################################################################

# Set simulation parameters
burnin <- 100                    # Burnin period for simulation
n <- 10000                       # Number of observations
Delta <- 0.1                     # Grid width
n_with_burnin <- n + burnin

# Create output directory for data
if (!dir.exists("output")) {
  dir.create("output")
}
if (!dir.exists("output/data")) {
  dir.create("output/data")
}

################################################################################
# 2. TRAWL FUNCTION DEFINITIONS
################################################################################

# Piecewise exponential trawl function with one breakpoint
a_trawl <- function(x, a0, lambda1, lambda2, c1) {
  a <- numeric(length(x))
  a_c1 <- a0 * exp(-lambda1 * c1)

  a[x < c1] <- a0 * exp(-lambda1 * (x[x < c1]))
  a[x >= c1] <- a_c1 * exp(-lambda2 * (x[x >= c1] - c1))

  return(a)
}

# Simple exponential trawl function
a_trawl_exp <- function(x, a0, lambda1) {
  a <- numeric(length(x))
  a[] <- a0 * exp(-lambda1 * x)
  return(a)
}

# Define the specific trawl functions to use
a <- function(x) {a_trawl(x, a0 = 5, lambda1 = 0.1, lambda2 = 0.5, c1 = 2)}
a2 <- function(x) {a_trawl_exp(x, a0 = 5, lambda1 = 0.3)}

################################################################################
# 3. SIMULATE TRAWL PROCESSES
################################################################################

cat("================================================================================\n")
cat("SIMULATING TRAWL PROCESSES\n")
cat("================================================================================\n")
cat("Parameters:\n")
cat("  n =", n, "\n")
cat("  Delta =", Delta, "\n")
cat("  Burnin =", burnin, "\n")
cat("\n")

# Reset seed for reproducible simulation
set.seed(123)

# --- POISSON TRAWL PROCESS ---
cat("Simulating Poisson process...\n")
distr <- "Poi"
distr_par <- 1 # rate of the Poisson process

path_Poi_with_burnin <- sim_weighted_trawl_gen(n_with_burnin, Delta,
                                               trawlfct_gen = a, distr, distr_par)$path

path_Poi <- path_Poi_with_burnin[(burnin+1):n_with_burnin]
cat("  ✓ Poisson process simulated (length:", length(path_Poi), ")\n")

# --- GAUSSIAN TRAWL PROCESS ---
cat("Simulating Gaussian process...\n")
set.seed(123)  # Reset seed for the Gaussian trawl
distr <- "Gauss"
distr_par <- c(0, 1)  # mean 0, std 1

path_Gauss_with_burnin <- sim_weighted_trawl_gen(n_with_burnin, Delta,
                                                 trawlfct_gen = a, distr, distr_par)$path
path_Gauss <- path_Gauss_with_burnin[(burnin+1):n_with_burnin]

cat("  ✓ Gaussian process simulated (length:", length(path_Gauss), ")\n")

################################################################################
# 4. SAVE SIMULATED DATA
################################################################################

cat("\nSaving simulated data...\n")

# Save paths as CSV
write.csv(data.frame(value = path_Poi),
          file = "output/data/path_poisson.csv",
          row.names = FALSE)
write.csv(data.frame(value = path_Gauss),
          file = "output/data/path_gaussian.csv",
          row.names = FALSE)

# Save simulation parameters and trawl functions
simulation_info <- list(
  n = n,
  Delta = Delta,
  burnin = burnin,
  seed = 123,
  trawl_params = list(
    a0 = 5,
    lambda1 = 0.1,
    lambda2 = 0.5,
    c1 = 2
  ),
  poisson_params = list(
    distr = "Poi",
    distr_par = 1
  ),
  gaussian_params = list(
    distr = "Gauss",
    distr_par = c(0, 1)
  )
)

saveRDS(simulation_info, file = "output/data/simulation_info.rds")

# Save trawl function definitions for analysis
saveRDS(a, file = "output/data/trawl_function_a.rds")
saveRDS(a2, file = "output/data/trawl_function_a2.rds")

cat("  ✓ Saved: output/data/path_poisson.csv\n")
cat("  ✓ Saved: output/data/path_gaussian.csv\n")
cat("  ✓ Saved: output/data/simulation_info.rds\n")
cat("  ✓ Saved: output/data/trawl_function_a.rds\n")
cat("  ✓ Saved: output/data/trawl_function_a2.rds\n")

################################################################################
# 5. COMPUTE AND SAVE THEORETICAL VALUES
################################################################################

cat("\nComputing theoretical values...\n")

# Define GammaInt(x) = ∫_x^∞ a(s) ds numerically
# Note: Named "GammaInt" to avoid conflict with stats::Gamma
GammaInt <- function(x) {
  integrate(a, lower = x, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
}

GammaInt2 <- function(x) {
  integrate(a2, lower = x, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
}

# ACF functions: ρ(x) = Γ(x) / Γ(0)
rho <- function(x) {
  GammaInt(x) / GammaInt(0)
}

rho2 <- function(x) {
  GammaInt2(x) / GammaInt2(0)
}

# Evaluate functions for plotting
x_vals <- seq(0, 10, by = 0.01)
theoretical_values <- list(
  x_vals = x_vals,
  a_vals = a_trawl(x_vals, a0 = 5, lambda1 = 0.1, lambda2 = 0.5, c1 = 2),
  a2_vals = a_trawl_exp(x_vals, a0 = 5, lambda1 = 0.3),
  Gamma_vals = sapply(x_vals, GammaInt),
  Gamma2_vals = sapply(x_vals, GammaInt2),
  rho_vals = sapply(x_vals, rho),
  rho2_vals = sapply(x_vals, rho2)
)

saveRDS(theoretical_values, file = "output/data/theoretical_values.rds")
saveRDS(rho, file = "output/data/rho_function.rds")

cat("  ✓ Saved: output/data/theoretical_values.rds\n")
cat("  ✓ Saved: output/data/rho_function.rds\n")

################################################################################
# 6. SUMMARY
################################################################################

cat("\n================================================================================\n")
cat("SIMULATION COMPLETE\n")
cat("================================================================================\n")
cat("\nSimulated data saved to output/data/:\n")
cat("  - path_poisson.csv (", length(path_Poi), " observations)\n", sep = "")
cat("  - path_gaussian.csv (", length(path_Gauss), " observations)\n", sep = "")
cat("  - simulation_info.rds\n")
cat("  - trawl_function_a.rds\n")
cat("  - trawl_function_a2.rds\n")
cat("  - theoretical_values.rds\n")
cat("  - rho_function.rds\n")
cat("\nYou can now run the analysis script: Step2_Analysis.R\n")
cat("\n")

################################################################################
# END OF SIMULATION SCRIPT
################################################################################
