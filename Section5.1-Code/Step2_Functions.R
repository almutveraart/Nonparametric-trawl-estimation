################################################################################
# STEP 2 ANALYSIS - FUNCTIONS FILE
#
# This script contains all functions for trawl process analysis including:
# - Trawl function definitions
# - Hypothesis testing
# - Plotting functions
# - Forecasting comparison
# - Main analysis workflow
################################################################################

################################################################################
# SECTION 0: TRAWL FUNCTION DEFINITIONS
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

################################################################################
# SECTION 1: HYPOTHESIS TESTING FUNCTIONS
################################################################################

#
#Trawl hypothesis testing using functions from the ambit package
perform_trawl_specification_test <- function(est_trawl, a_par_fitted, Delta, n,
                                             path_data, lag_values, c4_est,
                                             alpha = 0.05, parametric_name = "Parametric") {

  # Number of tests
  K <- length(lag_values)

  # Bonferroni correction
  bonf_alpha <- alpha / K
  z_crit <- qnorm(1 - alpha/2)
  z_crit_bonf <- qnorm(1 - bonf_alpha/2)

  # Variance of Levy seed (default 1)
  varlevyseed <- 1

  # Initialize results list
  test_results <- list()

  for (idx in 1:K) {
    lag_value <- lag_values[idx]

    # Array index: est_trawl[1] corresponds to lag 0, so lag L is at index L+1
    array_index <- lag_value + 1

    # Check bounds
    if (array_index > length(est_trawl)) {
      warning(paste("Lag", lag_value, "exceeds available data"))
      next
    }

    # Get estimates
    a_hat <- est_trawl[array_index]
    a_par <- a_par_fitted[array_index]

    # Compute asymptotic variance using user's function
    # asymptotic_variance_est(t, c4, varlevyseed, Delta, avector)
    t <- lag_value * Delta

    if (lag_value == 0) {
      # For lag 0, use RQ estimator
      av <- rq(path_data, Delta)
      std_error <- sqrt(av / (n * Delta))
    } else {
      # For lag > 0, use asymptotic_variance_est
      var_result <- asymptotic_variance_est(t, c4_est, varlevyseed, Delta, est_trawl)
      av <- var_result$v

      if (av <= 0) {
        warning(paste("Non-positive variance at lag", lag_value, "- setting to NA"))
        av <- NA
      }

      std_error <- sqrt(av / (n * Delta))
    }

    # Compute test statistic
    if (!is.na(av) && av > 0) {
      t_stat <- sqrt(n * Delta) * (a_hat - a_par) / sqrt(av)
      p_value <- 2 * (1 - pnorm(abs(t_stat)))
      reject_05 <- abs(t_stat) > z_crit
      reject_01 <- abs(t_stat) > qnorm(0.995)
      reject_bonf <- abs(t_stat) > z_crit_bonf
    } else {
      t_stat <- NA
      p_value <- NA
      reject_05 <- NA
      reject_01 <- NA
      reject_bonf <- NA
    }

    # Store results
    test_results[[idx]] <- list(
      Lag = lag_value,
      Time = lag_value * Delta,
      T_stat = t_stat,
      P_value = p_value,
      Reject_05 = reject_05,
      Reject_01 = reject_01,
      Reject_Bonf = reject_bonf,
      a_nonpar = a_hat,
      a_par = a_par,
      std_error = std_error
    )
  }

  # Convert to data frame
  results_df <- do.call(rbind, lapply(test_results, function(x) {
    data.frame(
      Lag = x$Lag,
      Time = x$Time,
      T_stat = x$T_stat,
      P_value = x$P_value,
      Reject_05 = x$Reject_05,
      Reject_01 = x$Reject_01,
      Reject_Bonf = x$Reject_Bonf,
      a_nonpar = x$a_nonpar,
      a_par = x$a_par,
      std_error = x$std_error
    )
  }))

  # Add attributes
  attr(results_df, "alpha") <- alpha
  attr(results_df, "n_tests") <- K
  attr(results_df, "bonf_alpha") <- bonf_alpha
  attr(results_df, "parametric_name") <- parametric_name

  return(list(
    test_results = results_df,
    parametric_name = parametric_name
  ))
}

#####################
###Creating testing plots
create_testing_plots <- function(test_results, analysis_type = "Analysis",
                                 parametric_name = "Parametric",
                                 save_path = NULL) {

  # Extract attributes
  alpha <- attr(test_results, "alpha")
  if (is.null(alpha)) alpha <- 0.05

  K <- attr(test_results, "n_tests")
  if (is.null(K)) K <- nrow(test_results)

  bonf_alpha <- attr(test_results, "bonf_alpha")
  if (is.null(bonf_alpha)) bonf_alpha <- alpha / K

  z_crit_bonf <- attr(test_results, "z_crit_bonf")
  if (is.null(z_crit_bonf)) z_crit_bonf <- qnorm(1 - bonf_alpha / 2)

  # 1. Test statistic plot - NO LEGEND TITLES, PROPER MARGINS
  p1 <- ggplot(test_results, aes(x = Lag, y = T_stat)) +
    geom_line(color = "black", linewidth = 1.5) +
    geom_point(aes(color = Reject_05), size = 4.5) +
    geom_hline(aes(yintercept = qnorm(0.975), linetype = "Alpha"),
               color = "red", linewidth = 1.2) +
    geom_hline(aes(yintercept = -qnorm(0.975), linetype = "Alpha"),
               color = "red", linewidth = 1.2) +
    geom_hline(aes(yintercept = z_crit_bonf, linetype = "Bonferroni"),
               color = "blue", linewidth = 1.2) +
    geom_hline(aes(yintercept = -z_crit_bonf, linetype = "Bonferroni"),
               color = "blue", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray",
               linewidth = 0.8) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                       labels = c("FALSE" = "Not rejected (α=0.05)",
                                  "TRUE" = "Rejected (α=0.05)"),
                       name = NULL) +  # CHANGED: NULL instead of "Decision"
    scale_linetype_manual(values = c("Alpha" = "dashed", "Bonferroni" = "dotted"),
                          labels = c("Alpha" = paste0("α=", alpha),
                                     "Bonferroni" = paste0("Bonferroni α=",
                                                           round(bonf_alpha, 4))),
                          name = NULL) +  # CHANGED: NULL instead of "Critical values"
    labs(x = "Lag",
         y = TeX("$T(t)$")) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.box.margin = margin(t = 10),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      plot.margin = margin(t = 10, r = 15, b = 40, l = 15, unit = "pt")
    )

  # 2. P-value plot - NO LEGEND TITLES, PROPER MARGINS
  p2 <- ggplot(test_results, aes(x = Lag, y = P_value)) +
    geom_line(color = "black", linewidth = 1.5) +
    geom_point(aes(color = Reject_05), size = 4.5) +
    geom_hline(aes(yintercept = alpha, linetype = "Alpha"),
               color = "red", linewidth = 1.2) +
    geom_hline(aes(yintercept = bonf_alpha, linetype = "Bonferroni"),
               color = "blue", linewidth = 1.2) +
    scale_y_log10() +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                       labels = c("FALSE" = "Not rejected",
                                  "TRUE" = "Rejected (α=0.05)"),
                       name = NULL) +  # CHANGED: NULL instead of "Decision"
    scale_linetype_manual(values = c("Alpha" = "dashed", "Bonferroni" = "dotted"),
                          labels = c("Alpha" = paste0("α=", alpha),
                                     "Bonferroni" = paste0("Bonferroni α=",
                                                           round(bonf_alpha, 4))),
                          name = NULL) +  # CHANGED: NULL instead of "Critical values"
    labs(x = "Lag",
         y = "P-value (log scale)") +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.box.margin = margin(t = 10),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      plot.margin = margin(t = 10, r = 15, b = 40, l = 15, unit = "pt")
    )

  # 3. Difference plot - NO LEGEND TITLES, PROPER MARGINS
  test_results$difference <- test_results$a_nonpar - test_results$a_par
  test_results$diff_se <- test_results$std_error / sqrt(K)
  test_results$diff_lower <- test_results$difference - qnorm(0.975) * test_results$diff_se
  test_results$diff_upper <- test_results$difference + qnorm(0.975) * test_results$diff_se

  p3 <- ggplot(test_results, aes(x = Lag, y = difference)) +
    geom_ribbon(aes(ymin = diff_lower, ymax = diff_upper),
                alpha = 0.2, fill = "blue") +
    geom_line(color = "black", linewidth = 1.5) +
    geom_point(aes(color = Reject_05), size = 4.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red",
               linewidth = 1.2) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                       labels = c("FALSE" = "Not rejected",
                                  "TRUE" = "Rejected (α=0.05)"),
                       name = NULL) +  # CHANGED: NULL instead of "Decision"
    labs(x = "Lag",
         y = TeX("$\\hat{a}(t) - a_{par}(t)$")) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.box.margin = margin(t = 10),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      plot.margin = margin(t = 10, r = 15, b = 40, l = 15, unit = "pt")
    )

  # Save plots if path provided
  if (!is.null(save_path)) {
    # Individual plots with appropriate height to accommodate legend
    ggsave(paste0(save_path, "_test_statistics.png"),
           plot = p1, width = 12, height = 7)
    ggsave(paste0(save_path, "_test_statistics.eps"),
           plot = p1, width = 12, height = 7, device = "eps")

    ggsave(paste0(save_path, "_p_values.png"),
           plot = p2, width = 12, height = 7)
    ggsave(paste0(save_path, "_p_values.eps"),
           plot = p2, width = 12, height = 7, device = "eps")

    ggsave(paste0(save_path, "_difference.png"),
           plot = p3, width = 12, height = 7)
    ggsave(paste0(save_path, "_difference.eps"),
           plot = p3, width = 12, height = 7, device = "eps")

    # Combined plot
    p_combined <- grid.arrange(p1, p2, p3, ncol = 1)
    ggsave(paste0(save_path, "_combined.png"),
           plot = p_combined, width = 12, height = 18)
    ggsave(paste0(save_path, "_combined.eps"),
           plot = p_combined, width = 12, height = 18, device = "eps")
  }

  return(list(p1 = p1, p2 = p2, p3 = p3))
}

################################################################################
# 4. ANALYSIS FUNCTION FOR BOTH PROCESS TYPES (WITH TESTING)
################################################################################

analyse_trawl_process <- function(path1, analysis_type, Delta, n) {

  cat("\n========================================\n")
  cat("Analyzing", analysis_type, "trawl process\n")
  cat("========================================\n")

  # Create subdirectory for this analysis type
  plot_dir <- paste0("output/plots/", tolower(analysis_type))
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }

  ################################################################################
  # 4.1 BASIC PLOTS
  ################################################################################

  cat("  Creating basic plots...\n")

  # Plot simulated path
  png(paste0(plot_dir, "/simulated_path_", tolower(analysis_type), ".png"),
      width = 1200, height = 600, res = 150)
  plot(1:length(path1), path1, type = 'l', col = 'blue',
       main = paste("Simulated", analysis_type, "Trawl Process Path"),
       xlab = "Time Index", ylab = "Value")
  dev.off()

  setEPS()
  postscript(paste0(plot_dir, "/simulated_path_", tolower(analysis_type), ".eps"),
             width = 8, height = 4)
  plot(1:length(path1), path1, type = 'l', col = 'blue',
       main = paste("Simulated", analysis_type, "Trawl Process Path"),
       xlab = "Time Index", ylab = "Value")
  dev.off()

  ################################################################################
  # 4.2 ACF ANALYSIS
  ################################################################################

  cat("  Computing ACF...\n")

  # Compute empirical ACF
  acf_result <- acf(path1, lag.max = 100, plot = FALSE)

  # Plot basic ACF
  png(paste0(plot_dir, "/acf_basic_", tolower(analysis_type), ".png"),
      width = 1200, height = 800, res = 150)
  plot(acf_result, main = paste("ACF:", analysis_type, "Trawl Process"))
  dev.off()

  setEPS()
  postscript(paste0(plot_dir, "/acf_basic_", tolower(analysis_type), ".eps"),
             width = 8, height = 6)
  plot(acf_result, main = paste("ACF:", analysis_type, "Trawl Process"))
  dev.off()

  # Theoretical ACF at observation lags
  lags_obs <- acf_result$lag * Delta
  rho_theoretical <- sapply(lags_obs, rho)

  # Plot ACF with theoretical overlay
  png(paste0(plot_dir, "/acf_with_theoretical_", tolower(analysis_type), ".png"),
      width = 1200, height = 800, res = 150)
  plot(acf_result$lag, acf_result$acf, type = 'h', lwd = 2, col = 'blue',
       main = paste("ACF with Theoretical Overlay:", analysis_type),
       xlab = "Lag", ylab = "ACF", ylim = c(0, 1))
  lines(acf_result$lag, rho_theoretical, col = 'red', lwd = 2)
  legend("topright", legend = c("Empirical", "Theoretical"),
         col = c("blue", "red"), lwd = 2)
  dev.off()

  setEPS()
  postscript(paste0(plot_dir, "/acf_with_theoretical_", tolower(analysis_type), ".eps"),
             width = 8, height = 6)
  plot(acf_result$lag, acf_result$acf, type = 'h', lwd = 2, col = 'blue',
       main = paste("ACF with Theoretical Overlay:", analysis_type),
       xlab = "Lag", ylab = "ACF", ylim = c(0, 1))
  lines(acf_result$lag, rho_theoretical, col = 'red', lwd = 2)
  legend("topright", legend = c("Empirical", "Theoretical"),
         col = c("blue", "red"), lwd = 2)
  dev.off()

  ################################################################################
  # 4.3 NONPARAMETRIC TRAWL ESTIMATION
  ################################################################################

  cat("  Estimating trawl function nonparametrically...\n")

  # Estimate trawl function
  my_lag <- 101
  est_trawl <- nonpar_trawlest(path1, Delta, lag = my_lag)$a_hat

  # Plot nonparametric estimate
  lags_to_plot <- 1:my_lag
  l_seq <- seq(from = 0, to = (my_lag - 1), by = 1)
  plot_length <- my_lag
  trawl_values <- est_trawl[lags_to_plot]

  png(paste0(plot_dir, "/trawl_nonparametric_", tolower(analysis_type), ".png"),
      width = 1200, height = 800, res = 150)
  plot(l_seq, trawl_values, type = 'p', pch = 19, col = 'black',
       main = paste("Nonparametric trawl function estimate:", analysis_type),
       xlab = "Lag", ylab = TeX("$\\hat{a}(\\cdot)$"))
  lines(l_seq, sapply(l_seq, function(l) a(l * Delta)),
        col = 'red', lwd = 2)
  legend("topright", legend = c("Estimated", "True"),
         col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
  dev.off()

  setEPS()
  postscript(paste0(plot_dir, "/trawl_nonparametric_", tolower(analysis_type), ".eps"),
             width = 8, height = 6)
  plot(l_seq, trawl_values, type = 'p', pch = 19, col = 'black',
       main = paste("Nonparametric trawl function estimate:", analysis_type),
       xlab = "Lag", ylab = TeX("$\\hat{a}(\\cdot)$"))
  lines(l_seq, sapply(l_seq, function(l) a(l * Delta)),
        col = 'red', lwd = 2)
  legend("topright", legend = c("Estimated", "True"),
         col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))
  dev.off()

  ################################################################################
  # 4.4 ASYMPTOTIC VARIANCE ESTIMATION
  ################################################################################

  cat("  Computing asymptotic variance...\n")

  # Estimate asymptotic variance using custom functions
  c4 <- c4est(path1, Delta)
  av_vector <- numeric(my_lag)

  for (i in 1:my_lag) {
    av_vector[i] <- asymptotic_variance_est((i - 1) * Delta, c4,
                                            varlevyseed = 1, Delta, est_trawl)$v
  }

  var_diag <- av_vector  # For compatibility with rest of code

  # Plot asymptotic standard errors
  plot_length_av <- min(31, my_lag)

  png(paste0(plot_dir, "/asymptotic_variance_", tolower(analysis_type), ".png"),
      width = 1200, height = 800, res = 150)
  plot(l_seq[1:plot_length_av], sqrt(var_diag[1:plot_length_av]),
       type = 'l', lwd = 2, col = 'blue',
       main = paste("Asymptotic Standard Errors:", analysis_type),
       xlab = "Lag", ylab = "Standard Error")
  dev.off()

  setEPS()
  postscript(paste0(plot_dir, "/asymptotic_variance_", tolower(analysis_type), ".eps"),
             width = 8, height = 6)
  plot(l_seq[1:plot_length_av], sqrt(var_diag[1:plot_length_av]),
       type = 'l', lwd = 2, col = 'blue',
       main = paste("Asymptotic Standard Errors:", analysis_type),
       xlab = "Lag", ylab = "Standard Error")
  dev.off()

  ################################################################################
  # 4.5 CONFIDENCE INTERVALS
  ################################################################################

  cat("  Constructing confidence intervals...\n")

  # Standard errors
  std_errors <- sqrt(var_diag / (n * Delta))

  # 95% confidence intervals
  z <- qnorm(0.975)
  lower <- est_trawl - z * std_errors
  upper <- est_trawl + z * std_errors

  # Plot with confidence intervals
  plot_length <- min(71, length(est_trawl))
  xx <- l_seq[1:plot_length]
  y1 <- upper[1:plot_length]
  y2 <- lower[1:plot_length]
  y3 <- est_trawl[1:plot_length]

  df <- data.frame(xx, y1, y2, y3)
  mdf <- melt(data = df, id.vars = "xx", measure.vars = c("y1", "y2", "y3"))

  true_trawl <- data.frame(
    Lag = l_seq[1:plot_length],
    TrueValue = sapply(l_seq[1:plot_length], function(l) a(l * Delta))
  )

  mdf$LegendGroup <- ifelse(mdf$variable == "y3", "Estimated function",
                            ifelse(mdf$variable %in% c("y1", "y2"), "Confidence bounds", NA))
  true_trawl$LegendGroup <- "True trawl function"

  g <- ggplot() +
    geom_point(data = mdf, aes(x = xx, y = value, colour = LegendGroup), size = 3,
               shape = ifelse(mdf$variable == "y3", 19, 1)) +
    geom_line(data = true_trawl, aes(x = Lag, y = TrueValue, colour = LegendGroup), linewidth = 1.5) +
    geom_point(data = true_trawl, aes(x = Lag, y = TrueValue, colour = LegendGroup),
               size = 2, shape = 19) +
    scale_color_manual(name = "Legend",
                       values = c("Confidence bounds" = "blue",
                                  "Estimated function" = "black",
                                  "True trawl function" = "red")) +
    guides(colour = guide_legend(override.aes = list(
      shape = c(1, 19, 19),
      linetype = c(0, 0, 1),
      size = c(3, 3, 1.5)
    ))) +
    xlab("l") +
    ylab(TeX("$\\hat{a}(\\cdot)$")) +
    ggtitle(paste("Nonparametric trawl estimate with 95% CI:", analysis_type)) +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(paste0(plot_dir, "/confidence_intervals_", tolower(analysis_type), ".png"),
         plot = g, width = 12, height = 8)
  ggsave(paste0(plot_dir, "/confidence_intervals_", tolower(analysis_type), ".eps"),
         plot = g, width = 14, height = 8, device = "eps")

  ################################################################################
  # 4.6 PARAMETRIC FITTING
  ################################################################################

  cat("  Fitting parametric models...\n")

  # Fit exponential trawl
  fit_exp <- nlsLM(y3 ~ a0 * exp(-lambda * xx),
                   data = df,
                   start = list(a0 = 5, lambda = 0.3),
                   control = nls.lm.control(maxiter = 200))

  # Fit piecewise exponential trawl
  fit_piecewise <- nlsLM(y3 ~ ifelse(xx < c1,
                                     a0 * exp(-lambda1 * xx),
                                     a0 * exp(-lambda1 * c1) * exp(-lambda2 * (xx - c1))),
                         data = df,
                         start = list(a0 = 5, lambda1 = 0.1, lambda2 = 0.5, c1 = 2),
                         control = nls.lm.control(maxiter = 200))

  # Fit piecewise exponential with TRUE breakpoint FIXED at lag 20
  df$c1_true <- 20  # TRUE breakpoint in LAG units
  fit_piecewise_true <- nlsLM(y3 ~ ifelse(xx <= c1_true,
                                          a0 * exp(-lambda1 * xx),
                                          a0 * exp(-lambda1 * c1_true) * exp(-lambda2 * (xx - c1_true))),
                              data = df,
                              start = list(a0 = 5, lambda1 = 0.1, lambda2 = 0.5),
                              lower = c(0.1, 0.01, 0.01),
                              upper = c(10, 2, 2),
                              control = nls.lm.control(maxiter = 200))

  # Extract parameters
  params_exp <- coef(fit_exp)
  params_piecewise <- coef(fit_piecewise)
  params_piecewise_true <- coef(fit_piecewise_true)

  cat("    Exponential fit parameters:\n")
  print(params_exp)
  cat("    Piecewise exponential fit parameters (optimised breakpoint):\n")
  print(params_piecewise)
  cat("    Piecewise exponential fit parameters (TRUE breakpoint at lag 20):\n")
  print(params_piecewise_true)

  # Plot parametric fits
  # Create predictions for the same range as df (plot_length points)
  pred_data_plot <- data.frame(xx = xx)
  df$pred_exp <- predict(fit_exp, newdata = pred_data_plot)
  df$pred_piecewise <- predict(fit_piecewise, newdata = pred_data_plot)

  g_param <- ggplot(df, aes(x = xx)) +
    geom_point(aes(y = y3), color = "black", size = 2) +
    geom_line(aes(y = pred_exp, color = "Exponential"), linewidth = 1) +
    geom_line(aes(y = pred_piecewise, color = "Piecewise exp."), linewidth = 1) +
    geom_line(data = true_trawl, aes(x = Lag, y = TrueValue, color = "True"), linewidth = 1.5) +
    scale_color_manual(name = "Model",
                       values = c("Exponential" = "orange",
                                  "Piecewise exp." = "purple",
                                  "True" = "red")) +
    xlab("l") +
    ylab(TeX("$\\hat{a}(\\cdot)$")) +
    ggtitle(paste("Parametric fits:", analysis_type)) +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(paste0(plot_dir, "/parametric_fits_", tolower(analysis_type), ".png"),
         plot = g_param, width = 12, height = 8)
  ggsave(paste0(plot_dir, "/parametric_fits_", tolower(analysis_type), ".eps"),
         plot = g_param, width = 14, height = 8, device = "eps")

  ################################################################################

  ################################################################################
  # 4.7 COMPREHENSIVE COMPARISON PLOT
  ################################################################################

  cat("  Creating comprehensive comparison plot...\n")

  plot_length_comp <- min(100, length(est_trawl))
  l_seq_comp <- seq(from = 0, to = 100, by = 1)

  est_trawl_method1 <- sapply(l_seq_comp[1:plot_length_comp], function(l) a(l * Delta))  # True
  est_trawl_method2 <- est_trawl[1:plot_length_comp]                   # Nonparametric

  # Create prediction data for parametric fits
  pred_data <- data.frame(xx = l_seq_comp[1:plot_length_comp])
  est_trawl_method3 <- predict(fit_piecewise_true, newdata = data.frame(xx = l_seq_comp[1:plot_length_comp], c1_true = 20))     # Piecewise exp
  est_trawl_method4 <- predict(fit_exp, newdata = pred_data)           # Simple exp
  # est_trawl_method5 <- trawl_hybrid[1:plot_length_comp]              # OLD: Hybrid with true breakpoint
  # Note: Method5 is now defined later as piecewise with optimised breakpoint

  # Prepare data for plotting
  xx <- l_seq_comp[1:plot_length_comp]
  y1 <- upper[1:plot_length_comp]
  y2 <- lower[1:plot_length_comp]
  y3 <- est_trawl_method2

  df_comp <- data.frame(xx, y1, y2, y3)
  df_comp$fit_exp <- est_trawl_method4
  df_comp$fit_piecewise <- est_trawl_method3

  mdf_comp <- melt(data = df_comp, id.vars = "xx", measure.vars = c("y1", "y2", "y3"))

  true_trawl_values <- data.frame(
    Lag = l_seq_comp[1:plot_length_comp],
    TrueValue = est_trawl_method1
  )

  # Assign legend groups
  mdf_comp$LegendGroup <- ifelse(mdf_comp$variable == "y3", "Estimated function",
                                 ifelse(mdf_comp$variable %in% c("y1", "y2"), "Confidence bounds", NA))
  true_trawl_values$LegendGroup <- "True trawl function"
  df_comp$LegendGroup_exp <- "Exponential fit"
  df_comp$LegendGroup_piecewise <- "Piecewise exponential fit"

  # Comprehensive comparison plot
  g_comparison <- ggplot() +
    # Confidence bounds and estimated function
    geom_point(data = mdf_comp, aes(x = xx, y = value, colour = LegendGroup), size = 3,
               shape = ifelse(mdf_comp$variable == "y3", 19, 1)) +

    # True trawl function
    geom_line(data = true_trawl_values, aes(x = Lag, y = TrueValue, colour = LegendGroup), linewidth = 1.5) +
    geom_point(data = true_trawl_values, aes(x = Lag, y = TrueValue, colour = LegendGroup),
               size = 2, shape = 19) +

    # Exponential fit
    geom_line(data = df_comp, aes(x = xx, y = fit_exp, colour = LegendGroup_exp), linewidth = 1) +

    # Piecewise exponential fit
    geom_line(data = df_comp, aes(x = xx, y = fit_piecewise, colour = LegendGroup_piecewise), linewidth = 1) +

    # Manual color legend
    scale_color_manual(name = "Legend",
                       values = c("Confidence bounds" = "blue",
                                  "Estimated function" = "black",
                                  "True trawl function" = "red",
                                  "Exponential fit" = "orange",
                                  "Piecewise exponential fit" = "purple")) +

    guides(colour = guide_legend(override.aes = list(
      shape = c(1, 19, 19, NA, NA),
      linetype = c(0, 0, 1, 1, 1),
      size = c(3, 3, 1.5, 1, 1)
    ))) +

    xlab("l") +
    ylab(TeX("$\\hat{a}(\\cdot)$")) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "right",
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          axis.title = element_text(size = 22),
          axis.text = element_text(size = 20))

  # Save PNG
  ggsave(paste0(plot_dir, "/comparison_all_methods_", tolower(analysis_type), ".png"),
         plot = g_comparison, width = 12, height = 8)

  # Save EPS
  ggsave(paste0(plot_dir, "/comparison_all_methods_", tolower(analysis_type), ".eps"),
         plot = g_comparison, width = 14, height = 8, device = "eps")

  ################################################################################
  # 4.8 FORMAL HYPOTHESIS TESTING - COMPREHENSIVE (LAG VALUES 0-60)
  ################################################################################

  cat("  Performing comprehensive formal hypothesis tests (lag values 0-60)...\n")

  # Select lag VALUES to test (ALL lags from 0 to 60)
  lag_values_to_test <- 0:60

  cat("    Testing at", length(lag_values_to_test), "lag values from", min(lag_values_to_test),
      "to", max(lag_values_to_test), "\n")

  # Estimate c4 using user's function
  cat("    Computing c4 estimate...\n")
  c4_estimate <- c4est(path1, Delta)
  cat("    c4 estimate:", c4_estimate, "\n")

  # Test exponential specification
  cat("    Testing exponential specification at all lags 0-60...\n")
  a_par_exp <- predict(fit_exp, newdata = data.frame(xx = l_seq))
  test_exp <- perform_trawl_specification_test(
    est_trawl = est_trawl,
    a_par_fitted = a_par_exp,
    Delta = Delta,
    n = n,
    path_data = path1,
    lag_values = lag_values_to_test,
    c4_est = c4_estimate,
    alpha = 0.05,
    parametric_name = "Exponential"
  )

  # Test against TRUE piecewise exponential
  cat("    Testing piecewise exponential specification (TRUE parameters) at all lags 0-60...\n")
  a_par_true_piecewise <- sapply(l_seq, function(x) a(x * Delta))
  test_piecewise <- perform_trawl_specification_test(
    est_trawl = est_trawl,
    a_par_fitted = a_par_true_piecewise,
    Delta = Delta,
    n = n,
    path_data = path1,
    lag_values = lag_values_to_test,
    c4_est = c4_estimate,
    alpha = 0.05,
    parametric_name = "Piecewise Exponential (True)"
  )

  # Save test results
  write.csv(test_exp$test_results,
            file = paste0("output/data/hypothesis_tests_comprehensive_exponential_",
                          tolower(analysis_type), ".csv"),
            row.names = FALSE)

  write.csv(test_piecewise$test_results,
            file = paste0("output/data/hypothesis_tests_comprehensive_piecewise_",
                          tolower(analysis_type), ".csv"),
            row.names = FALSE)

  # Create testing plots for exponential specification
  cat("    Creating comprehensive hypothesis testing plots for exponential...\n")
  testing_plots_exp <- create_testing_plots(
    test_results = test_exp$test_results,
    analysis_type = analysis_type,
    parametric_name = "Exponential",
    save_path = paste0(plot_dir, "/testing_comprehensive_exponential_", tolower(analysis_type))
  )

  # Create testing plots for piecewise exponential specification
  cat("    Creating comprehensive hypothesis testing plots for piecewise exponential...\n")
  testing_plots_piecewise <- create_testing_plots(
    test_results = test_piecewise$test_results,
    analysis_type = analysis_type,
    parametric_name = "Piecewise exp.",
    save_path = paste0(plot_dir, "/testing_comprehensive_piecewise_", tolower(analysis_type))
  )

  # # Create summary tables for the paper (representative subset of lags)
  # summary_table_exp <- create_summary_table(
  #   test_exp$test_results,
  #   lags_to_include = c(5, 10, 15, 20, 25, 30, 40, 50),
  #   digits = 2
  # )
  #
  # summary_table_piecewise <- create_summary_table(
  #   test_piecewise$test_results,
  #   lags_to_include = c(5, 10, 15, 20, 25, 30, 40, 50),
  #   digits = 2
  # )
  #
  # cat("\n    Summary Table - Exponential (for paper, representative lags):\n")
  # print(summary_table_exp)
  #
  # cat("\n    Summary Table - Piecewise Exponential (for paper, representative lags):\n")
  # print(summary_table_piecewise)
  #
  # write.csv(summary_table_exp,
  #           file = paste0("output/data/test_summary_table_comprehensive_exponential_",
  #                         tolower(analysis_type), ".csv"),
  #           row.names = FALSE)
  #
  # write.csv(summary_table_piecewise,
  #           file = paste0("output/data/test_summary_table_comprehensive_piecewise_",
  #                         tolower(analysis_type), ".csv"),
  #           row.names = FALSE)

  # Store test results for return
  test_results_exp <- test_exp
  test_results_piecewise <- test_piecewise

  ################################################################################
  # 4.9 DATA-DRIVEN BREAKPOINT IDENTIFICATION VIA HYPOTHESIS TESTS
  ################################################################################

  cat("  Using hypothesis tests to identify structural breakpoint...\n")

  # Use the comprehensive test results from exponential specification
  # Restrict search to lags 5-40 (breakpoint unlikely to be at very short or very long lags)
  search_range <- test_exp$test_results$Lag >= 5 & test_exp$test_results$Lag <= 40
  test_results_search <- test_exp$test_results[search_range, ]

  # Find lag where test statistic is largest (indicates strongest misspecification)
  max_t_idx <- which.max(abs(test_results_search$T_stat))
  estimated_breakpoint_lag <- test_results_search$Lag[max_t_idx]
  estimated_breakpoint_time <- estimated_breakpoint_lag * Delta

  cat(sprintf("    Identified breakpoint at lag %d (time %.2f)\n",
              estimated_breakpoint_lag, estimated_breakpoint_time))
  cat(sprintf("    True breakpoint is at time %.2f\n", 2.0))

  # Create HYBRID estimator with ESTIMATED breakpoint
  cat("    Creating hybrid estimator: nonparametric until breakpoint, then exponential...\n")

  # Use nonparametric estimates up to the breakpoint
  # Then fit exponential to data AFTER the breakpoint

  # Get indices for data after the estimated breakpoint
  # Array index: lag 0 is at index 1, so lag L is at index L+1
  breakpoint_array_idx <- estimated_breakpoint_lag + 1
  after_breakpoint_idx <- (breakpoint_array_idx + 1):min(length(est_trawl), plot_length)

  # Fit exponential to data AFTER breakpoint using log-linear regression
  # This ensures perfect continuity at the breakpoint

  # Get the nonparametric value AT the breakpoint for continuity
  a_at_breakpoint <- est_trawl[breakpoint_array_idx]

  # Fit exponential ONLY to estimate lambda
  # Use only first 30 lags after breakpoint for more stable fitting
  max_lags_for_fit <- min(30, length(after_breakpoint_idx))
  fit_idx <- after_breakpoint_idx[1:max_lags_for_fit]

  x_after <- (l_seq[fit_idx] * Delta) - estimated_breakpoint_time
  y_after <- est_trawl[fit_idx]

  # Filter out problematic values: keep only values > 10% of breakpoint value
  valid_mask <- y_after > (0.1 * a_at_breakpoint) & y_after > 0
  x_after_valid <- x_after[valid_mask]
  y_after_valid <- y_after[valid_mask]

  # Fit: y = a_at_breakpoint * exp(-lambda * (x - breakpoint_time))
  # Taking logs: log(y) = log(a_at_breakpoint) - lambda * x
  if (sum(valid_mask) >= 5) {  # Need at least 5 points
    log_y_after <- log(y_after_valid)
    log_a_bp <- log(a_at_breakpoint)

    # Fit: log_y = log_a_bp - lambda * x
    fit_log <- lm(log_y_after ~ x_after_valid)
    lambda_fitted <- -coef(fit_log)[2]  # Negative of slope

    # Bound lambda to reasonable range [0.01, 2.0]
    lambda_fitted <- pmax(0.01, pmin(2.0, lambda_fitted))
  } else {
    # Not enough valid points, use a default reasonable value
    lambda_fitted <- 0.5
    cat("    Warning: Not enough valid points for fitting, using default lambda=0.5\n")
  }

  # Store in a list for consistency
  fit_exp_after <- list(
    lambda = lambda_fitted,
    a_fixed = a_at_breakpoint,
    c_fixed = estimated_breakpoint_time
  )

  cat(sprintf("    Fitted exponential (after breakpoint): a0=%.3f, lambda=%.3f\n",
              fit_exp_after$a_fixed, fit_exp_after$lambda))

  # Fit PIECEWISE EXPONENTIAL with c1 as FREE parameter (for comparison)
  cat("    Fitting piecewise exponential with c1 as free parameter...\n")

  # Use better starting values
  a0_start <- max(df$y3, na.rm = TRUE)

  fit_piecewise_optimised <- nlsLM(
    y3 ~ ifelse(xx <= c1,
                a0 * exp(-lambda1 * xx),
                a0 * exp(-lambda1 * c1) * exp(-lambda2 * (xx - c1))),
    data = df,
    start = list(a0 = a0_start, lambda1 = 0.1, lambda2 = 0.5, c1 = 20),
    lower = c(0.1, 0.01, 0.01, 5),
    upper = c(10, 2, 2, 60),
    control = nls.lm.control(maxiter = 1000)
  )

  coefs_piecewise_opt <- coef(fit_piecewise_optimised)
  cat(sprintf("    Fitted piecewise (optimised): a0=%.3f, lambda1=%.3f, lambda2=%.3f, c1=%.3f (lag %.0f)\n",
              coefs_piecewise_opt["a0"], coefs_piecewise_opt["lambda1"],
              coefs_piecewise_opt["lambda2"], coefs_piecewise_opt["c1"],
              coefs_piecewise_opt["c1"]))

  # Create diagnostic plot showing the process
  cat("    Creating diagnostic plot...\n")

  # Panel 1: Test statistics showing where breakpoint is identified (subset for clarity)
  test_search_plot <- test_exp$test_results[test_exp$test_results$Lag >= 5 &
                                              test_exp$test_results$Lag <= 40, ]

  p_test_search <- ggplot(test_search_plot, aes(x = Lag, y = abs(T_stat))) +
    geom_line(color = "black", linewidth = 1.5) +
    geom_point(size = 3) +
    geom_vline(xintercept = estimated_breakpoint_lag,
               linetype = "dashed", color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = 20, linetype = "dotted", color = "red", linewidth = 1.2) +
    annotate("text", x = estimated_breakpoint_lag, y = max(abs(test_search_plot$T_stat)) * 0.9,
             label = paste0("Estimated\nbreakpoint\n(lag ", estimated_breakpoint_lag, ")"),
             color = "blue", size = 6, hjust = -0.1) +
    annotate("text", x = 20, y = max(abs(test_search_plot$T_stat)) * 0.7,
             label = "True\nbreakpoint\n(lag 20)",
             color = "red", size = 6, hjust = 1.1) +
    labs(x = "Lag",
         y = TeX("$|T(t)|$")) +
    theme_minimal(base_size = 20) +
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 20))

  # Panel 2: Fitted trawl functions comparison

  # Create hybrid estimator: nonparametric until breakpoint, exponential after
  xx_vals <- l_seq[1:plot_length]

  hybrid_estimated <- est_trawl[1:plot_length]  # Start with nonparametric

  # Replace values after breakpoint with exponential fit
  after_idx <- (breakpoint_array_idx + 1):plot_length
  if (length(after_idx) > 0) {
    hybrid_estimated[after_idx] <- fit_exp_after$a_fixed *
      exp(-fit_exp_after$lambda * (l_seq[after_idx] * Delta - fit_exp_after$c_fixed))
  }

  cat(sprintf("    Hybrid: using nonparametric for lags 0-%d, exponential for lags %d+\n",
              estimated_breakpoint_lag, estimated_breakpoint_lag + 1))

  df_compare <- data.frame(
    xx = xx_vals,
    y3 = est_trawl[1:plot_length],
    fit_exp = predict(fit_exp, newdata = data.frame(xx = xx_vals)),
    fit_piecewise_optimised = predict(fit_piecewise_optimised, newdata = data.frame(xx = xx_vals)),
    hybrid_estimated = hybrid_estimated
  )

  true_trawl_compare <- data.frame(
    Lag = l_seq[1:plot_length],
    TrueValue = sapply(l_seq[1:plot_length], function(l) a(l * Delta))
  )

  p_fits_comparison <- ggplot() +
    geom_point(data = df_compare, aes(x = xx, y = y3, color = "Nonparametric"),
               size = 2) +
    geom_line(data = true_trawl_compare, aes(x = Lag, y = TrueValue, color = "True trawl"),
              linewidth = 1.5) +
    geom_line(data = df_compare, aes(x = xx, y = fit_exp, color = "Exponential"),
              linewidth = 1.2) +
    geom_line(data = df_compare, aes(x = xx, y = hybrid_estimated,
                                     color = "Hybrid"),
              linewidth = 1.2, linetype = "dashed") +
    geom_line(data = df_compare, aes(x = xx, y = fit_piecewise_optimised,
                                     color = "Piecewise exp."),
              linewidth = 1.2) +
    geom_vline(xintercept = estimated_breakpoint_lag,
               linetype = "dashed", color = "blue", linewidth = 1.2) +
    geom_vline(xintercept = coefs_piecewise_opt["c1"],
               linetype = "dotted", color = "purple", linewidth = 1.2) +
    scale_color_manual(
      name = "Model",
      values = c(
        "Nonparametric" = "gray50",
        "True trawl" = "red",
        "Exponential" = "orange",
        "Hybrid" = "blue",
        "Piecewise exp." = "purple"
      )
    ) +
    labs(x = "Lag",
         y = TeX("$\\hat{a}(\\cdot)$")) +
    theme_minimal(base_size = 20) +
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 20),
          legend.position = "right",
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20))


  # Save individual panels separately
  ggsave(paste0(plot_dir, "/breakpoint_test_statistics_", tolower(analysis_type), ".png"),
         plot = p_test_search, width = 14, height = 6)
  ggsave(paste0(plot_dir, "/breakpoint_test_statistics_", tolower(analysis_type), ".eps"),
         plot = p_test_search, width = 14, height = 6, device = "eps")

  ggsave(paste0(plot_dir, "/breakpoint_model_comparison_", tolower(analysis_type), ".png"),
         plot = p_fits_comparison, width = 14, height = 6)
  ggsave(paste0(plot_dir, "/breakpoint_model_comparison_", tolower(analysis_type), ".eps"),
         plot = p_fits_comparison, width = 14, height = 6, device = "eps")

  # Combine panels
  g_breakpoint_diagnostic <- grid.arrange(
    p_test_search,
    p_fits_comparison,
    ncol = 1
  )

  # Save the diagnostic plot
  ggsave(paste0(plot_dir, "/breakpoint_identification_", tolower(analysis_type), ".png"),
         plot = g_breakpoint_diagnostic, width = 14, height = 12)
  ggsave(paste0(plot_dir, "/breakpoint_identification_", tolower(analysis_type), ".eps"),
         plot = g_breakpoint_diagnostic, width = 14, height = 12, device = "eps")

  cat(sprintf("    Saved breakpoint identification diagnostic plot (combined and separate panels)\n"))

  ################################################################################
  # CREATE STANDALONE HYBRID ESTIMATION PLOT
  ################################################################################

  cat("  Creating standalone hybrid estimation plot...\n")

  # Create simple hybrid estimation plot
  g_hybrid <- ggplot() +
    geom_point(data = df_compare, aes(x = xx, y = y3, color = "Nonparametric"), size = 2) +
    geom_line(aes(x = xx_vals, y = hybrid_estimated, color = "Hybrid"), linewidth = 1) +
    geom_vline(xintercept = estimated_breakpoint_lag, linetype = "dashed", color = "gray") +
    geom_line(data = true_trawl_compare, aes(x = Lag, y = TrueValue, color = "True"), linewidth = 1.5) +
    scale_color_manual(name = "Method",
                       values = c("Nonparametric" = "black",
                                  "Hybrid" = "green",
                                  "True" = "red")) +
    xlab("l") +
    ylab(TeX("$\\hat{a}(\\cdot)$")) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "right",
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          axis.title = element_text(size = 22),
          axis.text = element_text(size = 20))

  ggsave(paste0(plot_dir, "/hybrid_estimation_", tolower(analysis_type), ".png"),
         plot = g_hybrid, width = 12, height = 8)
  ggsave(paste0(plot_dir, "/hybrid_estimation_", tolower(analysis_type), ".eps"),
         plot = g_hybrid, width = 14, height = 8, device = "eps")

  cat(sprintf("    Saved standalone hybrid estimation plot\n"))

  # Print summary
  cat("\n  ========================================\n")
  cat("  Breakpoint Identification Summary\n")
  cat("  ========================================\n")
  cat(sprintf("  Estimated breakpoint: lag %d (time %.2f)\n",
              estimated_breakpoint_lag, estimated_breakpoint_time))
  cat(sprintf("  True breakpoint:      lag 20 (time %.2f)\n", 2.0))
  cat(sprintf("  Estimation error:     %.2f lags (%.3f time units)\n",
              abs(estimated_breakpoint_lag - 20),
              abs(estimated_breakpoint_time - 2.0)))
  cat("  ========================================\n\n")

  # Create FULL-LENGTH hybrid estimator for forecasting comparison
  # (the one created above was only for plotting)
  hybrid_estimated_full <- est_trawl  # Start with full nonparametric

  # Replace values after breakpoint with exponential fit (with perfect continuity)
  breakpoint_array_idx_full <- estimated_breakpoint_lag + 1
  after_idx_full <- (breakpoint_array_idx_full + 1):length(est_trawl)

  if (length(after_idx_full) > 0) {
    hybrid_estimated_full[after_idx_full] <- fit_exp_after$a_fixed *
      exp(-fit_exp_after$lambda * (l_seq[after_idx_full] * Delta - fit_exp_after$c_fixed))
  }

  # Create FULL-LENGTH piecewise exponential (optimised breakpoint) for forecasting
  piecewise_optimised_full <- predict(fit_piecewise_optimised,
                                      newdata = data.frame(xx = l_seq[1:length(est_trawl)]))

  # Save estimators for forecasting comparison
  est_trawl_method5 <- piecewise_optimised_full  # Piecewise with optimised breakpoint
  est_trawl_method6 <- hybrid_estimated_full     # Hybrid with data-driven breakpoint

  ################################################################################
  # 4.10 FORECASTING COMPARISON
  ################################################################################

  cat("  Performing forecasting analysis...\n")

  # Split data into training and test sets
  training_length <- floor(0.8 * length(path1))
  training <- path1[1:training_length]
  test <- path1[(training_length + 1):length(path1)]
  mean_training <- mean(training)

  # Create list of estimation methods (matching diagnostic plot)
  est_trawl_list <- list(
    True_trawl = est_trawl_method1,                      # Red
    Nonparametric = est_trawl_method2,                   # Gray
    Exponential = est_trawl_method4,                     # Orange (misspecified)
    Piecewise = est_trawl_method5,                       # Purple (optimised breakpoint)
    Hybrid = est_trawl_method6                           # Blue (test-based breakpoint)
  )

  # Compare forecasting performance across different horizons
  k_values <- 1:10  # Forecast horizons: h = k * Delta (focus on 10 lags)
  all_accuracy <- data.frame()

  for (k in k_values) {
    res <- compare_trawl_methods(est_trawl_list, Delta, k, training, test)

    # Extract accuracy measures (RMSE only)
    acc_df <- do.call(rbind, lapply(seq_along(res$accuracy), function(i) {
      data.frame(
        Method = names(est_trawl_list)[i],
        k = k,
        RMSE = res$accuracy[[i]]["Test set", "RMSE"]
      )
    }))

    all_accuracy <- rbind(all_accuracy, acc_df)
  }

  # Plot forecast accuracy comparison (RMSE only)
  # Use same colors as diagnostic plot
  p_forecast <- ggplot(all_accuracy, aes(x = factor(k), y = RMSE, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c(
      "True_trawl" = "red",
      "Nonparametric" = "gray50",
      "Exponential" = "orange",
      "Piecewise" = "purple",
      "Hybrid" = "blue"
    )) +
    labs(
      x = "Forecast horizon (k)",
      y = "RMSE",
      fill = "Estimation method") +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 22, ),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 20, ),
    )

  # Save PNG
  ggsave(paste0(plot_dir, "/forecast_accuracy_", tolower(analysis_type), ".png"),
         plot = p_forecast, width = 12, height = 8)

  # Save EPS
  ggsave(paste0(plot_dir, "/forecast_accuracy_", tolower(analysis_type), ".eps"),
         plot = p_forecast, width = 14, height = 8, device = "eps")

  # Create plot showing percentage increase in RMSE for Exponential vs alternatives
  # Reshape data to compute percentage increases
  rmse_wide <- all_accuracy %>%
    tidyr::pivot_wider(names_from = Method, values_from = RMSE)

  # Calculate percentage increase: ((Exponential - Alternative) / Alternative) * 100
  pct_increase <- data.frame(
    k = rmse_wide$k,
    vs_True = ((rmse_wide$Exponential - rmse_wide$True_trawl) / rmse_wide$True_trawl) * 100,
    vs_Nonparametric = ((rmse_wide$Exponential - rmse_wide$Nonparametric) / rmse_wide$Nonparametric) * 100,
    vs_Piecewise = ((rmse_wide$Exponential - rmse_wide$Piecewise) / rmse_wide$Piecewise) * 100,
    vs_Hybrid = ((rmse_wide$Exponential - rmse_wide$Hybrid) / rmse_wide$Hybrid) * 100
  )

  # Reshape for plotting
  pct_increase_long <- pct_increase %>%
    tidyr::pivot_longer(cols = starts_with("vs_"),
                        names_to = "Comparison",
                        values_to = "Percentage_Increase") %>%
    mutate(Comparison = case_when(
      Comparison == "vs_True" ~ "vs_True_trawl",
      TRUE ~ Comparison
    ))

  # Create plot
  p_pct_increase <- ggplot(pct_increase_long, aes(x = factor(k), y = Percentage_Increase,
                                                  color = Comparison, group = Comparison)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    scale_color_manual(
      values = c(
        "vs_True_trawl" = "red",
        "vs_Nonparametric" = "gray50",
        "vs_Piecewise" = "purple",
        "vs_Hybrid" = "blue"
      ),
      labels = c(
        "vs_True_trawl" = "vs True trawl",
        "vs_Nonparametric" = "vs Nonparametric",
        "vs_Piecewise" = "vs Piecewise",
        "vs_Hybrid" = "vs Hybrid"
      )
    ) +
    labs(x = "Forecast horizon (k)",
         y = "RMSE increase (%)",
         color = "Comparison") +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          legend.position = "right")

  # Save PNG
  ggsave(paste0(plot_dir, "/forecast_exponential_penalty_", tolower(analysis_type), ".png"),
         plot = p_pct_increase, width = 12, height = 8)

  # Save EPS
  ggsave(paste0(plot_dir, "/forecast_exponential_penalty_", tolower(analysis_type), ".eps"),
         plot = p_pct_increase, width = 14, height = 8, device = "eps")

  # Save accuracy results to CSV
  write.csv(all_accuracy,
            file = paste0("output/data/forecast_accuracy_", tolower(analysis_type), ".csv"),
            row.names = FALSE)

  cat("  Analysis complete for", analysis_type, "\n")

  # Return results for potential further analysis
  return(list(
    est_trawl = est_trawl,
    std_errors = std_errors,
    fit_exp = fit_exp,
    fit_piecewise = fit_piecewise,
    hybrid_estimated = hybrid_estimated_full,
    test_results_exp = test_results_exp,
    test_results_piecewise = test_results_piecewise,
    forecast_accuracy = all_accuracy
  ))
}


##################additional functio
compare_trawl_methods <- function(est_trawl_list, Delta, k, training, test) {
  n_methods <- length(est_trawl_list)
  training_length <- length(training)
  mean_training <- mean(training)

  accuracy_list <- list()
  forecast_list <- list()

  for (i in seq_len(n_methods)) {
    slices <- LebA_slice_est_approx(est_trawl_list[[i]], Delta, h = k * Delta)
    intersec <- slices$LebAintersection / slices$LebA
    differ <- slices$LebAsetdifference / slices$LebA

    forecast_vec <- numeric(length(test))
    forecast_vec[1] <- training[training_length] * intersec + mean_training * differ

    for (j in 2:length(test)) {
      forecast_vec[j] <- test[j-1] * intersec + mean_training * differ
    }

    forecast_list[[i]] <- forecast_vec
    accuracy_list[[i]] <- forecast::accuracy(forecast_vec, test)
  }

  # Prepare accuracy data frame for plotting
  accuracy_df <- do.call(rbind, lapply(seq_len(n_methods), function(i) {
    data.frame(
      Method = paste0("Method", i),
      RMSE = accuracy_list[[i]]["Test set", "RMSE"]#,
      #MAE = accuracy_list[[i]]["Test set", "MAE"]
    )
  }))

  accuracy_long <- tidyr::pivot_longer(accuracy_df, cols = c("RMSE"#, "MAE"
  ),
  names_to = "Metric", values_to = "Value")

  # Plot
  p <- ggplot(accuracy_long, aes(x = Method, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Forecast Accuracy Comparison",
         y = "Error", x = "Estimation Method") +
    theme_minimal()

  return(list(forecasts = forecast_list,
              accuracy = accuracy_list,
              accuracy_plot = p))
}
