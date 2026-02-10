# Load all packages
library(ggplot2)
library(forecast)
library(reshape2)
library(tidyr)
library(dplyr)
library(ambit)

source("AdditionalRFunctions.R")
# Create a file to write the results
file_connection <- file("output.txt", open = "wt")
on.exit(close(file_connection), add = TRUE)

# Increase font size in labels in ggplot
theme_update(text = element_text(size=30))


# List all the tickers used in the empirical study
all_tickers <- c("A","DFS","WAT","WM")



# Set the sample grid width to 5s=5/60 min
my_Delta <- 5/60

# For testing - uncomment the line below to limit days for all tickers
# testing_days_override <- 10

for(my_ticker in all_tickers){
  ticker = my_ticker
  print(paste("Ticker:", ticker))

  # Create the file name to read in
  file_name <- paste0(ticker, "_", "5s_data.txt")

  # Read in the data
  my_data <- as.matrix(read.table(file_name, sep=";"))

  ndays <- dim(my_data)[1]
  nobs <- dim(my_data)[2] #3961

  # Use all days by default, unless a global override was set for testing
  # Check if a global testing override exists (set before the ticker loop)
  if(exists("testing_days_override")) {
    used_days <- testing_days_override
  } else {
    used_days <- ndays  # Use actual number of days for this ticker
  }

  # Forecasting:
  n <- nobs #3961
  is_length <- floor(80/100*n) #in sample length

  h_vector <- 1:10
  hrange <- length(h_vector)
  noos <- n - is_length-max(h_vector)
  stopifnot(noos > 0L)

  # Estimate the trawl function for each file and store it
  # in a matrix

  my_lag <- 200  # For plotting/visualisation only
  my_lag_full <- is_length - 1  # For forecasting (full trawl returns is_length-1 when lag=is_length)

  all_trawlfcts <- matrix(nrow=used_days, ncol=my_lag)  # For plotting
  all_trawlfcts_full <- matrix(nrow=used_days, ncol=my_lag_full)  # For forecasting
  all_means <- numeric(used_days)
  all_lebA_fc <- numeric(used_days)

  #For nonparametric trawl
  all_Weight_Intersection <- matrix(0, nrow=used_days, ncol=hrange)
  all_Weight_SetDifference <- matrix(0, nrow=used_days, ncol=hrange)

  #For exponential fit
  all_Weight_Intersection_Exp <- matrix(0, nrow=used_days, ncol=hrange)
  all_Weight_SetDifference_Exp <- matrix(0, nrow=used_days, ncol=hrange)

  #For LM trawl fit
  all_Weight_Intersection_LM <- matrix(0, nrow=used_days, ncol=hrange)
  all_Weight_SetDifference_LM <- matrix(0, nrow=used_days, ncol=hrange)

  # For ACF-based slice estimation
  all_Weight_Intersection_acf <- matrix(0, nrow=used_days, ncol=hrange)
  all_Weight_SetDifference_acf <- matrix(0, nrow=used_days, ncol=hrange)

  # Create forecast matrices:
  CondMean <- array(0, dim = c(used_days, noos, hrange))
  CondMean_Exp <- array(0, dim = c(used_days, noos, hrange))
  CondMean_LM <- array(0, dim = c(used_days, noos, hrange))
  CondMean_acf <- array(0, dim = c(used_days, noos, hrange))

  NaiveFC <- array(0, dim = c(used_days, noos, hrange))
  ActualValues <- array(0, dim = c(used_days, noos, hrange))

  #Error measures
  CondMean_MSE <- matrix(0, nrow= used_days, ncol=hrange)
  CondMean_Exp_MSE <- matrix(0, nrow= used_days, ncol=hrange)
  CondMean_LM_MSE <- matrix(0, nrow= used_days, ncol=hrange)
  CondMean_acf_MSE <- matrix(0, nrow= used_days, ncol=hrange)

  Naive_MSE <- matrix(0, nrow= used_days, ncol=hrange)

  Ratio_Naive <- matrix(0, nrow= used_days, ncol=hrange)
  Ratio_Exp <- matrix(0, nrow= used_days, ncol=hrange)
  Ratio_LM <- matrix(0, nrow= used_days, ncol=hrange)
  Ratio_acf <- matrix(0, nrow= used_days, ncol=hrange)



  for(i in 1:used_days){#ndays){

    # Print the index for every 10th iteration in the loop
    if (i %% 10 == 0) {
      print(paste("Iteration:", i, ", Ticker:", ticker))
    }

    # Fitting the trawl function for plotting (lag=200)
    all_trawlfcts[i,] <- nonpar_trawlest(my_data[i,1:is_length], Delta=my_Delta, lag=my_lag)$a_hat

    # Fitting the full trawl function for forecasting
    # Note: when lag=is_length, nonpar_trawlest returns a vector of length is_length-1
    all_trawlfcts_full[i,] <- nonpar_trawlest(my_data[i,1:is_length], Delta=my_Delta, lag=is_length)$a_hat

    # Computing the training data mean
    all_means[i] <- mean(my_data[i,1:is_length])

    #Fit parametric exp trawl
    Expfit <- trawl::fit_Exptrawl(my_data[i,1:is_length], my_Delta)
    Exp_lambda <- Expfit$lambda

    #Fit parametric LM trawl
    LMfit <-trawl::fit_LMtrawl(my_data[i,1:is_length], my_Delta, GMMlag = 5)
    LM_alpha <-LMfit$alpha
    LM_H <-LMfit$H

    #Compute acf
    insampleacf <- acf(my_data[i,1:is_length], lag.max= (max(h_vector)+10),plot=FALSE)

    # Compute LebA from the full trawl function
    all_lebA_fc[i] <- sum(all_trawlfcts_full[i,]) * my_Delta

    #data <- my_data[i,1:is_length]
    for(h in 1:hrange){
      my_h <- h_vector[h]

      # Use pre-estimated full trawl function for slices
      slices <- LebA_slice_est_approx(all_trawlfcts_full[i,], my_Delta, my_h*my_Delta)

      # Normalise using the same LebA
      all_Weight_Intersection[i, h] <- slices$LebAintersection / slices$LebA
      all_Weight_SetDifference[i, h] <- slices$LebAsetdifference / slices$LebA

      #Compute ratios for parametric exp model:
      #Estimate Leb(A intersection A_h)/Leb(A)
      lebAintersection_ratio_exp <-exp(-my_h*my_Delta*Exp_lambda)

      #Estimate Leb(A \ A_h)/Leb(A)
      lebAsetdiff_ratio_exp <-1-lebAintersection_ratio_exp

      all_Weight_Intersection_Exp[i, h]<-lebAintersection_ratio_exp
      all_Weight_SetDifference_Exp[i,h]<-lebAsetdiff_ratio_exp

      #Compute ratios for parametric LM model:
      #Estimate Leb(A intersection A_h)/Leb(A)
      lebAintersection_ratio_LM <-(1+my_h*my_Delta/LM_alpha)^(1-LM_H)

      #Estimate Leb(A \ A_h)/Leb(A)
      lebAsetdiff_ratio_LM <-1-lebAintersection_ratio_LM

      all_Weight_Intersection_LM[i, h]<-lebAintersection_ratio_LM
      all_Weight_SetDifference_LM[i,h]<-lebAsetdiff_ratio_LM

      #Compute ratios acf based
      my_acf <- insampleacf$acf[my_h+1]

      # Components from trawl prediction formula:
      all_Weight_Intersection_acf[i, h]<-my_acf
      all_Weight_SetDifference_acf[i,h]<-1-my_acf


    }#end h loop (forecast horizon)
  } # End i loop (days in sample)

  # Next compute the forecast


  for(i in 1:used_days){#ndays){

    x <- my_data[i,]

    for(tau in 1:noos){

      # Compute the actual value used for making the forecast
      actualvalue <- x[(is_length+tau-1)]

      for(h in 1:hrange){
        my_h <- h_vector[h]

        # ActualValue
        ActualValues[i,tau,h]<- x[(is_length+tau-1+my_h)]

        # Compute forecast(s)
        # naive
        NaiveFC[i,tau,h] <- actualvalue

        # Trawl prediction formula:
        CondMean[i,tau,h] <- actualvalue*all_Weight_Intersection[i,h]
        +all_means[i]*all_Weight_SetDifference[i, h]
        CondMean_Exp[i,tau,h] <- actualvalue*all_Weight_Intersection_Exp[i,h]
        +all_means[i]*all_Weight_SetDifference_Exp[i, h]
        CondMean_LM[i,tau,h] <- actualvalue*all_Weight_Intersection_LM[i,h]
        +all_means[i]*all_Weight_SetDifference_LM[i, h]
        CondMean_acf[i,tau,h] <- actualvalue*all_Weight_Intersection_acf[i,h]
        +all_means[i]*all_Weight_SetDifference_acf[i, h]
      } #end for h

    }#end for tau

  }#end i for days



  # Computing the error measures


  for(i in 1:used_days){#ndays){

    for(h in 1:hrange){

      CondMean_MSE[i,h] <-my_mse(CondMean[i,,h], ActualValues[i,,h])
      Naive_MSE[i,h] <-my_mse(NaiveFC[i,,h], ActualValues[i,,h])
      Ratio_Naive[i,h] <-CondMean_MSE[i,h]/Naive_MSE[i,h]

      CondMean_Exp_MSE[i,h] <-my_mse(CondMean_Exp[i,,h], ActualValues[i,,h])
      Ratio_Exp[i,h] <-CondMean_MSE[i,h]/CondMean_Exp_MSE[i,h]

      CondMean_LM_MSE[i,h] <-my_mse(CondMean_LM[i,,h], ActualValues[i,,h])
      Ratio_LM[i,h] <-CondMean_MSE[i,h]/CondMean_LM_MSE[i,h]

      CondMean_acf_MSE[i,h] <-my_mse(CondMean_acf[i,,h], ActualValues[i,,h])
      Ratio_acf[i,h] <-CondMean_MSE[i,h]/CondMean_acf_MSE[i,h]
    }
  }





  #############
  # Create boxplots
  p <- plot_ratio_boxplot(Ratio_Naive)
  plot_file_name <- paste0(ticker, "_", "Ratio_Naive.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")

  #################
  p <- plot_ratio_boxplot(Ratio_Exp)
  plot_file_name <- paste0(ticker, "_", "Ratio_Exp.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")

  #################
  p<-plot_ratio_boxplot(Ratio_LM)
  plot_file_name <- paste0(ticker, "_", "Ratio_LM.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")

  #################
  p <- plot_ratio_boxplot(Ratio_acf)
  plot_file_name <- paste0(ticker, "_", "Ratio_acf.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")

  ###############################

  # Diebold-Mariano tests
  # Create error matrices:
  CondMean_Error <- CondMean-ActualValues
  # Naive_Error <- NaiveFC - ActualValues
  CondMean_acf_Error <- CondMean_acf - ActualValues
  CondMean_Exp_Error <- CondMean_Exp - ActualValues
  CondMean_LM_Error <- CondMean_LM - ActualValues

  # Nonparametric trawl versus acf
  CM_ACF_g <-matrix(0, nrow=used_days, ncol=hrange)
  CM_Exp_g <-matrix(0, nrow=used_days, ncol=hrange)
  CM_LM_g <-matrix(0, nrow=used_days, ncol=hrange)
  for(i in 1:used_days){#ndays){
    for(h in 1:hrange){
      my_h <- h_vector[h]
      CM_ACF_g[i,h]<-dm.test(CondMean_acf_Error[i,,h], CondMean_Error[i,,h] , alternative="greater", h=my_h, power =2)$p.value
      CM_Exp_g[i,h]<-dm.test(CondMean_Exp_Error[i,,h], CondMean_Error[i,,h] , alternative="greater", h=my_h, power =2)$p.value
      CM_LM_g[i,h]<-dm.test(CondMean_LM_Error[i,,h], CondMean_Error[i,,h] , alternative="greater", h=my_h, power =2)$p.value

    }
  }

  # Save the p-values
  if(ticker =="A"){
    A_pvalues <- as.vector(CM_ACF_g)
    A_pvalues_Exp <- as.vector(CM_Exp_g)
    A_pvalues_LM <- as.vector(CM_LM_g)
  }
  else if(ticker =="DFS"){
    DFS_pvalues <- as.vector(CM_ACF_g)
    DFS_pvalues_Exp <- as.vector(CM_Exp_g)
    DFS_pvalues_LM <- as.vector(CM_LM_g)
  }
  else if(ticker =="WAT"){
    WAT_pvalues <- as.vector(CM_ACF_g)
    WAT_pvalues_Exp <- as.vector(CM_Exp_g)
    WAT_pvalues_LM <- as.vector(CM_LM_g)
  }
  else if(ticker =="WM"){
    WM_pvalues <- as.vector(CM_ACF_g)
    WM_pvalues_Exp <- as.vector(CM_Exp_g)
    WM_pvalues_LM <- as.vector(CM_LM_g)
  }
  else{
    print("Other ticker")
  }

  #################

  p <- plot_DM_pvalues(CM_ACF_g)
  plot_file_name <- paste0(ticker, "_", "DM.eps")
  ggsave(plot_file_name,  plot=p, device = "eps", width = 20, height = 20, units = "cm")

  # Compute adjusted p-values using "BH"
  fdrs <- p.adjust(CM_ACF_g, method="BH")

  # Count the number of values <= 0.05
  count <- sum(fdrs <= 0.05)

  # Calculate the percentage
  percentage <- (count / length(fdrs)) * 100
  percentage


  print(percentage)

  # Write the result to the file
  cat("The 0.05 percentage for ticker", ticker, "is", percentage, "\n", file = file_connection)



  # Plot the histogram of the adjusted p-values
  data <- fdrs
  plot_adj_pvalues <- ggplot(data.frame(x = data), aes(x = x)) +
    geom_histogram(bins = 30, color = "black", fill = "lightblue") +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
    labs(x = "Values", y = "Frequency of adjusted p-values")

  plot_adj_pvalues
  # Save the boxplot
  plot_adj_pvalues_file <- paste0(ticker, "_", "AdjPvalues.eps")
  ggsave(plot_adj_pvalues_file, plot=plot_adj_pvalues, device = "eps", width = 20, height = 20, units = "cm")


  #Plot the trawl function (in sample)
  p <- plot_trawl_boxplots(all_trawlfcts[,1:10])

  plot_file_name <- paste0(ticker, "_", "TrawlFct_Boxplot.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")

  ######################
  #New plots: Compare the weights based on the four methods
  # Example: put all weight matrices into a list with names
  weights_list <- list(
    Nonpar.       = all_Weight_Intersection,
    ACF           = all_Weight_Intersection_acf,
    Exp           = all_Weight_Intersection_Exp,
    LM            = all_Weight_Intersection_LM
  )

  # Convert each to long format
  df_long <- bind_rows(
    lapply(names(weights_list), function(method) {
      mat <- weights_list[[method]]
      data.frame(
        day      = rep(1:nrow(mat), times = ncol(mat)),
        horizon  = rep(h_vector, each = nrow(mat)),
        method   = method,
        weight   = as.vector(mat)
      )
    }),
    .id = NULL
  )

  # Make method a factor for consistent ordering
  df_long$method <- factor(df_long$method,
                           levels = c("Nonpar.", "ACF","Exp", "LM"))

  # Boxplot comparison
  p<-ggplot(df_long, aes(x = factor(horizon), y = weight, fill = method)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
    labs(x = "Forecast horizon (k)", y = "Set intersection weight",
         title = "", fill="Method")
  p
  plot_file_name <- paste0(ticker, "_", "WeightComp_Intersec.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")

  # Boxplot comparison for set difference:
  weights_list <- list(
    Nonpar.       = all_Weight_SetDifference,
    ACF           = all_Weight_SetDifference_acf,
    Exp           = all_Weight_SetDifference_Exp,
    LM            = all_Weight_SetDifference_LM
  )

  # Convert each to long format
  df_long <- bind_rows(
    lapply(names(weights_list), function(method) {
      mat <- weights_list[[method]]
      data.frame(
        day      = rep(1:nrow(mat), times = ncol(mat)),
        horizon  = rep(h_vector, each = nrow(mat)),
        method   = method,
        weight   = as.vector(mat)
      )
    }),
    .id = NULL
  )

  # Make method a factor for consistent ordering
  df_long$method <- factor(df_long$method,
                           levels = c("Nonpar.", "ACF", "Exp", "LM"))

  p<-ggplot(df_long, aes(x = factor(horizon), y = weight, fill = method)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
    labs(x = "Forecast horizon lag (k)", y = "Set difference weight",
         fill="Method",
         title = "")
  p
  plot_file_name <- paste0(ticker, "_", "WeightComp_Setdiff.eps")
  ggsave(plot_file_name, plot=p, device = "eps", width = 20, height = 20, units = "cm")



} # end ticker for-loop


all_pvalues <-c(A_pvalues, DFS_pvalues,WAT_pvalues,WM_pvalues)

# Compute adjusted p-values using "BH"
fdrs_all <- p.adjust(all_pvalues, method="BH")

# Count the number of values <= 0.05
count_all <- sum(fdrs_all <= 0.05)

# Calculate the percentage
percentage_all <- (count_all / length(fdrs_all)) * 100
percentage_all
