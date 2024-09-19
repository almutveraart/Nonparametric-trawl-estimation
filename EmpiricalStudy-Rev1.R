# Load all packages
library(ggplot2)
library(forecast)
library(reshape2)
library(ambit)

# Create a file to write the results
file_connection <- file("output.txt", open = "wt")


# Increase font size in labels in ggplot
theme_update(text = element_text(size=30)) 


# List all the tickers used in the empirical study
all_tickers <- c("A","DFS","WAT","WM")



# Set the sample grid width to 5s=5/60 min
my_Delta <- 5/60 

for(my_ticker in all_tickers){
  ticker = my_ticker
  print(paste("Ticker:", ticker))
  
  # Create the file name to read in 
  file_name <- paste0(ticker, "_", "5s_data.txt")
  
  # Read in the data
  my_data <- as.matrix(read.table(file_name, sep=";"))
  
  ndays <- dim(my_data)[1]
  nobs <- dim(my_data)[2] #3961
  
  my_lag <- 20
  
  all_trawlfcts <- matrix(nrow=ndays, ncol=my_lag)
  
  for(i in 1:ndays){
    # Fitting the trawl function
    all_trawlfcts[i,] <- nonpar_trawlest(my_data[i,], Delta=my_Delta, lag=my_lag)$a_hat
  }
  
  # Convert the matrix to a data frame
  all_trawlfcts_df <- as.data.frame(all_trawlfcts)
  
  # Add column names "Lag" + number
  colnames(all_trawlfcts_df) <- 0:(ncol(all_trawlfcts)-1)
  
  # Reshape the data frame from wide to long format
  all_trawlfcts_long <- melt(all_trawlfcts_df, variable.name = "Lag",value.name = "Trawlfunction")
  
  # Create the boxplot
  boxplot <- ggplot(all_trawlfcts_long, aes(x = Lag, y = Trawlfunction)) +
    geom_boxplot(outlier.size = 0.75)+
    labs(x = "Lag", y = "Trawl function")+
    scale_x_discrete(breaks = levels(all_trawlfcts_long$Lag)[seq(1, length(levels(all_trawlfcts_long$Lag)), by = 2)]) 
  
  boxplot
  
  # Save the boxplot
  plot_file_name <- paste0(ticker, "_", "TrawlFct_Boxplot.eps")
  ggsave(plot_file_name, plot = boxplot, device = "eps", width = 20, height = 20, units = "cm")
  
  # Forecasting:
  n <- nobs #3961
  is_length <- floor(80/100*n) #in sample length
  
  h_vector <- 1:20 
  hrange <- length(h_vector)
  noos <- n - is_length-max(h_vector)
  
  # Estimate the trawl function for each file and store it
  # in a matrix
  
  my_lag <- 200
  
  all_trawlfcts <- matrix(nrow=ndays, ncol=my_lag)
  all_means <- numeric(ndays)
  all_lebA_fc <- numeric(ndays)
  for(i in 1:ndays){
    
    # Print the index for every 10th iteration in the loop
    if (i %% 10 == 0) {
      print(paste("LebA calculation - Iteration:", i, ", Ticker:", ticker))
    }
    
    # Fitting the trawl function
    all_trawlfcts[i,] <- nonpar_trawlest(my_data[i,1:is_length], Delta=my_Delta, lag=my_lag)$a_hat
    
    # Estimate the Lebesgue measure of the trawl set 
    all_lebA_fc[i] <- LebA_est(my_data[i,1:is_length], my_Delta)
    
    # Computing the training data mean
    all_means[i] <- mean(my_data[i,1:is_length])
  }
  
  all_Weight_Intersection <- matrix(0, nrow=ndays, ncol=hrange)
  all_Weight_SetDifference <- matrix(0, nrow=ndays, ncol=hrange)
  
  for(i in 1:ndays){
    
    # Print the index for every 10th iteration in the loop
    if (i %% 10 == 0) {
      print(paste("Slice calculation - Iteration:", i, ", Ticker:", ticker))
    }
    
    data <- my_data[i,1:is_length]
    for(h in 1:hrange){ 
      my_h <- h_vector[h]
      
      slices <- LebA_slice_est(data, my_Delta, my_h*my_Delta)
      
      #Estimate Leb(A intersection A_h)
      lebAintersection_fc <-slices$LebAintersection 
      
      #Estimate Leb(A \ A_h)
      lebAsetdiff_fc <-slices$LebAsetdifference 
      
      #Components from trawl prediction formula:
      all_Weight_Intersection[i, h]<-lebAintersection_fc/all_lebA_fc[i]
      all_Weight_SetDifference[i,h]<-lebAsetdiff_fc/all_lebA_fc[i]
      
      
    }
  } # End for loop for computing forecasting components
  
  # Next compute the forecast
  # Create forecast matrices:
  CondMean <- array(0, dim = c(ndays, noos, hrange))
  NaiveFC <- array(0, dim = c(ndays, noos, hrange))
  ActualValues <- array(0, dim = c(ndays, noos, hrange))
  
  for(i in 1:ndays){
    
    x <- my_data[i,]
    
    for(tau in 1:noos){
      
      # Compute the actual value used for making the forecast
      actualvalue <- x[(is_length+tau-1)]
      
      for(h in 1:hrange){ 
        
        # ActualValue
        ActualValues[i,tau,h]<- x[(is_length+tau-1+my_h)]
        
        # Compute forecast(s)
        # naive
        NaiveFC[i,tau,h] <- actualvalue
        
        # Trawl prediction formula:
        CondMean[i,tau,h] <- actualvalue*all_Weight_Intersection[i,h]+all_means[i]*all_Weight_SetDifference[i, h]
        
      } #end for h
      
    }#end for tau
    
  }#end i for days
  
  # Computing the error measures
  
  CondMean_MSE <- matrix(0, nrow= ndays, ncol=hrange)
  Naive_MSE <- matrix(0, nrow= ndays, ncol=hrange)
  Ratio1 <- matrix(0, nrow= ndays, ncol=hrange)
  for(i in 1:ndays){
    
    for(h in 1:hrange){
      
      CondMean_MSE[i,h] <-my_mse(CondMean[i,,h], ActualValues[i,,h])
      Naive_MSE[i,h] <-my_mse(NaiveFC[i,,h], ActualValues[i,,h])
      Ratio1[i,h] <-CondMean_MSE[i,h]/Naive_MSE[i,h] 
    }
  }
  
  
  
  ############################
  # Use the ACF to compute the slices
  all_Weight_Intersection_acf <- matrix(0, nrow=ndays, ncol=hrange)
  all_Weight_SetDifference_acf <- matrix(0, nrow=ndays, ncol=hrange)
  
  CondMean_acf <- matrix(0, nrow=ndays, ncol=hrange)
  
  
  for(i in 1:ndays){
    data <- my_data[i,1:is_length]
    insampleacf <- acf(data, lag= (max(h_vector)+10),plot=FALSE)
    
    for(h in 1:hrange){ # forecasting horizon loop in 1:hrange
      my_h <- h_vector[h]
      
      my_acf <- insampleacf$acf[my_h+1]
      
      # Components from trawl prediction formula:
      all_Weight_Intersection_acf[i, h]<-my_acf
      all_Weight_SetDifference[i,h]<-1-my_acf
      
    }
  }
  
  ###
  # Create forecast matrices:
  CondMean_acf <- array(0, dim = c(ndays, noos, hrange))
  
  for(i in 1:ndays){
    
    x <- my_data[i,]
    
    for(tau in 1:noos){
      
      # Compute the actual value used for making the forecast
      actualvalue <- x[(is_length+tau-1)]
      
      for(h in 1:hrange){ 
        
        # Trawl prediction formula:
        CondMean_acf[i,tau,h] <- actualvalue*all_Weight_Intersection_acf[i,h]+all_means[i]*all_Weight_SetDifference_acf[i, h]
        
      } #end for h
      
    }#end for tau
    
  }#end i for days
  
  # Computing the error measures
  
  CondMean_acf_MSE <- matrix(0, nrow= ndays, ncol=hrange)
  Ratio2 <- matrix(0, nrow= ndays, ncol=hrange)
  for(i in 1:ndays){
    
    for(h in 1:hrange){
      
      CondMean_acf_MSE[i,h] <-my_mse(CondMean_acf[i,,h], ActualValues[i,,h])
      
      Ratio2[i,h] <-CondMean_MSE[i,h]/CondMean_acf_MSE[i,h] 
    }
  }
  
  
  #############
  # Create boxplots
  # Convert the matrix to a data frame
  Ratio1_df <- as.data.frame(Ratio1)
  
  # Add column names "Lag" + number
  colnames(Ratio1_df) <- 1:ncol(Ratio1)
  
  # Reshape the data frame from wide to long format
  Ratio1_long <- melt(Ratio1_df, variable.name = "Lag",value.name = "MSEratio")
  
  # Create the boxplot
  boxplot <- ggplot(Ratio1_long, aes(x = Lag, y = MSEratio)) +
    geom_boxplot(outlier.size = 0.75)+
    geom_hline(yintercept = 1, color = "red", linewidth = 0.75) +
    labs(x = "Lag", y = "MSE Ratio")+
    scale_x_discrete(breaks = levels(Ratio1_long$Lag)[seq(1, length(levels(Ratio1_long$Lag)), by = 2)]) 
  
  
  boxplot
  # Save the boxplot
  plot_file_name <- paste0(ticker, "_", "Ratio1.eps")
  ggsave(plot_file_name, plot = boxplot, device = "eps", width = 20, height = 20, units = "cm")
  
  #################
  # Convert the matrix to a data frame
  Ratio2_df <- as.data.frame(Ratio2)
  
  # Add column names "Lag" + number
  colnames(Ratio2_df) <- 1:ncol(Ratio2)
  
  # Reshape the data frame from wide to long format
  Ratio2_long <- melt(Ratio2_df, variable.name = "Lag",value.name = "MSEratio")
  
  # Create the boxplot
  boxplot <- ggplot(Ratio2_long, aes(x = Lag, y = MSEratio)) +
    geom_boxplot(outlier.size = 0.75)+
    geom_hline(yintercept = 1, color = "red", linewidth = 0.75) +
    labs(x = "Lag", y = "MSE Ratio")+
    scale_x_discrete(breaks = levels(Ratio2_long$Lag)[seq(1, length(levels(Ratio2_long$Lag)), by = 2)]) 
  
  
  boxplot
  # Save the boxplot
  plot_file_name <- paste0(ticker, "_", "Ratio2.eps")
  ggsave(plot_file_name, plot = boxplot, device = "eps", width = 20, height = 20, units = "cm")
  
  
  # Diebold-Mariano tests
  # Create error matrices:
  CondMean_Error <- CondMean-ActualValues 
  # Naive_Error <- NaiveFC - ActualValues
  CondMean_acf_Error <- CondMean_acf - ActualValues
  
  # Nonparametric trawl versus naive
  CM_ACF_g <-matrix(0, nrow=ndays, ncol=hrange)
  for(i in 1:ndays){
    for(h in 1:hrange){
      CM_ACF_g[i,h]<-dm.test(CondMean_acf_Error[i,,h], CondMean_Error[i,,h] , alternative="greater", h=h, power =2)$p.value
    }
  }
  
  # Save the p-values
  as.vector(CM_ACF_g)
  all_tickers <- c("A","DFS","WAT","WM")
  if(ticker =="A"){
    A_pvalues <- as.vector(CM_ACF_g)
  }
  else if(ticker =="DFS"){
    DFS_pvalues <- as.vector(CM_ACF_g)
  }
  else if(ticker =="WAT"){
    WAT_pvalues <- as.vector(CM_ACF_g)
  }
  else if(ticker =="WM"){
    WM_pvalues <- as.vector(CM_ACF_g)
  }
  else{
    print("Other ticker")
  }
  
  #################
  # Convert the matrix to a data frame
  CM_ACF_g_df <- as.data.frame(CM_ACF_g)
  
  # Add column names 
  colnames(CM_ACF_g_df) <- 1:ncol(CM_ACF_g)
  
  # Reshape the data frame from wide to long format
  CM_ACF_g_long <- melt(CM_ACF_g_df, variable.name = "Lag",value.name = "pvalue")
  
  # Create the boxplot
  boxplot <- ggplot(CM_ACF_g_long, aes(x = Lag, y = pvalue)) +
    geom_boxplot(outlier.size = 0.75)+
    geom_hline(yintercept = 0.05, color = "red", linewidth = 0.75) +
    labs(x = "Lag", y = "p-value")+
    scale_x_discrete(breaks = levels(CM_ACF_g_long$Lag)[seq(1, length(levels(CM_ACF_g_long$Lag)), by = 2)]) 
  
  
  boxplot
  # Save the boxplot
  plot_file_name <- paste0(ticker, "_", "DM.eps")
  ggsave(plot_file_name, plot = boxplot, device = "eps", width = 20, height = 20, units = "cm")
  
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
  ggsave(plot_adj_pvalues_file, plot = plot_adj_pvalues, device = "eps", width = 20, height = 20, units = "cm")
  
} # end ticker for-loop


# Close the file connection and save results
close(file_connection)

all_pvalues <-c(A_pvalues, DFS_pvalues,WAT_pvalues,WM_pvalues)

# Compute adjusted p-values using "BH"
fdrs_all <- p.adjust(all_pvalues, method="BH")

# Count the number of values <= 0.05
count_all <- sum(fdrs_all <= 0.05)

# Calculate the percentage
percentage_all <- (count_all / length(fdrs_all)) * 100
percentage_all