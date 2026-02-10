# This file contains plotting routines for the empirical study
# in Section 5.2
# Convert the matrix to a data frame
plot_ratio_boxplot <-function(ratiomatrix){
  Ratio_df <- as.data.frame(ratiomatrix)

  # Add column names "Lag" + number
  colnames(Ratio_df) <- 1:ncol(ratiomatrix)

  # Reshape the data frame from wide to long format
  Ratio_long <- suppressMessages(melt(Ratio_df, variable.name = "Lag", value.name = "MSEratio"))

  theme_update(text = element_text(size=30))
  # Create the boxplot
  boxplot <- ggplot(Ratio_long, aes(x = Lag, y = MSEratio)) +
    geom_boxplot(outlier.size = 0.75)+
    geom_hline(yintercept = 1, color = "red", linewidth = 0.75) +
    labs(x = "Lag", y = "MSE Ratio")+
    scale_x_discrete(breaks = levels(Ratio_long$Lag)[seq(1, length(levels(Ratio_long$Lag)), by = 2)])


  return(boxplot)

}

plot_DM_pvalues <- function(pvalues)
{
  # Convert the matrix to a data frame
  df <- as.data.frame(pvalues)

  # Add column names
  colnames(df) <- 1:ncol(pvalues)

  # Reshape the data frame from wide to long format
  df_long <- suppressMessages(melt(df, variable.name = "Lag", value.name = "pvalue"))

  theme_update(text = element_text(size=30))
  # Create the boxplot
  boxplot <- ggplot(df_long, aes(x = Lag, y = pvalue)) +
    geom_boxplot(outlier.size = 0.75)+
    geom_hline(yintercept = 0.05, color = "red", linewidth = 0.75) +
    labs(x = "Lag", y = "p-value")+
    scale_x_discrete(breaks = levels(df_long$Lag)[seq(1, length(levels(df_long$Lag)), by = 2)])


  return(boxplot)

}


##Plot boxplots of trawl function estimates
plot_trawl_boxplots <- function(trawlfct_matrix){

  df <- as.data.frame(trawlfct_matrix)

  # Add column names "Lag" + number
  colnames(df) <- 0:(ncol(trawlfct_matrix)-1)

  # Reshape the data frame from wide to long format
  df_long <- suppressMessages(melt(df, variable.name = "Lag", value.name = "Trawlfunction"))

  theme_update(text = element_text(size=30))
  # Create the boxplot
  boxplot <- ggplot(df_long, aes(x = Lag, y = Trawlfunction)) +
    geom_boxplot(outlier.size = 0.75)+
    labs(x = "Lag", y = "Trawl function")+
    scale_x_discrete(breaks = levels(df_long$Lag)[seq(1, length(levels(df_long$Lag)), by = 2)])

  return(boxplot)

}
