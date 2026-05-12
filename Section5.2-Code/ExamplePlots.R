library(ggplot2)


# Read in the data
ticker <- "A"
file_name <- paste0(ticker, "_", "5s_data.txt")
my_data <- as.matrix(read.table(file_name, sep=";"))


# Data
x <- as.numeric(my_data[1, ])
df_ts <- data.frame(
  Time = seq_along(x),
  Value = x
)

# Common theme
my_theme <- theme_gray(base_size = 30) +
  theme(
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

# -----------------------------
# 1. One path
# -----------------------------
p1 <- ggplot(df_ts, aes(x = Time, y = Value)) +
  geom_line(color = "blue", linewidth = 0.4) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000, 4000)
  ) +
  labs(x = "Time step", y = "Spread price") +
  my_theme

ggsave(
  "A_OnePath.eps",
  plot = p1,
  device = cairo_ps,
  width = 7.5,
  height = 5
)

p1
# -----------------------------
# 2. Histogram
# -----------------------------
p2 <- ggplot(df_ts, aes(x = Value)) +
  geom_histogram(
    breaks = seq(-0.5, 27.5, by = 1),
    fill = "blue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0, 27, by = 5),
    limits = c(0, 27)
  ) +
  labs(x = "Value", y = "Count") +
  my_theme

ggsave(
  "A_Hist.eps",
  plot = p2,
  device = cairo_ps,
  width = 6,
  height = 5
)

p2

# -----------------------------
# 3. ACF
# -----------------------------
acf_obj <- acf(x, lag.max = 90, plot = FALSE)

df_acf <- data.frame(
  Lag = 1:90,
  ACF = as.numeric(acf_obj$acf[-1])
)

ci <- 1.96 / sqrt(length(x))

p3 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
  geom_col(fill = "blue", width = 0.8) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(
    yintercept = c(-ci, ci),
    color = "black",
    linetype = "dashed"
  ) +
  labs(x = "Lag", y = "ACF") +
  my_theme

ggsave(
  "A_ACF.eps",
  plot = p3,
  device = cairo_ps,
  width = 6,
  height = 5
)

p3
