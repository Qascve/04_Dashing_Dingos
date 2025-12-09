if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install with install.packages('here')")
}
library(here)

data_path <- here::here("groupwork2", "data", "KeyWestAnnualMeanTemperature.RData")
obj_name  <- load(data_path)
dat       <- get(obj_name)

temps <- dat$Temp

obs_cor <- cor(temps[-length(temps)], temps[-1])

set.seed(0)
n_perm   <- 10000

perm_cor <- numeric(n_perm)

for (i in 1:n_perm) {
  perm_temps   <- sample(temps)
  perm_cor[i]  <- cor(perm_temps[-length(perm_temps)],
                      perm_temps[-1])
}

p_value <- mean(perm_cor >= obs_cor)

cat("Observed correlation between successive years:",
    round(obs_cor, 3), "\n")

cat("Approximate permutation p-value (one-sided):",
    p_value, "\n")

if (p_value < 0.05) {
  cat("Conclusion: At the 0.05 level, successive years' ",
      "temperatures appear to be significantly positively ",
      "correlated.\n", sep = "")
} else {
  cat("Conclusion: At the 0.05 level, we do NOT find strong ",
      "evidence that successive years' temperatures are ",
      "positively correlated.\n", sep = "")
}

results_path <- here::here("groupwork2", "results", "TAutoCorr_results.txt")

out_lines <- c(
  "Autocorrelation in Key West annual mean temperature",
  paste("Observed correlation between successive years:",
        round(obs_cor, 3)),
  paste("Approximate permutation p-value (one-sided):",
        p_value)
)

if (p_value < 0.05) {
  out_lines <- c(out_lines,
                 "Conclusion (alpha = 0.05): successive years' temperatures are significantly positively correlated.")
} else {
  out_lines <- c(out_lines,
                 "Conclusion (alpha = 0.05): successive years' temperatures are not significantly positively correlated.")
}

writeLines(out_lines, con = results_path)

cat("Results written to", results_path, "\n")

png(filename = here::here("groupwork2", "results", "TAutoCorr_null_hist.png"))
hist(perm_cor,
     main = "Permutation distribution of correlation",
     xlab = "Correlation (permuted data)")
abline(v = obs_cor, lwd = 2)
dev.off()
