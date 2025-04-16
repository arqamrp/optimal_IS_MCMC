# Configuration
num_chains <- 20          # Number of independent chains per sampler
checkpoints <- c(1e2, 1e3, 1e4, 1e5)  # Sample sizes to evaluate
d <- ncol(X)             # Dimension of parameters

# Initialize storage arrays
results <- list(
  PMALA = list(
    estimates = array(NA, dim = c(num_chains, d, length(checkpoints))),
    times = numeric(num_chains)
  ),
  MNMALA = list(
    estimates = array(NA, dim = c(num_chains, d, length(checkpoints))),
    times = numeric(num_chains)
  ),
  MALA = list(
    estimates = array(NA, dim = c(num_chains, d, length(checkpoints))),
    times = numeric(num_chains)
  )
)

# Sampler configurations
config <- list(
  PMALA = list(
    step_size = 2.7e-3,
    lambda = 0.001,
    max_samples = max(checkpoints)
  ),
  MNMALA = list(
    step_size = 5.8e-4,
    max_samples = max(checkpoints)
  ),
  MALA = list(
    step_size = 5.7e-4,
    max_samples = max(checkpoints)
  )
)

# Run chains for PMALA
for(i in 1:num_chains) {
  cat("\nPMALA Chain", i, "/", num_chains)
  start_time <- Sys.time()
  
  # Run full chain
  result <- poisson_PMALA_sampler(
    X, y, 
    init = rep(0, d),
    num_samples = config$PMALA$max_samples,
    step_size = config$PMALA$step_size,
    lambda = config$PMALA$lambda,
    mu_hat = mle_coef_manual,
    NR = TRUE
  )
  
  # Store timing
  results$PMALA$times[i] <- difftime(Sys.time(), start_time, units = "secs")
  
  # Process checkpoints
  for(k in seq_along(checkpoints)) {
    n <- checkpoints[k]
    idx <- 1:n
    
    # Calculate weighted estimate
    log_post <- apply(result$chain[idx,], 1, log_post_poisson_simple)
    weights <- exp(log_post - result$log_opt_density[idx])
    results$PMALA$estimates[i,,k] <- colSums(result$chain[idx,] * weights) / sum(weights)
  }
}

# Run chains for MNMALA
for(i in 1:num_chains) {
  cat("\nMNMALA Chain", i, "/", num_chains)
  start_time <- Sys.time()
  
  # Run full chain
  chain <- MALA_sampler(
    log_opt_poisson_simple,
    grad_log_opt_poisson_simple,
    num_samples = config$MNMALA$max_samples,
    init = rep(0, d),
    step_size = config$MNMALA$step_size
  )
  
  # Store timing
  results$MNMALA$times[i] <- difftime(Sys.time(), start_time, units = "secs")
  
  # Process checkpoints
  for(k in seq_along(checkpoints)) {
    n <- checkpoints[k]
    idx <- 1:n
    
    # Calculate weighted estimate
    log_post <- apply(chain[idx,], 1, log_post_poisson_simple)
    log_opt <- apply(chain[idx,], 1, log_opt_poisson_simple)
    weights <- exp(log_post - log_opt)
    results$MNMALA$estimates[i,,k] <- colSums(chain[idx,] * weights) / sum(weights)
  }
}

# Run chains for MALA
for(i in 1:num_chains) {
  cat("\nMALA Chain", i, "/", num_chains)
  start_time <- Sys.time()
  
  # Run full chain
  chain <- MALA_sampler(
    log_post_poisson_simple,
    grad_log_post_poisson_simple,
    num_samples = config$MALA$max_samples,
    init = rep(0, d),
    step_size = config$MALA$step_size
  )
  
  # Store timing
  results$MALA$times[i] <- difftime(Sys.time(), start_time, units = "secs")
  
  # Process checkpoints
  for(k in seq_along(checkpoints)) {
    n <- checkpoints[k]
    results$MALA$estimates[i,,k] <- colMeans(chain[1:n,])
  }
}

# Calculate variances and timing statistics
analyze_results <- function(estimates) {
  apply(estimates, c(2,3), var)
}

variance_results <- list(
  PMALA = analyze_results(results$PMALA$estimates),
  MNMALA = analyze_results(results$MNMALA$estimates),
  MALA = analyze_results(results$MALA$estimates)
)

time_stats <- list(
  PMALA = c(
    mean = mean(results$PMALA$times),
    sd = sd(results$PMALA$times)
  ),
  MNMALA = c(
    mean = mean(results$MNMALA$times),
    sd = sd(results$MNMALA$times)
  ),
  MALA = c(
    mean = mean(results$MALA$times),
    sd = sd(results$MALA$times)
  )
)

# Print variance results
print_variance_table <- function(variances, checkpoint_idx) {
  data.frame(
    Parameter = colnames(X),
    PMALA = variances$PMALA[,checkpoint_idx],
    MNMALA = variances$MNMALA[,checkpoint_idx],
    MALA = variances$MALA[,checkpoint_idx]
  )
}

cat("\nVariance at 1e3 samples:")
print(print_variance_table(variance_results, 1))

cat("\nVariance at 1e4 samples:")
print(print_variance_table(variance_results, 2))

cat("\nVariance at 1e5 samples:")
print(print_variance_table(variance_results, 3))

cat("\nVariance at 1e6 samples:")
print(print_variance_table(variance_results, 4))


# Print timing statistics
cat("\nAverage Computation Times (seconds):")
print(data.frame(
  PMALA = time_stats$PMALA["mean"],
  MNMALA = time_stats$MNMALA["mean"],
  MALA = time_stats$MALA["mean"]
))



library(ggplot2)

# Calculate summed variances across parameters
summed_variances <- list(
  PMALA = colSums(variance_results$PMALA),
  MNMALA = colSums(variance_results$MNMALA),
  MALA = colSums(variance_results$MALA)
)

# Create plotting data frame
plot_data <- data.frame(
  log_n = rep(log10(checkpoints), 3),
  sum_var = c(summed_variances$PMALA, summed_variances$MNMALA, summed_variances$MALA),
  sampler = rep(c("PMALA", "MNMALA", "MALA"), each = length(checkpoints))
)

# Create plot
ggplot(plot_data, aes(x = log_n, y = sum_var, color = sampler)) +
  geom_line(linewidth = 1.2) + 
  geom_point(size = 3) +
  scale_x_continuous(breaks = log10(checkpoints), 
                     labels = expression(10^3, 10^4, 10^5, 10^6)) +
  labs(x = "Number of Samples (log scale)",
       y = "Total Variance Across Parameters",
       title = "Convergence Diagnostics by Sampler",
       color = "Sampler Type") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))
