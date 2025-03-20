HMC_sampler <- function(log_density_fn, grad_log_fn, init, num_samples = 1e4, step_size = 0.1, L = 10, verbose = FALSE) {
  num_accepts <- 0
  d <- length(init)  # Number of dimensions
  
  chain <- matrix(NA, nrow = num_samples, ncol = d)
  chain[1, ] <- init
  
  for(i in 2:num_samples) {
    current <- chain[i-1, ]
    momentum <- rnorm(d)  # Sample initial momentum from N(0, I)
    
    # Leapfrog integration
    prop <- current
    prop_momentum <- momentum + 0.5 * step_size * grad_log_fn(prop)  # Half step in momentum
    
    for(l in 1:L) {
      prop <- prop + step_size * prop_momentum  # Full step in position
      if(l != L) prop_momentum <- prop_momentum + step_size * grad_log_fn(prop)  # Full step in momentum
    }
    
    prop_momentum <- prop_momentum + 0.5 * step_size * grad_log_fn(prop)  # Final half step in momentum
    prop_momentum <- -prop_momentum
    
    log_MH_ratio <- log_density_fn(prop) - log_density_fn(current) + 
      sum(momentum^2 - prop_momentum^2) / 2  # Kinetic energy term
    
    if(log(runif(1)) < log_MH_ratio) {
      chain[i, ] <- prop
      num_accepts <- num_accepts + 1
    } else {
      chain[i, ] <- current
    }
  }
  
  if(verbose) print(paste("HMC Acceptance rate is", num_accepts / num_samples))
  return(chain)
}

rwMH_sampler <- function(log_density_fn, init,  prop_var = 1, num_samples = 1e4, verbose = FALSE) {
  num_accepts <- 0
  d <- length(init)  
  
  chain <- matrix(NA, nrow = num_samples, ncol = d)
  chain[1, ] <- init
  
  for(i in 2:num_samples) {
    prop <- chain[i-1, ] + sqrt(prop_var) * rnorm(d)
    log_MH_ratio <- log_density_fn(prop) - log_density_fn(chain[i-1, ])
    
    
    if(log(runif(1)) < log_MH_ratio) {
      chain[i, ] <- prop
      num_accepts <- num_accepts + 1
    } else {
      chain[i, ] <- chain[i-1, ]
    }
  }
  
  if(verbose) print(paste("RWMH Acceptance rate is", num_accepts / num_samples))
  return(chain)
}

MALA_sampler <- function(log_density_fn, grad_log_fn, init, num_samples = 1e4, step_size = 0.1, verbose = FALSE) {
  num_accepts <- 0
  d <- length(init) 
  
  chain <- matrix(NA, nrow = num_samples, ncol = d)
  chain[1, ] <- init
  
  for(i in 2:num_samples) {
    current <- chain[i-1, ]
    grad <- grad_log_fn(current)
    prop <- current + 0.5 * step_size * grad + sqrt(step_size) * rnorm(d)
    
    grad_prop <- grad_log_fn(prop)  # Gradient at proposed position
    
    # Compute log transition probabilities
    log_q_forward <- sum(dnorm(prop, mean = current + 0.5 * step_size * grad, sd = sqrt(step_size), log = TRUE))
    log_q_backward <- sum(dnorm(current, mean = prop + 0.5 * step_size * grad_prop, sd = sqrt(step_size), log = TRUE))
    
    log_MH_ratio <- log_density_fn(prop) - log_density_fn(current) + log_q_backward - log_q_forward
    
    if(log(runif(1)) < log_MH_ratio) {
      chain[i, ] <- prop
      num_accepts <- num_accepts + 1
    } else {
      chain[i, ] <- current
    }
  }
  
  if(verbose) print(paste("MALA Acceptance rate is", num_accepts / num_samples))
  return(chain)
}

plot_traces <- function(samples, title, true_values = NULL) {
  df <- as.data.frame(samples)
  df$Iteration <- 1:nrow(df)
  df_long <- reshape2::melt(df, id = "Iteration")
  
  p <- ggplot(df_long, aes(x = Iteration, y = value, color = variable)) +
    geom_line() +
    theme_minimal() +
    ggtitle(title)
  
  # Add horizontal reference lines if true values are provided
  if (!is.null(true_values)) {
    for (i in seq_along(true_values)) {
      p <- p + geom_hline(yintercept = true_values[i], linetype = "dashed", color = "black")
    }
  }
  
  print(p)
}

moreau_norm <- function(vec, lambda = 0, eps = 1e-5){
  vec_norm <- sqrt(sum(vec^2))
  if( vec_norm <= lambda & lambda != 0){
    return(vec_norm**2/ (2 * lambda) + eps) 
  }
  else{
    return(vec_norm - lambda/2 + eps)
  }
}

moreau_norm(0)

grad_log_moreau_norm <- function(vec, lambda = 0, eps = 1e-5) {
  vec_norm <- sqrt(sum(vec^2))  # L2 norm
  
  if (vec_norm <= lambda) {
    mureau_value <- (vec_norm^2) / (2 * lambda) + eps
    gradient_mureau <- vec / lambda
  } else {
    mureau_value <- vec_norm - (lambda / 2) + eps
    gradient_mureau <- vec / vec_norm
  }
  gradient_log <- gradient_mureau / mureau_value
  return(gradient_log)
}

############ Toy example 1

toy1_simulate_data <- function(n, d) {
  set.seed(42)
  X <- cbind(1, matrix(rnorm(n * (d - 1)), n, d - 1))  # design matrix with intercept
  beta_true <- rnorm(d)  # randomly genrate true coeffs
  prob <- 1 / (1 + exp(-X %*% beta_true)) 
  y <- rbinom(n, 1, prob)  # generate binary response
  return(list(X = X, y = y, beta_true = beta_true))
}

d <- 12
data <- toy1_simulate_data(n = 200, d = d)
X <- data$X
y <- data$y
beta_true <- data$beta_true
print(beta_true)

toy1_log_post <-  function(beta, prior_var = 1){
  #log_prior <- sum( dnorm( beta, mean = 0, sd = sqrt(prior_var), log = TRUE) ) 
  log_prior <- sum( log(0.5* dnorm( beta, mean = 0, sd = sqrt(prior_var))+ 0.5*dnorm( beta, mean = 2, sd = sqrt(prior_var)) ))
  eta <- X %*% beta
  log_lik <- sum(y * eta - log(1 + exp(eta))) 
  return(log_prior+log_lik)
}

toy1_grad_log_post <- function(beta,  prior_var = 1){
  eta <- X %*% beta
  p <- 1 / (1 + exp(-eta))
  grad <- t(X) %*% (y - p) - beta / prior_var
  return(grad)
}

toy1_log_opt_SIS <- function(beta, i){
  #return(toy1_log_post(beta) + log(moreau_norm(beta - toy1_mean_estimates_post[i, ])) )
  return(toy1_log_post(beta) + log(moreau_norm(beta)) )
}

toy1_grad_log_opt_SIS <- function(beta, i){
  #return(toy1_grad_log_post(beta) +  grad_log_moreau_norm(beta- toy1_mean_estimates_post[i, ]))
  return(toy1_grad_log_post(beta) +  grad_log_moreau_norm(beta))
}

toy1_log_opt <- function(beta, i){
  return(toy1_log_post(beta) + log(moreau_norm(beta - toy1_mean_estimates_post[i, ])) )

}

toy1_grad_log_opt <- function(beta, i){
  return(toy1_grad_log_post(beta) +  grad_log_moreau_norm(beta- toy1_mean_estimates_post[i, ]))
}



## prop_var = 0.04
# chain<- rwMH_sampler(toy1_log_post, init = rep(0, d), prop_var = 0.04, verbose = TRUE)
# plot_traces(chain, "1", beta_true)

## step-size = 0.045
#chain<- MALA_sampler(toy1_log_post, grad_log_fn = toy1_grad_log_post, init = rep(0, d), step_size = 0.085, verbose = TRUE)
#plot_traces(chain, "1", beta_true)


num_samples <- 1e4
burnin <- 500
num_chains <- 100
step_size <- 0.035

toy1_mean_estimates_post <- matrix(NA, nrow = num_chains, ncol = d)
toy1_mean_estimates_opt <- matrix(NA, nrow = num_chains, ncol = d)
toy1_mean_estimates_opt_SIS <- matrix(NA, nrow = num_chains, ncol = d)

chain<- rwMH_sampler(toy1_log_post, init = rep(0, d), prop_var =  0.04, verbose = TRUE)


for(i in 1:num_chains){
  if(i == 1){
    chain<- rwMH_sampler(function(beta){toy1_log_opt_SIS(beta, 1)}, init = rep(0, d), prop_var = step_size, verbose = TRUE)
    #plot_traces(chain, "SNIS", beta_true)
  }
  else chain<- rwMH_sampler(function(beta){toy1_log_opt_SIS(beta, i)}, init = rep(0, d), prop_var = step_size)
  weights <-exp( apply(chain[burnin:num_samples, , drop = FALSE], 1, toy1_log_post) - apply(chain[burnin:num_samples, , drop = FALSE], 1, toy1_log_opt_SIS))
  
  toy1_mean_estimates_opt_SIS[i, ] <- colSums(weights * chain[burnin:num_samples, , drop = FALSE]) / sum(weights)
}


# After computing weights:
ESS <- (sum(weights)^2) / sum(weights^2)
cat("Effective Sample Size:", ESS, "\n")

# Check for d = 12
toy1_mean_estimates_post2 <- matrix(NA, nrow = num_chains, ncol = d)
for(i in 1:num_chains){
  if(i == 1){
    #chain<- MALA_sampler(toy1_log_post, grad_log_fn = toy1_grad_log_post, init = rep(0, d), step_size = step_size, verbose = TRUE)
    chain<- rwMH_sampler(toy1_log_post, num_samples = 2e4, init = rep(0, d), prop_var = step_size, verbose = TRUE)
    #plot_traces(chain, "true posterior", beta_true)
  }
  else chain <-rwMH_sampler(toy1_log_post, num_samples = 2e4, init = rep(0, d), prop_var = step_size)
  toy1_mean_estimates_post2[i, ] <- colMeans(chain[burnin:num_samples, ])
}


colMeans(scale(toy1_mean_estimates_post2, center = TRUE, scale = FALSE)^2)

colMeans(scale(toy1_mean_estimates_opt, center = TRUE, scale = FALSE)^2)

sum(colMeans(scale(toy1_mean_estimates_opt_SIS, center = TRUE, scale = FALSE)^2))/sum(colMeans(scale(toy1_mean_estimates_post2, center = TRUE, scale = FALSE)^2))



for(i in 1:num_chains){
  if(i == 1){
    #chain<- MALA_sampler(toy1_log_post, grad_log_fn = toy1_grad_log_post, init = rep(0, d), step_size = step_size, verbose = TRUE)
    chain<- rwMH_sampler(toy1_log_post, init = rep(0, d), prop_var = step_size, verbose = TRUE)
    #plot_traces(chain, "true posterior", beta_true)
  }
  else chain <-rwMH_sampler(toy1_log_post, init = rep(0, d), prop_var = step_size)
  toy1_mean_estimates_post[i, ] <- colMeans(chain[burnin:num_samples, ])
}

chain<- rwMH_sampler(function(beta){toy1_log_opt(beta, 1)}, init = rep(0, d), prop_var = step_size, verbose = TRUE)

for(i in 1:num_chains){
  if(i == 1){
    chain<- rwMH_sampler(function(beta){toy1_log_opt(beta, 1)}, init = rep(0, d), prop_var = step_size, verbose = TRUE)
    #plot_traces(chain, "SNIS", beta_true)
  }
  else chain<- rwMH_sampler(function(beta){toy1_log_opt(beta, i)}, init = rep(0, d), prop_var = step_size)
  weights <-exp( apply(chain[burnin:num_samples, , drop = FALSE], 1, toy1_log_post) - apply(chain[burnin:num_samples, , drop = FALSE], 1, toy1_log_opt))
  
  toy1_mean_estimates_opt[i, ] <- colSums(weights * chain[burnin:num_samples, , drop = FALSE]) / sum(weights)
}


sum(colMeans(scale(toy1_mean_estimates_opt, center = TRUE, scale = FALSE)^2))/sum(colMeans(scale(toy1_mean_estimates_post2, center = TRUE, scale = FALSE)^2))
