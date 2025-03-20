
rwMH_sampler <- function(log_density_fn, init = 0, num_samples = 1e4, prop_var = 1, verbose = FALSE){
  num_accepts <- 0
  
  chain <- numeric(length = num_samples)
  chain[1] <- init
  
  for(i in 2:num_samples){
    prop <- chain[i-1] + sqrt(prop_var) * rnorm(n = 1, mean = 0, sd = 1)
    log_MH_ratio <- log_density_fn(prop) - log_density_fn(chain[i-1])
    
    
    if(log(runif(1)) < log_MH_ratio){
      chain[i] <- prop
      num_accepts <- num_accepts + 1
    }
    else{
      chain[i] <- chain[i-1]
    }
  }
  if(verbose) print(paste("Acceptance rate is", num_accepts/num_samples))
  return(chain)
}


MALA_sampler <- function(log_density_fn, grad_log_fn, init = 0, num_samples = 1e4, step_size = 0.1, verbose = FALSE) {
  num_accepts <- 0
  
  chain <- numeric(length = num_samples)
  chain[1] <- init
  
  for(i in 2:num_samples) {
    current <- chain[i-1]
    grad <- grad_log_fn(current)
    
    prop <- current + 0.5 * step_size * grad + sqrt(step_size) * rnorm(1)
    
    grad_prop <- grad_log_fn(prop)
    log_q_forward <- -((prop - current - 0.5 * step_size * grad)^2) / (2 * step_size)
    log_q_backward <- -((current - prop - 0.5 * step_size * grad_prop)^2) / (2 * step_size)
    log_MH_ratio <- log_density_fn(prop) - log_density_fn(current) + log_q_backward - log_q_forward
    
    if(log(runif(1)) < log_MH_ratio) {
      chain[i] <- prop
      num_accepts <- num_accepts + 1
    } else {
      chain[i] <- current
    }
  }
  
  if(verbose) print(paste("MALA Acceptance rate is", num_accepts / num_samples))
  return(chain)
}


HMC_sampler <- function(log_density_fn, grad_log_fn, init = 0, num_samples = 1e4, step_size = 0.1, L = 10, verbose = FALSE) {
  num_accepts <- 0
  
  chain <- numeric(length = num_samples)
  chain[1] <- init
  
  for(i in 2:num_samples) {
    current <- chain[i-1]
    momentum <- rnorm(1)  # Sample initial momentum
    
    # Leapfrog integration
    prop <- current
    prop_momentum <- momentum + 0.5 * step_size * grad_log_fn(prop)  # Half step in momentum
    
    for(l in 1:L) {
      prop <- prop + step_size * prop_momentum  # Full step in position
      if(l != L) prop_momentum <- prop_momentum + step_size * grad_log_fn(prop)  # Full step in momentum
    }
    
    prop_momentum <- prop_momentum + 0.5 * step_size * grad_log_fn(prop)  # Final half step in momentum
    prop_momentum <- -prop_momentum  # Negate momentum (reverse dynamics)
    
    log_MH_ratio <- log_density_fn(prop) - log_density_fn(current) + 
      (momentum^2 - prop_momentum^2) / 2  # Kinetic energy term
    
    if(log(runif(1)) < log_MH_ratio) {
      chain[i] <- prop
      num_accepts <- num_accepts + 1
    } else {
      chain[i] <- current
    }
  }
  
  if(verbose) print(paste("HMC Acceptance rate is", num_accepts / num_samples))
  return(chain)
}





true_mean <- 3
num_samples <- 10
sd <- 1

samples <- rnorm(num_samples, true_mean, sd)

analytical_posterior_mean <- sum(samples)/(num_samples +1)
analytical_posterior_var <- 1/(num_samples + 1)

prior_mean <- 0
prior_sd <- 1

toy1_log_prior <- function(mu){
  return(dnorm(mu, mean = prior_mean, sd = prior_sd, log = TRUE))
}

toy1_log_likelihood <- function(mu){
  return(sum( dnorm( samples, mean = mu, sd = 1, log = TRUE) ) )
}

toy1_moments <- function(mu){
  return(c(mu, mu^2))
}

toy1_log_opt_density <- function(mu){
  return(toy1_log_prior(mu) +toy1_log_likelihood(mu) + log(sqrt(sum(toy1_moments(mu)^2))))
}

toy1_log_opt_SNIS<- function(mu){
  true_moments <- c(analytical_posterior_mean, analytical_posterior_var+analytical_posterior_mean^2)
  return(toy1_log_prior(mu) +toy1_log_likelihood(mu) + log(sqrt(sum((toy1_moments(mu)-true_moments)^2))))
}

toy1_log_posterior <- function(mu){
  return( toy1_log_prior(mu) + toy1_log_likelihood(mu) ) 
}

true_mean <- 3
num_samples <- 10
sd <- 1

samples <- rnorm(num_samples, true_mean, sd)

analytical_posterior_mean <- sum(samples)/(num_samples +1)
analytical_posterior_var <- 1/(num_samples + 1)

prior_mean <- 0
prior_sd <- 1

toy1_log_prior <- function(mu){
  return(dnorm(mu, mean = prior_mean, sd = prior_sd, log = TRUE))
}

toy1_log_likelihood <- function(mu){
  return(sum( dnorm( samples, mean = mu, sd = 1, log = TRUE) ) )
}

toy1_moments <- function(mu){
  return(c(mu, mu^2))
}

toy1_log_opt_density <- function(mu){
  return(toy1_log_prior(mu) +toy1_log_likelihood(mu) + log(sqrt(sum(toy1_moments(mu)^2))))
}

toy1_log_opt_SNIS<- function(mu){
  true_moments <- c(analytical_posterior_mean, analytical_posterior_var+analytical_posterior_mean^2)
  return(toy1_log_prior(mu) +toy1_log_likelihood(mu) + log(sqrt(sum((toy1_moments(mu)-true_moments)^2))))
}

toy1_log_posterior <- function(mu){
  return( toy1_log_prior(mu) + toy1_log_likelihood(mu) ) 
}

mu_values <- seq(0, 6, length.out = 1000)
log_opt_density <- sapply(mu_values, toy1_log_opt_density)
log_opt_SNIS <- sapply(mu_values, toy1_log_opt_SNIS)
log_posterior_density <- sapply(mu_values, toy1_log_posterior)

opt_density <- exp(log_opt_density - max(log_opt_density))  
opt_density_SNIS <- exp(log_opt_SNIS - max(log_opt_SNIS))
posterior_density <- exp(log_posterior_density - max(log_posterior_density))


plot(mu_values, opt_density, type = "l", col = "blue", lwd = 2, 
     ylab = "Unnormalized Density", xlab = expression(mu),
     main = "Unnormalized Densities", ylim = c(0, max(opt_density, posterior_density)))

lines(mu_values, posterior_density, col = "red", lwd = 2)
lines(mu_values, opt_density_SNIS, col = "green", lwd = 2)


legend("topright", legend = c("SNIS optimal", "Pseud Optimal", "Posterior"),
       col = c("green", "blue", "red"), lwd = 2)

burnin <- 2000
num_chains <- 100

toy1_mean_estimates_posterior <- numeric(length = num_chains)
toy1_var_estimates_posterior <- numeric(length = num_chains)

toy1_mean_estimates_opt <- numeric(length = num_chains)
toy1_var_estimates_opt <- numeric(length = num_chains)

toy1_mean_estimates_snis_opt <- numeric(length = num_chains)
toy1_var_estimates_snis_opt <- numeric(length = num_chains)


for(i in 1:num_chains){
  if(i == 1){
    chain <- rwMH(toy1_log_posterior, prop_var = 0.49, num_samples = 1e5, verbose = TRUE)
    plot(chain)
  }
  else chain <- rwMH(toy1_log_posterior, prop_var = 0.49, num_samples = 1e5)
  toy1_mean_estimates_posterior[i] <- mean(chain[burnin:num_samples])
  toy1_var_estimates_posterior[i] <- mean(chain[burnin:num_samples]^2) - toy1_mean_estimates_posterior[i]^2
}

for(i in 1:num_chains){
  if(i == 1){
    chain <- rwMH(toy1_log_opt_density, prop_var = 0.49, num_samples = 1e5, verbose = TRUE)
    plot(chain)
  }
  else chain <- rwMH(toy1_log_opt_density, prop_var = 0.49, num_samples = 1e5)
  weights <-exp( sapply( chain[burnin:num_samples] , toy1_log_posterior) - sapply( chain[burnin:num_samples], toy1_log_opt_density))
  
  toy1_mean_estimates_opt[i] <- sum(weights*chain[burnin:num_samples])/sum(weights)
  toy1_var_estimates_opt[i] <- sum(weights*chain[burnin:num_samples]^2)/sum(weights)- toy1_mean_estimates_opt[i]^2
}

for(i in 1:num_chains){
  if(i == 1){
    chain <- rwMH(toy1_log_opt_SNIS, prop_var = 0.64, num_samples = 1e5, verbose = TRUE)
    plot(chain)
  }
  else chain <- rwMH(toy1_log_opt_SNIS, prop_var = 0.49, num_samples = 1e5)
  weights <-exp( sapply( chain[burnin:num_samples] , toy1_log_posterior) - sapply( chain[burnin:num_samples], toy1_log_opt_SNIS))
  
  toy1_mean_estimates_snis_opt[i] <- sum(weights*chain[burnin:num_samples])/sum(weights)
  toy1_var_estimates_snis_opt[i] <- sum(weights*chain[burnin:num_samples]^2)/sum(weights)- toy1_mean_estimates_snis_opt[i]^2
}

print(paste("Posterior mean estimator bias" , mean(toy1_mean_estimates_posterior) - analytical_posterior_mean)) 
print(paste("Pseudo-optimal mean estimator bias" , mean(toy1_mean_estimates_opt) - analytical_posterior_mean)) 
print(paste("SNIS optimal mean estimator bias" , mean(toy1_mean_estimates_snis_opt) - analytical_posterior_mean)) 

print(paste("Posterior mean estimator variance" , var(toy1_mean_estimates_posterior)))
print(paste("Pseudo-optimal mean estimator variance" , var(toy1_mean_estimates_opt)))
print(paste("SNIS optimal mean estimator variance" , var(toy1_mean_estimates_snis_opt)))
print(paste("Posterior variance estimator bias" , mean(toy1_var_estimates_posterior) - analytical_posterior_var)) 
print(paste("Pseuod-optimal variance estimator bias" , mean(toy1_var_estimates_opt) - analytical_posterior_var)) 
print(paste("SNIS optimal variance estimator bias" , mean(toy1_var_estimates_snis_opt) - analytical_posterior_var)) 

print(paste("Posterior estimator variance" , var(toy1_var_estimates_posterior)))
print(paste("Pseudo-optimal estimator variance" , var(toy1_var_estimates_opt)))
print(paste("SNIS optimal estimator variance" , var(toy1_var_estimates_snis_opt)))
