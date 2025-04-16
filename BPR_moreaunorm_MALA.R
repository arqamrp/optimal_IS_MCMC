# we use (-mu_hat)

library(MASS)
head(epil)
y <- epil$y
X <- model.matrix(~ lbase * trt + lage + V4, data = epil)
head(X)


### MLE est

neg_log_lik <- function(beta, X, y) {
  lambda <- exp(X %*% beta)
  -sum(dpois(y, lambda, log = TRUE))
}

init_beta <- rep(0, ncol(X))
mle_fit <- optim(init_beta, neg_log_lik, X = X, y = y, method = "BFGS", hessian = TRUE)

mle_coef_manual <- mle_fit$par
names(mle_coef_manual) <- colnames(X)
mle_coef_manual

se <- sqrt(diag(solve(mle_fit$hessian)))
se


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


moreau_norm <- function(vec, lambda = 0, eps = 1e-5){
  vec_norm <- sqrt(sum(vec^2))
  if( vec_norm <= lambda & lambda != 0){
    return(vec_norm**2/ (2 * lambda) + eps) 
  }
  else{
    return(vec_norm - lambda/2 + eps)
  }
}


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


log_likelihood_poisson <- function(beta, X, y) 
{
  eta <- X %*% beta
  log_lambda <- eta
  log_liks <- dpois(y, lambda = exp(log_lambda), log = TRUE)
  sum(log_liks)
}

log_prior_poisson <- function(beta, mu_prior, sigma_prior)
{
  log_priors <- dnorm(beta, mean = mu_prior, sd = sigma_prior, log = TRUE)
  sum(log_priors)
}

log_post_poisson <- function(beta, X, y, mu_prior= 0, sigma_prior=10){
  return(log_likelihood_poisson(beta, X, y) + log_prior_poisson(beta, mu_prior, sigma_prior))
}


log_opt_poisson <- function(beta, X, y, mu_hat, mu_prior= 0, sigma_prior=10)
{
  # return( 0.5*log(sum((beta-mu_hat)^2)) + log_prior_poisson(beta, mu_prior, sigma_prior) + log_likelihood_poisson(beta, X, y)) 
  return(log_post_poisson(beta, X, y) + log(moreau_norm(beta-mu_hat)) )
  
}


grad_log_post_poisson <- function(beta, X, y, mu_prior= 0, sigma_prior=10){
  eta <- X %*% beta 
  lambda_pois <- exp(eta)
  grad<- t(X) %*% (y-lambda_pois) - beta/sigma_prior^2
  return(grad)
}

grad_log_opt_poisson <- function(beta, X, y, mu_hat, mu_prior= 0, sigma_prior=10){
  return(grad_log_post_poisson(beta, X, y, mu_prior= 0, sigma_prior=10) + grad_log_moreau_norm(beta- mu_hat))
}

grad_log_opt_poisson_simple <- function(beta){
  grad_log_opt_poisson(beta, X, y, mle_coef_manual)
}

log_opt_poisson_simple <- function(beta){
  log_opt_poisson(beta, X, y, mle_coef_manual)
}

grad_log_post_poisson_simple <- function(beta){
  grad_log_post_poisson(beta, X, y)
}

log_post_poisson_simple <- function(beta){
  log_post_poisson(beta, X, y)
}


chain<- MALA_sampler(log_opt_poisson_simple, grad_log_fn = grad_log_opt_poisson_simple, num_samples = 1e5, init = rep(0, dim(X)[2]), step_size = 5.8e-4, verbose = TRUE)
chain2 <- MALA_sampler(log_post_poisson_simple, grad_log_fn = grad_log_post_poisson_simple, num_samples = 1e5, init = rep(0, dim(X)[2]), step_size = 5.7e-4, verbose = TRUE)

library(SimTools)
smc <- as.Smcmc(chain2)

plot(smc)

weights <- exp(apply(chain, 1, log_post_poisson, X, y) - apply(chain, 1, log_opt_poisson, X, y, mle_coef_manual))
est <- colSums(weights * chain.1$chain)/ sum(weights)

ESS <- (sum(weights)^2) / sum(weights^2)
ESS
length(weights)

colMeans(chain.1$chain)


weights
