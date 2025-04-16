# data <-read.csv("/Users/arqam/proj/bpr_sampler/betting.csv")
# head(data)
# 
# y<- data[,1]
# X <- as.matrix(data[,-1])
# head(X)

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


#### POISSON REGRESSION

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
  return( 0.5*log(sum((beta-mu_hat)^2)) + log_prior_poisson(beta, mu_prior, sigma_prior) + log_likelihood_poisson(beta, X, y)) 
}

grad_hess_f <- function(u, beta_current, lambda, X, y, mu_prior, sigma_prior, mu_hat, eps = 1e-8, hess = TRUE) {
  eta <- X %*% u
  lambda_pois <- as.vector(exp(eta))
  # Prevent numerical instability
  lambda_pois <- pmin(lambda_pois, 1e100)
  lambda_pois <- pmax(lambda_pois, 1e-100)
  
  n <- length(u)
  diff <- u - mu_hat
  norm_sq <- sum(diff^2) + eps
  inv_norm_sq <- 1/norm_sq
  inv_norm_sq_squared <- inv_norm_sq^2
  inv_lambda <- 1/lambda
  inv_sigma_prior_squared <- 1/sigma_prior^2
  
  # Gradient components
  grad_negloglik <- -t(X) %*% (y - lambda_pois)
  grad_neglogprior <- (u - mu_prior) * inv_sigma_prior_squared
  grad_neglognorm <- -diff * inv_norm_sq
  grad_prox <- (u - beta_current) * inv_lambda
  total_grad <- grad_negloglik + grad_neglogprior + grad_neglognorm + grad_prox
  
  if(hess) {
    identity_matrix <- diag(1, n)
    diff_outer <- tcrossprod(diff)
    X_weighted <- sweep(X, 1, sqrt(lambda_pois), "*")
    hess_negloglik <- crossprod(X_weighted)
    hess_neglogprior <- inv_sigma_prior_squared * identity_matrix
    hess_neglognorm <- -identity_matrix * inv_norm_sq + 2 * diff_outer * inv_norm_sq_squared
    hess_prox <- inv_lambda * identity_matrix
    total_hess <- hess_negloglik + hess_neglogprior + hess_neglognorm + hess_prox
  }
  
  list(
    gradient = total_grad,
    hessian = if(hess) total_hess else NULL,
    components = list(
      grad_negloglik = grad_negloglik,
      grad_neglogprior = grad_neglogprior,
      grad_neglognorm = grad_neglognorm,
      grad_prox = grad_prox
    ),
    diagnostics = list(
      lambda_pois = lambda_pois,
      norm_distance = sqrt(norm_sq)
    )
  )
}

moreau_opt_poisson <- function(beta, X, y, mu_prior, sigma_prior, lambda, mu_hat, NR = TRUE) {
  beta_tilde <- beta
  max_iter <- 100
  tol <- 1e-6
  learning_rate <- 0.01
  n <- length(beta)
  
  for (iter in 1:max_iter) {
    gh <- grad_hess_f(
      u = beta_tilde, 
      beta_current = beta,
      lambda = lambda,
      X = X, y = y,
      mu_prior = mu_prior,
      sigma_prior = sigma_prior,
      mu_hat = mu_hat,
      hess = NR
    )
    
    if(NR) {
      H_reg <- gh$hessian + 1e-4 * diag(n)
      
      delta <- tryCatch(
        qr.solve(H_reg, gh$gradient),
        error = function(e) tryCatch(
          MASS::ginv(H_reg) %*% gh$gradient,
          error = function(e) rep(NA, n)
        )
      )
      
      # Switch to GD if invalid NR step
      if(any(is.na(delta))) {
        delta <- learning_rate * gh$gradient
        print(paste("Switching to gradient descent for", iter, "th step"))
      }
    } else {
      delta <- learning_rate * gh$gradient
    }
    
    beta_tilde_new <- beta_tilde - delta
    
    if(any(is.na(beta_tilde_new))) {
      beta_tilde <- beta
      break
    }
    
    beta_tilde <- beta_tilde_new
    if(max(abs(gh$gradient)) < tol) break
  }
  
  list(
    log_density = log_opt_poisson(beta_tilde, X, y, mu_hat, mu_prior, sigma_prior) - 
      sum((beta_tilde - beta)^2)/(2*lambda),
    beta_tilde = beta_tilde,
    grad_log_density = (beta_tilde - beta)/lambda
  )
}




  

poisson_PMALA_sampler <- function(X, y, mu_prior = 0, sigma_prior = 10, lambda = 0.01,
                                  init = 0, num_samples = 1e4, step_size = 0.1, 
                                  verbose = TRUE, mu_hat = 0, NR = TRUE) 
{ 
  
  num_accepts <- 0
  d <- length(init)
  chain <- matrix(NA, num_samples, d)
  log_opt_density <- numeric(num_samples)
  chain[1, ] <- init
  
  MOP_current <- moreau_opt_poisson(init, X, y, mu_prior, sigma_prior, lambda, mu_hat, NR= NR)
  log_opt_density[1] <- MOP_current$log_density
  
  for (i in 2:num_samples) {
    current <- chain[i-1, ]
    
    prop <- current + step_size * MOP_current$grad_log_density/2 + sqrt(step_size) * rnorm(d)
    
    MOP_prop <- moreau_opt_poisson(prop, X, y, mu_prior, sigma_prior, lambda, mu_hat, NR = NR)
    
    log_target_ratio <- MOP_prop$log_density - MOP_current$log_density
    
    log_q_forward <- sum(dnorm(prop, current + step_size*MOP_current$grad_log_density/2, 
                               sqrt(step_size), log = TRUE))
    log_q_backward <- sum(dnorm(current, prop + step_size*MOP_prop$grad_log_density/2, 
                                sqrt(step_size), log = TRUE))
    
    log_alpha <- log_target_ratio + log_q_backward - log_q_forward
    log_alpha <- ifelse(is.na(log_alpha), -Inf, log_alpha)
    
    if (log(runif(1)) < log_alpha) {
      chain[i, ] <- prop
      log_opt_density[i] <- MOP_prop$log_density
      MOP_current <- MOP_prop
      num_accepts <- num_accepts + 1
    } else {
      chain[i, ] <- current
      log_opt_density[i] <- MOP_current$log_density
    }
    
    if(i%%1000 == 0){
      cat(i)
    }
  }
  
  
  if (verbose) message("Acceptance rate: ", round(num_accepts/(num_samples-1), 3))
  return(list(chain = chain, log_opt_density = log_opt_density))
}

library(profvis)
profvis({
  chain.1 <- poisson_PMALA_sampler(X, y, init= rep(0, dim(X)[2]) , num_samples = 1e2, step_size = 2e-3, lambda = 0.001, mu_hat = mle_coef_manual)
})


chain.1 <- poisson_PMALA_sampler(X, y, init= rep(0, dim(X)[2]) , num_samples = 1e4, step_size = 2.7e-3, lambda = 0.001, mu_hat = mle_coef_manual, NR = TRUE)


# library(tictoc)
# tic()
# # chain.1 <- poisson_PMALA_sampler(X, y, init= mle_coef_manual, num_samples = 1e5, step_size = 2e-3, lambda = 0.001)
# 
# toc()


dim(chain.1$chain)
# install.packages("fastmap")
# devtools::install_github("dvats/SimTools", ref = "Siddharth-Pathak")
library(SimTools)
smc <- as.Smcmc(chain.1$chain)

plot(smc)
weights <- exp(apply(chain.1$chain, 1, log_post_poisson, X, y) - chain.1$log_opt_density)
colSums(weights * chain.1$chain)/ sum(weights)
ESS <- (sum(weights)^2) / sum(weights^2)
ESS

colMeans(chain.1$chain)


mle_coef_manual
