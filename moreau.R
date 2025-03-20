moreau_norm <- function(vec, lambda){
  eps <- 1e-2
  vec_norm <- sqrt(sum(vec**2))
  if( vec_norm <= lambda){
    return(vec_norm**2/ (2 * lambda) + eps) 
  }
  else{
    return(vec_norm - lambda/2 + eps)
  }
}

moreau_norm_log_gradient <- function(vec, lambda) {
  eps <- 1e-2
  vec_norm <- sqrt(sum(vec**2))  # L2 norm
  
  if (vec_norm <= lambda) {
    mureau_value <- (vec_norm**2) / (2 * lambda) + eps
    gradient_mureau <- vec / lambda
  } else {
    mureau_value <- vec_norm - (lambda / 2) + eps
    gradient_mureau <- vec / vec_norm
  }
  
  # Gradient of log(mureau_norm)
  gradient_log <- gradient_mureau / mureau_value
  
  return(gradient_log)
}




vec <- c(3, 4)

moreau_norm_log_gradient(vec, 5)

