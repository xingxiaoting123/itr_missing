compute_c_star_value_star <- function(psix, psixy1, psixy2, muy1, muy2) {
  p <- length(psixy1)
  stopifnot(is.numeric(psixy1), is.numeric(psixy2), is.matrix(psix), ncol(psix) == p)
  
  delta_psixy <- psixy1 - psixy2
  alpha_unnorm <- solve(psix) %*% delta_psixy
  #norm_factor <- sqrt(sum(alpha_unnorm^2))  # Euclidean norm
  norm_factor <- as.numeric(sqrt(t(alpha_unnorm) %*% alpha_unnorm)) # use the matrix norm
  alpha_opt <- as.vector(alpha_unnorm) / norm_factor
  
  numerator <- as.numeric((muy2 - muy1) * sqrt(t(alpha_opt) %*% psix %*% alpha_opt))
  denominator <- as.numeric(abs(t(alpha_opt) %*% delta_psixy))
  c_star <- numerator / denominator
  
  value_star <- muy1 + (muy2 - muy1) * (pnorm(c_star) + dnorm(c_star) / c_star)
  
  return(list(
    c_star = c_star,
    value_star = value_star
  ))
}

