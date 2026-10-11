library(MASS)
# Run from the folder that contains this package (see README). The line
# below sources a helper that computes the population-level c* / value*.
source("compute_c_star_value_star.R")
generate_simulation_data <- function(n = 50, 
                                     p = 6, 
                                     sigma_noise1 =3.144761e-1, # control noise-rsquare 1--30, 3--4, 6---1, 9---0.5
                                     sigma_noise2 =3.097618e-1, 
                                     muy1 =0.73, #drg
                                     meandiff = 0.05, #muy2 pbo
                                     rho = 0.3,
                                     cov_type = c("equal", "ar1", "real"),
                                     beta1 = NULL,  
                                     beta2 = NULL,
                                     real_data_path = "cortical_thickness_data.csv" ) {
  cov_type <- match.arg(cov_type)
  
  # === Generate covariance matrix psix based on type ===
  if (cov_type == "equal") {
    seq_base_sd <- matrix(0.1, nrow = p, ncol = 1)
    corMat <- matrix(rho, nrow = p, ncol = p)
    diag(corMat) <- 1
    psix <- (seq_base_sd %*% t(seq_base_sd)) * corMat
    
  } else if (cov_type == "ar1") {
    seq_base_sd <- matrix(0.1, nrow = p, ncol = 1)
    corMat <- outer(1:p, 1:p, function(i, j) rho ^ abs(i - j))
    psix <- seq_base_sd %*% t(seq_base_sd) * corMat
    
  } else if (cov_type == "real") {
    dat0 <- read.csv(real_data_path, header = TRUE)
    for (j in 12:ncol(dat0)) {
      dat0[is.na(dat0[, j]), j] <- mean(dat0[, j], na.rm = TRUE)
    }
    dat1 <- dat0[, 12:(12 + p - 1)] # max p=214
    psix <- cov(dat1)
  }
  psixy1 <- psix %*% beta1
  psixy2 <- psix %*% beta2
  y1y1<-t(beta1)%*%psix%*%beta1+sigma_noise1^2 #sigma1
  y2y2<-t(beta2)%*%psix%*%beta2+sigma_noise2^2
  
  muy2 <- muy1-meandiff #y1 is drg y2 is pbo
  
  # Construct covariance matrices for DRG and PBO
  cov1 <- matrix(0, p+1, p+1)
  cov2 <- matrix(0, p+1, p+1)
  cov1[1:p, 1:p] <- psix
  cov2[1:p, 1:p] <- psix
  cov1[1:p, p+1] <- psixy1
  cov1[p+1, 1:p] <- psixy1
  cov2[1:p, p+1] <- psixy2
  cov2[p+1, 1:p] <- psixy2
  cov1[p+1, p+1] <- y1y1
  cov2[p+1, p+1] <- y2y2
  
  #beta1 <- solve(psix)%*%psixy1
  #beta2 <- solve(psix)%*%psixy2
  #bstar <- (beta1-beta2) #true contrast coefficients
  Rsquare1 <- t(beta1)%*%psix%*%beta1/sigma_noise1^2
  Rsquare2 <- t(beta2)%*%psix%*%beta2/sigma_noise2^2
  #browser()
  
  
  # === Compute theoretical c_star and value_star 30June ===###############################
  res_star <- compute_c_star_value_star(psix, psixy1, psixy2, muy1, muy2)
  c_star <- res_star$c_star
  value_star <- res_star$value_star
  # Eigen-decomposition and fix negative eigenvalues ######################################
  e1 <- eigen(cov1)
  e2 <- eigen(cov2)
  
  # 1 is drg, 2 is pbo
  #browser() N(0,cov1)
  simdat1 <- matrix(rnorm((p+1)*n), n, (p+1)) %*% sqrt(diag(e1$values)) %*% t(e1$vectors)
  simdat2 <- matrix(rnorm((p+1)*n), n, (p+1)) %*% sqrt(diag(e2$values)) %*% t(e2$vectors)

  #browser()
  sum(is.na(e1$values))
  x1 <- simdat1[, 1:p]
  y1 <- simdat1[, p+1] + muy1
  x2 <- simdat2[, 1:p]
  y2 <- simdat2[, p+1] + muy2
  
  # Calculate beta coefficients 30June
  #beta1 <- solve(psix) %*% psixy1
  #beta2 <- solve(psix) %*% psixy2
  
  # Expected potential outcomes
  W_drg <- x1 %*% beta1 + muy1
  W_drg_pbo <- x1 %*% beta2 + muy2
  W_pbo <- x2 %*% beta2 + muy2
  W_pbo_drg <- x2 %*% beta1 + muy1
  
  # Optimal treatment (true group) label
  true_group_drg <- ifelse(W_drg > W_drg_pbo, "drg", "pbo")
  true_group_pbo <- ifelse(W_pbo > W_pbo_drg, "pbo", "drg")
  
  # Final dataset
  # data <- data.frame(
  #  subject = 1:(2*n),
  #  group = rep(c("drg", "pbo"), each = n),
  #  true_group = c(true_group_drg, true_group_pbo),
  #  outcome = c(y1, y2),
  #  y_drg = c(y1, rep(NA, n)),
  #  y_pbo = c(rep(NA, n), y2)
  #)
  data <- data.frame(
    subject    = 1:(2*n),
    group      = rep(c("drg", "pbo"), each = n),
    true_group = c(true_group_drg, true_group_pbo),
    outcome    = c(y1, y2),
    y_drg      = c(y1, rep(NA, n)),
    y_pbo      = c(rep(NA, n), y2)
  )
  # old label X1 to XP
  # for (i in 1:p) {
  #  data[[paste0("X", i)]] <- c(x1[, i], x2[, i])
  # }
  # new label X1 to XP-03Dec
  idx_p <- sample.int(p, size = p, replace = FALSE)
  
  for (k in 1:p) {
    j <- idx_p[k]  # which original column goes to Xk
    data[[paste0("X", k)]] <- c(x1[, j], x2[, j])
  }
  
  true_params <- list(
    psix = psix, psixy1 = psixy1, psixy2 = psixy2,
    muy1 = muy1, muy2 = muy2,
    sigmay1 = sigma_noise1, sigmay2 = sigma_noise2,
    beta1 = as.vector(beta1), beta2 = as.vector(beta2)
  )
  
  list(data = data, c_star = c_star, value_star = value_star, true_params = true_params)
}



oracle_pcd_mc <- function(true_params, n_big = 1e6, seed = 20250813) {
  set.seed(seed)
  
  Σx <- true_params$psix
  μ1 <- true_params$muy1
  μ2 <- true_params$muy2
  s1 <- true_params$sigmay1
  s2 <- true_params$sigmay2
  β1 <- true_params$beta1
  β2 <- true_params$beta2
  
  p <- ncol(Σx)
  
  # We know all the parameter from the simulation data
  X <- MASS::mvrnorm(n = n_big, mu = rep(0, p), Sigma = Σx)
  
  #browser()
  mu1x <- μ1 + as.numeric(X %*% β1) #drg 
  mu2x <- μ2 + as.numeric(X %*% β2) #pbo 
  

  Y1 <- mu1x + rnorm(n_big, 0, s1)
  Y2 <- mu2x + rnorm(n_big, 0, s2)

  # if Y2>Y1, pbo group, other wise drg 
  #choose_arm2 <- (mu2x > mu1x)
  #mean(ifelse(choose_arm2, Y2 > Y1, Y1 > Y2))
  
  #modified on 16Dec2025
  # Oracle decision rule (based on means)
  group <- ifelse(mu2x > mu1x, "pbo", "drg")
  # True optimal treatment (based on realized outcomes)
  true_group <- ifelse(Y2 > Y1, "pbo", "drg")
  
  # Oracle PCD
  pcd <- mean(group == true_group)
}


