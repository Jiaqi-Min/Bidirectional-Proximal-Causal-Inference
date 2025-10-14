library(readr)
library(foreach)
library(doParallel)     
library(MASS)
library(numDeriv)

set.seed(123)

data <- read_csv("levitt20.csv")

selected_vars <- c("year", "fips", "murder_rate", "ear_murd", "afdc15", "prison", "population", "unemp", "income", "pover")
data_selected <- data[selected_vars]

new_colnames <- c("year", "fips", "Y", "X", "Z", "W", paste0("V", 1:4))
data_clean <- na.omit(data_selected)
colnames(data_clean) <- new_colnames
data <- data_clean

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

data$Z <- standardize(data$Z)
data$V1 <- standardize(data$V1)

GMMF <- function(para, mrf, data) {
  g0 <- mrf(para, data)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g^2)
  return(gmmf)
}

G <- function(para, bfun, data) {
  G1 <- function(para) {
    G1 <- apply(bfun(para, data), 2, mean)
    return(G1)
  }
  G <- jacobian(func=G1, x=para)
  return(G)
}

VAREST <- function(para, bfun, data) {
  bG <- solve(G(para, bfun, data))
  bg <- bfun(para, data)
  spsz <- dim(bg)[1]
  Omega <- t(bg) %*% bg / spsz
  Sigma <- bG %*% Omega %*% t(bG)
  return(Sigma / spsz)
}

Bi_TSLS <- function(data){
  V_names <- paste0("V", 1:4)
  V_terms <- paste(V_names, collapse = " + ")
  
  ZV_interactions <- paste("Z:", V_names, collapse = " + ")
  WV_interactions <- paste("W:", V_names, collapse = " + ")
  
  formula_W <- as.formula(paste("W ~ Z +", V_terms, "+", ZV_interactions))
  formula_Z <- as.formula(paste("Z ~ W +", V_terms, "+", WV_interactions))
  
  fit_W <- lm(formula_W, data = data)
  fit_Z <- lm(formula_Z, data = data)
  
  formula_X_Z <- as.formula(paste("X ~ Z + fit_W$fitted.values +", V_terms))
  formula_Y_Z <- as.formula(paste("Y ~ Z + fit_W$fitted.values +", V_terms))
  formula_X_W <- as.formula(paste("X ~ fit_Z$fitted.values + W +", V_terms))
  formula_Y_W <- as.formula(paste("Y ~ fit_Z$fitted.values + W +", V_terms))
  
  lm_X_Z <- lm(formula_X_Z, data = data)
  lm_Y_Z <- lm(formula_Y_Z, data = data)
  lm_X_W <- lm(formula_X_W, data = data)
  lm_Y_W <- lm(formula_Y_W, data = data)
  
  beta_yx <- coef(lm_X_W)["W"] / coef(lm_Y_W)["W"]
  beta_xy <- coef(lm_Y_Z)["Z"] / coef(lm_X_Z)["Z"]
  
  return(c(beta_xy = as.numeric(beta_xy), beta_yx = as.numeric(beta_yx)))
}

IV_estimation <- function(data){
  V_terms <- paste(paste0("V", 1:4), collapse = " + ")
  
  formula_x <- as.formula(paste("X ~ Z + W +", V_terms))
  formula_y <- as.formula(paste("Y ~ Z + W +", V_terms))
  
  model_x <- lm(formula_x, data = data)
  model_y <- lm(formula_y, data = data)
  
  beta_xy <- coef(model_y)["Z"] / coef(model_x)["Z"]
  beta_yx <- coef(model_x)["W"] / coef(model_y)["W"]
  
  return(c(beta_xy = as.numeric(beta_xy), beta_yx = as.numeric(beta_yx)))
}

OLS_estimation <- function(data){
  V_terms <- paste(paste0("V", 1:4), collapse = " + ")
  
  formula_x <- as.formula(paste("X ~ Y + Z + W +", V_terms))
  formula_y <- as.formula(paste("Y ~ X + Z + W +", V_terms))
  
  model_x <- lm(formula_x, data = data)
  model_y <- lm(formula_y, data = data)
  
  beta_xy <- coef(model_y)["X"]
  beta_yx <- coef(model_x)["Y"]
  
  return(c(beta_xy = as.numeric(beta_xy), beta_yx = as.numeric(beta_yx)))
}

calculate_traditional_se <- function(data) {
  V_terms <- paste(paste0("V", 1:4), collapse = " + ")
  
  formula_x_ols <- as.formula(paste("X ~ Y + Z + W +", V_terms))
  formula_y_ols <- as.formula(paste("Y ~ X + Z + W +", V_terms))
  model_x_ols <- lm(formula_x_ols, data = data)
  model_y_ols <- lm(formula_y_ols, data = data)
  ols_se_xy <- summary(model_y_ols)$coefficients["X", "Std. Error"]
  ols_se_yx <- summary(model_x_ols)$coefficients["Y", "Std. Error"]
  
  formula_x_iv <- as.formula(paste("X ~ Z + W +", V_terms))
  formula_y_iv <- as.formula(paste("Y ~ Z + W +", V_terms))
  model_x_iv <- lm(formula_x_iv, data = data)
  model_y_iv <- lm(formula_y_iv, data = data)
  iv_se_xy <- summary(model_y_iv)$coefficients["Z", "Std. Error"] / abs(coef(model_x_iv)["Z"])
  iv_se_yx <- summary(model_x_iv)$coefficients["W", "Std. Error"] / abs(coef(model_y_iv)["W"])
  
  return(list(
    ols_se = c(ols_se_xy, ols_se_yx),
    iv_se = c(iv_se_xy, iv_se_yx)
  ))
}

BiTSLS_moment_function_XY <- function(xi, data) {
  Z <- as.matrix(data$Z); W <- as.matrix(data$W) 
  V1 <- as.matrix(data$V1); V2 <- as.matrix(data$V2)
  V3 <- as.matrix(data$V3); V4 <- as.matrix(data$V4)
  X <- as.matrix(data$X); Y <- as.matrix(data$Y)
  
  delta_0 <- xi[1]; delta_z <- xi[2]; delta_v1 <- xi[3]; delta_v2 <- xi[4]
  delta_v3 <- xi[5]; delta_v4 <- xi[6]; delta_zv1 <- xi[7]; delta_zv2 <- xi[8]
  delta_zv3 <- xi[9]; delta_zv4 <- xi[10]
  
  theta_0 <- xi[11]; theta_z <- xi[12]; theta_w <- xi[13]; theta_v1 <- xi[14]
  theta_v2 <- xi[15]; theta_v3 <- xi[16]; theta_v4 <- xi[17]
  
  mu_0 <- xi[18]; mu_z <- xi[19]; mu_w <- xi[20]; mu_v1 <- xi[21]
  mu_v2 <- xi[22]; mu_v3 <- xi[23]; mu_v4 <- xi[24]
  
  W_tilde <- delta_0 + delta_z*Z + delta_v1*V1 + delta_v2*V2 + delta_v3*V3 + delta_v4*V4 +
    delta_zv1*Z*V1 + delta_zv2*Z*V2 + delta_zv3*Z*V3 + delta_zv4*Z*V4
  
  instruments_W <- cbind(1, Z, V1, V2, V3, V4, Z*V1, Z*V2, Z*V3, Z*V4)
  instruments_XY <- cbind(1, Z, W_tilde, V1, V2, V3, V4)
  
  residual_W <- W - W_tilde
  moment1 <- instruments_W * as.vector(residual_W)
  
  residual_X <- X - (theta_0 + theta_z*Z + theta_w*W_tilde + theta_v1*V1 + theta_v2*V2 + theta_v3*V3 + theta_v4*V4)
  moment2 <- instruments_XY * as.vector(residual_X)
  
  residual_Y <- Y - (mu_0 + mu_z*Z + mu_w*W_tilde + mu_v1*V1 + mu_v2*V2 + mu_v3*V3 + mu_v4*V4)
  moment3 <- instruments_XY * as.vector(residual_Y)
  
  G <- cbind(moment1, moment2, moment3)
  return(G)
}

BiTSLS_moment_function_YX <- function(xi, data) {
  Z <- as.matrix(data$Z); W <- as.matrix(data$W) 
  V1 <- as.matrix(data$V1); V2 <- as.matrix(data$V2)
  V3 <- as.matrix(data$V3); V4 <- as.matrix(data$V4)
  X <- as.matrix(data$X); Y <- as.matrix(data$Y)
  
  eta_0 <- xi[1]; eta_w <- xi[2]; eta_v1 <- xi[3]; eta_v2 <- xi[4]
  eta_v3 <- xi[5]; eta_v4 <- xi[6]; eta_wv1 <- xi[7]; eta_wv2 <- xi[8]
  eta_wv3 <- xi[9]; eta_wv4 <- xi[10]
  
  theta_0 <- xi[11]; theta_z <- xi[12]; theta_w <- xi[13]; theta_v1 <- xi[14]
  theta_v2 <- xi[15]; theta_v3 <- xi[16]; theta_v4 <- xi[17]
  
  mu_0 <- xi[18]; mu_z <- xi[19]; mu_w <- xi[20]; mu_v1 <- xi[21]
  mu_v2 <- xi[22]; mu_v3 <- xi[23]; mu_v4 <- xi[24]
  
  Z_tilde <- eta_0 + eta_w*W + eta_v1*V1 + eta_v2*V2 + eta_v3*V3 + eta_v4*V4 +
    eta_wv1*W*V1 + eta_wv2*W*V2 + eta_wv3*W*V3 + eta_wv4*W*V4
  
  instruments_Z <- cbind(1, W, V1, V2, V3, V4, W*V1, W*V2, W*V3, W*V4)
  instruments_XY <- cbind(1, Z_tilde, W, V1, V2, V3, V4)
  
  residual_Z <- Z - Z_tilde
  moment1 <- instruments_Z * as.vector(residual_Z)
  
  residual_X <- X - (theta_0 + theta_z*Z_tilde + theta_w*W + theta_v1*V1 + theta_v2*V2 + theta_v3*V3 + theta_v4*V4)
  moment2 <- instruments_XY * as.vector(residual_X)
  
  residual_Y <- Y - (mu_0 + mu_z*Z_tilde + mu_w*W + mu_v1*V1 + mu_v2*V2 + mu_v3*V3 + mu_v4*V4)
  moment3 <- instruments_XY * as.vector(residual_Y)
  
  G <- cbind(moment1, moment2, moment3)
  return(G)
}

construct_initial_values_XY <- function(data) {
  V_names <- paste0("V", 1:4)
  V_terms <- paste(V_names, collapse = " + ")
  ZV_interactions <- paste("Z:", V_names, collapse = " + ")
  
  formula_W <- as.formula(paste("W ~ Z +", V_terms, "+", ZV_interactions))
  fit_W <- lm(formula_W, data = data)
  formula_X_Z <- as.formula(paste("X ~ Z + fit_W$fitted.values +", V_terms))
  formula_Y_Z <- as.formula(paste("Y ~ Z + fit_W$fitted.values +", V_terms))
  lm_X_Z <- lm(formula_X_Z, data = data)
  lm_Y_Z <- lm(formula_Y_Z, data = data)
  inioptim <- c(coef(fit_W), coef(lm_X_Z), coef(lm_Y_Z))
  names(inioptim) <- NULL
  return(inioptim)
}

construct_initial_values_YX <- function(data) {
  V_names <- paste0("V", 1:4)
  V_terms <- paste(V_names, collapse = " + ")
  WV_interactions <- paste("W:", V_names, collapse = " + ")
  
  formula_Z <- as.formula(paste("Z ~ W +", V_terms, "+", WV_interactions))
  fit_Z <- lm(formula_Z, data = data)
  formula_X_W <- as.formula(paste("X ~ fit_Z$fitted.values + W +", V_terms))
  formula_Y_W <- as.formula(paste("Y ~ fit_Z$fitted.values + W +", V_terms))
  lm_X_W <- lm(formula_X_W, data = data)
  lm_Y_W <- lm(formula_Y_W, data = data)
  inioptim <- c(coef(fit_Z), coef(lm_X_W), coef(lm_Y_W))
  names(inioptim) <- NULL
  return(inioptim)
}

theoretical_method <- function(data) {
  bitsls_results <- Bi_TSLS(data)
  beta_xy_bitsls <- bitsls_results["beta_xy"]
  beta_yx_bitsls <- bitsls_results["beta_yx"]
  
  inioptim_xy <- construct_initial_values_XY(data)
  var_matrix_xy <- VAREST(inioptim_xy, BiTSLS_moment_function_XY, data)
  mu_z <- inioptim_xy[19]
  theta_z <- inioptim_xy[12]
  grad_xy <- rep(0, length(inioptim_xy))
  grad_xy[19] <- 1/theta_z
  grad_xy[12] <- -mu_z/(theta_z^2)
  var_beta_xy <- as.numeric(t(grad_xy) %*% var_matrix_xy %*% grad_xy)
  se_xy <- sqrt(var_beta_xy)
  
  inioptim_yx <- construct_initial_values_YX(data)
  var_matrix_yx <- VAREST(inioptim_yx, BiTSLS_moment_function_YX, data)
  theta_w <- inioptim_yx[13]
  mu_w <- inioptim_yx[20]
  grad_yx <- rep(0, length(inioptim_yx))
  grad_yx[13] <- 1/mu_w
  grad_yx[20] <- -theta_w/(mu_w^2)
  var_beta_yx <- as.numeric(t(grad_yx) %*% var_matrix_yx %*% grad_yx)
  se_yx <- sqrt(var_beta_yx)
  
  return(list(
    beta_xy = beta_xy_bitsls,
    beta_yx = beta_yx_bitsls,
    se_xy = se_xy,
    se_yx = se_yx
  ))
}

empirical_method <- function(data, n_bootstrap = 500) {
  original_results <- Bi_TSLS(data)
  beta_xy_original <- original_results["beta_xy"]
  beta_yx_original <- original_results["beta_yx"]
  
  n_obs <- nrow(data)
  bootstrap_estimates_xy <- numeric(n_bootstrap)
  bootstrap_estimates_yx <- numeric(n_bootstrap)
  
  for(b in 1:n_bootstrap) {
    set.seed(123 + b)
    boot_indices <- sample(1:n_obs, n_obs, replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    boot_results <- Bi_TSLS(boot_data)
    bootstrap_estimates_xy[b] <- boot_results["beta_xy"]
    bootstrap_estimates_yx[b] <- boot_results["beta_yx"]
  }
  
  se_xy <- sd(bootstrap_estimates_xy)
  se_yx <- sd(bootstrap_estimates_yx)
  
  return(list(
    beta_xy = beta_xy_original,
    beta_yx = beta_yx_original,
    se_xy = se_xy,
    se_yx = se_yx
  ))
}

run_bootstrap_parallel <- function(data, n_bootstrap = 500) {
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  clusterExport(cl, c("Bi_TSLS", "OLS_estimation", "IV_estimation", 
                      "theoretical_method", "empirical_method",
                      "GMMF", "G", "VAREST",
                      "BiTSLS_moment_function_XY", "BiTSLS_moment_function_YX",
                      "construct_initial_values_XY", "construct_initial_values_YX"))
  clusterEvalQ(cl, library(MASS))
  clusterEvalQ(cl, library(numDeriv))
  
  results <- foreach(i = 1:n_bootstrap, .combine = 'rbind') %dopar% {
    set.seed(123 + i)
    bootstrap_indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    bootstrap_data <- data[bootstrap_indices, ]
    
    ols_est <- OLS_estimation(bootstrap_data)
    iv_est <- IV_estimation(bootstrap_data)
    bitsls_est <- Bi_TSLS(bootstrap_data)
    theo_result <- theoretical_method(bootstrap_data)
    
    c(ols_est, iv_est, bitsls_est, theo_result$se_xy, theo_result$se_yx)
  }
  
  stopCluster(cl)
  return(results)
}

calculate_final_results <- function(bootstrap_results, original_data) {
  ols_orig <- OLS_estimation(original_data)
  iv_orig <- IV_estimation(original_data)
  bitsls_orig <- Bi_TSLS(original_data)
  bitsls_theo_orig <- theoretical_method(original_data)
  
  traditional_se <- calculate_traditional_se(original_data)
  
  ols_se_xy_emp <- sd(bootstrap_results[, 1])
  ols_se_yx_emp <- sd(bootstrap_results[, 2])
  
  iv_se_xy_emp <- sd(bootstrap_results[, 3])
  iv_se_yx_emp <- sd(bootstrap_results[, 4])
  
  bitsls_se_xy_emp <- sd(bootstrap_results[, 5])
  bitsls_se_yx_emp <- sd(bootstrap_results[, 6])
  
  final_results <- data.frame(
    Method = c("OLS", "IV", "Bi-TSLS (Empirical)", "Bi-TSLS (Theoretical)"),
    X_to_Y_Estimate = c(ols_orig["beta_xy"], iv_orig["beta_xy"], 
                        bitsls_orig["beta_xy"], bitsls_orig["beta_xy"]),
    X_to_Y_SE = c(traditional_se$ols_se[1], traditional_se$iv_se[1], 
                  bitsls_se_xy_emp, bitsls_theo_orig$se_xy),
    Y_to_X_Estimate = c(ols_orig["beta_yx"], iv_orig["beta_yx"], 
                        bitsls_orig["beta_yx"], bitsls_orig["beta_yx"]),
    Y_to_X_SE = c(traditional_se$ols_se[2], traditional_se$iv_se[2], 
                  bitsls_se_yx_emp, bitsls_theo_orig$se_yx),
    stringsAsFactors = FALSE
  )
  
  final_results$X_to_Y_CI_Lower <- final_results$X_to_Y_Estimate - 1.96 * final_results$X_to_Y_SE
  final_results$X_to_Y_CI_Upper <- final_results$X_to_Y_Estimate + 1.96 * final_results$X_to_Y_SE
  final_results$Y_to_X_CI_Lower <- final_results$Y_to_X_Estimate - 1.96 * final_results$Y_to_X_SE
  final_results$Y_to_X_CI_Upper <- final_results$Y_to_X_Estimate + 1.96 * final_results$Y_to_X_SE
  
  numeric_cols <- c("X_to_Y_Estimate", "X_to_Y_SE", "Y_to_X_Estimate", "Y_to_X_SE",
                    "X_to_Y_CI_Lower", "X_to_Y_CI_Upper", "Y_to_X_CI_Lower", "Y_to_X_CI_Upper")
  final_results[numeric_cols] <- lapply(final_results[numeric_cols], function(x) round(x, 4))
  
  return(final_results)
}

bootstrap_results <- run_bootstrap_parallel(data, n_bootstrap = 500)
final_results <- calculate_final_results(bootstrap_results, data)

print(final_results)
write.csv(final_results, "Real_Data_Results.csv", row.names = FALSE)