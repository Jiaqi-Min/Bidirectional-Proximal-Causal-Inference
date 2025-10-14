library(dplyr)
library(tidyr)  
library(stringr)
library(numDeriv)
library(foreach)
library(doParallel)
library(parallel)
library(openxlsx)

set.seed(123)

calculate_gmm_objective <- function(model_params, moment_func, dataset) {
  moment_matrix <- moment_func(model_params, dataset)
  mean_moments <- apply(moment_matrix, 2, mean)
  objective_value <- sum(mean_moments^2)
  return(objective_value)
}

calculate_gradient_matrix <- function(model_params, moment_func, dataset) {
  gradient_func <- function(params) {
    mean_moments <- apply(moment_func(params, dataset), 2, mean)
    return(mean_moments)
  }
  gradient_matrix <- jacobian(func = gradient_func, x = model_params)
  return(gradient_matrix)
}

estimate_variance_matrix <- function(model_params, moment_func, dataset) {
  grad_matrix <- calculate_gradient_matrix(model_params, moment_func, dataset)
  grad_inverse <- solve(grad_matrix)
  moment_matrix <- moment_func(model_params, dataset)
  sample_size <- nrow(moment_matrix)
  omega_matrix <- t(moment_matrix) %*% moment_matrix / sample_size
  variance_matrix <- grad_inverse %*% omega_matrix %*% t(grad_inverse)
  return(variance_matrix / sample_size)
}

moment_conditions_xy_direction <- function(model_params, dataset) {
  instrument_z <- as.matrix(dataset$Z)
  endogenous_w <- as.matrix(dataset$W)
  covariate_v <- as.matrix(dataset$V)
  outcome_x <- as.matrix(dataset$X)
  outcome_y <- as.matrix(dataset$Y)
  
  fs_intercept <- model_params[1]
  fs_coef_z <- model_params[2]
  fs_coef_v1 <- model_params[3]
  fs_coef_interaction <- model_params[4]
  
  x_intercept <- model_params[5]
  x_coef_z <- model_params[6]
  x_coef_v <- model_params[7]
  x_coef_w <- model_params[8]
  
  y_intercept <- model_params[9]
  y_coef_z <- model_params[10]
  y_coef_v <- model_params[11]
  y_coef_w <- model_params[12]
  
  interaction_zv <- instrument_z * covariate_v
  fitted_w <- fs_intercept + fs_coef_z * instrument_z + 
    fs_coef_v1 * covariate_v + fs_coef_interaction * interaction_zv
  
  instruments_first_stage <- cbind(1, instrument_z, covariate_v, interaction_zv)
  instruments_structural <- cbind(1, instrument_z, covariate_v, fitted_w)
  
  residual_w <- endogenous_w - fitted_w
  moment_first_stage <- instruments_first_stage * as.vector(residual_w)
  
  residual_x <- outcome_x - (x_intercept + x_coef_z * instrument_z + 
                               x_coef_v * covariate_v + x_coef_w * fitted_w)
  moment_x_equation <- instruments_structural * as.vector(residual_x)
  
  residual_y <- outcome_y - (y_intercept + y_coef_z * instrument_z + 
                               y_coef_v * covariate_v + y_coef_w * fitted_w)
  moment_y_equation <- instruments_structural * as.vector(residual_y)
  
  moment_matrix <- cbind(moment_first_stage, moment_x_equation, moment_y_equation)
  return(moment_matrix)
}

construct_initial_values_xy <- function(dataset) {
  first_stage_fit <- lm(W ~ Z + V + I(Z * V), data = dataset)
  x_equation_fit <- lm(X ~ Z + V + first_stage_fit$fitted.values, data = dataset)
  y_equation_fit <- lm(Y ~ Z + V + first_stage_fit$fitted.values, data = dataset)
  
  fs_init <- coef(first_stage_fit)
  x_init <- coef(x_equation_fit)
  y_init <- coef(y_equation_fit)
  
  names(fs_init) <- NULL
  names(x_init) <- NULL
  names(y_init) <- NULL
  
  initial_params <- c(fs_init, x_init, y_init)
  return(initial_params)
}

moment_conditions_yx_direction <- function(model_params, dataset) {
  instrument_z <- as.matrix(dataset$Z)
  endogenous_w <- as.matrix(dataset$W)
  covariate_v <- as.matrix(dataset$V)
  outcome_x <- as.matrix(dataset$X)
  outcome_y <- as.matrix(dataset$Y)
  
  z_inst_intercept <- model_params[1]
  z_inst_coef_w <- model_params[2]
  z_inst_coef_v1 <- model_params[3]
  z_inst_coef_interaction <- model_params[4]
  
  x_intercept <- model_params[5]
  x_coef_z <- model_params[6]
  x_coef_v <- model_params[7]
  x_coef_w <- model_params[8]
  
  y_intercept <- model_params[9]
  y_coef_z <- model_params[10]
  y_coef_v <- model_params[11]
  y_coef_w <- model_params[12]
  
  interaction_wv <- endogenous_w * covariate_v
  fitted_z <- z_inst_intercept + z_inst_coef_w * endogenous_w + 
    z_inst_coef_v1 * covariate_v + z_inst_coef_interaction * interaction_wv
  
  instruments_first_stage <- cbind(1, endogenous_w, covariate_v, interaction_wv)
  instruments_structural <- cbind(1, fitted_z, covariate_v, endogenous_w)
  
  residual_z <- instrument_z - fitted_z
  moment_first_stage <- instruments_first_stage * as.vector(residual_z)
  
  residual_x <- outcome_x - (x_intercept + x_coef_z * fitted_z + 
                               x_coef_v * covariate_v + x_coef_w * endogenous_w)
  moment_x_equation <- instruments_structural * as.vector(residual_x)
  
  residual_y <- outcome_y - (y_intercept + y_coef_z * fitted_z + 
                               y_coef_v * covariate_v + y_coef_w * endogenous_w)
  moment_y_equation <- instruments_structural * as.vector(residual_y)
  
  moment_matrix <- cbind(moment_first_stage, moment_x_equation, moment_y_equation)
  return(moment_matrix)
}

construct_initial_values_yx <- function(dataset) {
  first_stage_fit <- lm(Z ~ W + V + I(W * V), data = dataset)
  x_equation_fit <- lm(X ~ first_stage_fit$fitted.values + V + W, data = dataset)
  y_equation_fit <- lm(Y ~ first_stage_fit$fitted.values + V + W, data = dataset)
  
  z_inst_init <- coef(first_stage_fit)
  x_init <- coef(x_equation_fit)
  y_init <- coef(y_equation_fit)
  
  names(z_inst_init) <- NULL
  names(x_init) <- NULL
  names(y_init) <- NULL
  
  initial_params <- c(z_inst_init, x_init, y_init)
  return(initial_params)
}

bidirectional_tsls <- function(dataset) {
  first_stage_w <- lm(W ~ Z + V + I(Z * V), data = dataset)
  x_equation_xy <- lm(X ~ Z + V + first_stage_w$fitted.values, data = dataset)
  y_equation_xy <- lm(Y ~ Z + V + first_stage_w$fitted.values, data = dataset)
  
  coef_x_on_z <- coef(x_equation_xy)["Z"]
  coef_y_on_z <- coef(y_equation_xy)["Z"]
  causal_effect_xy <- coef_y_on_z / coef_x_on_z
  
  first_stage_z <- lm(Z ~ W + V + I(W * V), data = dataset)
  x_equation_yx <- lm(X ~ first_stage_z$fitted.values + V + W, data = dataset)
  y_equation_yx <- lm(Y ~ first_stage_z$fitted.values + V + W, data = dataset)
  
  coef_x_on_w <- coef(x_equation_yx)["W"]
  coef_y_on_w <- coef(y_equation_yx)["W"]
  causal_effect_yx <- coef_x_on_w / coef_y_on_w
  
  return(c(
    causal_effect_xy = as.numeric(causal_effect_xy),
    causal_effect_yx = as.numeric(causal_effect_yx)
  ))
}

theoretical_method <- function(data) {
  bitsls_results <- bidirectional_tsls(data)
  effect_xy_bitsls <- bitsls_results["causal_effect_xy"]
  effect_yx_bitsls <- bitsls_results["causal_effect_yx"]
  
  initial_params_xy <- construct_initial_values_xy(data)
  gmm_result_xy <- optim(
    par = initial_params_xy,
    fn = calculate_gmm_objective,
    moment_func = moment_conditions_xy_direction,
    dataset = data,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  estimated_params_xy <- gmm_result_xy$par
  variance_matrix_xy <- estimate_variance_matrix(
    estimated_params_xy,
    moment_conditions_xy_direction,
    data
  )
  
  coef_y_on_z <- estimated_params_xy[10]
  coef_x_on_z <- estimated_params_xy[6]
  gradient_xy <- rep(0, 12)
  gradient_xy[10] <- 1 / coef_x_on_z
  gradient_xy[6] <- -coef_y_on_z / (coef_x_on_z^2)
  
  variance_effect_xy <- as.numeric(t(gradient_xy) %*% variance_matrix_xy %*% gradient_xy)
  std_error_xy_gmm <- sqrt(variance_effect_xy)
  
  initial_params_yx <- construct_initial_values_yx(data)
  gmm_result_yx <- optim(
    par = initial_params_yx,
    fn = calculate_gmm_objective,
    moment_func = moment_conditions_yx_direction,
    dataset = data,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  estimated_params_yx <- gmm_result_yx$par
  variance_matrix_yx <- estimate_variance_matrix(
    estimated_params_yx,
    moment_conditions_yx_direction,
    data
  )
  
  coef_x_on_w <- estimated_params_yx[8]
  coef_y_on_w <- estimated_params_yx[12]
  gradient_yx <- rep(0, 12)
  gradient_yx[8] <- 1 / coef_y_on_w
  gradient_yx[12] <- -coef_x_on_w / (coef_y_on_w^2)
  
  variance_effect_yx <- as.numeric(t(gradient_yx) %*% variance_matrix_yx %*% gradient_yx)
  std_error_yx_gmm <- sqrt(variance_effect_yx)
  
  return(list(
    beta_xy = effect_xy_bitsls,
    beta_yx = effect_yx_bitsls,
    se_xy = std_error_xy_gmm,
    se_yx = std_error_yx_gmm
  ))
}

empirical_method <- function(data, n_bootstrap = 200) {
  original_results <- bidirectional_tsls(data)
  effect_xy_original <- original_results["causal_effect_xy"]
  effect_yx_original <- original_results["causal_effect_yx"]
  
  num_observations <- nrow(data)
  bootstrap_estimates_xy <- numeric(n_bootstrap)
  bootstrap_estimates_yx <- numeric(n_bootstrap)
  
  for (boot_iter in 1:n_bootstrap) {
    boot_indices <- sample(1:num_observations, num_observations, replace = TRUE)
    boot_dataset <- data[boot_indices, ]
    
    boot_results <- bidirectional_tsls(boot_dataset)
    bootstrap_estimates_xy[boot_iter] <- boot_results["causal_effect_xy"]
    bootstrap_estimates_yx[boot_iter] <- boot_results["causal_effect_yx"]
  }
  
  std_error_xy_bootstrap <- sd(bootstrap_estimates_xy)
  std_error_yx_bootstrap <- sd(bootstrap_estimates_yx)
  
  return(list(
    beta_xy = effect_xy_original,
    beta_yx = effect_yx_original,
    se_xy = std_error_xy_bootstrap,
    se_yx = std_error_yx_bootstrap
  ))
}

process_single_iteration <- function(iter_idx, data_file, scenario, combo_key, n_key) {
  all_data <- readRDS(data_file)
  current_data <- all_data[[iter_idx]][[scenario]][[combo_key]][[n_key]]
  
  rm(all_data)
  gc(verbose = FALSE)
  
  theo_result <- theoretical_method(current_data)
  emp_result <- empirical_method(current_data, n_bootstrap = 200)
  
  return(list(
    iter = iter_idx,
    xy_theo_estimate = theo_result$beta_xy,
    xy_theo_se = theo_result$se_xy,
    yx_theo_estimate = theo_result$beta_yx,
    yx_theo_se = theo_result$se_yx,
    xy_emp_estimate = emp_result$beta_xy,
    xy_emp_se = emp_result$se_xy,
    yx_emp_estimate = emp_result$beta_yx,
    yx_emp_se = emp_result$se_yx
  ))
}

calculate_method_results <- function(estimates, ses, true_value) {
  ci_lowers <- estimates - 1.96 * ses
  ci_uppers <- estimates + 1.96 * ses
  
  coverage_count <- sum(ci_lowers <= true_value & ci_uppers >= true_value)
  coverage_rate <- coverage_count / length(estimates)
  
  mean_estimate <- mean(estimates)
  mean_se <- mean(ses)
  mean_ci_lower <- mean(ci_lowers)
  mean_ci_upper <- mean(ci_uppers)
  
  return(list(
    mean_estimate = mean_estimate,
    mean_se = mean_se,
    mean_ci_lower = mean_ci_lower,
    mean_ci_upper = mean_ci_upper,
    coverage_rate = coverage_rate,
    n_valid = length(estimates)
  ))
}

parallel_GMM_analysis_batch <- function(data_file, n_cores = 50) {
  true_beta_xy <- 0.5   
  true_beta_yx <- -0.5  
  
  scenarios <- c("a", "b", "c")
  combo_key <- "alpha_u_0.5_gamma_u_-0.5"
  n_keys <- c("N_1000", "N_2000", "N_5000")
  total_iters <- 200
  batch_size <- n_cores
  n_batches <- ceiling(total_iters / batch_size)
  
  all_results <- list()
  
  for(scenario in scenarios) {
    for(n_key in n_keys) {
      xy_theo_estimates <- numeric(total_iters)
      xy_theo_ses <- numeric(total_iters)
      yx_theo_estimates <- numeric(total_iters)
      yx_theo_ses <- numeric(total_iters)
      xy_emp_estimates <- numeric(total_iters)
      xy_emp_ses <- numeric(total_iters)
      yx_emp_estimates <- numeric(total_iters)
      yx_emp_ses <- numeric(total_iters)
      
      for(batch in 1:n_batches) {
        start_idx <- (batch - 1) * batch_size + 1
        end_idx <- min(batch * batch_size, total_iters)
        batch_indices <- start_idx:end_idx
        current_batch_size <- length(batch_indices)
        
        cl <- makeCluster(min(current_batch_size, n_cores))
        registerDoParallel(cl)
        
        clusterSetRNGStream(cl, 2)
        
        clusterExport(cl, c("calculate_gmm_objective", "calculate_gradient_matrix", 
                            "estimate_variance_matrix", 
                            "moment_conditions_xy_direction", "moment_conditions_yx_direction",
                            "construct_initial_values_xy", "construct_initial_values_yx",
                            "bidirectional_tsls", "theoretical_method", "empirical_method",
                            "process_single_iteration"))
        
        clusterEvalQ(cl, {
          library(numDeriv)
        })
        
        batch_results <- foreach(iter = batch_indices, 
                                 .packages = c("numDeriv")) %dopar% {
                                   
                                   result <- process_single_iteration(iter, data_file, scenario, combo_key, n_key)
                                   
                                   list(iter = result$iter,
                                        xy_theo_estimate = result$xy_theo_estimate,
                                        xy_theo_se = result$xy_theo_se,
                                        yx_theo_estimate = result$yx_theo_estimate,
                                        yx_theo_se = result$yx_theo_se,
                                        xy_emp_estimate = result$xy_emp_estimate,
                                        xy_emp_se = result$xy_emp_se,
                                        yx_emp_estimate = result$yx_emp_estimate,
                                        yx_emp_se = result$yx_emp_se)
                                 }
        
        stopCluster(cl)
        
        for(j in 1:length(batch_results)) {
          result <- batch_results[[j]]
          idx <- result$iter
          xy_theo_estimates[idx] <- result$xy_theo_estimate
          xy_theo_ses[idx] <- result$xy_theo_se
          yx_theo_estimates[idx] <- result$yx_theo_estimate
          yx_theo_ses[idx] <- result$yx_theo_se
          xy_emp_estimates[idx] <- result$xy_emp_estimate
          xy_emp_ses[idx] <- result$xy_emp_se
          yx_emp_estimates[idx] <- result$yx_emp_estimate
          yx_emp_ses[idx] <- result$yx_emp_se
        }
        
        gc(verbose = FALSE)
      }
      
      xy_theo_results <- calculate_method_results(xy_theo_estimates, xy_theo_ses, true_beta_xy)
      yx_theo_results <- calculate_method_results(yx_theo_estimates, yx_theo_ses, true_beta_yx)
      xy_emp_results <- calculate_method_results(xy_emp_estimates, xy_emp_ses, true_beta_xy)
      yx_emp_results <- calculate_method_results(yx_emp_estimates, yx_emp_ses, true_beta_yx)
      
      result_key <- paste0(scenario, "_", n_key)
      all_results[[result_key]] <- list(
        XY_theoretical = xy_theo_results,
        XY_empirical = xy_emp_results,
        YX_theoretical = yx_theo_results,
        YX_empirical = yx_emp_results
      )
    }
  }
  
  return(all_results)
}

summarize_results_excel <- function(results) {
  summary_data <- list()
  
  for(key in names(results)) {
    parts <- strsplit(key, "_")[[1]]
    scenario <- parts[1]
    n_size_raw <- parts[3]
    n_size <- as.numeric(gsub("N", "", n_size_raw))
    
    xy_theo <- results[[key]]$XY_theoretical
    xy_emp <- results[[key]]$XY_empirical
    yx_theo <- results[[key]]$YX_theoretical
    yx_emp <- results[[key]]$YX_empirical
    
    row_data <- data.frame(
      Scenario = scenario,
      Sample_Size = n_size,
      Direction = c("X→Y", "X→Y", "Y→X", "Y→X"),
      Method = c("Theoretical", "Empirical", "Theoretical", "Empirical"),
      Estimate = c(
        round(xy_theo$mean_estimate, 4),
        round(xy_emp$mean_estimate, 4),
        round(yx_theo$mean_estimate, 4),
        round(yx_emp$mean_estimate, 4)
      ),
      SE = c(
        round(xy_theo$mean_se, 4),
        round(xy_emp$mean_se, 4),
        round(yx_theo$mean_se, 4),
        round(yx_emp$mean_se, 4)
      ),
      CI_Lower = c(
        round(xy_theo$mean_ci_lower, 4),
        round(xy_emp$mean_ci_lower, 4),
        round(yx_theo$mean_ci_lower, 4),
        round(yx_emp$mean_ci_lower, 4)
      ),
      CI_Upper = c(
        round(xy_theo$mean_ci_upper, 4),
        round(xy_emp$mean_ci_upper, 4),
        round(yx_theo$mean_ci_upper, 4),
        round(yx_emp$mean_ci_upper, 4)
      ),
      Coverage_Rate = c(
        round(xy_theo$coverage_rate, 4),
        round(xy_emp$coverage_rate, 4),
        round(yx_theo$coverage_rate, 4),
        round(yx_emp$coverage_rate, 4)
      ),
      Valid_Count = c(
        xy_theo$n_valid,
        xy_emp$n_valid,
        yx_theo$n_valid,
        yx_emp$n_valid
      )
    )
    
    summary_data[[key]] <- row_data
  }
  
  final_df <- do.call(rbind, summary_data)
  rownames(final_df) <- NULL
  
  final_df <- final_df[order(final_df$Scenario, final_df$Sample_Size, final_df$Direction, final_df$Method), ]
  
  return(final_df)
}

run_parallel_GMM_batch <- function(data_file, n_cores = 50) {
  test_data <- readRDS(data_file)
  rm(test_data)
  gc(verbose = FALSE)
  
  results <- parallel_GMM_analysis_batch(data_file, n_cores = n_cores)
  
  summary_table <- summarize_results_excel(results)
  
  write.csv(summary_table, "Sim_Main_Variance.csv", row.names = FALSE)
  
  return(list(results = results, summary = summary_table))
}

final_results <- run_parallel_GMM_batch("Data/Sim_Data_Main.rds", n_cores = 50)