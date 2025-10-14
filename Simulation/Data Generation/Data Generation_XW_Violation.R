library(foreach)
library(doParallel)

set.seed(123)

generate_data <- function(N, alpha_u, gamma_u, beta_xw, n_time = 5000) {
  chunk_size <- ceiling(N/getDoParWorkers())
  
  alpha_0 <- 1;      gamma_0 <- -1;      eta_0 <- 1
  beta_yx <- -0.5;   beta_xy <- 0.5;     
  alpha_v <- 1;      gamma_v <- -1;      eta_v <- -1
  alpha_z <- 1;      gamma_w <- 2;       eta_u <- -1
  delta_u <- 1;      delta_v <- -0.5
  
  v_vec <- rnorm(N)
  u_vec <- exp(v_vec) + 2 * rbinom(N, size = 1, prob = 0.5) - 1
  
  epsilon_z <- rnorm(N)
  z_vec <- 1 + delta_u * u_vec + delta_v * v_vec + epsilon_z
  
  epsilon_x <- rnorm(N)
  epsilon_y <- rnorm(N)
  epsilon_w <- rnorm(N)
  
  iota <- 1 / (1 - beta_yx * beta_xy - beta_yx * gamma_w * beta_xw)
  
  alpha_0_star <- (alpha_0 + beta_yx * gamma_0 + beta_yx * gamma_w * eta_0) * iota
  gamma_0_star <- ((beta_xy + gamma_w * beta_xw) * alpha_0 + gamma_0 + gamma_w * eta_0) * iota
  eta_0_star <- (beta_xw * alpha_0 + beta_yx * beta_xw * gamma_0 + (1 - beta_xy * beta_yx) * eta_0) * iota
  
  alpha_v_star <- (alpha_v + beta_yx * gamma_v + beta_yx * gamma_w * eta_v) * iota
  gamma_v_star <- ((beta_xy + gamma_w * beta_xw) * alpha_v + gamma_v + gamma_w * eta_v) * iota
  eta_v_star <- (beta_xw * alpha_v + beta_yx * beta_xw * gamma_v + (1 - beta_xy * beta_yx) * eta_v) * iota
  
  alpha_z_star <- alpha_z * iota
  gamma_z_star <- (beta_xy + gamma_w * beta_xw) * alpha_z * iota
  eta_z_star <- beta_xw * alpha_z * iota
  
  alpha_u_star <- (alpha_u + beta_yx * gamma_u + beta_yx * gamma_w * eta_u) * iota
  gamma_u_star <- ((beta_xy + gamma_w * beta_xw) * alpha_u + gamma_u + gamma_w * eta_u) * iota
  eta_u_star <- (beta_xw * alpha_u + beta_yx * beta_xw * gamma_u + (1 - beta_xy * beta_yx) * eta_u) * iota
  
  epsilon_x_star <- (epsilon_x + beta_yx * epsilon_y + beta_yx * gamma_w * epsilon_w) * iota
  epsilon_y_star <- ((beta_xy + gamma_w * beta_xw) * epsilon_x + epsilon_y + gamma_w * epsilon_w) * iota
  epsilon_w_star <- (beta_xw * epsilon_x + beta_yx * beta_xw * epsilon_y + (1 - beta_xy * beta_yx) * epsilon_w) * iota
  
  simulation_results <- foreach(i = seq(1, N, by = chunk_size), 
                                .combine = rbind, 
                                .packages = c("stats")) %dopar% {
                                  current_chunk <- min(chunk_size, N - i + 1)
                                  chunk_indices <- i:(i + current_chunk - 1)
                                  
                                  results <- matrix(NA, nrow = current_chunk, ncol = 6)
                                  colnames(results) <- c("X", "Y", "Z", "W", "V", "U")
                                  
                                  for (j in 1:current_chunk) {
                                    idx <- chunk_indices[j]
                                    
                                    x_val <- alpha_0_star + alpha_v_star * v_vec[idx] + alpha_z_star * z_vec[idx] + alpha_u_star * u_vec[idx] + epsilon_x_star[idx]
                                    y_val <- gamma_0_star + gamma_v_star * v_vec[idx] + gamma_z_star * z_vec[idx] + gamma_u_star * u_vec[idx] + epsilon_y_star[idx]
                                    w_val <- eta_0_star + eta_v_star * v_vec[idx] + eta_z_star * z_vec[idx] + eta_u_star * u_vec[idx] + epsilon_w_star[idx]
                                    
                                    results[j,] <- c(x_val, y_val, z_vec[idx], w_val, v_vec[idx], u_vec[idx])
                                  }
                                  
                                  results
                                }
  
  as.data.frame(simulation_results)
}

N <- c(5000)
n_iterations <- 200
beta_xw_values <- seq(-0.5, 0.5, by = 0.1)

alpha_u_fixed <- 0.5
gamma_u_fixed <- -0.5

all_iterations <- vector("list", n_iterations)

n_cores <- parallel::detectCores() - 1
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)
clusterExport(cluster, c("generate_data"))

for(i in 1:n_iterations) {
  current_data <- list()
  
  combo_key <- paste0("alpha_u_", alpha_u_fixed, "_gamma_u_", gamma_u_fixed)
  current_data[[combo_key]] <- list()
  
  for (beta_xw_val in beta_xw_values) {
    beta_key <- paste0("beta_XW_", beta_xw_val)
    current_data[[combo_key]][[beta_key]] <- list()
    
    for (n_val in N) {
      results <- generate_data(
        N = n_val,
        alpha_u = alpha_u_fixed,
        gamma_u = gamma_u_fixed,
        beta_xw = beta_xw_val
      )
      
      current_data[[combo_key]][[beta_key]][[paste0("N_", n_val)]] <- results
    }
  }
  
  all_iterations[[i]] <- current_data
}

stopCluster(cluster)
saveRDS(all_iterations, file = "Data/Sim_Data_Supp_XW_Violation.rds")