library(foreach)
library(doParallel)

set.seed(123)

generate_data <- function(n_sample, scenario, alpha_u, gamma_u, R_w, R_z, n_time = 5000) {
  chunk_size <- ceiling(n_sample/getDoParWorkers())  
  
  alpha_0 <- 1;      gamma_0 <- -1
  beta_YX <- -0.5;   beta_XY <- 0.5
  alpha_v <- 1;      gamma_v <- -1
  alpha_z <- 1;      gamma_w <- 2
  gamma_z <- alpha_z * R_z; alpha_w <- gamma_w * R_w
  
  delta_u <- 1;      eta_u <- -1
  delta_v <- -0.5;   eta_v  <- 0.5
  
  v_vec <- rnorm(n_sample + n_time)
  u_vec <- exp(v_vec) + 2 * rbinom(n_sample + n_time, size = 1, prob = 0.5) - 1
  
  if (scenario == "a") {
    epsilon_z <- rnorm(n_sample + n_time)
    epsilon_w <- rnorm(n_sample + n_time)
  } else if (scenario == "b") {
    epsilon_z <- runif(n_sample + n_time, -1, 1)
    epsilon_w <- runif(n_sample + n_time, -1, 1)
  } else if (scenario == "c") {
    epsilon_z <- 2 * rbinom(n_sample + n_time, size = 1, prob = 0.5) - 1
    epsilon_w <- 2 * rbinom(n_sample + n_time, size = 1, prob = 0.5) - 1
  }
  
  z_vec <- 1 + delta_u * u_vec + delta_v * v_vec + epsilon_z
  w_vec <- 1 + eta_u * u_vec + eta_v * v_vec + epsilon_w
  
  error_y <- 2 * rbinom(n_sample + n_time, size = 1, prob = 0.5) - 1
  error_x <- 2 * rbinom(n_sample + n_time, size = 1, prob = 0.5) - 1
  
  simulation_results <- foreach(i = seq(1, n_sample, by = chunk_size), 
                                .combine = rbind, 
                                .packages = c("stats")) %dopar% { 
                                  current_chunk <- min(chunk_size, n_sample - i + 1)
                                  chunk_indices <- i:(i + current_chunk - 1)
                                  
                                  results <- matrix(NA, nrow = current_chunk, ncol = 6)
                                  colnames(results) <- c("Y", "X", "Z", "W", "V", "U")
                                  
                                  for (j in 1:current_chunk) {
                                    idx <- chunk_indices[j]
                                    y_series <- numeric(n_time)
                                    x_series <- numeric(n_time)
                                    
                                    y_series[1] <- rnorm(1)
                                    x_series[1] <- rnorm(1)
                                    
                                    for (t in 2:n_time) {
                                      x_series[t] <- alpha_0 + beta_yx * y_series[t-1] + alpha_v * v_vec[idx] + alpha_z * z_vec[idx] + alpha_w * w_vec[idx] + alpha_u * u_vec[idx] + error_x[idx]
                                      y_series[t] <- gamma_0 + beta_xy * x_series[t-1] + gamma_v * v_vec[idx] + gamma_z * z_vec[idx] + gamma_w * w_vec[idx] + gamma_u * u_vec[idx] + error_y[idx]
                                    }
                                    
                                    results[j,] <- c(y_series[n_time], x_series[n_time], z_vec[idx], w_vec[idx], v_vec[idx], u_vec[idx])
                                  }
                                  
                                  results
                                }
  
  as.data.frame(simulation_results)
}

scenarios <- c("a", "b", "c")
N <- c(5000)
n_iterations <- 200
alpha_u_val <- 0.5
gamma_u_val <- -0.5

R_combinations <- list()
for(i in 0:10) {  
  R_w_val <- -0.5 + i * 0.1  
  R_z_val <- 0.5 - i * 0.1   
  R_combinations[[i+1]] <- c(R_w = R_w_val, R_z = R_z_val)
}

all_iterations <- vector("list", n_iterations)

n_cores <- parallel::detectCores() - 1
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)
clusterExport(cluster, c("generate_data"))

for(i in 1:n_iterations) {
  current_data <- list()
  
  for (scenario in scenarios) {
    current_data[[scenario]] <- list()
    combo_key <- paste0("alpha_u_", alpha_u_val, "_gamma_u_", gamma_u_val)
    current_data[[scenario]][[combo_key]] <- list()
    
    for (r_combo in R_combinations) {
      R_w_val <- r_combo["R_w"]
      R_z_val <- r_combo["R_z"]
      r_combo_key <- paste0("R_w_", R_w_val, "_R_z_", R_z_val)
      current_data[[scenario]][[combo_key]][[r_combo_key]] <- list()
      
      for (n_val in N) {
        results <- generate_data(
          n_sample = n_val,
          scenario = scenario,
          alpha_u = alpha_u_val,
          gamma_u = gamma_u_val,
          R_w = R_w_val,
          R_z = R_z_val
        )
        
        current_data[[scenario]][[combo_key]][[r_combo_key]][[paste0("N_", n_val)]] <- results
      }
    }
  }
  
  all_iterations[[i]] <- current_data
}

stopCluster(cluster)
saveRDS(all_iterations, file = "Data/Sim_Data_Sensitivity.rds")